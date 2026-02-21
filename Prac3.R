
setwd("C:/Users/User/OneDrive - The University of Manchester/GEOG71922/WEEK3/spatialecology")

#Install packages

install.packages(c("dismo","terra","spatstat","cowplot","ggplot2","precrec","glmnet","maxnet"))


#We also need to load the "MLR" package for this practical

install.packages("mlr")
#Load libraries 
library(terra) #for spatial data
library(sf) #data frames with geometry
library(spatstat) #for point process modelling and converting between raster and pixel image objects
#Function to convert raster to images for spatstat. Takes one argument "im" that should be a raster object

raster.as.im = function(im) {
  #get the resolution (cell size of the raster)
  r = raster::res(im)[1]
  #get the origin (bottom left corner of the raster/image)
  orig = ext(im)[c(1,3)]
  #set the coordinates of the columns which is just a series of number from zero increasing by 100 metres (the resolution of the raster) for every cell along the rows and columns.
  xx = orig[1] + seq(from=0,to=(ncol(im) - 1)*100,by=r)
  #set the coordinates of the columns
  yy = orig[2] + seq(from=0,to=(nrow(im) - 1)*100,by=r)
  
  #now build a single matrix with the cell values and dimension we want - note that we reverse the rows in the matrix by setting nrow(im):1. This is just because spatstat reads images in reverse order to how raster object are usually organised (so the below code ensures our image is not upside-down when we come to analyse it)
  mat=matrix(raster::values(im), ncol = ncol(im), 
             nrow = nrow(im), byrow = TRUE)[nrow(im):1, ]
  return(spatstat.geom::im(mat, 
                           xcol = xx, yrow = yy))
}

#load in the red squirrel data
sciurus=read.csv("Sciurus.csv")


#subset to cases with coordinate uncertainty <100 metres
sciurus=sciurus[sciurus$Coordinate.uncertainty_m<100,]

#make spatial points layer

#create crs object
sciurus.latlong=data.frame(x=sciurus$Longitude,y=sciurus$Latitude)

#Use coordinates object to create our spatial points object

sciurus.sp=st_as_sf(sciurus.latlong,coords=c("x","y"),crs="epsg:4326")
#First set the extent to the study area

scot=st_read('scotSamp.shp')

#load in the land cover map and then clip to the polygon

LCM=rast("LCMUK.tif")

#crop to the extent of the study area plus a little more (because we will lose a small amount of data in the next step)
LCM=crop(LCM,st_buffer(scot, dist= 1000))

#aggregate LCM raster
LCM=aggregate(LCM$LCMUK_1,fact=4,fun="modal")
#project squirrel data
sciurus.sp=st_transform(sciurus.sp,crs(LCM))
#now crop our points to the study area
sciurusFin=sciurus.sp[scot,]

#finally, mask the LCM to this boundary

LCM=crop(LCM,scot,mask=TRUE)

#inspect
plot(LCM)
plot(sciurusFin$geometry,add=T)
#access levels of the raster by treating them as categorical data ('factors' in R)

LCM=as.factor(LCM$LCMUK_1)


#create an vector object called reclass
reclass = c(0,1,rep(0,20))

# combine with the LCM categories into a matrix of old and new values.
RCmatrix=cbind(levels(LCM)[[1]],reclass)

RCmatrix=RCmatrix[,2:3]

#apply function to make sure new columns are numeric (here the "2" specifies that we want to apply the as.numeric function to columns, where "1" would have specified rows)
RCmatrix=apply(RCmatrix,2,FUN=as.numeric)
#Use the classify() function to asssign new values to LCM with our reclassification matrix

broadleaf=classify(LCM, RCmatrix)
#neighbourhood weights matrix to sum all available resources for each cell

#get number of pixels needed to cover the 1800 meter radius (round() just rounds this to the nearest interger - because you can't have 0.5 pixels)

nPix=round(1800/res(LCM)[1])

#next, you need to double this number (for the distance in two direction i.e. nPix is the radius and we want the diameter) and add one (because a weights matrix always needs to be odd so that you have a central cell and an equal number of grid cells either side).

nPix=(nPix*2)+1
#buiild weights matrix
weightsMatrix=matrix(1:nPix^2,nrow=nPix,ncol=nPix)
#get focal cell 
x=ceiling(ncol(weightsMatrix)/2)
y=ceiling(nrow(weightsMatrix)/2)


focalCell=weightsMatrix[x,y]
indFocal=which(weightsMatrix==focalCell,arr.ind = TRUE)
#compute distances
distances=list()

for(i in 1:nPix^2){
  ind.i=which(weightsMatrix==i,arr.ind=T)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distances[[i]]=dist.i
  
}
#add distance values to the weights matrix
weightsMatrix[]=unlist(distances)

#set cells outside search radius to NA
weightsMatrix[weightsMatrix>1800]=NA
#plot weights matrix
plot(rast(weightsMatrix))
#normalise the weights matrix by dividing all cell values by the number of cells. 

weightsMatrixNorm=weightsMatrix

weightsMatrixNorm[!is.na(weightsMatrixNorm)]=1/length(weightsMatrixNorm[!is.na(weightsMatrixNorm)])

#test to see for yourself
sum(weightsMatrixNorm,na.rm=T)
plot(rast(weightsMatrixNorm))
#sum neighbourhood values from all surrounding cells

lcm_wood_1800=focal(broadleaf,w=weightsMatrixNorm,fun="sum")

plot(lcm_wood_1800)
#create an vector object called reclassUrban which is zero for all classes except tghe two urban classes in the LCM
reclassUrban = c(rep(0,19),1,1)

# combine with the LCM categories into a matrix of old and new values.
RCmatrixUrban= cbind(levels(LCM)[[1]],reclass)

RCmatrixUrban=RCmatrixUrban[,2:3]

#apply function to make sure new columns are numeric (here the "2" specifies that we want to apply the as.numeric function to columns, where "1" would have specified rows)
RCmatrixUrban=apply(RCmatrixUrban,2,FUN=as.numeric)
#Use the reclassify() function to asssign new values to LCM with our reclassification matrix

urban = classify(LCM, RCmatrixUrban)

#neighbourhood weights matrix to sum all available resources for each cell

#get number of picels needed to cover the 2300 metre radius for the urban class (round() just rounds this to the nearest interger - because you can't have 0.5 pixels)

nPixUrban=round(2300/res(LCM)[1])

#next, you need to double this number (one for rows and one for columns) and add one (because a weights matrix always needs to be odd so that you have a central cell and an equal number of grid cells either side).

nPixUrban=(nPixUrban*2)+1

#buiild weights matrix
weightsMatrixUrban=matrix(1:nPixUrban^2,nrow=nPixUrban,ncol=nPixUrban)

#get focal cell 
x=ceiling(ncol(weightsMatrixUrban)/2)
y=ceiling(nrow(weightsMatrixUrban)/2)


focalCell=weightsMatrixUrban[x,y]


indFocal=which(weightsMatrixUrban==focalCell,arr.ind = TRUE)

#compute distances
distancesUrban=list()

for(i in 1:nPixUrban^2){
  ind.i=which(weightsMatrixUrban==i,arr.ind=T)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distancesUrban[[i]]=dist.i
  
}


#add distance values to the weights matrix
weightsMatrixUrban[]=unlist(distancesUrban)

#set cells outside search radius to NA
weightsMatrixUrban[weightsMatrixUrban>2300]=NA

#normalise the weights matrix by dividing all cell values by the number of cells. 
weightsMatrixUrban[!is.na(weightsMatrixUrban)]=1/length(weightsMatrixUrban[!is.na(weightsMatrixUrban)])

#sum urban class from all surrounding cells

lcm_urban_2300=focal(urban,w=weightsMatrixUrban,fun="sum")
demScot=rast('demScotland.tif')
demScot=terra::resample(demScot,lcm_wood_1800)

#inspect
plot(demScot)
#stack the covariate layers together
allEnv=c(lcm_wood_1800,lcm_urban_2300,demScot)

all.cov=st_drop_geometry(all.cov)

names(allEnv)=c("broadleaf","urban","elev")
set.seed(11)

#sample background - one point for every cell (9775)
back = spatSample(allEnv,size=2000,as.points=TRUE,method="random",na.rm=TRUE) 
back=back[!is.na(back$broadleaf),]

back=st_as_sf(back,crs="EPSG:27700")
# get environmental covariates at presence locations
eP=terra::extract(allEnv,sciurusFin)

#bind together the presence data using cbind() which binds together objects by column (i.e. with different columns but the same number of rows)
Pres.cov=st_as_sf(cbind(eP,sciurusFin))

Pres.cov$Pres=1


#Remove the first column which is just an ID field.
Pres.cov=Pres.cov[,-1]

#get coordinates for spatial cross-validation later
coordsPres=st_coordinates(Pres.cov)

#drop geometry column using st_drop_geometry()
Back.cov=st_as_sf(data.frame(back,Pres=0))


#get coordinates of background points for cross validation later
coordsBack=st_coordinates(back)

#combine
coords=data.frame(rbind(coordsPres,coordsBack))

#assign column names
colnames(coords)=c("x","y")
#combine pres and background
all.cov=rbind(Pres.cov,Back.cov)

#add coordinates
all.cov=cbind(all.cov,coords)


#remove any NAs
all.cov=na.omit(all.cov)
all.cov=st_drop_geometry(all.cov)
library(mlr)
task=all.cov
head(all.cov)
task$Pres=as.factor(task$Pres)

task = makeClassifTask(data = task[,c(1:4)], target = "Pres",
                       positive = "1", coordinates = task[,5:6])
################ Binomial (logistic regression)

#use the make learner function to build the model approach. Fix factors prediction is set to TRUE here because our outcome is a factor (i.e. categorical: 0-1)
lrnBinomial = makeLearner("classif.binomial",
                          predict.type = "prob",
                          fix.factors.prediction = TRUE)



##########Random Forest

install.packages("ranger")
library(ranger)
lrnRF = makeLearner("classif.ranger",
                    predict.type = "prob",
                    fix.factors.prediction = TRUE)
#set up resampling strategy (non-spatial cross-validation)
perf_levelCV = makeResampleDesc(method = "RepCV", predict = "test", folds = 5, reps = 5)

#set up resampling strategy for spatial cross-validation

perf_level_spCV = makeResampleDesc(method = "SpRepCV", folds = 5, reps = 5) #sampling strategy to run five fold re-sampling five times
#Binomial conventional cross validation (K fold)


cvBinomial = resample(learner = lrnBinomial, task =task,
                      resampling = perf_levelCV, 
                      measures = mlr::auc,
                      show.info = FALSE)

print(cvBinomial)
#create Spatial Resampling Plots
plots = createSpatialResamplingPlots(task,resample=cvBinomial,
                                     crs=crs(allEnv),datum=crs(allEnv),color.test = "red",point.size = 1)

library(cowplot)
#use the cowplot function to plot all folds out in a grid
cowplot::plot_grid(plotlist = plots[["Plots"]], ncol = 3, nrow = 2,
                   labels = plots[["Labels"]])
##Binomial spatial cross validation

sp_cvBinomial = resample(learner = lrnBinomial, task =task,
                         resampling = perf_level_spCV, 
                         measures = mlr::auc,
                         show.info = FALSE)

print(sp_cvBinomial)
#make partition plots
plotsSP = createSpatialResamplingPlots(task,resample=sp_cvBinomial,
                                       crs=crs(allEnv),datum=crs(allEnv),color.test = "red",point.size = 1)

#use the cowplot function to plot all folds out in a grid
cowplot::plot_grid(plotlist = plotsSP[["Plots"]], ncol = 3, nrow = 2,
                   labels = plotsSP[["Labels"]])
###### Random Forest evaluation


##random sampling cross-validation

cvRF = resample(learner = lrnRF, task =task,
                resampling = perf_levelCV, 
                measures = mlr::auc,
                show.info = FALSE)

print(cvRF)
#spatial partitioning cross-validation

sp_cvRF = resample(learner = lrnRF, task =task,
                   resampling = perf_level_spCV, 
                   measures = mlr::auc,
                   show.info = FALSE)

print(sp_cvRF)
getParamSet(lrnRF)
paramsRF = makeParamSet(
  makeIntegerParam("mtry",lower = 1,upper = 3),
  makeIntegerParam("min.node.size",lower = 1,upper = 20),
  makeIntegerParam("num.trees",lower = 100,upper = 500)
)
# specifying random parameter value search

tune_level = makeResampleDesc(method = "SpCV", iters = 5)

ctrl = makeTuneControlRandom(maxit = 50)
tuned_RF = tuneParams(learner = lrnRF,
                      task = task,
                      resampling = tune_level,
                      measures = mlr::auc,
                      par.set = paramsRF,
                      control = ctrl,
                      show.info = FALSE)

print(tuned_RF)
# make a grid for spatial blocking

area_grid = st_make_grid(sciurusFin, c(50000, 50000), what = "polygons", square = T)
area_grid_sf=st_as_sf(area_grid)

#the original grid from st_make_grid() is contained in a list so we can use the position of each square of that list as its ID.
area_grid_sf$grid_id=1:length(lengths(area_grid))
# plot
plot(area_grid_sf$x)

#check all points in correct place
plot(sciurusFin$geometry,add=T)
# set folds
folds=area_grid_sf$grid_id

dataPoints=st_as_sf(all.cov,coords = c("x","y"))

st_crs(dataPoints)=crs(area_grid_sf)
library(dismo)
library(maxnet)
library(glmnet)

#set number of folds to use
folds=5

#partition presence and background data and assign to folds using the kfold() function.

Pres.cov=all.cov[all.cov$Pres==1,]
Back.cov=all.cov[all.cov$Pres==0,]


kfold_pres = kfold(Pres.cov, folds)
kfold_back = kfold(Back.cov, folds)

eMax=list()

for (i in 1:folds) {
  train = Pres.cov[kfold_pres!= i,]
  test = Pres.cov[kfold_pres == i,]
  backTrain=Back.cov[kfold_back!=i,]
  backTest=Back.cov[kfold_back==i,]
  dataTrain=rbind(train,backTrain)
  dataTest=rbind(test,backTest)
  maxnetMod=maxnet(dataTrain$Pres, dataTrain[,1:3])
  
  eMax[[i]] = evaluate(p=dataTest[ which(dataTest$Pres==1),],a=dataTest[which(dataTest$Pres==0),],maxnetMod)
}

#print the result

aucMax = sapply(eMax, function(x){slot(x, 'auc')} )

print(mean(aucMax))
#load precrec for more versatile model evaluation methods
library(precrec)

folds=area_grid_sf$grid_id

maxEvalList=list()
for (i in folds) {
  
  
  
  #select fold - select all folds which are not 'i' to train the model
  gridTrain=subset(area_grid_sf,area_grid_sf$grid_id!=i)
  
  #spatially subset the training data
  train=data.frame(st_drop_geometry(st_intersection(gridTrain, dataPoints)))
  
  #spatially subset the test data 
  gridTest=subset(area_grid_sf,area_grid_sf$grid_id==i)
  
  
  test=data.frame(st_drop_geometry(st_intersection(gridTest, dataPoints)))  
  
  
  #build model
  train
  maxnetMod=maxnet(train$Pres, train[1:3],
                   classes="lq")   
  
  
  #make prediction
  pred=predict(maxnetMod, test,type="cloglog")
  
  #use precrec to get roc curves
  precrec_proc=evalmod(scores = pred,labels = test$Pres,mode = "prcroc")
  
  
  #get auc value
  modauc=precrec::auc(precrec::evalmod(scores = pred, 
                                       labels = test$Pres))
  
  
  #add to list
  maxEvalList[[i]]=modauc$aucs[1]
  
  #track
  print(i)
  
}
mean(unlist(maxEvalList))
#Specify the model
glm.sciurus=glm(Pres~broadleaf+urban+elev,binomial(link='logit'),
                data=all.cov)

#predict and inspect the output
prGLM=predict(allEnv,glm.sciurus,type="response")

#plot
plot(prGLM)
#build new data frame based on mean of elev and urban but varying values for broadleaf. 
glmNew=data.frame(broadleaf=seq(0,max(all.cov$broadleaf),length=1000),
                  elev=mean(all.cov$elev),
                  urban=mean(all.cov$urban))


# use type = "response" for probability-scale predictions and chose to return the standard error of the prediction (se.fit=TRUE)   
preds = predict(glm.sciurus, newdata = glmNew, type = "response", se.fit = TRUE)
glmNew$fit = preds$fit
glmNew$se = preds$se.fit

head(glmNew)
library(ggplot2)
#use ggplot to plot the fitted values against the broadleaf variable. Use the ribbon() function in ggplot to add a 95% confidence interval.


ggplot(glmNew, aes(x = broadleaf, y = fit)) +
  
  geom_ribbon(data = glmNew, aes(y = fit, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "blue", alpha = 0.3) +
  geom_line(data = glmNew, aes(y = fit)) 
#build new data frame based on mean of elev and broadleaf but varying values for urban. 
glmNewUrban=data.frame(urban=seq(0,max(all.cov$urban),length=1000),
                       elev=mean(all.cov$elev),
                       broadleaf=mean(all.cov$broadleaf))


# use type = "response" for probability-scale predictions    
predUrban = predict(glm.sciurus, newdata = glmNewUrban, type = "response", se.fit = TRUE)
glmNewUrban$fit = predUrban$fit
glmNewUrban$se = predUrban$se.fit
ggplot(glmNewUrban, aes(x = urban, y = fit)) +
  
  geom_ribbon(data = glmNewUrban, aes(y = fit, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "blue", alpha = 0.3) +
  geom_line(data = glmNewUrban, aes(y = fit)) 
#build new data frame based on mean of broaleaf and urban but varying values for elev. 
glmNewElev=data.frame(elev=seq(0,max(all.cov$elev),length=1000),
                      urban=mean(all.cov$urban),
                      broadleaf=mean(all.cov$broadleaf))


# use type = "response" for probability-scale predictions    
predElev = predict(glm.sciurus, newdata = glmNewElev, type = "response", se.fit = TRUE)
glmNewElev$fit = predElev$fit
glmNewElev$se = predElev$se.fit
ggplot(glmNewElev, aes(x = elev, y = fit)) +
  
  geom_ribbon(data = glmNewElev, aes(y = fit, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "blue", alpha = 0.3) +
  geom_line(data = glmNewElev, aes(y = fit)) 
#Specify the glm model using the appropriate polynomial term for boradleaf

glm.sciurus=glm(Pres~poly(broadleaf,3)+urban+elev,binomial(link='logit'),
                data=all.cov)
#load spatstat
library(spatstat)

#convert broadleaf object
broadleafIm=raster.as.im(raster(allEnv$broadleaf))

#convert urban object
urbanIm=raster.as.im(raster(allEnv$urban))

#convert elevation object
elevIm=raster.as.im(raster(allEnv$elev))
#create study window using the broadleafIm image object as a template
window.poly=as.owin(urbanIm)

plot(window.poly)#inspect
#get coordinates from our sciurusFin object to create point pattern
sciurusCoords=st_coordinates(sciurusFin) 

#create point pattern ('ppp') object with the ppp() function
pppSciurus=ppp(sciurusCoords[,1],sciurusCoords[,2], window = window.poly)
#inspect
plot(allEnv$broadleaf)

plot(pppSciurus,add=T)
#use as.ppp() to remove points outside the window
pppSciurus=as.ppp(pppSciurus)
#rescale from m to km
pppSciurus=spatstat.geom::rescale(pppSciurus, 1000)

#rescale image objects to conform with the point pattern 

broadleafIm=spatstat.geom::rescale(broadleafIm,1000)
elevIm=spatstat.geom::rescale(elevIm,1000)
urbanIm=spatstat.geom::rescale(urbanIm,1000)
Kcsr=envelope(pppSciurus,Kest,nsim=39,VARIANCE=T,nSD=1,global =TRUE)
plot(Kcsr,shade=c("hi","lo"),legend=T) #check whether our data fall within this envelope
ndTry=seq(100,1000,by=100)


for(i in ndTry){
  #build quadrature scheme
  Q.i=quadscheme(pppSciurus,method = "grid",nd=i)
  #run a simple model
  fit.i=ppm(Q.i~broadleafIm+elevIm+urbanIm)
  print(i)
  #print the AIC statistic
  print(AIC(fit.i))
}
Q=quadscheme(pppSciurus,method = "grid",nd=900)
plot(rhohat(pppSciurus,broadleafIm)) # broadleaf
plot(rhohat(pppSciurus,elevIm)) # elevation
plot(rhohat(pppSciurus,urbanIm)) # urban
firstPPMod=ppm(Q~poly(broadleafIm,3)+poly(elevIm,2)+poly(urbanIm,2)+x+y) #point process model based on polynomial terms for the three covariates and the projected coordinates
firstModEnv=envelope(firstPPMod,Kest,nsim=39,VARIANCE=TRUE,nSD=1,global =TRUE)
plot(firstModEnv)
thomasMod=kppm(Q~poly(broadleafIm,3)+poly(elevIm,2)+poly(urbanIm,2)+x+y,"Thomas") # apply Thomas cluster process to the data.
thomasEnv=envelope(thomasMod,Kinhom,nsim=39,VARIANCE=TRUE,nSD=1,global=TRUE)
plot(thomasEnv)
plot(roc(thomasMod))
auc.kppm(thomasMod)
prPPMod=predict(thomasMod)

plot(prPPMod)
plot(rast(prPPMod))
ppmRast=rast(prPPMod,)

writeRaster(ppmRast,"ppmSciurus.tif")
