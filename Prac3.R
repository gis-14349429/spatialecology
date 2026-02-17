
setwd("C:/Users/s16188ss/OneDrive - The University of Manchester/GEOG71922/WEEK3/spatialecology")

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

names(allEnv)=c("broadleaf","urban","elev")
