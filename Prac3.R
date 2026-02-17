
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
