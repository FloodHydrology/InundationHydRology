###################################################################################
#Name: Synthetic Landscape
#Coder: C. Nathan Jones
#Date: 4/11/2016
#Purpose: Develop synthetic landscape (z) for PC Model Comparison
##################################################################################

####################################################################################
#Step 1: Setup Workspace
####################################################################################
#Clear Memory
rm(list=ls(all=TRUE))

#Define parent directory
wd<-"C:\\Users/cnjones/Google Drive/Existing Conditions/DEM_Inundate/"

#Set Working Directory
setwd("M:\\McLaughlin_Lab/Jones/workspace/ConnModWorkshop/Model_Comparison/")

#Download packages (use install.packages function if not already downloaded)
library(raster)
library(sp)
library(gstat)
library(rgdal)
library(maptools)
library(foreign)
library(rgeos)

####################################################################################
#Step 2: Create Function to Create Synthetic DEM
####################################################################################
#Model Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#lanscape size (Note, it will be a square)
area<-10000 #m^2 

#Size of grid cell (e.g. length of one side)
grid<-1 #m

#landscape slope
slope<-0.01 #m/m

#number of wetlands
n.wetlands<-16

#proportion of wetlands on landscape
wetland_upland_proportion<-0.2

#wetland depth
invert<-1 #m

#Watershed depth
watershed_depth<-0.5#m

#Create Raster~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x<-matrix(0, nrow=sqrt(area)/grid, ncol=sqrt(area)/grid)
y<-matrix(0, nrow=sqrt(area)/grid, ncol=sqrt(area)/grid)
z<-matrix(0, nrow=sqrt(area)/grid, ncol=sqrt(area)/grid)

#Add coordinates
for(i in 1:(sqrt(area)/grid)){
  x[,i]<-i*grid-(grid/2)
  y[i,]<-i*grid-(grid/2)
  }

#Add elevation data (valley slope)
interp<-approxfun(data.frame(c(1,sqrt(area)/grid), c(100, 100-slope*sqrt(area))))
for(i in 1:(sqrt(area)/grid)){z[i,]<-interp(i)}

#Create Raster
dem <-raster(
              z,
              xmn=range(x)[1], xmx=range(x)[2],
              ymn=range(y)[1], ymx=range(y)[2]
)

#Create Wetands shapes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Define wetland and upland areas
watershed_area<-area/n.wetlands
wetland_area<-wetland_upland_proportion*watershed_area

#Define coordinates of wetlands (center)
pnts<-matrix(seq(watershed_area^0.5,area^0.5,watershed_area^0.5))-(watershed_area^.5)/2
pnts<-data.frame(rep(pnts,4), rep(pnts, each=4))
  colnames(pnts)<-c("x","y")
pnts$z<-extract(dem, pnts)
pnts.shp<-SpatialPoints(pnts)

#Create function to define the coordinates of the wetland perimeter
circle.fun<-function(df, wetid){
  radius<-(wetland_area/pi)^.5
  circle<-seq(0, 2 * pi, length.out = 100)
  circle<-cbind(df$x + radius * sin(circle), df$y + radius * cos(circle))
  circle
}

#Create Wetlands Polygon
wetlands.shp<-data.frame(circle.fun(pnts[1,]))
wetlands.shp<-Polygon(wetlands.shp)
wetlands.shp<-Polygons(list(wetlands.shp),1)
wetlands.shp<-SpatialPolygons(list(wetlands.shp))
for(i in 2:n.wetlands){
  df<-data.frame(circle.fun(pnts[i,]))
  df<-Polygon(df)
  df<-Polygons(list(df),i)
  df<-SpatialPolygons(list(df))
  wetlands.shp<-gUnion(wetlands.shp, df)
}
wetlands.shp<-disaggregate(wetlands.shp)

#Add z value to wetlands.shp
wetlands.shp$z<-pnts$z[over(pnts.shp,wetlands.shp)]

#Convert wetlands.shp to raster
wetlands.grd <- dem*0
extent(wetlands.grd) <- extent(dem)
wetlands.grd<- rasterize(wetlands.shp, wetlands.grd, wetlands.shp$z)

#Incoorporate wetlands into raster~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Add wetland invert to landscape
delete.grd<-wetlands.grd*0
delete.grd[is.na(delete.grd)] <- 1
invert.grd<-wetlands.grd-invert
invert.grd[is.na(invert.grd)] <- 0
dem_nowatershed<-dem*delete.grd+invert.grd

#Export Raster
writeRaster(dem_nowatershed,"dem_nowatershed.asc" )

#Incooporate Watershed Boundaries into DEM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#number of lines in one direction
boundaries<-c(0,seq(watershed_area^0.5,area^0.5,watershed_area^0.5))
n.lines<-length(boundaries)

#Create horizontal lines
boundaries.shp<-matrix(c(0, boundaries[1], sqrt(area), boundaries[1]), nrow=2, ncol=2, byrow=T)
boundaries.shp<-Line(boundaries.shp)
boundaries.shp<-Lines(list(boundaries.shp), ID="a")
boundaries.shp<-SpatialLines(list(boundaries.shp))
for(i in 2:n.lines){
  df<-matrix(c(0, boundaries[i], sqrt(area), boundaries[i]), nrow=2, ncol=2, byrow=T)
  df<-Line(df)
  df<-Lines(list(df), ID="a")
  df<-SpatialLines(list(df))
  boundaries.shp<-gUnion(boundaries.shp, df)
}

#Create vertical lines
vert.shp<-matrix(c(boundaries[1], 0,boundaries[1],sqrt(area)), nrow=2, ncol=2, byrow=T)
vert.shp<-Line(vert.shp)
vert.shp<-Lines(list(vert.shp), ID="a")
vert.shp<-SpatialLines(list(vert.shp))
boundaries.shp<-gUnion(boundaries.shp, vert.shp)
for(i in 2:n.lines){
  df<-matrix(c(boundaries[i],0,boundaries[i], sqrt(area)), nrow=2, ncol=2, byrow=T)
  df<-Line(df)
  df<-Lines(list(df), ID="a")
  df<-SpatialLines(list(df))
  boundaries.shp<-gUnion(boundaries.shp, df)
}

#Convert watershed bounderies to raster
extent(dem)<-extent(boundaries.shp)
boundaries.grd <- dem*0
extent(boundaries.grd) <- extent(boundaries.shp)
boundaries.grd<- rasterize(boundaries.shp, boundaries.grd)

#Obtain elevation along boundaries from DEM
bound_ele.grd<-boundaries.grd
bound_ele.grd<-bound_ele.grd*dem

#Create wetland invert grid
ws_invert.grd<-wetlands.grd-watershed_depth

#Convert boundary and wetland invert to points
bound_ele.pnt<-rasterToPoints(bound_ele.grd)
ws_invert.pnt<-rasterToPoints(ws_invert.grd)
interp.pnts<-data.frame(rbind(ws_invert.pnt,bound_ele.pnt))
  colnames(interp.pnts)<-c("x","y","z")

#Interpolate to watershed invert
blank.grd<-wetlands.grd*0
IDW.grd<-gstat(id="layer", formula=z~1, locations=~x+y, data=interp.pnts, nmax=7, set=list(idp=4.2))
IDW.grd<-interpolate(blank.grd,IDW.grd)

#Add wetland
delete.grd<-wetlands.grd*0
delete.grd[is.na(delete.grd)] <- 1
invert.grd<-wetlands.grd-invert
invert.grd[is.na(invert.grd)] <- 0
dem_watershed<-IDW.grd*delete.grd+invert.grd

#Export Raster
writeRaster(dem_watershed,"dem_watershed.asc" )




