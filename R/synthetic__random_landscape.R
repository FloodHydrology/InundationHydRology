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
library(actuar)

####################################################################################
#Step 2: Define Inputs
####################################################################################
#Model Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#lanscape size (Note, it will be a square)
area<-1000 #ha  

#Size of cell_size cell (e.g. length of one side)
cell_size<-3 #m

#wetland depth
invert<-1 #m

#landscape slope
slope<-0.001 #m/m

####################################################################################
#Step 3: Create Raster
####################################################################################
#Convert units to meters
area<-area*10000 #convert to m^2

#Define coordinates
x<-matrix(0, nrow=sqrt(area)/cell_size, ncol=sqrt(area)/cell_size)
y<-matrix(0, nrow=sqrt(area)/cell_size, ncol=sqrt(area)/cell_size)
z<-matrix(0, nrow=sqrt(area)/cell_size, ncol=sqrt(area)/cell_size)

#Add coordinates
for(i in 1:(sqrt(area)/cell_size)){
  x[,i]<-i*cell_size-(cell_size/2)
  y[i,]<-i*cell_size-(cell_size/2)
  }

#Add elevation data (valley slope)
interp<-approxfun(data.frame(c(1,sqrt(area)/cell_size), c(100, 100-slope*sqrt(area))))
for(i in 1:(sqrt(area)/cell_size)){z[i,]<-interp(i)}

#Create Raster
dem <-raster(
              z,
              xmn=range(x)[1], xmx=range(x)[2],
              ymn=range(y)[1], ymx=range(y)[2]
)

#clean up workspace
remove(list=c("x","y","z","i","slope"))

####################################################################################
#Step 4: Create Wetlands
####################################################################################
#Determine locationa of wetlands~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#determine number and size of wet
wetland_area<-area*0.12 #From Van Meter and Basu 2015 (Des Moines Lobe)

#Iterate through wetland scenarios until correct area
set.seed(1)
dif<-area
while(abs(dif)>(50)){
  n.wetlands<-round(runif(1, 5, 500), digits = 0)
  pnts<-rplcon(n.wetlands, 10^3, 1.67)
  dif<-wetland_area-sum(pnts)
}

#Randomly select location of each wetland
pnts<-data.frame(coordinates(sampleRandom(dem,n.wetlands, sp=T)), pnts)
  colnames(pnts)<-c("x","y", "area_m2")

#add z data and create shp
pnts$z<-extract(dem, pnts[,c("x","y")])
pnts.shp<-SpatialPoints(pnts)

#Create Wetland shapes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create function to define the coordinates of the wetland perimeter
circle.fun<-function(df){
  radius<-(df$area_m2/pi)^.5
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

#Convert wetlands.shp to raster
wetlands.grd <- dem*0
extent(wetlands.grd) <- extent(dem)
wetlands.grd<- rasterize(wetlands.shp, wetlands.grd, pnts.shp@coords[over(wetlands.shp, pnts.shp),4])

#Incoorporate wetlands into raster~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Add wetland invert to landscape
delete.grd<-wetlands.grd*0
delete.grd[is.na(delete.grd)] <- 1
invert.grd<-wetlands.grd-invert
invert.grd[is.na(invert.grd)] <- 0
dem<-dem*delete.grd+invert.grd

#Smooth Raster
dem_filter<- focal(dem, w=matrix(1/225,nrow=15,ncol=15))

#Export Raster
writeRaster(dem,"dem_random.asc" )

#plot
plot(dem_filter, col=colorRampPalette(c("#0571b0","#92c5de","#f7f7f7","#f4a582","#ca0020"))(50))