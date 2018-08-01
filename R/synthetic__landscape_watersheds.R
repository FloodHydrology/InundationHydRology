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
setwd("M:\\McLaughlin_Lab/Jones/workspace/PC Modeling Workshop/Model_Comparison/DEM/")

#Download packages (use install.packages function if not already downloaded)
library(raster)
library(sp)
library(gstat)
library(rgdal)
library(maptools)
library(foreign)
library(rgeos)
library(actuar)
library(poweRlaw)

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
#Determine number, size, and location of wetlands~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#determine number and size of wet
wetland_area<-area*0.12 #From Van Meter and Basu 2015 (Des Moines Lobe)

#Iterate through wetland scenarios until correct area
set.seed(111)
dif<-area
while(abs(dif)>(10)){
  n.wetlands<-round(runif(1, 5, 500), digits = 0)
  pnts<-rplcon(n.wetlands, 10^3, 1.67)
  pnts<-pnts[pnts<(0.1*wetland_area)]
  dif<-wetland_area-sum(pnts)
}
n.wetlands<-length(pnts)

#Randomly select location of each wetland
pnts<-data.frame(coordinates(sampleRandom(dem,n.wetlands, sp=T)), pnts)
  colnames(pnts)<-c("x","y", "area_m2")

#add wetid data
pnts<-pnts[order(-pnts$area_m2),]
pnts$WetID<-seq(1,length(pnts[,1]))

#calculate wetland volume and depth
p<-4 #shape factor from Hayashi and Kamp [200]
pnts$volume_m3<-(0.25*(pnts$area_m2/10000)^1.4742)*10000 #Equation  from Wu and Lane [2016]
pnts$max_depth_m<-(pnts$volume_m3*(1+2/p))/pnts$area_m2
  
#add z data and create shp
pnts$z<-extract(dem, pnts[,c("x","y")])

#add data about the watershed depth
pnts$ws_depth<-abs(rnorm(nrow(pnts),2,0.5))

#Create Function to add wetland shape (with bathymetry)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fun<-function(WetID){
  #Setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #define variables
  area_max<-pnts$area_m2[pnts$WetID==WetID]
  x<-pnts$x[pnts$WetID==WetID]
  y<-pnts$y[pnts$WetID==WetID]
  z<-pnts$z[pnts$WetID==WetID]
  ws_depth<-pnts$ws_depth[pnts$WetID==WetID]
  depth<-pnts$max_depth_m[pnts$WetID==WetID]
  
  #create circle function to define circle
  circle.fun<-function(area, x, y){
    radius<-(area/pi)^.5
    circle<-seq(0, 2 * pi, length.out = 2*pi*sqrt(area/pi))
    circle<-cbind(x + radius * sin(circle), y + radius * cos(circle))
    circle
  }
  
  #Represent bathymetry with points~~~~~~~~~~~~~~~~~~~~~~~~
  #Calculate variables
  radius<-(area_max/pi)^0.5
  n<-round(radius/cell_size)
  
  #create dataframe to house points
  bath<-data.frame(matrix(0, ncol=3))
    colnames(bath)<-c("x","y","z")
  bath$x[1]<-x
  bath$y[1]<-y
  bath$z[1]<-z-depth
    
  #Add "rings" to bath df
  for(i in 1:(n+1)){
    area_ring<-pi*(i*cell_size)^2
    df<-data.frame(circle.fun(area_ring, x,y),depth*(i*cell_size/radius)^p+(z-depth))
      colnames(df)<-c("x","y","z")
    bath<-rbind(bath,df)
  }
  
  #Adjust elevation for watershed invert elevation
  bath$z<-bath$z-ws_depth
  
  #Represent watershed surface~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Calculate variables
  ws_radius<-(area_max*10/pi)^0.5
  m<-round((ws_radius-radius)/cell_size)
  
  #create dataframe to house points
  ws<-data.frame(matrix(0, ncol=3))
  colnames(ws)<-c("x","y","z")

  #Add "rings" to ws df
  for(i in 0:m){
    area_ring<-pi*(i*cell_size+radius)^2
    df<-data.frame(circle.fun(area_ring, x,y),ws_depth/(ws_radius-radius)*(i*cell_size)+(z-ws_depth))
    colnames(df)<-c("x","y","z")
    ws<-rbind(ws,df)
  }
  ws<-ws[-1,]
  
  #Create raster of Wetland Bathymetry~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Conbime bathymetry and ws surface
  bath<-rbind(bath,ws)
  
  #Remove points that are > current dem surface
  bath$dem<-extract(dem,bath[,1:2])
  bath<-bath[bath$z<bath$dem,]
  
  #If any points are left:
  if(nrow(bath)>1){
  
    #Create gird
    wetland.grd<-rasterize(bath[,1:2], dem, bath[,3])
    
    #Create Clip
    clip<-wetland.grd*0
    clip[is.na(clip)] <- 1
    wetland.grd[is.na(wetland.grd)] <- 0
    
    #append dem
    dem<-dem*clip+wetland.grd
    
    #assign dem to the global environment
    assign('dem', dem, envir = .GlobalEnv)
  }
}

#Run function
for(i in 1:length(pnts[,1])){fun(i)}

#Smooth with raster filter
dem<-focal(dem, w=matrix(1/25,nrow=5,ncol=5)) 

#Export Raster
writeRaster(dem,"dem_ws_bathymetry.asc",overwrite=TRUE)

#plot
nbcol <- 99
color <- c("grey",terrain.colors(nbcol))
plot(dem, col=color)

############################################################################################
#Plot in 3D (https://chitchatr.wordpress.com/2014/04/25/plotting-dems-in-3d-using-persp-part-2/)
#create zDatamatrix
zData<-as.matrix(dem)
x = (cell_size * (1:nrow(zData)))    
y = (cell_size * (1:ncol(zData)))
nrzmat <- nrow(zData)
nczmat <- ncol(zData)
facetValues <- (zData[-1, -1] + zData[-1, -nczmat] + zData[-nrzmat, -1] + zData[-nrzmat, -nczmat])/4
facetcol <- cut(facetValues, nbcol+1)

jpeg("C:\\Users/cnjones/Desktop/test.jpg", height=6, width=6, units="in", res=300)
res = persp(x, y, z = zData*50, theta = 120, phi = 45,
            col = color[facetcol],
            scale = FALSE, expand = 0.75, 
            ltheta = 75, shade = 0.75, border = NA,
            box = F, ticktype = "detailed")
dev.off()