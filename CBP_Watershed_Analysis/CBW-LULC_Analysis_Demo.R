####################################################################################
#Name: CBW-LULC Analysis Demo
#Coder: C. Nathan Jones
#Date: 6/6/2018
#Purpose: Demo watershed delineation and LULC change analysis for USGS gages accross
#         the Chesepeake Bay Watershed
#         # Helpful link -- http://ibis.geog.ubc.ca/~rdmoore/rcode/ShannonFallsMap.r
####################################################################################

#To do: 
#Create "scripting" setup to initiate multiple instances of R
#For each instance of R, copy paste SAGA cmd executable so you don't execute at the same time 
#         as another script.

####################################################################################
#Step 1: Setup Workspace------------------------------------------------------------
####################################################################################
#Clear Memory
rm(list=ls(all=TRUE))

#Set working directory
setwd("//storage.research.sesync.org/njones-data/Research Projects/LULC_CBW/")

#Add appropriate libararies
#install.packages(c('RSAGA', 'rgdal', 'rgeos', 'raster', 'sp', 'dplyr', 'magrittr'))
library(RSAGA)
library(rgdal)
library(rgeos)
library(raster)
library(sp)
library(dplyr)
library(magrittr)

#Download required data
dem<-raster("Spatial_Data/NHDPlus02/Elev_Unit_a/elev_cm")
gages<-readOGR("Spatial_Data/NHDPlus02/StreamGageEvent.shp")
  gages<-gages[gages$SUBREGION %in% c("0205","0206","0207","0208"),]
  gages<-spTransform(gages, dem@crs)
HUC<-readOGR("Spatial_Data/NHDPlus02/WBD/WBD_Subwatershed.shp")
  HUC<-spTransform(HUC, dem@crs)
LULC2001<-raster("Spatial_Data/nlcd_2001_landcover_2011_edition_2014_10_10/nlcd_2001_landcover_2011_edition_2014_10_10/nlcd_2001_landcover_2011_edition_2014_10_10.img")
LULC2006<-raster("Spatial_Data/nlcd_2006_landcover_2011_edition_2014_10_10/nlcd_2006_landcover_2011_edition_2014_10_10.img")
LULC2011<-raster("Spatial_Data/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img")

####################################################################################
#Step 2: Delineation Function-------------------------------------------------------
####################################################################################
#Create function to delineate watershed and complete LULC change anlayis~~~~~~~~~~~~
fun<-function(ID){

    #Setup work space~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Identify gage
    gage<-gages[ID,]
    staid<-gage$SOURCE_FEA
    
    #Create working directory
    dir.create(paste0("//storage.research.sesync.org/njones-data/Workspace/",staid,"/"))
    setwd(paste0("//storage.research.sesync.org/njones-data/Workspace/",staid,"/"))
    
    #Setup SAGA environment
    myenv = rsaga.env(workspace = getwd(), 
                      path = 'C:\\Program Files/QGIS 2.16.0/apps/saga', 
                      modules = 'C:\\Program Files/QGIS 2.16.0/apps/saga/modules'
    )
    
    #Define HUC of interest
    HUC<-HUC[gage,]
    
    #crop dem to huc extent
    res<-res(dem)[1]
    ext<-raster(extent(HUC)*1.1, res=res)
    dem<-crop(dem, ext)
    
    #Delineate watershed using RSAGA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Write outputs
    writeRaster(dem, file="dem.sgrd", format="SAGA",overwrite=T)
    
    #Fill sinks
    rsaga.fill.sinks(
      in.dem = "dem.sgrd", out.dem = "dem_fill.sgrd", 
      method = "planchon.darboux.2001", minslope = 0.1,
      env = myenv
    )
    
    #Flow Accumulation
    rsaga.parallel.processing("dem_fill.sgrd", 
                              out.carea = "fac.sgrd", 
                              env = myenv)
    
    #Snap pour point
    fac<-raster("fac.sdat")
    gage<-data.frame(extract(fac, gBuffer(gage, width=90), cellnumbers=T))
    gage<-gage$cell[which.max(gage$value)]
    gage<-data.frame(xyFromCell(dem,gage))
    coordinates(gage) = ~ x + y
    
    #Watershed Delineation
    rsaga.geoprocessor(
      lib = 'ta_hydrology', module = 4, env = myenv,
      param = list(TARGET_PT_X = gage@coords[1],
                   TARGET_PT_Y = gage@coords[2],
                   ELEVATION = 'dem_fill.sgrd',
                   AREA = 'basin.sgrd',
                   METHOD = 0
      )					 
    )
    
    #Convert to polygon
    rsaga.geoprocessor(
      lib = 'shapes_grid', module = 6, env = myenv,
      param = list(GRID = 'basin.sgrd',
                   POLYGONS = 'basin.shp',
                   CLASS_ALL = 0,
                   CLASS_ID = 100,
                   SPLIT = 0
      )
    )
    
    #Pull back into workspace
    basin.shp<-readOGR("basin.shp")
    
    #Tabulate LULC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #crop lulc2001 to huc extent
    res<-res(LULC2001)[1]
    ext<-raster(extent(basin.shp)*1.1, res=res)
    LULC2001<-crop(LULC2001, ext)
    LULC2001<-mask(LULC2001, basin.shp)
    LULC2001<-extract(LULC2001)
    LULC2001<-LULC2001[!(is.na(LULC2001))]
    
    #crop lulc2006 to huc extent
    res<-res(LULC2006)[1]
    ext<-raster(extent(basin.shp)*1.1, res=res)
    LULC2006<-crop(LULC2006, ext)
    LULC2006<-mask(LULC2006, basin.shp)
    LULC2006<-extract(LULC2006)
    LULC2006<-LULC2006[!(is.na(LULC2006))]
    
    #crop lulc2011 to huc extent
    res<-res(LULC2011)[1]
    ext<-raster(extent(basin.shp)*1.1, res=res)
    LULC2011<-crop(LULC2011, ext)
    LULC2011<-mask(LULC2011, basin.shp)
    LULC2011<-extract(LULC2011)
    LULC2011<-LULC2011[!(is.na(LULC2011))]
    
    #Create output data.frame
    output<-data.frame(gage=staid, 
                       area_ha = gArea(basin.shp)/10000)
    
    #Caluclate water area
    output$urban2001<-length(LULC2001[LULC2001>=20 & LULC2001<30])*900/10000
    output$urban2006<-length(LULC2006[LULC2006>=20 & LULC2006<30])*900/10000
    output$urban2011<-length(LULC2011[LULC2011>=20 & LULC2011<30])*900/10000
    
    #Calculate barren area
    output$barren2001<-length(LULC2001[LULC2001>=30 & LULC2001<40])*900/10000
    output$barren2006<-length(LULC2006[LULC2006>=30 & LULC2006<40])*900/10000
    output$barren2011<-length(LULC2011[LULC2011>=30 & LULC2011<40])*900/10000
    
    #Calculate forest area
    output$forest2001<-length(LULC2001[LULC2001>=40 & LULC2001<50])*900/10000
    output$forest2006<-length(LULC2006[LULC2006>=40 & LULC2006<50])*900/10000
    output$forest2011<-length(LULC2011[LULC2011>=40 & LULC2011<50])*900/10000
    
    #Calculate shrub area
    output$shrub2001<-length(LULC2001[LULC2001>=50 & LULC2001<60])*900/10000
    output$shrub2006<-length(LULC2006[LULC2006>=50 & LULC2006<60])*900/10000
    output$shrub2011<-length(LULC2011[LULC2011>=50 & LULC2011<60])*900/10000
    
    #Calculate Herbaceous area
    output$herb2001<-length(LULC2001[LULC2001>=70 & LULC2001<80])*900/10000
    output$herb2006<-length(LULC2006[LULC2006>=70 & LULC2006<80])*900/10000
    output$herb2011<-length(LULC2011[LULC2011>=70 & LULC2011<80])*900/10000
    
    #Calculate Crop area
    output$crop2001<-length(LULC2001[LULC2001>=80 & LULC2001<90])*900/10000
    output$crop2006<-length(LULC2006[LULC2006>=80 & LULC2006<90])*900/10000
    output$crop2011<-length(LULC2011[LULC2011>=80 & LULC2011<90])*900/10000
    
    #Calculate Wetlands area
    output$wetland2001<-length(LULC2001[LULC2001>=90])*900/10000
    output$wetland2006<-length(LULC2006[LULC2006>=90])*900/10000
    output$wetland2011<-length(LULC2011[LULC2011>=90])*900/10000
    
    #Clean up workspace~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Kill directory
    unlink(paste0("//storage.research.sesync.org/njones-data/Workspace/",staid), recursive = T)

    #export output data.frame
    output
}

#Create function to catch errors
execute<-function(n){tryCatch(fun(n), error=function(e) c(gages$SOURCE_FEA[n], rep(0,22)))}

#Run function!
t0<-Sys.time()
output<-lapply(seq(1,778), execute)
tf<-Sys.time()

#unlist
setwd("//storage.research.sesync.org/njones-data/Research Projects/LULC_CBW/")
output<-do.call(rbind, output)
write.csv(output, "output.csv")
