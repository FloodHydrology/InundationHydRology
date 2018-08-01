GIW_subshed_delineation<-function(
  workspace="C:\\ScratchWorkspace\\",
  wbt_path="C:/WBT/whitebox_tools",
  dem,
  depressions,
  wetland){
  
  #Set Working Directory
  setwd(paste(workspace))  
  
  #Export rater to workspace
  writeRaster(dem, paste0(workspace,"dem.tif"), overwrite=T)
  
  #Breach single-cell pits
  system(paste(paste(wbt_path), 
               "-r=BreachSingleCellPits", 
               paste0("--wd=",workspace),
               "--dem='dem.tif'", 
               "-o='dem_breachedsinglecells.tif'"))
  
  #Filter the DEM
  system(paste(paste(wbt_path), 
               "-r=EdgePreservingMeanFilter", 
               paste0("--wd=",workspace),
               "-i='dem_breachedsinglecells.tif'", 
               "-o='dem_edgepreservingfilter.tif'",
               "filter=10", 
               "threshold=100"))
  
  #Breach Analysis of the DEM
  system(paste(paste(wbt_path), 
               "-r=BreachDepressions", 
               paste0("--wd=",workspace),
               "--dem=dem_edgepreservingfilter.tif", 
               "-o=dem_breach.tif"))
  
  #Flow Direction
  system(paste(paste(wbt_path), 
               "-r=D8Pointer", 
               paste0("--wd=",workspace),
               "--dem='dem_breach.tif'", 
               "-o='fdr.tif'",
               "--out_type=sca"))
  
  #Flow Accumulation
  system(paste(paste(wbt_path), 
               "-r=DInfFlowAccumulation", 
               paste0("--wd=",workspace),
               "--dem='dem_breach.tif'", 
               "-o='fac.tif'",
               "--out_type=sca"))
  
  #Read fac and fdr rasters into R environment
  fdr<-raster(paste0(workspace,"fdr.tif"))
  fac<-raster(paste0(workspace,"fac.tif"))
  
  #Extract depression of interest
  wetland_dep<-depressions
  wetland_dep[wetland_dep!=raster::extract(depressions, wetland)]<-NA
  wetland_dep<-wetland_dep*0+1
  
  #Define watershed pour point (max fac in depression)  
  wetland_fac<-wetland_dep*fac
  max_wetland_fac<-cellStats(wetland_fac, max)
  pp<-rasterToPoints(wetland_fac, fun=function(x){x==max_wetland_fac})
  pp<-SpatialPointsDataFrame(pp, data.frame(x=1))
  pp@proj4string<-dem@crs
  writeOGR(pp,paste0(workspace,"."),"pp", drive="ESRI Shapefile", overwrite=T)
  
  #Watershed Delineation
  system(paste(paste(wbt_path), 
               "-r=Watershed", 
               paste0("--wd=",workspace),
               "--d8_pntr='fdr.tif'", 
               "--pour_pts=pp.shp",
               "-o=watershed.tif"))
  
  #Import watershed shape
  watershed<-raster(paste0(workspace,"watershed.tif"))
  
  #Clip relevant layers to watershed
  fac<-fac*watershed
  fdr<-fac*watershed
  depressions<-depressions*watershed
  
  #If there are other basins complete the analysis again to remove
  depressions[depressions==raster::extract(depressions, wetland)]<-NA
  n_depression<-length(unique(depressions))
  while(n_depression>0){
    #Extract depression of interest
    wetland_dep<-depressions
    wetland_dep[wetland_dep!=unique(wetland_dep)[1]]<-NA
    wetland_dep<-wetland_dep*0+1
    
    #Define watershed pour point (max fac in depression)  
    wetland_fac<-wetland_dep*fac
    max_wetland_fac<-cellStats(wetland_fac, max)
    pp<-rasterToPoints(wetland_fac, fun=function(x){x==max_wetland_fac})
    pp<-matrix(pp[1,], ncol=3)
    pp<-SpatialPointsDataFrame(pp, data.frame(x=1))
    pp@proj4string<-dem@crs
    writeOGR(pp,paste0(workspace,"."),"pp", drive="ESRI Shapefile", overwrite=T)
    
    #Watershed Delineation
    system(paste(paste(wbt_path), 
                 "-r=Watershed", 
                 paste0("--wd=",workspace),
                 "--d8_pntr='fdr.tif'", 
                 "--pour_pts=pp.shp",
                 "-o=subshed.tif"))
    subshed<-raster(paste0(workspace,"subshed.tif"))
    
    #Remove subshed from watershed
    subshed[is.na(subshed)]<-0
    watershed<-watershed-subshed
    watershed[watershed!=1]<-NA
    
    #Recalculate nubmer of depressions
    depressions<-depressions*watershed
    n_depression<-length(unique(depressions))
  }
  
  #Export Watershed
  writeRaster(watershed, paste0(paste(wetland$Name),"_subshed.asc"), overwrite=T)
  watershed
}
