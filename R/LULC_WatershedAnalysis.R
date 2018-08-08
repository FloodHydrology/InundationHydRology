LULC_WatershedAnalysis<-function(
  dem=dem,  
  pnts=pnts,
  unique_id="SOURCE_FEA",
  streams=streams,
  mask=mask, 
  threshold=1000,
  LULC2001=LULC2001,
  LULC2006=LULC2006,
  LULC2011=LULC2011, 
  data_dir=data_dir){
  
  #Mask dem and gages
  dem<-crop(dem, mask)
  pnts<-pnts[mask,]
  
  #Add UID to pnts
  pnts$ID<-seq(1,length(pnts))
  
  #Export DEM and stream layer to local working directory
  writeRaster(dem, paste0(scratch_dir,"dem.tif"), overwrite=T)
  writeOGR(pnts,paste0(scratch_dir,"."),"gages", driver = "ESRI Shapefile", overwrite_layer = T)
  
  #Fill "single cell" depressions
  system(paste(paste(wbt_dir),
               "-r=FillSingleCellPits",
               paste0("--wd=",scratch_dir),
               "--dem='dem.tif'",
               "-o='dem_breach_minor.tif'"))
  
  #Gaussian Filter
  system(paste(paste(wbt_dir), 
               "-r=GaussianFilter", 
               paste0("--wd=",scratch_dir),
               "-i='dem_breach_minor.tif'", 
               "-o='dem_filter.tif'",
               "--sigma=3"))
  
  #Breach larger depressions
  system(paste(paste(wbt_dir), 
               "-r=BreachDepressions", 
               paste0("--wd=",scratch_dir),
               "--dem='dem_filter.tif'", 
               "-o='dem_breach_major.tif'"))
  
  #Create Flow Accumulation Raster
  system(paste(paste(wbt_dir), 
               "-r=D8FlowAccumulation", 
               "--out_type='cells'",
               paste0("--wd=",scratch_dir),
               "--dem='dem_breach_major.tif'", 
               "-o='fac.tif'"))
  
  #Create Stream Raster [fac>1000]
  fac<-raster(paste0(scratch_dir,"fac.tif"))
  fac[fac<threshold]<-NA
  fac<-fac*0+1
  fac@crs<-p
  writeRaster(fac,paste0(scratch_dir,"flowgrid.tiff"), overwrite=T)
  
  #Run flow direction [note we can skip breaching and/or filling sinks b/c we are using NHD data
  system(paste(paste(wbt_dir), 
               "-r=D8Pointer", 
               paste0("--wd=",scratch_dir),
               "--dem='dem_breach_major.tif'", 
               "-o='fdr.tif'",
               "--out_type=sca"))
  
  #Create pour pnt raster
  system(paste(paste(wbt_dir), 
               "-r=VectorPointsToRaster", 
               paste0("--wd=",scratch_dir),
               "-i='gages.shp'", 
               "--field=ID",
               "-o=pp.tif",
               "--assign=min",
               "--nodata",
               "--base=dem.tif"))
  
  #Jenson Snap Pour point
  system(paste(paste(wbt_dir),
               "-r=JensonSnapPourPoints", 
               paste0("--wd=",scratch_dir),
               "--pour_pts='pp.tif'", 
               "--streams='flowgrid.tif'",
               "-o='pp_snap.tif",
               "--snap_dist=1000"
  ))
  
  #Convert back to point file
  snapgrid<-raster(paste0(scratch_dir,"pp_snap.tif"))
  snappnts<-rasterToPoints(snapgrid, fun=function(x){x>0})
  snappnts<-SpatialPointsDataFrame(snappnts[,1:2], data.frame(ID=snappnts))
  
  #Create function to delineate watershed and tabulate LULC
  fun<-function(ID){
    
    #Select Pour Point and Export
    pnt<-snappnts[snappnts$ID.pp_snap==ID,]
    pnt@proj4string<-dem@crs
    writeOGR(pnt,paste0(scratch_dir,"."),"snap", driver = "ESRI Shapefile", overwrite_layer=T)
    
    #Watershed Tool
    system(paste(paste(wbt_dir), 
                 "-r=Watershed", 
                 paste0("--wd=",scratch_dir),
                 "--d8_pntr='fdr.tif'",
                 "--pour_pts='snap.shp'", 
                 paste0("-o='watershed.tif'")))
    
    #Import watershed raster into R environment
    ws_grid<-raster(paste0(scratch_dir,"watershed.tif"), overwrite=T)
    
    #Write copy of raster to data directory
    if(!dir.exists(paste0(data_dir,"watershed"))){
      dir.create(paste0(data_dir,"watershed"))
    }
    writeRaster(ws_grid, 
                paste0(data_dir,
                       "watershed/watershed_",
                       pnts@data[pnts$ID==ID,paste(unique_id)], 
                       ".tif"), 
                overwrite=T)
    
    #Lmit LULC grids to watershed
    LULC2001<-LULC2001*ws_grid
    LULC2006<-LULC2006*ws_grid
    LULC2011<-LULC2011*ws_grid
    
    #Calculate Watershed Area
    cell_size<-res(ws_grid)[1]*res(ws_grid)[2]
    area<-cellStats(ws_grid, sum)*cell_size
    
    #Calculate LULC area
    #Create output data.frame
    output<-data.frame(ID=ID,ws_area = area)
    
    #Caluclate water area
    output$urban2001<-length(LULC2001[LULC2001>=20 & LULC2001<30])*cell_size
    output$urban2006<-length(LULC2006[LULC2006>=20 & LULC2006<30])*cell_size
    output$urban2011<-length(LULC2011[LULC2011>=20 & LULC2011<30])*cell_size
    
    #Calculate barren area
    output$barren2001<-length(LULC2001[LULC2001>=30 & LULC2001<40])*cell_size
    output$barren2006<-length(LULC2006[LULC2006>=30 & LULC2006<40])*cell_size
    output$barren2011<-length(LULC2011[LULC2011>=30 & LULC2011<40])*cell_size
    
    #Calculate forest area
    output$forest2001<-length(LULC2001[LULC2001>=40 & LULC2001<50])*cell_size
    output$forest2006<-length(LULC2006[LULC2006>=40 & LULC2006<50])*cell_size
    output$forest2011<-length(LULC2011[LULC2011>=40 & LULC2011<50])*cell_size
    
    #Calculate shrub area
    output$shrub2001<-length(LULC2001[LULC2001>=50 & LULC2001<60])*cell_size
    output$shrub2006<-length(LULC2006[LULC2006>=50 & LULC2006<60])*cell_size
    output$shrub2011<-length(LULC2011[LULC2011>=50 & LULC2011<60])*cell_size
    
    #Calculate Herbaceous area
    output$herb2001<-length(LULC2001[LULC2001>=70 & LULC2001<80])*cell_size
    output$herb2006<-length(LULC2006[LULC2006>=70 & LULC2006<80])*cell_size
    output$herb2011<-length(LULC2011[LULC2011>=70 & LULC2011<80])*cell_size
    
    #Calculate Crop area
    output$crop2001<-length(LULC2001[LULC2001>=80 & LULC2001<90])*cell_size
    output$crop2006<-length(LULC2006[LULC2006>=80 & LULC2006<90])*cell_size
    output$crop2011<-length(LULC2011[LULC2011>=80 & LULC2011<90])*cell_size
    
    #Calculate Wetlands area
    output$wetland2001<-length(LULC2001[LULC2001>=90])*cell_size
    output$wetland2006<-length(LULC2006[LULC2006>=90])*cell_size
    output$wetland2011<-length(LULC2011[LULC2011>=90])*cell_size
    
    #Export output
    output
  }
  
  #Run function
  output<-lapply(X = snappnts$ID.pp_snap, FUN = fun)
  output<-do.call(rbind, output)
  
  #merge output
  pnts<-pnts@data
  pnts<-merge(pnts, output, by.x="ID", by.y="ID")
  
  #Export output
  pnts
}