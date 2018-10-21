WatershedAnalysis<-function(
  dem=dem,  
  pnts=pnts,
  unique_id="FldActI",
  streams=streams,
  mask=mask, 
  threshold=300,
  scratch_dir=scratch_dir,
  wbt_dir=wbt_dir,
  data_dir=data_dir, 
  output_dir=output_dir){
  
  #Mask dem and gages
  dem<-crop(dem, as.vector(st_bbox(mask))[c(1, 3, 2, 4)])
  pnts<-pnts[mask,]
  
  #Export DEM and stream layer to local working directory
  writeRaster(dem, 
              paste0(scratch_dir,"dem.tif"), 
              overwrite=T)
  st_write(pnts, paste0(scratch_dir,"pnts.shp"), delete_layer = T)
           
  #Gaussian Filter
  system(paste(paste(wbt_dir), 
               "-r=GaussianFilter", 
               paste0("--wd=",scratch_dir),
               "-i='dem.tif'", 
               "-o='dem_filter.tif'",
               "--sigma=3"))
  
  #Fill "single cell" depressions
  system(paste(paste(wbt_dir),
                "-r=FillSingleCellPits",
                paste0("--wd=",scratch_dir),
                "--dem='dem_filter.tif'",
                "-o='dem_breach_minor.tif'"))
  
  #Breach larger depressions
  system(paste(paste(wbt_dir), 
               "-r=BreachDepressions", 
               paste0("--wd=",scratch_dir),
               "--dem='dem_breach_minor.tif'", 
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
               "-i='pnts.shp'", 
               paste0("--field=",unique_id),
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
               paste0("--snap_dist=",threshold)
  ))
  
  #Convert back to point file
  snapgrid<-raster(paste0(scratch_dir,"pp_snap.tif"))
  snappnts<-data.frame(rasterToPoints(snapgrid, fun=function(x){x>0}))
    colnames(snappnts)<-c('x','y',paste(unique_id))
  snappnts$watershedID<-seq(1, length(snappnts[,1]))
  snappnts<-st_as_sf(snappnts, 
              coords=c("x","y"), 
              crs=paste(dem@crs))
  st_write(snappnts, paste0(scratch_dir,"snap.shp"), delete_layer = T)
  snappnts<-data.frame(snappnts)
  
  #Create folder for watershed files
  dir.create(paste0(scratch_dir,"temp"))
  
  #Watershed Analysis
  system(paste(paste(wbt_dir), 
               "-r=UnnestBasins", 
               paste0("--wd=",scratch_dir),
               "--d8_pntr='fdr.tif'",
               "--pour_pts='snap.shp'", 
               paste0("-o='",scratch_dir,"temp/watershed.tif'")))

  #Create list of raster files
  watersheds <- list.files(paste0(scratch_dir,"temp")) 
  watersheds <- paste0(scratch_dir, "temp\\", watersheds)
  
  #Create function to print individual watershed files 
  fun<-function(i){
    r<-raster(watersheds[i])
    l<-unique(r)
    for(j in 1:length(l)){
      #Define unique id
      uid<-snappnts[,c(unique_id,"watershedID")]
      uid<-uid[,unique_id][uid$watershedID==l[j]]
      
      #Isolate watershed points
      n<-r
      n[n!=l[j]]<-NA
      
      #print to output location
      writeRaster(n, paste0(output_dir,uid,".tif"))
    }
  }
  
  #Apply function
  sapply(seq(1, length(watersheds)), fun)
  
  #Remove files
  unlink(paste0(scratch_dir,"temp"), recursive=T)
  file.remove(paste0(scratch_dir,list.files(scratch_dir)))
}