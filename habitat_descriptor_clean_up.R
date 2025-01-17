library(terra)
library(sf)


## clean up memory
rm(list= ls())
gc()

#Author: Mickael Vasquez
#Version: 20250115_01

# Cleans up a habitat descriptor raster
#A shapefile is provided as input, containing a series of polygon within which the raster has to be cleaned up.
#3 types of action are possible:
#  1) within the polygon the cell values are removed (i.e. replaced by NAs)
#  2) within the polygon the cell values of a given value are replaced with another value
#  3) within the polygon a sieve filter is applied (to reduce the ‘salt-and-pepper’ effect).
#The attributes of the shapefile are the following:
#  - TYPE: describe the action to be done (-1: remove values inside the polygon, 1: replace a pixel value with another value; 2: apply a sieve filter)
#  - IN_VAL: for TYPE=1 only. Value to be replaced
#  - OUT_VAL: for TYPE=1 only. new value
#  - MAX_SIZE: for TYPE=2 only. Corresponds to the 'threshold' parameter of the function sieve




#-------------------------- script parametrisation ---------------------------------------------------
#max memory to use (in GBytes)
maxmemory<-20
#working directory
workingDirectory<-"D:/temp/EUSeaMap_exercises/cleanup_example"
clean_up_shape_file<-"input/cleanup_biozones.shp"
input_raster_file<-"output/biozones.tif"
output_raster_file<-"output/biozones_corrected.tif"

#-----------------------------------------------------------------------------------------------------


########################################### FUNCTIONS ##############################################################################

check_clean_up_inputs<-function (sdf_clean_up) {
  
  sdf_clean_up_copy<-sdf_clean_up
  sdf_clean_up_copy$cpt<-seq(0,(nrow(sdf_clean_up_copy)-1))
  hasError<-FALSE
  msg<-""
  sdf_test<-sdf_clean_up_copy[!sdf_clean_up_copy$TYPE %in% c(-1,1,2),]
  if (nrow(sdf_test) != 0) {
    hasError<-TRUE
    msg<-sprintf("record number %s: TYPE=%s. This is not an appropriate value for TYPE (should be -1, 1, or 2)",sdf_test$cpt,sdf_test$TYPE)
  }
  sdf_test<-sdf_clean_up_copy[sdf_clean_up_copy$TYPE == 1 & sdf_clean_up_copy$IN_VAL==0,]
  if (nrow(sdf_test) != 0) {
    hasError<-TRUE
    msg<-sprintf("record number %s: TYPE=1 and IN_VAL=0. 0 is not an appropriate value for IN_VAL when TYPE=1",sdf_test$cpt)
  }
  sdf_test<-sdf_clean_up_copy[sdf_clean_up_copy$TYPE == 1 & sdf_clean_up_copy$OUT_VAL==0,]
  if (nrow(sdf_test) != 0) {
    hasError<-TRUE
    msg<-sprintf("record number %s: TYPE=1 and IN_VAL=0. 0 is not an appropriate value for IN_VAL when TYPE=1",sdf_test$cpt)
  }
  sdf_test<-sdf_clean_up_copy[sdf_clean_up_copy$TYPE == 2 & sdf_clean_up_copy$MAX_SIZE==0,]
  if (nrow(sdf_test) != 0) {
    hasError<-TRUE
    msg<-sprintf("record number %s: TYPE=2 and MAX_SIZE=0. 0 is not an appropriate value for MAX_SIZE when TYPE=2",sdf_test$cpt)
  }
  lstCheckRes<-list(hasError,msg)
}


clean_up_raster<-function (sdf,input_raster) {
  #sdf is a dataframe with one single record
  rst<-crop(input_raster,ext(sdf))
  rst_clean_up<-rasterize(vect(sdf), rst, field = "TYPE")
  rst<-lapp(sds(rst,rst_clean_up),
            fun=function(r,rc){ifelse(!is.na(rc),r,NA)})
  if (sdf$TYPE==1) rst<-ifel(rst==sdf$IN_VAL,sdf$OUT_VAL,NA)
  else if (sdf$TYPE==2) {
    rst<-sieve(rst, sdf$MAX_SIZE)
    rst<-ifel(rst==0,NA,rst)
  }
  rst
}


########################################### MAIN SCRIPT ##############################################################################
start <- Sys.time()

setwd(workingDirectory)




#temp folder will be used by R if requires to create temp files
if (!file.exists("temp")){
  dir.create("temp")
  
}

terraOptions(tempdir=paste(getwd(),"/temp",sep=""),progress=10,
             memmax=20)

sdf_clean_up <- st_read(clean_up_shape_file)

suppressWarnings(checkList<-check_clean_up_inputs(sdf_clean_up))
hasError<-checkList[[1]]


if (! hasError) {

  input_raster<-rast(input_raster_file)
  output_raster<-input_raster
  sdf_clean_up_tmp<-sdf_clean_up[sdf_clean_up$TYPE==-1,]
  if (nrow(sdf_clean_up_tmp)!=0) { 
    print ("erasing pixels inside eraser polygons (TYPE=-1)")
    rst_eraser<-rasterize(vect(sdf_clean_up_tmp), input_raster, field = "TYPE")
    output_raster<-lapp(sds(output_raster,rst_eraser),
                          fun=function(rst1,rst2){ifelse(!is.na(rst2),NA,rst1)})
  }
  sdf_clean_up_tmp<-sdf_clean_up[sdf_clean_up$TYPE == 2,]
  if (nrow(sdf_clean_up_tmp)!=0) { 
    lst_cleaned_up_rasters<-list()
    print ("replacing existing pixel values by other ones inside sieve polygons (TYPE=2)")
    for (i in 1:nrow(sdf_clean_up_tmp)) {
      print (sprintf("      polygon %s out of %s",i,nrow(sdf_clean_up_tmp)))
      sdf_tmp<-sdf_clean_up_tmp[i,]
      rst<-clean_up_raster(sdf_tmp,output_raster)
      lst_cleaned_up_rasters[[i]]<-rst
    }
    
    print ("merging...")
    lst_cleaned_up_rasters[[nrow(sdf_clean_up_tmp)+1]]<-output_raster
    output_raster<-do.call(merge,lst_cleaned_up_rasters)
  }
  sdf_clean_up_tmp<-sdf_clean_up[sdf_clean_up$TYPE == 1,]
  if (nrow(sdf_clean_up_tmp)!=0) { 
    lst_cleaned_up_rasters<-list()
    print ("replacing existing pixel values by other ones inside replace polygons (TYPE=1)")
    for (i in 1:nrow(sdf_clean_up_tmp)) {
      print (sprintf("      polygon %s out of %s",i,nrow(sdf_clean_up_tmp)))
      sdf_tmp<-sdf_clean_up_tmp[i,]
      rst<-clean_up_raster(sdf_tmp,output_raster)
      lst_cleaned_up_rasters[[i]]<-rst
    }
    
    print ("merging...")
    lst_cleaned_up_rasters[[nrow(sdf_clean_up_tmp)+1]]<-output_raster
    rst<-do.call(merge,lst_cleaned_up_rasters)
  }
  print (sprintf("writing %s raster",output_raster_file))
  writeRaster(rst, filename=output_raster_file, overwrite=TRUE,  datatype='INT1U')

} else print (checkList[[2]])

tmpFiles(remove=T)

end <- Sys.time()

gc()

print (paste("Completed in",difftime(end,start),"!"))

