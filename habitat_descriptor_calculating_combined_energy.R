library(raster)


#install.packages('sdm')

#clear away all variables
rm(list= ls())
gc()



#Author: Mickaël Vasquez
#Version: 20210427_01
#Date: 27/04/2021
#Calculates the layer on overall energy by combining the layer on current-induces energy and the one
#on wave-induced energy. The rule is that the maximum value of the two layers is kept 
#e.g. if for a cell the current is moderate and the wave is low, then the combined energy is moderate.  
#When it comes to the confidence, for each cell the confidence is that of the habitat descriptor that contributes to the calculation of the overall energy, 
#ie if for a cell the combined energy is that of wave (because wave > current), then the combined confidence is that of wave. 
#If for a cell the combined energy is that of wave and current (because wave = current), the confidence is the average confidence of
#wave and current, rounded up


#-------------------------- script parametrisation ---------------------------------------------------
#chunk size and max memory to load for processings
chunksize<-1e+09
maxmemory<-1e+09

#working directory
workingDirectory<-"E:/travail/EUSeaMap4/WP1/modeling/script_publication/hands_on_dataset/partialBoB"

#path to the input raster file name corresponding to the current-induced energy
current_fileName<-"output/Currents.tif"
#path to the input confidence raster file name corresponding to the current-induced energy
current_confidence_fileName<-"output/Currents_confidence_overall.tif"

#path to the input raster file name corresponding to the current-induced energy
wave_fileName<-"output/Wave.tif"
#path to the input confidence raster file name corresponding to the current-induced energy
wave_confidence_fileName<-"output/Wave_confidence_overall.tif"

#path  to the output energy file name
ouput_fileName<-"output/energy.tif"

#path  to the output confidence file name
ouput_confidence_fileName<-"output/energy_confidence_overall.tif"

#Is it required to calculate confidence?
calculateConfidence<-TRUE


########################################## FUNCTIONS ##############################################################################
check_Inputs<-function () {
  #checks that all the input files exist
  
  hasError<-FALSE
  msg<-""
  
  if (calculateConfidence) {
    inputChecks<-c(current_fileName,
                   current_confidence_fileName,
                   wave_fileName,
                   wave_confidence_fileName
                  )
  } else {
    inputChecks<-c(current_fileName,
                   wave_fileName
                  )
  }
  for (i in 1:length(inputChecks)) {
    file<-inputChecks[i]
    if (!file.exists(file)) {
      msg<-sprintf("script aborted: file %s does not exist",
                   file)
      hasError<-TRUE
      break
    }
  }
    
  lstCheckRes<-list(hasError,msg)
}


########################################### MAIN SCRIPT ##############################################################################

start <- Sys.time()

setwd(workingDirectory)

#temp folder will be used by R if requires to create temp files
if (!file.exists("temp")){
  dir.create("temp")
  
}
#set raster options (most important are chunksize and maxmemory)
rasterOptions(tmpdir=paste(getwd(),"/temp",sep=""),progress="text",format="GTiff",
              chunksize=chunksize,
              maxmemory=maxmemory)

suppressWarnings(checkList<-check_Inputs())
hasError<-checkList[[1]]


if (! hasError) {


  print("Combining wave and current rasters")
  
  
  current_raster<-raster(current_fileName)
  
  wave_raster<-raster(wave_fileName)
  
  stack_energy<-stack (current_raster,wave_raster)
  
  writeRaster(max(stack_energy),ouput_fileName,
              datatype='INT1U',overwrite=TRUE)
  
  
  rm(stack_energy)
  
  if (calculateConfidence) {
    
    print("Combining wave and current confidence rasters")
    
    current_confidence_raster<-raster(current_confidence_fileName)
    
    wave_confidence_raster<-raster(wave_confidence_fileName)
    
    raster_combined_confidence<-overlay(current_raster,wave_raster,
                                        current_confidence_raster,wave_confidence_raster,
                                        fun=function(cur,wave,
                                                     cur_conf,wave_conf){ifelse(cur>wave,cur_conf,
                                                                                ifelse(wave>cur,wave_conf,
                                                                                       floor(((cur_conf+wave_conf)/2)+0.5)))})
    writeRaster(raster_combined_confidence,ouput_confidence_fileName,
                datatype='INT1U',overwrite=TRUE)
  }
  
  
  #remove all tmp file created by the raster package while overlaying ou writing chunk by chunk
  file.remove(Sys.glob(file.path("temp", "r_tmp*")))
  end <- Sys.time()
  print (paste("Completed in",difftime(end,start),"!"))
} else print (checkList[[2]])

gc()





