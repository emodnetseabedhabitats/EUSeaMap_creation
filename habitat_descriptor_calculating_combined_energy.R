library(terra)

#clear away all variables
rm(list= ls())
gc()



#Author: Mickael Vasquez
#Version: 20230503_01
#Date: 03/05/2023
#N.B.: main update compared to previous version: the use of terra package instead of raster package
#Calculates the layer on overall energy by combining the layer on current-induces energy and the one
#on wave-induced energy. The rule is that the maximum value of the two layers is kept 
#e.g. if for a cell the current is moderate and the wave is low, then the combined energy is moderate.  
#When it comes to the confidence, for each cell the confidence is that of the habitat descriptor that contributes to the calculation of the overall energy, 
#ie if for a cell the combined energy is that of wave (because wave > current), then the combined confidence is that of wave. 
#If for a cell the combined energy is that of wave and current (because wave = current), the confidence is the average confidence of
#wave and current, rounded up
#Additionally to the full-coverage energy confidence raster, a raster dedicated to EUNIS 2007-11 is created; In that raster
#only cells where seabed substrate is rock and biozone is infralittoral or circalittoral will have confidence values


#-------------------------- script parametrisation ---------------------------------------------------
#max memory to use (in GBytes)
maxmemory<-20

#working directory
workingDirectory<-"E:/travail/euseamap_phase5/WP1/modeling/models/Atlantic_SE"

#path to the input raster file name corresponding to the current-induced energy
current_fileName<-"output/currents.tif"
#path to the input confidence raster file name corresponding to the current-induced energy
current_confidence_fileName<-"output/currents_confidence_overall.tif"

#path to the input raster file name corresponding to the current-induced energy
wave_fileName<-"output/waves.tif"
#path to the input confidence raster file name corresponding to the current-induced energy
wave_confidence_fileName<-"output/waves_confidence_overall.tif"

#path to the input raster file name corresponding to the biozones. Will be used for the creation of the EUNIS2007-11 confidence raster 
biozone_raster_fileName<-"output/biozones.tif"

#path to the input raster file name corresponding to the seabed substrate. Will be used for the creation of the EUNIS2007-11 confidence raster 
substrate_raster_fileName<-"input/seabed_substrate_2021_raster.tif"

#path  to the output energy file name
ouput_fileName<-"output/energy.tif"

#path  to the output confidence file name
ouput_confidence_fileName<-"output/energy_confidence_overall.tif"

#path  to the output EUNIS2007-11 confidence file name
ouput_EUNIS_confidence_fileName<-"output/energy_confidence_overall_EUNIS2007-11.tif"

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
                   wave_confidence_fileName,
                   biozone_raster_fileName,
                   substrate_raster_fileName
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

#set terra options (most important is memmax)
terraOptions(tempdir=paste(getwd(),"/temp",sep=""),progress=10,
             memmax=maxmemory)

suppressWarnings(checkList<-check_Inputs())
hasError<-checkList[[1]]


if (! hasError) {


  print("Combining wave and current rasters")
  
  
  current_raster<-rast(current_fileName)
  
  wave_raster<-rast(wave_fileName)
  
  stack_energy<-rast(list(current_raster,wave_raster))
  
  writeRaster(max(stack_energy),ouput_fileName,
              datatype='INT1U',overwrite=TRUE)
  
  
  rm(stack_energy)
  #remove all tmp file created by the terra package
  tmpFiles(remove=T)
  
  if (calculateConfidence) {
    
    print("Combining wave and current confidence rasters")
    
    current_confidence_raster<-rast(current_confidence_fileName)
    
    wave_confidence_raster<-rast(wave_confidence_fileName)
    
    raster_combined_confidence<-lapp(sds(current_raster,wave_raster,
                                         current_confidence_raster,wave_confidence_raster),
                                         fun=function(cur,wave,
                                                      cur_conf,wave_conf){ifelse(cur>wave,cur_conf,
                                                                                ifelse(wave>cur,wave_conf,
                                                                                       floor(((cur_conf+wave_conf)/2)+0.5)))})
    writeRaster(raster_combined_confidence,ouput_confidence_fileName,
                datatype='INT1U',overwrite=TRUE)
    
    #Only where substrate=rock and biozone in (infra,circa,deep circa)
    print("Creating EUNIS2007-11 confidence raster")
    biozone_raster<-rast(biozone_raster_fileName)
    substrate_raster<-rast(substrate_raster_fileName)
    raster_combined_confidence<-lapp(sds(raster_combined_confidence,substrate_raster,biozone_raster),
                                     fun=function(e,s,b){ifelse(s==70 & b <= 30, e, NA)})
                                       
    writeRaster(raster_combined_confidence,ouput_EUNIS_confidence_fileName,
                datatype='INT1U',overwrite=TRUE)
    
    rm(biozone_raster,substrate_raster,raster_combined_confidence)
    #remove all tmp file created by the terra package
    tmpFiles(remove=T)
  }
  
  
  
  end <- Sys.time()
  print (paste("Completed in",difftime(end,start),"!"))
} else print (checkList[[2]])

gc()





