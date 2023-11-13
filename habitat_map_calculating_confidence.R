library(terra)


#clear away all variables
rm(list= ls())
gc()



#Author: Mickael Vasquez
#Version: 20230405_01
#Date: 05/04/2023
#N.B.: main update compared to previous version: the use of terra package instead of raster package
#For one or several habitat classification, creation of the habitat confidence layer, which is the the combination of the individual habitat descriptor confidence layers
#For each cell, the minimum confidence of the input confidence layers is kept
#For the individual habitat classifications, the confidence layers of individual habitat descriptor that are involved in the classification
#are provided as input to the script via a csv file, where each rows describes, for each habitat descriptor:
#   - the name of the habitat classification for which the habitat descriptor is used
#   - The name of the confidence raster file
#   - The path to the folder that contains the file, relatively to the working directory


#-------------------------- script parametrisation ---------------------------------------------------
#max memory to use (in GBytes)
maxmemory<-20

#working directory
workingDirectory<-"E:/travail/euseamap_phase5/WP1/modeling/models/Atlantic_SE"

#the config csv file that describes all the bits that have to be merged
config_csvFileName<-"Atlantic_SE_habitat_calculating_habitat_confidence.csv"


#-----------------------------------------------------------------------------------------------------


########################################### FUNCTIONS ##############################################################################
check_Inputs<-function () {
  #checks that
  #i) the configuration file exists
  #ii) all the raster files that are indicated in the configuration file exist
  
  hasError<-FALSE
  msg<-""
  #Does the configuration table exist?
  if (!file.exists(file.path("config_files",config_csvFileName))) {
    msg<-sprintf("script aborted: No file %s was found in config_files folder",
                 config_csvFileName)
    hasError<-TRUE
  } else {
    #the configuration table exists. Open it and check that the raster files that are indicated exist
    units_to_be_merged <- read.table(file=file.path("config_files",config_csvFileName), header=TRUE, sep=",", as.is = TRUE)
    #columns to check because they may contain raster file names. Here there is only one
    columnsToCheck<-c("confidence_fileName")
    for (i in 1:nrow(units_to_be_merged)) {
      folder<-units_to_be_merged$folder[i]
      for (j in 1:length(columnsToCheck)) {
        file<-units_to_be_merged[[columnsToCheck[j]]][i]
        if (!file.exists(file.path(folder,file))) {
          msg<-sprintf("script aborted: file %s, indicated in the column %s of the csv file, does not exist in folder %s",
                       file,
                       columnsToCheck[j],
                       folder)
          hasError<-TRUE
          break
        }
      }
      if (hasError) break
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
  #as.is=TRUE because otherwise string columns are considered as factors, and this make the script fail
  units_to_be_merged <- read.table(file=file.path("config_files",config_csvFileName), header=TRUE, sep=",", as.is = TRUE)
  
  habitatMapNames<-unique(units_to_be_merged$habitat_classification_name)
  
  for (i in 1:length(habitatMapNames)) {
    print (sprintf ("creating the confidence layer for %s...",habitatMapNames[i]))
    habitat_map_units_to_be_merged<-units_to_be_merged[units_to_be_merged$habitat_classification_name==habitatMapNames[i],]
    filesToCombine<-file.path (habitat_map_units_to_be_merged$folder,habitat_map_units_to_be_merged$confidence_fileName)
    
    if (length(filesToCombine) > 1) {
      rastersToCombine<-rast(filesToCombine)
      final_confidence_raster<-min(rastersToCombine,na.rm=TRUE)
    } else {
      #in case the habitat map comprises one habitat descriptor only (not the case currently, but one day maybe...)
      final_confidence_raster<-rast(filesToCombine[1])
    }
    writeRaster(final_confidence_raster,file.path("output",sprintf("%s_confidence_overall.tif",habitatMapNames[i])),
                datatype='INT1U',overwrite=TRUE)
    
    rm(final_confidence_raster,rastersToCombine)
    #remove all tmp file created by the terra package
    tmpFiles(remove=T)
  
  }
  
  end <- Sys.time()
  print (paste("Completed in",difftime(end,start),"!"))
} else print (checkList[[2]])



gc()





