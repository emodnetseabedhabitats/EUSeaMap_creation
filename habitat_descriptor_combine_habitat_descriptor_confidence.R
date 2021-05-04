library(raster)


#install.packages('sdm')

#clear away all variables
rm(list= ls())
gc()



#Author: Mickaël Vasquez
#Version: 20210423_01
#Date: 23/04/2021
#Creation of habitat descriptor layers according to a set of inputs (all described in a standardised csv configuration file)
#For each habitat descriptor class , the csv file provides the following
#   - The name of the habitat descriptor to which the class is attached
#   - The name of the input raster file that characterises the habitat descriptor class spatial distribution
#   - The name of the input raster file that characterises the habitat descriptor class presence probability
#   - The name of the input raster file that characterises the confidence in the habitat descriptor class 
#     occurrence based on the probability 
#   - The name of the input raster file that characterises the overall confidence in the habitat descriptor class occurrence
#   - The name of the folder that contains the input raster file. Classically the output folder
#As outputs the script creates for each habitat descriptor:
#   - A spatial distribution raster, result of all the class spatial distribution rasters merging, 
#     the name of which is what is indicated in the csv file in the column "habitat_descriptor_shortName"
#   - A probability raster, result of all the class probability rasters merging,
#     the name of which is what is indicated in the csv file in the column  "habitat_descriptor_shortName" + "_proba"
#   - A confidence based on probability raster, result of all the class confidence based on probability rasters merging, 
#     the name of which is what is indicated in the csv file in the column "habitat_descriptor_shortName" + "_proba_confidence"
#   - An overall confidence raster, result of all the class overall confidence raster rasters merging, 
#     the name of which is what is indicated in the csv file in the column 
#     "habitat_descriptor_shortName" + "_combined_confidence"




#-------------------------- script parametrisation ---------------------------------------------------
#chunk size and max memory to load for processings
chunksize<-1e+09
maxmemory<-1e+09

#working directory
workingDirectory<-"E:/travail/EUSeaMap4/WP1/modeling/script_publication/hands_on_dataset/partialBoB"

#the config csv file that describes all the bits that have to be merged
config_csvFileName<-"habitat_descriptor_combining_confidence.csv"

outputConfidenceraster<-"overall_confidence.tif"


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

#set raster options (most important are chunksize and maxmemory)
rasterOptions(tmpdir=paste(getwd(),"/temp",sep=""),progress="text",format="GTiff",
              chunksize=chunksize,
              maxmemory=maxmemory)

suppressWarnings(checkList<-check_Inputs())
hasError<-checkList[[1]]


if (! hasError) {
  print("combining the confidence rasters...")
  #as.is=TRUE because otherwise string columns are considered as factors, and this make the script fail
  units_to_be_merged <- read.table(file=file.path("config_files",config_csvFileName), header=TRUE, sep=",", as.is = TRUE)
  
  filesToCombine<-file.path (units_to_be_merged$folder,units_to_be_merged$confidence_fileName)
  rastersToCombine<-stack(filesToCombine)
  final_confidence_raster<-min(rastersToCombine)
  
  writeRaster(final_confidence_raster,file.path("output",outputConfidenceraster),datatype='INT1U',overwrite=TRUE)
  
  rm(final_confidence_raster,rastersToCombine)
  
  #remove all tmp file created by the raster package while overlaying ou writing chunk by chunk
  file.remove(Sys.glob(file.path("temp", "r_tmp*")))
  end <- Sys.time()
  print (paste("Completed in",difftime(end,start),"!"))
} else print (checkList[[2]])



gc()





