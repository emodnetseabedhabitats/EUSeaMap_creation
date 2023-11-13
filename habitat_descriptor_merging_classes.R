library(terra)


#clear away all variables
rm(list= ls())
gc()



#Author: Mickael Vasquez
#Version: 20230414_01
#Date: 14/04/2023
#N.B.: main update compared to previous version: the use of terra package instead of raster package
# 14/04/2023: a bug was corrected
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
#     the name of which is what is indicated in the csv file in the column "habitat_descriptor_shortName" + "_confidence_based_on_proba"
#   - An overall confidence raster, result of all the class overall confidence raster rasters merging, 
#     the name of which is what is indicated in the csv file in the column 
#     "habitat_descriptor_shortName" + "_confidence_overall"




#-------------------------- script parametrisation ---------------------------------------------------
#max memory to use (in GBytes)
maxmemory<-20

#working directory
workingDirectory<-"E:/travail/euseamap_phase5/WP1/modeling/models/Caribbean"

#the config csv file that describes all the bits that have to be merged
config_csvFileName<-"habitat_descriptor_merging_inputs.csv"

#Which outputs are required?
#The habitat descriptor spatial distribution raster?
output_habitat_descriptor_raster<-TRUE
#The overall confidence raster?
output_overall_confidence<-TRUE
#The confidence raster based on probability?
output_confidence_based_on_proba<-TRUE
#The probability raster?
output_probability_raster<-FALSE

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
    #columns to check because they may contain raster file names
    ColumnsNames<-c("class_fileName",
                    "based_on_proba_confidence_fileName",
                    "overall_confidence_fileName",
                    "proba_fileName"
                    )
    #The file names contained in the above columns will be checked only the corresponding output is required
    #(see booleans in script parametrisation)
    outputBooleans<-c(output_habitat_descriptor_raster,
                      output_confidence_based_on_proba,
                      output_overall_confidence,
                      output_probability_raster
    )
    columnsToCheck<-c()
    cpt<-0
    for (i in 1:length(outputBooleans)) {
      if (outputBooleans[i]) {
        cpt<-cpt+1
        columnsToCheck[cpt]<-ColumnsNames[i]
      }
    }
    for (i in 1:nrow(units_to_be_merged)) {
      folder<-units_to_be_merged$folder[i]
      
      for (j in 1:length(columnsToCheck)) {
        file<-units_to_be_merged[[columnsToCheck[j]]][i]
        #if no file indicated, continue
        if (file=="") next
        
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


merge_habitat_descriptor_units<-function(units_to_be_merged,habitatDescriptorShortName) {
  
  
  units_to_be_merged<-subset(units_to_be_merged,units_to_be_merged$habitat_descriptor_shortName==habitatDescriptorShortName)
  number_of_units<-nrow(units_to_be_merged) 
  
  lst_class_raster_to_merge<-list()
  lst_confidence_in_threshold_raster_to_merge<-list()
  lst_combined_confidence_raster_to_merge<-list()
  lst_proba_raster_to_merge<-list()
  for (i in 1:number_of_units) {
    folderName<-units_to_be_merged[i,"folder"]
    if (output_habitat_descriptor_raster) {
      fileName<-units_to_be_merged[i,"class_fileName"]
      if (fileName!="") {
        if (!file.exists(file.path(folderName,fileName))){
          print(paste("No file",fileName,"in",folderName))
        }
        rst<-rast(file.path(folderName,fileName))
        lst_class_raster_to_merge[[length(lst_class_raster_to_merge)+1]]<-rst
      }
    }
    if (output_confidence_based_on_proba) {
      fileName<-units_to_be_merged[i,"based_on_proba_confidence_fileName"]
      if (fileName!="") {
        if (!file.exists(file.path(folderName,fileName))){
          print(paste("No file",fileName,"in",folderName))
        }
        rst<-rast(file.path(folderName,fileName))
        lst_confidence_in_threshold_raster_to_merge[[length(lst_confidence_in_threshold_raster_to_merge)+1]]<-rst
      }
    } 
    if (output_overall_confidence) {
      fileName<-units_to_be_merged[i,"overall_confidence_fileName"]
      if (fileName!="") {
        if (!file.exists(file.path(folderName,fileName))){
          print(paste("No file",fileName,"in",folderName))
        }
        rst<-rast(file.path(folderName,fileName))
        lst_combined_confidence_raster_to_merge[[length(lst_combined_confidence_raster_to_merge)+1]]<-rst
      }
    } 
    if (output_probability_raster) {
      fileName<-units_to_be_merged[i,"proba_fileName"]
      if (fileName!="") {
        if (!file.exists(file.path(folderName,fileName))){
          print(paste("No file",fileName,"in",folderName))
        }
        rst<-rast(file.path(folderName,fileName))
        lst_proba_raster_to_merge[[length(lst_proba_raster_to_merge)+1]]<-rst
      }
    }
  }
  #Merge class files
  if (output_habitat_descriptor_raster) {
    print ("   Merging class files")
    raster_merged_units <- do.call(merge, lst_class_raster_to_merge)
    outFile<-sprintf("%s.tif",habitatDescriptorShortName)
    print (sprintf("   writing %s",outFile))
    writeRaster(raster_merged_units,file.path("output",outFile),datatype='INT1U',overwrite=TRUE)
    rm(raster_merged_units,lst_class_raster_to_merge)
    gc()
  }
  #Merge confidence in threshold files
  if (output_confidence_based_on_proba) {
    print ("   Merging confidence based on probability")
    raster_merged_units <- do.call(merge, lst_confidence_in_threshold_raster_to_merge)
    outFile<-sprintf("%s_confidence_based_on_proba.tif",habitatDescriptorShortName)
    print (sprintf("   writing %s",outFile))
    writeRaster(raster_merged_units,file.path("output",outFile),datatype='INT1U',overwrite=TRUE)
    rm(raster_merged_units,lst_confidence_in_threshold_raster_to_merge)
    gc()
  }
  #Merge combined confidence files
  if (output_overall_confidence) {
    print ("   Merging overall confidence files")
    raster_merged_units <- do.call(merge, lst_combined_confidence_raster_to_merge)
    outFile<-sprintf("%s_confidence_overall.tif",habitatDescriptorShortName)
    print (sprintf("   writing %s",outFile))
    writeRaster(raster_merged_units,file.path("output",outFile),datatype='INT1U',overwrite=TRUE)
    rm(raster_merged_units,lst_combined_confidence_raster_to_merge)
    gc()
  }
  #Merge probability files
  if (output_probability_raster) {
    number_of_proba_rasters<-length(lst_proba_raster_to_merge)
    if (number_of_proba_rasters!=0) {
      print ("   Merging probability files")
      raster_proba<-lst_proba_raster_to_merge[[1]]
      if (number_of_proba_rasters>1) {
        for (i in 2:length(lst_proba_raster_to_merge)) {
          r<-lst_proba_raster_to_merge[[i]]
          raster_proba<-lapp(sds(raster_proba,r),
                             fun=function(r1,r2){ifelse(r1>r2,r1,r2)})
        }
      }
      outFile<-sprintf("%s_proba.tif",habitatDescriptorShortName)
      print (sprintf("   writing %s",outFile))
      writeRaster(raster_proba,file.path("output",outFile),datatype='FLT4S',overwrite=TRUE)
      rm(raster_proba,lst_proba_raster_to_merge)
      gc()
    }
  }
}


########################################### MAIN SCRIPT ##############################################################################

start <- Sys.time()

setwd(workingDirectory)

if (!file.exists("output")){
  dir.create("output")
}

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
  
  habitatDescriptorShortNames<-unique(units_to_be_merged$habitat_descriptor_shortName)
  
  for (i in 1:length(habitatDescriptorShortNames)) {
    habitatDescriptorShortName<-habitatDescriptorShortNames[i]
    print(paste("Merging",habitatDescriptorShortName,"..."))
    merge_habitat_descriptor_units (units_to_be_merged,habitatDescriptorShortName)
    #remove all tmp file created by the terra package
    tmpFiles(remove=T)
  }
  
  end <- Sys.time()
  print (paste("Completed in",difftime(end,start),"!"))
} else print (checkList[[2]])



gc()





