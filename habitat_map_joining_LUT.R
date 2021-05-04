library(rgdal)
library(readxl)


#clear away all variables
rm(list= ls())
gc()

#Author: Mickaël Vasquez
#Version: 20210417_01
#Date: 17/04/2021
#Joins the habitat shapefile where habitats are described by a code to a look-up table that crosswalks each habitat
#code with various habitat classification
#the script assumes that the input habitat shapefile is in the folder 'output' of the working directoty

#-------------------------- script parametrisation ---------------------------------------------------
#working directory (i.e. the directory that contains the habitat shapefile)
workingDirectory<-"E:/travail/EUSeaMap4/WP1/modeling/script_publication/hands_on_dataset/partialBoB"


#Name of habitat shapefile
habitat_shapefile_name<-"habitat_polygon"
#name of the shapefile column that contains the habitat codes
modelled_habitat_code_column<-"hab_code"

#path for the look-up table (the table that crosswalks EUSeaMap Numeric codes and classification classes such as EUNIS)
LUT_table_excel_file<-"input/LUT_ATLPOL20190426.xlsx"
#name of the look-up table column that describes the habitat codes
LUT_habitat_code<-"ModelCode"
#name of the shapefile that will be created
output_habitat_shapefile_name<-"final_habitat_polygon"
#-------------------------- end script parametrisation ---------------------------------------------------

########################################## FUNCTIONS ##############################################################################
check_Inputs<-function () {
  #checks that all the input files exist
  
  hasError<-FALSE
  msg<-""
  file<-sprintf("output/%s.shp",habitat_shapefile_name)
  if (!file.exists(file)) {
    msg<-sprintf("script aborted: file %s does not exist",
                 file)
    hasError<-TRUE
  } else {
    if (!file.exists(LUT_table_excel_file)) {
      msg<-sprintf("script aborted: file %s does not exist",
                   LUT_table_excel_file)
      hasError<-TRUE
    }
  }
  
  lstCheckRes<-list(hasError,msg)
}

########################################### MAIN SCRIPT ##############################################################################

start <- Sys.time()

setwd(workingDirectory)

suppressWarnings(checkList<-check_Inputs())
hasError<-checkList[[1]]


if (! hasError) {

  habitat_sdf <-readOGR(file.path(getwd(),"output"),habitat_shapefile_name,stringsAsFactors=FALSE)
  lut_df<-read_excel(LUT_table_excel_file)
  new_habitat_sdf<-merge(x = habitat_sdf, y = lut_df, by.x = modelled_habitat_code_column, by.y = LUT_habitat_code)
  names(new_habitat_sdf)[names(new_habitat_sdf)==modelled_habitat_code_column]<-LUT_habitat_code
  #only the fields that are in the LUT are kept in the final table
  new_habitat_sdf<-new_habitat_sdf[,(names(new_habitat_sdf) %in% names(lut_df))]
  
  writeOGR(new_habitat_sdf,file.path(getwd(),"output"),output_habitat_shapefile_name,driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  end <- Sys.time()
  print (paste("Completed in",difftime(end,start),"!"))
} else print (checkList[[2]])
