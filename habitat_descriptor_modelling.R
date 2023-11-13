library(terra)


#clear away all variables
rm(list= ls())
gc()


#Author: Mickael Vasquez
#Version: 202305126_01
#Date: 26/05/2023
#N.B.: main update compared to previous version: 
#- use of terra package instead of raster
#- integration of distance from manually drawn boundary equation type
#Creation of habitat descriptor classes according to a set of inputs (all described in a standardised csv configuration file)
#For each boundary of the habitat descriptor class, the csv file provides the following
#   - the name of the environmental variable raster(e.g. "temperature.tif")
#   - Equation type  (1=fuzzy, 2=glm) in the form of a constant value or, if vary spatially, of a raster name
#   - Fuzzy/GLM functions slope and intercept, in the form of a constant value or, if these parameters vary spatially, of a raster name
#   - probability threshold in the form of a constant value or, if this threshold varies spatially, of a raster name
#   - probability max probability in the form of a constant value or, if this max probability varies spatially, of a raster name
#   - the name of the confidence in the variable values raster(e.g. "confidence_temperature.tif"), where pixels have a value in (1,2,3), ie (low, moderate, high)
#As outputs the script creates for each habitat descriptor class (i.e. each row of the csv file):
#   - a probability raster, where the file has the name provided in the column "shortname" of the csv file
#   - a spatial distribution raster, where pixels have the value provided in the column "code" of the csv file
#   - a confidence based on probability raster, where pixels have a value in (1,2,3), ie (low, moderate, high)
#   - an overall confidence, combination of the confeidence in the threshold values and the variable values raster, 
#     where pixels have a value in (1,2,3), ie (low, moderate, high)


#-------------------------- script parametrisation ---------------------------------------------------
#max memory to use (in GBytes)
maxmemory<-20

#working directory

workingDirectory<-"E:/travail/euseamap_phase5/WP1/modeling/models/Baltic"

#the config csv file that describes all the habitat descriptor classes that have to be created
config_csvFileName<-"Baltic_habitat_descriptor_modelling_biozones.csv"

habitat_descriptor_probability_rasters_as_output<-TRUE
#-------------------------- end script parametrisation ---------------------------------------------------


########################################### FUNCTIONS ##############################################################################

check_habitat_descriptor_Inputs<-function (habitatDescriptors) {
  
  columnsToCheck<-c("upper_boundary_variable",
                    "upper_boundary_equation_type",
                    "upper_boundary_confidence_in_variable",
                    "upper_boundary_slope",
                    "upper_boundary_intercept",
                    "upper_boundary_threshold",
                    "lower_boundary_variable",
                    "lower_boundary_equation_type",
                    "lower_boundary_confidence_in_variable",
                    "lower_boundary_slope",
                    "lower_boundary_intercept",
                    "lower_boundary_threshold"
                    )
  hasError<-FALSE
  msg<-""
  for (i in 1:length(columnsToCheck)) {
    files<-unique(habitatDescriptors[[columnsToCheck[i]]])
    files<-files[files!="" & is.na(as.numeric(files)) & !is.na(files)]
    if (length(files)!=0) {
      for (j in 1:length(files)) {
        if (!file.exists(file.path("input",files[j]))) {
          msg<-sprintf("Error:  The script had to be terminated due to file %s, indicated in the column %s of the csv file, not existing in the input folder",
                       files[j],
                       columnsToCheck[i])
          hasError<-TRUE
          break
        }
      }
      if (hasError) break
    }
  }
  lstCheckRes<-list(hasError,msg)
}

calculate_probability<-function(boundary_equation_type,boundary_variable_file,
                                boundary_slope,boundary_intercept) {
  
  equation_type_vary_spatially<-is.na(as.numeric(boundary_equation_type))
  equation_parameters_vary_spatially<-is.na(as.numeric(boundary_slope))
  
  boundary_threshold<-habitatDescriptors[i,paste(strColumnSuffix,"threshold",sep="_")]
  threshold_vary_spatially<-is.na(as.numeric(boundary_threshold))
  
  #make a raster with the raster file
  raster_variable<-rast(file.path("input",boundary_variable_file))
  #calculation of the probability
  #first calculate ax+b
  if (equation_parameters_vary_spatially) {
    #the upper boundary parameters (slope, intercept, threshold) vary spatially
    #so raster files are defined in the config file and we use the lapp function for calculation
    boundary_slope<-rast(file.path("input",boundary_slope))
    boundary_intercept<-rast(file.path("input",boundary_intercept))
    raster_ax_plus_b<-lapp(sds(raster_variable,boundary_slope,boundary_intercept),
                           fun=function(x, a, b){((a *x ) + b)})
  } else {
    #the upper boundary parameters (slope, intercept, threshold) are constant values
    #so we use app function (instead of lapp) for calculation
    boundary_slope<-as.numeric(boundary_slope)
    boundary_intercept<-as.numeric(boundary_intercept)
    raster_ax_plus_b<-app(raster_variable,
                          fun=function(x){((boundary_slope * x) + boundary_intercept)})
  }
  rm(boundary_slope,boundary_intercept,raster_variable)
  gc()
  #then if it's a fuzzy rule, values of ax+b that are <0 and >1 have to be set to resp. 0 and 1
  #if it's a glm, the probablity is exp(ax+b)/(1+exp(ax+b))
  if (equation_type_vary_spatially) {
    #the equation type vary spatially, so a raster file is defined in the config file. The lapp function is used
    boundary_equation_type<-rast(file.path("input",boundary_equation_type))
    raster_proba<-lapp(sds(raster_ax_plus_b,boundary_equation_type),
                       fun=function(x,type){ifelse(type==1,ifelse(x<0,0,ifelse(x>1,1,x)),ifelse(type==3,ifelse(x>1,1,x),exp(x)/(1+exp(x))))})
  } else {
    #the equation type does not vary spatially; so a constant value is defined in the config file. The app function is used
    if (boundary_equation_type==1) {
      #fuzzy
      raster_proba<-app(raster_ax_plus_b,
                        fun=function(x){ifelse(x<0,0,ifelse(x>1,1,x))})
    } else if (boundary_equation_type==3) {
      #distance to manually drawn boundary
      raster_proba<-app(raster_ax_plus_b,
                        fun=function(x){ifelse(x>1,1,x)})
    } else {
      #GLM
      raster_proba<-app(raster_ax_plus_b,
                        fun=function(x){exp(x)/(1+exp(x))})
    }
  }
  #masking probability values that are not within the class probability range
  if (threshold_vary_spatially) {
    boundary_threshold<-rast(file.path("input",boundary_threshold))
    raster_proba<-lapp(sds(raster_proba,boundary_threshold),
                       fun=function(p,threshold){ifelse(p < threshold,NA,p)})
  }else {
    boundary_threshold<-as.numeric(boundary_threshold)
    if (boundary_equation_type==3) {
      #distance to manually drawn boundary: set null when = threshold
      raster_proba<-app(raster_proba,
                        fun=function(p){ifelse(p == boundary_threshold,NA,p)})
    } else {
      #fuzzy or GLM: set null when < threshold
      raster_proba<-app(raster_proba,
                        fun=function(p){ifelse(p < boundary_threshold,NA,p)})
    }
  }
  rm (raster_ax_plus_b,boundary_equation_type,boundary_threshold)
  return (raster_proba)
  
  
}


create_habitat_descriptor_class_raster<-function(raster_proba,HD_code) {
  
  raster_boundary_class<-rast(extent=ext(raster_proba),resolution=res(raster_proba),crs=crs(raster_proba))
  values(raster_boundary_class)<-HD_code
  raster_boundary_class<-mask(raster_boundary_class,raster_proba)
  
  return(raster_boundary_class)
}


create_Confidence_rasters<-function(raster_proba,strColumnSuffix) {
  
  boundary_threshold<-habitatDescriptors[i,paste(strColumnSuffix,"threshold",sep="_")]
  threshold_vary_spatially<-is.na(as.numeric(boundary_threshold))
  boundary_max_probability<-habitatDescriptors[i,paste(strColumnSuffix,"max_probability",sep="_")]
  max_probability_vary_spatially<-is.na(as.numeric(boundary_max_probability))

  print("         Classifying confidence based on probability")
  
  if (threshold_vary_spatially) {
    boundary_threshold<-rast(file.path("input",boundary_threshold))
    if (max_probability_vary_spatially) {
      boundary_max_probability<-rast(file.path("input",boundary_max_probability))
      probaRange<-lapp(sds(boundary_threshold,boundary_max_probability),
                       fun=function(threshold,max_proba){max_proba-threshold})
    } else {
      boundary_max_probability<-as.numeric(boundary_max_probability)
      probaRange<-app(boundary_threshold,
                      fun=function(threshold){boundary_max_probability-threshold})
    }
    raster_boundary_confidence<-lapp(sds(raster_proba,probaRange,boundary_threshold),
                                     fun=function(p,probaRange,threshold){ifelse(p>=threshold+0.6*probaRange,3,ifelse(p<threshold+0.2*probaRange,1,2))})
  }else {
    boundary_threshold<-as.numeric(boundary_threshold)
    if (max_probability_vary_spatially) {
      boundary_max_probability<-rast(file.path("input",boundary_max_probability))
      probaRange<-app(boundary_max_probability,
                      fun=function(max_proba){max_proba-boundary_threshold})
      raster_boundary_confidence<-lapp(sds(raster_proba,probaRange),
                                       fun=function(p,probaRange){ifelse(p>=boundary_threshold+0.6*probaRange,3,ifelse(p<boundary_threshold+0.2*probaRange,1,2))})
    } else{
      boundary_max_probability<-as.numeric(boundary_max_probability)
      probaRange<-boundary_max_probability-boundary_threshold
      highConf_bound<-boundary_threshold+0.6*probaRange
      lowConf_bound<-boundary_threshold+0.2*probaRange
      if (boundary_equation_type==3) {
        #distance to manually drawn boundary: set null when = threshold
        raster_boundary_confidence<-app(raster_proba,
                                        fun=function(p){ifelse(p>=highConf_bound,3,ifelse(p<lowConf_bound,1,2))})
      } else {
        #fuzzy or GLM: set null when < threshold
        raster_boundary_confidence<-app(raster_proba,
                                        fun=function(p){ifelse(p>=highConf_bound,3,ifelse(p<lowConf_bound,1,2))})
      }
    }
  }
  rm(probaRange,boundary_threshold,boundary_max_probability)
  gc()
  print("         intersecting confidence based on probability and confidence in variable")
  raster_confidence_in_variable<-rast(file.path("input",
                                                habitatDescriptors[i,paste(strColumnSuffix,"confidence_in_variable",sep="_")]))
  
  raster_combined_confidence<-lapp(sds(raster_boundary_confidence,raster_confidence_in_variable),
                                   fun=function(ct,cv){ifelse(ct==2 & cv ==1,1,
                                                       ifelse(ct==3 & cv ==1,2,
                                                       ct))})
 
 
  lstRasters<-list("raster_boundary_confidence"=raster_boundary_confidence,
                   "raster_combined_confidence"=raster_combined_confidence)
  
  rm(raster_confidence_in_variable)
  gc()
  
  return(lstRasters)
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

#open the configuration table
if (!file.exists(file.path("config_files",config_csvFileName))) print(paste("No file",config_csvFileName,"in config_files folder"))

#as.is=TRUE because otherwise string columns are considered as factors, and this make the script fail
habitatDescriptors <- read.table(file=file.path("config_files",config_csvFileName), header=TRUE, sep=",", as.is = TRUE)

suppressWarnings(checkList<-check_habitat_descriptor_Inputs(habitatDescriptors))
hasError<-checkList[[1]]
if (! hasError) {

  number_of_HDs<-nrow(habitatDescriptors)
  for (i in 1:number_of_HDs) {
    full_HD_name<-habitatDescriptors[i,"longName"]
    HD_code<-habitatDescriptors[i,"code"]
    if (is.character(habitatDescriptors[i,"upper_boundary_equation_type"])) {
      hasUpperBoundary<-(habitatDescriptors[i,"upper_boundary_equation_type"]!="")
    } else {
      hasUpperBoundary<-!is.na(as.numeric(habitatDescriptors[i,"upper_boundary_equation_type"]))
    }
    if (is.character(habitatDescriptors[i,"lower_boundary_equation_type"])) {
      hasLowerBoundary<-(habitatDescriptors[i,"lower_boundary_equation_type"]!="")
    } else {
      hasLowerBoundary<-!is.na(as.numeric(habitatDescriptors[i,"lower_boundary_equation_type"]))
    }
    
    #Processing upper boundary
    if (hasUpperBoundary) {
      #The habitat descriptor has a upper boundary
      strBoundaryType<-"upper boundary"
      strColumnSuffix<-"upper_boundary"
      ############## PROBABILITY CALCULATION ###################################################################################
      print(paste("Calculating probability raster for",full_HD_name, strBoundaryType))
      #get upper boundary equation type  (1=fuzzy, 2=glm)
      boundary_equation_type<-habitatDescriptors[i,paste(strColumnSuffix,"equation_type",sep="_")]
      #get the boundary parameters (equation slope, intercept, and threshold)
      boundary_variable_file<-habitatDescriptors[i,paste(strColumnSuffix,"variable",sep="_")]
      boundary_slope<-habitatDescriptors[i,paste(strColumnSuffix,"slope",sep="_")]
      boundary_intercept<-habitatDescriptors[i,paste(strColumnSuffix,"intercept",sep="_")]
      #calculate probability
      raster_proba<-calculate_probability (boundary_equation_type,boundary_variable_file,
                                           boundary_slope,boundary_intercept)
      gc()
      ############## CLASSIFICATION & CONFIDENCE ##############################################################################################
      print(paste("Creating habitat descriptor spatial distribution raster for",full_HD_name, strBoundaryType))
      raster_boundary_class<-create_habitat_descriptor_class_raster(raster_proba,HD_code)
      if (!(hasLowerBoundary & hasUpperBoundary)) {
        #has only a upper boundary. Write raster
        writeRaster(raster_boundary_class,file.path("output",paste(habitatDescriptors[i,"shortName"],".tif",sep="")),datatype='INT1U', overwrite=TRUE)
      } else {
        raster_upper_boundary_class<-raster_boundary_class
      } 
      rm(raster_boundary_class)
      print(paste("Creating confidence rasters for",full_HD_name, strBoundaryType))
      lstClassifiedRasters<-create_Confidence_rasters(raster_proba,strColumnSuffix)
      raster_boundary_confidence<-lstClassifiedRasters$raster_boundary_confidence
      raster_combined_confidence<-lstClassifiedRasters$raster_combined_confidence
      if (!(hasLowerBoundary & hasUpperBoundary)) {
        #has only a upper boundary. Write rasters
        writeRaster(raster_boundary_confidence,file.path("output",paste(habitatDescriptors[i,"shortName"],"confidence_based_on_proba.tif",sep="_")),datatype='INT1U',overwrite=TRUE)
        writeRaster(raster_combined_confidence,file.path("output",paste(habitatDescriptors[i,"shortName"],"confidence_overall.tif",sep="_")),datatype='INT1U',overwrite=TRUE)
      }else{
        #writeRaster(raster_combined_confidence,file.path("temp",paste(habitatDescriptors[i,"shortName"],"upper_combined_confidence.tif",sep="_")),datatype='INT1U',overwrite=TRUE)
        raster_upper_boundary_conf<-raster_boundary_confidence
        raster_upper_combined_conf<-raster_combined_confidence
      }
      rm(raster_boundary_confidence,raster_combined_confidence)
      if (habitat_descriptor_probability_rasters_as_output) {
        if (!(hasLowerBoundary & hasUpperBoundary)) {
          #has only a upper boundary. Write raster
          writeRaster(raster_proba,file.path("output",paste(habitatDescriptors[i,"shortName"],"proba.tif",sep="_")),overwrite=TRUE)
        } else {
          raster_upper_boundary_proba<-raster_proba
        }
        rm(raster_proba)
      }
      gc()
    }
    
    if (hasLowerBoundary) {
      #The habitat descriptor has a upper boundary
      strBoundaryType<-"lower boundary"
      strColumnSuffix<-"lower_boundary"
      ############## PROBABILITY CALCULATION ###################################################################################
      print(paste("Calculating probability raster for",full_HD_name, strBoundaryType))
      #get boundary equation type  (1=fuzzy, 2=glm)
      boundary_equation_type<-habitatDescriptors[i,paste(strColumnSuffix,"equation_type",sep="_")]
      #get the boundary parameters (equation slope, intercept, and threshold)
      boundary_variable_file<-habitatDescriptors[i,paste(strColumnSuffix,"variable",sep="_")]
      boundary_slope<-habitatDescriptors[i,paste(strColumnSuffix,"slope",sep="_")]
      boundary_intercept<-habitatDescriptors[i,paste(strColumnSuffix,"intercept",sep="_")]
      #calculate probability
      raster_proba<-calculate_probability (boundary_equation_type,boundary_variable_file,
                                           boundary_slope,boundary_intercept)
      gc()
      ############## CLASSIFICATION & CONFIDENCE ###############################################################################################
      print(paste("Creating habitat descriptor spatial distribution raster for",full_HD_name, strBoundaryType))
      raster_boundary_class<-create_habitat_descriptor_class_raster(raster_proba,HD_code)
      if (!(hasLowerBoundary & hasUpperBoundary)) {
        #has only a upper boundary. Write raster
        writeRaster(raster_boundary_class,file.path("output",paste(habitatDescriptors[i,"shortName"],".tif",sep="")),datatype='INT1U', overwrite=TRUE)
      } else {
        raster_lower_boundary_class<-raster_boundary_class
      } 
      rm(raster_boundary_class)
      print(paste("Creating confidence rasters for",full_HD_name, strBoundaryType))
      lstClassifiedRasters<-create_Confidence_rasters(raster_proba,strColumnSuffix)
      raster_boundary_confidence<-lstClassifiedRasters$raster_boundary_confidence
      raster_combined_confidence<-lstClassifiedRasters$raster_combined_confidence
      if (!(hasLowerBoundary & hasUpperBoundary)) {
        #has only a upper boundary. Write rasters
        writeRaster(raster_boundary_confidence,file.path("output",paste(habitatDescriptors[i,"shortName"],"confidence_based_on_proba.tif",sep="_")),datatype='INT1U',overwrite=TRUE)
        writeRaster(raster_combined_confidence,file.path("output",paste(habitatDescriptors[i,"shortName"],"confidence_overall.tif",sep="_")),datatype='INT1U',overwrite=TRUE)
      }else{
        #writeRaster(raster_combined_confidence,file.path("temp",paste(habitatDescriptors[i,"shortName"],"lower_combined_confidence.tif",sep="_")),datatype='INT1U',overwrite=TRUE)
        raster_lower_boundary_conf<-raster_boundary_confidence
        raster_lower_combined_conf<-raster_combined_confidence
      }
      rm(raster_boundary_confidence,raster_combined_confidence)
      if (habitat_descriptor_probability_rasters_as_output) {
        if (!(hasLowerBoundary & hasUpperBoundary)) {
          #has only a upper boundary. Write raster
          writeRaster(raster_proba,file.path("output",paste(habitatDescriptors[i,"shortName"],"proba.tif",sep="_")),overwrite=TRUE)
        } else {
          raster_lower_boundary_proba<-raster_proba
        }
        rm(raster_proba)
      }
      gc()
    }
    gc()
    
    if (hasLowerBoundary & hasUpperBoundary) {
      #union of the two rasters
      print(paste("Merging",full_HD_name,"upper and lower boundary rasters"))
      #spatial distribution rasters. The 2 rasters are merged
      print("         merging upper and lower boundary class spatial distribution raster")
      raster_boundary_class<-min(raster_upper_boundary_class,raster_lower_boundary_class)
      writeRaster(raster_boundary_class,
                  file.path("output",paste(habitatDescriptors[i,"shortName"],".tif",sep="")),
                  datatype='INT1U',overwrite=TRUE)
      rm(raster_upper_boundary_class,raster_lower_boundary_class,raster_boundary_class)
      
      #confidence rasters
      #Confidence based on proba The min probability is kept
      print("         merging upper and lower boundary confidence based on probability raster")
      raster_boundary_conf<-min(raster_lower_boundary_conf,raster_upper_boundary_conf)
      #Merge overall confidence of lower and upper boundary. Complex rule
      #if the confidence based on proba of the upper boundary is 1 or 2, then the overall confidence of the lower boundary is kept
      #else (ie the confidence based on proba of the upper boundary is 3)
      #       if the confidence based on proba of the lower boundary is 1 or 2, then the overall confidence of the lower boundary is kept
      #       else (ie the confidence based on proba of the lower boundary is 3), then the overall confidence of the lower boundary and the overall confidence of the upper boundary are averaged and rounded (with 2.5=3)
      print("         merging upper and lower boundary overall confidence raster")
      raster_combined_conf<-lapp(sds(raster_lower_combined_conf,raster_upper_combined_conf,raster_lower_boundary_conf,raster_upper_boundary_conf),
                                 fun=function(r1,r2,r3,r4){ifelse(r3 %in% c(1,2),r1,ifelse(r4 %in% c(1,2),r2,floor(((r1+r2)/2)+0.5)))})
      writeRaster(raster_boundary_conf,file.path("output",paste(habitatDescriptors[i,"shortName"],"confidence_based_on_proba.tif",sep="_")),datatype='INT1U',overwrite=TRUE)
      writeRaster(raster_combined_conf,file.path("output",paste(habitatDescriptors[i,"shortName"],"confidence_overall.tif",sep="_")),datatype='INT1U',overwrite=TRUE)
      rm(raster_lower_combined_conf,raster_upper_combined_conf,raster_combined_conf,
         raster_lower_boundary_conf,raster_upper_boundary_conf,raster_boundary_conf)
      if (habitat_descriptor_probability_rasters_as_output) {
        #probability. The min probability is kept
        print("         merging upper and lower boundary probability raster")
        raster_boundary_proba<-min(raster_lower_boundary_proba,raster_upper_boundary_proba)
        writeRaster(raster_boundary_proba,file.path("output",paste(habitatDescriptors[i,"shortName"],"_proba.tif",sep="")),overwrite=TRUE)
        rm(raster_lower_boundary_proba,raster_upper_boundary_proba,raster_boundary_proba)
      }
      
      
    }
    #remove all tmp file created by the terra package
    tmpFiles(remove=T)
    gc()
  }
} else print (checkList[[2]])
end <- Sys.time()

print (paste("Completed in",difftime(end,start),"!"))

