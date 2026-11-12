
######################################################################################################################################################
# R script for the analyses of the paper entitled
# "Telomeres in a spatial context: a tool to understand ageing rate variation in wild populations" 
#
# Authors (The Pablo's Team):
#   Pablo Burraco
#   Pablo M. Lucas
#   Pablo Salmon
#
# For doubts/comments/questions about the code, please contact to Pablo M. Lucas at lucas.pablo.2010@gmail.com
######################################################################################################################################################

#Schema

#0 DATA PREPARATION
  #0.1  Read and transform Climatic database CHELSA  
  #0.2 Obtaintion of coordinates of our raster cells 
#1 MODELLING TELOMERE LONGITUDE
  #1.1 Wee load the sampling data
  #1.2 We cross the vector of sampling with the environmental variables for Europe
  #1.3 Descriptive analysis and variable selection
    #1.3.1 Descriptive analysis of the data
    #1.3.2 Calculation of univariable models for variable selection
  #1.4 Calculation of a null model with only the intercept
  #1.5 We compare univariable models with the null model
  #1.6 Check the correlation of variables with VIF 
  #1.7 Average modelling of telomere lenght 
    #1.7.1 Model evaluation
    #1.7.2 Prediction of telomere  
    #1.7.3 Spatial prediction of telomere lenght over Europe. Close R and open a new one in version 3.53
#2 HABITAT MODEL USING GBIF DATA
  #2.1 Download of occurrences from GBIF
  #2.2 We create a polygon and a raster with the study area/area of sampling from GBIF
  #2.3 We add to the data of the presences from GBIF with the study area, the data of ID of pixels and the environmetnal variables
  #2.4 Calculation of univariable models of habitat for variable selection  
  #2.5 Calculation of a null model with only the intercept
  #2.6 We compare univariable models with the null model
  #2.7 Check the correlation of variables
 #2.8 We are going to  use biomod2 
   #2.8.1 We formating the data for biomod
   #2.8.2 Modeling options 
   #2.8.3 Building models 
   #2.8.4 Ensemble Modeling
   #2.8.5 Ensemble Modeling Evaluation and variable importance  
   #2.8.6 Ensemble Forecasting of Current Scenario 
    ##2.9 We see the results of the models
      #2.9.1 We see the respose curve for each variable (Thanks for the code in this section to Prof. Antoine Guisan and Dr. Olivier Broennimann)   
      #2.9.2 We create a map with the predcition of habitat (Figure 2b)
#3 CORRELATION TELOMERE LENGHT WITH HABITAT
 #3.1 Validation of telomere length as presence/absences for species distribution

###################################################################################################################################################
#0 DATA PREPARATION
###################################################################################################################################################
  #We use R version 3.5.3
  rm(list=ls()) 
  ls()
  your_wd<-"WRITE HERE YOUR WORKING DIRECTORY FOLDER"          
  setwd(your_wd)
  #Load the required libraries
  install.packages("reshape")  
  install.packages("measurements")
  install.packages("questionr")
  install.packages("lme4")
  install.packages("lmerTest")
  install.packages("usdm")
  install.packages("MuMIn")
  install.packages("rlang")
  install.packages("dplyr")
  install.packages("akima")
  install.packages("rgl")
  install.packages("readxl")
  library(rgdal) 
  library(raster) 
  library(readxl)  
  library(measurements)  
  library(questionr) 
  library(lme4) 
  library(lmerTest) 
  library(usdm)
  library(MuMIn)
  library(rlang) 
  library(dplyr) 
  library(akima) 
  library(rgl) 
  library(boot) 

  #0.1  Read and transform Climatic database CHELSA  
    #We have previously downloaded the data from CHELSA database to a folder called "Original"
    folder_with_CHELSA_original<- paste0(your_wd,"/Original")
    list_files<-list.files(path=folder_with_CHELSA_original)
    DF_spatial_data<-c()
    files_list<-c()
     for (p in 1:19){ 
        file_in_loop<-list_files[p]
        print(paste0("########## Starting projection of Bioclimatic variable ",file_in_loop," ##########"))       
        r1 <- raster(paste(folder_with_CHELSA_original,"/",file_in_loop,sep="")) 
        crs(r1) <- CRS('+init=EPSG:4326')
        #We are going to cut the map for reduce to Europe
        e <- extent(-20,105,15,75) 
        m1_crop <- crop(r1,e) 
        #Raster of reference: 
        new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
        pr <- projectRaster(m1_crop, new_raster, method='bilinear')
        e2<- extent(-2800000,3100000,-800000,4800000) 
        m1_crope2 <- crop(pr,e2)
        #We save the raster:
          folder_to_save<-paste0(your_wd,"/Projected/")
          files_list<-cbind(files_list,file_in_loop)
          #We are going to create a new folder if it does not exist
          if (file.exists(folder_to_save)){
            assign(file_in_loop,pr)
            save(list=file_in_loop, file=paste(folder_to_save, file_in_loop, ".Rdata",sep=""))         
          } else {
            dir.create(file.path(your_wd, "Projected"))
            assign(file_in_loop,pr)
            save(list=file_in_loop, file=paste(folder_to_save, file_in_loop, ".Rdata",sep=""))       
          }
          vec_pr<-extract(pr,e2)
          DF_spatial_data<-cbind(DF_spatial_data,vec_pr) 
     }
    colnames(DF_spatial_data)<-files_list
    #We reorder the columns
    DF_spatial_data_2 <- DF_spatial_data[,c("CHELSA_bio10_1.tif","CHELSA_bio10_2.tif","CHELSA_bio10_3.tif",
                            "CHELSA_bio10_4.tif","CHELSA_bio10_5.tif","CHELSA_bio10_6.tif","CHELSA_bio10_7.tif",
                            "CHELSA_bio10_8.tif","CHELSA_bio10_9.tif","CHELSA_bio10_10.tif","CHELSA_bio10_11.tif",
                            "CHELSA_bio10_12.tif","CHELSA_bio10_13.tif","CHELSA_bio10_14.tif","CHELSA_bio10_15.tif",
                            "CHELSA_bio10_16.tif","CHELSA_bio10_17.tif","CHELSA_bio10_18.tif","CHELSA_bio10_19.tif")]
       
    Scenario_current<-as.data.frame(DF_spatial_data_2)
    colnames(Scenario_current)<-c("Clim_1","Clim_2","Clim_3","Clim_4","Clim_5","Clim_6","Clim_7","Clim_8","Clim_9","Clim_10","Clim_11","Clim_12","Clim_13","Clim_14","Clim_15","Clim_16","Clim_17","Clim_18","Clim_19")
        save(Scenario_current, file="Scenario_current.RData")

  #0.2 Obtain of coordinates of our raster cells 
        #We create a new raster (with the same resolution of spatial extent as variables over Europe)
        #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
        new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
        #We are going to use Europe Albers Equal Area Conic  
        newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
        #We defined the projection of our raster:
        projection(new_raster) <- newproj
        vec_values<-seq(1,33040000,1)
        values(new_raster)<-vec_values
        #We obtain a matrix with three columns (x, y and an ID for the pixel)
        data_matrix <- rasterToPoints(new_raster)
        colnames(data_matrix)<-c("x","y","ID_pixel")
        save(data_matrix, file="data_matrix.RData")

###################################################################################################################################################
#1 MODELLING TELOMERE LONGITUDE
###################################################################################################################################################
  #1.1 Wee load the sampling data. 
  #We load the data of telomere lenght for each subpopulation. 
    #We used information contained in Dupoué, A.et al. 2017. https://doi.org/10.1038/s41598-017-17323-z 
    #Coordinates from these locations can be extracted from Supplementary Table S1 in Dupoué, A. et al. 2018.  https://doi.org/10.1111/1365-2656.12843
    DF_pop <- read.csv(file="Coordinates_Dupoue_pop.csv",head=TRUE,sep=",")
    dms2dec <- function(dms, separators = c("º", "°", "\'", "\"")) {
    dms <- as.character(dms)
    dms <- gsub(pattern = " ", replacement = "", x = dms)
    for (s in separators) dms <- gsub(pattern = s, replacement = "_splitHere_", x = dms)
    splits <- strsplit(dms, split = "_splitHere_")
    n <- length(dms)
    deg <- min <- sec <- hem <- vector("character", n)
    for (i in 1:n) {
      deg[i] <- splits[[i]][1]
      min[i] <- splits[[i]][2]
      sec[i] <- splits[[i]][3]
      hem[i] <- splits[[i]][4]
    }
    dec <- as.numeric(deg) + (as.numeric(min) / 60) + (as.numeric(sec) / 3600)
    sign <- ifelse (hem %in% c("N", "E"), 1, -1)
    dec <- sign * dec
    return(dec)
    }  # 
    DF_pop$latitude.decimal <- dms2dec(DF_pop$Latitude)
    DF_pop$longitude.decimal <- dms2dec(DF_pop$Longitude)
    DF_coordenates<-DF_pop[,c("longitude.decimal","latitude.decimal")]     
    DF_coordenates<-rbind(BAS,BEL,CAR,CHA,COM,JOC,JON,LAJ,MON,PUY)    
    df_Populations<-data.frame(Populations=c("Barnesac","Belair","Caroux","Chalet","Col_du_cheval_mort","Source_Loire","Gerbier_de_Jonc","Lajo","Montselgues","Puy_Mary"))
    DF_coordenates<-cbind(DF_coordenates2,df_Populations,id=1:10)  
    colnames(DF_coordenates)<-c("Coor_X","Coor_Y","Populations","id")
    save(DF_coordenates, file="DF_coordenates.RData")
  
    #We transform the sampling points to spatial data
    pre_coor_Y<-as.numeric(DF_coordenates$Coor_Y)
    pre_coor_X<-as.numeric(DF_coordenates$Coor_X)
    mapa_presences<-as.data.frame(cbind(pre_coor_X,pre_coor_Y))
    class(mapa_presences)
    coordinates(mapa_presences)=~pre_coor_X+pre_coor_Y
    class(mapa_presences)
    #Metadata indicate that are in WGS84(epsg:4326)
    proj4string(mapa_presences) = CRS("+init=epsg:4326")
    mapa_mapa_presences_df_WGS=SpatialPointsDataFrame(mapa_presences, data.frame(id=1:length(mapa_presences)))
    #We have created a spatial point object with the projection  WGS84: mapa_sub_coord_uncert_564_df and we are going to define other projection and to project in this new one  
    #We are going to use Europe Albers Equal Area Conic  
    newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
    #We change the projection
    mapa_mapa_presences_aea <- spTransform(mapa_mapa_presences_df_WGS, newproj)
    #We save the raw occurrences proyected in AEA in a point shapefile       
    writeOGR(obj=mapa_mapa_presences_aea, "Sampling_pop",layer="Sampling_pop", "ESRI Shapefile",overwrite_layer=TRUE)
          #Raster of reference: 
          new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
          #We defined the projection of our raster:
          projection(new_raster) <- newproj
          #Rasterize the data
          sps_raster<-rasterize(mapa_mapa_presences_aea, new_raster, field="id")
          e2<- extent(-2800000,3100000,-800000,4800000) 
          vec_sampling<-extract(sps_raster,e2)
          #We exported the raster to ascii format. The we can convert it in Arcgis 9.4 using converts form Ascii to raster
          writeRaster(sps_raster, 'Raster_sampling.img',overwrite=TRUE) 
          save(vec_sampling,file="vec_sampling.RData")

    #Raster of reference for extrapolate: 
        x_coord <- c(2,  2,  5, 5)
        y_coord <- c(42, 46, 46, 42)
        xym <- cbind(x_coord, y_coord)
        p = Polygon(xym)
        ps = Polygons(list(p),1)
        sps = SpatialPolygons(list(ps))
        proj4string(sps) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
        sps_EAEA <- spTransform(sps, newproj)
        new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
        #We defined the projection of our raster:
        projection(new_raster) <- newproj
        #Rasterize the data
        gbif_download_area_raster<-rasterize(sps_EAEA, new_raster)#, field="id"
        e2<- extent(-2800000,3100000,-800000,4800000) 
        vec_gbif_download_area_raster<-extract(gbif_download_area_raster,e2)
        sum(na.omit(vec_gbif_download_area_raster))#
        #We exported the raster to ascii format. The we can convert it in Arcgis 9.4 using converts form Ascii to raster
        writeRaster(vec_gbif_download_area_raster, 'vec_gbif_download_area_raster.img',overwrite=TRUE) 
        save(vec_gbif_download_area_raster,file="vec_gbif_download_area_raster.RData")

  #1.2 We cross the vector of sampling with the environmental variables for Europe
  #We use R version 4.0.2
    rm(list=ls()) 
    ls()
    setwd("WRITE HERE YOUR WORKING DIRECTORY FOLDER")
    
    #Load the required libraries
    install.packages("reshape")  
    install.packages("measurements")
    install.packages("questionr")
    install.packages("lme4")
    install.packages("lmerTest")
    install.packages("usdm")
    install.packages("MuMIn")
    install.packages("rlang")
    install.packages("dplyr")
    install.packages("akima")
    install.packages("rgl")
    install.packages("readxl")
    library(rgdal) 
    library(raster) 
    library(readxl)  
    library(measurements)  
    library(questionr) 
    library(lme4) 
    library(lmerTest) 
    library(usdm)
    library(MuMIn)
    library(rlang) 
    library(dplyr) 
    library(akima) 
    library(rgl) 
    library(shiny) 
    library(boot) 
    #Vector with the sampling areas 
    load("vec_sampling.RData")# 
    #Data of coordinates and Ids of pixels of the raster of reference Europe Albers Equal Area Conic 
    load("data_matrix.RData")# 
    #Scenario current (Current Climate variables download from CHELSA and projected to Europe Albers Equal Area Conic 
    load("Scenario_current.RData")# 
    #Vector with the area of prediction 
    load("vec_gbif_download_area_raster.RData")# 
    load("DF_coordenates.RData")
    DF_LC_current_data_presence<-as.data.frame(cbind(vec_sampling,data_matrix,Scenario_current,vec_gbif_download_area_raster))#,vec_density_with_NAs_scaled))
    save(DF_LC_current_data_presence, file="DF_LC_current_data_presence.RData")

    #We select the study area. Here We are going to apply the scale 
    DF_LC_current_data_presence_study_area <- subset(DF_LC_current_data_presence, vec_gbif_download_area_raster>0)
    DF_LC_current_data_presence_study_area$vec_sampling[is.na(DF_LC_current_data_presence_study_area$vec_sampling)]<-0
    DF_LC_current_data_presence_study_area_no_na<-na.omit(DF_LC_current_data_presence_study_area)

    DF_LC_current_data_presence_study_area_no_na$Clim_1_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_1/max(DF_LC_current_data_presence_study_area_no_na$Clim_1)
    DF_LC_current_data_presence_study_area_no_na$Clim_2_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_2/max(DF_LC_current_data_presence_study_area_no_na$Clim_2)
    DF_LC_current_data_presence_study_area_no_na$Clim_3_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_3/max(DF_LC_current_data_presence_study_area_no_na$Clim_3)
    DF_LC_current_data_presence_study_area_no_na$Clim_4_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_4/max(DF_LC_current_data_presence_study_area_no_na$Clim_4)
    DF_LC_current_data_presence_study_area_no_na$Clim_5_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_5/max(DF_LC_current_data_presence_study_area_no_na$Clim_5)
    DF_LC_current_data_presence_study_area_no_na$Clim_6_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_6/max(DF_LC_current_data_presence_study_area_no_na$Clim_6)
    DF_LC_current_data_presence_study_area_no_na$Clim_7_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_7/max(DF_LC_current_data_presence_study_area_no_na$Clim_7)
    DF_LC_current_data_presence_study_area_no_na$Clim_8_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_8/max(DF_LC_current_data_presence_study_area_no_na$Clim_8)
    DF_LC_current_data_presence_study_area_no_na$Clim_9_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_9/max(DF_LC_current_data_presence_study_area_no_na$Clim_9)
    DF_LC_current_data_presence_study_area_no_na$Clim_10_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_10/max(DF_LC_current_data_presence_study_area_no_na$Clim_10)
    DF_LC_current_data_presence_study_area_no_na$Clim_11_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_11/max(DF_LC_current_data_presence_study_area_no_na$Clim_11)
    DF_LC_current_data_presence_study_area_no_na$Clim_12_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_12/max(DF_LC_current_data_presence_study_area_no_na$Clim_12)
    DF_LC_current_data_presence_study_area_no_na$Clim_13_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_13/max(DF_LC_current_data_presence_study_area_no_na$Clim_13)
    DF_LC_current_data_presence_study_area_no_na$Clim_14_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_14/max(DF_LC_current_data_presence_study_area_no_na$Clim_14)
    DF_LC_current_data_presence_study_area_no_na$Clim_15_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_15/max(DF_LC_current_data_presence_study_area_no_na$Clim_15)
    DF_LC_current_data_presence_study_area_no_na$Clim_16_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_16/max(DF_LC_current_data_presence_study_area_no_na$Clim_16)
    DF_LC_current_data_presence_study_area_no_na$Clim_17_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_17/max(DF_LC_current_data_presence_study_area_no_na$Clim_17)
    DF_LC_current_data_presence_study_area_no_na$Clim_18_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_18/max(DF_LC_current_data_presence_study_area_no_na$Clim_18)
    DF_LC_current_data_presence_study_area_no_na$Clim_19_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_19/max(DF_LC_current_data_presence_study_area_no_na$Clim_19)

    DF_LC_current_data_presence_study_area_no_na$Clim_1_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_1_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_2_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_2_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_3_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_3_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_4_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_4_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_5_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_5_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_6_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_6_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_7_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_7_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_8_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_8_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_9_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_9_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_10_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_10_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_11_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_11_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_12_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_12_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_13_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_13_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_14_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_14_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_15_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_15_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_16_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_16_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_17_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_17_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_18_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_18_scal^2
    DF_LC_current_data_presence_study_area_no_na$Clim_19_c_scal<-DF_LC_current_data_presence_study_area_no_na$Clim_19_scal^2
    save(DF_LC_current_data_presence_study_area_no_na, file="DF_LC_current_data_presence_study_area_no_na.RData")
    #We select the areas of sampling the lizards
    data_sampling_areas <- subset(DF_LC_current_data_presence_study_area_no_na, vec_sampling>0)
    save(data_sampling_areas, file="data_sampling_areas.RData")
    #We load the data of telomere lenght for each subpopulation. 
    #We used information contained in Dupoué, A. et al. 2017. https://doi.org/10.1038/s41598-017-17323-z
    DF0 <- read.csv(file="DupoueData.csv",head=TRUE,sep=",")
    #We merge coordinates of each subpopulation with the data of telomere lenght for each subpopulation  
      #We merge the environmental variables of the sampled areas with DF_coordenates
      merge_sampling_coordinates<-merge(data_sampling_areas,DF_coordenates, by.x="vec_sampling", by.y = "id",all.x = T)
      #We merge the environmental variables of each subpopulation with the information of telomere lenght
      merge_sampling_coordinates_telomeres<-merge(DF0,merge_sampling_coordinates, by.x="Population", by.y = "Populations",all.x = T)
      merge_sampling_coordinates_telomeres$telomere_scal<-(merge_sampling_coordinates_telomeres$telomere..kb.-min(merge_sampling_coordinates_telomeres$telomere..kb.)+0.0001)/(max(merge_sampling_coordinates_telomeres$telomere..kb.)-min(merge_sampling_coordinates_telomeres$telomere..kb.)+0.001)
      save(merge_sampling_coordinates_telomeres, file="merge_sampling_coordinates_telomeres.RData")

  #1.3 Descriptive analysis and variable selection
    #1.3.1 Descriptive analysis of the data
      library(lme4)  
      library(lawstat)  
      library(e1071)  
      library(glmmTMB)
      library(boot)
      library(MuMIn)
      hist(merge_sampling_coordinates_telomeres$telomere..kb.,breaks = 11)
      #We chec the simetry of he data  
      symmetry.test(
        merge_sampling_coordinates_telomeres$telomere..kb.,
        option = c("MGG", "CM", "M"),
        side = c("both", "left", "right"),
        boot = TRUE,
        B = 1000,
        q = 8/9
      ) 
      #We asses the skewness and we seee that has negative skewness which indicates that the mean of the data values is less
      #than the median, and the data distribution is left-skewed.     
      skewness(merge_sampling_coordinates_telomeres$telomere..kb.)#-0.7754334 so it will be better fit to a beta distribution  
  
    #1.3.2 Calculation of univariable models for variable selection
      var_names<-colnames(merge_sampling_coordinates_telomeres)
      names_var_clim<-var_names[34:52]
      DF_AIC_univar_models_telomere<-c()
      
          for (v in 1:19){
            var_in_loop_clim<-var_names[c(42+v)]
            write_model_prefix_clim<-c("(telomere_scal)~")
            fixed_term_clim<-var_in_loop_clim
            write_model_var_in_loop_clim_1<-paste(var_in_loop_clim,sep="")# lineal effect
            write_model_var_in_loop_clim_2<-paste("I(",var_in_loop_clim,"^2",")",sep="")# quadratic effect
            write_model_var_in_loop_clim_3<-paste(var_in_loop_clim," + ","I(",var_in_loop_clim,"^2",")",sep="")# lineal effect and quadratic effect
        
            write_model_1_clim<-paste(write_model_prefix_clim,write_model_var_in_loop_clim_1,"+(1|Population)")# lineal effect
            write_model_2_clim<-paste(write_model_prefix_clim,write_model_var_in_loop_clim_2,"+(1|Population)")# quadratic effect
            write_model_3_clim<-paste(write_model_prefix_clim,write_model_var_in_loop_clim_3,"+(1|Population)")# lineal effect and quadratic effect
          
            from_mod_1_clim<-as.formula(write_model_1_clim)
            from_mod_2_clim<-as.formula(write_model_2_clim)
            from_mod_3_clim<-as.formula(write_model_3_clim)
            
            mod_1_clim<- glmmTMB(formula = from_mod_1_clim , family=beta_family(link="logit"), data = merge_sampling_coordinates_telomeres)  
            mod_2_clim<- glmmTMB(formula = from_mod_2_clim , family=beta_family(link="logit"), data = merge_sampling_coordinates_telomeres)  
            mod_3_clim<- glmmTMB(formula = from_mod_3_clim , family=beta_family(link="logit"), data = merge_sampling_coordinates_telomeres)  
            
            #mod_null<- glm(formula = (sp_1_europa)~1, data = data_for_mod, family = "binomial", na.action = "na.fail") 
            AICc_mod_1_clim<-AICc(mod_1_clim)
            AICc_mod_2_clim<-AICc(mod_2_clim)
            AICc_mod_3_clim<-AICc(mod_3_clim)
            AICc_mod_clim_123<-cbind(AICc_mod_1_clim,AICc_mod_2_clim,AICc_mod_3_clim)
            AICc_mod_clim_123_min<-min(AICc_mod_clim_123[1,])
            AICc_in_loop<-cbind(AICc_mod_clim_123_min,AICc_mod_clim_123)
            DF_AIC_univar_models_telomere<-rbind(DF_AIC_univar_models_telomere,AICc_in_loop)
          }
          colnames(DF_AIC_univar_models_telomere)<-c("AICc_min","AICc_lin","AICc_quadr","AICc_quadr_lin")
          rownames(DF_AIC_univar_models_telomere)<-names_var_clim
          sorted_variables_clim<-sort(DF_AIC_univar_models_telomere[,1])
          DF_AIC_univar_models_telomere<-as.data.frame(DF_AIC_univar_models_telomere)
          #We select the 3 better
          save(sorted_variables_clim, file="sorted_variables_clim.RData")
          save(DF_AIC_univar_models_telomere, file="DF_AIC_univar_models_telomere.RData")# Supplementary Table A2
      

 #1.4 Calculation of a null model with only the intercept
    mod_null<-glmmTMB(telomere_scal~ (1 | Population), family=beta_family(link="logit"), data = merge_sampling_coordinates_telomeres)
    AICc(mod_null)
    
 #1.5 We compare univariable models with the null model 
    DF_AIC_univar_models_telomere$Delta_AIC_from_null<-DF_AIC_univar_models_telomere$AICc_min-AICc(mod_null)
    save(DF_AIC_univar_models_telomere, file="DF_AIC_univar_models_telomere.RData")
    write.csv(DF_AIC_univar_models_telomere, file="DF_AIC_univar_models_telomere.csv")
    #We select the best three variables based on AICc
    DF_AIC_univar_models_telomere_best_models <- subset(DF_AIC_univar_models_telomere, Delta_AIC_from_null<(-2.5))
    save(DF_AIC_univar_models_telomere_best_models, file="DF_AIC_univar_models_telomere_best_models.RData")

 #1.6 Check the correlation of variables
    cor_merge_sampling_coordinates_telomeres<-cor(merge_sampling_coordinates_telomeres[,c(34:52)])#,"Clim_11","Clim_7"
    write.csv(cor_merge_sampling_coordinates_telomeres, file="cor_merge_sampling_coordinates_telomeres.csv")# Supplementary Table A3
    #We retaain variables Clim_6_scal (whioh is the best based in AICc) and with Clim_7_scal (which is the second best uncorrelated variable). We include quadratic effects of boths

  #1.7 Average modelling of telomere lenght 
    #We write a full model which include all variables. #### el previo era de "Clim_2","Clim_11"
    fullmodel_no_NA<-glmmTMB(telomere_scal~ Clim_6_scal+Clim_6_c_scal+ Clim_7_scal+Clim_7_c_scal+(1 | Population), family=beta_family(link="logit"), data = merge_sampling_coordinates_telomeres)
    #fullmodel_no_NA<-glmmTMB(telomere_scal~ Clim_6_scal+Clim_6_c_scal+ (1 | Population), family=beta_family(link="logit"), data = merge_sampling_coordinates_telomeres)
    summary(fullmodel_no_NA)
    AICc(fullmodel_no_NA)
    #With dredge funtion we are going to obtain all possible models combining the variables included in the full model 
    dr_reduced_model_glmm<- dredge(fullmodel_no_NA, rank= "AICc")
    #load("dr_reduced_model_glmm.RData")
    save (dr_reduced_model_glmm, file='dr_reduced_model_glmm.RData')
    write.csv(dr_reduced_model_glmm, file="dr_reduced_model_glmm.csv")# Supplementary Table A4

    #We get the models subsetting the models which are equally good (delta<3)
    get_models_glmm<-get.models(dr_reduced_model_glmm, subset=delta < 5)
    save (get_models_glmm, file='get_models_glmm.RData')
    #WE calculate the average model# S3 method for model.selection
    avg_model<-model.avg(get_models_glmm)
    avg_model$ msTable
    avg_model$ coefficients
    summary(avg_model)
    save (avg_model, file='avg_model.RData')# Supplementary Table A5

    #1.7.1 Model evaluation
      pred_avg_model<-predict(avg_model,newdata=NULL)
      rsq <- function (x, y) cor(x, y) ^ 2
      rsq(pred_avg_model,merge_sampling_coordinates_telomeres$telomere_scal)# R2 included in Supplementary Table A5 0.7330269

    #1.7.2 Prediction of telomere  
      #We apply the average model obtained before using the mean values for all populations. We choose subset option
      # Prediction of random effects in linear and generalized linear models under model misspecification https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3066692/
      prediction_tel_study_area<-inv.logit(-10.116 -27.764*DF_LC_current_data_presence_study_area_no_na$Clim_6_c_scal -33.242*DF_LC_current_data_presence_study_area_no_na$Clim_6_scal+2.850*DF_LC_current_data_presence_study_area_no_na$Clim_7_scal+2.609*DF_LC_current_data_presence_study_area_no_na$Clim_7_c_scal)

      DF_LC_current_data_presence_study_area_no_na2<-cbind(DF_LC_current_data_presence_study_area_no_na,prediction_tel_study_area)
      ID_pixel_pred<-DF_LC_current_data_presence_study_area_no_na2[,c("ID_pixel","prediction_tel_study_area")]
      
      merge_sampling_coordinates_telomeres2<-merge(data_matrix,ID_pixel_pred, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
      vec_pred_over_Europe<-merge_sampling_coordinates_telomeres2[,c("prediction_tel_study_area")]
      save (vec_pred_over_Europe, file='vec_pred_over_Europe.RData')
      
    #1.7.3 Spatial prediction of telomere lenght over Europe. Close R and open a new one in version 3.53
        #We are going to create a map with the prediction
        #We use R version 3.5.3
        rm(list=ls()) 
        setwd("WRITE HERE YOUR WORKING DIRECTORY FOLDER")
        ls()
        library(rgdal) 
        library(raster) 
        library(readxl)  
        library(questionr) 
        library(lme4) 
        library(lmerTest) 
        library(usdm)
        library(MuMIn)
        library(rlang) 
        library(dplyr) 
        library(akima) 
        library(rgl) 
        library(boot) 
 
        setwd("WRITE HERE YOUR WORKING DIRECTORY FOLDER")
        load("vec_pred_over_Europe.RData")
           #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #Raster of reference: 
            telomer_predict <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
            #We defined the projection of our raster:
            projection(telomer_predict) <- newproj
            #Rasterize the data
            values(telomer_predict)<-vec_pred_over_Europe
            #We exported the raster 
            writeRaster(telomer_predict, 'telomer_predict.img',overwrite=TRUE) 
            save(telomer_predict,file="telomer_predict.RData")


###################################################################################################################################################
#2 HABITAT MODEL USING GBIF DATA
###################################################################################################################################################
  #We use R version 3.5.3
  rm(list=ls()) 
  ls()
  your_wd<-"WRITE HERE YOUR WORKING DIRECTORY FOLDER"          
  setwd(your_wd)
    #install.packages("rgbif")   
    install.packages("ggmap")
    install.packages("maps")
    install.packages("mapdata")
    install.packages("maptools")
    install.packages("readxl")
    #install_github("ropensci/CoordinateCleaner")
    #install.packages("CoordinateCleaner")
    #install.packages("devtools")
    #install.packages("countrycode")
    install.packages("rlang")
    packageurl <- "http://cran.r-project.org/src/contrib/Archive/rgbif/rgbif_2.0.0.tar.gz"
    install.packages(packageurl, repos=NULL, type="source")
    packageurl <- "http://cran.r-project.org/src/contrib/Archive/CoordinateCleaner/CoordinateCleaner_2.0-9.tar.gz"
    install.packages(packageurl, repos=NULL, type="source")

    library(ggmap)
    library(maps)
    library(ggplot2)
    library(mapdata)
    library(rgbif)     
    library(grid)  
    library(maptools)
    library(rgdal)
    library(sp)
    library(rgdal) 
    library(raster) 
    library(readxl) 
    library(devtools)
    library(dplyr)
    library(CoordinateCleaner)    
    library(MuMIn)    
  #2.1 Download of occurrences from GBIF
    DF_spatial_data_sps<-c()
    DF_n_data_by_sps<-c()
    GBIF_data_by_sps<-c()
    list_sps_in_loop<-c()
    sps_in_loop<-c()
    sps_in_loop$Species<-"Zootoca vivipara"
    vec_sps_name_in_loop<-as.vector(sps_in_loop$Species)
    nam_sps<-as.data.frame(sps_in_loop$Species)
    key<-c(NA)
    num_occur<-c(NA)
    data_by_sps_in_loop<-c()
    name_list<-c(NA)
    occ_search<-c(NA)
    sub_occ_search<-c(NA)
    sub_year<-c(NA)
    dat_clean<-c(NA)
    n_dat_clean<-c(NA)
    i=1
    key<-9702100
    num_occur<-NA
    #Get number of occurrences records:  
    num_occur<-occ_count(taxonKey = key)
    n_dat_clean<-NA
    sum_flags_data<-t(as.data.frame(c(NA,NA,NA,NA,NA,NA,NA,NA)))
    colnames(sum_flags_data)<-c(".val", ".equ",".zer",".cap",".cen",".gbf",".inst",".summary")
    rownames(sum_flags_data)<-c("summary(flags)")
    occ_search_sps_1<-occ_search(taxonKey = key, return='data', limit = 50000, hasCoordinate=TRUE,decimalLatitude='43,46',decimalLongitude='2,5') #,fields=c("species","coordinateUncertaintyInMeters","coordinatePrecision","year","decimalLongitude","decimalLatitude","basisOfRecord","countryCode")
    save(occ_search_sps_1, file="occ_search_sps_1.Rdata")  
    data_by_sps_in_loop<-c(NA,NA,NA)
    sub_year<-subset(occ_search_sps_1,year>1989)
    sub_europe_lat<-subset(sub_year,decimalLongitude>2&decimalLongitude<5)
    sub_europe_lon<-subset(sub_europe_lat,decimalLatitude>43&decimalLatitude<46)
    #We subset by coordinate uncertainity https://terms.tdwg.org/wiki/dwc:coordinateUncertaintyInMeters    
    # or by coordinate precission https://terms.tdwg.org/wiki/dwc:coordinatePrecision
    attach(sub_europe_lon)
    dat_sub_uncert<-subset(sub_europe_lon,coordinateUncertaintyInMeters<564.1)
    detach(sub_europe_lon)
    #Basis of record https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/BasisOfRecord.html
    sub_basis<-subset(dat_sub_uncert, basisOfRecord == "HUMAN_OBSERVATION" |  basisOfRecord == "LIVING_SPECIMEN" |  basisOfRecord == "MACHINE_OBSERVATION" |  basisOfRecord == "OBSERVATION")
    ####  COORDINATE CLEANER  ######
    #We are going to use coordinate cleaner, We have followed the turorial https://ropensci.github.io/CoordinateCleaner/articles/Tutorial_Cleaning_GBIF_data_with_CoordinateCleaner.html
    dat<-sub_basis
    dat <- data.frame(dat)
    #funtion https://ropensci.github.io/CoordinateCleaner/reference/clean_coordinates.html
    flags <- clean_coordinates(x = dat, lon = "decimalLongitude", lat = "decimalLatitude",
      countries = "countryCode", 
      species = "species",
      tests = c("capitals", "centroids", "equal","gbif", "institutions",
      "zeros")) # most test are on by default
    sum_flags_data<-t(as.data.frame(summary(flags)))
    #We take the clean data:
    dat_clean <- dat[flags$.summary,]
    #We are going to create a shapefile of points for 564 uncert:
      n_dat_clean<-NA
      Y564<-dat_clean$decimalLatitude
      X564<-dat_clean$decimalLongitude
      mapa_dat_clean<-as.data.frame(cbind(X564,Y564))
      class(mapa_dat_clean)
      coordinates(mapa_dat_clean)=~X564+Y564
      class(mapa_dat_clean)
      proj4string(mapa_dat_clean) = CRS("+init=epsg:4326")
      #We have created a spatial point object with the projection  WGS84: mapa_sub_coord_uncert_564_df and we are going to define other projection and to project in this new one  
      #We are going to use Europe Albers Equal Area Conic  
      newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
      #We change the projection of the presences
      mapa_dat_clean_df_EAEAC <- spTransform(mapa_dat_clean, newproj)
      mapa_dat_clean_df_EAEAC_df=SpatialPointsDataFrame(mapa_dat_clean_df_EAEAC, data.frame(id=1:length(mapa_dat_clean_df_EAEAC)))
      writeOGR(obj=mapa_dat_clean_df_EAEAC_df, your_wd,layer="Occurrences_sps_1", "ESRI Shapefile",overwrite_layer=TRUE)
      df_mapa_dat_clean_df_EAEAC<-as.data.frame(mapa_dat_clean_df_EAEAC)
      n_dat_clean<-nrow(df_mapa_dat_clean_df_EAEAC)
    #We rasterize the data of presences from GBIF         
      #Raster of reference: 
      new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
      #We defined the projection of our raster:
      projection(new_raster) <- newproj
      #Rasterize the points of species selected
      sps_raster<-rasterize(mapa_dat_clean_df_EAEAC, new_raster, field=1)
      #pr <- projectRaster(m1_crop, new_raster, method='bilinear')
      e2<- extent(-2800000,3100000,-800000,4800000) 
      Vector_sps_1<-extract(sps_raster,e2)
      #We save the raster:
      save(Vector_sps_1, file="Vector_sps_1.Rdata")  
      load("Vector_sps_1.Rdata")
      writeRaster(sps_raster, 'sps_raster.img',overwrite=TRUE) 

    vec_sps_name_in_loop<-as.vector(sps_in_loop$Species)
    nam_sps<-as.data.frame(sps_in_loop$Species)
    #We count the number of presences:
    n_presences_sps<-NA
    no_na_vec<-na.omit(Vector_sps_1)
    n_presences_sps<-length(no_na_vec)
    DF_data_sps<-as.data.frame(cbind(i,vec_sps_name_in_loop,n_presences_sps,sum_flags_data))
    data_by_sps_in_loop<-cbind(key,num_occur,n_dat_clean)
    DF_n_data_by_sps<-rbind(DF_n_data_by_sps, DF_data_sps)
    GBIF_data_by_sps<-rbind(GBIF_data_by_sps,data_by_sps_in_loop)
    list_sps_in_loop<-rbind(list_sps_in_loop,"Species"=vec_sps_name_in_loop)
        
    GBIF_data_by_sps<-as.data.frame(GBIF_data_by_sps)
    df_datos<-cbind(DF_n_data_by_sps,GBIF_data_by_sps)
    df_datos_numeric<-df_datos

    df_datos_numeric$i <- as.numeric(as.character(df_datos_numeric$i))
    df_datos_numeric$n_presences_sps <- as.numeric(as.character(df_datos_numeric$n_presences_sps))
    df_datos_numeric$.val <- as.numeric(as.character(df_datos_numeric$.val))
    df_datos_numeric$.equ <- as.numeric(as.character(df_datos_numeric$.equ))
    df_datos_numeric$.zer <- as.numeric(as.character(df_datos_numeric$.zer))
    df_datos_numeric$.cap  <- as.numeric(as.character(df_datos_numeric$.cap))
    df_datos_numeric$.cen <- as.numeric(as.character(df_datos_numeric$.cen))
    df_datos_numeric$.gbf <- as.numeric(as.character(df_datos_numeric$.gbf))
    df_datos_numeric$.inst <- as.numeric(as.character(df_datos_numeric$.inst))
    df_datos_numeric$.summary <- as.numeric(as.character(df_datos_numeric$.summary))
    df_datos_numeric$key <- as.numeric(as.character(df_datos_numeric$key))
    df_datos_numeric$num_occur <- as.numeric(as.character(df_datos_numeric$num_occur))
    df_datos_numeric$n_dat_clean <- as.numeric(as.character(df_datos_numeric$n_dat_clean))
    rownames(df_datos_numeric)<-df_datos_numeric$i
    df_datos_numeric<-df_datos_numeric[,c(1,2,12,13,4,5,6,7,8,9,10,11,14,3)]
    colnames(df_datos_numeric)<-c("i","Species","GBIF_key","Total_occur","Coord_clean.val","Coord_clean.equ","Coord_clean.zer","Coord_clean.cap","Coord_clean.cen","Coord_clean.gbf","Coord_clean.inst","Coord_clean.summary","Filter_occur","N_pixels_presence")
    #save(df_datos_numeric,file="df_datos_numeric.RData")
      
  #2.2 We create a polygon and a raster with the study area/area of sampling from GBIF
    x_coord <- c(2,  2,  5, 5)
    y_coord <- c(43, 46, 46, 43)
    xym <- cbind(x_coord, y_coord)
    p = Polygon(xym)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    proj4string(sps) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    sps_EAEA <- spTransform(sps, newproj)
    sps_EAEA_df=SpatialPolygonsDataFrame(sps_EAEA, data.frame(id=1:length(sps_EAEA)))
    writeOGR(obj=sps_EAEA_df, your_wd ,layer="Study_area", "ESRI Shapefile",overwrite_layer=TRUE)
    
      #Raster of reference: 
      new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
      #We defined the projection of our raster:
      projection(new_raster) <- newproj
      #Rasterize the data
      gbif_download_area_raster<-rasterize(sps_EAEA, new_raster)
      e2<- extent(-2800000,3100000,-800000,4800000) 
      vec_gbif_download_area_raster<-extract(gbif_download_area_raster,e2)
      #We exported the raster to ascii format. The we can convert it in Arcgis 9.4 using converts form Ascii to raster
      writeRaster(gbif_download_area_raster, 'gbif_download_area_raster.img',overwrite=TRUE) 
      #save(vec_gbif_download_area_raster,file="vec_gbif_download_area_raster.RData")
    
  #2.3 We add to the data of the presences from GBIF with the study area, the data of ID of pixels and the environmetnal variables
    #Data of coordinates and Ids of pixels of the raster
    load("data_matrix.RData")# 
    #Scenario current (Current Climate variables)
    load("Scenario_current.RData")# 
    #We add all data  
    DF_LC_current_data_presence_habitat<-as.data.frame(cbind(Vector_sps_1,vec_gbif_download_area_raster,data_matrix,Scenario_current))
    #save(DF_LC_current_data_presence_habitat, file="DF_LC_current_data_presence_habitat.RData")
    #We select the study area
    data_study_area <- subset(DF_LC_current_data_presence_habitat, vec_gbif_download_area_raster>0)
    data_study_area$Vector_sps_1[is.na(data_study_area$Vector_sps_1)]<-0
    #We select terrestrial areas:
    data_terrestrial_area <- na.omit(data_study_area)
    #save(data_terrestrial_area, file="data_terrestrial_area.RData")
    #We select the presences
    data_presences <- subset(data_terrestrial_area, Vector_sps_1==1)
    n_presences<-nrow(data_presences)#55 obs
    #We select the absences
    data_pre_absences <- subset(data_terrestrial_area, Vector_sps_1==0)
    #Random selection:
    #Same number as presences for univariable models or average modelling
    data_absences<-data_pre_absences[sample(nrow(data_pre_absences), n_presences), ]
    #20 times the number of presences for biomod2
    data_absences_biomod<-data_pre_absences[sample(nrow(data_pre_absences), 20*n_presences), ]
    #We merge the presences and pseudoabsences calculated:
    presences_absences<-rbind(data_presences,data_absences)
    presences_absences_biomod<-rbind(data_presences,data_absences_biomod)
    save(presences_absences, file="presences_absences.RData")
    save(presences_absences_biomod, file="presences_absences_biomod.RData")
  
  #2.4 Calculation of univariable models of habitat for variable selection  
    var_names_hab<-colnames(presences_absences)
    names_var_hab<-var_names_hab[6:24]

        DF_AIC_univar_models_habitat<-c()
        for (v in 1:19){
          var_in_loop_habitat<-names_var_hab[v]
          write_model_prefix_habitat<-c("(Vector_sps_1)~")
          fixed_term_habitat<-var_in_loop_habitat
          write_model_var_in_loop_habitat_1<-paste(var_in_loop_habitat,sep="")# lineal effect
          write_model_var_in_loop_habitat_2<-paste("I(",var_in_loop_habitat,"^2",")",sep="")# quadratic effect
          write_model_var_in_loop_habitat_3<-paste(var_in_loop_habitat," + ","I(",var_in_loop_habitat,"^2",")",sep="")# lineal effect and quadratic effect,
      
          write_model_1_habitat<-paste(write_model_prefix_habitat,write_model_var_in_loop_habitat_1)# lineal effect
          write_model_2_habitat<-paste(write_model_prefix_habitat,write_model_var_in_loop_habitat_2)# quadratic effect
          write_model_3_habitat<-paste(write_model_prefix_habitat,write_model_var_in_loop_habitat_3)# lineal effect and quadratic effect
          
          from_mod_1_habitat<-as.formula(write_model_1_habitat)
          from_mod_2_habitat<-as.formula(write_model_2_habitat)
          from_mod_3_habitat<-as.formula(write_model_3_habitat)
          
          mod_1_habitat<- glm(formula = from_mod_1_habitat, data = presences_absences, family = "binomial", na.action = "na.fail") 
          mod_2_habitat<- glm(formula = from_mod_2_habitat, data = presences_absences, family = "binomial", na.action = "na.fail") 
          mod_3_habitat<- glm(formula = from_mod_3_habitat, data = presences_absences, family = "binomial", na.action = "na.fail")  
          
          AICc_mod_1_habitat<-AICc(mod_1_habitat)
          AICc_mod_2_habitat<-AICc(mod_2_habitat)
          AICc_mod_3_habitat<-AICc(mod_3_habitat)
          AICc_mod_habitat_123<-cbind(AICc_mod_1_habitat,AICc_mod_2_habitat,AICc_mod_3_habitat)
          AICc_mod_habitat_123_min<-min(AICc_mod_habitat_123[1,])
          AICc_in_loop_habitat<-cbind(AICc_mod_habitat_123_min,AICc_mod_habitat_123)
          DF_AIC_univar_models_habitat<-rbind(DF_AIC_univar_models_habitat,AICc_in_loop_habitat)
        }
        colnames(DF_AIC_univar_models_habitat)<-c("AICc_min","AICc_lin","AICc_quadr","AICc_quadr_lin")
        rownames(DF_AIC_univar_models_habitat)<-names_var_hab
        DF_AIC_univar_models_habitat<-as.data.frame(DF_AIC_univar_models_habitat)
        sorted_variables_habitat<-sort(DF_AIC_univar_models_habitat[,1])
        save(sorted_variables_habitat, file="sorted_variables_habitat.RData")# Supplementary Table A6
        save(DF_AIC_univar_models_habitat, file="DF_AIC_univar_models_habitat.RData")

 #2.5 Calculation of a null model with only the intercept
    mod_null_habitat <- glm(Vector_sps_1 ~  1, data=presences_absences)
    summary(mod_null_habitat)
    AICc(mod_null_habitat)
 #2.6 We compare univariable models with the null model
    DF_AIC_univar_models_habitat$Delta_AIC_from_null<-DF_AIC_univar_models_habitat$AICc_min-AICc(mod_null_habitat)
    save(DF_AIC_univar_models_habitat, file="DF_AIC_univar_models_habitat.RData")
    write.csv(DF_AIC_univar_models_habitat, file="DF_AIC_univar_models_habitat.csv")
    #We select the best three variables based on AICc
    DF_AIC_univar_models_habitat_best_models <- subset(DF_AIC_univar_models_habitat, Delta_AIC_from_null<(-3))
    DF_AIC_univar_models_habitat_best_models
    save(DF_AIC_univar_models_habitat_best_models, file="DF_AIC_univar_models_habitat_best_models.RData")
        
  #2.7 Check the correlation of variables 
        cor_habitat<-cor(presences_absences[,c(6:24)])#
        write.csv(cor_habitat, file="cor_habitat.csv")# Supplementary Table A7
    
  #2.8 We are going to  use biomod2 to fit habitat models
    #We use R version 4.0.2  
    rm(list=ls()) 
    ls()
    your_wd<-"WRITE HERE YOUR WORKING DIRECTORY FOLDER"          
    setwd(your_wd)
    library(biomod2) 
    library(slam)
    library(MASS)
    library(usdm)
    library(rgdal) 
    library(biomod2) 
    library(caret)
    library(MuMIn)
    library(raster)
    load("presences_absences_biomod.Rdata")    
    table(presences_absences_biomod$Vector_sps_1)
    #We convert the "0" of the column of presences of the species to NA
    sp_1<-presences_absences_biomod[,c("Vector_sps_1")]
    sp_1[sp_1=="0"] <- NA
    sum_presences<-sum(sp_1,na.rm =T)                
  
      #2.8.1 We formating the data for biomod before was Clim_3","Clim_5"
        format_data<-BIOMOD_FormatingData(resp.var=sp_1,expl.var= presences_absences_biomod[,c("Clim_2","Clim_5")],resp.xy=presences_absences_biomod[,c("x","y")],resp.name = "resp_name",
                                            PA.nb.rep=8,PA.nb.absences=sum_presences, PA.strategy="random")       
        save(format_data, file="format_data.RData")
      #2.8.2 Modeling options 
          #We customize the set of parameters and options of the models that are going to be build   
          ptm_Model_opt <- proc.time()  
          Model_opt<-BIOMOD_ModelingOptions(
            GLM=list(type="quadratic",interaction.level=0),
            GBM=list(n.trees=3000),
            RF = list(ntree = 750),
          )
          save(Model_opt, file="Model_opt.RData")
  
      #2.8.3 Building models 
          Sps_models<-BIOMOD_Modeling(
            data=format_data,
            models=c("GLM","GBM","RF"),
            models.options = Model_opt,
            NbRunEval = 2,
            DataSplit = 70,
            VarImport = 4,#VarImport = 0,
            models.eval.meth = c("TSS"),
            do.full.models=F,
            modeling.id = "ex2")
          save(Sps_models, file="Sps_models.RData")
  
      #2.8.4 Ensemble Modeling
          Sps_essemble_models<-BIOMOD_EnsembleModeling(modeling.output=Sps_models,
                                                       em.by = "all", 
                                                       eval.metric = 'TSS',
                                                       eval.metric.quality.threshold = 0.3,
                                                       models.eval.meth = c('TSS'),
                                                       prob.mean = FALSE,
                                                       prob.cv = FALSE,
                                                       committee.averaging = TRUE,
                                                       prob.mean.weight = FALSE,
                                                       VarImport = 4)
          save(Sps_essemble_models, file="Sps_essemble_models.RData")
        
      #2.8.5 Ensemble Modeling Evaluation and variable importance  
        #Evaluation
            eval<-get_evaluations(Sps_essemble_models,as.data.frame=T)
            save(eval, file="eval.RData")
            write.csv(eval, file='eval.csv')
          #Variable importance
            var_imp_df<- get_variables_importance(Sps_essemble_models,as.data.frame=T) 
            var_imp_df_mean<-apply(var_imp_df,1,mean)
            save (var_imp_df, file='var_imp_df.RData')
            write.csv(var_imp_df, file='var_imp_df.csv')
            save (var_imp_df_mean, file='var_imp_df_mean.RData')
            write.csv(var_imp_df_mean, file='var_imp_df_mean.csv')

      #2.8.6 Ensemble Forecasting of Current Scenario 
        #Data of coordinates and Ids of pixels of the raster
        load("data_matrix.RData")# 
        #Scenario current (Current Climate variables)
        load("Scenario_current.RData")# 
        #This is the vector with the study area (The area of the population where we extrapolate models)
        load("vec_gbif_download_area_raster.RData")
        DF_prediction_0<-as.data.frame(cbind(data_matrix,Scenario_current,vec_gbif_download_area_raster))
        #We select the study area
        DF_prediction_1 <- subset(DF_prediction_0, vec_gbif_download_area_raster>0)
        DF_prediction_2<-na.omit(DF_prediction_1)

          load("Sps_essemble_models.RData")
          variables_names<-Sps_essemble_models@ expl.var.names
          data_pre_for_current_projection<-DF_prediction_2[,variables_names]

          Sps_esemble_curr_proj<-BIOMOD_EnsembleForecasting(EM.output=Sps_essemble_models,
            new.env=data_pre_for_current_projection,
            proj.name="Current",
            selected.models = 'all', 
            binary.meth = "TSS", 
            do.stack=FALSE,
            build.clamping.mask = T)
          save(Sps_esemble_curr_proj, file="Sps_esemble_curr_proj.RData")
          habitat_biomod<-(Sps_esemble_curr_proj@ proj@ val[,1])           
          save (habitat_biomod, file='habitat_biomod.RData')
          merge_pred_biomod_pixels<-cbind(DF_prediction_2,habitat_biomod)
          merge_Europe_habitat_biomod<-merge(data_matrix,merge_pred_biomod_pixels, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          habitat_biomod_Europe<-merge_Europe_habitat_biomod$habitat_biomod
          save (habitat_biomod_Europe, file='habitat_biomod_Europe.RData')
      
    ##2.9 We see the results of the models
      install.packages("ggplot2")
      library(gbm) 
      library(ggplot2)

      #2.9.1 We see the respose curve for each variable (THANKS for the code on this section to Prof. Antoine Guisan and Dr. Olivier Broennimann)   
        #This code was facilated by Prof. Antoine Guisan and Dr. Olivier Broennimann
    
   #2.9.2 We create a map with the predcition of habitat (Figure 2b)
      #We use R version 3.5.3
      rm(list=ls()) 
      ls()
      your_wd<-"WRITE HERE YOUR WORKING DIRECTORY FOLDER"          
      setwd(your_wd)
      library(rgdal) 
      library(raster)
      load("habitat_biomod_Europe.RData")
          habitat_biomod_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
          #We are going to use Europe Albers Equal Area Conic  
          newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
          #We defined the projection of our raster:
          projection(habitat_biomod_raster) <- newproj
          values(habitat_biomod_raster)<-habitat_biomod_Europe 
          writeRaster(habitat_biomod_raster, "habitat_biomod_raster.img", overwrite=TRUE) 

          
###################################################################################################################################################
#3 CORRELATION TELOMERE LENGHT WITH HABITAT
###################################################################################################################################################
  #Load the data of the telomere prediction scaled
  load("vec_pred_over_Europe.RData")
  load("vec_gbif_download_area_raster.RData")
  
  merge_habitat_telomer_predict<-as.data.frame(cbind(habitat_biomod_Europe,vec_pred_over_Europe,vec_gbif_download_area_raster))
  colnames(merge_habitat_telomer_predict)<-c("Habitat_biomod","Telomere_lenght","Study_area")
  save(merge_habitat_telomer_predict, file="merge_habitat_telomer_predict.RData")
  #We exclude area outside study area and the sea of the study area
  terrestrial_study_Area_predictions<-na.omit(merge_habitat_telomer_predict)
  save(terrestrial_study_Area_predictions, file="terrestrial_study_Area_predictions.RData")

  # Table A9. Linear model for predicted values of telomere length index (ranging from 0 to 1) and using scaled values of estimated values of habitat suitability as predictors (ranging from 0 to 1).
  terrestrial_study_Area_predictions$Habitat_biomod_scaled<- terrestrial_study_Area_predictions$Habitat_biomod/1000
  lm_telomere_lineal_biomod<-lm(Telomere_lenght ~ Habitat_biomod_scaled, data = terrestrial_study_Area_predictions)
  summary(lm_telomere_lineal_biomod)#Table A9
  
 #3.1 Validation of telomere length as presence/absences for species distribution
  #We use R version 4.0.2  
  rm(list=ls()) 
  ls()
  your_wd<-"WRITE HERE YOUR WORKING DIRECTORY FOLDER"          
  setwd(your_wd)
  install.packages("data.table")
  install.packages("caret")
      library(data.table) 
      library(caret) 
  load("merge_habitat_telomer_predict.RData")
  DF0 <- read.csv(file="DupoueData.csv",head=TRUE,sep=",")

  quantile_cut<-quantile(DF0$telomere..kb., probs = c(0, 0.10, 0.95, 1)) # quartile 10% = 11.887
  scaled_quantile<-(quantile_cut[2]-min(DF0$telomere..kb.)+0.0001)/(max(DF0$telomere..kb.)-min(DF0$telomere..kb.)+0.001)
  merge_habitat_telomer_predict$Telomere_lenght_binomial = ifelse(merge_habitat_telomer_predict$Telomere_lenght < scaled_quantile,0,1)
  merge_habitat_telomer_predict2<-merge_habitat_telomer_predict
  save(merge_habitat_telomer_predict2, file="merge_habitat_telomer_predict2.RData")
  Telomere_lenght_binomial<-merge_habitat_telomer_predict2$Telomere_lenght_binomial

  load("presences_absences.RData")
  load("data_matrix.RData")# 
  merge_habitat_telomer_predict3<-cbind(data_matrix,merge_habitat_telomer_predict2)
  merge_habitat_telomer_predict4<-merge(merge_habitat_telomer_predict3,presences_absences, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
  save (merge_habitat_telomer_predict4, file='merge_habitat_telomer_predict4.RData')

  merge_habitat_telomer_predict5<-na.omit(merge_habitat_telomer_predict4)
  table_observed_predicted<- table(merge_habitat_telomer_predict5$Vector_sps_1,merge_habitat_telomer_predict5$Telomere_lenght_binomial)
  conf_matrix <-table(merge_habitat_telomer_predict5$Vector_sps_1,merge_habitat_telomer_predict5$Telomere_lenght_binomial)
  statistics<-confusionMatrix(conf_matrix)
  statistics#Values of sensivity, specifity and accuracy
  




        