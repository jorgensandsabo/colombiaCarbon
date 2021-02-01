### Code to add spatial raster data (e.g. climate / topography) to plot data

# Use machine identifier to automatically set correct working path
if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}

setwd(dir.path)

### Point files
Spatialdata <- read.csv("Data\\vegetation\\Plots_prep.csv", stringsAsFactors = FALSE)

### SpatialPoints object
plotcoords <- sp::SpatialPoints(data.frame(cbind(Spatialdata$long,Spatialdata$lat)),proj4string=sp::CRS("+init=epsg:4326"))

#################################
######## WORLDCLIM DATA #########
#################################
# Worldclim data + PET (CGIAR-CSI DATABASE) + AET 
nrasters <- c("01","02","03","04","05","06","07","08","09","10","11","12")
precrasters <- list();tminrasters <- list();tmaxrasters <- list();tmeanrasters <- list();petrasters <- list();aetrasters <- list();biovariables <- list();tavgs <- list();pets <- list();precs <- list();aets <- list();sradrasters <- list();srads <- list()

for(i in 1:length(nrasters)){precrasters[[i]] <- raster::crop(raster::raster(paste("GIS\\WorldClim_PrecMonth\\wc2.0_30s_prec_",nrasters[i],".tif",sep="")),raster::extent(c(-81.3,-66.2,-4.53,13.03)))
  tminrasters[[i]] <- raster::crop(raster::raster(paste("GIS\\WorldClim_TempMinMonth\\wc2.0_30s_tmin_",nrasters[i],".tif",sep="")),raster::extent(c(-81.3,-66.2,-4.53,13.03)))
  tmaxrasters[[i]] <- raster::crop(raster::raster(paste("GIS\\WorldClim_TempMaxMonth\\wc2.0_30s_tmax_",nrasters[i],".tif",sep="")),raster::extent(c(-81.3,-66.2,-4.53,13.03)))
  tmeanrasters[[i]] <- raster::crop(raster::raster(paste("GIS\\WorldClim_TempAvgMonth\\wc2.0_30s_tavg_",nrasters[i],".tif",sep="")),raster::extent(c(-81.3,-66.2,-4.53,13.03)))
  petrasters[[i]] <- raster::crop(raster::raster(paste("GIS\\PET_CGIAR-CSI\\et0_month\\et0_",nrasters[i],".tif",sep="")),raster::extent(c(-81.3,-66.2,-4.53,13.03)))
  aetrasters[[i]] <- raster::crop(raster::raster(paste("GIS\\SoilWaterBalance\\aet_monthly\\aet_",i,".tif",sep="")),raster::extent(c(-81.3,-66.2,-4.53,13.03)))
  sradrasters[[i]] <- raster::crop(raster::raster(paste("GIS\\WorldClim_SolarRadiation\\wc2.0_30s_srad_",nrasters[i],".tif",sep="")),raster::extent(c(-81.3,-66.2,-4.53,13.03)))
}

for(i in 1:length(plotcoords)){
  print(paste("starting plot ",i," of ",length(plotcoords),sep=""))
  prec <- vector();tmin <- vector();tmax <- vector();tavg <- vector();pet <- vector();aet <- vector();srad <- vector()
  for(j in 1:length(nrasters)){
    prec[j] <- raster::extract(precrasters[[j]],sp::spTransform(plotcoords[i],raster::projection(precrasters[[j]])))
    tmin[j] <- raster::extract(tminrasters[[j]],sp::spTransform(plotcoords[i],raster::projection(tminrasters[[j]])))
    tmax[j] <- raster::extract(tmaxrasters[[j]],sp::spTransform(plotcoords[i],raster::projection(tmaxrasters[[j]])))
    tavg[j] <- raster::extract(tmeanrasters[[j]],sp::spTransform(plotcoords[i],raster::projection(tmeanrasters[[j]])))
    pet[j] <- raster::extract(petrasters[[j]],sp::spTransform(plotcoords[i],raster::projection(petrasters[[j]])))
    aet[j] <- raster::extract(aetrasters[[j]],sp::spTransform(plotcoords[i],raster::projection(aetrasters[[j]])))
    srad[j] <- raster::extract(sradrasters[[j]],sp::spTransform(plotcoords[i],raster::projection(sradrasters[[j]])))
  }
  
  biovariables[[i]] <- dismo::biovars(prec,tmin,tmax)
  tavgs[[i]] <- tavg;pets[[i]] <- pet;precs[[i]] <- prec;aets[[i]] <- aet;srads[[i]] <- srad
}

# Add WorldClim data to Spatialdata
#Spatialdata$TotPrec <- NA;Spatialdata$MeanTemp <- NA;Spatialdata$PrecVar <- NA;Spatialdata$TempVar <- NA
Spatialdata$TotPrec <- NA
Spatialdata$PrecVar <- NA
Spatialdata$PrecMax <- NA
Spatialdata$PrecMin <- NA
Spatialdata$PrecWetQ <- NA
Spatialdata$PrecDryQ <- NA
Spatialdata$PrecWarmQ <- NA
Spatialdata$MeanTemp <- NA
Spatialdata$TempVar <- NA
Spatialdata$TempMax <- NA
Spatialdata$TempMin <- NA
Spatialdata$TempRange <- NA
Spatialdata$MeanDiurR <- NA
Spatialdata$Isotherm <- NA
Spatialdata$TempWetQ <- NA
Spatialdata$TempDryQ <- NA
Spatialdata$TempWarmQ <- NA
Spatialdata$TempColdQ <- NA
Spatialdata$MeanArid <- NA
Spatialdata$AridMax <- NA
Spatialdata$AridMin <- NA
Spatialdata$AridVar <- NA
Spatialdata$DryMonths <- NA
Spatialdata$CWDmax <- NA
Spatialdata$SolRad <- NA
for(i in 1:nrow(Spatialdata)){
  WD <- precs[[i]]-pets[[i]]
  ordWD <- c(seq(first(which(WD == max(WD))), 12, 1), seq(1, first(which(WD == max(WD)))-1, 1))
  WD <- WD[order(ordWD)]
  WD[1] <- 0
  for(j in 2:length(WD)){
    ifelse(WD[j-1] + WD[j] < 0, WD[j] <- WD[j-1] + WD[j], WD[j] <- 0)
  }
  
  Spatialdata$TotPrec[i] <- biovariables[[i]][12]
  Spatialdata$PrecVar[i] <- biovariables[[i]][15]
  Spatialdata$PrecMax[i] <- biovariables[[i]][13]
  Spatialdata$PrecMin[i] <- biovariables[[i]][14]
  Spatialdata$PrecWetQ[i] <- biovariables[[i]][16]
  Spatialdata$PrecDryQ[i] <- biovariables[[i]][17]
  Spatialdata$PrecWarmQ[i] <- biovariables[[i]][18]
  Spatialdata$MeanTemp[i] <- biovariables[[i]][1]
  Spatialdata$TempVar[i] <- biovariables[[i]][4]
  Spatialdata$TempMax[i] <- biovariables[[i]][5]
  Spatialdata$TempMin[i] <- biovariables[[i]][6]
  Spatialdata$TempRange[i] <- biovariables[[i]][7]
  Spatialdata$MeanDiurR[i] <- biovariables[[i]][2]
  Spatialdata$Isotherm[i] <- biovariables[[i]][3]
  Spatialdata$TempWetQ[i] <- biovariables[[i]][8]
  Spatialdata$TempDryQ[i] <- biovariables[[i]][9]
  Spatialdata$TempWarmQ[i] <- biovariables[[i]][10]
  Spatialdata$TempColdQ[i] <- biovariables[[i]][11]
  Spatialdata$MeanArid[i] <- mean(precs[[i]])/mean(pets[[i]])
  Spatialdata$AridMax[i] <- min(precs[[i]]/pets[[i]])
  Spatialdata$AridMin[i] <- max(precs[[i]]/pets[[i]])
  Spatialdata$AridVar[i] <- var(precs[[i]]/pets[[i]])/Spatialdata$MeanArid[i]
  Spatialdata$DryMonths[i] <- length(which(precs[[2]]/pets[[2]] < 1))
  Spatialdata$CWDmax[i] <- min(WD)
  Spatialdata$SolRad[i] <- mean(srads[[i]])
}

########################
###### ELEVATIONS ######
########################
library(reticulate)
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()            

# Raster
ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
ALOS_elev <- ALOS$select('AVE_DSM')

# Points
geompts <- list()
for(i in 1:nrow(Spatialdata)){geompts[[i]] <- ee$Geometry$Point(c(Spatialdata$long[i],Spatialdata$lat[i]))}
geompts <- ee$FeatureCollection(c(unlist(geompts)))

# Elevations
pts_elev <- ALOS$reduceRegions(geompts, ee$Reducer$mean(),scale=30.922080775909325)$getInfo()
ALOSelev <- vector()
for(i in 1:length(pts_elev$features)){
  ALOSelev[i] <- pts_elev$features[[i]]$properties$AVE_DSM
}
spatialdata <- cbind.data.frame(Spatialdata,ALOSelev)

# Slope
ALOS_slope <- ee$Terrain$slope(ALOS)
pts_slope <- ALOS_slope$reduceRegions(geompts, ee$Reducer$mean(), 30)$getInfo()
ALOSslope <- vector()
for(i in 1:length(pts_slope$features)){
  ALOSslope[i] <- pts_slope$features[[i]]$properties$mean
}
spatialdata <- cbind.data.frame(spatialdata,ALOSslope)

# Aspect
ALOS_aspect <- ee$Terrain$aspect(ALOS)
pts_aspect <- ALOS_aspect$reduceRegions(geompts, ee$Reducer$mean(), 30)$getInfo()
ALOSaspect <- vector()
for(i in 1:length(pts_aspect$features)){
  ALOSaspect[i] <- pts_aspect$features[[i]]$properties$mean
}
spatialdata <- cbind.data.frame(spatialdata,ALOSaspect)

# Continuous heat-insolation load index
ALOS_chili <- ee$Image("CSP/ERGo/1_0/Global/ALOS_CHILI")
pts_chili <- ALOS_chili$reduceRegions(geompts, ee$Reducer$mean(), 30)$getInfo()
ALOSchili <- vector()
for(i in 1:length(pts_chili$features)){
  ALOSchili[i] <- pts_chili$features[[i]]$properties$mean
}
spatialdata <- cbind.data.frame(spatialdata,ALOSchili)

########################
###### TPI300  #########
########################
tpi300 <- raster::raster("GIS\\ALOS_DEM30\\tpi300.tif")

tpi300s <- vector()
for(i in 1:length(plotcoords)){
  tpi300s[i] <- raster::extract(tpi300,sp::spTransform(plotcoords[i],raster::projection(tpi300)))
}

tpi300_stdi <- (((tpi300s-mean(tpi300s))/sd(tpi300s))*100)+0.5
tpi300_sp <- vector()
for(i in 1:length(tpi300_stdi)){
   if(tpi300_stdi[i] > 100){tpi300_sp[i] <- "ridge"}
   else if(tpi300_stdi[i] <= 100 & tpi300_stdi[i] > 50){tpi300_sp[i] <- "upper slope"}
   else if(tpi300_stdi[i] <= 50 & tpi300_stdi[i] > -50 & spatialdata$ALOSslope[i] > 10){tpi300_sp[i] <- "middle slope"}
   else if(tpi300_stdi[i] <= 50 & tpi300_stdi[i] > -50 & spatialdata$ALOSslope[i] <= 10){tpi300_sp[i] <- "flats slope"}
   else if(tpi300_stdi[i] < 50 & tpi300_stdi[i] >= -100){tpi300_sp[i] <- "lower slope"}
   else if(tpi300_stdi[i] < -100){tpi300_sp[i] <- "valley"}
}

tpi300c <- tpi300_sp
for(i in 1:length(tpi300_sp)){
  if(tpi300c[i] %in% c("middle slope", "upper slope", "lower slope")){tpi300c[i] <- "slope"}  
  if(tpi300c[i] == "flats slope"){tpi300c[i] <- "flat"} 
}

spatialdata <- cbind(spatialdata,tpi300s)
#spatialdata <- cbind(spatialdata,tpi300_sp)
spatialdata <- cbind(spatialdata, tpi300c)

###########################
###### SOIL DATA  ######### Currently only dominant soil type ("soildata" not in use)
###########################
# soildata <- read.csv("GIS\\SoilDatabase1_21\\HWSD_DATA.csv", stringsAsFactors = F)
# 
# soiltype <- read.csv("GIS\\SoilDatabase1_21\\HWSD_SMU.csv", stringsAsFactors = F)
# soilraster <- raster::raster("GIS\\SoilDatabase1_21\\hwsd.tif")
# MU_GLOBAL <- raster::extract(soilraster, sp::spTransform(plotcoords, raster::projection(soilraster)))
# spatialdata <- cbind(spatialdata, MU_GLOBAL)
# spatialdata <- dplyr::left_join(spatialdata,soiltype[,c("MU_GLOBAL","SU_SYMBOL")])
# spatialdata$MU_GLOBAL <- NULL
# names(spatialdata)[which(names(spatialdata) == "SU_SYMBOL")] <- "SoilType"

################
#### SAVE ######
################
write.csv(spatialdata,"Data\\vegetation\\Spatialdata.csv",row.names=F)








#### TESTS
ALOShli <- raster::raster("ALOShli.tif")

ALOShlis <- vector()
for(i in 1:length(plotcoords)){
  ALOShlis[i] <- raster::extract(ALOShli,sp::spTransform(plotcoords[i],raster::projection(ALOShli)))
}
spatialdata <- cbind(spatialdata,ALOShlis)

ALOSsra <- raster::raster("ALOSsra.tif")
ALOSsras <- vector()
for(i in 1:length(plotcoords)){
  ALOSsras[i] <- raster::extract(ALOSsra,sp::spTransform(plotcoords[i],raster::projection(ALOSsra)))
}
spatialdata <- cbind(spatialdata,ALOSsras)

ALOS_mtpi <- ee$Image("CSP/ERGo/1_0/Global/ALOS_mTPI")
pts_mtpi <- ALOS_mtpi$reduceRegions(geompts, ee$Reducer$mean(), 30)$getInfo()
ALOSmtpi <- vector()
for(i in 1:length(pts_mtpi$features)){
  ALOSmtpi[i] <- pts_mtpi$features[[i]]$properties$mean
}
spatialdata <- cbind.data.frame(spatialdata,ALOSmtpi)















## REMOVE

########################
###### ELEVATIONS ######
########################
library(reticulate)
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()            

# Raster
ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
ALOS_elev <- ALOS$select('AVE_DSM')

# Points
geompts <- list()
for(i in 1:nrow(Spatialdata)){geompts[[i]] <- ee$Geometry$Point(c(Spatialdata$long[i],Spatialdata$lat[i]))}
geompts <- ee$FeatureCollection(c(unlist(geompts)))

# Elevations
pts_elev <- ALOS$reduceRegions(geompts, ee$Reducer$mean(),scale=30.922080775909325)$getInfo()
ALOSelev <- vector()
for(i in 1:length(pts_elev$features)){
  ALOSelev[i] <- pts_elev$features[[i]]$properties$AVE_DSM
}
spatialdata <- cbind.data.frame(Spatialdata,ALOSelev)







######### 
dataset <- ee$ImageCollection('COPERNICUS/S2_SR')$filterDate('2015-01-01', '2020-01-30')$filter(ee$Filter$lt('CLOUDY_PIXEL_PERCENTAGE',5))

pts_elev <- dataset$reduceRegions(geompts, ee$Reducer$mean(),scale=30.922080775909325)$getInfo()
ALOSelev <- vector()
for(i in 1:length(pts_elev$features)){
  ALOSelev[i] <- pts_elev$features[[i]]$properties$AVE_DSM
}
spatialdata <- cbind.data.frame(Spatialdata,ALOSelev)


