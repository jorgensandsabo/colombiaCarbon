# Experimenting with getting data analyzed via GEE from the R console and bringing the result directly back into R.

# IMPORTANT: In order for this to work at all, it is necessary to first create a GEE account, download Anaconda, and 
# run GEE_setup.sh from a bash shell. I make no guarantees about whether GEE_setup.sh works properly on Windows.
# I experimented with running the entire bash script from R via system2(), but ultimately we'll have to manually input 
# our GEE authentication keys either way, so it's quite cumbersome to do entirely via R, and we don't get the benefit of
# complete automation anyway, since we'll need to manually input the GEE key no matter what.

# Code adapted from https://philippgaertner.github.io/2019/12/earth-engine-rstudio-reticulate/

library(reticulate)

##### Set up GEE session #####
#Sys.setenv(PATH = paste(c("C:/ProgramData/Anaconda3/bin", Sys.getenv("PATH")), collapse = .Platform$path.sep))

use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication
#np <- import("numpy")       # Import Numpy        not sure this is necessary; it's a holdover from an example I borrowed from
#pd <- import("pandas")      # Import Pandas       ditto the above

##### Extract ALOS elevation from a given point (can be looped over our points if necessary) #####
#point_of_interest <- ee$Geometry$Point(-72.37150102, -0.631402982)
point_of_interest <- ee$Geometry$Point(-72.3715, -0.631403)
ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
ALOS_elev <- ALOS$select('AVE_DSM')
poi_elev <- ALOS$reduceRegions(point_of_interest, ee$Reducer$mean(), 30)$getInfo()
poi_elev$features[[1]]$properties$AVE_DSM

##### Extract ALOS elevations from many points #####
# (e.g. from shapefile or dataframe); prefer solution that doesn't require staging shapefile through Google Cloud Storage
# Write shapefile for point PSP6 as test
#library(sf)
#PSP6 <- st_sf(st_sfc(st_point(c(-72.37150102, -0.631402982), "PSP6")), crs = 4326)
#st_write(PSP6, dsn = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Points/PSP6/PSP6.shp", 
#         layer = "/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/Points/PSP6/PSP6.shp", 
#         driver = "ESRI Shapefile")

ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
ALOS_elev <- ALOS$select('AVE_DSM')
pts <- read.csv("C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD\\Data\\elevations\\CO_sampling_points_metafile.csv")
pts <- pts[-which(pts$long > 0),]

# Featurecollection of point geometries
geompts <- list()
for(i in 1:nrow(pts)){geompts[[i]] <- ee$Geometry$Point(c(pts$long[i],pts$lat[i]))}
geompts <- ee$FeatureCollection(c(unlist(geompts)))

# Extract ALOS elevations for all points - combine into dataframe
pts_elev <- ALOS$reduceRegions(geompts, ee$Reducer$mean(),scale=30.922080775909325)$getInfo()
ALOSelev <- vector()
for(i in 1:length(pts_elev$features)){
  ALOSelev[i] <- pts_elev$features[[i]]$properties$AVE_DSM
}
ALOSelev <- as.data.frame(cbind(point_id=as.character(pts$point_id),ALOSelev))

##### Extract mean elevation from 100-m buffer around point (or many points) #####


##### Do some fancier calculation on raster at a point (or many points) in GEE; import result to R #####
# For example, get slope and aspect data from ALOS
ALOS_slope <- ee$Terrain$slope(ALOS)

# Extract ALOS elevations for all points - combine into dataframe
pts_slope <- ALOS_slope$reduceRegions(geompts, ee$Reducer$mean(), 30)$getInfo()
ALOSslope <- vector()
for(i in 1:length(pts_slope$features)){
  ALOSslope[i] <- pts_slope$features[[i]]$properties$mean
}
ALOSslope <- cbind(point_id=as.character(pts$point_id),ALOSslope)

##### Extract actual elevation raster over some defined polygon (e.g. a 10 km buffer around a point), import to R as raster object #####
# (This might be more relevant if we want to automatically pull in GFC rasters at some future date to derive complicated fragstat-style measures that we can't figure out how to code directly in GEE)
