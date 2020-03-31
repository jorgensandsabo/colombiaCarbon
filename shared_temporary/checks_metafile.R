### Temporary scrips - check the raw project metafile 
setwd("C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD")

meta <- read.csv("Data\\elevations\\CO_sampling_points_metafile.csv",na.strings=c("NA","na"))

### A few checks
unique(meta$site)
unique(meta$cluster)
unique(meta$region)
unique(meta$mountain_slope)
unique(meta$ecosystem)
unique(meta$biome)
unique(meta$state)
unique(meta$county)
unique(meta$locality)
which(is.na(meta$lat));which(is.na(meta$long))

#############################
## CHECK ALOS ELEVATION DATA

### AW3D30 tiles were downloaded from https://www.eorc.jaxa.jp/ALOS/en/aw3d30/index.htm
### Tiles were merged into a composite covering Colombia using gdal_merge tool in QGIS

### Plot coordinates
metacoords <- sp::SpatialPoints(data.frame(cbind(meta$long,meta$lat)),proj4string=sp::CRS("+init=epsg:4326"))
raster::plot(metacoords)

### Extract AW3D30 point elevations
rasterALOS <- raster::raster("GIS\\ALOS_DEM30\\ALOS_DEM_merged.tif")
ALOSelev2 <- raster::extract(rasterALOS,sp::spTransform(metacoords,raster::projection(rasterALOS)))
meta <- cbind(meta,ALOSelev2)

### Plot elevation differences
plot(meta$ALOSelev,meta$elev_ALOS30m);abline(0,1)
plot(meta$ALOSelev,meta$elev_strn90m);abline(0,1)
plot(meta$ALOSelev,meta$elev_gps);abline(0,1)

### Add difference columns, check manually
meta$diffALOS_ALOS <- abs(meta$ALOSelev - meta$elev_ALOS30m)
meta$diffALOS_SRTM <- abs(meta$ALOSelev - meta$elev_strn90m)
meta$diffALOS_GPS <- abs(meta$ALOSelev - meta$elev_gps)
meta$diffSRTM_GPS <- abs(meta$elev_strn90m - meta$elev_gps)

###############################
## CHECK AGAINST PREVIOUS DATA
plots <- read.csv("Data_clean\\Clean_raw\\Plotdata_clean.csv")
spatial <- read.csv("Data_clean\\Clean_raw\\Spatialdata.csv")

### Same lat/long?
interm <- meta
names(interm)[1] <- "SiteCode"
comb <- dplyr::left_join(spatial,interm,by="SiteCode")
rm(interm)

plot(comb$lat.x,comb$lat.y);abline(0,1)
plot(comb$long.x,comb$long.y);abline(0,1)
diffs <- comb[which((comb$lat.x != comb$lat.y) & (comb$long.x != comb$long.y)),]
diffs$latdiffs <- diffs$lat.x - diffs$lat.y
diffs$longdiffs <- diffs$long.x - diffs$long.y

### Missing coords?
test <- plots[which(!plots$SiteCode %in% meta$point_id),]
test2 <- meta[which(!meta$point_id %in% plots$SiteCode),]

