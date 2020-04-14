
#### Function to download ALOS (AW3D30_v2.2) elevation and cloud-reduced median landsat imagery for an ROI.
#### Accepts single coordinates or vectors of coordinates
#### Square/circle around single point (length = buffer.width)
#### Rectangular boundary around vectors of points (min distance to edge = buffer.width)
dl_alos_ls <- function(roi_name = "MtEverest", roi_long = 86.922623, roi_lat = 27.986065,
                       buffer.width = 1000, max.error = 1, pixscale = 30.922080775909325,
                       cloud_percentile = 30, cloud_score = 1, max.pixels = 1e20){

  # Get geometry of roi
  coords <- lapply(1:length(roi_lat),function(x)cbind(roi_long,roi_lat)[x,])
  
  if(length(coords) == 1){
    roi <- ee$Geometry$Point(coords[[1]])$buffer(buffer.width,max.error)$bounds()
  }else if(length(coords) > 1){
    roi <- ee$Geometry$MultiPoint(coords)$buffer(buffer.width,max.error)$bounds()
  }
  
  # Get ALOS elevation raster from Earth Engine
  ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
  ALOS_elev <- ALOS$select("AVE_DSM")

  # Get cloud-free landsat raster from Earth Engine
  LSraw <- ee$ImageCollection("LANDSAT/LC08/C01/T1")
  LScf_allbands <- ee$Algorithms$Landsat$simpleComposite(
    collection = LSraw,
    percentile = cloud_percentile,
    cloudScoreRange = cloud_score);
  LScf <- LScf_allbands$select(list("B4","B3","B2"))
  
  # Reduce rasters to roi
  latlng <- ee$Image$pixelLonLat()$addBands(ALOS_elev)$addBands(LScf)
  latlng <- latlng$reduceRegion(reducer = ee$Reducer$toList(),
                                geometry = roi,
                                maxPixels = max.pixels,
                                scale=pixscale)

  # Convert to arrays
  lats <- np$array((ee$Array(latlng$get("latitude"))$getInfo()))
  lngs <- np$array((ee$Array(latlng$get("longitude"))$getInfo()))
  ras_vals <- np$array((ee$Array(latlng$get("AVE_DSM"))$getInfo()))

  red <- np$array((ee$Array(latlng$get("B4"))$getInfo()))
  green <- np$array((ee$Array(latlng$get("B3"))$getInfo()))
  blue <- np$array((ee$Array(latlng$get("B2"))$getInfo()))

  # Convert to elevation raster
  ras <- data.frame(x = lngs, y = lats, ras = ras_vals)
  rasterElev <- raster::rasterFromXYZ(ras,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

  # Convert landsat to rasterBrick
  rasterRGB <- raster::brick(list(raster::rasterFromXYZ(data.frame(x = lngs, y = lats, red = red),crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"),
                                  raster::rasterFromXYZ(data.frame(x = lngs, y = lats, green = green),crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"),
                                  raster::rasterFromXYZ(data.frame(x = lngs, y = lats, blue = blue),crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")))
  
  # Return elevation and landsat rasters
  return(list(rasterElev = rasterElev, rasterRGB = rasterRGB,
              info = list(roi_name = roi_name, 
                          roi_long = roi_long,
                          roi_lat = roi_lat,
                          pixscale = pixscale)))
}



########## Stretch functions, stolen from RStoolbox
RGBstretch <- function (x, method = "lin", quantiles = c(0.02,0.98), band = NULL) {
  if(method == "lin"){
    if(length(quantiles) == 1) quantiles <- c(0,1) + c(quantiles, -quantiles)/100
    v <- quantile(x, quantiles, na.rm = TRUE)
    if(diff(v)==0) {
      ## sometimes lowr and upr quantile can be identical, which breaks the color calculation --> enforce a minimal distance by adding ~0
      v[2] <- v[2] + 1e-9
    }
    temp <-  (x - v[1])/(v[2] - v[1])
    temp[temp < 0] <- 0
    temp[temp > 1] <- 1
    return(temp)
  } 
  
  if(method == "hist"){
    ecdfun <- ecdf(x)
    return(ecdfun(x))
  } 
  
  if(method == "log"){
    x <- log(x + 1)
    x <- x - min(x,na.rm=TRUE)
    return(x / max(x,na.rm=TRUE))         
  }
  
  if(method == "sqrt"){
    x <- sqrt(x)
    x <- x - min(x,na.rm=TRUE)
    return(x /max(x,na.rm=TRUE))
  }
}


########## Function to prepare the rasters for rayshader
########## Change lightness and stretch values of RGB raster
########## Reproject to a squared metric crs
########## Convert to matrices
prep_elev_RGB <- function(rasterElev, rasterRGB, roi_lat, roi_long,
                          lab_longs = NULL, lab_lats = NULL, lab_names = NULL,
                          projection = "albersequalarea", pixscale = 30.922080775909325,
                          lightness = 1, stretch.method = "lin", stretch.quantiles = c(0.02,0.98)){
  
  # Change visuals of RGB raster (light + stretch)
  raster::values(rasterRGB) <- raster::values(rasterRGB)^(1/lightness)
  
  if(stretch.method %in% c("lin","hist","log","sqrt")){
    for(i in 1:3){
      raster::values(rasterRGB[[i]]) <- RGBstretch(raster::values(rasterRGB[[i]]), method = stretch.method, quantiles = stretch.quantiles)
    }
  }
  
  # Choose projection
  if(projection == "albersequalarea"){
    proj <- paste("+proj=aea +lat_1=",roi_lat[1]-5," +lat_2=",roi_lat[1]+5,sep="")
  }
  
  # Labels
  if(!(is.null(lab_longs) | is.null(lab_lats) | is.null(lab_names))){
    lpts <- sp::SpatialPoints(cbind(lab_longs,lab_lats),proj4string=sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  }else{
    lpts <- "no labels provided"}
  
  # Reproject
  if(projection %in% c("albersequalarea")){
    rasterElev <- raster::projectRaster(rasterElev,crs=proj,res=pixscale)
    rasterRGB <- raster::projectRaster(rasterRGB,crs=proj,res=pixscale)
    
    if(!(is.null(lab_longs) | is.null(lab_lats) | is.null(lab_names))){
      lpts_reproj <- sp::spTransform(lpts,proj)
      lpts <- data.frame(lpts_reproj,row.names = lab_names)
    }
  }
  
  
  # Convert to matrices
  elevvals_m = rayshader::raster_to_matrix(rasterElev$ras)
  
  RGBvals_m <- array(NA,dim=c(dim(rasterRGB[[i]])[1],dim(rasterRGB[[i]])[2],3))
  for(i in 1:3){
    RGBvals_m[,,i] <- matrix(raster::values(rasterRGB[[i]]),nrow=dim(rasterRGB[[i]])[1],ncol=dim(rasterRGB[[i]])[2],byrow=TRUE)
  }
  
  # Return objects
  return(list(elevvals_m = elevvals_m, RGBvals_m = RGBvals_m, rasterElev = rasterElev, rasterRGB = rasterRGB, lab_pts = lpts,
              info = list(pixscale = pixscale)))
}









### Load required package and GEE
library(reticulate)
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication
np <- import("numpy")       # Import Numpy        needed for converting gee raster to R raster object
pd <- import("pandas")      # Import Pandas       ditto the above

# Read coordinate file
dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD\\Data_raw\\elevations"
pts <- read.csv(paste(dir.path,"\\CO_sampling_points_metafile.csv",sep=""),stringsAsFactors = F)
pts <- pts[-which(pts$long > 0),]   # Needed due to a cluster with error in the current metafile

# Examples
  # Point
    region_name <- "Hytte" ; lat <- 58.698105 ; long <- 6.820840
    region_name <- "Preikestolen" ; lat <- 58.986239 ; long <- 6.224986
    region_name <- "Moi" ; lat <- 58.456153 ; long <- 6.552710
    region_name <- "SMF1" ; lat <- pts[which(pts$point_id==point_name),]$lat ; long <- pts[which(pts$point_id==point_name),]$long
  # Multiple points
    region_name <- c("CC","RA","PP") ; lat <- pts[which(pts$site %in% region_name),]$lat ; long <- pts[which(pts$site %in% region_name),]$long ; lab_longs <- c(-73.21665,-73.44127,-73.37779) ; lab_lats <- c(6.824891,6.969141,6.848298) ; lab_names <- c("Chicamocha Canyon","Pauxi Pauxi","Reinita Azul")
    region_name <- c("CC","RA") ; lat <- pts[which(pts$site %in% region_name),]$lat ; long <- pts[which(pts$site %in% region_name),]$long ; lab_longs <- c(-73.21665,-73.37779) ; lab_lats <- c(6.824891,6.848298) ; lab_names <- c("Chicamocha Canyon","Reinita Azul")
    region_name <- "IG" ; lat <- pts[which(pts$site==region_name),]$lat ; long <- pts[which(pts$site==region_name),]$long
  # Boundary rectangle
    region_name <- "Jotunheimen" ; lat <- c(61.75,61.75,61.5,61.5) ; long <- c(8.13,8.5,8.13,8.5) 
    
# Run functions
alos_ls <- dl_alos_ls(roi_name = region_name, roi_lat = lat, roi_long = long, buffer.width = 5000)
alos_ls_prep <- prep_elev_RGB(rasterElev = alos_ls$rasterElev, rasterRGB = alos_ls$rasterRGB,
                              roi_lat = alos_ls$info$roi_lat, roi_long = alos_ls$info$roi_long,
                              lab_longs = lab_longs, lab_lats = lab_lats, lab_names = lab_names,
                              pixscale = alos_ls$info$pixscale, stretch.method = "lin",
                              lightness = 5)

# Rayshader
library(rayshader)

zscale <- alos_ls_prep$info$pixscale
heightmap <- alos_ls_prep$elevvals_m
overlay <- alos_ls_prep$RGBvals_m
lab_pts <- alos_ls_prep$lab_pts

rgl::clear3d()
heightmap %>%
  sphere_shade(texture = "desert",zscale=zscale) %>%
  add_overlay(overlay, alphacolor = NULL, alphalayer = 0.99) %>%
  add_shadow(ray_shade(heightmap, zscale = 3), 0.5) %>%
  add_shadow(ambient_shade(heightmap),0.5) %>%
  plot_3d(heightmap, zscale = zscale, fov = 30, lineantialias = TRUE, theta = 45, phi = 35, zoom = 0.8)

for(i in 1:nrow(lab_pts)){
  render_label(heightmap = heightmap, 
               x = dim(heightmap)[1]*(lab_pts[i,1] - raster::extent(alos_ls_prep$rasterElev)[1])/(raster::extent(alos_ls_prep$rasterElev)[2] - raster::extent(alos_ls_prep$rasterElev)[1]),
               y = dim(heightmap)[2]*(lab_pts[i,2] - raster::extent(alos_ls_prep$rasterElev)[3])/(raster::extent(alos_ls_prep$rasterElev)[4] - raster::extent(alos_ls_prep$rasterElev)[3]), 
               z = 1500, linewidth = 1, zscale = zscale, text = rownames(lab_pts)[i],relativez = T,
               linecolor="red",textcolor="red",freetype=F)
}

#render_depth(focus = 0.5, focallength = 10, clear = TRUE,filename="test.png")
render_snapshot(filename = "test.png")
#render_highquality(filename = "test.png")
rgl::clear3d()


