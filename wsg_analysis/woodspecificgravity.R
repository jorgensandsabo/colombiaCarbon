#################
################# Analysis of forest structure data
#################
if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}

setwd(dir.path)

library(dplyr)
library(brms)
library(ggplot2)

###############
## Read data ##
###############
Treedata <- read.csv("Data\\vegetation\\Trees_prep.csv", header = T, stringsAsFactors = F)
Spatial <- read.csv("Data\\vegetation\\Spatialdata.csv", header = T, stringsAsFactors = F)

Treedata <- dplyr::filter(Treedata, SiteCode %in% Spatial$SiteCode)      # Remove trees from plots not in spatial data (lacking latitude/longitude data)
Spatial[which(Spatial$SiteCode == "01-Sep"),]$SiteCode <- "SEP1"
Spatial[which(Spatial$SiteCode == "02-Sep"),]$SiteCode <- "SEP2"
Spatial[which(Spatial$SiteCode == "03-Sep"),]$SiteCode <- "SEP3"

# Filter forest data only
treesF <- Treedata
treesF <- dplyr::left_join(treesF, Spatial, by = "SiteCode")
treesF <- dplyr::filter(treesF, HabitatP == "Forest")

treesF <- dplyr::filter(treesF, !is.na(lat))              
treesF <- dplyr::filter(treesF, is.na(FAge) | FAge > 20)  # Exclude secondary young forest from Choco data
treesF <- dplyr::filter(treesF,AreaCode != "CC")          # Dry forest site
treesF <- dplyr::filter(treesF,AreaCode != "TU")          # High elevation polylepis forest

# Remove plots with no cores
treesF <- dplyr::filter(treesF, !SiteCode %in% c("C15_P2","PSF11","PSF12","TAF1","TAF2","TSF5"))

# Remove palms and ferns
treesF <- dplyr::filter(treesF, !Species %in% c("Palm","Fern"))

# Remove species with less than two individuals or no cores
coredSp <- treesF %>% mutate(cored = WSG/WSG)
coredSp <- coredSp %>% group_by(Species) %>% summarise(total = n(), cored = sum(cored,na.rm=T))
coredSp <- dplyr::filter(coredSp, cored > 1)
coredSp <- dplyr::filter(coredSp, total > cored)

## Add "speciesid" column (1 = identified, 0 = not identified)
treesF[which(!treesF$Species %in% coredSp$Species),]$Species <- NA
treesF[which(is.na(treesF$Species)),]$Species <- "nospec"
treesF$speciesid <- ifelse(treesF$Species == "nospec", 0, 1)

# Calculate volume and plot volume-weights
treesF <- dplyr::mutate(treesF, vol = exp(2.789 - 1.414 * log(DBH_used) + 1.178  * (log(DBH_used)^2) - 0.118 * (log(DBH_used)^3)))
treesF$vol <- treesF$vol / 1000     # From dm^3 to m^3 for readable effect sizes. Back-transform before AGB calculation
treesF <- left_join(treesF, treesF %>% group_by(SiteCode) %>% summarise(totVOL = sum(vol)))
treesF$VOLw <- treesF$vol / treesF$totVOL

### Remove plots - only for repeating the analysis without weird plots
# Plots dominated by oak
#left_join(treesF %>% group_by(SiteCode) %>% summarise("ntrees" = n(), "AreaCode" = first(AreaCode)),treesF %>% count(SiteCode,Species)) %>% filter(Species == "Quercus Humboldtii") %>% mutate(perc = 100*n/ntrees)
#left_join(treesF %>% group_by(SiteCode) %>% summarise("ntrees" = n(), "AreaCode" = first(AreaCode)),treesF %>% count(SiteCode,Species)) %>% filter(Species == "Quercus Humboldtii") %>% mutate(perc = 100*n/ntrees) %>% group_by(AreaCode) %>% summarise("ntrees" = sum(ntrees), "noak" = sum(n)) %>% mutate(perc = 100 * noak / ntrees) %>% print()
#treesF <- filter(treesF, !AreaCode %in% c("AT","GA","SE","SO","TA","TS","VC","VI"))

# Secondary: western andes, ALF, BRF, PPF, RBF
#treesF <- filter(treesF, !AreaCode %in% c("AL","BR","PP","RB"))
#treesF <- dplyr::filter(treesF, is.na(HabitatChoco) | HabitatChoco == "P")

# Choco
#treesF <- filter(treesF, Dataset != "Chocó")

## Cluster points by distance
plotdata <- as.data.frame(treesF %>% group_by(SiteCode) %>% summarise(AreaCode = first(AreaCode), Cluster = first(Cluster), Region = first(region), lat = first(lat), long = first(long)))
coords <- sp::SpatialPoints(cbind(plotdata$long, plotdata$lat))
distmat <- sp::spDists(cbind(plotdata$long, plotdata$lat),longlat = TRUE)
hc <- hclust(as.dist(distmat), "ward.D2")

  # Find smallest within-area distance variance that ensures larger within-area maximum distance than smallest between-area minimum distance
  for(i in 1:100){
    treeheight <- i
    if(fpc::cluster.stats(as.dist(distmat), clustering = cutree(hc, h = i))$max.diameter - fpc::cluster.stats(as.dist(distmat), clustering = cutree(hc, h = i))$min.separation < 0){break}}
  # Plot dendrogram / map
  plot(hc, labels = plotdata$SiteCode) ; abline(h = seq(0,50,10), col = "red") ; abline(h = c(15,25), col = "blue")
  raster::plot(coords, col = cutree(hc, h = treeheight), pch = 18)

# Add area numbers
plotdata$AreaNum <- cutree(hc, h = treeheight)      

data.frame(plotdata %>% group_by(AreaNum) %>% summarise(unique(AreaCode))) # Merges AL/SM, AT/IG/GA, PP/RA, TA/TS/VI/VC
fpc::cluster.stats(as.dist(distmat), clustering = cutree(hc, h = treeheight))$max.diameter
fpc::cluster.stats(as.dist(distmat), clustering = cutree(hc, h = treeheight))$min.separation

## Slight changes to regions
plotdata$Region[which(plotdata$Region == "oriental na")] <- "oriental western slope"
plotdata$Region[which(plotdata$AreaCode == "CQ")] <- "oriental eastern slope"
plotdata$RegionNum <- as.numeric(as.factor(plotdata$Region))

### Add grouping levels as numbers
treesF$SiteNum <- as.numeric(as.factor(treesF$SiteCode))
treesF$ClusNum <- as.numeric(as.factor(treesF$Cluster))
treesF <- left_join(treesF, plotdata[,c("SiteCode","AreaNum")])
treesF <- left_join(treesF, plotdata[,c("SiteCode","RegionNum")])
treesF <- left_join(treesF, plotdata[,c("SiteCode","Region")])

### Colour and region names for plots and legends
unique(treesF$Region)
treesF$Region[which(treesF$Region == "oriental eastern slope")] <- "Oriental - east"
treesF$Region[which(treesF$Region == "oriental western slope")] <- "Oriental - west"
treesF$Region[which(treesF$Region == "central-eastern")] <- "Central - south"
treesF$Region[which(treesF$Region == "santa_marta")] <- "Santa Marta"
treesF$Region[which(treesF$Region == "central")] <- "Central - north"
treesF$Region[which(treesF$Region == "amazonia")] <- "Amazon"
treesF$Region[which(treesF$Region == "occidental")] <- "Occidental"

unique(treesF[order(treesF$RegionNum),]$Region)
regcol <- c("red","green","blue","yellow","grey","black","white")
regcol <- c("#267234","#656150","#b5b09e","#84a955","#84a955","#6f5a67","#44437e")
regcol <- c("#189a39","#55928f","#b5b09e","#84a955","#decd60","#6f5a67","#44437e")


### Plot dataset
plotsF <- treesF %>% group_by(SiteCode) %>% summarise(SiteNum = first(SiteNum), ClusNum = first(ClusNum), AreaNum = first(AreaNum), RegionNum = first(RegionNum),
                                                      ntree = n(), ncore = sum(WSG/WSG,na.rm=T), vol = sum(vol), Region = first(Region),
                                                      nQuercus = sum(Species == "Quercus Humboldtii"), nSpec = sum(!Species %in% c("Quercus Humboldtii","nospec")))
plotsF <- as.data.frame(left_join(plotsF,Spatial))
plotsF$region <- NULL
plotsF <- data.frame(plotsF, "ElevS" = plotsF$ALOSelev / 1000, "TotPrecS" = plotsF$TotPrec/1000, "TempVarS" = scale(plotsF$TempVar),
                     "PrecVarS" = scale(plotsF$PrecVar), "SlopeS" = plotsF$ALOSslope/10)
plotsF$volha <- (plotsF$vol / plotsF$Size) * 10000

### Remove temporary objects
rm(coords, coredSp, distmat, hc, plotdata, Spatial, Treedata, i, treeheight)

###########################################
## Predictive model for individual trees ## 
###########################################
moddat <- list("Nobs" = length(treesF$WSG),                       # Number of observations
               "WSG" = treesF$WSG, "vol" = treesF$vol,            # WSG data
               "Species" = as.numeric(as.factor(treesF$Species)), # Vector of species
               "nspec" = length(unique(treesF$Species)),          # Number of species
               "speciesid" = treesF$speciesid,                    # 0/1 - species determined or not
               "SiteNum" = treesF$SiteNum,
               "ClusNumS" = data.frame(treesF %>% group_by(SiteNum) %>% summarise(first(ClusNum)))[,2],
               "AreaNum" = treesF$AreaNum,
               "AreaNumC" = data.frame(treesF %>% group_by(ClusNum) %>% summarise(first(AreaNum)))[,2],
               "nsites" = length(unique(treesF$SiteNum)),
               "ncluster" = length(unique(treesF$ClusNum)),
               "narea" = length(unique(treesF$AreaNum)))

modinits <- function(){list(alpha = rnorm(1,0.5,0.1), muspec = rnorm(1,0,0.001), sdspec = rnorm(1,0.036,0.01), 
                            sdsite = rnorm(1,0.03,0.01), sdcluster = rnorm(1,0.017,0.01), sdarea = rnorm(1,0.05,0.01))}

n.chains <- 4
n.burnin <- 10000
n.iter <- n.burnin + 50000
n.thin <- 50
n.cores <- n.chains

start.time <- Sys.time()
WSGtreemod <- jagsUI::jags(model.file = "JAGS//WSG_trees.txt", data = moddat, inits = modinits,
                           parameters.to.save = c("alpha","var1","var0","b_vol","muspec","varspec","vararea","varsite","varcluster","loglik","WSG.pred","res","b_area","b_cluster","b_site","mu"),
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                           parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res","mu"))
time.used <- Sys.time()-start.time

# saveRDS(WSGtreemod, "Output\\WSG\\WSGtreemod.RDS")

WSGdraws <- WSGtreemod$sims.list$WSG.pred
WSGdraws[, which(!is.na(treesF$WSG))] <- matrix(data = rep(treesF$WSG[which(!is.na(treesF$WSG))], nrow(WSGdraws)), nrow = nrow(WSGdraws), ncol = length(treesF$WSG[which(!is.na(treesF$WSG))]), byrow = TRUE)
WSGdraws <- as.data.frame(t(WSGdraws) * treesF$VOLw)
WSGdraws <- data.frame(WSGdraws) %>% mutate(SiteNum = treesF$SiteNum) %>% reshape2::melt(id.vars = "SiteNum") %>% reshape2::dcast(SiteNum ~ variable, value.var="value", fun.aggregate = sum)

# saveRDS(WSGdraws, "Output\\WSG\\WSGdraws.RDS")

#############################
## Plot-average WSG models ##
#############################
## Extract posterior draws of volume-weighted plot averages
# WSGdraws <- readRDS("Output\\WSG\\WSGdraws.RDS")
WSGdraws_short <- WSGdraws[,seq(1, ncol(WSGdraws), 40)]

## Run model
n.chains1 <- 4                ; n.chains2 <- 4
n.burnin1 <- 1000             ; n.burnin2 <- 1000
n.iter1 <- n.burnin1 + 10000  ; n.iter2 <- n.burnin2 + 50000
n.thin1 <- 100                ; n.thin2 <- 500
n.cores1 <- n.chains1         ; n.cores2 <- n.chains2

plotmods <- list("nullmods" = list("null" = list(), "nullclus" = list(), "nullarea" = list()), 
                 "climmods" = list("clim" = list(), "climclus" = list(), "climarea" = list()))
Rhats <- list("nullmods" = list("null" = list(), "nullclus" = list(), "nullarea" = list()), 
                 "climmods" = list("clim" = list(), "climclus" = list(), "climarea" = list()))
predictors <- as.matrix(plotsF[order(plotsF$SiteCode), c("ElevS","TotPrecS","TempVarS","PrecVarS","SlopeS")])

start.time <- Sys.time()
for(i in 2:ncol(WSGdraws_short)){
  loop.time <- Sys.time()
  
  mod <- jagsUI::jags(model.file = "JAGS//WSG_plots_null.txt", inits = NULL,
                      data = list("Nobs" = nrow(WSGdraws_short),
                                  "WSG" = WSGdraws_short[1:nrow(WSGdraws_short), i]),
                      parameters.to.save = c("alpha","var0",
                                             "WSG.pred","loglik","res","mu"),
                      n.chains = n.chains1, n.iter = n.iter1, n.burnin = n.burnin1, n.thin = n.thin1,
                      parallel = T, n.cores = n.cores1, codaOnly = c("WSG.pred","loglik","res","mu"), DIC = T)
  
  Rhats$nullmods$null$alpha <- append(Rhats$nullmods$null$alpha, rstan::Rhat(matrix(mod$sims.list$alpha, ncol = n.chains1)))
  Rhats$nullmods$null$var0 <- append(Rhats$nullmods$null$var0, rstan::Rhat(matrix(mod$sims.list$var0, ncol = n.chains1)))
  Rhats$nullmods$null$mu <- rbind(Rhats$nullmods$null$mu, apply(mod$sims.list$mu, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains1))))
  
  plotmods$nullmods$null$alpha <- append(plotmods$nullmods$null$alpha, mod$sims.list$alpha)
  #plotmods$nullmods$null$var0 <- append(plotmods$nullmods$null$var0, mod$sims.list$var0)
  plotmods$nullmods$null$WSGpred <- rbind(plotmods$nullmods$null$WSGpred, mod$sims.list$WSG.pred)
  #plotmods$nullmods$null$loglik <- rbind(plotmods$nullmods$null$loglik, mod$sims.list$loglik)
  plotmods$nullmods$null$res <- rbind(plotmods$nullmods$null$res, mod$sims.list$res)
  #plotmods$nullmods$null$mu <- rbind(plotmods$nullmods$null$mu, mod$sims.list$mu)
  plotmods$nullmods$null$DIC <- append(plotmods$nullmods$null$DIC, mod$DIC)
  plotmods$nullmods$null$R2 <- append(plotmods$nullmods$null$R2, apply(mod$sims.list$mu, 1, function(x) cor(x, WSGdraws_short[,i])^2))
  plotmods$nullmods$null$RMSE <- append(plotmods$nullmods$null$RMSE, apply(mod$sims.list$res, 1, function(x) sqrt(mean(x^2))))
  
  print(paste("Loop number ", i-1, " of ", ncol(WSGdraws_short)-1, " took ", round(as.numeric(Sys.time() - loop.time, "mins"),2), " minutes", sep = ""))
}
mod1.time <- Sys.time() - start.time ; start.time <- Sys.time()
for(i in 2:ncol(WSGdraws_short)){
  start.time <- Sys.time()
  
  mod <- jagsUI::jags(model.file = "JAGS//WSG_plots_null_spatial.txt", inits = NULL,
                      data = list("Nobs" = nrow(WSGdraws_short),
                                  "WSG" = WSGdraws_short[1:nrow(WSGdraws_short), i],
                                  "SpatialNum" = plotsF[order(plotsF$SiteCode),]$ClusNum,
                                  "nspatial" = length(unique(plotsF$ClusNum))),
                      parameters.to.save = c("alpha","varspatial","var0","b_spatial",
                                             "WSG.pred","loglik","res","mu"),
                      n.chains = n.chains2, n.iter = n.iter2, n.burnin = n.burnin2, n.thin = n.thin2,
                      parallel = T, n.cores = n.cores2, codaOnly = c("WSG.pred","loglik","res","mu","b_spatial"), DIC = T)
  
  Rhats$nullmods$nullclus$alpha <- append(Rhats$nullmods$nullclus$alpha, rstan::Rhat(matrix(mod$sims.list$alpha, ncol = n.chains2)))
  Rhats$nullmods$nullclus$var0 <- append(Rhats$nullmods$nullclus$var0, rstan::Rhat(matrix(mod$sims.list$var0, ncol = n.chains2)))
  Rhats$nullmods$nullclus$mu <- rbind(Rhats$nullmods$nullclus$mu, apply(mod$sims.list$mu, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains2))))
  Rhats$nullmods$nullclus$varspatial <- append(Rhats$nullmods$nullclus$var0, rstan::Rhat(matrix(mod$sims.list$var0, ncol = n.chains2)))
  Rhats$nullmods$nullclus$b_spatial <- rbind(Rhats$nullmods$nullclus$b_spatial, apply(mod$sims.list$b_spatial, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains2))))
  
  plotmods$nullmods$nullclus$alpha <- append(plotmods$nullmods$nullclus$alpha, mod$sims.list$alpha)
  #plotmods$nullmods$nullclus$var0 <- append(plotmods$nullmods$nullclus$var0, mod$sims.list$var0)
  #plotmods$nullmods$nullclus$varspatial <- append(plotmods$nullmods$nullclus$varspatial, mod$sims.list$varspatial)
  #plotmods$nullmods$nullclus$b_spatial <- rbind(plotmods$nullmods$nullclus$b_spatial, mod$sims.list$b_spatial)
  plotmods$nullmods$nullclus$WSGpred <- rbind(plotmods$nullmods$nullclus$WSGpred, mod$sims.list$WSG.pred)
  #plotmods$nullmods$nullclus$loglik <- rbind(plotmods$nullmods$nullclus$loglik, mod$sims.list$loglik)
  plotmods$nullmods$nullclus$res <- rbind(plotmods$nullmods$nullclus$res, mod$sims.list$res)
  #plotmods$nullmods$nullclus$mu <- rbind(plotmods$nullmods$nullclus$mu, mod$sims.list$mu)
  plotmods$nullmods$nullclus$DIC <- append(plotmods$nullmods$nullclus$DIC, mod$DIC)
  plotmods$nullmods$nullclus$R2 <- append(plotmods$nullmods$nullclus$R2, apply(mod$sims.list$mu, 1, function(x) cor(x, WSGdraws_short[,i])^2))
  plotmods$nullmods$nullclus$RMSE <- append(plotmods$nullmods$nullclus$RMSE, apply(mod$sims.list$res, 1, function(x) sqrt(mean(x^2))))
  
  print(paste("Loop number ", i-1, " of ", ncol(WSGdraws_short)-1, " took ", round(as.numeric(Sys.time() - start.time, "mins"),2), " minutes", sep = ""))
}
mod2.time <- Sys.time() - start.time ; start.time <- Sys.time()
for(i in 2:ncol(WSGdraws_short)){
  start.time <- Sys.time()
  
  mod <- jagsUI::jags(model.file = "JAGS//WSG_plots_null_spatial.txt", inits = NULL,
                      data = list("Nobs" = nrow(WSGdraws_short),
                                  "WSG" = WSGdraws_short[1:nrow(WSGdraws_short), i],
                                  "SpatialNum" = plotsF[order(plotsF$SiteCode),]$AreaNum,
                                  "nspatial" = length(unique(plotsF$AreaNum))),
                      parameters.to.save = c("alpha","varspatial","var0","b_spatial",
                                             "WSG.pred","loglik","res","mu"),
                      n.chains = n.chains2, n.iter = n.iter2, n.burnin = n.burnin2, n.thin = n.thin2,
                      parallel = T, n.cores = n.cores2, codaOnly = c("WSG.pred","loglik","res","mu","b_spatial"), DIC = T)
  
  Rhats$nullmods$nullarea$alpha <- append(Rhats$nullmods$nullarea$alpha, rstan::Rhat(matrix(mod$sims.list$alpha, ncol = n.chains2)))
  Rhats$nullmods$nullarea$var0 <- append(Rhats$nullmods$nullarea$var0, rstan::Rhat(matrix(mod$sims.list$var0, ncol = n.chains2)))
  Rhats$nullmods$nullarea$mu <- rbind(Rhats$nullmods$nullarea$mu, apply(mod$sims.list$mu, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains2))))
  Rhats$nullmods$nullarea$varspatial <- append(Rhats$nullmods$nullarea$var0, rstan::Rhat(matrix(mod$sims.list$var0, ncol = n.chains2)))
  Rhats$nullmods$nullarea$b_spatial <- rbind(Rhats$nullmods$nullarea$b_spatial, apply(mod$sims.list$b_spatial, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains2))))
  
  plotmods$nullmods$nullarea$alpha <- append(plotmods$nullmods$nullarea$alpha, mod$sims.list$alpha)
  #plotmods$nullmods$nullarea$var0 <- append(plotmods$nullmods$nullarea$var0, mod$sims.list$var0)
  #plotmods$nullmods$nullarea$varspatial <- append(plotmods$nullmods$nullarea$varspatial, mod$sims.list$varspatial)
  #plotmods$nullmods$nullarea$b_spatial <- rbind(plotmods$nullmods$nullarea$b_spatial, mod$sims.list$b_spatial)
  plotmods$nullmods$nullarea$WSGpred <- rbind(plotmods$nullmods$nullarea$WSGpred, mod$sims.list$WSG.pred)
  #plotmods$nullmods$nullarea$loglik <- rbind(plotmods$nullmods$nullarea$loglik, mod$sims.list$loglik)
  plotmods$nullmods$nullarea$res <- rbind(plotmods$nullmods$nullarea$res, mod$sims.list$res)
  #plotmods$nullmods$nullarea$mu <- rbind(plotmods$nullmods$nullarea$mu, mod$sims.list$mu)
  plotmods$nullmods$nullarea$DIC <- append(plotmods$nullmods$nullarea$DIC, mod$DIC)
  plotmods$nullmods$nullarea$R2 <- append(plotmods$nullmods$nullarea$R2, apply(mod$sims.list$mu, 1, function(x) cor(x, WSGdraws_short[,i])^2))
  plotmods$nullmods$nullarea$RMSE <- append(plotmods$nullmods$nullarea$RMSE, apply(mod$sims.list$res, 1, function(x) sqrt(mean(x^2))))
  
  print(paste("Loop number ", i-1, " of ", ncol(WSGdraws_short)-1, " took ", round(as.numeric(Sys.time() - start.time, "mins"),2), " minutes", sep = ""))
}
mod3.time <- Sys.time() - start.time ; start.time <- Sys.time()
for(i in 2:ncol(WSGdraws_short)){
  start.time <- Sys.time()
  
  mod <- jagsUI::jags(model.file = "JAGS//WSG_plots_cov.txt", inits = NULL,
                      data = list("Nobs" = nrow(WSGdraws_short),
                                  "WSG" = WSGdraws_short[1:nrow(WSGdraws_short), i],
                                  "cov" = predictors,
                                  "ncov" = ncol(predictors)),
                      parameters.to.save = c("alpha","b_cov","var0",
                                             "WSG.pred","loglik","res","mu"),
                      n.chains = n.chains1, n.iter = n.iter1, n.burnin = n.burnin1, n.thin = n.thin1,
                      parallel = T, n.cores = n.cores1, codaOnly = c("WSG.pred","loglik","res","mu","b_cov"), DIC = T)
  
  Rhats$climmods$clim$alpha <- append(Rhats$climmods$clim$alpha, rstan::Rhat(matrix(mod$sims.list$alpha, ncol = n.chains1)))
  Rhats$climmods$clim$var0 <- append(Rhats$climmods$clim$var0, rstan::Rhat(matrix(mod$sims.list$var0, ncol = n.chains1)))
  Rhats$climmods$clim$mu <- rbind(Rhats$climmods$clim$mu, apply(mod$sims.list$mu, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains1))))
  Rhats$climmods$clim$b_cov <- rbind(Rhats$climmods$clim$b_cov, apply(mod$sims.list$b_cov, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains1))))
  
  plotmods$climmods$clim$alpha <- append(plotmods$climmods$clim$alpha, mod$sims.list$alpha)
  #plotmods$climmods$clim$var0 <- append(plotmods$climmods$clim$var0, mod$sims.list$var0)
  plotmods$climmods$clim$b_cov <- rbind(plotmods$climmods$clim$b_cov, mod$sims.list$b_cov)
  plotmods$climmods$clim$WSGpred <- rbind(plotmods$climmods$clim$WSGpred, mod$sims.list$WSG.pred)
  #plotmods$climmods$clim$loglik <- rbind(plotmods$climmods$clim$loglik, mod$sims.list$loglik)
  plotmods$climmods$clim$res <- rbind(plotmods$climmods$clim$res, mod$sims.list$res)
  #plotmods$climmods$clim$mu <- rbind(plotmods$climmods$clim$mu, mod$sims.list$mu)
  plotmods$climmods$clim$DIC <- append(plotmods$climmods$clim$DIC, mod$DIC)
  plotmods$climmods$clim$R2 <- append(plotmods$climmods$clim$R2, apply(mod$sims.list$mu, 1, function(x) cor(x, WSGdraws_short[,i])^2))
  plotmods$climmods$clim$RMSE <- append(plotmods$climmods$clim$RMSE, apply(mod$sims.list$res, 1, function(x) sqrt(mean(x^2))))
  
  print(paste("Loop number ", i-1, " of ", ncol(WSGdraws_short)-1, " took ", round(as.numeric(Sys.time() - start.time, "mins"),2), " minutes", sep = ""))
}
mod4.time <- Sys.time() - start.time ; start.time <- Sys.time()
for(i in 2:ncol(WSGdraws_short)){
  start.time <- Sys.time()
  
  mod <- jagsUI::jags(model.file = "JAGS//WSG_plots_spatial_cov.txt", inits = NULL,
                      data = list("Nobs" = nrow(WSGdraws_short),
                                  "WSG" = WSGdraws_short[1:nrow(WSGdraws_short), i],
                                  "SpatialNum" = plotsF[order(plotsF$SiteCode),]$ClusNum,
                                  "nspatial" = length(unique(plotsF$ClusNum)),
                                  "cov" = predictors,
                                  "ncov" = ncol(predictors)),
                      parameters.to.save = c("alpha","b_cov","varspatial","var0","b_spatial",
                                             "WSG.pred","loglik","res","mu"),
                      n.chains = n.chains2, n.iter = n.iter2, n.burnin = n.burnin2, n.thin = n.thin2,
                      parallel = T, n.cores = n.cores2, codaOnly = c("WSG.pred","loglik","res","mu","b_cov","b_spatial"), DIC = T)
  
  Rhats$climmods$climclus$alpha <- append(Rhats$climmods$climclus$alpha, rstan::Rhat(matrix(mod$sims.list$alpha, ncol = n.chains2)))
  Rhats$climmods$climclus$var0 <- append(Rhats$climmods$climclus$var0, rstan::Rhat(matrix(mod$sims.list$var0, ncol = n.chains2)))
  Rhats$climmods$climclus$mu <- rbind(Rhats$climmods$climclus$mu, apply(mod$sims.list$mu, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains2))))
  Rhats$climmods$climclus$varspatial <- append(Rhats$climmods$climclus$var0, rstan::Rhat(matrix(mod$sims.list$var0, ncol = n.chains2)))
  Rhats$climmods$climclus$b_spatial <- rbind(Rhats$climmods$climclus$b_spatial, apply(mod$sims.list$b_spatial, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains2))))
  Rhats$climmods$climclus$b_cov <- rbind(Rhats$climmods$climclus$b_cov, apply(mod$sims.list$b_cov, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains2))))
  
  plotmods$climmods$climclus$alpha <- append(plotmods$climmods$climclus$alpha, mod$sims.list$alpha)
  #plotmods$climmods$climclus$var0 <- append(plotmods$climmods$climclus$var0, mod$sims.list$var0)
  #plotmods$climmods$climclus$varspatial <- append(plotmods$climmods$climclus$varspatial, mod$sims.list$varspatial)
  #plotmods$climmods$climclus$b_spatial <- rbind(plotmods$climmods$climclus$b_spatial, mod$sims.list$b_spatial)
  plotmods$climmods$climclus$b_cov <- rbind(plotmods$climmods$climclus$b_cov, mod$sims.list$b_cov)
  plotmods$climmods$climclus$WSGpred <- rbind(plotmods$climmods$climclus$WSGpred, mod$sims.list$WSG.pred)
  #plotmods$climmods$climclus$loglik <- rbind(plotmods$climmods$climclus$loglik, mod$sims.list$loglik)
  plotmods$climmods$climclus$res <- rbind(plotmods$climmods$climclus$res, mod$sims.list$res)
  #plotmods$climmods$climclus$mu <- rbind(plotmods$climmods$climclus$mu, mod$sims.list$mu)
  plotmods$climmods$climclus$DIC <- append(plotmods$climmods$climclus$DIC, mod$DIC)
  plotmods$climmods$climclus$R2 <- append(plotmods$climmods$climclus$R2, apply(mod$sims.list$mu, 1, function(x) cor(x, WSGdraws_short[,i])^2))
  plotmods$climmods$climclus$RMSE <- append(plotmods$climmods$climclus$RMSE, apply(mod$sims.list$res, 1, function(x) sqrt(mean(x^2))))
  
  print(paste("Loop number ", i-1, " of ", ncol(WSGdraws_short)-1, " took ", round(as.numeric(Sys.time() - start.time, "mins"),2), " minutes", sep = ""))
}
mod5.time <- Sys.time() - start.time ; start.time <- Sys.time()
for(i in 2:ncol(WSGdraws_short)){
  start.time <- Sys.time()
  
  mod <- jagsUI::jags(model.file = "JAGS//WSG_plots_spatial_cov.txt", inits = NULL,
                      data = list("Nobs" = nrow(WSGdraws_short),
                                  "WSG" = WSGdraws_short[1:nrow(WSGdraws_short), i],
                                  "SpatialNum" = plotsF[order(plotsF$SiteCode),]$AreaNum,
                                  "nspatial" = length(unique(plotsF$AreaNum)),
                                  "cov" = predictors,
                                  "ncov" = ncol(predictors)),
                      parameters.to.save = c("alpha","b_cov","varspatial","var0","b_spatial",
                                             "WSG.pred","loglik","res","mu"),
                      n.chains = n.chains2, n.iter = n.iter2, n.burnin = n.burnin2, n.thin = n.thin2,
                      parallel = T, n.cores = n.cores2, codaOnly = c("WSG.pred","loglik","res","mu","b_cov","b_spatial"), DIC = T)
  
  Rhats$climmods$climarea$alpha <- append(Rhats$climmods$climarea$alpha, rstan::Rhat(matrix(mod$sims.list$alpha, ncol = n.chains2)))
  Rhats$climmods$climarea$var0 <- append(Rhats$climmods$climarea$var0, rstan::Rhat(matrix(mod$sims.list$var0, ncol = n.chains2)))
  Rhats$climmods$climarea$mu <- rbind(Rhats$climmods$climarea$mu, apply(mod$sims.list$mu, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains2))))
  Rhats$climmods$climarea$varspatial <- append(Rhats$climmods$climarea$var0, rstan::Rhat(matrix(mod$sims.list$var0, ncol = n.chains2)))
  Rhats$climmods$climarea$b_spatial <- rbind(Rhats$climmods$climarea$b_spatial, apply(mod$sims.list$b_spatial, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains2))))
  Rhats$climmods$climarea$b_cov <- rbind(Rhats$climmods$climarea$b_cov, apply(mod$sims.list$b_cov, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains2))))
  
  plotmods$climmods$climarea$alpha <- append(plotmods$climmods$climarea$alpha, mod$sims.list$alpha)
  #plotmods$climmods$climarea$var0 <- append(plotmods$climmods$climarea$var0, mod$sims.list$var0)
  #plotmods$climmods$climarea$varspatial <- append(plotmods$climmods$climarea$varspatial, mod$sims.list$varspatial)
  #plotmods$climmods$climarea$b_spatial <- rbind(plotmods$climmods$climarea$b_spatial, mod$sims.list$b_spatial)
  plotmods$climmods$climarea$b_cov <- rbind(plotmods$climmods$climarea$b_cov, mod$sims.list$b_cov)
  plotmods$climmods$climarea$WSGpred <- rbind(plotmods$climmods$climarea$WSGpred, mod$sims.list$WSG.pred)
  #plotmods$climmods$climarea$loglik <- rbind(plotmods$climmods$climarea$loglik, mod$sims.list$loglik)
  plotmods$climmods$climarea$res <- rbind(plotmods$climmods$climarea$res, mod$sims.list$res)
  #plotmods$climmods$climarea$mu <- rbind(plotmods$climmods$climarea$mu, mod$sims.list$mu)
  plotmods$climmods$climarea$DIC <- append(plotmods$climmods$climarea$DIC, mod$DIC)
  plotmods$climmods$climarea$R2 <- append(plotmods$climmods$climarea$R2, apply(mod$sims.list$mu, 1, function(x) cor(x, WSGdraws_short[,i])^2))
  plotmods$climmods$climarea$RMSE <- append(plotmods$climmods$climarea$RMSE, apply(mod$sims.list$res, 1, function(x) sqrt(mean(x^2))))
  
  print(paste("Loop number ", i-1, " of ", ncol(WSGdraws_short)-1, " took ", round(as.numeric(Sys.time() - start.time, "secs"),2), " seconds", sep = ""))
}
mod6.time <- Sys.time() - start.time

# saveRDS(plotmods, "Output\\WSG\\plotmods.RDS")
# saveRDS(Rhats, "Output\\WSG\\plotmods_rhats.RDS")

##############################
## Spatial cross-validation ##
##############################
kfolds <- list("plot" = list(), "cluster" = list(), "area" = list(), "region" = list())

## Extract posterior draws of volume-weighted plot averages
#WSGdraws <- readRDS("Output\\WSG\\WSGdraws.RDS")
WSGdraws_short <- WSGdraws[,seq(1, ncol(WSGdraws), 40)]

n.chains <- 4   
n.burnin <- 1000
n.iter <- n.burnin + 10000
n.thin <- 100
n.cores <- n.chains

# Select spatial level
spat_cv <- "plot" # region, area, cluster, plot

if(spat_cv == "region"){spatlevel <- plotsF$RegionNum} 
if(spat_cv == "area"){spatlevel <- plotsF$AreaNum}
if(spat_cv == "cluster"){spatlevel <- plotsF$ClusNum}
if(spat_cv == "plot"){spatlevel <- plotsF$SiteNum}
nfolds <- length(unique(spatlevel))
randomplots <- sample(1:nrow(plotsF), nrow(plotsF), replace = FALSE, prob = NULL)
sortedplots <- plotsF[order(spatlevel),]$SiteNum
foldlength <- as.vector(table(spatlevel[order(spatlevel)]))
last <- cumsum(foldlength)
first <- last - foldlength + 1
kfoldlist <- list("elpd" = vector(), "res" = vector(), "mu" = vector(), "WSGpred" = vector(), "prederror" = vector(), "sitenum" = vector(), "foldid" = vector(),
                  "rhatmax_alpha" = vector(), "rhatmax_var0" = vector(), "rhatmax_bcov" = vector())
kfold <- list("spatialfold" = kfoldlist, "spatialfoldNull" = kfoldlist, "randomfold" = kfoldlist,"randomfoldNull" = kfoldlist)
predictors <- as.matrix(plotsF[order(plotsF$SiteCode), c("ElevS","TotPrecS","TempVarS","PrecVarS","SlopeS")])

start.time <- Sys.time()
for(k in 1:4){
  ifelse(k %in% c(1,2), plotorder <- sortedplots, plotorder <- randomplots)
  ifelse(k %in% c(1,3), model.file <- "JAGS//WSG_plots_CV.txt", model.file <- "JAGS//WSG_plots_null_CV.txt")
  
  for(i in 1:length(foldlength)){
    pred <- plotorder[first[i]:last[i]]
    train <- plotorder[-which(plotorder %in% pred)]
    
    elpd.iter <- vector() ; res.iter <- vector() ; mu.iter <- vector() ; WSGpred.iter <- vector() ; prederror.iter <- vector()
    rhat_alpha <- vector() ; rhat_var0 <- vector() ; rhat_bcov <- vector()
    for(j in 2:ncol(WSGdraws_short)){
      start.loop.time <- Sys.time()
      
      mod <- jagsUI::jags(model.file = model.file, inits = NULL,
                          data = list("Nobst" = length(WSGdraws_short[train, j]),
                                      "WSGt" = WSGdraws_short[train, j],
                                      "covt" = matrix(as.numeric(predictors[train,]), ncol = ncol(predictors)),
                                      "Nobsp" = length(WSGdraws_short[pred, j]),
                                      "WSGp" = WSGdraws_short[pred, j],
                                      "covp" = matrix(as.numeric(predictors[pred,]), ncol = ncol(predictors)),
                                      "ncov" = ncol(predictors)),
                          parameters.to.save = c("alpha","b_cov","var0","WSG.pred","loglik","res","mup","prederror"),
                          n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                          parallel = T, n.cores = n.cores, codaOnly = c("b_cov","var0","WSG.pred","loglik","res","mup","b_cov","prederror"), DIC = F)

      elpd.iter <- do.call("rbind", list(elpd.iter, as.matrix(mod$sims.list$loglik)))
      res.iter <- do.call("rbind", list(res.iter, as.matrix(mod$sims.list$res^2)))
      mu.iter <- do.call("rbind", list(mu.iter, as.matrix(mod$sims.list$mup)))
      WSGpred.iter <- do.call("rbind", list(WSGpred.iter, as.matrix(mod$sims.list$WSG.pred)))
      prederror.iter <- do.call("rbind", list(prederror.iter, as.matrix(mod$sims.list$prederror)))

      rhat_alpha <- c(rhat_alpha, rstan::Rhat(matrix(mod$sims.list$alpha, ncol = n.chains)))     # NOT WORKING
      rhat_var0 <- c(rhat_var0, rstan::Rhat(matrix(mod$sims.list$var0, ncol = n.chains)))
      if(k %in% c(1,3)){rhat_bcov <- rbind(rhat_bcov, apply(mod$sims.list$b_cov, 2, function(x) rstan::Rhat(matrix(x, ncol = n.chains))))}else{rhat_bcov <- rbind(rhat_bcov, rep(NA,5))}

      print(paste(spat_cv," level: ",ifelse(k %in% c(1,2), "Spatial","Random"),ifelse(k %in% c(1,3), " environmental"," null")," fold ",i," of ", nfolds,sep=""))
      print(paste("Loop number ", j-1, " of ", ncol(WSGdraws_short)-1, " took ", round(as.numeric(Sys.time() - start.loop.time, "secs"),2), " seconds", sep = ""))
    }

    kfold[[k]]$elpd <- cbind(kfold[[k]]$elpd, elpd.iter)
    kfold[[k]]$res <- cbind(kfold[[k]]$res, res.iter)
    kfold[[k]]$mu <- cbind(kfold[[k]]$mu, mu.iter)
    kfold[[k]]$WSGpred <- cbind(kfold[[k]]$WSGpred, WSGpred.iter)
    kfold[[k]]$prederror <- cbind(kfold[[k]]$prederror, prederror.iter)
    kfold[[k]]$sitenum <- c(kfold[[k]]$sitenum, pred)
    kfold[[k]]$foldid <- c(kfold[[k]]$foldid, rep(i, length(pred)))
    
    kfold[[k]]$rhatmax_alpha <- c(kfold[[k]]$rhatmax_alpha, max(rhat_alpha))
    kfold[[k]]$rhatmax_var0 <- c(kfold[[k]]$rhatmax_var0, max(rhat_var0))
    kfold[[k]]$rhatmax_bcov <- rbind(kfold[[k]]$rhatmax_bcov, apply(rhat_bcov,2,max))
  }

  #kfold[[k]][c("elpd","res","mu","WSGpred")] <- lapply(kfold[[k]][c("elpd","res","mu","WSGpred")], function(x) x[,order(kfold[[k]]$sitenum)])
  #kfold[[k]][c("sitenum","foldid")] <- lapply(kfold[[k]][c("sitenum","foldid")], function(x) x[order(kfold[[k]]$sitenum)])
  
  if(spat_cv == "plot"){if(k == 2){break}}
}

if(spat_cv == "region"){kfolds$region <- kfold ; regmod.time <- Sys.time() - start.time}
if(spat_cv == "area"){kfolds$area <- kfold ; areamod.time <- Sys.time() - start.time}
if(spat_cv == "cluster"){kfolds$cluster <- kfold; clusmod.time <- Sys.time() - start.time}
if(spat_cv == "plot"){kfolds$plot <- kfold ; plotmod.time <- Sys.time() - start.time
                      kfolds$plot$randomfold <- NULL
                      kfolds$plot$randomfoldNull <- NULL}

# saveRDS(kfolds, file = "Output\\WSG\\kfolds.RDS")

##################
## Model output ##
##################
#### Individual tree model
#WSGtreemod <- readRDS("Output\\WSG\\WSGtreemod.RDS")

## Model estimate table
R2 <- apply(WSGtreemod$sims.list$mu[,which(!is.na(treesF$WSG))], 1, function(x) cor(x, treesF$WSG[which(!is.na(treesF$WSG))])^2)
RMSE <- apply(WSGtreemod$sims.list$res[,which(!is.na(treesF$WSG))], 1, function(x) sqrt(sum(x^2)/length(x)))

treemod_esttab <- WSGtreemod$summary[c("alpha","b_vol","vararea","varcluster","varsite","varspec","var0","var1"),c(1,3,7)]
rownames(treemod_esttab) <- c("Intercept","Volume (m3)","Area","Cluster","Site","Species","Residual (no species)","Residual (species)")
treemod_esttab <- rbind(treemod_esttab,
                        "RMSE" = c(mean(RMSE), RMSE[order(RMSE)][length(RMSE)*0.025], RMSE[order(RMSE)][length(RMSE)*0.975]),
                        "R2" = c(mean(R2), R2[order(R2)][length(R2)*0.025], R2[order(R2)][length(R2)*0.975]))

#write.csv(round(treemod_esttab,4) ,"Output\\WSG\\treemod_esttab.csv")

## Range of WSG at different scales
MeanPred <- data.frame(WSG = WSGtreemod$mean$WSG.pred, SiteNum = treesF$SiteNum, AreaNum = treesF$AreaNum, ClusNum = treesF$ClusNum)
range(data.frame(MeanPred %>% group_by(AreaNum) %>% summarise(WSGm = mean(WSG)))[,2]) ; max(data.frame(MeanPred %>% group_by(AreaNum) %>% summarise(WSGm = mean(WSG)))[,2]) - min(data.frame(test %>% group_by(AreaNum) %>% summarise(WSGm = mean(WSG)))[,2])
mean(data.frame(MeanPred %>% group_by(AreaNum, ClusNum) %>% summarise(WSGm = mean(WSG)) %>% group_by(AreaNum) %>% summarise(max(WSGm)-min(WSGm)))[,2]) ; max(data.frame(MeanPred %>% group_by(AreaNum, ClusNum) %>% summarise(WSGm = mean(WSG)) %>% group_by(AreaNum) %>% summarise(max(WSGm)-min(WSGm)))[,2])
mean(data.frame(MeanPred %>% group_by(ClusNum, SiteNum) %>% summarise(WSGm = mean(WSG)) %>% group_by(ClusNum) %>% summarise(max(WSGm)-min(WSGm)))[,2]) ; max(data.frame(MeanPred %>% group_by(ClusNum, SiteNum) %>% summarise(WSGm = mean(WSG)) %>% group_by(ClusNum) %>% summarise(max(WSGm)-min(WSGm)))[,2])
                                                  
## Whisker plot - random effects
# Spatial random effects output dataframes, sorted by mean WSG in regions
b_area <- data.frame(AreaNum = 1:length(unique(treesF$AreaNum)), bArea = WSGtreemod$mean$b_area, bArea95l = WSGtreemod$q2.5$b_area, bArea95u = WSGtreemod$q97.5$b_area, RegionNum = data.frame(treesF %>% group_by(AreaNum) %>% summarise(first(RegionNum)))[,2])
  b_area <- left_join(b_area, b_area %>% group_by(RegionNum) %>% summarise(RegOrd = mean(bArea)))
  b_area <- b_area[order(b_area$RegOrd, b_area$bArea),]
b_cluster <- data.frame(ClusNum = 1:length(unique(treesF$Cluster)), bCluster = WSGtreemod$mean$b_cluster, bCluster95l = WSGtreemod$q2.5$b_cluster, bCluster95u = WSGtreemod$q97.5$b_cluster, RegionNum = data.frame(treesF %>% group_by(ClusNum) %>% summarise(first(RegionNum)))[,2])
  b_cluster <- left_join(b_cluster, b_cluster %>% group_by(RegionNum) %>% summarise(RegOrd = mean(bCluster)))
  b_cluster <- b_cluster[order(b_cluster$RegOrd, b_cluster$bCluster),]
b_site <- data.frame(SiteNum = 1:length(unique(treesF$SiteCode)), bSite = WSGtreemod$mean$b_site, bSite95l = WSGtreemod$q2.5$b_site, bSite95u = WSGtreemod$q97.5$b_site, AreaNum = data.frame(treesF %>% group_by(SiteNum) %>% summarise(first(AreaNum)))[,2], ClusNum = data.frame(treesF %>% group_by(SiteNum) %>% summarise(first(ClusNum)))[,2], RegionNum = data.frame(treesF %>% group_by(SiteNum) %>% summarise(first(RegionNum)))[,2])
  b_site <- left_join(b_site, b_site %>% group_by(RegionNum) %>% summarise(RegOrd = mean(bSite)))
  b_site <- b_site[order(b_site$RegOrd, b_site$bSite),]


tiff(file = "Output\\WSG\\plotwhisker.tiff", width = 6000, height = 3000, res = 900)

  layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE), widths=c(0.6,0.2,0.2), heights=c(3,2.5))
  par(mar = c(0.4,0.4,0.1,0.1), oma = c(0,3,0,0))
  
  plot(NA, ylim = c(min(b_cluster$bCluster95l, b_area$bArea95l, b_site$bSite95l), max(b_cluster$bCluster95u, b_area$bArea95u, b_site$bSite95u)), 
       xlim = c(0,nrow(b_site)), xaxt='n', yaxt = 'n', fg = "gray60")
  for(i in 1:nrow(b_site)){segments(x0=i, y0=b_site$bSite95l[i], x1=i, y1=b_site$bSite95u[i])}
  points(b_site$bSite, pch = 19, cex = 0.7, col = regcol[b_site$RegionNum])
  axis(2)
  text(x = 5, y = 0.2, "Plot", cex = 1)
  abline(h = 0)
  
  plot(NA, ylim = c(min(b_cluster$bCluster95l, b_area$bArea95l, b_site$bSite95l), max(b_cluster$bCluster95u, b_area$bArea95u, b_site$bSite95u)), 
       xlim = c(0,nrow(b_cluster)), xaxt='n', yaxt = 'n', fg = "gray60")
  for(i in 1:nrow(b_cluster)){segments(x0=i, y0=b_cluster$bCluster95l[i], x1=i, y1=b_cluster$bCluster95u[i])}
  points(b_cluster$bCluster, pch = 19, cex = 0.7, col = regcol[b_cluster$RegionNum])
  axis(2)
  text(x = 5, y = 0.2, "Cluster", cex = 1)
  abline(h = 0)
  
  plot(NA, ylim = c(min(b_cluster$bCluster95l, b_area$bArea95l, b_site$bSite95l), max(b_cluster$bCluster95u, b_area$bArea95u, b_site$bSite95u)), 
       xlim = c(0,nrow(b_area)), yaxt='n', xaxt='n', fg = "gray60")#, mar = c(5,0.5,0.1,2), )
  for(i in 1:nrow(b_area)){segments(x0=i, y0=b_area$bArea95l[i], x1=i, y1=b_area$bArea95u[i])}
  points(b_area$bArea, pch = 19, cex = 0.7, col = regcol[b_area$RegionNum])
  text(x = 5, y = 0.2, "Area", cex = 1)
  abline(h = 0)
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("center", legend =  unique(treesF[order(treesF$RegionNum),]$Region), pch=16, pt.cex=1.5, cex=1, bty='n',
         col = regcol, title = "Region", title.adj = -0.005)
  #legend(x = 0.1, y = 0.99, legend =  "Region", cex=2, bty='n', y.intersp = 0, x.intersp = 0)
  
  mtext("Spatial effect on WSG", at = 0.5, side = 2, outer = TRUE, line = 2, cex = 0.7)

dev.off()

######################################################################
#### Plot models
#plotmods <- readRDS("Output\\WSG\\plotmods.RDS")
#WSGdraws <- readRDS("Output\\WSG\\WSGdraws.RDS")
#WSGdraws_short <- WSGdraws[,seq(1, ncol(WSGdraws), 40)]

## Summary table
plotmod_out <- rbind(
  do.call("cbind",lapply(plotmods$climmods, function(x) rbind(
    cbind(mean(x$alpha), x$alpha[order(x$alpha)][length(x$alpha) * 0.025], x$alpha[order(x$alpha)][length(x$alpha) * 0.975]),
    #cbind(cbind(apply(plotmods$climmods$climarea$b_cov, 2, mean)),cbind(apply(plotmods$climmods$climarea$b_cov, 2, function(y) y[order(y)][length(y)*0.025])),cbind(apply(plotmods$climmods$climarea$b_cov, 2, function(y) y[order(y)][length(y)*0.975]))),
    cbind(cbind(apply(x$b_cov, 2, mean)),cbind(apply(x$b_cov, 2, function(y) y[order(y)][length(y)*0.025])),cbind(apply(x$b_cov, 2, function(y) y[order(y)][length(y)*0.975]))),
    cbind(mean(x$DIC), x$DIC[order(x$DIC)][length(x$DIC) * 0.025], x$DIC[order(x$DIC)][length(x$DIC) * 0.975]),
    cbind(mean(x$R2), x$R2[order(x$R2)][length(x$R2) * 0.025], x$R2[order(x$R2)][length(x$R2) * 0.975]),
    cbind(mean(x$RMSE), x$RMSE[order(x$RMSE)][length(x$RMSE) * 0.025], x$RMSE[order(x$RMSE)][length(x$RMSE) * 0.975])))),
  do.call("cbind",lapply(plotmods$nullmods, function(x) rbind(
    cbind(mean(x$DIC), x$DIC[order(x$DIC)][length(x$DIC) * 0.025], x$DIC[order(x$DIC)][length(x$DIC) * 0.975]),
    cbind(mean(x$R2), x$R2[order(x$R2)][length(x$R2) * 0.025], x$R2[order(x$R2)][length(x$R2) * 0.975]),
    cbind(mean(x$RMSE), x$RMSE[order(x$RMSE)][length(x$RMSE) * 0.025], x$RMSE[order(x$RMSE)][length(x$RMSE) * 0.975])))))
row.names(plotmod_out) <- c("Intercept", "Elevation", "Precipitation", "Temperature variation", "Precipitation variation", "Slope",
                            "DIC","R2","RMSE","DIC (null)","R2 (null)","RMSE (null)")
colnames(plotmod_out) <- c("mean (no spatial effect)","CI-l","CI-u","mean (cluster random effect)","CI-l","CI-u","mean (area random effect)","CI-l","CI-u")

#write.csv(round(plotmod_out, 3), "Output\\WSG\\plotmod_out.csv")

## Range of WSG along environmental predictors
mean(plotmods$climmods$clim$b_cov[1]) * (max(range(plotsF$ElevS))-min(range(plotsF$ElevS))) ; print("WSGv - Elevation (nonspatial)")
mean(plotmods$climmods$climclus$b_cov[1]) * (max(range(plotsF$ElevS))-min(range(plotsF$ElevS))) ; print("WSGv - Elevation (cluster)")
mean(plotmods$climmods$clim$b_cov[2]) * (max(range(plotsF$TotPrecS))-min(range(plotsF$TotPrecS))) ; print("WSGv - TotPrec (nonspatial)")
mean(plotmods$climmods$climclus$b_cov[2]) * (max(range(plotsF$TotPrecS))-min(range(plotsF$TotPrecS))) ; print("WSGv - TotPrec (cluster)")

## Plot effects
MeanEnv <- colMeans(plotsF[,c("ElevS","TotPrecS","TempVarS","PrecVarS","SlopeS")])
meanb <- cbind(mean(plotmods$climmods$clim$b_cov) * MeanEnv, mean(plotmods$climmods$climclus$b_cov) * MeanEnv, mean(plotmods$climmods$climarea$b_cov) * MeanEnv)
meanl <- cbind(apply(plotmods$climmods$clim$b_cov, 2, function(x) sort(x)[length(x) * 0.025]) * MeanEnv, apply(plotmods$climmods$climclus$b_cov, 2, function(x) sort(x)[length(x) * 0.025]) * MeanEnv, apply(plotmods$climmods$climarea$b_cov, 2, function(x) sort(x)[length(x) * 0.025]) * MeanEnv)
meanu <- cbind(apply(plotmods$climmods$clim$b_cov, 2, function(x) sort(x)[length(x) * 0.975]) * MeanEnv, apply(plotmods$climmods$climclus$b_cov, 2, function(x) sort(x)[length(x) * 0.975]) * MeanEnv, apply(plotmods$climmods$climarea$b_cov, 2, function(x) sort(x)[length(x) * 0.975]) * MeanEnv)

amean_elev <- sapply(plotmods$climmods, function(x) mean(x$alpha)) + colSums(meanb[-1,])
al_elev <- sapply(plotmods$climmods, function(x) sort(x$alpha)[length(x$alpha) * 0.025]) + colSums(meanb[-1,])
au_elev <- sapply(plotmods$climmods, function(x) sort(x$alpha)[length(x$alpha) * 0.975]) + colSums(meanb[-1,])

plot(rowMeans(WSGdraws_short[,-1]) ~ plotsF$ElevS,
     #col = as.numeric(as.factor(plotsF$AreaNum)),
     xlab = "Elevation", ylab = "Volume-weighted WSG")
polygon(x = c(0, 10, 10, 0), y = c(al_elev[1], al_elev[1] + sort(plotmods$climmods$clim$b_cov[,1])[length(plotmods$climmods$clim$b_cov[,1])*0.025] * 10, au_elev[1] + sort(plotmods$climmods$clim$b_cov[,1])[length(plotmods$climmods$clim$b_cov[,1])*0.975] * 10, au_elev[1]),
        border = "blue", lty = 3)
polygon(x = c(0, 10, 10, 0), y = c(al_elev[2], al_elev[2] + sort(plotmods$climmods$climclus$b_cov[,1])[length(plotmods$climmods$climclus$b_cov[,1])*0.025] * 10, au_elev[2] + sort(plotmods$climmods$climclus$b_cov[,1])[length(plotmods$climmods$climclus$b_cov[,1])*0.975] * 10, au_elev[2]),
        border = "green", lty = 3)
polygon(x = c(0, 10, 10, 0), y = c(al_elev[3], al_elev[3] + sort(plotmods$climmods$climarea$b_cov[,1])[length(plotmods$climmods$climarea$b_cov[,1])*0.025] * 10, au_elev[3] + sort(plotmods$climmods$climarea$b_cov[,1])[length(plotmods$climmods$climarea$b_cov[,1])*0.975] * 10, au_elev[3]),
        border = "red", lty = 3)

abline(amean_elev[1], mean(plotmods$climmods$clim$b_cov[,1]), col = "blue")
abline(amean_elev[2], mean(plotmods$climmods$climclus$b_cov[,1]), col = "green")
abline(amean_elev[3], mean(plotmods$climmods$climarea$b_cov[,1]), col = "red")

amean_prec <- sapply(plotmods$climmods, function(x) mean(x$alpha)) + colSums(meanb[-2,])
al_prec <- sapply(plotmods$climmods, function(x) sort(x$alpha)[length(x$alpha) * 0.025]) + colSums(meanb[-2,])
au_prec <- sapply(plotmods$climmods, function(x) sort(x$alpha)[length(x$alpha) * 0.975]) + colSums(meanb[-2,])

plot(rowMeans(WSGdraws_short[,-1]) ~ plotsF$TotPrecS,
     #col = as.numeric(as.factor(plotsF$AreaNum)),
     xlab = "Total annual precipitation (m)", ylab = "Volume-weighted WSG")
polygon(x = c(0, 10, 10, 0), y = c(al_prec[1], al_prec[1] + sort(plotmods$climmods$clim$b_cov[,2])[length(plotmods$climmods$clim$b_cov[,2])*0.025] * 10, au_prec[1] + sort(plotmods$climmods$clim$b_cov[,2])[length(plotmods$climmods$clim$b_cov[,2])*0.975] * 10, au_prec[1]),
        border = "blue", lty = 3)
polygon(x = c(0, 10, 10, 0), y = c(al_prec[2], al_prec[2] + sort(plotmods$climmods$climclus$b_cov[,2])[length(plotmods$climmods$climclus$b_cov[,2])*0.025] * 10, au_prec[2] + sort(plotmods$climmods$climclus$b_cov[,2])[length(plotmods$climmods$climclus$b_cov[,2])*0.975] * 10, au_prec[2]),
        border = "green", lty = 3)
polygon(x = c(0, 10, 10, 0), y = c(al_prec[3], al_prec[3] + sort(plotmods$climmods$climarea$b_cov[,2])[length(plotmods$climmods$climarea$b_cov[,2])*0.025] * 10, au_prec[3] + sort(plotmods$climmods$climarea$b_cov[,2])[length(plotmods$climmods$climarea$b_cov[,2])*0.975] * 10, au_prec[3]),
        border = "red", lty = 3)

abline(amean_prec[1], mean(plotmods$climmods$clim$b_cov[,2]), col = "blue")
abline(amean_prec[2], mean(plotmods$climmods$climclus$b_cov[,2]), col = "green")
abline(amean_prec[3], mean(plotmods$climmods$climarea$b_cov[,2]), col = "red")


######################################################################
#### K-fold CV models
#kfolds <- readRDS("Output\\WSG\\kfolds.RDS")
#WSGdraws <- readRDS("Output\\WSG\\WSGdraws.RDS")
#WSGdraws_short <- WSGdraws[,seq(1, ncol(WSGdraws), 40)]

## Cross-validation mean table
elpd <- lapply(kfolds, function(z) lapply(z, function(x) apply(x$elpd, 1, sum)))
RMSE <- lapply(kfolds, function(z) lapply(z, function(x) apply(x$res, 1, function(y) sqrt(sum(y) / length(y)))))
WSGdraws_expanded <- do.call("cbind",replicate(nrow(kfolds[[1]][[1]]$mu) / ncol(WSGdraws_short[,-1]), WSGdraws_short[,-1]))
R2 <- lapply(kfolds, function(y) lapply(y, function(x) sapply(1:nrow(x$mu), function(i) cor(x$mu[i, order(x$sitenum)], WSGdraws_expanded[,i])^2)))

elpdmean <- unlist(lapply(elpd, function(x) lapply(x, mean)))
RMSEmean <- unlist(lapply(RMSE, function(x) lapply(x, mean)))
R2mean <- unlist(lapply(R2, function(x) lapply(x, mean)))

CVtab <- rbind("leave-one-out" = c(elpdmean[1:2],NA,NA,RMSEmean[1:2],NA,NA,R2mean[1:2],NA,NA),
               "108-fold (cluster)" = c(elpdmean[3:6],RMSEmean[3:6],R2mean[3:6]),
               "23-fold (area)" = c(elpdmean[7:10],RMSEmean[7:10],R2mean[7:10]),
               "7-fold (region)" = c(elpdmean[11:14],RMSEmean[11:14],R2mean[11:14]))
colnames(CVtab) <- c("elpdS","elpdSN","elpdR","elpdRN","RMSPES","RMSPESN","RMSPER","RMSPERN","pR2S","pR2SN","pR2R","pR2RN")

CVtab <- CVtab[, c("elpdRN","elpdR","RMSPER","pR2R","elpdSN","elpdS","RMSPES","pR2S")]

#write.csv(CVtab,"Output\\WSG\\plotmod_CVtab.csv")

######################
## AGB and CV table ## 
######################
plotsF <- plotsF[order(plotsF$SiteNum),]
#WSGdraws <- readRDS("Output\\WSG\\WSGdraws.RDS")

WSGdraws_plot <- WSGdraws[,-1]
WSGdraws_area <- lapply(split(WSGdraws_plot, plotsF$AreaNum), function(x) do.call("rbind",x))
WSGdraws_cluster <- lapply(split(WSGdraws_plot, plotsF$ClusNum), function(x) do.call("rbind",x))
WSGdraws_region <- lapply(split(WSGdraws_plot, plotsF$RegionNum), function(x) do.call("rbind",x))

AGBdraws_plot <- WSGdraws_plot * plotsF$volha
AGBdraws_cluster <- lapply(split(AGBdraws_plot, plotsF$ClusNum), function(x) do.call("rbind",x))
AGBdraws_area <- lapply(split(AGBdraws_plot, plotsF$AreaNum), function(x) do.call("rbind",x))
AGBdraws_region <- lapply(split(AGBdraws_plot, plotsF$RegionNum), function(x) do.call("rbind",x))

# AGB by region
AGBtab <- rbind(
  cbind(
    cbind(
      do.call("rbind",lapply(AGBdraws_region, function(x) mean(apply(x, 1, mean)))),
      do.call("rbind",lapply(AGBdraws_region, function(x) sort(apply(x, 1, mean))[nrow(x)*0.025])),
      do.call("rbind",lapply(AGBdraws_region, function(x) sort(apply(x, 1, mean))[nrow(x)*0.975]))),
    cbind(
      do.call("rbind",lapply(AGBdraws_region, function(x) mean(apply(x, 1, function(y) sd(y)/mean(y))))),
      do.call("rbind",lapply(AGBdraws_region, function(x) sort(apply(x, 1, function(y) sd(y)/mean(y)))[nrow(x)*0.025])),
      do.call("rbind",lapply(AGBdraws_region, function(x) sort(apply(x, 1, function(y) sd(y)/mean(y)))[nrow(x)*0.975]))),
    as.vector(as.data.frame(plotsF %>% group_by(RegionNum) %>% summarise(CVvol = sd(volha)/mean(volha)))$CVvol),
    cbind(
      do.call("rbind",lapply(WSGdraws_region, function(x) mean(apply(x, 1, function(y) sd(y)/mean(y))))),
      do.call("rbind",lapply(WSGdraws_region, function(x) sort(apply(x, 1, function(y) sd(y)/mean(y)))[nrow(x)*0.025])),
      do.call("rbind",lapply(WSGdraws_region, function(x) sort(apply(x, 1, function(y) sd(y)/mean(y)))[nrow(x)*0.975])))
  ), as.vector(c(
    mean(apply(AGBdraws_plot, 2, function(x) mean(x))),
    sort(apply(AGBdraws_plot, 2, function(x) mean(x)))[ncol(AGBdraws_plot)*0.025],
    sort(apply(AGBdraws_plot, 2, function(x) mean(x)))[ncol(AGBdraws_plot)*0.975],
    mean(apply(AGBdraws_plot, 2, function(x) sd(x) / mean(x))),
    sort(apply(AGBdraws_plot, 2, function(x) sd(x) / mean(x)))[ncol(AGBdraws_plot)*0.025],
    sort(apply(AGBdraws_plot, 2, function(x) sd(x) / mean(x)))[ncol(AGBdraws_plot)*0.975],
    sd(plotsF$volha) / mean(plotsF$volha),
    mean(apply(WSGdraws_plot, 2, function(x) sd(x) / mean(x))),
    sort(apply(WSGdraws_plot, 2, function(x) sd(x) / mean(x)))[ncol(WSGdraws_plot)*0.025],
    sort(apply(WSGdraws_plot, 2, function(x) sd(x) / mean(x)))[ncol(WSGdraws_plot)*0.975])
  ), as.vector(c(
    NA,NA,NA,
    mean(apply(do.call("cbind",lapply(AGBdraws_cluster, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))),
    sort(apply(do.call("cbind",lapply(AGBdraws_cluster, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y)))[length(sort(apply(do.call("cbind",lapply(AGBdraws_cluster, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))))*0.025],
    sort(apply(do.call("cbind",lapply(AGBdraws_cluster, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y)))[length(sort(apply(do.call("cbind",lapply(AGBdraws_cluster, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))))*0.975],
    sd(as.data.frame(plotsF %>% group_by(ClusNum) %>% summarise(meanvol = mean(volha)))$meanvol) / mean(as.data.frame(plotsF %>% group_by(ClusNum) %>% summarise(meanvol = mean(volha)))$meanvol),
    mean(apply(do.call("cbind",lapply(WSGdraws_cluster, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))),
    sort(apply(do.call("cbind",lapply(WSGdraws_cluster, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y)))[length(sort(apply(do.call("cbind",lapply(WSGdraws_cluster, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))))*0.025],
    sort(apply(do.call("cbind",lapply(WSGdraws_cluster, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y)))[length(sort(apply(do.call("cbind",lapply(WSGdraws_cluster, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))))*0.975])
  ), as.vector(c(
    NA,NA,NA,
    mean(apply(do.call("cbind",lapply(AGBdraws_area, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))),
    sort(apply(do.call("cbind",lapply(AGBdraws_area, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y)))[length(sort(apply(do.call("cbind",lapply(AGBdraws_area, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))))*0.025],
    sort(apply(do.call("cbind",lapply(AGBdraws_area, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y)))[length(sort(apply(do.call("cbind",lapply(AGBdraws_area, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))))*0.975],
    sd(as.data.frame(plotsF %>% group_by(AreaNum) %>% summarise(meanvol = mean(volha)))$meanvol) / mean(as.data.frame(plotsF %>% group_by(AreaNum) %>% summarise(meanvol = mean(volha)))$meanvol),
    mean(apply(do.call("cbind",lapply(WSGdraws_area, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))),
    sort(apply(do.call("cbind",lapply(WSGdraws_area, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y)))[length(sort(apply(do.call("cbind",lapply(WSGdraws_area, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))))*0.025],
    sort(apply(do.call("cbind",lapply(WSGdraws_area, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y)))[length(sort(apply(do.call("cbind",lapply(WSGdraws_area, function(x) apply(x, 1, mean))), 1, function(y) sd(y)/mean(y))))*0.975])
  ))

colnames(AGBtab) <- c("AGB","AGB-l","AGB-u","CV (AGB)","CV-l (AGB)","CV-u (AGB)","CV (vol)","CV (WSG)","CV-l (WSG)","CV-u (WSG)")
rownames(AGBtab) <- c(unique(plotsF$Region)[order(unique(plotsF$RegionNum))], "Overall (between-plot)","Overall (between-cluster)","Overall (between-area)")

#write.csv(round(AGBtab,2) ,"Output\\WSG\\AGBtab.csv")

#############################
## Error/bias in AGB table ##
#############################
plotsF <- plotsF[order(plotsF$SiteNum),]
#WSGdraws <- readRDS("Output\\WSG\\WSGdraws.RDS")
#kfolds <- readRDS("Output\\WSG\\kfolds.RDS")

# Extract relative errors from k-fold CV of plot models
WSGmoderror <- lapply(lapply(kfolds, function(x)x$spatialfold), function(z)  z$prederror[,order(z$sitenum)])
AGBmoderror <- lapply(WSGmoderror, function(x) apply(x, 1, function(y) y * plotsF$volha))
AGBmoderror <- lapply(AGBmoderror, function(x) as.data.frame(x))

# Calculate errors from using spatial average
AGBdraws_plot <- WSGdraws[,-1] * plotsF$volha
AGBdraws_cluster <- lapply(split(AGBdraws_plot, plotsF$ClusNum), function(x) do.call("rbind",x))
AGBdraws_area <- lapply(split(AGBdraws_plot, plotsF$AreaNum), function(x) do.call("rbind",x))
AGBdraws_region <- lapply(split(AGBdraws_plot, plotsF$RegionNum), function(x) do.call("rbind",x))

AGBmeanerror <- list()
#AGBmeanerror$cluster <- t(do.call("cbind", lapply(AGBdraws_cluster, function(x) rowMeans(x) - x)))[order(plotsF$ClusNum)]      # PLOT-WISE DIFFERENCE FROM CLUSTER MEAN
AGBmeanerror$area <- t(do.call("cbind", lapply(AGBdraws_area, function(x) rowMeans(x) - x)))[order(plotsF$AreaNum),]          # PLOT-WISE DIFFERENCE FROM AREA MEAN
AGBmeanerror$region <- t(do.call("cbind", lapply(AGBdraws_region, function(x) rowMeans(x) - x)))[order(plotsF$RegionNum),]      # PLOT-WISE DIFFERENCE FROM REGION MEAN
AGBmeanerror$overall <- apply(AGBdraws_plot, 2, function(x) mean(x) - x)                               # PLOT-WISE DIFFERENCE FROM OVERALL MEAN
AGBmeanerror <- lapply(AGBmeanerror, as.data.frame)

# Table
AGBmoderror_byreg <- lapply(AGBmoderror, function(x) lapply(split(x, plotsF$RegionNum), function(y) do.call("rbind",y)))
AGBmeanerror_byreg <- lapply(AGBmeanerror, function(x) lapply(split(x, plotsF$RegionNum), function(y) do.call("rbind",y)))

MSDtab <- rbind(
  cbind(
    do.call("cbind", lapply(lapply(AGBmoderror_byreg, function(x) lapply(x, function(y) mean(rowMeans(y)))), unlist)),
    do.call("cbind", lapply(lapply(AGBmeanerror_byreg, function(x) lapply(x, function(y) mean(rowMeans(y)))), unlist))
  ),
  cbind(
    do.call("cbind", lapply(AGBmoderror, function(x) mean(as.matrix(x)))),
    do.call("cbind", lapply(AGBmeanerror, function(x) mean(as.matrix(x))))))

MSDtab_CIl <- rbind(
  cbind(
    do.call("cbind", lapply(lapply(AGBmoderror_byreg, function(x) lapply(x, function(y) rowMeans(y)[order(rowMeans(y))][nrow(y)*0.025])), unlist)),
    do.call("cbind", lapply(lapply(AGBmeanerror_byreg, function(x) lapply(x, function(y) rowMeans(y)[order(rowMeans(y))][nrow(y)*0.025])), unlist))
  ),
  cbind(
    do.call("cbind", lapply(AGBmoderror, function(x) colMeans(x)[order(colMeans(x))][ncol(x)*0.025])),
    do.call("cbind", lapply(AGBmeanerror, function(x) colMeans(x)[order(colMeans(x))][ncol(x)*0.025]))))

MSDtab_CIu <- rbind(
  cbind(
    do.call("cbind", lapply(lapply(AGBmoderror_byreg, function(x) lapply(x, function(y) rowMeans(y)[order(rowMeans(y))][nrow(y)*0.975])), unlist)),
    do.call("cbind", lapply(lapply(AGBmeanerror_byreg, function(x) lapply(x, function(y) rowMeans(y)[order(rowMeans(y))][nrow(y)*0.975])), unlist))
  ),
  cbind(
    do.call("cbind", lapply(AGBmoderror, function(x) colMeans(x)[order(colMeans(x))][ncol(x)*0.975])),
    do.call("cbind", lapply(AGBmeanerror, function(x) colMeans(x)[order(colMeans(x))][ncol(x)*0.975]))))

MSDtab_full <- data.frame(NA, nrow = nrow(MSDtab))
cnames <- vector()
for(i in 1:ncol(MSDtab)){
  MSDtab_full <- cbind(MSDtab_full,cbind(MSDtab[,i],MSDtab_CIl[,i],MSDtab_CIu[,i]))
  #MSDtab_full <- cbind(MSDtab_full, paste(round(MSDtab[,i], 2), " (", round(MSDtab_CIl[,i], 2), " - ", round(MSDtab_CIu[,i], 2), ")", sep=""))
  cnames <- append(cnames,c(colnames(MSDtab)[i],"CI-l","CI-u"))
}
MSDtab_full[,1:2] <- NULL
colnames(MSDtab_full) <- cnames
#colnames(MSDtab_full) <- colnames(MSDtab)
rownames(MSDtab_full) <- c(unique(plotsF$Region[order(plotsF$RegionNum)]),"Overall")

#write.csv(round(MSDtab_full, 2), "Output\\WSG\\MSDtab.csv")

##############################
## Predict plot-average WSG ## 
##############################
WSGdraws <- readRDS("Output\\WSG\\WSGdraws.RDS")
WSGtreemod <- WSGtreemod <- readRDS("Output\\WSG\\WSGtreemod.RDS")
plotsF <- plotsF[order(plotsF$SiteNum),]

## Predict values across dataset
plotsF$WSGpred <- rowMeans(WSGdraws[,-1])
treesF$WSGused <- WSGtreemod$mean$WSG.pred
# treesF$WSGused[which(!is.na(treesF$WSG))] <- treesF$WSG[which(!is.na(treesF$WSG))]

# ## Weigh WSG by volume
# treesF$WSGv <- treesF$WSGused * treesF$VOLw

plot(treesF$WSG[which(!is.na(treesF$WSG))],treesF$WSGused[which(!is.na(treesF$WSG))])

## Calculate plot averages
#plotsF <- treesF %>% group_by(SiteCode) %>% summarise(WSGn = mean(WSG,na.rm=T), WSGp = mean(WSGused), WSGv = sum(WSGv), WSGvm = sum(WSGvm,na.rm=T))
#plotsF <- left_join(plotsF, treesF %>% group_by(SiteCode) %>% summarise(WSGn = mean(WSG,na.rm=T), WSGp = mean(WSGused), WSGv = sum(WSGv), WSGvm = sum(WSGvm,na.rm=T)))

## Figures
temp <- plotsF[-which(plotsF$ncore == 0),]
cols <- RColorBrewer::brewer.pal(9, "Blues")[ceiling(plotsF$ncore/7) + 1]
plot((plotsF$WSGv-plotsF$WSGvm) ~ plotsF$ncore, ylab = "Estimated - measured plot average WSGv", xlab = "Number of cores") ; abline(h = 0)
plot(plotsF$WSGvm~plotsF$WSGv, col = cols, pch=16) ; abline(0,1)

#########################################
## Raw data overview table and figures ## 
#########################################
WSGdraws <- readRDS("Output\\WSG\\WSGdraws.RDS")
WSGtreemod <- WSGtreemod <- readRDS("Output\\WSG\\WSGtreemod.RDS")

plotsF <- plotsF[order(plotsF$SiteNum),]
plotsF$WSGpred <- rowMeans(WSGdraws[,-1])

treesF$WSGused <- WSGtreemod$mean$WSG.pred
treesF$WSGused[which(!is.na(treesF$WSG))] <- treesF$WSG[which(!is.na(treesF$WSG))]

## This part adds volume-weighted WSG for cored trees only
treesFm <- filter(treesF, !is.na(WSG))
treesFm <- left_join(treesFm, treesFm %>% group_by(SiteCode) %>% summarise(totDBHm = sum(DBH_used), totVOLm = sum(vol)))
treesFm$VOLw <- treesFm$vol / treesFm$totVOLm
treesFm$WSGvm <- treesFm$WSG * treesFm$VOLw
treesF <- left_join(treesF, treesFm[,c("SiteCode","TreeN","StemN","totVOLm","WSGvm")])
plotsF <- left_join(plotsF, treesF %>% group_by(SiteCode) %>% summarise(WSGm = mean(WSG, na.rm=T), WSGvm = sum(WSGvm, na.rm = T)))

## Summary table
datatab <- rbind(
           as.data.frame(plotsF %>% group_by(Region, AreaNum) %>% summarise(Lat = mean(lat), Long = mean(long), Plots = n(),
                         Volume = mean(volha), WSG = mean(WSGm, na.rm = T), WSGv = mean(WSGvm),
                         ntrees = sum(ntree), cores = sum(ncore)/sum(ntree), Quercus = sum(nQuercus)/sum(ntree), Species = sum(nSpec)/sum(ntree),
                         Elevation = mean(ALOSelev), elevmin = min(ALOSelev), elevmax = max(ALOSelev),
                         Precipitation = mean(TotPrec), precmin = min(TotPrec), precmax = max(TotPrec),
                         Tempvar = mean(TempVar), tempvarmin = min(TempVar), tempvarmax = max(TempVar),
                         Precvar = mean(PrecVar), precvarmin = min(PrecVar), precvarmax = max(PrecVar))),
           as.data.frame(plotsF %>% group_by(Region) %>% summarise(AreaNum = NA, Lat = mean(lat), Long = mean(long), Plots = n(),
                         Volume = mean(volha), WSG = mean(WSGm, na.rm = T), WSGv = mean(WSGvm),
                         ntrees = sum(ntree), cores = sum(ncore)/sum(ntree), Quercus = sum(nQuercus)/sum(ntree), Species = sum(nSpec)/sum(ntree),
                         Elevation = mean(ALOSelev), elevmin = min(ALOSelev), elevmax = max(ALOSelev),
                         Precipitation = mean(TotPrec), precmin = min(TotPrec), precmax = max(TotPrec),
                         Tempvar = mean(TempVar), tempvarmin = min(TempVar), tempvarmax = max(TempVar),
                         Precvar = mean(PrecVar), precvarmin = min(PrecVar), precvarmax = max(PrecVar))))
datatab <- datatab[order(datatab$Region, datatab$AreaNum, na.last = FALSE),]
datatab[,-c(1,2)] <- round(datatab[,-c(1,2)],2)
datatab[,-c(1:12)] <- round(datatab[,-c(1:12)])
datatab$Elevation <- paste(datatab$Elevation, " (", datatab$elevmin, " - ", datatab$elevmax, ")",sep="")
datatab$Precipitation <- paste(datatab$Precipitation, " (", datatab$precmin, " - ", datatab$precmax, ")",sep="")
datatab$Tempvar <- paste(datatab$Tempvar, " (", datatab$tempvarmin, " - ", datatab$tempvarmax, ")",sep="")
datatab$Precvar <- paste(datatab$Precvar, " (", datatab$precvarmin, " - ", datatab$precvarmax, ")",sep="")
datatab[,c("elevmin","elevmax","precmin","precmax","tempvarmin","tempvarmax","precvarmin","precvarmax")] <- NULL

#write.csv(datatab, "Output\\WSG\\datatab.csv", row.names = F)

## Histograms and CVs of WSG - noID, ID, Quercus
treesQuercus <- dplyr::filter(treesF, Species == "Quercus Humboldtii")
treesID <- dplyr::filter(treesF, !Species %in% c("Quercus Humboldtii","nospec"))
treesNOID <- dplyr::filter(treesF, Species == "nospec")

tiff(file = "Output\\WSG\\WSGdistribution.tiff", width = 10000, height = 5000, res = 900)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), widths = c(0.5,0.5), height = c(0.8, 1.2))
par(mar = c(4, 4, 0.5, 0.5))
  hist(treesNOID$WSG, xlab = "WSG", main = NA, xlim = c(0.1, 1.4))
    legend("topright", legend = c(paste("Unidentified trees (n = ", length(which(!is.na(treesNOID$WSG))), ")", sep = ""), 
                                 paste("Coefficient of variation: ", round(sd(treesNOID$WSG, na.rm=T)/mean(treesNOID$WSG, na.rm=T), 2), sep = "")),
           ncol = 1, bty = "n")
    abline(v = mean(treesNOID$WSG, na.rm = T), col = "red")
  hist(treesQuercus$WSG, xlab = "WSG", main = NA, xlim = c(0.1, 1.4))
    legend("topright", legend = c(paste("Quercus Humboldtii (n = ", length(which(!is.na(treesQuercus$WSG))), ")", sep = ""), 
                                  paste("Coefficient of variation: ", round(sd(treesQuercus$WSG, na.rm=T)/mean(treesQuercus$WSG, na.rm=T), 2), sep = "")),
           ncol = 1, bty = "n")
    abline(v = mean(treesQuercus$WSG, na.rm = T), col = "red")
  boxplot(WSG ~ Species, treesID, xlab = "Identifier", xaxt = "n")
    legend("topleft", legend = c(paste("Trees with plot-specific identifier (n = ", length(which(!is.na(treesID$WSG))), ")", sep = ""), 
                                  paste("Coefficient of variation: ", round(mean(data.frame(treesID %>% group_by(Species) %>% summarise(CV = sd(WSG, na.rm = T) / mean(WSG, na.rm = T)))$CV), 2), sep = "")),
           ncol = 1, bty = "n")

dev.off()

## Tree size
treesF$DBH2 <- ifelse(treesF$DBH > 50, 50, treesF$DBH)
treesF$DBHwsg <- ifelse(is.na(treesF$WSG), NA, treesF$DBH2)

ggplot(treesF)+
  theme_classic()+
  scale_x_continuous(breaks = seq(5, 50, 5), expand = c(0,0), 
                     labels = c("5-10","10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50","50+"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0,0))+
  geom_histogram(aes(x = DBH2), fill = "darkgreen", binwidth = 5)+
  geom_histogram(aes(x = DBHwsg), fill = "darkorange", binwidth = 5)+
  labs(x = "DBH range", y = "Count")

ggsave("Output\\WSG\\tree_size_WSG.tiff", width = 11, height = 10, units = "cm")

## Observed vs. estimates WSG figure
temp <- plotsF[-which(plotsF$ncore == 0),]
cols <- RColorBrewer::brewer.pal(9, "Greens")[ceiling(plotsF$ncore/7) + 1]

tiff(file = "Output\\WSG\\WSG_scatterplots.tiff", width = 10000, height = 5000, res = 900)

layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths = c(0.5,0.5))
par(mar = c(4, 4, 0.5, 0.5))
    ## Scatter plot: observed vs. estimated individual WSG
    plot(colMeans(WSGtreemod$sims.list$WSG.pred[,which(!is.na(treesF$WSG))]) ~ treesF$WSG[which(!is.na(treesF$WSG))],
         xlab = "Measured WSG", ylab = "Estimated WSG")
    abline(0,1)
    legend("topright", bty = "n",
           paste("Pearson's r =", round(cor(treesF$WSG[which(!is.na(treesF$WSG))], colMeans(WSGtreemod$sims.list$WSG.pred[,which(!is.na(treesF$WSG))])), 2)))
    ## Scatter plot: observed vs. estimated plot WSGv
    plot(plotsF$WSGpred ~ plotsF$WSGvm, col = cols, pch=16,
         xlab = "Measured WSGv (average of cores)", ylab = "Estimated WSGv",
         xlim = c(min(plotsF$WSGpred, plotsF$WSGvm),max(plotsF$WSGpred, plotsF$WSGvm)), ylim = c(min(plotsF$WSGpred, plotsF$WSGvm),max(plotsF$WSGpred, plotsF$WSGvm)))
    abline(0,1)
    legend("bottomright", legend = lapply(split(plotsF$ncore,ceiling(plotsF$ncore/7)+1),min), pch=16, pt.cex=1.5, cex=1, bty='n',
           col = unique(cols)[order(unique(ceiling(plotsF$ncore/7) + 1))], title = "Number of cores", title.adj = -0.005)
    legend("topleft", bty = "n",
           paste("Pearson's r =", round(cor(plotsF$WSGpred, plotsF$WSGvm), 2)))

dev.off()

## WSG vs. size (trees and plots) NOT FINISHED
plot(WSGtreemod$mean$WSG.pred[which(!is.na(treesF$WSG))] ~ treesF$vol[which(!is.na(treesF$WSG))])
abline(WSGtreemod$q50$alpha,WSGtreemod$q50$b_vol)
abline(WSGtreemod$q2.5$alpha,WSGtreemod$q2.5$b_vol)
abline(WSGtreemod$q97.5$alpha,WSGtreemod$q97.5$b_vol)

plot(rowMeans(WSGdraws[,-1]) ~ plotsF$vol)
abline(lm(rowMeans(WSGdraws[,-1]) ~ plotsF$volha))

##################
## R-hat checks ##
##################
#WSGtreemod <- readRDS("Output\\WSG\\WSGtreemod.RDS")
lapply(WSGtreemod$Rhat, max)

#plotmods_rhats <- readRDS("Output\\WSG\\plotmods_rhats.RDS")
lapply(plotmods_rhats, function(x) lapply(x, function(y) lapply(y, max)))

#kfolds <- readRDS("Output\\WSG\\kfolds.RDS")
lapply(kfolds, function(x) lapply(x,function(y) max(y$rhatmax_alpha)))
lapply(kfolds, function(x) lapply(x,function(y) max(y$rhatmax_var0)))
lapply(kfolds, function(x) lapply(x,function(y) apply(y$rhatmax_bcov, 2, max)))

############################
## Posterior model checks ##
############################
#WSGtreemod <- readRDS("Output\\WSG\\WSGtreemod.RDS")
#plotmods <- readRDS("Output\\WSG\\plotmods.RDS")
#WSGdraws <- readRDS("Output\\WSG\\WSGdraws.RDS")

y <- treesF$WSG[which(!is.na(treesF$WSG))]
yrep <- WSGtreemod$sims.list$WSG.pred[,which(!is.na(treesF$WSG))]
group <- treesF$speciesid[which(!is.na(treesF$WSG))]

y <- rowMeans(WSGdraws[,-1])
yrep <- plotmods$climmods$clim$WSGpred
group <- plotsF$RegionNum

bayesplot::ppc_dens_overlay(y = y, yrep = yrep[1:100, ])
bayesplot::ppc_hist(y, yrep[1:10, ])
bayesplot::ppc_stat(y, yrep, stat = "max")
bayesplot::ppc_stat(y, yrep, stat = "min")
bayesplot::ppc_stat(y, yrep, stat = "mean")
bayesplot::ppc_stat(y, yrep, stat = "median")
bayesplot::ppc_stat_grouped(y, yrep, group = group, stat = "mean")

bayesplot::ppc_error_scatter(y, yrep[1:10,])
bayesplot::ppc_error_scatter_avg(y, yrep)
bayesplot::ppc_scatter_avg(y, yrep)
bayesplot::ppc_scatter_avg_grouped(y, yrep, group = group)

## Moran's I
residuals <- colMeans(plotmods$climmods$clim$res)
residuals <- WSGtreemod$mean$res

coordmat <- as.matrix(dist(cbind(plotsF$long, plotsF$lat)))
coordmat <- as.matrix(dist(cbind(treesF$long, treesF$lat)))

dists <- coordmat + 1
dists.inv <- 1/dists
diag(dists.inv) <- 0
dists.inv[is.infinite(dists.inv)] <- 0

as.data.frame(ape::Moran.I(residuals, dists.inv, na.rm = T))

spatialEco::morans.plot(residuals,coords = cbind(plotsF$long,plotsF$lat))
spatialEco::morans.plot(residuals,coords = cbind(treesF$long,treesF$lat))

