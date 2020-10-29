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
hc <- hclust(as.dist(distmat), "single")
plot(hc, labels = plotdata$SiteCode) ; abline(h = seq(0,50,10), col = "red") ; abline(h = c(15,25), col = "blue")
raster::plot(coords, col = cutree(hc, h = 10), pch = 18)

plotdata$AreaNum <- cutree(hc, h = 10)      # Merges AL/SM, AT/IG/GA, TA/TS/VI/VC    # Splits PS

## Slight changes to regions
plotdata$Region[which(plotdata$Region == "oriental na")] <- "oriental western slope"
plotdata$Region[which(plotdata$AreaCode == "CQ")] <- "oriental eastern slope"
plotdata$RegionNum <- as.numeric(as.factor(plotdata$Region))

### Add grouping levels as numbers
treesF$SiteNum <- as.numeric(as.factor(treesF$SiteCode))
treesF$ClusNum <- as.numeric(as.factor(treesF$Cluster))
treesF <- left_join(treesF, plotdata[,c("SiteCode","AreaNum")])
treesF <- left_join(treesF, plotdata[,c("SiteCode","RegionNum")])

### Plot dataset
plotsF <- treesF %>% group_by(SiteCode) %>% summarise(SiteNum = first(SiteNum), ClusNum = first(ClusNum), AreaNum = first(AreaNum), RegionNum = first(RegionNum),
                                                      ntree = n(), ncore = sum(WSG/WSG,na.rm=T), vol = sum(vol))
plotsF <- as.data.frame(left_join(plotsF,Spatial))
plotsF <- data.frame(plotsF, "ElevS" = plotsF$ALOSelev / 1000, "TotPrecS" = plotsF$TotPrec/1000, "TempVarS" = scale(plotsF$TempVar),
                     "PrecVarS" = scale(plotsF$PrecVar), "SlopeS" = plotsF$ALOSslope/10)
plotsF$volha <- (plotsF$vol / plotsF$Size) * 10000

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

# saveRDS(WSGtreemod, "JAGS\\WSGtreemod.RDS")
# WSGtreemod <- readRDS("JAGS\\WSGtreemod.RDS")

## Posterior checks
y <- treesF$WSG[which(!is.na(treesF$WSG))]
yrep <- WSGtreemod$sims.list$WSG.pred[,which(!is.na(treesF$WSG))]
bayesplot::ppc_dens_overlay(y = y, yrep = yrep[1:100, ])
bayesplot::ppc_hist(y, yrep[1:10, ])
bayesplot::ppc_stat(y, yrep, stat = "max")
bayesplot::ppc_stat(y, yrep, stat = "min")
bayesplot::ppc_stat(y, yrep, stat = "mean")
bayesplot::ppc_stat(y, yrep, stat = "median")
bayesplot::ppc_stat_grouped(y, yrep, group = treesF$speciesid[which(!is.na(treesF$WSG))], stat = "mean")

bayesplot::ppc_error_scatter(y, yrep[1:10,])
bayesplot::ppc_error_scatter_avg(y, yrep)
bayesplot::ppc_scatter_avg(y, yrep)

## Scatter plot
plot(y ~ colMeans(yrep),
     xlab = "Measured WSG", ylab = "Estimated WSG")
abline(0,1)

## Whisker plot
# Spatial random effects output dataframes
b_area <- data.frame(AreaNum = 1:length(unique(treesF$AreaNum)), bArea = WSGtreemod$mean$b_area, bArea95l = WSGtreemod$q2.5$b_area, barea95u = WSGtreemod$q97.5$b_area)
b_cluster <- data.frame(ClusNum = 1:length(unique(treesF$Cluster)), bCluster = WSGtreemod$mean$b_cluster, bCluster95l = WSGtreemod$q2.5$b_cluster, bCluster95u = WSGtreemod$q97.5$b_cluster)
b_site <- data.frame(SiteNum = 1:length(unique(treesF$SiteCode)), bSite = WSGtreemod$mean$b_site, bSite95l = WSGtreemod$q2.5$b_site, bSite95u = WSGtreemod$q97.5$b_site, AreaNum = data.frame(treesF %>% group_by(SiteNum) %>% summarise(first(AreaNum)))[,2], ClusNum = data.frame(treesF %>% group_by(SiteNum) %>% summarise(first(ClusNum)))[,2], RegionNum = data.frame(treesF %>% group_by(SiteNum) %>% summarise(first(RegionNum)))[,2])

# Sort dataframes by increasing WSG
b_site <- left_join(b_site, b_site %>% group_by(RegionNum) %>% summarise(RegOrd = mean(bSite)))
b_site <- left_join(b_site, b_site %>% group_by(AreaNum) %>% summarise(AreaOrd = mean(bSite)))
b_site <- left_join(b_site, b_site %>% group_by(ClusNum) %>% summarise(ClusOrd = mean(bSite)))
b_site <- b_site[order(b_site$RegOrd, b_site$AreaOrd, b_site$ClusOrd, b_site$bSite),]

# Plot by site
#par(mfrow=c(3,1))
#tiff(file = "Output\\WSG\\plotwhisker.tiff", width = 6000, height = 3000, res = 900)
  plot(NA, ylim = c(min(b_site$bSite95l), max(b_site$bSite95u)), xlim = c(0,nrow(b_site)), xlab = "Site index", ylab = "Random effect")
  for(i in 1:nrow(b_site)){segments(x0=i, y0=b_site$bSite95l[i], x1=i, y1=b_site$bSite95u[i])}
  points(b_site$bSite, pch = 19, cex = 0.7, col = b_site$AreaNum)
  abline(h = 0)
#dev.off()
  
#plot(NA, ylim = c(min(b_cluster$bCluster95l), max(b_cluster$bCluster95u)), xlim = c(0,nrow(b_cluster)), xlab = "cluster index", ylab = "Random effect")
#for(i in 1:nrow(b_cluster)){segments(x0=i, y0=b_cluster$bCluster95l[i], x1=i, y1=b_cluster$bCluster95u[i])}
#points(b_cluster$bCluster, pch = 19, cex = 0.7, col = b_cluster$AreaNum)
#abline(h = 0)

## Model fit
R2 <- apply(WSGtreemod$sims.list$mu[,which(!is.na(treesF$WSG))], 1, function(x) cor(x, treesF$WSG[which(!is.na(treesF$WSG))])^2)
RMSE <- apply(WSGtreemod$sims.list$res[,which(!is.na(treesF$WSG))], 1, function(x) sqrt(sum(x^2)/length(x)))

## Add cluster/site ranges
treesF <- cbind(treesF, WSGpred = WSGtreemod$mean$WSG.pred)
areaavg <- as.data.frame(treesF %>% group_by(AreaCode) %>% summarise(WSGpred = mean(WSGpred)))$WSGpred
clusteravg <- as.data.frame(treesF %>% group_by(AreaCode,Cluster) %>% summarise(WSGpred = mean(WSGpred)))$WSGpred

## Model estimate table
treemod_esttab <- WSGtreemod$summary[c("alpha","b_vol","vararea","varcluster","varsite","varspec","var0","var1"),c(1,3,7)]
rownames(treemod_esttab) <- c("Intercept","Volume (m3)","Area","Cluster","Site","Species","Residual (no species)","Residual (species)")
treemod_esttab <- rbind(treemod_esttab,
                        "RMSE" = c(mean(RMSE), RMSE[order(RMSE)][length(RMSE)*0.025], RMSE[order(RMSE)][length(RMSE)*0.975]),
                        "R2" = c(mean(R2), R2[order(R2)][length(R2)*0.025], R2[order(R2)][length(R2)*0.975]))

treemod_esttab <- round(treemod_esttab,4)
#write.csv(treemod_esttab,"Output\\WSG\\treemod_esttab.csv")

## Moran's I
tree.dists <- as.matrix(dist(cbind(treesF$long, treesF$lat)))
tree.dists <- tree.dists + 1
tree.dists.inv <- 1/tree.dists
diag(tree.dists.inv) <- 0
tree.dists.inv[is.infinite(tree.dists.inv)] <- 0
ape::Moran.I(WSGtreemod$mean$res, tree.dists.inv)

#############################
## Plot-average WSG models ##
#############################
## Extract posterior draws of volume-weighted plot averages
WSGdraws <- WSGtreemod$sims.list$WSG.pred
WSGdraws[, which(!is.na(treesF$WSG))] <- matrix(data = rep(treesF$WSG[which(!is.na(treesF$WSG))], nrow(WSGdraws)), nrow = nrow(WSGdraws), ncol = length(treesF$WSG[which(!is.na(treesF$WSG))]), byrow = TRUE)
WSGdraws <- as.data.frame(t(WSGdraws) * treesF$VOLw)
WSGdraws <- data.frame(WSGdraws) %>% mutate(SiteNum = treesF$SiteNum) %>% reshape2::melt(id.vars = "SiteNum") %>% reshape2::dcast(SiteNum ~ variable, value.var="value", fun.aggregate = sum)
WSGdraws_short <- WSGdraws[,seq(1, ncol(WSGdraws), 40)]

## Function to run models
WSGmod_multicov <- function(covariates = NULL, spatscale = NULL, modelfile = NA){
    moddat <- list("Nobs" = nrow(WSGdraws_short),
                   "WSG" = WSGdraws_short[order(WSGdraws_short$SiteNum),-1],
                   "ndraws" = ncol(WSGdraws_short) - 1)
  #moddat <- list("Nobs" = nrow(plotsF),
  #               "WSG" = plotsF$WSGv)
  
  if(!is.null(spatscale)){
    moddat <- append(moddat, list("SpatialNum" = spatscale,
                                  "nspatial" = length(unique(spatscale))))}
  
  if(!is.null(covariates)){
    moddat <- append(moddat, list("cov" = covariates,
                                  "ncov" = ncol(covariates)))}
  
  mod <- jagsUI::jags(model.file = modelfile, data = moddat, inits = modinits,#inits = NULL,
                      parameters.to.save = c("alpha","b_cov","varspatial","var0","b_spatial",
                                             "WSG.pred","loglik","MSE","mu"),
                      n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                      parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res"))
  return(mod)
}

## Run model
n.chains <- 4
n.burnin <- 5000
n.iter <- n.burnin + 50000
n.thin <- 50
n.cores <- n.chains

#area
modinits <- function(){list(alpha = rnorm(1, 0.5, 0.001), tauspatial = 1/rnorm(1, 0.004, 0.0001), tau = 1/rnorm(1, 0.004, 0.0001), 
                            b_spatial = sapply(plotmods$simpmod$area$mean$b_spatial, function(x) rnorm(1, x, 0.01)))}
#cluster
modinits <- function(){list(alpha = rnorm(1, 0.507, 0.001), tauspatial = 1/rnorm(1, 0.003, 0.0001), tau = 1/rnorm(1, 0.003, 0.0001), 
                            b_spatial = sapply(plotmods$simpmod$clus$mean$b_spatial, function(x) rnorm(1, x, 0.01)))}

plotmods$simpmod$area <- WSGmod_multicov(modelfile = "JAGS//WSG_plots_test.txt", spatscale = plotsF$AreaNum)
plotmods$simpmod$clus <- WSGmod_multicov(modelfile = "JAGS//WSG_plots_test.txt", spatscale = plotsF$ClusNum)

plotmods <- list()
plotmods$nullmods$null <- WSGmod_multicov(modelfile = "JAGS//WSG_plots_null.txt")
plotmods$nullmods$nullclus <- WSGmod_multicov(modelfile = "JAGS//WSG_plots_null_spatial.txt", spatscale = plotsF$ClusNum)
plotmods$nullmods$nullarea <- WSGmod_multicov(modelfile = "JAGS//WSG_plots_null_spatial.txt", spatscale = plotsF$AreaNum)
plotmods$climmods$clim <- WSGmod_multicov(covariates = plotsF[order(plotsF$SiteCode),c("ElevS","TotPrecS","TempVarS","PrecVarS","SlopeS")],
                                  modelfile = "JAGS//WSG_plots_cov.txt")
plotmods$climmods$climclus <- WSGmod_multicov(covariates = plotsF[order(plotsF$SiteCode),c("ElevS","TotPrecS","TempVarS","PrecVarS","SlopeS")],
                                     modelfile = "JAGS//WSG_plots_spatial_cov.txt", spatscale = plotsF$ClusNum)
start.time <- Sys.time()
plotmods$climmods$climarea <- WSGmod_multicov(covariates = plotsF[order(plotsF$SiteCode),c("ElevS","TotPrecS","TempVarS","PrecVarS","SlopeS")],
                                     modelfile = "JAGS//WSG_plots_spatial_cov.txt", spatscale = plotsF$AreaNum)
time.used <- Sys.time()-start.time
#saveRDS(plotmods, "JAGS\\WSGplotmods.RDS")
#plotmods <- readRDS("JAGS\\WSGplotmods.RDS")

## Summary table
plotmod_out <- lapply(plotmods$climmods, function(x) x$summary[1:8,c(1,3,7)])
for(i in 1:length(plotmod_out)){
  DIC <- plotmods$climmods[[i]]$DIC
  R2 <- apply(plotmods$climmods[[i]]$sims.list$mu, 1, function(x) cor(x, rowMeans(WSGdraws_short[,-1]))^2)
  RMSE <- apply(plotmods$climmods[[i]]$sims.list$MSE,1,function(x) sqrt(sum(x)/length(x)))
  #R2 <- apply(plotmods$climmods[[i]]$sims.list$mu, 1, function(x) cor(x, plotsF$WSGv)^2)
  #RMSE <- apply(plotmods$climmods[[i]]$sims.list$res,1,function(x) sqrt(sum(x^2)/length(x)))
  
  DICnull <- plotmods$nullmods[[i]]$DIC
  R2null <- apply(plotmods$nullmods[[i]]$sims.list$mu, 1, function(x) cor(x, rowMeans(WSGdraws_short[,-1]))^2)
  RMSEnull <- apply(plotmods$nullmods[[i]]$sims.list$MSE,1,function(x) sqrt(sum(x)/length(x)))
  #R2null <- apply(plotmods$nullmods[[i]]$sims.list$mu, 1, function(x) cor(x, plotsF$WSGv)^2)
  #RMSEnull <- apply(plotmods$nullmods[[i]]$sims.list$res,1,function(x) sqrt(sum(x^2)/length(x)))
  
  plotmod_out[[i]] <- rbind(plotmod_out[[i]], 
                     cbind(DIC, NA, NA),
                     cbind(mean(R2), R2[order(R2)][length(R2)*0.025], R2[order(R2)][length(R2)*0.975]),
                     cbind(mean(RMSE), RMSE[order(RMSE)][length(RMSE)*0.025], RMSE[order(RMSE)][length(RMSE)*0.975]), 
                     cbind(DICnull, NA, NA),
                     cbind(mean(R2null), R2null[order(R2null)][length(R2null)*0.025], R2null[order(R2null)][length(R2null)*0.975]),
                     cbind(mean(RMSEnull), RMSEnull[order(RMSEnull)][length(RMSEnull)*0.025], RMSEnull[order(RMSEnull)][length(RMSEnull)*0.975]))
  
  rownames(plotmod_out[[i]]) <- c("Intercept","Elevation","Total precipitation","Temperature variation","Precitipation variation","Slope","Residuals spatial","Residuals","DIC","R2","RMSE","DICnull","R2null","RMSEnull")
  }

plotmod_out <- do.call(cbind,plotmod_out)
#write.csv(round(plotmod_out,4),"Output\\WSG\\plotmod_out.csv")

## Range of WSG along environmental predictors
plotmods$climmods$clim$summary[2,1] * (max(range(plotsF$ElevS))-min(range(plotsF$ElevS))) ; print("WSGv - Elevation (nonspatial)")
plotmods$climmods$climclus$summary[2,1] * (max(range(plotsF$ElevS))-min(range(plotsF$ElevS))) ; print("WSGv - Elevation (cluster)")
plotmods$climmods$clim$summary[3,1] * (max(range(plotsF$TotPrecS))-min(range(plotsF$TotPrecS))) ; print("WSGv - TotPrec (nonspatial)")
plotmods$climmods$climclus$summary[3,1] * (max(range(plotsF$TotPrecS))-min(range(plotsF$TotPrecS))) ; print("WSGv - TotPrec (cluster)")

## Plot effects
MeanEnv <- colMeans(plotsF[,c("ElevS","TotPrecS","TempVarS","PrecVarS","SlopeS")])
meanb <- cbind(plotmods$climmods$clim$mean$b_cov * MeanEnv, plotmods$climmods$climclus$mean$b_cov * MeanEnv, plotmods$climmods$climarea$mean$b_cov * MeanEnv)
meanl <- cbind(plotmods$climmods$clim$q2.5$b_cov * MeanEnv, plotmods$climmods$climclus$q2.5$b_cov * MeanEnv, plotmods$climmods$climarea$q2.5$b_cov * MeanEnv)
meanu <- cbind(plotmods$climmods$clim$q97.5$b_cov * MeanEnv, plotmods$climmods$climclus$q97.5$b_cov * MeanEnv, plotmods$climmods$climarea$q97.5$b_cov * MeanEnv)

amean_elev <- sapply(plotmods$climmods, function(x) x$mean$alpha) + colSums(meanb[-1,])
al_elev <- sapply(plotmods$climmods, function(x) x$q2.5$alpha) + colSums(meanb[-1,])
au_elev <- sapply(plotmods$climmods, function(x) x$q97.5$alpha) + colSums(meanb[-1,])

plot(rowMeans(WSGdraws_short[,-1]) ~ plotsF$ElevS,
     #col = as.numeric(as.factor(plotsF$AreaNum)),
     xlab = "Elevation", ylab = "Volume-weighted WSG")
  polygon(x = c(0, 10, 10, 0), y = c(al_elev[1], al_elev[1] + plotmods$climmods$clim$q2.5$b_cov[1] * 10, au_elev[1] + plotmods$climmods$clim$q97.5$b_cov[1] * 10, au_elev[1]),
          border = "blue", lty = 3)
  polygon(x = c(0, 10, 10, 0), y = c(al_elev[2], al_elev[2] + plotmods$climmods$climclus$q2.5$b_cov[1] * 10, au_elev[2] + plotmods$climmods$climclus$q97.5$b_cov[1] * 10, au_elev[2]),
          border = "green", lty = 3)
  polygon(x = c(0, 10, 10, 0), y = c(al_elev[3], al_elev[3] + plotmods$climmods$climarea$q2.5$b_cov[1] * 10, au_elev[3] + plotmods$climmods$climarea$q97.5$b_cov[1] * 10, au_elev[3]),
          border = "red", lty = 3)

  abline(amean_elev[1], plotmods$climmods$clim$mean$b_cov[1], col = "blue")
  abline(amean_elev[2], plotmods$climmods$climclus$mean$b_cov[1], col = "green")
  abline(amean_elev[3], plotmods$climmods$climarea$mean$b_cov[1], col = "red")

amean_prec <- sapply(plotmods$climmods, function(x) x$mean$alpha) + colSums(meanb[-2,])
al_prec <- sapply(plotmods$climmods, function(x) x$q2.5$alpha) + colSums(meanb[-2,])
au_prec <- sapply(plotmods$climmods, function(x) x$q97.5$alpha) + colSums(meanb[-2,])
  
plot(rowMeans(WSGdraws_short[,-1]) ~ plotsF$TotPrecS,
     #col = as.numeric(as.factor(plotsF$AreaNum)),
     xlab = "Total annual precipitation (m)", ylab = "Volume-weighted WSG")
  polygon(x = c(0, 10, 10, 0), y = c(al_prec[1], al_prec[1] + plotmods$climmods$clim$q2.5$b_cov[2] * 10, au_prec[1] + plotmods$climmods$clim$q97.5$b_cov[2] * 10, au_prec[1]),
          border = "blue", lty = 3)
  polygon(x = c(0, 10, 10, 0), y = c(al_prec[2], al_prec[2] + plotmods$climmods$climclus$q2.5$b_cov[2] * 10, au_prec[2] + plotmods$climmods$climclus$q97.5$b_cov[2] * 10, au_prec[2]),
          border = "green", lty = 3)
  polygon(x = c(0, 10, 10, 0), y = c(al_prec[3], al_prec[3] + plotmods$climmods$climarea$q2.5$b_cov[2] * 10, au_prec[3] + plotmods$climmods$climarea$q97.5$b_cov[2] * 10, au_prec[3]),
          border = "red", lty = 3)
  
  abline(amean_prec[1], plotmods$climmods$clim$mean$b_cov[2], col = "blue")
  abline(amean_prec[2], plotmods$climmods$climclus$mean$b_cov[2], col = "green")
  abline(amean_prec[3], plotmods$climmods$climarea$mean$b_cov[2], col = "red")
  
##############################
## Spatial cross-validation ##
##############################
kfolds <- list("plot" = list(), "cluster" = list(), "area" = list(), "region" = list())

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
kfold <- list("spatialfold" = list(), "spatialfoldNull" = list(), "randomfold" = list(),"randomfoldNull" = list())

start.time <- Sys.time()

for(j in 1:4){
  ifelse(j %in% c(1,2), plotorder <- sortedplots, plotorder <- randomplots)
  ifelse(j %in% c(1,3), model.file <- "JAGS//WSG_plots_CV.txt", model.file <- "JAGS//WSG_plots_null_CV.txt")
  
  elpd <- list() ; MSE <- list(); WSGpred <- list(); sitenum <- list(); foldid <- list() ; mu <- list(); WSGpredM <- list(); WSGpred <- list() ; rhatmax <- list()
  for(i in 1:length(foldlength)){
    print(paste(ifelse(j %in% c(1,2), "Spatial","Random"),ifelse(j %in% c(1,3), " environmental"," null")," fold ",i," of ", nfolds,sep=""))
    
    pred <- plotorder[first[i]:last[i]]
    train <- plotorder[-which(plotorder %in% pred)]
    plotsF_p <- plotsF[pred,]
    plotsF_t <- plotsF[train,]
    WSGdraws_shortp <- WSGdraws_short[plotsF_p$SiteNum,-1]
    WSGdraws_shortt <- WSGdraws_short[plotsF_t$SiteNum,-1]
    
    moddat <- list("Nobst" = nrow(WSGdraws_shortt),
                   "WSGt" = WSGdraws_shortt,
                   "covt" = plotsF_t[,c("MeanTempS","TotPrecS","TempVarS","PrecVarS","SlopeS")],
                   "Nobsp" = nrow(WSGdraws_shortp),
                   "WSGp" = WSGdraws_shortp,
                   "ndraws" = ncol(WSGdraws_shortp),
                   "covp" = plotsF_p[,c("MeanTempS","TotPrecS","TempVarS","PrecVarS","SlopeS")],
                   "ncov" = 5)
    
    mod <- jagsUI::jags(model.file = model.file, data = moddat, inits = NULL,
                        parameters.to.save = c("alpha","var0","b_cov",
                                               "WSG.pred","elpd","MSE","mup"),
                        n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                        parallel = T, n.cores = n.cores)
    
    elpd[[i]] <- mod$sims.list$elpd
    MSE[[i]] <- mod$sims.list$MSE
    mu[[i]] <- mod$sims.list$mup
    WSGpred[[i]] <- mod$sims.list$WSG.pred
    sitenum[[i]] <- pred
    foldid[[i]] <- rep(i, length(pred))
    rhatmax[[i]] <- max(unlist(mod$Rhat))
  }
  
  kfold[[j]]$elpd <- do.call(cbind, elpd)
  kfold[[j]]$MSE <- do.call(cbind, MSE)
  kfold[[j]]$mu <- do.call(cbind, mu)
  kfold[[j]]$WSGpred <- do.call(cbind, WSGpred)
  kfold[[j]]$sitenum <- unlist(sitenum)
  kfold[[j]]$foldid <- unlist(foldid)
  kfold[[j]]$rhatmax <- unlist(rhatmax)
  
  kfold[[j]][c("elpd","MSE","mu","WSGpred")] <- lapply(kfold[[j]][c("elpd","MSE","mu","WSGpred")], function(x) x[,order(kfold[[j]]$sitenum)])
  kfold[[j]][c("sitenum","foldid","rhatmax")] <- lapply(kfold[[j]][c("sitenum","foldid","rhatmax")], function(x) x[order(kfold[[j]]$sitenum)])
  
  if(spat_cv == "plot"){if(j == 2){break}}
}
if(spat_cv == "region"){kfolds$region <- kfold}
if(spat_cv == "area"){kfolds$area <- kfold}
if(spat_cv == "cluster"){kfolds$cluster <- kfold}
if(spat_cv == "plot"){kfolds$plot <- kfold}

time.used[length(time.used) + 1] <- Sys.time() - start.time

# saveRDS(kfolds, file = "JAGS\\kfoldoutput.RDS")
# kfolds <- readRDS("JAGS\\kfoldoutput.RDS")

## Cross-validation mean table
elpd <- lapply(kfolds, function(x) lapply(x, function(x) apply(x$elpd, 1, sum)))
RMSE <- lapply(kfolds, function(x) lapply(x, function(x) apply(x$MSE, 1, function(y) sqrt(sum(y) / length(y)))))
R2 <- lapply(kfolds, function(x) lapply(x, function(x) apply(x$mu, 1, function(y) cor(y, rowMeans(WSGdraws_short)))))

elpdmean <- unlist(lapply(elpd, function(x) lapply(x, mean)))
RMSEmean <- unlist(lapply(RMSE, function(x) lapply(x, mean)))
R2mean <- unlist(lapply(R2, function(x) lapply(x, mean)))

CVtab <- rbind("leave-one-out" = c(elpdmean[1:4],MSEmean[1:4],R2mean[1:4]),
               "108-fold (cluster)" = c(elpdmean[5:8],MSEmean[5:8],R2mean[5:8]),
               "23-fold (area)" = c(elpdmean[9:12],MSEmean[9:12],R2mean[9:12]),
               "7-fold (region)" = c(elpdmean[13:16],MSEmean[13:16],R2mean[13:16]))
colnames(CVtab) <- c("elpdS","elpdSN","elpdR","elpdRN","RMSPES","RMSPESN","RMSPER","RMSPERN","pR2S","pR2SN","pR2R","pR2RN")

#write.csv(CVtab,"Output\\WSG\\plotmod_CVtab.csv")

##############################
## Predict plot-average WSG ## 
##############################
plotsF <- plotsF[order(plotsF$SiteNum),]
WSGdraws_plot <- WSGdraws[,-1]
WSGdraws_area <- lapply(split(WSGdraws_plot, plotsF$AreaNum), function(x) do.call("rbind",x))
WSGdraws_areareg <- lapply(split(t(do.call("cbind", lapply(WSGdraws_area, rowMeans))), as.data.frame(plotsF %>% group_by(AreaNum) %>% summarise(RegionNum = first(RegionNum)))$RegionNum), function(x) do.call("matrix",list(x, ncol = 4000)))
WSGdraws_region <- lapply(split(WSGdraws_plot, plotsF$RegionNum), function(x) do.call("rbind",x))

AGBdraws_plot <- WSGdraws_plot * plotsF$volha
AGBdraws_area <- lapply(split(AGBdraws_plot, plotsF$AreaNum), function(x) do.call("rbind",x))
AGBdraws_areareg <- lapply(split(t(do.call("cbind", lapply(AGBdraws_area, rowMeans))), as.data.frame(plotsF %>% group_by(AreaNum) %>% summarise(RegionNum = first(RegionNum)))$RegionNum), function(x) do.call("matrix",list(x, ncol = 4000)))
AGBdraws_region <- lapply(split(AGBdraws_plot, plotsF$RegionNum), function(x) do.call("rbind",x))

# AGB by region
test <- rbind(
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
    sort(apply(WSGdraws_plot, 2, function(x) sd(x) / mean(x)))[ncol(WSGdraws_plot)*0.975])))

write.csv(test,"Output\\WSG\\test.csv")

# Test CV between areas within region
test2 <- rbind(
  cbind(do.call("rbind",lapply(AGBdraws_region, function(x) mean(apply(x, 1, mean)))),
        do.call("rbind",lapply(AGBdraws_region, function(x) mean(apply(x, 1, function(y) sd(y)/mean(y))))),
        as.vector(as.data.frame(plotsF %>% group_by(RegionNum) %>% summarise(CVvol = sd(volha)/mean(volha)))$CVvol),
        do.call("rbind",lapply(WSGdraws_region, function(x) mean(apply(x, 1, function(y) sd(y)/mean(y))))),
        # Area level
        do.call("rbind",lapply(AGBdraws_areareg, function(x) mean(apply(x, 2, function(y) sd(y)/mean(y))))),
        as.vector(as.data.frame(plotsF %>% group_by(AreaNum) %>% summarise(vol = mean(volha), reg = first(RegionNum)) %>% group_by(reg) %>% summarise(CVvol = sd(vol)/mean(vol)))$CVvol),
        do.call("rbind",lapply(WSGdraws_areareg, function(x) mean(apply(x, 2, function(y) sd(y)/mean(y)))))),
        as.vector(c(
            mean(apply(AGBdraws_plot, 2, function(x) mean(x))),
            mean(apply(AGBdraws_plot, 2, function(x) sd(x) / mean(x))),
            sd(plotsF$volha) / mean(plotsF$volha),
            mean(apply(WSGdraws_plot, 2, function(x) sd(x) / mean(x))),
            sd(unlist(lapply(lapply(AGBdraws_area, rowMeans), mean))) / mean(unlist(lapply(lapply(AGBdraws_area, rowMeans), mean))),
            sd(as.data.frame(plotsF %>% group_by(AreaNum) %>% summarise(meanvol = mean(volha)))$meanvol) / mean(as.data.frame(plotsF %>% group_by(AreaNum) %>% summarise(meanvol = mean(volha)))$meanvol),
            sd(unlist(lapply(lapply(WSGdraws_area, rowMeans), mean))) / mean(unlist(lapply(lapply(WSGdraws_area, rowMeans), mean))))))

##############################
## Predict plot-average WSG ## 
##############################
## Predict values across dataset
treesF$WSGpred <- WSGtreemod$mean$WSG.pred
treesF$WSGused <- WSGtreemod$mean$WSG.pred
treesF$WSGused[which(!is.na(treesF$WSG))] <- treesF$WSG[which(!is.na(treesF$WSG))]

## Weigh WSG by volume
treesF$WSGv <- treesF$WSGused * treesF$VOLw

## This part adds volume-weighted WSG for cored trees only
treesFm <- filter(treesF, !is.na(WSG))
treesFm <- left_join(treesFm, treesFm %>% group_by(SiteCode) %>% summarise(totDBHm = sum(DBH_used), totVOLm = sum(vol)))
treesFm$VOLw <- treesFm$vol / treesFm$totVOLm
treesFm$WSGvm <- treesFm$WSG * treesFm$VOLw
treesF <- left_join(treesF, treesFm[,c("SiteCode","TreeN","StemN","totVOLm","WSGvm")])

## Calculate plot averages
#plotsF <- treesF %>% group_by(SiteCode) %>% summarise(WSGn = mean(WSG,na.rm=T), WSGp = mean(WSGused), WSGv = sum(WSGv), WSGvm = sum(WSGvm,na.rm=T))
plotsF <- left_join(plotsF, treesF %>% group_by(SiteCode) %>% summarise(WSGn = mean(WSG,na.rm=T), WSGp = mean(WSGused), WSGv = sum(WSGv), WSGvm = sum(WSGvm,na.rm=T)))

## Figures
temp <- plotsF[-which(plotsF$ncore == 0),]
cols <- RColorBrewer::brewer.pal(9, "Blues")[ceiling(plotsF$ncore/7) + 1]
plot((plotsF$WSGv-plotsF$WSGvm) ~ plotsF$ncore, ylab = "Estimated - measured plot average WSGv", xlab = "Number of cores") ; abline(h = 0)
plot(plotsF$WSGvm~plotsF$WSGv, col = cols, pch=16) ; abline(0,1)

###############################
## Compare biomass estimates ##
###############################
plotsF <- left_join(plotsF,plotdata[,c("SiteCode","Region")])

plotsF <- left_join(plotsF, as.data.frame(plotsF %>% group_by(AreaNum) %>% summarise(WSGareaMu = mean(WSGv))))
plotsF <- left_join(plotsF, as.data.frame(plotsF %>% group_by(RegionNum) %>% summarise(WSGregMu = mean(WSGv))))
plotsF <- mutate(plotsF, WSGallMu = mean(WSGv))

plotsF <- left_join(plotsF, CVout_region$CVoutput$spatialfold[,c("SiteNum","WSGpred")]) ; plotsF <- dplyr::rename(plotsF, c("WSGregEnv" = "WSGpred"))
plotsF <- left_join(plotsF, CVout_area$CVoutput$spatialfold[,c("SiteNum","WSGpred")]) ; plotsF <- dplyr::rename(plotsF, c("WSGareaEnv" = "WSGpred"))
plotsF <- left_join(plotsF, CVout_cluster$CVoutput$spatialfold[,c("SiteNum","WSGpred")]) ; plotsF <- dplyr::rename(plotsF, c("WSGclusEnv" = "WSGpred"))
plotsF <- left_join(plotsF, CVout_loo$CVoutput[,c("SiteNum","WSGpred")]) ; plotsF <- dplyr::rename(plotsF, c("WSGsiteEnv" = "WSGpred"))

plotsF$volha <- ((plotsF$vol * 1000) / plotsF$Size) * 10000     # Volume back to dm^3
plotsF$AGBind <- (plotsF$volha * plotsF$WSGv) / 1000

plotsF$AGBareaMu <- (plotsF$volha * plotsF$WSGareaMu) / 1000
plotsF$AGBregMu <- (plotsF$volha * plotsF$WSGregMu) / 1000
plotsF$AGBallMu <- (plotsF$volha * plotsF$WSGallMu) / 1000

plotsF$AGBregEnv <- (plotsF$volha * plotsF$WSGregEnv) / 1000
plotsF$AGBareaEnv <- (plotsF$volha * plotsF$WSGareaEnv) / 1000
plotsF$AGBclusEnv <- (plotsF$volha * plotsF$WSGclusEnv) / 1000
plotsF$AGBsiteEnv <- (plotsF$volha * plotsF$WSGsiteEnv) / 1000

AGBest <- plotsF[,c("AGBsiteEnv","AGBclusEnv","AGBareaEnv","AGBregEnv","AGBareaMu","AGBregMu","AGBallMu")]
AGBdev <- apply(AGBest, 2, function(x) x - plotsF$AGBind)

BMcomp <- data.frame("MAD_AGB" = apply(AGBdev, 2, function(x) median(abs(x))),
                     "MSD_AGB" = apply(AGBdev, 2, function(x) mean(x)),
                     "mMAD_plot" = apply(AGBdev, 2, function(x) max(tapply(x, plotsF$SiteNum, function(y) median(abs(y))))),
                     "mMSD_plot" = apply(apply(AGBdev, 2, function(x) range(tapply(x, plotsF$SiteNum, function(y) mean(y)))), 2, function(z) z[which(abs(z) == max(abs(z)))]),
                     "mMAD_cluster" = apply(AGBdev, 2, function(x) max(tapply(x, plotsF$ClusNum, function(y) median(abs(y))))),
                     "mMSD_cluster" = apply(apply(AGBdev, 2, function(x) range(tapply(x, plotsF$ClusNum, function(y) mean(y)))), 2, function(z) z[which(abs(z) == max(abs(z)))]),
                     "mMAD_area" = apply(AGBdev, 2, function(x) max(tapply(x, plotsF$AreaNum, function(y) median(abs(y))))),
                     "mMSD_area" = apply(apply(AGBdev, 2, function(x) range(tapply(x, plotsF$AreaNum, function(y) mean(y)))), 2, function(z) z[which(abs(z) == max(abs(z)))]),
                     "mMAD_region" = apply(AGBdev, 2, function(x) max(tapply(x, plotsF$RegionNum, function(y) median(abs(y))))),
                     "mMSD_region" = apply(apply(AGBdev, 2, function(x) range(tapply(x, plotsF$RegionNum, function(y) mean(y)))), 2, function(z) z[which(abs(z) == max(abs(z)))]))
row.names(BMcomp) <- c("Site","Cluster","Area","Region","Area mean", "Region mean", "Overall mean")

#write.csv(BMcomp,"Output//WSG//BMcomp.csv")

data.frame("MAD_AGB" = apply(AGBdev, 2, function(x) median(abs(x))),
           "MSD_AGB" = apply(AGBdev, 2, function(x) mean(x)),
           "MAD_cluster" = apply(AGBdev, 2, function(x) median(tapply(x, plotsF$ClusNum, function(y) median(abs(y))))),
           "MSD_cluster" = apply(AGBdev, 2, function(x) mean(tapply(x, plotsF$ClusNum, function(y) mean(y)))),
           "MAD_area" = apply(AGBdev, 2, function(x) median(tapply(x, plotsF$AreaNum, function(y) median(abs(y))))),
           "MSD_area" = apply(AGBdev, 2, function(x) mean(tapply(x, plotsF$AreaNum, function(y) mean(y)))),
           "MAD_region" = apply(AGBdev, 2, function(x) median(tapply(x, plotsF$RegionNum, function(y) median(abs(y))))),
           "MSD_region" = apply(AGBdev, 2, function(x) mean(tapply(x, plotsF$RegionNum, function(y) mean(y)))))

AGBesttab <- as.data.frame(plotsF %>% group_by(Region) %>% summarise(AGBregion = mean(AGBind)))
AGBesttab <- left_join(AGBesttab, as.data.frame(plotsF %>% group_by(Region,AreaNum) %>% summarise(AGBarea= mean(AGBind)) %>% group_by(Region) %>% summarise(AGBareaMin = min(AGBarea), AGBareaMax = max(AGBarea))))
AGBesttab <- left_join(AGBesttab, as.data.frame(plotsF %>% group_by(Region,ClusNum) %>% summarise(AGBcluster= mean(AGBind)) %>% group_by(Region) %>% summarise(AGBclusterMin = min(AGBcluster), AGBclusterMax = max(AGBcluster))))
AGBesttab <- left_join(AGBesttab, as.data.frame(plotsF %>% group_by(Region,SiteNum) %>% summarise(AGBsite= mean(AGBind)) %>% group_by(Region) %>% summarise(AGBplotMin = min(AGBsite), AGBplotMax = max(AGBsite))))

#write.csv(AGBesttab, "Output//WSG//AGBesttab.csv")









WSGest <- plotsF[,c("WSGsiteEnv","WSGclusEnv","WSGareaEnv","WSGregEnv","WSGareaMu","WSGregMu","WSGallMu")]
WSGdev <- apply(WSGest, 2, function(x) x - plotsF$WSGv)

BMcomp <- data.frame("MAD_WSG" = apply(WSGdev, 2, function(x) median(abs(x))),
                     "Bias_WSG" = apply(WSGdev, 2, function(x) mean(x)),
                     "MAD_AGB" = apply(AGBdev, 2, function(x) median(abs(x))),
                     "Bias_AGB" = apply(AGBdev, 2, function(x) mean(x)),
                     "MAD_Perc" = apply(AGBdev, 2, function(x) 100 * median(abs(x)) / median(plotsF$AGBind)),
                     "Bias_Perc" = apply(AGBdev, 2, function(x) 100 * mean(x) / mean(plotsF$AGBind)))

row.names(BMcomp) <- c("Site","Cluster","Area","Region","Area mean", "Region mean", "Overall mean")





#plotmods$nullunc <- WSGmod_multicov(modelfile = "JAGS//WSG_plots_null_unc.txt")

### OLD

## Retrieve DIC estimates
DIC_comp <- as.data.frame(matrix(NA, nrow = length(plotmods), ncol = 2))
for(i in 1:length(plotmods)){
  DIC_comp[i,1] <- plotmods[[i]]$DIC
  DIC_comp[i,2] <- plotmods[[i]]$pD
}
row.names(DIC_comp) <- names(plotmods) ; names(DIC_comp) <- c("DIC","pD") ; DIC_comp$model <- row.names(DIC_comp)

## Retrieve Bayesian R2 estimates
R2_comp <- as.data.frame(matrix(NA, nrow = length(plotmods), ncol = 4))
for(i in 1:length(plotmods)){
  R2 <- apply(plotmods[[i]]$sims.list$mu, 1, function(x) cor(x, plotsF$WSGv)^2)
  #var.mu <- apply(plotmods[[i]]$sims.list$mu,1,var)
  #var.res <- apply(plotmods[[i]]$sims.list$res,1,var)
  #R2 <- var.mu / (var.mu + var.res)
  R2_comp[i,1] <- mean(R2)
  R2_comp[i,2] <- sd(R2)
  R2_comp[i,3] <- R2[order(R2)][length(R2)*0.025]
  R2_comp[i,4] <- R2[order(R2)][length(R2)*0.975]
}
row.names(R2_comp) <- names(plotmods) ; names(R2_comp) <- c("R2","R2sd","R2q2.5","R2q97.5") ; R2_comp$model <- row.names(R2_comp)

## Retrieve RMSE estimates
RMSE_complist <- list()
for(j in 1:length(plotmods)){
  RMSE_comp <- as.data.frame(matrix(NA, nrow = length(plotmods), ncol = 4))
  for(i in 1:length(plotmods[[j]])){
    RMSE <- apply(plotmods$climmods[[i]]$sims.list$res,1,function(x) sqrt(sum(x^2)/length(x)))
    RMSE_comp[i,1] <- mean(RMSE)
    RMSE_comp[i,2] <- sd(RMSE)
    RMSE_comp[i,3] <- RMSE[order(RMSE)][length(RMSE)*0.025]
    RMSE_comp[i,4] <- RMSE[order(RMSE)][length(RMSE)*0.975]
  }
  row.names(RMSE_comp) <- names(plotmods[[j]]) ; names(RMSE_comp) <- c("mean","RMSEsd","2.5%","97.5%") ; RMSE_comp$model <- row.names(RMSE_comp)
  RMSE_complist[[j]] <- RMSE_comp
}

test <- t(WSGtreemod$sims.list$WSG.pred)
test2 <- matrix(treesF$WSG[which(!is.na(treesF$WSG))], nrow=length(treesF$WSG[which(!is.na(treesF$WSG))]),ncol=4000)
test[which(!is.na(treesF$WSG)),] <- test2
test <- apply(test, 2, function(x) x * treesF$VOLw)

test5 <- rowsum(test, group = treesF$SiteNum)

start.time <- Sys.time()
plotmods$null2 <- jagsUI::jags(model.file = "JAGS//WSG_plots_null.txt", data = list("Nobs" = length(plotsF$WSGv), "WSG" = test5, ndraws = ncol(test5)), 
                              inits = NULL, parameters.to.save = c("alpha","var0","WSG.pred","loglik","res"),
                              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                              parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res"))
time.used <- Sys.time()-start.time

plotmods$area <- jagsUI::jags(model.file = "JAGS//WSG_plots_area.txt", data = list("Nobs" = length(plotsF$WSGv), "WSG" = plotsF$WSGv, "AreaNum" = as.numeric(as.factor(plotsF$Cluster)), "narea" = length(unique(plotsF$Cluster))), 
                              inits = NULL, parameters.to.save = c("alpha","var0","vararea","b_cov","b_area","fit","fit.pred","WSG.pred","loglik","res"),
                              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                              parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res"))
plotmods$clim <- WSGmod_multicov(covariatesN = cbind(plotsF$MeanTemp, plotsF$MeanPrec, plotsF$TempVar,plotsF$PrecVar,plotsF$ALOSslope,plotsF$ALOSchili,plotsF$CWDmax))

plotmods$elev <- WSGmod_multicov(covariatesN = cbind(plotsF$elevs))
plotmods$MeanArid <- WSGmod_multicov(covariatesN = cbind(plotsF$arids))
plotmods$AridVar <- WSGmod_multicov(covariatesN = cbind(plotsF$aridvars))
plotmods$TempVar <- WSGmod_multicov(covariatesN = cbind(plotsF$tempvars))
plotmods$slope <- WSGmod_multicov(covariatesN = cbind(plotsF$slopes))
plotmods$aspect <- WSGmod_multicov(covariatesN = cbind(plotsF$aspects))
plotmods$tpi <- WSGmod_multicov(covariatesC = cbind(plotsF$tpi300c))
plotmods$soil <- WSGmod_multicov(covariatesC = cbind(plotsF$SoilType))
plotmods$clim <- WSGmod_multicov(covariatesN = cbind(plotsF$MeanTemp, plotsF$MeanPrec, plotsF$TempVar,plotsF$PrecVar,plotsF$ALOSslope,plotsF$ALOSchili,plotsF$CWDmax))
plotmods$local <- WSGmod_multicov(covariatesN = cbind(plotsF$elevs, plotsF$arids, plotsF$aridvars,plotsF$tempvars,plotsF$slopes,plotsF$aspects), covariatesC = cbind(plotsF$tpi300c, plotsF$SoilType))

plotmods$full <- WSGmod_multicov(covariatesN = cbind(plotsF$slopes,plotsF$aspects), covariatesC = cbind(plotsF$tpi300c))

test <- plotsF[,c("MeanTemp","TotPrec","PrecVar","TempVar","ALOSslope","PrecDryQ","CWDmax")]
PerformanceAnalytics::chart.Correlation(cbind(plotsF$WSGv,test), histogram=TRUE, pch=19)

wtest <- cbind(plotsF$MeanTemp,plotsF$MeanArid,plotsF$MeanPrec,plotsF$AridVar,plotsF$TempVar,plotsF$ALOSslope,cos(plotsF$ALOSaspect))
test2 <- prcomp(test, scale = T, center = T)
plotmods$PCA <- WSGmod_multicov(covariatesN = cbind(test2$x[,1], test2$x[,2], test2$x[,3], test2$x[,4], test2$x[,5],test2$x[,6]))

factoextra::fviz_eig(test2)
factoextra::fviz_pca_ind(test2)
factoextra::fviz_pca_var(test2)

## Retrieve loo estimates
loos <- list()
loomod <- plotmods$null;loos$null <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$area;loos$area <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$elev;loos$elev <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$MeanArid;loos$MeanArid <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$AridVar;loos$AridVar <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$TempVar;loos$TempVar <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$clim;loos$clim <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$slope;loos$slope <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$aspect;loos$aspect <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$local; loos$local <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loos_comp <- as.data.frame(loo_compare(loos)) ; loos_comp$model <- row.names(loos_comp)

loomod <- plotmods$clim;loos$clim <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$climclus;loos$climclus <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))

#loomod <- plotmods$PCA;loos$PCA <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))

## Retrieve DIC estimates
DIC_comp <- as.data.frame(matrix(NA, nrow = length(plotmods), ncol = 2))
for(i in 1:length(plotmods)){
  DIC_comp[i,1] <- plotmods[[i]]$DIC
  DIC_comp[i,2] <- plotmods[[i]]$pD
}
row.names(DIC_comp) <- names(plotmods) ; names(DIC_comp) <- c("DIC","pD") ; DIC_comp$model <- row.names(DIC_comp)

## Retrieve Bayesian R2 estimates
R2_comp <- as.data.frame(matrix(NA, nrow = length(plotmods), ncol = 4))
for(i in 1:length(plotmods)){
  var.mu <- apply(plotmods[[i]]$sims.list$mu,1,var)
  #var.res <- (plotmods[[i]]$sims.list$sigma)^2
  var.res <- apply(plotmods[[i]]$sims.list$res,1,var)
  R2 <- var.mu / (var.mu + var.res)
  R2_comp[i,1] <- mean(R2)
  R2_comp[i,2] <- sd(R2)
  R2_comp[i,3] <- R2[order(R2)][length(R2)*0.025]
  R2_comp[i,4] <- R2[order(R2)][length(R2)*0.975]
}
row.names(R2_comp) <- names(plotmods) ; names(R2_comp) <- c("R2","R2sd","R2q2.5","R2q97.5") ; R2_comp$model <- row.names(R2_comp)

apply(plotmods$clim$sims.list$WSG.pred,1,function(x) cor(x,plotsF$WSGv))

## Full comparison table
comptab <- left_join(left_join(DIC_comp, loos_comp), R2_comp)

## Summary table
plotmod_out <- round(plotmods$clim$summary[1:5,c(1,3,6)],3)
rownames(plotmod_out) <- c("Intercept","Elevation","Aridity","Aridity variability","Temperature variability")
plotmod_out2 <- round(plotmods$clim2$summary[c(1,3:7),c(1,3,6)],3)
rownames(plotmod_out2) <- c("Intercept","AreaVar","Elevation","Aridity","Aridity variability","Temperature variability")
plotmod_out3 <- round(plotmods$clim3$summary[c(1,3:7),c(1,3,6)],3)
rownames(plotmod_out3) <- c("Intercept","AreaVar","Elevation","Aridity","Aridity variability","Temperature variability")

RMSE <- sqrt(sum((plotmods$clim2$mean$res)^2)/length(plotmods$clim2$mean$res))
R2_bayes <- var(plotmods$clim3$mean$mu) / (var(plotmods$clim3$mean$mu) + var(plotmods$clim3$mean$res))
R2_cor <- (cor(plotmods$clim2$mean$WSG.pred, plotsF$WSGv))^2

#write.csv(plotmod_out,"Output\\WSG\\plotmod_out.csv")

## 
checkmod <- plotmods$climmods$clim
checkmod$summary[1:9,c(1:3,7:10)]
plot(plotsF$WSGv ~ plotsF$MeanTempS, col = as.numeric(as.factor(plotsF$AreaNum)))
intercept <- checkmod$summary[1,1] + checkmod$summary[which(rownames(checkmod$summary)=="b_cov[2]"),1] * mean(plotsF$arids) + 
             checkmod$summary[which(rownames(checkmod$summary)=="b_cov[3]"),1] * mean(plotsF$aridvars) + 
             checkmod$summary[which(rownames(checkmod$summary)=="b_cov[4]"),1] * mean(plotsF$tempvars)
abline(intercept,checkmod$summary[which(rownames(checkmod$summary)=="b_cov[1]"),1] ,col="black")

plot(plotsF$WSGv ~ plotsF$arids, col = as.numeric(as.factor(plotsF$AreaNum)))
intercept <- checkmod$summary[1,1] + checkmod$summary[which(rownames(checkmod$summary)=="b_cov[1]"),1] * mean(plotsF$elevs) + 
  checkmod$summary[which(rownames(checkmod$summary)=="b_cov[3]"),1] * mean(plotsF$aridvars) + 
  checkmod$summary[which(rownames(checkmod$summary)=="b_cov[4]"),1] * mean(plotsF$tempvars)
abline(intercept,checkmod$summary[which(rownames(checkmod$summary)=="b_cov[2]"),1] ,col="red")

plot(plotsF$WSGv ~ plotsF$aridvars, col = as.numeric(as.factor(plotsF$AreaNum)))
intercept <- checkmod$summary[1,1] + checkmod$summary[which(rownames(checkmod$summary)=="b_cov[1]"),1] * mean(plotsF$elevs) + 
  checkmod$summary[which(rownames(checkmod$summary)=="b_cov[2]"),1] * mean(plotsF$arids) + 
  checkmod$summary[which(rownames(checkmod$summary)=="b_cov[4]"),1] * mean(plotsF$tempvars)
abline(intercept,checkmod$summary[which(rownames(checkmod$summary)=="b_cov[3]"),1] ,col="blue")

plot(plotsF$WSGv ~ plotsF$tempvars, col = as.numeric(as.factor(plotsF$AreaNum)))
intercept <- checkmod$summary[1,1] + checkmod$summary[which(rownames(checkmod$summary)=="b_cov[1]"),1] * mean(plotsF$elevs) + 
  checkmod$summary[which(rownames(checkmod$summary)=="b_cov[2]"),1] * mean(plotsF$arids) + 
  checkmod$summary[which(rownames(checkmod$summary)=="b_cov[3]"),1] * mean(plotsF$aridvars)
abline(intercept,checkmod$summary[which(rownames(checkmod$summary)=="b_cov[4]"),1] ,col="blue")

## Posterior checks
y <- plotsF$WSGv
yrep <- plotmods$soil$sims.list$WSG.pred
bayesplot::ppc_dens_overlay(y = y, yrep = yrep[1:100, ])
bayesplot::ppc_hist(y, yrep[1:10, ])
bayesplot::ppc_stat(y, yrep, stat = "max")
bayesplot::ppc_stat(y, yrep, stat = "min")
bayesplot::ppc_stat(y, yrep, stat = "mean")
bayesplot::ppc_stat(y, yrep, stat = "median")
bayesplot::ppc_stat_grouped(y, yrep, group = plotsF$RegionNum, stat = "mean")

bayesplot::ppc_scatter(y, yrep[1:10, ])
bayesplot::ppc_scatter_avg(y, yrep)
bayesplot::ppc_scatter_avg_grouped(y, yrep, group = plotsF$RegionNum)

## Moran's I
plot.dists <- as.matrix(dist(cbind(plotsF$long, plotsF$lat)))
plot.dists <- plot.dists + 1
plot.dists.inv <- 1/plot.dists
diag(plot.dists.inv) <- 0
plot.dists.inv[is.infinite(plot.dists.inv)] <- 0

res_intercept <- plotmods$null$mean$res; MI_intercept <- as.data.frame(ape::Moran.I(res_intercept, plot.dists.inv))
res_area <- plotmods$area$mean$res; MI_area <- as.data.frame(ape::Moran.I(res_area, plot.dists.inv))
res_elev <- plotmods$elev$mean$res; MI_elev <- as.data.frame(ape::Moran.I(res_elev, plot.dists.inv))
res_tempvar <- plotmods$TempVar$mean$res; MI_tempvar <- as.data.frame(ape::Moran.I(res_tempvar, plot.dists.inv))
res_clim <- plotmods$clim$mean$res; MI_clim <- as.data.frame(ape::Moran.I(res_clim, plot.dists.inv))
cbind("Model" = c("intercept","area","elev","tempvar","clim"),
      rbind(MI_intercept,MI_area,MI_elev,MI_tempvar,MI_clim))

spatialEco::morans.plot(plotsF$WSGv,coords = cbind(plotsF$long,plotsF$lat))
spatialEco::morans.plot(plotmods$clim$mean$res,coords = cbind(plotsF$long,plotsF$lat))
spatialEco::morans.plot(plotmods$clim2$mean$res,coords = cbind(plotsF$long,plotsF$lat))
spatialEco::morans.plot(plotmods$clim3$mean$res,coords = cbind(plotsF$long,plotsF$lat))
spatialEco::morans.plot(plotmods$null$mean$res,coords = cbind(plotsF$long,plotsF$lat))
spatialEco::morans.plot(plotmods$clim4$mean$res,coords = cbind(plotsF$long,plotsF$lat))

res_PCA <- plotmods$PCA$mean$res; MI_clim <- as.data.frame(ape::Moran.I(res_PCA, plot.dists.inv))
