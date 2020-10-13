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

# Calculate volume
treesF <- dplyr::mutate(treesF, vol = exp(2.789 - 1.414 * log(DBH_used) + 1.178  * (log(DBH_used)^2) - 0.118 * (log(DBH_used)^3)))
treesF$vol <- treesF$vol / 1000     # From dm^3 to m^3 for readable effect sizes. Back-transform before AGB calculation

# Remove plots with no cores
treesF <- dplyr::filter(treesF, !SiteCode %in% c("C15_P2","PSF11","PSF12","TAF1","TAF2","TSF5"))

# Remove species with less than two individuals or no cores
coredSp <- treesF %>% mutate(cored = WSG/WSG)
coredSp <- coredSp %>% group_by(Species) %>% summarise(total = n(), cored = sum(cored,na.rm=T))
coredSp <- dplyr::filter(coredSp, cored > 1)
coredSp <- dplyr::filter(coredSp, total > cored)

# Add 
treesF[which(!treesF$Species %in% coredSp$Species),]$Species <- NA
treesF[which(is.na(treesF$Species)),]$Species <- "nospec"
treesF$speciesid <- ifelse(treesF$Species == "nospec", 0, 1)

# Add grouping levels as numbers
treesF$SiteNum <- as.numeric(as.factor(treesF$SiteCode))
treesF$ClusNum <- as.numeric(as.factor(treesF$Cluster))
treesF$AreaNum <- as.numeric(as.factor(treesF$AreaCode))
treesF$RegionNum <- as.numeric(as.factor(treesF$region))

# Plot distances
test <- sp::spDists(cbind(plotsF$long,plotsF$lat),longlat = TRUE)
test[lower.tri(test)] <- NA ; diag(test) <- NA

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
               #"Regionnum" = treesF$RegionNum,
               #"RegionNumA" = data.frame(treesF %>% group_by(AreaNum) %>% summarise(first(RegionNum)))[,2],
               "nsites" = length(unique(treesF$SiteCode)),
               "ncluster" = length(unique(treesF$Cluster)),
               "narea" = length(unique(treesF$AreaCode)),
               #"nregion" = length(unique(treesF$region)),
               "sitecores" = max(as.data.frame(treesF %>% group_by(SiteNum) %>% summarise(ncores = sum(!is.na(WSG))))$ncores)- as.data.frame(treesF %>% group_by(SiteNum) %>% summarise(ncores = sum(!is.na(WSG))))$ncores)

modinits <- function(){list(alpha = rnorm(1,0.5,0.1), muspec = rnorm(1,0,0.01), sdspec = rnorm(1,0.025,0.005))}

n.chains <- 4
n.burnin <- 20000
n.iter <- n.burnin + 10000
n.thin <- 10
n.cores <- n.chains

start.time <- Sys.time()
WSGtreemod <- jagsUI::jags(model.file = "JAGS//WSG_trees_varstruct.txt", data = moddat, inits = modinits,
                           parameters.to.save = c("alpha","var1","var0","b_vol","dsite","muspec","varspec","vararea","varsite","varcluster","RMSE","RMSE.pred","loglik","WSG.pred","res","b_area","b_cluster","b_site","mu"),
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                           parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res","mu"))
time.used <- Sys.time()-start.time

save(WSGtreemod, file = "JAGS\\WSGtreemod_out.Rdata")

## Model performance criteria
WSGtreeloo <- loo::loo(WSGtreemod$sims.list$loglik, r_eff = loo::relative_eff(exp(WSGtreemod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(WSGtreemod$sims.list$loglik) / n.chains)))))
WSGtreeDIC <- data.frame(cbind(WSGtreemod$DIC, WSGtreemod$pD)) ; colnames(WSGtreeDIC) <- c("DIC","pD")

## Posterior checks
y <- treesF$WSG[which(!is.na(treesF$WSG))]
yrep <- WSGtreemod$sims.list$WSG.pred[,which(!is.na(treesF$WSG))]
bayesplot::ppc_dens_overlay(y = y, yrep = yrep[1:100, ])
#bayesplot::ppc_dens(y, yrep)
bayesplot::ppc_hist(y, yrep[1:10, ])
bayesplot::ppc_stat(y, yrep, stat = "max")
bayesplot::ppc_stat(y, yrep, stat = "min")
bayesplot::ppc_stat(y, yrep, stat = "mean")
bayesplot::ppc_stat(y, yrep, stat = "median")
bayesplot::ppc_stat_grouped(y, yrep, group = treesF$speciesid[which(!is.na(treesF$WSG))], stat = "mean")

bayesplot::ppc_error_scatter(y, yrep)

bayesplot::ppc_error_scatter_avg(y, yrep)
bayesplot::ppc_scatter_avg(y, yrep)

## Scatter plot
plot(y~colMeans(yrep),
     xlab = "Estimated WSG", ylab = "Measured WSG")
abline(0,1)

## Whisker plot
b_area <- data.frame(AreaNum = 1:length(unique(treesF$AreaCode)), bArea = WSGtreemod$mean$b_area)
b_cluster <- data.frame(ClusNum = 1:length(unique(treesF$Cluster)), bCluster = WSGtreemod$mean$b_cluster, AreaNum = data.frame(treesF %>% group_by(ClusNum) %>% summarise(first(AreaNum)))[,2])
b_site <- data.frame(SiteNum = 1:length(unique(treesF$SiteCode)), bSite = WSGtreemod$mean$b_site, bSite95l = WSGtreemod$q2.5$b_site, bSite95u = WSGtreemod$q97.5$b_site, AreaNum = data.frame(treesF %>% group_by(SiteNum) %>% summarise(first(AreaNum)))[,2], ClusNum = data.frame(treesF %>% group_by(SiteNum) %>% summarise(first(ClusNum)))[,2], RegionNum = data.frame(treesF %>% group_by(SiteNum) %>% summarise(first(RegionNum)))[,2])
b_site2 <- left_join(b_site,b_cluster[,c("ClusNum","bCluster")])
b_site2 <- left_join(b_site2,b_area[,c("AreaNum","bArea")])
b_site2$bFull <- b_site2$bArea + b_site2$bCluster + b_site2$bSite

rect(152, -1, 250, 1, col="gray", lwd = 0, border = NA)

b_site2 <- left_join(b_site2, b_site2 %>% group_by(RegionNum) %>% summarise(RegOrd = mean(bSite)))
b_site2 <- left_join(b_site2, b_site2 %>% group_by(AreaNum) %>% summarise(AreaOrd = mean(bSite)))
b_site2 <- left_join(b_site2, b_site2 %>% group_by(ClusNum) %>% summarise(ClusOrd = mean(bSite)))
b_site2 <- b_site2[order(b_site2$RegOrd,b_site2$AreaOrd,b_site2$ClusOrd,b_site2$bSite),]
plot(NA, ylim = c(min(b_site2$bSite95l), max(b_site2$bSite95u)), xlim = c(0,350), xlab = "Site index", ylab = "Random effect")
#for(i in 1:nrow(b_site2)){rect(i-0.5, -1, i+0.5, 1, lwd = 0, border = NA, col = col2rgb(as.factor(b_site2$RegionNum)[i], alpha = TRUE))}
for(i in 1:nrow(b_site2)){segments(x0=i, y0=b_site2$bSite95l[i], x1=i, y1=b_site2$bSite95u[i])}
points(b_site2$bSite, pch = 19, col = b_site2$AreaNum)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
abline(h = 0)

test <- b_site %>% group_by(ClusNum) %>% filter(n() > 2) %>% summarise(diff = max(bSite) - min(bSite))
test2 <- b_cluster %>% group_by(AreaNum) %>% filter(n() > 2) %>% summarise(diff = max(bCluster) - min(bCluster))

## Add cluster/site ranges
## Model estimate table
treemod_esttab <- WSGtreemod$summary[c("alpha","b_vol","varspec","vararea","varcluster","varsite","dsite"),c(1:3,7)]
rownames(treemod_esttab) <- c("mean","volume","species","area","cluster","site","dsite")
treemod_esttab <- round(treemod_esttab,4)
#write.csv(treemod_esttab,"Output\\WSG\\treemod_esttab.csv")

## Bayesian R2
varres <- apply(WSGtreemod$sims.list$res[,which(!is.na(treesF$WSG))], 1, var)
varmu <- apply(WSGtreemod$sims.list$mu[,which(!is.na(treesF$WSG))], 1, var)
R2 <- varmu / (varmu + varres)
mean(R2)
R2[order(R2)][length(R2)*0.025]
R2[order(R2)][length(R2)*0.975]

## Residual variation
mean(varres)
varres[order(varres)][length(varres)*0.025]
varres[order(varres)][length(varres)*0.975]

## Moran's I
tree.dists <- as.matrix(dist(cbind(treesF$long, treesF$lat)))
tree.dists <- tree.dists + 1
tree.dists.inv <- 1/tree.dists
diag(tree.dists.inv) <- 0
tree.dists.inv[is.infinite(tree.dists.inv)] <- 0
ape::Moran.I(WSGtreemod$mean$res, tree.dists.inv)

##############################
## Predict plot-average WSG ## 
##############################
## Predict values across dataset
treesF$WSGpred <- WSGtreemod$mean$WSG.pred
treesF$WSGused <- WSGtreemod$mean$WSG.pred
treesF$WSGused[which(!is.na(treesF$WSG))] <- treesF$WSG[which(!is.na(treesF$WSG))]

## Weigh WSG by volume/DBH 
treesF <- left_join(treesF, treesF %>% group_by(SiteCode) %>% summarise(totDBH = sum(DBH_used), totVOL = sum(vol)))
treesF$VOLw <- treesF$vol / treesF$totVOL
treesF$DBHw <- treesF$DBH_used / treesF$totDBH
treesF$WSGv <- treesF$WSGused * treesF$VOLw
treesF$WSGd <- treesF$WSGused * treesF$DBHw

## This part adds volume-weighted WSG for cored trees only
treesFm <- filter(treesF, !is.na(WSG))
treesFm <- left_join(treesFm, treesFm %>% group_by(SiteCode) %>% summarise(totDBHm = sum(DBH_used), totVOLm = sum(vol)))
treesFm$VOLw <- treesFm$vol / treesFm$totVOLm
treesFm$DBHw <- treesFm$DBH_used / treesFm$totDBHm
treesFm$WSGvm <- treesFm$WSG * treesFm$VOLw
treesFm$WSGdm <- treesFm$WSG * treesFm$DBHw
treesF <- left_join(treesF, treesFm[,c("SiteCode","TreeN","StemN","totDBHm","totVOLm","WSGvm","WSGdm")])

## Calculate plot averages
plotsF <- treesF %>% group_by(SiteCode) %>% summarise(WSGn = mean(WSG,na.rm=T), WSGp = mean(WSGused), WSGv = sum(WSGv), WSGd = sum(WSGd), WSGvm = sum(WSGvm,na.rm=T), WSGdm = sum(WSGdm,na.rm=T), ntree = n(), ncore = sum(WSG/WSG,na.rm=T), vol = sum(vol))
plotsF <- as.data.frame(left_join(plotsF,Spatial))

## Figures
temp <- plotsF[-which(plotsF$ncore == 0),]
cols <- RColorBrewer::brewer.pal(9, "Blues")[ceiling(plotsF$ncore/7) + 1]
plot((plotsF$WSGv-plotsF$WSGvm) ~ plotsF$ncore, ylab = "Estimated - measured plot average WSGv", xlab = "Number of cores") ; abline(h = 0)
plot(plotsF$WSGvm~plotsF$WSGv, col = cols, pch=16) ; abline(0,1)



#############################
## Plot-average WSG models ##
#############################
# Scale covariates
plotsF$elevs <- plotsF$ALOSelev / 1000
plotsF$tempvars <- scale(plotsF$TempVar)
plotsF$arids <- plotsF$MeanArid
plotsF$aridvars <- scale(plotsF$AridVar)
plotsF$slopes <- scale(plotsF$ALOSslope)
plotsF$aspects <- scale(cos(plotsF$ALOSaspect))

plotsF$SiteNum <- as.numeric(as.factor(plotsF$SiteCode))
plotsF$ClusNum <- as.numeric(as.factor(plotsF$Cluster))
plotsF$AreaNum <- as.numeric(as.factor(plotsF$AreaCode))
plotsF$regionNum <- as.numeric(as.factor(plotsF$region))

## Function to run models
WSGmod_multicov <- function(covariatesN = NULL, covariatesC = NULL){
  moddat <- list("Nobs" = length(plotsF$WSGv),
                 "WSG" = plotsF$WSGv)
                 #"AreaNum" = plotsF$ClusNum,
                 #"narea" = length(unique(plotsF$Cluster)))
  
  moddat <- append(moddat, list("SiteNum" = plotsF$SiteNum,
                                "ClusNum" = plotsF$ClusNum,
                                "ClusNumS" = data.frame(plotsF %>% group_by(SiteNum) %>% summarise(first(ClusNum)))[,2],
                                "AreaNum" = plotsF$AreaNum,
                                "AreaNumC" = data.frame(plotsF %>% group_by(ClusNum) %>% summarise(first(AreaNum)))[,2],
                                "nsites" = length(unique(plotsF$SiteCode)),
                                "ncluster" = length(unique(plotsF$Cluster)),
                                "narea" = length(unique(plotsF$AreaCode))))
                              
  if(!is.null(covariatesC)){
    modelfile <- "JAGS//WSG_plots_covC.txt"
  }
  if(!is.null(covariatesN)){
    modelfile <- "JAGS//WSG_plots_covN.txt"
    if(!is.null(covariatesC)){
      modelfile <- "JAGS//WSG_plots_covNC.txt"
    }
  }
 
 #modelfile <-  "JAGS//WSG_plots_noarea.txt"
  modelfile <- "JAGS//WSG_plots_fullspatial.txt"
 
  if(!is.null(covariatesN)){
    moddat <- append(moddat, list("cov" = covariatesN,
                                  "ncov" = ncol(covariatesN)))}
  if(!is.null(covariatesC)){
    covCat <- lapply(as.data.frame(covariatesC), as.factor)
    covCat <- lapply(covCat[], as.numeric)
    covCat <- as.data.frame(covCat)
    moddat <- append(moddat, list("covC" = covCat,
                                  "ncovC" = ncol(covCat),
                                  "ncovCN" = as.vector(apply(as.data.frame(covCat),2, function(x) length(unique(x))))))}

  mod <- jagsUI::jags(model.file = modelfile, data = moddat, inits = NULL,
                      parameters.to.save = c("alpha","var0","vararea","varcovC","b_cov","b_covC","b_area","fit","fit.pred",
                                             "WSG.pred","loglik","res","mu"),
                      n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                      parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res"))
  return(mod)
}

## Run model
n.chains <- 4
n.burnin <- 10000
n.iter <- n.burnin + 1000
n.thin <- 1
n.cores <- n.chains

plotmods <- list()
plotmods$null <- jagsUI::jags(model.file = "JAGS//WSG_plots_null.txt", data = list("Nobs" = length(plotsF$WSGv), "WSG" = plotsF$WSGv), 
                              inits = NULL, parameters.to.save = c("alpha","var0","fit","fit.pred","WSG.pred","loglik","res"),
                              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                              parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res"))
plotmods$area <- jagsUI::jags(model.file = "JAGS//WSG_plots_area.txt", data = list("Nobs" = length(plotsF$WSGv), "WSG" = plotsF$WSGv, "AreaNum" = as.numeric(as.factor(plotsF$Cluster)), "narea" = length(unique(plotsF$Cluster))), 
                              inits = NULL, parameters.to.save = c("alpha","var0","vararea","b_cov","b_area","fit","fit.pred","WSG.pred","loglik","res"),
                              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                              parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res"))
plotmods$elev <- WSGmod_multicov(covariatesN = cbind(plotsF$elevs))
plotmods$MeanArid <- WSGmod_multicov(covariatesN = cbind(plotsF$arids))
plotmods$AridVar <- WSGmod_multicov(covariatesN = cbind(plotsF$aridvars))
plotmods$TempVar <- WSGmod_multicov(covariatesN = cbind(plotsF$tempvars))

plotmods$slope <- WSGmod_multicov(covariatesN = cbind(plotsF$slopes))
plotmods$aspect <- WSGmod_multicov(covariatesN = cbind(plotsF$aspects))
plotmods$tpi <- WSGmod_multicov(covariatesC = cbind(plotsF$tpi300c))

plotmods$soil <- WSGmod_multicov(covariatesC = cbind(plotsF$SoilType))

plotmods$clim4 <- WSGmod_multicov(covariatesN = cbind(plotsF$elevs, plotsF$arids, plotsF$aridvars,plotsF$tempvars))
plotmods$local <- WSGmod_multicov(covariatesN = cbind(plotsF$elevs, plotsF$arids, plotsF$aridvars,plotsF$tempvars,plotsF$slopes,plotsF$aspects), covariatesC = cbind(plotsF$tpi300c, plotsF$SoilType))

plotmods$full <- WSGmod_multicov(covariatesN = cbind(plotsF$slopes,plotsF$aspects), covariatesC = cbind(plotsF$tpi300c))

test <- cbind(plotsF$ALOSelev,plotsF$MeanArid,plotsF$AridVar,plotsF$TempVar) #,plotsF$ALOSslope,cos(plotsF$ALOSaspect))
test2 <- prcomp(test, scale = T, center = T)
plotmods$PCA <- WSGmod_multicov(covariatesN = cbind(test2$x[,1], test2$x[,2], test2$x[,3], test2$x[,4]))

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
plotmod_out3 <- round(plotmods$clim4$summary[c(1,3:7),c(1,3,6)],3)
rownames(plotmod_out2) <- c("Intercept","AreaVar","Elevation","Aridity","Aridity variability","Temperature variability")

RMSE <- sqrt(sum((plotmods$clim$mean$res)^2)/length(plotmods$clim$mean$res))
R2_bayes <- var(plotmods$clim$mean$mu) / (var(plotmods$clim$mean$mu) + var(plotmods$clim$mean$res))
R2_cor <- (cor(plotmods$clim$mean$WSG.pred, plotsF$WSGv))^2

#write.csv(plotmod_out,"Output\\WSG\\plotmod_out.csv")

## 
checkmod <- plotmods$clim
checkmod$summary[1:9,c(1:3,7:10)]
plot(plotsF$WSGv ~ plotsF$elevs, col = as.numeric(as.factor(plotsF$AreaCode)))
intercept <- checkmod$summary[1,1] + checkmod$summary[which(rownames(checkmod$summary)=="b_cov[2]"),1] * mean(plotsF$arids) + 
             checkmod$summary[which(rownames(checkmod$summary)=="b_cov[3]"),1] * mean(plotsF$aridvars) + 
             checkmod$summary[which(rownames(checkmod$summary)=="b_cov[4]"),1] * mean(plotsF$tempvars)
abline(intercept,checkmod$summary[which(rownames(checkmod$summary)=="b_cov[1]"),1] ,col="blue")

plot(plotsF$WSGv ~ plotsF$arids, col = as.numeric(as.factor(plotsF$AreaCode)))
intercept <- checkmod$summary[1,1] + checkmod$summary[which(rownames(checkmod$summary)=="b_cov[1]"),1] * mean(plotsF$elevs) + 
  checkmod$summary[which(rownames(checkmod$summary)=="b_cov[3]"),1] * mean(plotsF$aridvars) + 
  checkmod$summary[which(rownames(checkmod$summary)=="b_cov[4]"),1] * mean(plotsF$tempvars)
abline(intercept,checkmod$summary[which(rownames(checkmod$summary)=="b_cov[2]"),1] ,col="blue")

plot(plotsF$WSGv ~ plotsF$aridvars, col = as.numeric(as.factor(plotsF$AreaCode)))
intercept <- checkmod$summary[1,1] + checkmod$summary[which(rownames(checkmod$summary)=="b_cov[1]"),1] * mean(plotsF$elevs) + 
  checkmod$summary[which(rownames(checkmod$summary)=="b_cov[2]"),1] * mean(plotsF$arids) + 
  checkmod$summary[which(rownames(checkmod$summary)=="b_cov[4]"),1] * mean(plotsF$tempvars)
abline(intercept,checkmod$summary[which(rownames(checkmod$summary)=="b_cov[3]"),1] ,col="blue")

plot(plotsF$WSGv ~ plotsF$tempvars, col = as.numeric(as.factor(plotsF$AreaCode)))
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
bayesplot::ppc_stat_grouped(y, yrep, group = plotsF$region, stat = "mean")

bayesplot::ppc_scatter(y, yrep[1:10, ])
bayesplot::ppc_scatter_avg(y, yrep)
bayesplot::ppc_scatter_avg_grouped(y, yrep, group = plotsF$region)

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

## LOO CV (ELPPD, RMSE)
loodf <- data.frame(matrix(NA, nrow = length(plotsF$SiteCode), ncol = 4, dimnames = list(rownames=NULL,colnames = c("elpd","res","WSGpred","SiteNum"))))
loosummary <- data.frame(matrix(NA, nrow = 1, ncol = 3, dimnames = list(rownames = NULL, colnames = c("elpd","RMSE","R2"))))
elpd <- list() ; res <- list(); WSGpred <- list(); sitenum <- list() ; mu <- list()
for(i in 1:length(plotsF$SiteCode)){
  print(paste("loo fold ",i," of ", length(plotsF$SiteCode),sep=""))
  
  plotsF_p <- plotsF[i,]
  plotsF_t <- plotsF[-i,]
  
  moddat <- list("Nobst" = length(plotsF_t$WSGv),
                 "WSGt" = plotsF_t$WSGv,
                 "covt" = cbind(plotsF_t$elevs, plotsF_t$arids, plotsF_t$aridvars,plotsF_t$tempvars),
                 "Nobsp" = length(plotsF_p$WSGv),
                 "WSGp" = plotsF_p$WSGv,
                 "covp" = cbind(plotsF_p$elevs, plotsF_p$arids, plotsF_p$aridvars,plotsF_p$tempvars),
                 "ncov" = 4)
  
  mod <- jagsUI::jags(model.file = "JAGS//WSG_plots_noarea_CV.txt", data = moddat, inits = NULL,
                      parameters.to.save = c("alpha","var0","b_cov",
                                             "WSG.pred","elpd","res","mup"),
                      n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                      parallel = T, n.cores = n.cores)
  
  elpd[[i]] <- mod$mean$elpd
  res[[i]] <- mod$mean$res
  mu[[i]] <- mod$mean$mup
  WSGpred[[i]] <- mod$mean$WSG.pred
  sitenum[[i]] <- plotsF[i,]$SiteNum
}

loodf$elpd <- unlist(elpd) 
loodf$res <- unlist(res)
loodf$mu <- unlist(mu)
loodf$WSGpred <- unlist(WSGpred) 
loodf$SiteNum <- unlist(sitenum)

loosummary$elpd <- sum(unlist(loodf$elpd))
loosummary$RMSE <- sqrt(sum((loodf$res)^2)/length(loodf$res))
#loosummary$R2 <- var(loodf$mu) / (var(loodf$mu) + var(loodf$res))
loosummary$R2 <- cor(loodf[order(loodf$SiteNum),]$WSGpred,plotsF[order(plotsF$SiteNum),]$WSGv)

CVout_loo <- list("CVsummary" = loosummary, "CVoutput" = loodf)

## K-fold CV (ELPPD, RMSE)
spat_cv <- "area" # region, area, cluster

if(spat_cv == "region"){spatlevel <- plotsF$region} 
  if(spat_cv == "area"){spatlevel <- plotsF$AreaCode}
    if(spat_cv == "cluster"){spatlevel <- plotsF$Cluster}
nfolds <- length(unique(spatlevel))
randomplots <- sample(1:nrow(plotsF), nrow(plotsF), replace = FALSE, prob = NULL)
sortedplots <- plotsF[order(spatlevel),]$SiteNum
foldlength <- as.vector(table(spatlevel[order(spatlevel)]))
last <- cumsum(foldlength)
first <- last - foldlength + 1
kfolddf <- data.frame(matrix(NA,nrow = length(plotsF$SiteCode), ncol=6, dimnames = list(rownames=NULL,colnames = c("elpd","res","mu","WSGpred","SiteNum","foldid"))))
summarydf <- data.frame(matrix(NA, nrow = 1, ncol = 3, dimnames = list(rownames = NULL, colnames = c("elpd","RMSE","R2"))))
kfold <- list("spatialfold" = kfolddf, "randomfold" = kfolddf)
CVsummary <- list("spatialfold" = summarydf, "randomfold" = summarydf)

for(j in 1:2){
  ifelse(j == 1, plotorder <- sortedplots, plotorder <- randomplots)
  
  elpd <- list() ; res <- list(); WSGpred <- list(); sitenum <- list(); foldid <- list() ; mu <- list()
  for(i in 1:nfolds){
    print(paste(ifelse(j == 1, "Spatial","Random")," fold ",i," of ", nfolds,sep=""))
    
    pred <- plotorder[first[i]:last[i]]
    train <- plotorder[-which(plotorder %in% pred)]
    plotsF_p <- plotsF[pred,]
    plotsF_t <- plotsF[train,]
    
    moddat <- list("Nobst" = length(plotsF_t$WSGv),
                   "WSGt" = plotsF_t$WSGv,
                   "covt" = cbind(plotsF_t$elevs, plotsF_t$arids, plotsF_t$aridvars,plotsF_t$tempvars),
                   "Nobsp" = length(plotsF_p$WSGv),
                   "WSGp" = plotsF_p$WSGv,
                   "covp" = cbind(plotsF_p$elevs, plotsF_p$arids, plotsF_p$aridvars,plotsF_p$tempvars),
                   "ncov" = 4)
    
    mod <- jagsUI::jags(model.file = "JAGS//WSG_plots_noarea_CV.txt", data = moddat, inits = NULL,
                        parameters.to.save = c("alpha","var0","b_cov",
                                               "WSG.pred","elpd","res","mup"),
                        n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                        parallel = T, n.cores = n.cores)
    
    elpd[[i]] <- mod$mean$elpd
    res[[i]] <- mod$mean$res
    mu[[i]] <- mod$mean$mup
    WSGpred[[i]] <- mod$mean$WSG.pred
    sitenum[[i]] <- pred
    foldid[[i]] <- rep(i, length(pred))
  }
  
  kfold[[j]]$elpd <- unlist(elpd) 
  kfold[[j]]$res <- unlist(res)
  kfold[[j]]$mu <- unlist(mu)
  kfold[[j]]$WSGpred <- unlist(WSGpred) 
  kfold[[j]]$SiteNum <- unlist(sitenum) 
  kfold[[j]]$foldid <- unlist(foldid)
  
  kfold[[j]] <- kfold[[j]][order(kfold[[j]]$SiteNum),]
  
  CVsummary[[j]]$elpd <- sum(unlist(kfold[[j]]$elpd))
  CVsummary[[j]]$RMSE <- sqrt(sum((kfold[[j]]$res)^2)/length(kfold[[j]]$res)) 
  #CVsummary[[j]]$R2 <- var(kfold[[j]]$mu) / (var(kfold[[j]]$mu) + var(kfold[[j]]$res))
  CVsummary[[j]]$R2 <- (cor(kfold[[j]][order(kfold[[j]]$SiteNum),]$WSGpred,plotsF[order(plotsF$SiteNum),]$WSGv))^2
}
if(spat_cv == "region"){CVout_region <- list("CVsummary" = CVsummary, "CVoutput" = kfold)}
if(spat_cv == "area"){CVout_area <- list("CVsummary" = CVsummary, "CVoutput" = kfold)}
if(spat_cv == "cluster"){CVout_cluster <- list("CVsummary" = CVsummary, "CVoutput" = kfold)}

## Cross-validation table
clusfolds <- max(CVout_cluster$CVoutput$randomfold$foldid)
areafolds <- max(CVout_area$CVoutput$randomfold$foldid)
regfolds <- max(CVout_region$CVoutput$randomfold$foldid)
CVtab <- data.frame(cbind(rbind(CVout_loo$CVsummary, CVout_cluster$CVsummary$randomfold, CVout_area$CVsummary$randomfold, CVout_region$CVsummary$randomfold),
                          rbind(c("elpd" = NA, "RMSE" = NA, "R2" = NA), CVout_cluster$CVsummary$spatialfold, CVout_area$CVsummary$spatialfold, CVout_region$CVsummary$spatialfold)),
                    row.names = c("leave-one-out", paste(clusfolds, "-fold (cluster)", sep = ""), paste(areafolds , "-fold (area)", sep = ""), paste(regfolds , "-fold (region)", sep = "")))
write.csv(CVtab,"Output\\WSG\\plotmod_CVtab.csv")

###############################
## Compare biomass estimates ##
###############################
plotsF <- left_join(plotsF, as.data.frame(plotsF %>% group_by(AreaCode) %>% summarise(WSGareaMu = mean(WSGv))))
plotsF <- left_join(plotsF, as.data.frame(plotsF %>% group_by(region) %>% summarise(WSGregMu = mean(WSGv))))
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

mean(abs(plotsF$AGBind - plotsF$AGBareaMu)) ; mean(plotsF$AGBareaMu) - mean(plotsF$AGBind)
mean(abs(plotsF$AGBind - plotsF$AGBallMu)) ; mean(plotsF$AGBallMu) - mean(plotsF$AGBind)
mean(abs(plotsF$AGBind - plotsF$AGBareaEnv)) ; mean(plotsF$AGBareaEnv) - mean(plotsF$AGBind)
mean(abs(plotsF$AGBind - plotsF$AGBsiteEnv)) ; mean(plotsF$AGBsiteEnv) - mean(plotsF$AGBind)

BMcomp <- as.data.frame(cbind(rbind(mean(abs(plotsF$WSGv - plotsF$WSGsiteEnv)),
                                    mean(abs(plotsF$WSGv - plotsF$WSGclusEnv)),
                                    mean(abs(plotsF$WSGv - plotsF$WSGareaEnv)),
                                    mean(abs(plotsF$WSGv - plotsF$WSGregEnv)),
                                    mean(abs(plotsF$WSGv - plotsF$WSGareaMu)),
                                    mean(abs(plotsF$WSGv - plotsF$WSGregMu)),
                                    mean(abs(plotsF$WSGv - plotsF$WSGallMu))),
                              rbind(mean(plotsF$WSGsiteEnv) - mean(plotsF$WSGv),
                                    mean(plotsF$WSGclusEnv) - mean(plotsF$WSGv),
                                    mean(plotsF$WSGareaEnv) - mean(plotsF$WSGv),
                                    mean(plotsF$WSGregEnv) - mean(plotsF$WSGv),
                                    mean(plotsF$WSGareaMu) - mean(plotsF$WSGv),
                                    mean(plotsF$WSGregMu) - mean(plotsF$WSGv),
                                    mean(plotsF$WSGallMu) - mean(plotsF$WSGv)),
                              rbind(mean(abs(plotsF$AGBind - plotsF$AGBsiteEnv)),
                                    mean(abs(plotsF$AGBind - plotsF$AGBclusEnv)),
                                    mean(abs(plotsF$AGBind - plotsF$AGBareaEnv)),
                                    mean(abs(plotsF$AGBind - plotsF$AGBregEnv)),
                                    mean(abs(plotsF$AGBind - plotsF$AGBareaMu)),
                                    mean(abs(plotsF$AGBind - plotsF$AGBregMu)),
                                    mean(abs(plotsF$AGBind - plotsF$AGBallMu))),
                              rbind(mean(plotsF$AGBsiteEnv) - mean(plotsF$AGBind),
                                    mean(plotsF$AGBclusEnv) - mean(plotsF$AGBind),
                                    mean(plotsF$AGBareaEnv) - mean(plotsF$AGBind),
                                    mean(plotsF$AGBregEnv) - mean(plotsF$AGBind),
                                    mean(plotsF$AGBareaMu) - mean(plotsF$AGBind),
                                    mean(plotsF$AGBregMu) - mean(plotsF$AGBind),
                                    mean(plotsF$AGBallMu) - mean(plotsF$AGBind))),
                        row.names = c("Site","Cluster","Area","Region","Area mean", "Region mean", "Overall mean"))
names(BMcomp) <- c("RMSD","Bias","RMSD","Bias")
