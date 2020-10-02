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

# Scale covariates
#Plotdata$elevs <- Plotdata$ALOSelev / 1000
#Plotdata$tempvars <- scale(Plotdata$TempVar)
#Plotdata$arids <- Plotdata$MeanArid
#Plotdata$aridvars <- scale(Plotdata$AridVar)
#Plotdata$slopes <- scale(Plotdata$ALOSslope)
#Plotdata$aspects <- scale(cos(Plotdata$ALOSaspect))

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
               "nsites" = length(unique(treesF$SiteCode)),
               "ncluster" = length(unique(treesF$Cluster)),
               "narea" = length(unique(treesF$AreaCode)),
               "sitecores" = max(as.data.frame(treesF %>% group_by(SiteNum) %>% summarise(ncores = sum(!is.na(WSG))))$ncores)- as.data.frame(treesF %>% group_by(SiteNum) %>% summarise(ncores = sum(!is.na(WSG))))$ncores)

modinits <- function(){list(alpha = rnorm(1,0.5,0.1), muspec = rnorm(1,0,0.01), sdspec = rnorm(1,0.025,0.005))}

n.chains <- 4
n.burnin <- 10000
n.iter <- n.burnin + 2000
n.thin <- 1
n.cores <- n.chains

WSGtreemod <- jagsUI::jags(model.file = "JAGS//WSG_trees_varstruct.txt", data = moddat, inits = modinits,
                           parameters.to.save = c("alpha","var1","var0","b_vol","dsite","muspec","varspec","vararea","varsite","varcluster","fit","fit.pred","loglik","WSG.pred"),
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                           parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik"))

save(WSGtreemod, file = "JAGS\\WSGtreemod_out.Rdata")

## Model performance criteria
WSGtreeloo <- loo::loo(WSGtreemod$sims.list$loglik, r_eff = loo::relative_eff(exp(WSGtreemod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(WSGtreemod$sims.list$loglik) / n.chains)))))
WSGtreeDIC <- data.frame(cbind(WSGtreemod$DIC, WSGtreemod$pD)) ; colnames(WSGtreeDIC) <- c("DIC","pD")

## Model estimate table
treemod_esttab <- WSGtreemod$summary[c("alpha","b_vol","varspec","vararea","varcluster","varsite","dsite"),c(1:3,7)]
rownames(treemod_esttab) <- c("mean","volume","species","area","cluster","site")
treemod_esttab <- round(treemod_esttab,4)
#write.csv(treemod_esttab,"Output\\WSG\\treemod_esttab.csv")

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

bayesplot::ppc_scatter(y, yrep[1:10, ])
bayesplot::ppc_scatter_avg(y, yrep)
bayesplot::ppc_scatter_avg_grouped(y, yrep, group = treesF$speciesid[which(!is.na(treesF$WSG))])

## Moran's I
tree.dists <- as.matrix(dist(cbind(treesF$long, treesF$lat)))
tree.dists <- tree.dists + 1
tree.dists.inv <- 1/tree.dists
diag(tree.dists.inv) <- 0
tree.dists.inv[is.infinite(tree.dists.inv)] <- 0
ape::Moran.I(WSGtreemod$mean$res, plot.dists.inv)

#############################
## Plot-average WSG models ##  ADD CATEGORICAL VARIABLES
#############################
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

## Function to run models
WSGmod_multicov <- function(covariates){
  moddat <- list("Nobs" = length(plotsF$WSGv),
                 "WSG" = plotsF$WSGv,
                 "cov" = covariates,
                 "ncov" = ncol(covariates),
                 "AreaNum" = as.numeric(as.factor(plotsF$AreaCode)),
                 "narea" = length(unique(plotsF$AreaCode)))
  
  mod <- jagsUI::jags(model.file = "JAGS//WSG_plots_cov.txt", data = moddat, inits = NULL,
                      parameters.to.save = c("alpha","sigma","sdarea","b_cov","b_area","fit","fit.pred",
                                             "WSG.pred","loglik","res"),
                      n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                      parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res"))
  return(mod)
}

## Run model
plotmods <- list()
plotmods$null <- jagsUI::jags(model.file = "JAGS//WSG_plots_null.txt", data = list("Nobs" = length(plotsF$WSGv), "WSG" = plotsF$WSGv, "AreaNum" = as.numeric(as.factor(plotsF$AreaCode)), "narea" = length(unique(plotsF$AreaCode))), 
                              inits = NULL, parameters.to.save = c("alpha","sigma","sdarea","b_cov","b_area","fit","fit.pred","WSG.pred","loglik","res"),
                              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                              parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res"))
plotmods$elev <- WSGmod_multicov(covariates = cbind(plotsF$elevs))
plotmods$MeanArid <- WSGmod_multicov(covariates = cbind(plotsF$arids))
plotmods$AridVar <- WSGmod_multicov(covariates = cbind(plotsF$aridvars))
plotmods$TempVar <- WSGmod_multicov(covariates = cbind(plotsF$tempvars))
plotmods$clim <- WSGmod_multicov(covariates = cbind(plotsF$elevs, plotsF$arids, plotsF$aridvars,plotsF$tempvars))

# plotmods$slope <- WSGmod_multicov(covariates = cbind(plotsF$slopes))
# plotmods$aspect <- WSGmod_multicov(covariates = cbind(plotsF$aspects))
# plotmods$tpi <- WSGmod_multicov(covariates = cbind(plotsF$tpis))
# plotmods$local <- WSGmod_multicov(covariates = cbind(plotsF$slopes,plotsF$aspects,plotsF$tpis))

## Retrieve loo estimates
loos <- list()
loomod <- plotmods$null;loos$null <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$elev;loos$elev <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$MeanArid;loos$MeanArid <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$AridVar;loos$AridVar <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$TempVar;loos$TempVar <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loomod <- plotmods$clim;loos$clim <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
# loomod <- plotmods$slope;loos$slope <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
# loomod <- plotmods$aspect;loos$aspect <- loo::loo(loomod$sims.list$loglik, r_eff = loo::relative_eff(exp(loomod$sims.list$loglik), unlist(lapply(c(1:n.chains), function(x) rep(x, nrow(loomod$sims.list$loglik) / n.chains)))))
loos_comp <- as.data.frame(loo_compare(loos)) ; loos_comp$model <- row.names(loos_comp)

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
  var.fit <- apply(plotmods[[i]]$sims.list$WSG.pred,1,var)
  var.res <- (plotmods[[i]]$sims.list$sigma)^2
  R2 <- var.fit / (var.fit + var.res)
  R2_comp[i,1] <- mean(R2)
  R2_comp[i,2] <- sd(R2)
  R2_comp[i,3] <- sort(R2)[(length(R2)/100)*2.5]
  R2_comp[i,4] <- sort(R2)[(length(R2)/100)*97.5]
}
row.names(R2_comp) <- names(plotmods) ; names(R2_comp) <- c("R2","R2sd","R2q2.5","R2q97.5") ; R2_comp$model <- row.names(R2_comp)

## Full comparison table
comptab <- left_join(left_join(DIC_comp, loos_comp), R2_comp)

# write.csv(plotmods[[6]]$summary,"C:\\Users\\jorgesan\\Desktop\\clim.csv")

## 
checkmod <- plotmods$TempVar
checkmod$summary[1:8,c(1:3,7:10)]
plot(plotsF$WSGv ~ plotsF$tempvars, col = as.numeric(as.factor(plotsF$AreaCode)))
abline(checkmod$summary[1,1],checkmod$summary[4,1])

## Posterior checks
y <- plotsF$WSGv
yrep <- plotmods$clim$sims.list$WSG.pred
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
res_clim <- plotmods$elev$mean$res; MI_clim <- as.data.frame(ape::Moran.I(res_clim, plot.dists.inv))
cbind("Model" = c("intercept","clim"),
      rbind(MI_intercept,MI_clim))

############################# 
## Plot-average WSG models ## 
#############################
plotsF <- left_join(plotsF, as.data.frame(plotsF %>% group_by(AreaCode) %>% summarise(WSGarea = mean(WSGv))))
plotsF <- left_join(plotsF, as.data.frame(plotsF %>% group_by(region) %>% summarise(WSGregion = mean(WSGv))))
plotsF <- mutate(plotsF, WSGoverall = mean(WSGv))

plotsF$volha <- ((plotsF$vol * 1000) / plotsF$Size) * 10000     # Volume back to dm^3
plotsF$AGBind <- (plotsF$volha * plotsF$WSGv) / 1000
plotsF$AGBarea <- (plotsF$volha * plotsF$WSGarea) / 1000
plotsF$AGBregion <- (plotsF$volha * plotsF$WSGregion) / 1000
plotsF$AGBoverall <- (plotsF$volha * plotsF$WSGoverall) / 1000

plot(AGBind ~ AGBarea, plotsF) ; abline(0,1)
plot(AGBind ~ AGBregion, plotsF) ; abline(0,1)
plot(AGBind ~ AGBoverall, plotsF) ; abline(0,1)

# Function to run models
WSGmod_multicov <- function(covariates){
  moddat <- list("Nobs" = length(treesF$WSG),
                 "WSG" = treesF$WSG,
                 "cov" = covariates,
                 "ncov" = ncol(covariates),
                 "AreaNum" = as.numeric(as.factor(treesF$AreaCode)),
                 "narea" = length(unique(treesF$AreaCode)))
  
  mod <- jagsUI::jags(model.file = "JAGS//WSG_plots_cov.txt", data = moddat, inits = NULL,
                      parameters.to.save = c("alpha","sigma","sdarea","b_cov","b_area","fit","fit.pred",
                                             "WSG.pred","loglik","res"),
                      n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                      parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res"))
  return(mod)
}

testtreemodelev <- WSGmod_multicov(covariates = cbind(treesF$ALOSelev/1000))



