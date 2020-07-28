#################
################# Analysis of forest structure data
#################
setwd("C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD")

library(dplyr)
library(brms)

########################
## Read R environment ##
########################
load("Output\\models_AGB.Rdata")

###############
## Read data ##
###############
Plotdata <- read.csv("Data\\raw\\Plots_BM.csv",header=T,stringsAsFactors = F)
Treedata <- read.csv("Data\\raw\\Trees_BM.csv",header=T,stringsAsFactors = F)
Deaddata <- read.csv("Data\\raw\\Dead_BM.csv",header=T,stringsAsFactors = F)
PTrees <- read.csv("Data\\raw\\PTrees_BM.csv",header=T,stringsAsFactors = F)
Paramo <- read.csv("Data\\raw\\Paramo_BM.csv",header=T,stringsAsFactors = F)
Grass <- read.csv("Data\\raw\\Grass_BM.csv",header=T,stringsAsFactors = F)
Spatial <- read.csv("Data\\raw\\Spatialdata2.csv",header=T,stringsAsFactors = F)

Plotdata$AreaCode <- NULL;Plotdata$Cluster <- NULL
Plotdata <- dplyr::left_join(Plotdata,Spatial,by="SiteCode")
Plotdata$Cluster <- paste(Plotdata$AreaCode,Plotdata$Cluster,sep="")
Plotdata$sizeclass <- ifelse(Plotdata$Size > 100,1,0)

### Reduce dataset to forest points
plotF <- dplyr::filter(Plotdata,HabitatP=="Forest")
plotF <- dplyr::filter(plotF,AreaCode != "CC")
plotF <- dplyr::filter(plotF,AreaCode != "TU")
#plotF <- dplyr::filter(plotF,Dataset != "Chocó")
plotF <- dplyr::filter(plotF,is.na(FAge) | FAge > 20)

## Summarise structure variables
plotF <- dplyr::left_join(plotF,as.data.frame(dplyr::summarise(dplyr::group_by(Treedata,SiteCode),BA=sum(DBH,na.rm=TRUE))),)
plotF <- dplyr::left_join(plotF,as.data.frame(dplyr::summarise(dplyr::group_by(Treedata,SiteCode),WSG=mean(WSG,na.rm=TRUE))))
plotF <- dplyr::left_join(plotF,as.data.frame(dplyr::summarise(dplyr::group_by(Treedata,SiteCode),WSGused=mean(WSG_used,na.rm=TRUE))))
plotF <- dplyr::left_join(plotF,as.data.frame(dplyr::summarise(dplyr::group_by(dplyr::filter(Treedata,is.na(StemN) | StemN==1),SiteCode),Stems=dplyr::n())))
plotF <- dplyr::left_join(plotF,as.data.frame(dplyr::summarise(dplyr::group_by(dplyr::filter(Treedata,!is.na(WSG)),SiteCode),NWSG=dplyr::n())))
plotF <- dplyr::left_join(plotF,as.data.frame(dplyr::summarise(dplyr::group_by(Treedata,SiteCode),BAmax=max(DBH,na.rm=TRUE))),)

### Weighted average WSG
# Add plot average values for species
tree_wsg <- dplyr::filter(Treedata, !Species %in% c("Palm","Fern"));tree_wsg$Species <- as.character(tree_wsg$Species)
tree_wsg <- tree_wsg %>% group_by(SiteCode,Species) %>% mutate(WSG = ifelse(is.na(WSG),mean(WSG,na.rm=TRUE),WSG))
# Add plot average values (excluding oak/polylepis)
tree_wsg[which(is.na(tree_wsg$Species) | !tree_wsg$Species %in% c("Quercus Humboldtii","Polylepis")),]$Species <- "nospec"
tree_wsg <- tree_wsg %>% group_by(SiteCode,Species) %>% mutate(WSG = ifelse(is.na(WSG),mean(WSG,na.rm=TRUE),WSG))
# Add species average to plots without cores for Quercus Humboldtii and Polylepis
tree_wsg[which(tree_wsg$Species == "Quercus Humboldtii" & is.na(tree_wsg$WSG)),]$WSG <- mean(dplyr::filter(tree_wsg,Species=="Quercus Humboldtii")$WSG,na.rm=TRUE)
tree_wsg[which(tree_wsg$Species == "Polylepis" & is.na(tree_wsg$WSG)),]$WSG <- mean(dplyr::filter(tree_wsg,Species=="Polylepis")$WSG,na.rm=TRUE)
# Plot average weighted WSG
plotF <- dplyr::left_join(plotF,as.data.frame(dplyr::summarise(dplyr::group_by(tree_wsg,SiteCode),WSGw=mean(WSG,na.rm=TRUE))))
rm(tree_wsg)

plotF$Stems <- (plotF$Stems/plotF$Size)*10000
plotF$BA <- (plotF$BA/plotF$Size)*10000

## Total aboveground forest carbon
plotF <- dplyr::mutate(plotF,TotAGBha = TreeAGB + DeadAGB)
plotF <- dplyr::filter(plotF,!is.na(plotF$lat))

# Covariate correlations
test <- as.data.frame(cbind(TotAGBha=plotF$TotAGBha,Elevation=plotF$ALOSelev,MeanTemp=plotF$MeanTemp,TempMax=plotF$TempMax,TempMin=plotF$TempMin,TotPrec=plotF$TotPrec,PrecVar=plotF$PrecVar,TempVar=plotF$TempVar,PrecMax=plotF$PrecMax,PrecMin=plotF$PrecMin,CMI=plotF$CMI,MeanArid=plotF$MeanArid,AridMax=plotF$AridMax,AridMin=plotF$AridMin,AridN=plotF$AridN,AridVar=plotF$AridVar,Slope=plotF$ALOSslope,Aspect=plotF$ALOSaspect,SolRad=plotF$SolRad))
test <- as.data.frame(cbind(TotAGBha=plotF$TotAGBha,Elevation=plotF$ALOSelev,MeanTemp=plotF$MeanTemp,TempVar=plotF$TempVar,MeanArid=plotF$MeanArid,AridVar=plotF$AridVar,tpi300=plotF$tpi300,Slope=plotF$ALOSslope,Aspect=plotF$ALOSaspect,SolRad=plotF$SolRad))
test <- as.data.frame(cbind(TotAGBha=plotF$TotAGBha,Elevation=plotF$ALOSelev,MeanTemp=plotF$MeanTemp,TempVar=plotF$TempVar,MeanArid=plotF$MeanArid,AridVar=plotF$AridVar,tpi300=plotF$tpi300))
PerformanceAnalytics::chart.Correlation(test)

# Scale covariates
plotF$elevs <- scale(plotF$ALOSelev)
plotF$tempvars <- scale(plotF$TempVar)
plotF$arids <- scale(plotF$MeanArid)
plotF$aridvars <- scale(plotF$AridVar)
plotF$slopes <- scale(plotF$ALOSslope)
plotF$aspects <- scale(cos(plotF$ALOSaspect))

plotF$tpi300_sp[which(plotF$tpi300_sp %in% c("middle slope","upper slope", "lower slope"))] <- "slope"
plotF$tpi300_sp[which(plotF$tpi300_sp == "flats slope")] <- "flat"
plotF$tpis <- as.factor(plotF$tpi300_sp)

###########################################
## Run models - environmental covariates ##
###########################################
# Choose variable
svar <- "AGB"          #AGB, WSG, STEMS, BA, BAmax

# Define model terms
if(svar == "AGB"){
  plotF_svar <- dplyr::filter(plotF,!is.na(TotAGBha))
  names(plotF_svar)[which(names(plotF_svar) == "TotAGBha")] <- "svar"
  names(plotF_svar)[which(names(plotF_svar) == "sizeclass")] <- "sdvar"
  family <- brmsfamily("lognormal")
  nullstruct <- "(1|AreaCode)"
  fullformula <- "svar ~ elevs + aridvars + tempvars + arids"
  sparseformula <- "svar ~ elevs + tempvars"
}else if(svar == "WSG"){
  plotF_svar <- dplyr::filter(plotF,!is.na(WSGw))
  plotF_svar <- dplyr::filter(plotF_svar,NWSG > 0)
  names(plotF_svar)[which(names(plotF_svar) == "WSGw")] <- "svar"
  names(plotF_svar)[which(names(plotF_svar) == "NWSG")] <- "sdvar"
  family <- brmsfamily("gaussian")
  nullstruct <- "tpis + (1|AreaCode)"
  fullformula <- "svar ~ elevs + I(elevs^2) + aridvars + tempvars + arids"
  sparseformula <- "svar ~ elevs + I(elevs^2) + arids"
}else if(svar == "BA"){
  plotF_svar <- dplyr::filter(plotF,!is.na(BA))
  names(plotF_svar)[which(names(plotF_svar)=="BA")] <- "svar"
  names(plotF_svar)[which(names(plotF_svar) == "sizeclass")] <- "sdvar"
  family <- brmsfamily("gaussian")
  nullstruct <- "(1|AreaCode)"
  fullformula <- "svar ~ elevs + aridvars + tempvars + arids"
  sparseformula <- "svar ~ elevs + tempvars + arids"
}else if(svar == "STEMS"){
  plotF_svar <- dplyr::filter(plotF,!is.na(Stems))
  names(plotF_svar)[which(names(plotF_svar)=="Stems")] <- "svar"
  names(plotF_svar)[which(names(plotF_svar) == "sizeclass")] <- "sdvar"
  family <- brmsfamily("gaussian")
  nullstruct <- "(1|AreaCode)"
  fullformula <- "svar ~ elevs + aridvars + tempvars + arids"
  sparseformula <- "svar ~ elevs + arids + aridvars + tempvars"
}else if(svar == "BAmax"){
  plotF_svar <- dplyr::filter(plotF,!is.na(BAmax))
  names(plotF_svar)[which(names(plotF_svar)=="BAmax")] <- "svar"
  names(plotF_svar)[which(names(plotF_svar) == "sizeclass")] <- "sdvar"
  family <- brmsfamily("lognormal")
  nullstruct <- "(1|AreaCode)"
  fullformula <- "svar ~ elevs + aridvars + tempvars + arids"
  sparseformula <- "svar ~ elevs + arids"
}

# Other model specifications
chains <- 4
cores <- 4
warmup <- 5000
iter <- 10000 + warmup

# Fit models
models <- list()

# Null models
models$intercept_nos <- brms::brm(svar ~ 1,
                                  data = plotF_svar, family = family,
                                  chains = chains, cores = cores, warmup = warmup, iter = iter)
models$intercept <- brms::brm(bf(svar ~ 1,
                                 sigma ~ sdvar),
                              data = plotF_svar, family = family,
                              chains = chains, cores = cores, warmup = warmup, iter = iter)
models$area <- brms::brm(bf(svar ~ (1|AreaCode),
                            sigma ~ sdvar),
                         data = plotF_svar, family = family,
                         chains = chains, cores = cores, warmup = warmup, iter = iter)

# Local covariates
models$tpi <- brms::brm(bf(svar ~ tpis + (1|AreaCode),
                           sigma ~ sdvar),
                        data = plotF_svar, family = family,
                        chains = 4, cores = 4)
models$slope <- brms::brm(bf(svar ~ slopes + (1|AreaCode),
                             sigma ~ sdvar),
                          data = plotF_svar, family = family,
                          chains = 4, cores = 4)
models$aspect <- brms::brm(bf(svar ~ aspects + (1|AreaCode),
                              sigma ~ sdvar),
                           data = plotF_svar, family = family,
                           chains = 4, cores = 4)

# Single parameter models
models$elev <- brms::brm(bf(paste("svar ~ elevs + ", nullstruct, sep=""),
                            sigma ~ sdvar),
                         data = plotF_svar, family = family,
                         chains = 4, cores = 4)
models$elev2 <- brms::brm(bf(paste("svar ~ elevs + I(elevs^2) + ", nullstruct, sep=""),
                             sigma ~ sdvar),
                          data = plotF_svar, family = family,
                          chains = 4, cores = 4)
models$tempvar <- brms::brm(bf(paste("svar ~ tempvars + ", nullstruct, sep=""),
                            sigma ~ sdvar),
                         data = plotF_svar, family = family,
                         chains = 4, cores = 4)
models$tempvar2 <- brms::brm(bf(paste("svar ~ tempvars + I(tempvars^2) + ", nullstruct, sep=""),
                             sigma ~ sdvar),
                          data = plotF_svar, family = family,
                          chains = 4, cores = 4)
models$meanarid <- brms::brm(bf(paste("svar ~ arids + ", nullstruct, sep=""),
                               sigma ~ sdvar),
                            data = plotF_svar, family = family,
                            chains = 4, cores = 4)
models$meanarid2 <- brms::brm(bf(paste("svar ~ arids + I(arids^2) + ", nullstruct, sep=""),
                                sigma ~ sdvar),
                             data = plotF_svar, family = family,
                             chains = 4, cores = 4)
models$aridvar <- brms::brm(bf(paste("svar ~ aridvars + ", nullstruct, sep=""),
                               sigma ~ sdvar),
                            data = plotF_svar, family = family,
                            chains = 4, cores = 4)
models$aridvar2 <- brms::brm(bf(paste("svar ~ aridvars + I(aridvars^2) + ", nullstruct, sep=""),
                                sigma ~ sdvar),
                             data = plotF_svar, family = family,
                             chains = 4, cores = 4)

# Full models
models$full <- brms::brm(bf(paste(fullformula, " + ",nullstruct, sep = ""),
                            sigma ~ sdvar),
                         data = plotF_svar, family = family,
                         chains = 4, cores = 4)
models$sparse <- brms::brm(bf(paste(sparseformula, " + ", nullstruct, sep = ""),
                              sigma ~ sdvar),
                           data = plotF_svar, family = family,
                           chains = 4, cores = 4)

# Estimate loo
loos <- list()
loos$intercept_nos <- loo::loo(models$intercept_nos, reloo=T)
loos$intercept <- loo::loo(models$intercept, reloo=T)
loos$area <- loo::loo(models$area, reloo=T)
loos$tpi <- loo::loo(models$tpi, reloo=T)
loos$slope <- loo::loo(models$slope, reloo=T)
loos$aspect <- loo::loo(models$aspect, reloo=T)

loos$elev <- loo::loo(models$elev, reloo=T)
loos$elev2 <- loo::loo(models$elev2, reloo=T)
loos$tempvar <- loo::loo(models$tempvar, reloo=T)
loos$tempvar2 <- loo::loo(models$tempvar2, reloo=T)
loos$meanarid <- loo::loo(models$meanarid, reloo=T)
loos$meanarid2 <- loo::loo(models$meanarid2, reloo=T)
loos$aridvar <- loo::loo(models$aridvar, reloo=T)
loos$aridvar2 <- loo::loo(models$aridvar2, reloo=T)

loos$full <- loo::loo(models$full, reloo=T)
loos$sparse <- loo::loo(models$sparse, reloo=T)

# Bayesian R^2
rsq <- list()
rsq$intercept_nos <- brms::bayes_R2(models$intercept_nos)
rsq$intercept <- brms::bayes_R2(models$intercept)
rsq$area <- brms::bayes_R2(models$area)
rsq$tpi <- brms::bayes_R2(models$tpi)
rsq$slope <- brms::bayes_R2(models$slope)
rsq$aspect <- brms::bayes_R2(models$aspect)

rsq$elev <- brms::bayes_R2(models$elev)
rsq$elev2 <- brms::bayes_R2(models$elev)
rsq$tempvar <- brms::bayes_R2(models$tempvar)
rsq$tempvar2 <- brms::bayes_R2(models$tempvar2)
rsq$meanarid <- brms::bayes_R2(models$meanarid)
rsq$meanarid2 <- brms::bayes_R2(models$meanarid2)
rsq$aridvar <- brms::bayes_R2(models$aridvar)
rsq$aridvar2 <- brms::bayes_R2(models$aridvar2)

rsq$full <- brms::bayes_R2(models$full)
rsq$sparse <- brms::bayes_R2(models$sparse)

# Bayesian R^2
loo_rsq <- list()
loo_rsq$intercept_nos <- brms::loo_R2(models$intercept_nos, reloo=T)
loo_rsq$intercept <- brms::loo_R2(models$intercept, reloo=T)
loo_rsq$area <- brms::loo_R2(models$area, reloo=T)
loo_rsq$tpi <- brms::loo_R2(models$tpi, reloo=T)
loo_rsq$slope <- brms::loo_R2(models$slope, reloo=T)
loo_rsq$aspect <- brms::loo_R2(models$aspect, reloo=T)

loo_rsq$elev <- brms::loo_R2(models$elev, reloo=T)
loo_rsq$elev2 <- brms::loo_R2(models$elev2, reloo=T)
loo_rsq$tempvar <- brms::loo_R2(models$tempvar, reloo=T)
loo_rsq$tempvar2 <- brms::loo_R2(models$tempvar2, reloo=T)
loo_rsq$meanarid <- brms::loo_R2(models$meanarid, reloo=T)
loo_rsq$meanarid2 <- brms::loo_R2(models$meanarid2, reloo=T)
loo_rsq$aridvar <- brms::loo_R2(models$aridvar, reloo=T)
loo_rsq$aridvar2 <- brms::loo_R2(models$aridvar2, reloo=T)

loo_rsq$full <- brms::loo_R2(models$full, reloo=T)
loo_rsq$sparse <- brms::loo_R2(models$sparse, reloo=T)

# Model comparisons
loo_compare <- loo::loo_compare(loos)

rsq_compare <- sapply(rsq, function(x) x[,1:2])
rsq_compare <- t(round(rsq_compare[order(rsq_compare[,1],decreasing=T),],3))

loo_rsq_compare <- sapply(loo_rsq, function(x) x[1])
loo_rsq_compare <- round(loo_rsq_compare[order(loo_rsq_compare,decreasing=T)],3)
names(loo_rsq_compare) <- stringr::str_remove(names(loo_rsq_compare),".R2")

## Moran's I
plot.dists <- as.matrix(dist(cbind(plotF_svar$long, plotF_svar$lat)))
plot.dists.inv <- 1/plot.dists
diag(plot.dists.inv) <- 0
res_intercept <- brms::predictive_error(models$intercept) ; MI_intercept <- as.data.frame(ape::Moran.I(colMeans(res_intercept),plot.dists.inv))
res_area <- brms::predictive_error(models$area) ; MI_area <- as.data.frame(ape::Moran.I(colMeans(res_area),plot.dists.inv))

##################
## Model checks ##
##################
modcheck <- models$full
parcheck <- names(modcheck$fit)[startsWith(names(modcheck$fit),"b")]
parcheck2 <- parcheck[-which(endsWith(parcheck,"t"))]
parcheck2 <- parcheck2[-which(parcheck2 == "b_sigma_sdvar")]

yrep <- brms::posterior_predict(modcheck)
y <- plotF_svar$svar
group <- plotF_svar$AreaCode

# launch_shinystan(modcheck)       # Launches browser window with diagnostics plots

# Random effects
randcheck <- names(modcheck$fit)[startsWith(names(modcheck$fit),"r")]
#randcheck2 <- randcheck[-which(endsWith(randcheck,"t]"))]
bayesplot::mcmc_intervals(modcheck, pars = randcheck, prob = 0.75)

# Posterior predictions
bayesplot::ppc_dens_overlay(y,yrep[1:50,])
bayesplot::ppc_hist(y,yrep[1:5,])
bayesplot::ppc_stat_grouped(y = y, yrep = yrep, group = group, stat = "median")
bayesplot::ppc_intervals(y = y, yrep = yrep, x = plotF_svar$elevs[,1])

bayesplot::ppc_stat(y = y, yrep = yrep, stat = "median")
bayesplot::ppc_stat_grouped(y = y, yrep = yrep, group = group, stat = "median")

plot(brms::conditional_effects(modcheck, "elevs"), points = T)

# Parameter estimates
bayesplot::mcmc_areas(modcheck, pars = parcheck, prob = 0.75)
bayesplot::mcmc_intervals(modcheck, pars = parcheck, prob = 0.75)
bayesplot::mcmc_hist(modcheck, pars = parcheck)
bayesplot::mcmc_dens_overlay(modcheck, pars = parcheck)
bayesplot::mcmc_violin(modcheck, pars = parcheck)
bayesplot::mcmc_trace(modcheck, pars = parcheck)
bayesplot::mcmc_trace_highlight(modcheck, pars = parcheck, highlight = 4)
bayesplot::mcmc_pairs(modcheck, pars = parcheck)
bayesplot::mcmc_acf(modcheck, pars = parcheck)

# NUTS diagnostics
lp <- bayesplot::log_posterior(modcheck)
nutspar <- bayesplot::nuts_params(modcheck)
rhat <- bayesplot::rhat(modcheck)
neff <- bayesplot::neff_ratio(modcheck)
pardraw <- as.array(modcheck)

range(rhat)
range(neff)

bayesplot::mcmc_parcoord(pardraw, np = nutspar, pars = parcheck)
bayesplot::mcmc_pairs(pardraw, np = nutspar)
bayesplot::mcmc_trace(pardraw, np = nutspar)

bayesplot::mcmc_nuts_acceptance(nutspar,lp)
bayesplot::mcmc_nuts_divergence(nutspar,lp)
bayesplot::mcmc_nuts_stepsize(nutspar,lp)
bayesplot::mcmc_nuts_energy(nutspar,lp)
bayesplot::mcmc_nuts_treedepth(nutspar,lp)

bayesplot::mcmc_rhat(rhat) + bayesplot::yaxis_text(hjust=0)
bayesplot::mcmc_neff(neff) + bayesplot::yaxis_text(hjust=0)

###############
##  OUTPUTS  ##
###############
## Model comparison table
comparetab <- merge(loo_rsq_compare, rsq_compare, by = "row.names") ; row.names(comparetab) <- comparetab$Row.names ; comparetab$Row.names <- NULL
comparetab <- merge(comparetab, loo_compare, by = "row.names")
comparetab <- comparetab[c("Row.names","x","Estimate","Est.Error","elpd_loo","se_elpd_loo","elpd_diff","se_diff")]
modorder <- c("intercept_nos","intercept","area","slope","aspect","tpi","elev","elev2","tempvar","tempvar2","meanarid","meanarid2","aridvar","aridvar2","full","sparse")
comparetab <- comparetab[order(match(comparetab$Row.names,modorder)),]
names(comparetab) <- c("model","RsqAdj","Rsq","RsqErr","elpd","elpdSD","elpdDiff","elpdDiffSD")

## Moran's I table
MItab <- data.frame(cbind(MI_intercept, MI_area), row.names = svar)

###############
#### SAVE  ####
###############
## Save R environment
save.image(file = paste("Output\\models_", svar, ".Rdata", sep = ""))
rm(list=ls())

## Save tables
write.csv(comparetab, file = paste("Output\\mod_compare_", svar, ".csv", sep = ""), row.names = F)
write.table(MItab, file = "Output\\MoransI.csv", sep = ",", col.names = NA, row.names = T, append = T)

# Summary table
if(svar == "AGB"){
  AGB_pred <- cbind(plotF_svar, "AGB_pred" = exp(colMeans(log(brms::posterior_predict(models$full)))))
  sum_tab <- AGB_pred %>% group_by(AreaCode) %>% summarise(Lat = median(lat, na.rm = T),
                                                        Long = median(long, na.rm = T),
                                                        AGB = median(svar, na.rm = T),
                                                        AGBpred = median(AGB_pred, na.rm = T),
                                                        WSG = median(WSG, na.rm = T),
                                                        BA = median(BA, na.rm = T),
                                                        Stems = median(Stems, na.rm = T),
                                                        Elevation = median(ALOSelev, na.rm = T),
                                                        Aridity = median(MeanArid, na.rm = T),
                                                        TempVar = median(TempVar, na.rm = T),
                                                        AridVar = median(AridVar, na.rm = T))
  write.csv(sum_tab,"Output\\SummaryTable.csv", row.names = F)
}

## Save plots
parcheck <- names(models$full$fit)[startsWith(names(models$full$fit),"b")]
parcheck2 <- parcheck[-which(endsWith(parcheck,"t"))]
parcheck2 <- parcheck2[-which(parcheck2 == "b_sigma_sdvar")]
covars <- sub("b_","",parcheck2)
covars <- covars[-grep(c("tpisslope"), covars)]
covars <- covars[-grep(c("tpisridge"), covars)]
covars <- covars[-grep(c("IelevsE2"), covars)]
covars <- sub("valley","",covars)

# Posterior intervals
png(file = paste("Output\\Intervals_", svar, ".png", sep = ""), width = 1024, height = 128)
bayesplot::mcmc_intervals(models$full, pars = parcheck2, prob = 0.75) + ggplot2::ggtitle(svar)
dev.off()

# Conditional effects
for(i in 1:length(covars)){
  png(file = paste("Output\\Effect_", svar, "_", covars[i], ".png",sep = ""), width = 1024, height = 512)
  plot(brms::conditional_effects(models$full,covars[i]),points=T)
  dev.off()
}





######################################
## Run models - structure variables ##
######################################
# Define model terms
plotF_AGB <- dplyr::filter(plotF,!is.na(TotAGBha))
plotF_AGB <- dplyr::filter(plotF_AGB,!is.na(WSG))

# Other model specifications
chains <- 4
cores <- 4
warmup <- 5000
iter <- 10000 + warmup

AGBmods <- list()
AGBmods$intercept_nos <- brms::brm(bf(TotAGBha ~ 1),
                         data = plotF_AGB, family = "lognormal",
                         chains = chains, cores = cores, warmup = warmup, iter = iter)
AGBmods$intercept <- brms::brm(bf(TotAGBha ~ 1,
                            sigma ~ sizeclass),
                         data = plotF_AGB, family = "lognormal",
                         chains = chains, cores = cores, warmup = warmup, iter = iter)
AGBmods$area <- brms::brm(bf(TotAGBha ~ (1 | AreaCode),
                            sigma ~ sizeclass),
                         data = plotF_AGB, family = "lognormal",
                         chains = chains, cores = cores, warmup = warmup, iter = iter)
AGBmods$WSG <- brms::brm(bf(TotAGBha ~ WSG + (WSG | AreaCode),
                           sigma ~ sizeclass),
                           data = plotF_AGB, family = "lognormal",
                           chains = chains, cores = cores, warmup = warmup, iter = iter)
AGBmods$BA <- brms::brm(bf(TotAGBha ~ BA + (BA | AreaCode),
                            sigma ~ sizeclass),
                         data = plotF_AGB, family = "lognormal",
                         chains = chains, cores = cores, warmup = warmup, iter = iter)
AGBmods$BAmax <- brms::brm(bf(TotAGBha ~ BAmax + (BAmax | AreaCode),
                            sigma ~ sizeclass),
                         data = plotF_AGB, family = "lognormal",
                         chains = chains, cores = cores, warmup = warmup, iter = iter)
AGBmods$Stems <- brms::brm(bf(TotAGBha ~ Stems + (Stems | AreaCode),
                            sigma ~ sizeclass),
                         data = plotF_AGB, family = "lognormal",
                         chains = chains, cores = cores, warmup = warmup, iter = iter)

plot(brms::conditional_effects(test5, "BA"), points = T)

# Estimate loo
AGBloos <- list()
AGBloos$intercept_nos <- loo::loo(AGBmods$intercept_nos, reloo=T)
AGBloos$intercept <- loo::loo(AGBmods$intercept, reloo=T)
AGBloos$area <- loo::loo(AGBmods$area, reloo=T)

AGBloos$WSG <- loo::loo(AGBmods$WSG, reloo=T)
AGBloos$BA <- loo::loo(AGBmods$BA, reloo=T)
AGBloos$BAmax <- loo::loo(AGBmods$BAmax, reloo=T)
AGBloos$Stems <- loo::loo(AGBmods$Stems, reloo=T)

# Bayesian R^2
AGBrsq <- list()
AGBrsq$intercept_nos <- brms::bayes_R2(AGBmods$intercept_nos)
AGBrsq$intercept <- brms::bayes_R2(AGBmods$intercept)
AGBrsq$area <- brms::bayes_R2(AGBmods$area)

AGBrsq$WSG <- brms::bayes_R2(AGBmods$WSG)
AGBrsq$BA <- brms::bayes_R2(AGBmods$BA)
AGBrsq$BAmax <- brms::bayes_R2(AGBmods$BAmax)
AGBrsq$Stems <- brms::bayes_R2(AGBmods$Stems)

# Bayesian R^2
AGBloo_rsq <- list()
AGBloo_rsq$intercept_nos <- brms::loo_R2(AGBmods$intercept_nos, reloo=T)
AGBloo_rsq$intercept <- brms::loo_R2(AGBmods$intercept, reloo=T)
AGBloo_rsq$area <- brms::loo_R2(AGBmods$area, reloo=T)

AGBloo_rsq$WSG <- brms::loo_R2(AGBmods$WSG, reloo=T)
AGBloo_rsq$BA <- brms::loo_R2(AGBmods$BA, reloo=T)
AGBloo_rsq$BAmax <- brms::loo_R2(AGBmods$BAmax, reloo=T)
AGBloo_rsq$Stems <- brms::loo_R2(AGBmods$Stems, reloo=T)

# Model comparisons
AGBloo_compare <- loo::loo_compare(AGBloos)

AGBrsq_compare <- sapply(AGBrsq, function(x) x[,1:2])
AGBrsq_compare <- t(round(AGBrsq_compare[order(AGBrsq_compare[,1],decreasing=T),],3))

AGBloo_rsq_compare <- sapply(AGBloo_rsq, function(x) x[1])
AGBloo_rsq_compare <- round(AGBloo_rsq_compare[order(AGBloo_rsq_compare,decreasing=T)],3)
names(AGBloo_rsq_compare) <- stringr::str_remove(names(AGBloo_rsq_compare),".R2")

## Moran's I
plot.dists <- as.matrix(dist(cbind(plotF_AGB$long, plotF_AGB$lat)))
plot.dists.inv <- 1/plot.dists
diag(plot.dists.inv) <- 0
res_AGBintercept <- brms::predictive_error(AGBmods$intercept) ; MI_AGBintercept <- as.data.frame(ape::Moran.I(colMeans(res_AGBintercept),plot.dists.inv))
res_AGBarea <- brms::predictive_error(AGBmods$area) ; MI_AGBarea <- as.data.frame(ape::Moran.I(colMeans(res_AGBarea),plot.dists.inv))

###############
##  OUTPUTS  ##
###############
## Model comparison table
comparetab <- merge(AGBloo_rsq_compare, AGBrsq_compare, by = "row.names") ; row.names(comparetab) <- comparetab$Row.names ; comparetab$Row.names <- NULL
comparetab <- merge(comparetab, AGBloo_compare, by = "row.names")
comparetab <- comparetab[c("Row.names","x","Estimate","Est.Error","elpd_loo","se_elpd_loo","elpd_diff","se_diff")]
modorder <- c("intercept_nos","intercept","area","WSG","Stems","BA","BAmax")
comparetab <- comparetab[order(match(comparetab$Row.names,modorder)),]
names(comparetab) <- c("model","RsqAdj","Rsq","RsqErr","elpd","elpdSD","elpdDiff","elpdDiffSD")

comparetab$estimate <- NA
for(i in 4:nrow(comparetab)){
  comparetab$estimate[i] <- round(summary(AGBmods[comparetab$model[i]][[1]])$fixed[3,1],4)
}

## Moran's I table
MItab <- data.frame(cbind(MI_AGBintercept, MI_AGBarea))

###############
#### SAVE  ####
###############
## Save R environment
save.image(file = paste("Output\\AGBmodels", ".Rdata", sep = ""))
rm(list=ls())

## Save tables
write.csv(comparetab, file = paste("Output\\AGBmod_compare", ".csv", sep = ""), row.names = F)
write.table(MItab, file = "Output\\AGB_MoransI.csv", sep = ",", col.names = NA, row.names = T, append = T)

# Conditional effects
test <- c("WSG","BA","BAmax","Stems")
for(i in 1:length(test)){
  modelcheck <- AGBmods[[i+3]]$fit
  parcheck <- names(modelcheck)[startsWith(names(modelcheck),"b")]
  parcheck2 <- parcheck[-which(endsWith(parcheck,"t"))]
  parcheck2 <- parcheck2[-which(parcheck2 == "b_sigma_sizeclass")]
  covars <- sub("b_","",parcheck2)
  
  png(file = paste("Output\\AGBintervals_", covars, ".png",sep = ""), width = 1024, height = 100)
  bayesplot::mcmc_intervals(modelcheck, pars = parcheck2, prob = 0.75) + ggplot2::ggtitle(test[i])
  dev.off()
  
  png(file = paste("Output\\AGBeffect_", covars, ".png",sep = ""), width = 1024, height = 512)
  plot(brms::conditional_effects(AGBmods[[i+3]],covars),points=T)
  dev.off()
  
  randcheck <- names(modelcheck)[startsWith(names(modelcheck),"r")]
  randcheck2 <- randcheck[-which(endsWith(randcheck,"t]"))]
  
  png(file = paste("Output\\AGB_randomeffects_", covars, ".png",sep = ""), width = 1024, height = 512)
  bayesplot::mcmc_intervals(modelcheck, pars = randcheck2, prob = 0.75)
  dev.off()
}





