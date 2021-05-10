#################
################# Prepare carbon data
#################
if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}
setwd(dir.path)
library(dplyr)

###############
## Read data ##
###############
Treedata <- read.csv("Data\\vegetation\\Trees_prep.csv", header = T, stringsAsFactors = F)
Spatial <- read.csv("Data\\vegetation\\Spatialdata.csv", header = T, stringsAsFactors = F)
Grassdata <- read.csv("Data\\vegetation\\Grass_prep.csv", header = T, stringsAsFactors = F)
Paramodata <- read.csv("Data\\vegetation\\Paramo_prep.csv", header = T, stringsAsFactors = F)
Deaddata <- read.csv("Data\\vegetation\\Dead_prep.csv", header = T, stringsAsFactors = F)

# Remove trees from plots not in spatial data (lacking latitude/longitude data)
Treedata <- dplyr::filter(Treedata, SiteCode %in% Spatial$SiteCode)      
Grassdata <- dplyr::filter(Grassdata, SiteCode %in% Spatial$SiteCode)
Paramodata <- dplyr::filter(Paramodata, SiteCode %in% Spatial$SiteCode)
Deaddata <- dplyr::filter(Deaddata, SiteCode %in% Spatial$SiteCode)

# Quickfix
Spatial[which(Spatial$SiteCode == "01-Sep"),]$SiteCode <- "SEP1"
Spatial[which(Spatial$SiteCode == "02-Sep"),]$SiteCode <- "SEP2"
Spatial[which(Spatial$SiteCode == "03-Sep"),]$SiteCode <- "SEP3"
Spatial[which(Spatial$SiteCode == "sep.01"),]$SiteCode <- "SEP1"
Spatial[which(Spatial$SiteCode == "sep.02"),]$SiteCode <- "SEP2"
Spatial[which(Spatial$SiteCode == "sep.03"),]$SiteCode <- "SEP3"

# Plotdata
plotdata <- dplyr::filter(Spatial, HabitatP %in% c("Forest","Paramo"))
plotdata <- dplyr::filter(plotdata, is.na(FAge) | FAge > 20)
plotdata <- plotdata[,c("SiteCode","Cluster","AreaCode","region","lat","long","Size","HabitatP")]
plotdata <- dplyr::filter(plotdata, AreaCode != "CC")
plotdata$SiteNum <- as.numeric(as.factor(plotdata$SiteCode))
plotdata$ClusNum <- as.numeric(as.factor(plotdata$Cluster))

# Filter data
Deaddata <- Deaddata[which(Deaddata$SiteCode %in% plotdata$SiteCode),]
Grassdata <- Grassdata[which(Grassdata$SiteCode %in% plotdata$SiteCode),]
Treedata <- Treedata[which(Treedata$SiteCode %in% plotdata$SiteCode),]
Paramodata <- Paramodata[which(Paramodata$SiteCode %in% plotdata$SiteCode),]

##########################
## Add spatial grouping ##   Currently: max 20.53 km diameter, min 22.6 km separation - WITHIN habitat
##########################
coords_forest <- sp::SpatialPoints(cbind(plotdata$long[which(plotdata$HabitatP == "Forest")], plotdata$lat[which(plotdata$HabitatP == "Forest")]))
distmat_forest <- sp::spDists(cbind(plotdata$long[which(plotdata$HabitatP == "Forest")], plotdata$lat[which(plotdata$HabitatP == "Forest")]),longlat = TRUE)
hc_forest <- hclust(as.dist(distmat_forest), "ward.D2")
areanum_forest <- cutree(hc_forest, h = 65)

coords_paramo <- sp::SpatialPoints(cbind(plotdata$long[which(plotdata$HabitatP == "Paramo")], plotdata$lat[which(plotdata$HabitatP == "Paramo")]))
distmat_paramo <- sp::spDists(cbind(plotdata$long[which(plotdata$HabitatP == "Paramo")], plotdata$lat[which(plotdata$HabitatP == "Paramo")]),longlat = TRUE)
hc_paramo <- hclust(as.dist(distmat_paramo), "ward.D2")
areanum_paramo <- max(areanum_forest) + cutree(hc_paramo, h = 65)

plotdata$AreaNum <- NA
plotdata[which(plotdata$HabitatP == "Forest"),]$AreaNum <- areanum_forest
plotdata[which(plotdata$HabitatP == "Paramo"),]$AreaNum <- areanum_paramo

fpc::cluster.stats(as.dist(distmat_forest), clustering = cutree(hc_forest, h = 65))$max.diameter
fpc::cluster.stats(as.dist(distmat_forest), clustering = cutree(hc_forest, h = 65))$min.separation
fpc::cluster.stats(as.dist(distmat_paramo), clustering = cutree(hc_paramo, h = 65))$max.diameter
fpc::cluster.stats(as.dist(distmat_paramo), clustering = cutree(hc_paramo, h = 65))$min.separation

# Add biogeography
biogeo <- as.data.frame(plotdata %>% group_by(AreaCode) %>% summarise(region = first(region)))
biogeo <- biogeo %>% mutate(biogeo = case_when(
  region == "central" ~ "Central Andes",
  region == "central-eastern" ~ "Central Andes",
  region == "santa_marta" ~ "Santa Marta",
  region == "amazonia" ~ "Amazonia",
  region == "oriental eastern slope" ~ "East Andes",
  region == "oriental western slope" ~ "East Andes",
  region == "oriental na" ~ "East Andes",
  region == "occidental" ~ "West Andes"))
biogeo <- biogeo %>% mutate(biogeo = case_when(
  AreaCode == "EP" ~ "Magdalena valley",
  AreaCode == "PP" ~ "Magdalena valley",
  TRUE ~ biogeo
))
biogeo$bioNum <- as.numeric(as.factor(biogeo$biogeo))
plotdata <- dplyr::left_join(plotdata, biogeo[,c("AreaCode","biogeo","bioNum")])
plotdata$region <- NULL

######################
## Prepare datasets ##
######################
Treedata <- dplyr::filter(Treedata, !Species %in% c("Palm","Fern"))
AddPoints <- unique(plotdata[which(!plotdata$SiteCode %in% unique(Treedata$SiteCode)),]$SiteCode)
Treedata[(nrow(Treedata)+1):(nrow(Treedata)+length(AddPoints)),]$SiteCode <- AddPoints
Treedata <- dplyr::left_join(Treedata, plotdata, by = "SiteCode")

# Remove species with less than two individuals or no cores
coredSp <- Treedata %>% mutate(cored = WSG/WSG)
coredSp <- coredSp %>% group_by(Species) %>% summarise(total = n(), cored = sum(cored,na.rm=T))
coredSp <- dplyr::filter(coredSp, cored > 1)
coredSp <- dplyr::filter(coredSp, total > cored)

## Add "speciesid" column (1 = identified, 0 = not identified)
Treedata[which(!Treedata$Species %in% coredSp$Species),]$Species <- NA
Treedata[which(is.na(Treedata$Species)),]$Species <- "nospec"
Treedata$speciesid <- ifelse(Treedata$Species == "nospec", 0, 1)

## Calculate tree volume
Treedata <- dplyr::mutate(Treedata, vol = case_when(
  AreaCode == "CC" ~ exp(2.183 - 0.665 * log(DBH_used) + 0.892 * (log(DBH_used)^2) - 0.097 * (log(DBH_used)^3)),  # BMAlvII2.3 - Dry forest
  Species == "Polylepis" ~ 0.0694*(DBH_used^2.35996),                                                             # Polylepis equation (Espinoza & Quispe 2005)
  TRUE ~ exp(2.789 - 1.414 * log(DBH_used) + 1.178 * (log(DBH_used)^2) - 0.118 * (log(DBH_used)^3)),              # BMAlvII2.4 - All forest types
))
Treedata <- left_join(Treedata, Treedata %>% group_by(SiteCode) %>% summarise(totVOL = sum(vol)))
Treedata$VOLw <- Treedata$vol / Treedata$totVOL
Treedata[which(is.na(Treedata$VOLw)),]$VOLw <- 1

## Add to plotdata
plotdata <- left_join(plotdata, Treedata %>% group_by(SiteCode) %>% summarise(nQuercus = sum(Species == "Quercus Humboldtii"), nSpec = sum(!Species %in% c("Quercus Humboldtii","nospec")), ntree = sum(TreeN/TreeN, na.rm=T), ncore = sum(WSG/WSG,na.rm=T), Treevol = sum(vol, na.rm=T)))
plotdata$Treevol <- (plotdata$Treevol / plotdata$Size) * 10000
plotdata <- dplyr::left_join(plotdata, dplyr::filter(Treedata, !is.na(WSG)) %>% group_by(SiteNum) %>% mutate(wsgw = WSG * (vol / sum(vol))) %>% summarise(wsg_w = sum(wsgw)))

## Other datasets
AddPoints <- unique(plotdata[which(!plotdata$SiteCode %in% unique(Deaddata$SiteCode)),]$SiteCode)
Deaddata[(nrow(Deaddata)+1):(nrow(Deaddata)+length(AddPoints)),]$SiteCode <- AddPoints
Deaddata <- dplyr::filter(Deaddata, SiteCode %in% plotdata$SiteCode)
Deaddata <- dplyr::left_join(Deaddata, plotdata[,c("SiteCode","SiteNum","ClusNum","AreaNum")])

AddPoints <- unique(plotdata[which(!plotdata$SiteCode %in% unique(Paramodata$SiteCode)),]$SiteCode)
Paramodata[(nrow(Paramodata)+1):(nrow(Paramodata)+length(AddPoints)),]$SiteCode <- AddPoints
Paramodata <- dplyr::filter(Paramodata, SiteCode %in% plotdata$SiteCode)
Paramodata <- dplyr::left_join(Paramodata, plotdata[,c("SiteCode","SiteNum","ClusNum","AreaNum","Size")])

AddPoints <- unique(plotdata[which(!plotdata$SiteCode %in% unique(Grassdata$SiteCode)),]$SiteCode)
Grassdata[(nrow(Grassdata)+1):(nrow(Grassdata)+length(AddPoints)),]$SiteCode <- AddPoints
Grassdata <- dplyr::filter(Grassdata, SiteCode %in% plotdata$SiteCode)
Grassdata <- dplyr::left_join(Grassdata, plotdata[,c("SiteCode","SiteNum","ClusNum","AreaNum","HabitatP")])

### Remove temporary objects
rm(coords_forest, coords_paramo, distmat_forest, distmat_paramo, hc_forest, hc_paramo, areanum_forest, areanum_paramo)
rm(coredSp, Spatial, AddPoints)

####################################
## WSG model for individual trees ## 
####################################
moddat <- list("n_trees" = length(Treedata$WSG),                      # Number of observations
               "n_cores" = length(which(!is.na(Treedata$WSG))),       # Number of cores
               "n_spec" = length(unique(Treedata$Species)),           # Number of species
               "n_area" = length(unique(plotdata$AreaNum)),           # Number of areas
               "n_cluster" = length(unique(plotdata$ClusNum)),        # Number of clusters
               "n_site" = length(unique(plotdata$SiteNum)),           # Number of sites
               "WSG" = Treedata$WSG[which(!is.na(Treedata$WSG))],     # WSG data
               "core_treenum" = which(!is.na(Treedata$WSG)),          # Cored trees
               "species" = as.numeric(as.factor(Treedata$Species)),   # Vector of species
               "speciesid" = Treedata$speciesid,                      # 0/1 - species determined or not
               "site" = Treedata$SiteNum,                                                                        # Site
               "cluster_site" = data.frame(plotdata %>% group_by(SiteNum) %>% summarise(first(ClusNum)))[,2],    # Cluster for each site
               "area_cluster" = data.frame(plotdata %>% group_by(ClusNum) %>% summarise(first(AreaNum)))[,2])    # Area for each cluster

wsgtreemod <- rstan::stan(file = "STAN\\carbdiv_wsgtreemodel.stan", data = moddat, chains = 3, iter = 2000, cores = 3)

bayesplot::mcmc_trace(wsgtreemod, pars = c("alpha","sigma1","sigma2","mu_spec","sigma_spec","sigma_area","sigma_cluster","sigma_site"))
pairs(wsgtreemod, pars = c("alpha","sigma1","sigma2","mu_spec","sigma_spec","sigma_spec_raw[1]","b_spec[1]"), condition = "energy__")
pairs(wsgtreemod, pars = c("b_area[1]","sigma_area","b_cluster[1]","sigma_cluster","sigma_cluster_raw[1]","b_site[1]","sigma_site","sigma_site_raw[1]"), condition = "energy__")

# Extract WSG and tree AGB estimates
wsgmu <- rstan::extract(wsgtreemod, "mu")[[1]]
wsgsigma <- rstan::extract(wsgtreemod, "sigma")[[1]]
wsgdraws <- sapply(1:nrow(wsgmu), function(x) rnorm(ncol(wsgmu), wsgmu[x,], wsgsigma[x,]))
wsgdraws[which(!is.na(Treedata$WSG)),] <- t(matrix(data = rep(Treedata$WSG[which(!is.na(Treedata$WSG))], ncol(wsgdraws)), ncol = length(Treedata$WSG[which(!is.na(Treedata$WSG))]), nrow = ncol(wsgdraws), byrow = TRUE))
wsgdraws_volw <- wsgdraws * Treedata$VOLw 

wsgdraws_plot <- data.frame(wsgdraws_volw) %>% dplyr::mutate(SiteNum = Treedata$SiteNum) %>% reshape2::melt(id.vars = "SiteNum") %>% reshape2::dcast(SiteNum ~ variable, value.var="value", fun.aggregate = sum)
treedraws_plot <- dplyr::left_join(plotdata[which(!is.na(plotdata$Treevol)), c("SiteNum","Treevol")], wsgdraws_plot)
treedraws_plot[,-c(1:2)] <- (treedraws_plot$Treevol * treedraws_plot[,-c(1:2)]) / 1000
treedraws_plot <- treedraws_plot[,-2]

# Posterior parameter table
summary_wsgtreemod <- rstan::summary(wsgtreemod, pars = c("alpha","sigma1","sigma2","sigma_area","sigma_cluster","sigma_site","mu_spec","sigma_spec"))$summary[,c("mean","2.5%","97.5%","n_eff","Rhat")]

# Scatter plots measured/estimated (stem-level, plot-level)
cols <- RColorBrewer::brewer.pal(9, "Greens")[ceiling(plotdata[order(plotdata$SiteNum),]$ncore/7) + 1]
tiff(file = "Output\\CarbDiv\\WSG_scatterplots.tiff", width = 10000, height = 5000, res = 900)
  layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths = c(0.5,0.5))
  par(mar = c(4, 4, 0.5, 0.5))
  ## Scatter plot: observed vs. estimated individual WSG
  plot(colMeans(rstan::extract(wsgtreemod, pars = "mu")[[1]])[which(!is.na(Treedata$WSG))] ~ Treedata$WSG[which(!is.na(Treedata$WSG))],
       xlab = "Measured WSG", ylab = "Estimated WSG")
  abline(0,1)
  legend("topright", bty = "n",
         paste("Pearson's r =", round(cor(Treedata$WSG[which(!is.na(Treedata$WSG))], colMeans(rstan::extract(wsgtreemod, pars = "mu")[[1]])[which(!is.na(Treedata$WSG))]), 2)))
  ## Scatter plot: observed vs. estimated plot WSG
  plot(rowMeans(wsgdraws_plot[,-1]) ~ plotdata[order(plotdata$SiteNum),]$wsg_w, col = cols, pch=16,
       xlab = "Measured WSGv (average of cores)", ylab = "Estimated WSGv",
       xlim = c(min(rowMeans(wsgdraws_plot[,-1]), plotdata[order(plotdata$SiteNum),]$wsg_w, na.rm = T), max(rowMeans(wsgdraws_plot[,-1]), plotdata[order(plotdata$SiteNum),]$wsg_w, na.rm = T)), ylim = c(min(rowMeans(wsgdraws_plot[,-1]), plotdata[order(plotdata$SiteNum),]$wsg_w, na.rm = T), max(rowMeans(wsgdraws_plot[,-1]), plotdata[order(plotdata$SiteNum),]$wsg_w, na.rm = T)))
  abline(0,1)
  legend("bottomright", legend = lapply(split(plotdata[order(plotdata$SiteNum),]$ncore,ceiling(plotdata[order(plotdata$SiteNum),]$ncore/7)+1),min), pch=16, pt.cex=1.5, cex=1, bty='n',
         col = unique(cols)[order(unique(ceiling(plotdata[order(plotdata$SiteNum),]$ncore/7) + 1))], title = "Number of cores", title.adj = -0.005)
  legend("topleft", bty = "n",
         paste("Pearson's r =", round(cor(rowMeans(wsgdraws_plot[,-1]), plotdata[order(plotdata$SiteNum),]$wsg_w, use = "complete.obs"), 2)))
dev.off()

# Save model object
wsgtreemod@stanmodel@dso <- new("cxxdso")
saveRDS(wsgtreemod, file = "Output\\CarbDiv\\wsgtreemod.rds")

#########################
## Grass biomass model ## 
#########################
moddat <- list("n_site" = length(Grassdata$GrassAGB),                                  # Number of observations
               "n_grass" = length(which(!is.na(Grassdata$GrassAGB))),                  # Number of grass samples
               "n_area" = length(unique(plotdata$AreaNum)),                            # Number of areas
               "n_cluster" = length(unique(plotdata$ClusNum)),
               "n_habitat" = length(unique(Grassdata$HabitatP)),                       # Number of habitats
               "grass" = Grassdata$GrassAGB[which(!is.na(Grassdata$GrassAGB))],        # Grass biomass data
               "habitat" = as.numeric(as.factor(Grassdata %>% group_by(SiteNum) %>% summarise(hab = first(HabitatP)) %>% pull(hab))),                  # Vector of habitat                                           # Area of each observation
               "cluster" = Grassdata %>% group_by(SiteNum) %>% summarise(clus = first(ClusNum)) %>% pull(clus),
               "area_cluster" = plotdata %>% group_by(ClusNum) %>% summarise(area = first(AreaNum)) %>% pull(area),
               "grass_sitenum" = Grassdata$SiteNum[which(!is.na(Grassdata$GrassAGB))])                    # Site of each grass sample

grassmod <- rstan::stan(file = "STAN\\carbdiv_grassmodel.stan", data = moddat, chains = 3, iter = 2000, cores = 3)

bayesplot::mcmc_trace(grassmod, pars = c("alpha[1]","alpha[2]","sigma_hab[1]","sigma_hab[2]","sigma_cluster","sigma_area"))
pairs(grassmod, pars = c("alpha[1]","alpha[2]","sigma_hab[1]","sigma_hab[2]","sigma_cluster","sigma_area"), condition = "energy__")

grassmu <- rstan::extract(grassmod, "mu")[[1]]
grassigma <- rstan::extract(grassmod, "sigma")[[1]]
grassdraws <- sapply(1:nrow(grassmu), function(x) rnorm(ncol(grassmu), grassmu[x,], grassigma[x,]))
grassdraws[which(!is.na(Grassdata$GrassAGB)),] <- t(matrix(data = rep(Grassdata$GrassAGB[which(!is.na(Grassdata$GrassAGB))], ncol(grassdraws)), ncol = length(Grassdata$GrassAGB[which(!is.na(Grassdata$GrassAGB))]), nrow = ncol(grassdraws), byrow = TRUE))
grassdraws <- data.frame(SiteNum = Grassdata$SiteNum, grassdraws)

# Posterior parameter table

# Save model object
grassmod@stanmodel@dso <- new("cxxdso")
saveRDS(grassmod, file = "Output\\CarbDiv\\grassmod.rds")

#############################################
## Predict deadwood biomass across dataset ##
#############################################
Deaddata <- left_join(Deaddata, plotdata[,c("SiteCode","Size")])

# Create DecayFactor column - values from Chao et. al. 2008
Deaddata <- dplyr::mutate(Deaddata, DecayFactor = case_when(
                          DecayClass %in% c(1,2) ~ 0.82,
                          DecayClass %in% c(3,4) ~ 0.66,
                          DecayClass == 5 ~ 0.46))

# Create volume * decayfactor column (DecayVol)
Deaddata <- dplyr::mutate(Deaddata, DecayVol = DecayFactor * Volume)
Decayvolha <- as.data.frame(Deaddata %>% group_by(SiteNum) %>% summarise(TotVol = ((sum(DecayVol, na.rm=T) / first(Size)) * 10000) / 1000000)) # t/ha after multiplying with WSG (g/cm3)

# Calculate biomass
deaddraws_plot <- Decayvolha$TotVol * wsgdraws_plot[,-1]
deaddraws_plot <- data.frame(SiteNum = Decayvolha$SiteNum, deaddraws_plot)

############
## Paramo ##
############
paramodraws <- as.data.frame(Paramodata %>% group_by(SiteNum) %>% summarise(ParamoAGB = (sum(AGB, na.rm=T) / first(Size)) * 10000 / 1000))

##################
## Save objects ##
##################
treedraws_plot <- treedraws_plot[order(treedraws_plot$SiteNum),]
deaddraws_plot <- deaddraws_plot[order(deaddraws_plot$SiteNum),]
grassdraws <- grassdraws[order(grassdraws$SiteNum),]
paramodraws <- paramodraws[order(paramodraws$SiteNum),]

AGBdraws_plot <- replace(treedraws_plot, is.na(treedraws_plot), 0)[,-1] + 
                 replace(deaddraws_plot, is.na(deaddraws_plot), 0)[,-1] + 
                 replace(grassdraws, is.na(grassdraws), 0)[,-1] + 
                 as.data.frame(matrix(rep(paramodraws$ParamoAGB, ncol(treedraws_plot[,-1])), nrow = nrow(treedraws_plot[,-1]), ncol = ncol(treedraws_plot[,-1])))
AGBdraws_plot <- cbind(SiteNum = treedraws_plot[,1], AGBdraws_plot)
AGBlist <- list(AGBdraws = AGBdraws_plot, treedraws = treedraws_plot, deaddraws = deaddraws_plot, grassdraws = grassdraws, paramodraws = paramodraws, wsgdraws = wsgdraws_plot, plotdata = plotdata[order(plotdata$SiteNum),])

saveRDS(AGBlist, "Output\\CarbDiv\\AGBlist.RDS")


##############################
## Upscale to cluster level ##
##############################
#AGBlist <- readRDS("Output\\CarbDiv\\AGBlist.RDS")
spatial <- read.csv("Data\\vegetation\\Spatialdata.csv")
AGBlist$plotdata <- dplyr::left_join(AGBlist$plotdata, spatial[,c("SiteCode","ALOSelev")])   ## Add environmental covariates

AGBlist_forest <- lapply(AGBlist, function(x) x[which(x$SiteNum %in% AGBlist$plotdata[AGBlist$plotdata$HabitatP == "Forest",]$SiteNum),])
AGBlist_paramo <- lapply(AGBlist, function(x) x[which(x$SiteNum %in% AGBlist$plotdata[AGBlist$plotdata$HabitatP == "Paramo",]$SiteNum),])

forestdata <- list(n_point = nrow(AGBlist_forest$plotdata),
                   n_cluster = length(unique((AGBlist_forest$plotdata$ClusNum))),
                   lcarbon_point = log(as.numeric(rowMeans(AGBlist_forest$AGBdraws[order(AGBlist_forest$plotdata$ClusNum),-1]))),
                   elevation = scale(as.data.frame(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(mean(ALOSelev)))[,2])[,1],
                   cluster = as.numeric(as.factor(AGBlist_forest$plotdata[order(AGBlist_forest$plotdata$ClusNum),]$ClusNum)),
                   plotsize = AGBlist_forest$plotdata %>% mutate(plotsize = case_when(Size < 200 ~ 1, Size >= 200 ~ 0)) %>% pull(plotsize))

paramodata <- list(n_point = nrow(AGBlist_paramo$plotdata),
                   n_cluster = length(unique((AGBlist_paramo$plotdata$ClusNum))),
                   lcarbon_point = log(as.numeric(rowMeans(AGBlist_paramo$AGBdraws[order(AGBlist_paramo$plotdata$ClusNum),-1]))),
                   elevation = scale(as.data.frame(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(mean(ALOSelev)))[,2])[,1],
                   cluster = as.numeric(as.factor(AGBlist_paramo$plotdata[order(AGBlist_paramo$plotdata$ClusNum),]$ClusNum)))

clusmod_forest1 <- rstan::stan(file = "STAN\\carbdiv_clus_forest1.stan", data = forestdata, chains = 3, iter = 5000)
clusmod_paramo1 <- rstan::stan(file = "STAN\\carbdiv_clus_paramo1.stan", data = paramodata, chains = 3, iter = 5000, control = list(adapt_delta=0.99))
  # 6 divergent transitions after warmup

clusmod_forest2 <- rstan::stan(file = "STAN\\carbdiv_clus_forest2.stan", data = forestdata, chains = 3, iter = 5000)
clusmod_paramo2 <- rstan::stan(file = "STAN\\carbdiv_clus_paramo2.stan", data = paramodata, chains = 3, iter = 5000)

## CHECKS TEMP
pairs(clusmod_forest1, pars=c("lp__","alpha","sigma_cluster","lsigma_point","beta_elev"))
pairs(clusmod_forest1, pars=c("alpha","sigma_cluster","lsigma_point","beta_elev","carbon_cluster[6]","lsigma_point","sigma_c_raw[6]"))

cluscarb <- rstan::extract(clusmod_forest1,"carbon_cluster")[[1]]
cluscarb2 <- rstan::extract(clusmod_forest2,"carbon_cluster")[[1]]

plot(colMeans(cluscarb)[forestdata$cluster],exp(colMeans(cluscarb2)))

#################
## Save output ##
#################
AGBlist_clus <- list()
AGBlist_clus$AGBdraws <- cbind(as.data.frame(rstan::extract(clusmod_forest1,"carbon_cluster")),as.data.frame(rstan::extract(clusmod_paramo1,"carbon_cluster")))
AGBlist_clus$clusdata <- as.data.frame(AGBlist$plotdata %>% group_by(HabitatP,ClusNum) %>% summarise(AreaCode = first(AreaCode), Cluster = first(Cluster),
                                                                                                         ClusNum = first(ClusNum), AreaNum = first(AreaNum),
                                                                                                         lat = mean(lat), long = mean(long), Size = mean(Size), Habitat = first(HabitatP),
                                                                                                         nQuercus = sum(nQuercus), nSpec = sum(nSpec), ntree = sum(ntree), ncore = sum(ncore)))
AGBlist_clus$plotdata <- AGBlist$plotdata

saveRDS(AGBlist_clus, "Output//CarbDiv//AGBlist_clus2.RDS")





