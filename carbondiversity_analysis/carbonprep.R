#################
################# Prepare carbon data
#################
if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}
setwd(dir.path)
library(dplyr)

bayes_R2_res <- function(y, ypred) {
  e <- ypred - y
  var_ypred <- apply(ypred, 1, var)
  var_e <- apply(e, 1, var)
  var_ypred / (var_ypred + var_e)
}

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
plotdata <- dplyr::filter(plotdata, region != "occidental")
plotdata <- dplyr::filter(plotdata, region != "santa_marta")
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
  #region == "santa_marta" ~ "Santa Marta",
  region == "amazonia" ~ "Amazonia",
  AreaCode %in% c("AL","SM","SE","JU","CH","MO") ~ "East - east slope",
  AreaCode %in% c("CQ","EP","TN","IG","GA","AT","TS","TA","VI","VC","BE","TU","SO","CC","RA","PP","HT") ~ "East - west slope"))
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
rm(coredSp, Spatial, AddPoints, biogeo)

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
               "site" = Treedata$SiteNum,                                                                             # Site for each observation
               "cluster_site" = plotdata %>% group_by(SiteNum) %>% summarise(clus = first(ClusNum)) %>% pull(clus),   # Cluster for each site
               "area_cluster" = plotdata %>% group_by(ClusNum) %>% summarise(area = first(AreaNum)) %>% pull(area))   # Area for each cluster

wsgtreemod <- rstan::stan(file = "STAN\\carbdiv_wsgtreemodel.stan", data = moddat, chains = 4, iter = 5000, cores = 4, thin = 1)

# Extract WSG and tree AGB estimates
wsgmu <- rstan::extract(wsgtreemod, "mu")[[1]]
wsgsigma <- rstan::extract(wsgtreemod, "sigma")[[1]]
wsgmoddraws <- sapply(1:nrow(wsgmu), function(x) rnorm(ncol(wsgmu), wsgmu[x,], wsgsigma[x,]))
#wsgmoddraws <- sapply(seq(1, nrow(wsgmu), 10), function(x) rnorm(ncol(wsgmu), wsgmu[x,], wsgsigma[x,]))
wsgdraws <- wsgmoddraws
wsgdraws[which(!is.na(Treedata$WSG)),] <- t(matrix(data = rep(Treedata$WSG[which(!is.na(Treedata$WSG))], ncol(wsgmoddraws)), ncol = length(Treedata$WSG[which(!is.na(Treedata$WSG))]), nrow = ncol(wsgmoddraws), byrow = TRUE))
wsgdraws_volw <- wsgdraws * Treedata$VOLw 

wsgdraws_plot <- data.frame(wsgdraws_volw) %>% dplyr::mutate(SiteNum = Treedata$SiteNum) %>% reshape2::melt(id.vars = "SiteNum") %>% reshape2::dcast(SiteNum ~ variable, value.var="value", fun.aggregate = sum)
treedraws_plot <- dplyr::left_join(plotdata[which(!is.na(plotdata$Treevol)), c("SiteNum","Treevol")], wsgdraws_plot)
treedraws_plot[,-c(1:2)] <- (treedraws_plot$Treevol * treedraws_plot[,-c(1:2)]) / 1000
treedraws_plot <- treedraws_plot[,-2]

# Posterior parameter table
R2 <- bayes_R2_res(y = t(wsgmoddraws), ypred = wsgmu)
rmse <- sqrt(rowMeans((wsgmoddraws[which(!is.na(Treedata$WSG)),] - moddat$WSG)^2))
postest_wsgtreemod <- rstan::summary(wsgtreemod, pars = c("alpha","sigma1","sigma2","sigma_area","sigma_cluster","sigma_site","mu_spec","sigma_spec"))$summary[,c("mean","2.5%","97.5%")]

postest_wsgtreemod <- rbind(data.frame(postest_wsgtreemod),
                            "RMSE" = data.frame("mean" = mean(rmse), "2.5%" = sort(rmse)[length(rmse)*0.025], "97.5%" = sort(rmse)[length(rmse)*0.975]),
                            "R2" = data.frame("mean" = mean(R2), "2.5%" = sort(R2)[length(R2)*0.025], "97.5%" = sort(R2)[length(R2)*0.975]))
write.csv(postest_wsgtreemod, "Output//CarbDiv//postest_wsgtreemod.csv")

# Posterior predictive checks and nuts 
#bayesplot::mcmc_trace(wsgtreemod, pars = c("alpha","sigma1","sigma2","mu_spec","sigma_spec","sigma_area","sigma_cluster","sigma_site"))
#bayesplot::mcmc_parcoord(wsgtreemod, np = bayesplot::nuts_params(wsgtreemod,div_size = 3), pars = c("alpha","sigma1","sigma2","sigma_area","sigma_cluster","sigma_site","mu_spec","sigma_spec"))
#bayesplot::mcmc_parcoord(wsgtreemod, np = bayesplot::nuts_params(wsgtreemod,div_size = 3), pars = c("b_area[1]","b_cluster[1]","b_site[1]","sigma_area","sigma_cluster_raw[1]","sigma_site_raw[1]"))

pwsg_pairs1 <- bayesplot::mcmc_pairs(wsgtreemod, np = bayesplot::nuts_params(wsgtreemod), pars = c("alpha","sigma1","sigma2","mu_spec","sigma_spec","b_spec[1]","sigma_spec_raw[1]"))
pwsg_pairs2 <- bayesplot::mcmc_pairs(wsgtreemod, np = bayesplot::nuts_params(wsgtreemod), pars = c("b_area[1]","b_cluster[1]","b_site[1]","sigma_area","sigma_cluster","sigma_site","sigma_cluster_raw[1]","sigma_site_raw[1]"))

pwsg_grid <- bayesplot::bayesplot_grid(
                bayesplot::ppc_dens_overlay(moddat$WSG, t(wsgmoddraws)[1:100,which(!is.na(Treedata$WSG))]),
                bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(wsgtreemod)),
                bayesplot::mcmc_rhat(bayesplot::rhat(wsgtreemod)),
                bayesplot::mcmc_neff(bayesplot::neff_ratio(wsgtreemod), size = 2)
              )

ggplot2::ggsave(filename = "Output\\CarbDiv\\wsg_pairs1.png", plot = pwsg_pairs1, dpi = 900, width = 20, height = 20)
ggplot2::ggsave(filename = "Output\\CarbDiv\\wsg_pairs2.png", plot = pwsg_pairs2, dpi = 900, width = 20, height = 20)
ggplot2::ggsave(filename = "Output\\CarbDiv\\wsg_grid.png", plot = pwsg_grid, dpi = 900, width = 20, height = 20)

# Scatter plots measured/estimated (stem-level, plot-level)
cols <- RColorBrewer::brewer.pal(9, "Greens")[ceiling(plotdata[order(plotdata$SiteNum),]$ncore/7) + 1]
tiff(file = "Output\\CarbDiv\\WSG_scatterplots.tiff", width = 10000, height = 5000, res = 900)
  layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths = c(0.5,0.5))
  par(mar = c(4, 4, 0.5, 0.5))
  ## Scatter plot: observed vs. estimated individual WSG
  plot(colMeans(wsgmu)[which(!is.na(Treedata$WSG))] ~ Treedata$WSG[which(!is.na(Treedata$WSG))],
       xlab = "Measured WSG", ylab = "Estimated WSG")
  abline(0,1)
  legend("topright", bty = "n",
         paste("Pearson's r =", round(cor(Treedata$WSG[which(!is.na(Treedata$WSG))], colMeans(wsgmu)[which(!is.na(Treedata$WSG))]), 2)))
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

# WSG by taxonomy plot
treesQuercus <- dplyr::filter(Treedata, Species == "Quercus Humboldtii")
treesID <- dplyr::filter(Treedata, !Species %in% c("Quercus Humboldtii","nospec"))
treesNOID <- dplyr::filter(Treedata, Species == "nospec")

tiff(file = "Output\\CarbDiv\\WSG_taxonomy.tiff", width = 10000, height = 5000, res = 900)
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

# Tree size distribution plot
Treedata$DBH2 <- ifelse(Treedata$DBH > 50, 50, Treedata$DBH)
Treedata$DBHwsg <- ifelse(is.na(Treedata$WSG), NA, Treedata$DBH2)

library(ggplot2)
ggplot(Treedata)+
  theme_classic()+
  scale_x_continuous(breaks = seq(5, 50, 5), expand = c(0,0), 
                     labels = c("5-10","10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50","50+"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0,0))+
  geom_histogram(aes(x = DBH2), fill = "darkgreen", binwidth = 5)+
  geom_histogram(aes(x = DBHwsg), fill = "darkorange", binwidth = 5)+
  labs(x = "DBH range", y = "Count")

ggsave("Output\\CarbDiv\\WSG_treesize.tiff", width = 11, height = 10, units = "cm")

# Moran's I
residuals <- moddat$WSG - colMeans(wsgmu[,moddat$core_treenum])
coordmat <- as.matrix(dist(cbind(Treedata$long[which(!is.na(Treedata$WSG))], Treedata$lat[which(!is.na(Treedata$WSG))])))

dists <- coordmat + 1
dists.inv <- 1/dists
diag(dists.inv) <- 0
dists.inv[is.infinite(dists.inv)] <- 0

as.data.frame(ape::Moran.I(residuals, dists.inv, na.rm = T))

# Save model object
wsgtreemod@stanmodel@dso <- new("cxxdso")
saveRDS(wsgtreemod, file = "Output\\CarbDiv\\wsgtreemod.rds")

#########################
## Grass biomass model ## 
#########################
moddat <- list("n_site" = length(Grassdata$GrassAGB),                                  # Number of observations
               "n_grass" = length(which(!is.na(Grassdata$GrassAGB))),                  # Number of grass samples
               "n_habitat" = length(unique(Grassdata$HabitatP)),                       # Number of habitats
               "n_cluster" = length(unique(plotdata$ClusNum)),                         # Number of clusters
               "grass" = Grassdata$GrassAGB[which(!is.na(Grassdata$GrassAGB))],        # Grass biomass data
               "habitat" = as.numeric(as.factor(Grassdata %>% group_by(SiteNum) %>% summarise(hab = first(HabitatP)) %>% pull(hab))),     # Habitat at each site
               "hab_cluster" = as.numeric(as.factor(Grassdata %>% group_by(ClusNum) %>% summarise(hab = first(HabitatP)) %>% pull(hab))), # Habitat at each cluster
               "cluster" = Grassdata %>% group_by(SiteNum) %>% summarise(clus = first(ClusNum)) %>% pull(clus),                           # Cluster at each site
               "grass_sitenum" = Grassdata$SiteNum[which(!is.na(Grassdata$GrassAGB))])                                                    # Site of each grass sample

grassmod <- rstan::stan(file = "STAN\\carbdiv_grassmodel.stan", data = moddat, chains = 4, iter = 5000, cores = 4, thin = 1)

pairs(grassmod, pars = c("alpha_hab","alpha_clus[1]","sigma_hab","sigma_cluster_hab_raw[1]"))

grassmu <- rstan::extract(grassmod, "mu")[[1]]
grassigma <- rstan::extract(grassmod, "sigma")[[1]]
grassmoddraws <- sapply(1:nrow(grassmu), function(x) rnorm(ncol(grassmu), grassmu[x,], grassigma[x,]))
grassdraws <- grassmoddraws
grassdraws[which(!is.na(Grassdata$GrassAGB)),] <- t(matrix(data = rep(Grassdata$GrassAGB[which(!is.na(Grassdata$GrassAGB))], ncol(grassmoddraws)), ncol = length(Grassdata$GrassAGB[which(!is.na(Grassdata$GrassAGB))]), nrow = ncol(grassmoddraws), byrow = TRUE))
grassdraws <- data.frame(SiteNum = Grassdata$SiteNum, grassdraws)

# Posterior parameter table
R2 <- bayes_R2_res(y = t(grassmoddraws), ypred = grassmu)
rmse <- sqrt(rowMeans((grassmoddraws[which(!is.na(Grassdata$GrassAGB)),] - moddat$grass)^2))
postest_grassmod <- rstan::summary(grassmod, pars = c("alpha_hab","sigma_hab","sigma_cluster_hab"))$summary[,c("mean","2.5%","97.5%")]

postest_grassmod <- rbind(data.frame(postest_grassmod),
                          "RMSE" = data.frame("mean" = mean(rmse), "2.5%" = sort(rmse)[length(rmse)*0.025], "97.5%" = sort(rmse)[length(rmse)*0.975]),
                          "R2" = data.frame("mean" = mean(R2), "2.5%" = sort(R2)[length(R2)*0.025], "97.5%" = sort(R2)[length(R2)*0.975]))
write.csv(postest_grassmod, "Output//CarbDiv//postest_grassmod.csv")

# Posterior predictive checks and nuts 
#bayesplot::mcmc_trace(grassmod, pars = c("alpha_hab[1]","alpha_hab[2]","sigma_hab[1]","sigma_hab[2]","sigma_cluster_hab[1]","sigma_cluster_hab[2]"))
#bayesplot::mcmc_parcoord(grassmod, np = bayesplot::nuts_params(grassmod,div_size = 3), pars = c("alpha_hab[1]","alpha_hab[2]","sigma_hab[1]","sigma_hab[2]","sigma_cluster_hab[1]","sigma_cluster_hab[2]"))

pgrass_pairs <- bayesplot::mcmc_pairs(grassmod, np = bayesplot::nuts_params(grassmod), pars = c("alpha_hab[1]","alpha_hab[2]","sigma_hab[1]","sigma_hab[2]","sigma_cluster_hab[1]","sigma_cluster_hab[2]","sigma_cluster_hab_raw[1]"))

pgrass_grid <- bayesplot::bayesplot_grid(
  bayesplot::ppc_dens_overlay(moddat$grass, t(grassmoddraws)[1:100,which(!is.na(Grassdata$GrassAGB))]),
  bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(grassmod)),
  bayesplot::mcmc_rhat(bayesplot::rhat(grassmod)),
  bayesplot::mcmc_neff(bayesplot::neff_ratio(grassmod), size = 2)
)

ggplot2::ggsave(filename = "Output\\CarbDiv\\grass_pairs.png", plot = pgrass_pairs, dpi = 900, width = 20, height = 20)
ggplot2::ggsave(filename = "Output\\CarbDiv\\grass_grid.png", plot = pgrass_grid, dpi = 900, width = 20, height = 20)

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

AGBlist$plotdata <- dplyr::left_join(AGBlist$plotdata,data.frame(SiteNum = AGBlist$AGBdraws[,1], meanAGB = rowMeans(AGBlist$AGBdraws[,-1])))

saveRDS(AGBlist, "Output\\CarbDiv\\AGBlist.RDS")

##############################
## Upscale to cluster level ##
##############################
AGBlist <- readRDS("Output\\CarbDiv\\AGBlist.RDS")
spatial <- read.csv("Data\\vegetation\\Spatialdata.csv")
AGBlist$plotdata <- dplyr::left_join(AGBlist$plotdata, spatial[,c("SiteCode","ALOSelev","TotPrec","TempVar","PrecVar")])

AGBlist_forest <- lapply(AGBlist, function(x) x[which(x$SiteNum %in% AGBlist$plotdata[AGBlist$plotdata$HabitatP == "Forest",]$SiteNum),])
AGBlist_paramo <- lapply(AGBlist, function(x) x[which(x$SiteNum %in% AGBlist$plotdata[AGBlist$plotdata$HabitatP == "Paramo",]$SiteNum),])

###############################
## Run forest cluster models ##
###############################
# Prepare data
forestpredictors <- list(# Single model
                         "elevation" =     cbind(scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))[,1]),
                         "precipitation" = cbind(scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(prec = mean(TotPrec)) %>% pull(prec))[,1]),
                         "tempvar" =       cbind(scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(tempvar = mean(TempVar)) %>% pull(tempvar))[,1]),
                         "precvar" =       cbind(scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(precvar = mean(PrecVar)) %>% pull(precvar))[,1]),
                         # Elevation + covariate models
                         "elev+prec" =     cbind(scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))[,1],
                                                 scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(prec = mean(TotPrec)) %>% pull(prec))[,1]),
                         "elev+tempvar" =  cbind(scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))[,1],
                                                 scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(tempvar = mean(TempVar)) %>% pull(tempvar))[,1]),
                         "elev+precvar" =  cbind(scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))[,1],
                                                 scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(precvar = mean(PrecVar)) %>% pull(precvar))[,1]))

AGBdraws_short <- AGBlist_forest$AGBdraws[order(AGBlist_forest$AGBdraws$SiteNum), seq(2,dim(AGBlist_forest$AGBdraws[,-1])[2], 20)]
AGBdraws_short <- AGBdraws_short[,1:5]
warnings_forest <- list()
postest_forest <- list()
totaldraws <- 10000

# Run intercept models
modcode <- rstan::stan_model("STAN\\carbdiv_clus_forest_intercept.stan")
warnings <- list(low_bfmi = NULL, rhat = NULL, ess = NULL, divergences = NULL, treedepth = NULL)
  for(i in 1:dim(AGBdraws_short)[2]){
    print(paste("Model ",i," of ", dim(AGBdraws_short)[2], sep = ""))
    # Run model
    forestdata <- list(n_point = dim(AGBdraws_short)[1],                                                                          # Number of points
                       n_cluster = length(unique(AGBlist_forest$plotdata$ClusNum)),                                               # Number of clusters
                       lcarbon_point = log(AGBdraws_short[,i]),                                                                   # Carbon at each point
                       cluster = as.numeric(as.factor(AGBlist_forest$plotdata[order(AGBlist_forest$plotdata$SiteNum),]$ClusNum))) # Cluster at each point
    clusmod_forest <- rstan::sampling(object = modcode, data = forestdata, init = function(){list(alpha = rnorm(1, 200, 20))}, 
                                      chains = 4, cores = 4, iter = 5000, thin = 1, control = list(adapt_delta = 0.99))
    
    # Extract and combine posterior estimates
    postest <- rstan::extract(clusmod_forest)
    postest <- lapply(postest, function(x) if(length(dim(x)) > 1){x[seq(1, dim(postest[[1]]), (dim(postest[[1]]) * dim(AGBdraws_short)[2]) / totaldraws),]}
                                           else{x[seq(1, dim(postest[[1]]), (dim(postest[[1]]) * dim(AGBdraws_short)[2]) / totaldraws)]})
    ifelse(i == 1, postcomb <- postest,
                   postcomb <- lapply(seq_along(postest), function(x) if(length(dim(postest[[x]])) > 1){rbind(postcomb[[x]], postest[[x]])}
                                                                      else{append(postcomb[[x]], postest[[x]])}))

    # Save warnings
    warnings$low_bfmi <- append(warnings$low_bfmi, length(rstan::get_low_bfmi_chains(clusmod_forest)))
    warnings$rhat <- append(warnings$rhat,  length(which(rstan::summary(clusmod_forest)[[1]][,10] > 1.01)))
    warnings$ess <- append(warnings$ess, length(which(rstan::summary(clusmod_forest)[[1]][,9] < 400)))
    warnings$divergences <- append(warnings$divergence, rstan::get_num_divergent(clusmod_forest))
    warnings$treedepth <- append(warnings$treedepth, rstan::get_num_max_treedepth(clusmod_forest))
  }
names(postcomb) <- names(postest)
warnings_forest[[1]] <- warnings
postest_forest[[1]] <- postcomb
names(postest_forest)[[1]] <- "intercept"
names(warnings_forest)[[1]] <- "intercept"

pairs(clusmod_forest, pars = c("alpha","cluster_mean","carbon_cluster[1]","cluster_lmean[1]","lsigma_point","sigma_cluster"))

# Run covariate models
modcode <- rstan::stan_model("STAN\\carbdiv_clus_forest.stan")
modnum <- 0
for(j in 1:(length(forestpredictors))){
  warnings <- list(low_bfmi = NULL, rhat = NULL, ess = NULL, divergences = NULL, treedepth = NULL)
  for(i in 1:dim(AGBdraws_short)[2]){
    modnum <- modnum + 1
    print(paste("Model ", modnum," of ", length(forestpredictors) * dim(AGBdraws_short)[2], sep = ""))
    # Run model
    forestdata <- list(n_point = dim(AGBdraws_short)[1],                                                                          # Number of points
                       n_cluster = length(unique(AGBlist_forest$plotdata$ClusNum)),                                               # Number of clusters
                       n_pred = ncol(forestpredictors[[j]]),                                                                      # Number of predictors
                       lcarbon_point = log(AGBdraws_short[,i]),                                                                   # Carbon at each point
                       cluster = as.numeric(as.factor(AGBlist_forest$plotdata[order(AGBlist_forest$plotdata$SiteNum),]$ClusNum)), # Cluster at each point
                       predictor = forestpredictors[[j]])                                                                         # Predictor at each
    clusmod_forest <- rstan::sampling(object = modcode, data = forestdata, init = function(){list(alpha = rnorm(1, 200, 20), meanfraction = rnorm(1,0.5,0.001))}, 
                                      chains = 4, cores = 4, iter = 5000, thin = 1, control = list(adapt_delta = 0.99))
    
    # Extract and combine posterior estimates
    postest <- rstan::extract(clusmod_forest)
    postest <- lapply(postest, function(x) if(length(dim(x)) > 1){x[seq(1, dim(postest[[1]]), (dim(postest[[1]]) * dim(AGBdraws_short)[2]) / totaldraws),]}
                                           else{x[seq(1, dim(postest[[1]]), (dim(postest[[1]]) * dim(AGBdraws_short)[2]) / totaldraws)]})
    ifelse(i == 1, postcomb <- postest,
                   postcomb <- lapply(seq_along(postest), function(x) if(length(dim(postest[[x]])) > 1){rbind(postcomb[[x]], postest[[x]])}
                                                                      else{append(postcomb[[x]], postest[[x]])}))
    
    # Save warnings
    warnings$low_bfmi <- append(warnings$low_bfmi, length(rstan::get_low_bfmi_chains(clusmod_forest)))
    warnings$rhat <- append(warnings$rhat,  length(which(rstan::summary(clusmod_forest)[[1]][,10] > 1.01)))
    warnings$ess <- append(warnings$ess, length(which(rstan::summary(clusmod_forest)[[1]][,9] < 400)))
    warnings$divergences <- append(warnings$divergence, rstan::get_num_divergent(clusmod_forest))
    warnings$treedepth <- append(warnings$treedepth, rstan::get_num_max_treedepth(clusmod_forest))
  }
  names(postcomb) <- names(postest)
  warnings_forest[[j+1]] <- warnings
  postest_forest[[j+1]] <- postcomb
}
names(postest_forest)[-1] <- names(forestpredictors)
names(warnings_forest)[-1] <- names(forestpredictors)

pairs(clusmod_forest, pars = c("alpha","cluster_mean[1]","carbon_cluster[1]","cluster_lmean[1]","lsigma_point","meanfraction","sigma_cluster[1]","beta"))

# Save output
saveRDS(postest_forest, "Output\\CarbDiv\\postest_clusmod_forest.R")
saveRDS(warnings_forest, "Output\\CarbDivwarnings_clusmod_forest.R")

###############################
## Run paramo cluster models ##
###############################
# Prepare data
paramopredictors <- list(# Single model
  "elevation" =     cbind(scale(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))[,1]),
  "precipitation" = cbind(scale(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(prec = mean(TotPrec)) %>% pull(prec))[,1]),
  "tempvar" =       cbind(scale(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(tempvar = mean(TempVar)) %>% pull(tempvar))[,1]),
  "precvar" =       cbind(scale(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(precvar = mean(PrecVar)) %>% pull(precvar))[,1]),
  # Elevation + covariate models
  "elev+prec" =     cbind(scale(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))[,1],
                          scale(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(prec = mean(TotPrec)) %>% pull(prec))[,1]),
  "elev+tempvar" =  cbind(scale(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))[,1],
                          scale(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(tempvar = mean(TempVar)) %>% pull(tempvar))[,1]),
  "elev+precvar" =  cbind(scale(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))[,1],
                          scale(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(precvar = mean(PrecVar)) %>% pull(precvar))[,1]))

AGBdraws_short <- AGBlist_paramo$AGBdraws[order(AGBlist_paramo$AGBdraws$SiteNum), seq(2,dim(AGBlist_paramo$AGBdraws[,-1])[2], 20)]
AGBdraws_short <- AGBdraws_short[,1:2]
AGBdraws_short[AGBdraws_short < 1] <- 0.01  ## TEMPORARY FIX
warnings_paramo <- list()
postest_paramo <- list()
totaldraws <- 10000

# Run intercept models
modcode <- rstan::stan_model("STAN\\carbdiv_clus_paramo_intercept.stan")
warnings <- list(low_bfmi = NULL, rhat = NULL, ess = NULL, divergences = NULL, treedepth = NULL)
for(i in 1:dim(AGBdraws_short)[2]){
  print(paste("Model ",i," of ", dim(AGBdraws_short)[2], sep = ""))
  # Run model
  paramodata <- list(n_point = dim(AGBdraws_short)[1],                                                                          # Number of points
                     n_cluster = length(unique(AGBlist_paramo$plotdata$ClusNum)),                                               # Number of clusters
                     lcarbon_point = log(AGBdraws_short[,i]),                                                                   # Carbon at each point
                     cluster = as.numeric(as.factor(AGBlist_paramo$plotdata[order(AGBlist_paramo$plotdata$SiteNum),]$ClusNum))) # Cluster at each point
  clusmod_paramo <- rstan::sampling(object = modcode, data = paramodata, init = function(){list(alpha = rnorm(1, 30, 5))}, 
                                    chains = 4, cores = 4, iter = 5000, thin = 1, control = list(adapt_delta = 0.99, max_treedepth = 15))
  
  # Extract and combine posterior estimates
  postest <- rstan::extract(clusmod_paramo)
  postest <- lapply(postest, function(x) if(length(dim(x)) > 1){x[seq(1, dim(postest[[1]]), (dim(postest[[1]]) * dim(AGBdraws_short)[2]) / totaldraws),]}
                    else{x[seq(1, dim(postest[[1]]), (dim(postest[[1]]) * dim(AGBdraws_short)[2]) / totaldraws)]})
  ifelse(i == 1, postcomb <- postest,
         postcomb <- lapply(seq_along(postest), function(x) if(length(dim(postest[[x]])) > 1){rbind(postcomb[[x]], postest[[x]])}
                            else{append(postcomb[[x]], postest[[x]])}))
  
  # Save warnings
  warnings$low_bfmi <- append(warnings$low_bfmi, length(rstan::get_low_bfmi_chains(clusmod_paramo)))
  warnings$rhat <- append(warnings$rhat,  length(which(rstan::summary(clusmod_paramo)[[1]][,10] > 1.01)))
  warnings$ess <- append(warnings$ess, length(which(rstan::summary(clusmod_paramo)[[1]][,9] < 400)))
  warnings$divergences <- append(warnings$divergence, rstan::get_num_divergent(clusmod_paramo))
  warnings$treedepth <- append(warnings$treedepth, rstan::get_num_max_treedepth(clusmod_paramo))
}
names(postcomb) <- names(postest)
warnings_paramo[[1]] <- warnings
postest_paramo[[1]] <- postcomb
names(postest_paramo)[[1]] <- "intercept"
names(warnings_paramo[[1]]) <- "intercept"

pairs(clusmod_paramo, pars = c("alpha","cluster_mean","carbon_cluster[1]","cluster_lmean[1]","lsigma_point","sigma_cluster","sigma_cluster_raw[1]"))

# Run covariate models
modcode <- rstan::stan_model("STAN\\carbdiv_clus_paramo.stan")
modnum <- 0
for(j in 1:(length(paramopredictors))){
  warnings <- list(low_bfmi = NULL, rhat = NULL, ess = NULL, divergences = NULL, treedepth = NULL)
  for(i in 1:dim(AGBdraws_short)[2]){
    modnum <- modnum + 1
    print(paste("Model ", modnum," of ", length(forestpredictors) * dim(AGBdraws_short)[2], sep = ""))
    # Run model
    paramodata <- list(n_point = dim(AGBdraws_short)[1],                                                                          # Number of points
                       n_cluster = length(unique(AGBlist_paramo$plotdata$ClusNum)),                                               # Number of clusters
                       n_pred = ncol(paramopredictors[[j]]),                                                                      # Number of predictors
                       lcarbon_point = log(AGBdraws_short[,i]),                                                                   # Carbon at each point
                       cluster = as.numeric(as.factor(AGBlist_paramo$plotdata[order(AGBlist_paramo$plotdata$SiteNum),]$ClusNum)), # Cluster at each point
                       predictor = paramopredictors[[j]])                                                                         # Predictor at each
    clusmod_paramo <- rstan::sampling(object = modcode, data = paramodata, init = function(){list(alpha = rnorm(1, 30, 5))}, 
                                      chains = 4, cores = 4, iter = 5000, thin = 1, control = list(adapt_delta = 0.99, max_treedepth = 15))
    
    # Extract and combine posterior estimates
    postest <- rstan::extract(clusmod_paramo)
    postest <- lapply(postest, function(x) if(length(dim(x)) > 1){x[seq(1, dim(postest[[1]]), (dim(postest[[1]]) * dim(AGBdraws_short)[2]) / totaldraws),]}
                      else{x[seq(1, dim(postest[[1]]), (dim(postest[[1]]) * dim(AGBdraws_short)[2]) / totaldraws)]})
    ifelse(i == 1, postcomb <- postest,
           postcomb <- lapply(seq_along(postest), function(x) if(length(dim(postest[[x]])) > 1){rbind(postcomb[[x]], postest[[x]])}
                              else{append(postcomb[[x]], postest[[x]])}))
    
    # Save warnings
    warnings$low_bfmi <- append(warnings$low_bfmi, length(rstan::get_low_bfmi_chains(clusmod_paramo)))
    warnings$rhat <- append(warnings$rhat,  length(which(rstan::summary(clusmod_paramo)[[1]][,10] > 1.01)))
    warnings$ess <- append(warnings$ess, length(which(rstan::summary(clusmod_paramo)[[1]][,9] < 400)))
    warnings$divergences <- append(warnings$divergence, rstan::get_num_divergent(clusmod_paramo))
    warnings$treedepth <- append(warnings$treedepth, rstan::get_num_max_treedepth(clusmod_paramo))
  }
  names(postcomb) <- names(postest)
  warnings_paramo[[j+1]] <- warnings
  postest_paramo[[j+1]] <- postcomb
}
names(postest_paramo)[-1] <- names(paramopredictors)
names(warnings_paramo)[-1] <- names(paramopredictors)

pairs(clusmod_paramo, pars = c("alpha","cluster_mean[1]","carbon_cluster[1]","cluster_lmean[1]","lsigma_point","sigma_cluster","sigma_cluster_raw[1]","beta"))

# Save output
saveRDS(postest_paramo, "Output\\CarbDiv\\postest_clusmod_paramo.R")
saveRDS(warnings_paramo, "Output\\CarbDiv\\warnings_clusmod_paramo.R")

#################################
## Output from selected models ##
#################################
AGBlist <- readRDS("Output\\CarbDiv\\AGBlist.RDS")
spatial <- read.csv("Data\\vegetation\\Spatialdata.csv")
AGBlist$plotdata <- dplyr::left_join(AGBlist$plotdata, spatial[,c("SiteCode","ALOSelev","TotPrec","TempVar","PrecVar")])

AGBlist_forest <- lapply(AGBlist, function(x) x[which(x$SiteNum %in% AGBlist$plotdata[AGBlist$plotdata$HabitatP == "Forest",]$SiteNum),])
AGBlist_paramo <- lapply(AGBlist, function(x) x[which(x$SiteNum %in% AGBlist$plotdata[AGBlist$plotdata$HabitatP == "Paramo",]$SiteNum),])

postest_paramo <- readRDS("Output\\CarbDiv\\postest_clusmod_paramo.R")
postest_forest <- readRDS("Output\\CarbDiv\\postest_clusmod_forest.R")

forestdata <- list(n_point = dim(AGBdraws_short)[1],                                                                          # Number of points
                   n_cluster = length(unique(AGBlist_forest$plotdata$ClusNum)),                                               # Number of clusters
                   n_pred = ncol(forestpredictors[[j]]),                                                                      # Number of predictors
                   lcarbon_point = log(AGBdraws_short[,i]),                                                                   # Carbon at each point
                   cluster = as.numeric(as.factor(AGBlist_forest$plotdata[order(AGBlist_forest$plotdata$SiteNum),]$ClusNum)), # Cluster at each point
                   predictor = forestpredictors[[j]])   
paramodata <- list(n_point = dim(AGBdraws_short)[1],                                                                          # Number of points
                   n_cluster = length(unique(AGBlist_paramo$plotdata$ClusNum)),                                               # Number of clusters
                   n_pred = ncol(paramopredictors[[j]]),                                                                      # Number of predictors
                   lcarbon_point = log(AGBdraws_short[,i]),                                                                   # Carbon at each point
                   cluster = as.numeric(as.factor(AGBlist_paramo$plotdata[order(AGBlist_paramo$plotdata$SiteNum),]$ClusNum)), # Cluster at each point
                   predictor = paramopredictors[[j]])        

forest_cluster_point <-  as.numeric(as.factor(AGBlist_forest$plotdata[order(AGBlist_forest$plotdata$SiteNum),]$ClusNum))
paramo_cluster_point <- as.numeric(as.factor(AGBlist_paramo$plotdata[order(AGBlist_paramo$plotdata$SiteNum),]$ClusNum))

# Parameter estimate table - forest
posttab_forest <- list()
for(i in 1:length(postest_forest)){
  AGBmu <- postest_forest[[i]]$cluster_lmean
  AGBsigma <- postest_forest[[i]]$lsigma_point
  AGBpoint <-  sapply(1:nrow(AGBmu), function(x) rnorm(ncol(AGBmu), AGBmu[x,], AGBsigma[x]))
  R2 <- bayes_R2_res(y = t(AGBpoint), ypred = AGBmu)
  rmse <- sqrt(rowMeans((AGBpoint[forest_cluster_point,] - log(rowMeans(AGBlist_forest$AGBdraws)))^2))
  elpd <- loo::loo(postest_forest[[i]]$log_lik, moment_match = T, k_threshold = 0.5)
  
  posttab_forest[[i]] <- rbind(data.frame(mean = unlist(lapply(postest_forest[[i]][names(postest_forest[[i]][which(names(postest_forest[[i]]) %in% c("alpha","beta","lsigma_point","meanfraction"))])],
                                                          function(x) if(is.null(dim(x))){mean(x)}else{colMeans(x)})),
                                     CI2.5 = unlist(lapply(postest_forest[[i]][names(postest_forest[[i]][which(names(postest_forest[[i]]) %in% c("alpha","beta","lsigma_point","meanfraction"))])],
                                                           function(x) if(is.null(dim(x))){x[order(x)][length(x)*0.025]}else{apply(x, 2, function(y) y[order(y)][length(y) * 0.025])})),
                                     CI97.5 = unlist(lapply(postest_forest[[i]][names(postest_forest[[i]][which(names(postest_forest[[i]]) %in% c("alpha","beta","lsigma_point","meanfraction"))])],
                                                            function(x) if(is.null(dim(x))){x[order(x)][length(x)*0.975]}else{apply(x, 2, function(y) y[order(y)][length(y) * 0.975])}))),
                          "RMSE" = data.frame("mean" = mean(rmse), "CI2.5" = sort(rmse)[length(rmse)*0.025], "CI97.5" = sort(rmse)[length(rmse)*0.975]),
                          "R2" = data.frame("mean" = mean(R2), "CI2.5" = sort(R2)[length(R2)*0.025], "CI97.5" = sort(R2)[length(R2)*0.975]),
                          "elpd" = data.frame("mean" = elpd$estimates[1,1], "CI2.5" = elpd$estimates[1,2], "CI97.5" = NA))
}
names(posttab_forest) <- names(postest_forest)

write.csv(do.call(rbind, posttab_forest), "Output\\CarbDiv\\postest_clusmod_forest_temp.csv")

# Parameter estimate table - forest
posttab_paramo <- list()
for(i in 1:length(postest_paramo)){
  AGBmu <- postest_paramo[[i]]$cluster_lmean
  AGBsigma <- postest_paramo[[i]]$lsigma_point
  AGBpoint <-  sapply(1:nrow(AGBmu), function(x) rnorm(ncol(AGBmu), AGBmu[x,], AGBsigma[x]))
  R2 <- bayes_R2_res(y = t(AGBpoint), ypred = AGBmu)
  rmse <- sqrt(rowMeans((AGBpoint[paramo_cluster_point,] - log(rowMeans(AGBlist_paramo$AGBdraws)))^2))
  elpd <- loo::loo(postest_paramo[[i]]$log_lik, moment_match = T, k_threshold = 0.5)
  
  posttab_paramo[[i]] <- rbind(data.frame(mean = unlist(lapply(postest_paramo[[i]][names(postest_paramo[[i]][which(names(postest_paramo[[i]]) %in% c("alpha","beta","lsigma_point","meanfraction"))])],
                                                               function(x) if(is.null(dim(x))){mean(x)}else{colMeans(x)})),
                                          CI2.5 = unlist(lapply(postest_paramo[[i]][names(postest_paramo[[i]][which(names(postest_paramo[[i]]) %in% c("alpha","beta","lsigma_point","meanfraction"))])],
                                                                function(x) if(is.null(dim(x))){x[order(x)][length(x)*0.025]}else{apply(x, 2, function(y) y[order(y)][length(y) * 0.025])})),
                                          CI97.5 = unlist(lapply(postest_paramo[[i]][names(postest_paramo[[i]][which(names(postest_paramo[[i]]) %in% c("alpha","beta","lsigma_point","meanfraction"))])],
                                                                 function(x) if(is.null(dim(x))){x[order(x)][length(x)*0.975]}else{apply(x, 2, function(y) y[order(y)][length(y) * 0.975])}))),
                               "RMSE" = data.frame("mean" = mean(rmse), "CI2.5" = sort(rmse)[length(rmse)*0.025], "CI97.5" = sort(rmse)[length(rmse)*0.975]),
                               "R2" = data.frame("mean" = mean(R2), "CI2.5" = sort(R2)[length(R2)*0.025], "CI97.5" = sort(R2)[length(R2)*0.975]),
                               "elpd" = data.frame("mean" = elpd$estimates[1,1], "CI2.5" = elpd$estimates[1,2], "CI97.5" = NA))
}
names(posttab_paramo) <- names(postest_paramo)

write.csv(do.call(rbind, posttab_paramo), "Output\\CarbDiv\\postest_clusmod_paramo_temp.csv")

### REMOVE:
      # LOO comparison table
      loos_forest <- lapply(postest_forest, function(x) loo::loo(x$log_lik, moment_match = T, k_threshold = 0.5))
      as.table(loo::loo_compare(loos_forest))
      
      # LOO comparison table
      loos_paramo <- lapply(postest_paramo, function(x) loo::loo(x$log_lik, moment_match = T, k_threshold = 0.5))
      as.table(loo::loo_compare(loos_paramo))
      
      # Extract AGB estimates - forest model
      AGBmu_forest <- posterior_forest$cluster_lmean
      AGBsigma_forest <- posterior_forest$lsigma_point
      AGBpoint_forest <-  sapply(1:nrow(AGBmu_forest), function(x) rnorm(ncol(AGBmu_forest), AGBmu_forest[x,], AGBsigma_forest[x]))
      R2 <- bayes_R2_res(y = t(AGBpoint_forest), ypred = AGBmu_forest)
      rmse <- sqrt(rowMeans((AGBpoint_forest[forestdata$cluster,] - forestdata$lcarbon_point)^2))
      posttab_forest <- rbind(data.frame(mean = unlist(lapply(posterior_forest[names(posterior_forest[which(names(posterior_forest) %in% c("alpha","beta","lsigma_point","meanfraction"))])],
                                                             function(x) if(is.null(dim(x))){mean(x)}else{colMeans(x)})),
                                         CI2.5 = unlist(lapply(posterior_forest[names(posterior_forest[which(names(posterior_forest) %in% c("alpha","beta","lsigma_point","meanfraction"))])],
                                                               function(x) if(is.null(dim(x))){x[order(x)][length(x)*0.025]}else{apply(posterior_forest$sigma_cluster_raw, 2, function(y) y[order(y)][length(y) * 0.025])})),
                                         CI97.5 = unlist(lapply(posterior_forest[names(posterior_forest[which(names(posterior_forest) %in% c("alpha","beta","lsigma_point","meanfraction"))])],
                                                                function(x) if(is.null(dim(x))){x[order(x)][length(x)*0.975]}else{apply(posterior_forest$sigma_cluster_raw, 2, function(y) y[order(y)][length(y) * 0.975])}))),
                             "RMSE" = data.frame("mean" = mean(rmse), "CI2.5" = sort(rmse)[length(rmse)*0.025], "CI97.5" = sort(rmse)[length(rmse)*0.975]),
                             "R2" = data.frame("mean" = mean(R2), "CI2.5" = sort(R2)[length(R2)*0.025], "CI97.5" = sort(R2)[length(R2)*0.975]))
      write.csv(posttab_forest, "Output//CarbDiv//postest_clusmod_forest.csv")
      
      # Extract AGB estimates - paramo model
      AGBmu_paramo <- posterior_paramo$cluster_lmean
      AGBsigma_paramo <- posterior_paramo$lsigma_point
      AGBpoint_paramo <-  sapply(1:nrow(AGBmu_paramo), function(x) rnorm(ncol(AGBmu_paramo), AGBmu_paramo[x,], AGBsigma_paramo[x]))
      R2 <- bayes_R2_res(y = t(AGBpoint_paramo), ypred = AGBmu_paramo)
      rmse <- sqrt(rowMeans((AGBpoint_paramo[paramodata$cluster,] - paramodata$lcarbon_point)^2))
      posttab_paramo <- rbind(data.frame(mean = unlist(lapply(posterior_paramo[names(posterior_paramo[which(names(posterior_paramo) %in% c("alpha","beta","lsigma_point"))])],
                                                              function(x) if(is.null(dim(x))){mean(x)}else{colMeans(x)})),
                                         CI2.5 = unlist(lapply(posterior_paramo[names(posterior_paramo[which(names(posterior_paramo) %in% c("alpha","beta","lsigma_point"))])],
                                                               function(x) if(is.null(dim(x))){x[order(x)][length(x)*0.025]}else{apply(posterior_paramo$sigma_cluster_raw, 2, function(y) y[order(y)][length(y) * 0.025])})),
                                         CI97.5 = unlist(lapply(posterior_paramo[names(posterior_paramo[which(names(posterior_paramo) %in% c("alpha","beta","lsigma_point"))])],
                                                                function(x) if(is.null(dim(x))){x[order(x)][length(x)*0.975]}else{apply(posterior_paramo$sigma_cluster_raw, 2, function(y) y[order(y)][length(y) * 0.975])}))),
                             "RMSE" = data.frame("mean" = mean(rmse), "CI2.5" = sort(rmse)[length(rmse)*0.025], "CI97.5" = sort(rmse)[length(rmse)*0.975]),
                             "R2" = data.frame("mean" = mean(R2), "CI2.5" = sort(R2)[length(R2)*0.025], "CI97.5" = sort(R2)[length(R2)*0.975]))
      write.csv(posttab_paramo, "Output//CarbDiv//postest_clusmod_paramo.csv")
      
# Selected models
posterior_forest <- postest_forest$`elev+prec`
posterior_paramo <- postest_paramo$intercept
      
# Errors plot forest model
point_lmean <- posterior_forest$cluster_lmean[,forestdata$cluster]
tiff(file = "Output\\CarbDiv\\forestresiduals.tiff", width = 10000, height = 5000, res = 900)
  par(mfrow = c(2,3))
  plot(rowMeans(forestdata$lcarbon_point - t(point_lmean)) ~ rowMeans(t(point_lmean)),
       xlab = "Estimated expected point carbon", ylab = "Observed minus expected point carbon") ; abline(h = 0)
  legend("topleft", bty = "n",
         paste("Pearson's r =", round(cor(forestdata$lcarbon_point, colMeans(point_lmean)), 2)))
  plot(aggregate(exp(forestdata$lcarbon_point), list(forestdata$cluster), mean)[,2] - colMeans(posterior_forest$cluster_mean) ~ colMeans(posterior_forest$cluster_mean), 
       xlab = "Estimated expected cluster carbon", ylab = "Observed minus expected cluster carbon") ; abline(h = 0)
  plot(colMeans(posterior_forest$carbon_cluster) - colMeans(posterior_forest$cluster_mean) ~ colMeans(posterior_forest$cluster_mean), 
       xlab = "Estimated expected cluster carbon", ylab = "Estimated minus expected cluster carbon") ; abline(h = 0)
  plot(aggregate(exp(forestdata$lcarbon_point), list(forestdata$cluster), mean)[,2] ~ colMeans(posterior_forest$carbon_cluster), 
       xlab = "Estimated cluster carbon", ylab = "Observed cluster carbon") ; abline(0,1)
  legend("topleft", bty = "n",
         paste("Pearson's r =", round(cor(aggregate(exp(forestdata$lcarbon_point), list(forestdata$cluster), mean)[,2], colMeans(posterior_forest$carbon_cluster)), 2)))
  plot(aggregate(exp(forestdata$lcarbon_point), list(forestdata$cluster), mean)[,2] ~ forestdata$predictor[,1],
       xlab = "Elevation", ylab = "Observed average cluster value")
  points(mean(posterior_forest$alpha) + mean(posterior_forest$beta) * forestdata$predictor[,1] ~ forestdata$predictor[,1], col = "red")
  plot(colMeans(posterior_forest$carbon_cluster) ~ forestdata$predictor[,1],
       xlab = "Elevation", ylab = "Estimated cluster value")
  points(mean(posterior_forest$alpha) + mean(posterior_forest$beta) * forestdata$predictor[,1] ~ forestdata$predictor[,1], col = "red")
dev.off()

# Errors plot paramo model
point_lmean <- posterior_paramo$cluster_lmean[,paramodata$cluster]
cluster_mean <- rep(mean(posterior_paramo$cluster_mean), length(unique(paramodata$cluster)))
tiff(file = "Output\\CarbDiv\\paramoresiduals.tiff", width = 10000, height = 5000, res = 900)
  par(mfrow = c(2,3))
  plot(rowMeans(paramodata$lcarbon_point - t(point_lmean)) ~ colMeans(point_lmean),
       xlab = "Estimated expected point carbon", ylab = "Observed minus expected point carbon") ; abline(h = 0)
  legend("topleft", bty = "n",
         paste("Pearson's r =", round(cor(paramodata$lcarbon_point, colMeans(point_lmean)), 2)))
  plot(aggregate(exp(paramodata$lcarbon_point), list(paramodata$cluster), mean)[,2] - cluster_mean ~ cluster_mean, 
       xlab = "Estimated expected cluster carbon", ylab = "Observed minus expected cluster carbon") ; abline(h = 0)
  plot(colMeans(posterior_paramo$carbon_cluster) - cluster_mean ~ cluster_mean, 
       xlab = "Estimated expected cluster carbon", ylab = "Estimated minus expected cluster carbon") ; abline(h = 0)
  plot(aggregate(exp(paramodata$lcarbon_point), list(paramodata$cluster), mean)[,2] ~ colMeans(posterior_paramo$carbon_cluster), 
       xlab = "Estimated cluster carbon", ylab = "Observed cluster carbon") ; abline(0,1)
  legend("topleft", bty = "n",
         paste("Pearson's r =", round(cor(aggregate(exp(paramodata$lcarbon_point), list(paramodata$cluster), mean)[,2], colMeans(posterior_paramo$carbon_cluster)), 2)))
  plot(aggregate(exp(paramodata$lcarbon_point), list(paramodata$cluster), mean)[,2] ~ paramodata$predictor[,1],
       xlab = "Elevation", ylab = "Observed average cluster value")
  plot(colMeans(posterior_paramo$carbon_cluster) ~ paramodata$predictor[,1],
       xlab = "Elevation", ylab = "Estimated cluster value")
dev.off()

# Moran's I
residuals_forest <- forestdata$lcarbon_point - colMeans(posterior_forest$cluster_lmean)[forestdata$cluster]
coordmat_forest <- 1 / (as.matrix(dist(cbind(AGBlist_forest$plotdata[order(AGBlist_forest$plotdata$SiteNum),]$long, AGBlist_forest$plotdata[order(AGBlist_forest$plotdata$SiteNum),]$lat))) + 1)
diag(coordmat_forest) <- 0
coordmat_forest[is.infinite(coordmat_forest)] <- 0
as.data.frame(ape::Moran.I(residuals_forest, coordmat_forest, na.rm = T))

residuals_paramo <- paramodata$lcarbon_point - colMeans(posterior_paramo$cluster_lmean)[paramodata$cluster]
coordmat_paramo <- 1 / (as.matrix(dist(cbind(AGBlist_paramo$plotdata[order(AGBlist_paramo$plotdata$SiteNum),]$long, AGBlist_paramo$plotdata[order(AGBlist_paramo$plotdata$SiteNum),]$lat))) + 1)
diag(coordmat_paramo) <- 0
coordmat_paramo[is.infinite(coordmat_paramo)] <- 0
as.data.frame(ape::Moran.I(residuals_paramo, coordmat_paramo, na.rm = T))


#################
## Save output ##
#################
AGBlist_clus <- list()
AGBlist_clus$AGBdraws <- as.data.frame(t(cbind(posterior_forest$carbon_cluster, posterior_paramo$carbon_cluster)))
AGBlist_clus$AGBdraws <- data.frame(ClusNum = c(sort(unique(AGBlist_forest$plotdata[order(AGBlist_forest$plotdata$SiteNum),]$ClusNum)), sort(unique(AGBlist_paramo$plotdata[order(AGBlist_paramo$plotdata$SiteNum),]$ClusNum))), AGBlist_clus$AGBdraws)
AGBlist_clus$clusdata <- as.data.frame(AGBlist$plotdata %>% group_by(HabitatP,ClusNum) %>% summarise(AreaCode = first(AreaCode), Cluster = first(Cluster), biogeo = first(biogeo),
                                                                                                     ClusNum = first(ClusNum), AreaNum = first(AreaNum), bioNum = first(bioNum),
                                                                                                     lat = mean(lat), long = mean(long), Size = mean(Size), Habitat = first(HabitatP),
                                                                                                     nQuercus = sum(nQuercus), nSpec = sum(nSpec), ntree = sum(ntree), ncore = sum(ncore),
                                                                                                     meanAGB = mean(meanAGB),
                                                                                                     ALOSelev = mean(ALOSelev), TotPrec = mean(TotPrec), TempVar = mean(TempVar), PrecVar = mean(PrecVar)))
AGBlist_clus$plotdata <- AGBlist$plotdata

saveRDS(AGBlist_clus, "Output//CarbDiv//AGBlist_clus.RDS")











### UNUSED - FIX


# Posterior predictive checks and nuts - forest model
bayesplot::mcmc_trace(clusmod_forest, pars = c("alpha","cluster_mean[1]","a","b","lsigma_small","lsigma_offset","sigma_area","alpha_area[1]","beta[1]","beta[2]","beta[3]","beta[4]","carbon_cluster[1]","sigma_cluster[1]"))
bayesplot::mcmc_parcoord(clusmod_forest, np = bayesplot::nuts_params(clusmod_forest,div_size = 3), pars = c("alpha","cluster_mean[1]","lsigma_small","lsigma_offset","sigma_area","alpha_area[1]","beta[1]","beta[2]","beta[3]","beta[4]","carbon_cluster[1]","sigma_cluster"))

pclusforest_pairs1 <- bayesplot::mcmc_pairs(clusmod_forest, np = bayesplot::nuts_params(clusmod_forest), pars = c("alpha","sigma_area","alpha_area[1]","beta[1]","beta[2]","beta[3]","beta[4]"))
pclusforest_pairs2 <- bayesplot::mcmc_pairs(clusmod_forest, np = bayesplot::nuts_params(clusmod_forest), pars = c("lsigma_small","lsigma_offset","carbon_cluster[1]","cluster_mean[1]","sigma_cluster"))

pclusforest_grid <- bayesplot::bayesplot_grid(
  bayesplot::ppc_dens_overlay(forestdata$lcarbon_point, t(AGBpoint_forest)[1:100,]),
  bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(clusmod_forest)),
  bayesplot::mcmc_rhat(bayesplot::rhat(clusmod_forest)),
  bayesplot::mcmc_neff(bayesplot::neff_ratio(clusmod_forest), size = 2)
)

ggplot2::ggsave(filename = "Output\\CarbDiv\\clusforest_pairs1.png", plot = pclusforest_pairs1, dpi = 900, width = 20, height = 20)
ggplot2::ggsave(filename = "Output\\CarbDiv\\clusforest_pairs2.png", plot = pclusforest_pairs2, dpi = 900, width = 20, height = 20)
ggplot2::ggsave(filename = "Output\\CarbDiv\\clusforest_grid.png", plot = pclusforest_grid, dpi = 900, width = 20, height = 20)

# Posterior predictive checks and nuts- paramo model
bayesplot::mcmc_trace(clusmod_paramo, pars = c("alpha","cluster_mean[1]","lsigma_point","beta[1]","carbon_cluster[1]","sigma_cluster"))
bayesplot::mcmc_parcoord(clusmod_paramo, np = bayesplot::nuts_params(clusmod_paramo,div_size = 3), pars = c("alpha","cluster_mean[1]","lsigma_point","beta[1]","carbon_cluster[1]","sigma_cluster"))

pclusparamo_pairs1 <- bayesplot::mcmc_pairs(clusmod_paramo, np = bayesplot::nuts_params(clusmod_paramo), pars = c("alpha","sigma_area","alpha_area[1]","beta[1]","lsigma_point","carbon_cluster[1]","cluster_mean[1]","sigma_cluster"))
pclusparamo_pairs2 <- bayesplot::mcmc_pairs(clusmod_paramo, np = bayesplot::nuts_params(clusmod_paramo), pars = c("carbon_cluster[1]","cluster_mean[1]","sigma_cluster","sigma_c_raw[1]","alpha_area[1]","sigma_area","sigma_area_raw[1]"))

pclusparamo_grid <- bayesplot::bayesplot_grid(
  bayesplot::ppc_dens_overlay(paramodata$lcarbon_point, t(AGBpoint_paramo)[1:100,]),
  bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(clusmod_paramo)),
  bayesplot::mcmc_rhat(bayesplot::rhat(clusmod_paramo)),
  bayesplot::mcmc_neff(bayesplot::neff_ratio(clusmod_paramo), size = 2)
)

ggplot2::ggsave(filename = "Output\\CarbDiv\\clusparamo_pairs1.png", plot = pclusparamo_pairs1, dpi = 900, width = 20, height = 20)
ggplot2::ggsave(filename = "Output\\CarbDiv\\clusparamo_pairs2.png", plot = pclusparamo_pairs2, dpi = 900, width = 20, height = 20)
ggplot2::ggsave(filename = "Output\\CarbDiv\\clusparamo_grid.png", plot = pclusparamo_grid, dpi = 900, width = 20, height = 20)

