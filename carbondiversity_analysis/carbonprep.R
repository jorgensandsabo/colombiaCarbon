#################
################# Prepare carbon data
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
Grassdata <- read.csv("Data\\vegetation\\Grass_prep.csv", header = T, stringsAsFactors = F)
Paramodata <- read.csv("Data\\vegetation\\Paramo_prep.csv", header = T, stringsAsFactors = F)
Deaddata <- read.csv("Data\\vegetation\\Dead_prep.csv", header = T, stringsAsFactors = F)
PTrees <- read.csv("Data\\vegetation\\Ptrees_prep.csv", header = T, stringsAsFactors = F)

Treedata <- dplyr::filter(Treedata, SiteCode %in% Spatial$SiteCode)      # Remove trees from plots not in spatial data (lacking latitude/longitude data)
Grassdata <- dplyr::filter(Grassdata, SiteCode %in% Spatial$SiteCode)
Paramodata <- dplyr::filter(Paramodata, SiteCode %in% Spatial$SiteCode)
Deaddata <- dplyr::filter(Deaddata, SiteCode %in% Spatial$SiteCode)
PTrees <- dplyr::filter(PTrees, SiteCode %in% Spatial$SiteCode)

Spatial[which(Spatial$SiteCode == "01-Sep"),]$SiteCode <- "SEP1"
Spatial[which(Spatial$SiteCode == "02-Sep"),]$SiteCode <- "SEP2"
Spatial[which(Spatial$SiteCode == "03-Sep"),]$SiteCode <- "SEP3"

######################################
## Plot dataset with spatial groups ##
######################################
plotdata <- dplyr::filter(Spatial, HabitatP %in% c("Forest","Paramo","Pasture","Riparian","Fragment") | Ptrees == "X")
plotdata <- dplyr::filter(plotdata, is.na(FAge) | FAge > 20)
plotdata <- plotdata[,-c(1:2,17:35,37:49)]

## Pasture data changes
## If radius and plot data: remove plot data
## Additional rows for riparian and fragment data - separate sitecode
plotdata <- rbind(cbind(plotdata[which(plotdata$HabitatP %in% c("Forest","Paramo")),], "method" = "plot"),
                  cbind(plotdata %>% filter((Ptrees == "" | is.na(Ptrees)) & plotdata$HabitatP == "Pasture"), "method" = "plot"),
                  cbind(plotdata[which(plotdata$Ptrees == "X"),], "method" = "radius") %>% mutate(HabitatP = "Pasture", Plot = "", Size = NA),
                  cbind(plotdata[which(plotdata$HabitatP %in% c("Riparian","Fragment")),], "method" = "plot") %>% mutate(SiteCode = paste(SiteCode, "_F", sep = ""), Ptrees = ""))

Treedata[which(Treedata$SiteCode %in% gsub("_F","", dplyr::filter(plotdata, HabitatP %in% c("Riparian","Fragment"))$SiteCode)),]$SiteCode <- paste(Treedata[which(Treedata$SiteCode %in% gsub("_F","", dplyr::filter(plotdata, HabitatP %in% c("Riparian","Fragment"))$SiteCode)),]$SiteCode, "_F", sep="")
Treedata <- Treedata[-which(!Treedata$SiteCode %in% dplyr::filter(plotdata, method == "plot")$SiteCode),]
PTrees <- PTrees[-which(!PTrees$SiteCode %in% dplyr::filter(plotdata, method == "radius")$SiteCode),]

## Add site and cluster
plotdata$SiteNum <- as.numeric(as.factor(plotdata$SiteCode))
plotdata$ClusNum <- as.numeric(as.factor(plotdata$Cluster))

## Cluster points by distance
coords <- sp::SpatialPoints(cbind(plotdata$long, plotdata$lat))
distmat <- sp::spDists(cbind(plotdata$long, plotdata$lat),longlat = TRUE)
hc <- hclust(as.dist(distmat), "ward.D2")

# Find smallest within-area distance variance that ensures larger within-area maximum distance than smallest between-area minimum distance
# treeheight <- vector() ; distsep <- vector()
# for(i in 1:100){
#   treeheight[i] <- i
#   distsep[i] <- fpc::cluster.stats(as.dist(distmat), clustering = cutree(hc, h = i))$max.diameter - fpc::cluster.stats(as.dist(distmat), clustering = cutree(hc, h = i))$min.separation
# }
treeheight <- 65      ## Relatively random at the moment

  # Plot dendrogram / map
  plot(hc, labels = plotdata$SiteCode) ; abline(h = seq(0,50,10), col = "red") ; abline(h = c(15,25), col = "blue")
  raster::plot(coords, col = cutree(hc, h = treeheight), pch = 18)

# Add area numbers
plotdata$AreaNum <- cutree(hc, h = treeheight)

data.frame(plotdata %>% group_by(AreaNum) %>% summarise(unique(AreaCode)))
fpc::cluster.stats(as.dist(distmat), clustering = cutree(hc, h = treeheight))$max.diameter
fpc::cluster.stats(as.dist(distmat), clustering = cutree(hc, h = treeheight))$min.separation

## Add region numbers - and slight changes to regions FIX!!
plotdata$region[which(plotdata$region == "oriental na")] <- "oriental western slope"
plotdata$region[which(plotdata$AreaCode == "CQ")] <- "oriental eastern slope"
plotdata$region[which(plotdata$AreaCode == "CH")] <- "oriental eastern slope"
plotdata$RegionNum <- as.numeric(as.factor(plotdata$region))

######################
## Prepare datasets ##
######################
alltrees <- rbind(cbind(Treedata, "Dist" = NA), PTrees)

# Remove palms and ferns
alltrees <- dplyr::filter(alltrees, !Species %in% c("Palm","Fern"))

# Add points without trees for prediction of average WSG
AddPoints <- unique(plotdata[which(!plotdata$SiteCode %in% unique(alltrees$SiteCode)),]$SiteCode)
#AddPoints <- unique(dplyr::filter(plotdata, Plot == "X" & !SiteCode %in% unique(alltrees$SiteCode))$SiteCode)
alltrees[(nrow(alltrees)+1):(nrow(alltrees)+length(AddPoints)),]$SiteCode <- AddPoints

# Combine tree data with plot data
alltrees <- dplyr::filter(alltrees, SiteCode %in% plotdata$SiteCode)
alltrees <- dplyr::left_join(alltrees, plotdata, by = "SiteCode")

# Remove species with less than two individuals or no cores
coredSp <- alltrees %>% mutate(cored = WSG/WSG)
coredSp <- coredSp %>% group_by(Species) %>% summarise(total = n(), cored = sum(cored,na.rm=T))
coredSp <- dplyr::filter(coredSp, cored > 1)
coredSp <- dplyr::filter(coredSp, total > cored)

## Add "speciesid" column (1 = identified, 0 = not identified)
alltrees[which(!alltrees$Species %in% coredSp$Species),]$Species <- NA
alltrees[which(is.na(alltrees$Species)),]$Species <- "nospec"
alltrees$speciesid <- ifelse(alltrees$Species == "nospec", 0, 1)

## Calculate tree volume
alltrees <- dplyr::mutate(alltrees, vol = case_when(
  AreaCode == "CC" ~ exp(2.183 - 0.665 * log(DBH_used) + 0.892 * (log(DBH_used)^2) - 0.097 * (log(DBH_used)^3)),  # BMAlvII2.3 - Dry forest
  Species == "Polylepis" ~ 0.0694*(DBH_used^2.35996),                                                             # Polylepis equation (Espinoza & Quispe 2005)
  TRUE ~ exp(2.789 - 1.414 * log(DBH_used) + 1.178 * (log(DBH_used)^2) - 0.118 * (log(DBH_used)^3)),              # BMAlvII2.4 - All forest types
))
alltrees <- left_join(alltrees, alltrees %>% group_by(SiteCode) %>% summarise(totVOL = sum(vol)))
alltrees$VOLw <- alltrees$vol / alltrees$totVOL
alltrees[which(is.na(alltrees$VOLw)),]$VOLw <- 1

## Add to plotdata
plotdata <- left_join(plotdata,
                      alltrees %>% filter(method == "radius") %>% group_by(SiteCode) %>% summarise(MaxDist = max(Dist), MaxNum = max(TreeN)) %>% mutate(dist = case_when(
                        MaxNum > 9 ~ MaxDist,
                        MaxNum < 10 ~ 50)) %>% mutate(Psize = case_when(
                          is.na(dist) ~ pi * (50^2),
                          !is.na(dist) ~ pi * (dist^2))) %>% select(c(SiteCode, Psize)))
plotdata <- mutate(plotdata, Size = case_when(
                      is.na(Psize) ~ as.numeric(Size),
                      !is.na(Psize) ~ Psize
                    )) %>% select(-Psize)

plotdata <- left_join(plotdata, alltrees %>% group_by(SiteCode) %>% summarise(nQuercus = sum(Species == "Quercus Humboldtii"), nSpec = sum(!Species %in% c("Quercus Humboldtii","nospec")), ntree = sum(TreeN/TreeN, na.rm=T), ncore = sum(WSG/WSG,na.rm=T), Treevol = sum(vol, na.rm=T)))

plotdata$Treevol <- (plotdata$Treevol / plotdata$Size) * 10000

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
rm(coords, coredSp, distmat, hc, Spatial, treeheight, AddPoints)

####################################
## WSG model for individual trees ## 
####################################
#PTrees$HabitatP <- "Pasture"
#alltrees <- rbind(cbind(Treedata[,c("SiteNum","WSG","Species","speciesid","HabitatP","ClusNum","AreaNum","VOLw")], "method" = "plot"),
                 # cbind(PTrees[,c("SiteNum","WSG","Species","speciesid","HabitatP","ClusNum","AreaNum","VOLw")], "method" = "radius"))

moddat <- list("Nobs" = length(alltrees$WSG),                        # Number of observations
               "WSG" = alltrees$WSG,                                 # WSG data
               "Species" = as.numeric(as.factor(alltrees$Species)),  # Vector of species
               "nspec" = length(unique(alltrees$Species)),           # Number of species
               "speciesid" = alltrees$speciesid,                     # 0/1 - species determined or not
               "habitat" = as.numeric(as.factor(alltrees$HabitatP)), # Vector of habitat     
               "nhabitat" = length(unique(alltrees$HabitatP)),       # Number of habitats
               "SiteNum" = alltrees$SiteNum,
               "ClusNumS" = data.frame(plotdata %>% group_by(SiteNum) %>% summarise(first(ClusNum)))[,2],
               "AreaNumC" = data.frame(plotdata %>% group_by(ClusNum) %>% summarise(first(AreaNum)))[,2],
               "nsites" = length(unique(plotdata$SiteNum)),
               "ncluster" = length(unique(plotdata$ClusNum)),
               "narea" = length(unique(plotdata$AreaNum)))

modinits <- function(){list(alpha = rnorm(1,0.5,0.1), muspec = rnorm(1,0,0.001), sdspec = rnorm(1,0.036,0.01),
                            sdsite = rnorm(1,0.03,0.01), sdcluster = rnorm(1,0.017,0.01), sdarea = rnorm(1,0.05,0.01))}

n.chains <- 4
n.burnin <- 10000
n.iter <- n.burnin + 50000
n.thin <- 50
n.cores <- n.chains

start.time <- Sys.time()
WSGtreemod <- jagsUI::jags(model.file = "JAGS//CarbDiv_WSGmod.txt", data = moddat, inits = modinits,
                           parameters.to.save = c("alpha","var1","var0","muspec","varspec","vararea","varsite","varcluster","loglik","WSG.pred","res","b_area","b_cluster","b_site","b_habitat","mu"),
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                           parallel = T, n.cores = n.cores, codaOnly = c("WSG.pred","loglik","res","mu"))
time.used.WSG <- Sys.time()-start.time
#saveRDS(WSGtreemod, "Output\\CarbonDiversity\\WSGtreemod.RDS")
#WSGtreemod <- readRDS("Output\\CarbonDiversity\\WSGtreemod.RDS")

#########################
## Grass biomass model ## 
#########################
#Grassdata[which(Grassdata$HabitatP == ""),]$HabitatP <- "Pasture"
Grassdata[which(Grassdata$HabitatP == "Riparian"),]$HabitatP <- "Forest"
Grassdata[which(Grassdata$HabitatP == "Fragment"),]$HabitatP <- "Forest"

moddat <- list("Nobs" = length(Grassdata$GrassAGB),                   # Number of observations
               "Grass" = Grassdata$GrassAGB,                          # Grass biomass data
               "habitat" = as.numeric(as.factor(Grassdata$HabitatP)), # Vector of habitat     
               "nhabitat" = length(unique(Grassdata$HabitatP)),       # Number of habitats
               "AreaNum" = Grassdata$AreaNum,
               "narea" = length(unique(plotdata$AreaNum)))

n.chains <- 4
n.burnin <- 10000
n.iter <- n.burnin + 50000
n.thin <- 50
n.cores <- n.chains

modinits <- function(){list(b_area = rnorm(length(unique(plotdata$AreaNum)),1,0.1))}

start.time <- Sys.time()
GrassAGBmod <- jagsUI::jags(model.file = "JAGS//CarbDiv_grassBM.txt", data = moddat, inits = modinits,
                            parameters.to.save = c("alpha","var1","tau","vararea","Grass.pred","res","b_habitat","b_area","mu","loglik"),
                            n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                            parallel = T, n.cores = n.cores, codaOnly = c("Grass.pred","loglik","res","mu"))
time.used.grass <- Sys.time()-start.time
#saveRDS(GrassAGBmod, "Output\\CarbonDiversity\\GrassAGBmod.RDS")
#GrassAGBmod <- readRDS("Output\\CarbonDiversity\\GrassAGBmod.RDS")

#############################################
## Predict WSG and tree AGB across dataset ##
#############################################
## Predict WSG for individual trees
WSGdraws <- WSGtreemod$sims.list$WSG.pred
WSGdraws[, which(!is.na(Treedata$WSG))] <- matrix(data = rep(Treedata$WSG[which(!is.na(Treedata$WSG))], nrow(WSGdraws)), nrow = nrow(WSGdraws), ncol = length(Treedata$WSG[which(!is.na(Treedata$WSG))]), byrow = TRUE)

WSGdraws_volw <- as.data.frame(t(WSGdraws) * alltrees$VOLw)

## Calculate plot WSG and AGB
WSGdraws_plot <- data.frame(WSGdraws_volw) %>% dplyr::mutate(SiteNum = alltrees$SiteNum) %>% reshape2::melt(id.vars = "SiteNum") %>% reshape2::dcast(SiteNum ~ variable, value.var="value", fun.aggregate = sum)

Treedraws_plot <- dplyr::left_join(plotdata[which(!is.na(plotdata$Treevol)), c("SiteNum","Treevol")], WSGdraws_plot)
Treedraws_plot[,-c(1:2)] <- (Treedraws_plot$Treevol * Treedraws_plot[,-c(1:2)]) / 1000
Treedraws_plot <- Treedraws_plot[,-2]

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
Deaddraws_plot <- Decayvolha$TotVol * WSGdraws_plot[,-1]
Deaddraws_plot <- cbind(SiteNum = Decayvolha$SiteNum, Deaddraws_plot)

###########
## Grass ##
###########
Grassdata <- Grassdata[order(Grassdata$SiteNum),]
Grassdraws <- GrassAGBmod$sims.list$Grass.pred
Grassdraws[, which(!is.na(Grassdata$GrassAGB))] <- matrix(data = rep(Grassdata$GrassAGB[which(!is.na(Grassdata$GrassAGB))], nrow(Grassdraws)), nrow = nrow(Grassdraws), ncol = length(Grassdata$GrassAGB[which(!is.na(Grassdata$GrassAGB))]), byrow = TRUE)
Grassdraws <- cbind(Grassdata$SiteNum,as.data.frame(t(Grassdraws)))

############
## Paramo ##
############
Paramodraws <- as.data.frame(Paramodata %>% group_by(SiteNum) %>% summarise(ParamoAGB = (sum(AGB, na.rm=T) / first(Size)) * 10000 / 1000))

##################
## Save objects ##
##################
AGBdraws_plot <- replace(Treedraws_plot, is.na(Treedraws_plot), 0)[,-1] + 
                 replace(Deaddraws_plot, is.na(Deaddraws_plot), 0)[,-1] + 
                 replace(Grassdraws, is.na(Grassdraws), 0)[,-1] + 
                 as.data.frame(matrix(rep(Paramodraws$ParamoAGB, ncol(Treedraws_plot[,-1])), nrow = nrow(Treedraws_plot[,-1]), ncol = ncol(Treedraws_plot[,-1])))
AGBdraws_plot <- cbind(SiteNum = Treedraws_plot[,1], AGBdraws_plot)
AGBlist <- list(AGBdraws = AGBdraws_plot, Treedraws = Treedraws_plot, Deaddraws = Deaddraws_plot, Grassdraws = Grassdraws, Paramodraws = Paramodraws, plotdata = plotdata)

saveRDS(AGBlist, "Output\\CarbonDiversity\\AGBlist.RDS")

