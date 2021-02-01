#################
################# Carbon and biodiversity relationships
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
AGBlist <- readRDS("Output\\CarbonDiversity\\AGBlist.RDS")

#v5 <- cmdstanr::read_cmdstan_csv("occupancy_v5_threads-202012282018-1-261afe.csv")
#draws <- posterior::as_draws_df(v5$post_warmup_draws[1:2000,,])
draws <- readRDS("Diversity_code_data\\v5draws.RDS")
bird_data <- readRDS("Diversity_code_data\\bird_stan_data4_package.RDS")
birds <- readRDS("Diversity_code_data\\birds.RDS")

# Generate z_info dataset
z_info <- data.frame(bird_data$data[8:43])
z_info$point <- birds$point
z_info$species <- birds$species

# Random species for test analyses
#z_info <- dplyr::filter(z_info, z_info$species %in% sample(unique(z_info$species), 25))

# Remove points
#AGBlist <- lapply(AGBlist, function(x) x[which(x$SiteNum %in% subset(AGBlist$plotdata, Dataset != "Chocó")$SiteNum),])
#AGBlist <- lapply(AGBlist, function(x) x[which(x$SiteNum %in% subset(AGBlist$plotdata, HabitatP != "Paramo")$SiteNum),])

# Remove point without carbon data
z_info <- dplyr::filter(z_info, point %in% AGBlist$plotdata$SiteCode)

# Remove point without bird data 
AGBlist_birds <- lapply(AGBlist, function(x) x[-which(x$SiteNum %in% AGBlist$plotdata[which(!AGBlist$plotdata$SiteCode %in% z_info$point),]$SiteNum),])

# Additional bird data
birdlife_extents <- read.csv("C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD\\Diversity_code_data\\birdlife_scraper\\birdlife_extents.csv")
birdlife_extents <- birdlife_extents[which(!is.na(birdlife_extents$extent_breeding_resident)),]
birdlife_extents$invrange <- 1/scale(birdlife_extents$extent_breeding_resident, center = F)
birdlife_extents$X <- gsub("[ ]", "_", birdlife_extents$X)

redlist_status <- read.csv("C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD\\Diversity_code_data\\birdlife_scraper\\redlist_status.csv")
redlist_status <- dplyr::filter(redlist_status, !V1  %in% c("Extinct","Extinct in the Wild","Data Deficient"))
redlist_status$X <- gsub("[ ]", "_", redlist_status$X)

################################
## Draw cluster carbon values ##
################################
spatial <- read.csv("Data\\vegetation\\Spatialdata.csv")
plotdata <- left_join(AGBlist_birds$plotdata, spatial)

carbdat <- list("nsite" = nrow(AGBlist_birds$plotdata),
                "nclus" = max(AGBlist_birds$plotdata$ClusNum),
                "carb" = as.numeric(rowMeans(AGBlist_birds$AGBdraws[order(AGBlist_birds$AGBdraws$SiteNum),-1])),
                "clus" = AGBlist_birds$plotdata[order(AGBlist_birds$plotdata$SiteNum),]$ClusNum,
                "mediancarb" = mutate_all(left_join(data.frame(cluster = c(seq(1, max(AGBlist_birds$plotdata$ClusNum), 1))),
                                                    aggregate(rowMeans(AGBlist_birds$AGBdraws[order(AGBlist_birds$AGBdraws$SiteNum),-1]), list(cluster = AGBlist_birds$plotdata[order(AGBlist_birds$plotdata$SiteNum), c("ClusNum")]), mean)),
                                                    ~replace(., is.na(.), 0))[,2],
                "size" = mutate_all(left_join(data.frame(cluster = c(seq(1, max(AGBlist_birds$plotdata$ClusNum), 1))),
                                              aggregate(AGBlist_birds$plotdata[order(AGBlist_birds$plotdata$SiteNum),"Size"], list(cluster = AGBlist_birds$plotdata[order(AGBlist_birds$plotdata$SiteNum), c("ClusNum")]), mean)),
                                    ~replace(., is.na(.), 0))[,2],
                "elev" = (data.frame(left_join(data.frame(ClusNum = c(1:max(AGBlist_birds$plotdata$ClusNum))), plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev))) %>% mutate(elev = case_when(is.na(elev) ~ 0, !is.na(elev) ~ elev)))$elev)/1000,
                "habitat" = as.numeric(as.factor(data.frame(left_join(data.frame(ClusNum = c(1:max(AGBlist_birds$plotdata$ClusNum))), plotdata %>% group_by(ClusNum) %>% summarise(habitat = first(HabitatP))) %>% mutate(habitat = case_when(is.na(habitat) ~ "Forest", !is.na(habitat) ~ habitat)))$habitat)),
                "area" = as.numeric(as.factor(data.frame(left_join(data.frame(ClusNum = c(1:max(AGBlist_birds$plotdata$ClusNum))), plotdata %>% group_by(ClusNum) %>% summarise(area = first(AreaNum))) %>% mutate(area = case_when(is.na(area) ~ as.integer(1), !is.na(area) ~ area)))$area)),
                "narea" = max(plotdata$AreaNum))

n.chains <- 2
n.burnin <- 5000
n.iter <- n.burnin + 5000
n.thin <- 1
n.cores <- n.chains

start.time <- Sys.time()
carbmod_clus <- jagsUI::jags(model.file = "JAGS//CarbDiv_plot2clus.txt", data = carbdat, #inits = modinits,
                             parameters.to.save = c("alpha","mc","bhab","belev","barea","carb.clus","tau.clus","mu.clus","d","sd.clus"),
                             n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                             parallel = T, n.cores = n.cores, codaOnly = c("loglik"))
t.carbmod <- Sys.time()-start.time

cAGB_field <- carbdat$mediancarb[-27]
cAGB_est <- carbmod_clus$mean$carb.clus[-27]
cAGB_hab <- carbdat$habitat[-27]
cAGB_elev <- carbdat$elev[-27]

plot(cAGB_est ~ cAGB_field)
abline(0,1) ; abline(lm(cAGB_est ~ cAGB_field), col = "red")

par(mfrow = c(1,2))
plot(cAGB_est ~ cAGB_elev, col = cAGB_hab, pch = 16)
abline(lm(cAGB_est[which(cAGB_hab == 1)] ~ cAGB_elev[which(cAGB_hab == 1)]), col = 1)
abline(lm(cAGB_est[which(cAGB_hab == 2)] ~ cAGB_elev[which(cAGB_hab == 2)]), col = 2)
abline(lm(cAGB_est[which(cAGB_hab == 3)] ~ cAGB_elev[which(cAGB_hab == 3)]), col = 3)
plot(cAGB_field ~ cAGB_elev, col = cAGB_hab, pch = 16)
abline(lm(cAGB_field[which(cAGB_hab == 1)] ~ cAGB_elev[which(cAGB_hab == 1)]), col = 1)
abline(lm(cAGB_field[which(cAGB_hab == 2)] ~ cAGB_elev[which(cAGB_hab == 2)]), col = 2)
abline(lm(cAGB_field[which(cAGB_hab == 3)] ~ cAGB_elev[which(cAGB_hab == 3)]), col = 3)

dataset <- as.numeric(as.factor(data.frame(left_join(data.frame(ClusNum = c(1:max(AGBlist_birds$plotdata$ClusNum))), plotdata %>% group_by(ClusNum) %>% summarise(dataset = first(Dataset))))$dataset))[-27]
plot(cAGB_est[which(cAGB_hab == 1)] ~ cAGB_field[which(cAGB_hab == 1)], pch = 16, col = dataset[which(cAGB_hab == 1)])
abline(0,1) ; abline(lm(cAGB_est[which(cAGB_hab == 1)] ~ cAGB_field[which(cAGB_hab == 1)]), col = "red")

plot(cAGB_est[which(cAGB_hab == 1)] ~ cAGB_elev[which(cAGB_hab == 1)], pch = 16, col = dataset[which(cAGB_hab == 1)])
abline(lm(cAGB_est[which(cAGB_hab == 1)] ~ cAGB_elev[which(cAGB_hab == 1)]), col = "red")
plot(cAGB_field[which(cAGB_hab == 1)] ~ cAGB_elev[which(cAGB_hab == 1)], pch = 16, col = dataset[which(cAGB_hab == 1)])
abline(lm(cAGB_field[which(cAGB_hab == 1)] ~ cAGB_elev[which(cAGB_hab == 1)]), col = "red")

#lines(smooth.spline(carbmod_clus$mean$carb.clus[which(carbdat$habitat == 1)] ~ carbdat$elev[which(carbdat$habitat == 1)], spar = 1.3), col = 1)
#lines(smooth.spline(carbmod_clus$mean$carb.clus[which(carbdat$habitat == 2)] ~ carbdat$elev[which(carbdat$habitat == 2)], spar = 1.3), col = 2)
#lines(smooth.spline(carbmod_clus$mean$carb.clus[which(carbdat$habitat == 3)] ~ carbdat$elev[which(carbdat$habitat == 3)], spar = 1.3), col = 3)

#############################
## Draw area carbon values ##
#############################
carbdat <- list("nsite" = nrow(AGBlist_birds$plotdata),
                "narea" = max(AGBlist_birds$plotdata$AreaNum),
                "carb" = as.numeric(rowMeans(AGBlist_birds$AGBdraws[order(AGBlist_birds$AGBdraws$SiteNum),-1])),
                "area" = AGBlist_birds$plotdata[order(AGBlist_birds$plotdata$SiteNum),]$AreaNum,
                "mediancarb" = mutate_all(left_join(data.frame(area = c(seq(1, max(AGBlist_birds$plotdata$AreaNum), 1))),
                                                    aggregate(rowMeans(AGBlist_birds$AGBdraws[order(AGBlist_birds$AGBdraws$SiteNum),-1]), list(area = AGBlist_birds$plotdata[order(AGBlist_birds$plotdata$SiteNum), c("AreaNum")]), median)),
                                          ~replace(., is.na(.), 0))[,2],
                "size" = mutate_all(left_join(data.frame(area = c(seq(1, max(AGBlist_birds$plotdata$AreaNum), 1))),
                                              aggregate(AGBlist_birds$plotdata[order(AGBlist_birds$plotdata$SiteNum),"Size"], list(area = AGBlist_birds$plotdata[order(AGBlist_birds$plotdata$SiteNum), c("AreaNum")]), mean)),
                                    ~replace(., is.na(.), 0))[,2])

n.chains <- 4
n.burnin <- 1000
n.iter <- n.burnin + 50000
n.thin <- 50
n.cores <- n.chains

start.time <- Sys.time()
carbmod_area <- jagsUI::jags(model.file = "JAGS//CarbDiv_plot2area.txt", data = carbdat, #inits = modinits,
                             parameters.to.save = c("carb.area","tau.area","mu.area","d","sd.area"),
                             n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                             parallel = T, n.cores = n.cores, codaOnly = c("loglik"))
t.carbmod <- Sys.time()-start.time

###########################
## Draw alpha model data ##
###########################
# Draw bird data
Zdraws <- get_specZ_draws_expand(10, draws, z_info, pointid = AGBlist_birds$plotdata[,c("SiteCode","SiteNum")])

# Restricted-range species
Zdraws_rr <- dplyr::filter(Zdraws, species %in% dplyr::filter(birdlife_extents, (!migratory_status  %in% c("full migrant","unknown")) & extent_breeding_resident < 50000)$X)

# Redlist status
Zdraws_endangered <- dplyr::filter(Zdraws, species %in% filter(redlist_status, V1 %in% c("Critically Endangered","Endangered"))$X)
Zdraws_endangered <- dplyr::filter(Zdraws, species %in% filter(redlist_status, V1 %in% c("Least Concern"))$X)

# Inverse range weighted species richness
Zdraws_rangew <- dplyr::filter(Zdraws, species %in% dplyr::filter(birdlife_extents, !migratory_status  %in% c("full migrant","unknown"))$X)
invrange <- as.numeric(dplyr::left_join(Zdraws_rangew, dplyr::filter(birdlife_extents, !migratory_status  %in% c("full migrant","unknown"))[,c("X","invrange")], by = c("species" = "X"))$invrange)
Zdraws_rangew[,3:12] <- Zdraws_rangew[,3:12] * invrange

# Choose spatial index
Ztemp <- Zdraws
spatIndex <- left_join(Ztemp, AGBlist_birds$plotdata[,c("SiteNum","AreaNum")])$AreaNum ; carb <- carbmod_area$mean$carb.area[-which(carbmod_area$mean$carb.area == 0)]
spatIndex <- left_join(Ztemp, AGBlist_birds$plotdata[,c("SiteNum","ClusNum")])$ClusNum ; carb <- carbmod_clus$mean$carb.clus[-which(carbmod_clus$mean$carb.clus == 0)]
spatIndex <- Ztemp$SiteNum ; carb <- rowMeans(AGBlist_birds$AGBdraws[order(AGBlist_birds$AGBdraws$SiteNum),-1])

moddat <- list("nsite" = length(unique(spatIndex)),
               "niter" = dim(Ztemp[,-c(1,2)])[2],
               "div" = divdraws(Ztemp, spatIndex),              ### Sorted by plot/cluster/area number
               "carb" = carb) 

#####################################
## Linear model diversity ~ carbon ##
#####################################
#modinits <- function(){list(alpha = )

n.chains <- 4
n.burnin <- 1000
n.iter <- n.burnin + 1000
n.thin <- 1
n.cores <- n.chains

start.time <- Sys.time()
divmod_linear <- jagsUI::jags(model.file = "JAGS//CarbDiv_linear.txt", data = moddat, #inits = modinits,
                              parameters.to.save = c("alpha","b1","sd","mu","div_pred"),
                              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                              parallel = T, n.cores = n.cores, codaOnly = c("loglik"))
t.linmod <- Sys.time()-start.time

plot(rowMeans(moddat$div[,-1]) ~ moddat$carb)#, col = as.factor(AGBlist_birds$plotdata$HabitatP))
points(divmod_linear$mean$div_pred ~ moddat$carb, col = "red")

#####################################
## Spline model diversity ~ carbon ##
#####################################
#modinits <- function(){list(alpha = )

n.chains <- 4
n.burnin <- 5000
n.iter <- n.burnin + 1000
n.thin <- 1
n.cores <- n.chains

start.time <- Sys.time()
divmod_spline <- jagsUI::jags(model.file = "JAGS//CarbDiv_spline.txt", data = moddat, #inits = modinits,
                              parameters.to.save = c("alpha","kz","nz","kprob","theta","b","sd","mu","carbspline","carbfunc","div_pred"),
                              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                              parallel = T, n.cores = n.cores, codaOnly = c("loglik"))
t.splinemod <- Sys.time()-start.time

plot(rowMeans(moddat$div[,-1]) ~ moddat$carb)#, col = as.factor(AGBlist_birds$plotdata$RegionNum), pch = 18)
points(divmod_spline$mean$div_pred ~ moddat$carb, col = "blue", pch = 18)

#########################################
## Polynomial model diversity ~ carbon ##
#########################################
#modinits <- function(){list(alpha = )

n.chains <- 4
n.burnin <- 1000
n.iter <- n.burnin + 1000
n.thin <- 1
n.cores <- n.chains

start.time <- Sys.time()
divmod_polynomial <- jagsUI::jags(model.file = "JAGS//CarbDiv_polynomial.txt", data = moddat, #inits = modinits,
                                  parameters.to.save = c("alpha","b1","b2","sd","mu","div_pred"),
                                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                                  parallel = T, n.cores = n.cores, codaOnly = c("loglik"))
t.polymod <- Sys.time()-start.time

plot(rowMeans(moddat$div[,-1]) ~ moddat$carb, col = as.factor(AGBlist_birds$plotdata$region))
points(divmod_polynomial$mean$div_pred ~ moddat$carb, col = "blue")

##########################
## Draw beta model data ##
##########################

## Remove upper/lower triangle - don't include self-distances in vector - permute

test <- simplify2array(by(Zdraws[,-c(1:2)], Zdraws$species, as.matrix))

bdiv_sor <- sapply(1:ncol(Zdraws[,-c(1:2)]), function(x) as.matrix(vegan::betadiver(cbind(test[,x,]), "sor")), simplify = "array")
carb_diff <- as.matrix(dist(AGBlist_birds$AGBdraws[,2], diag = NA))
geo_diff <- as.matrix(dist(cbind(AGBlist_birds$plotdata$long,AGBlist_birds$plotdata$lat)))



rowMeans(test)

ape::mantel.test(bdiv_sor[,,1],carb_diff)

vegan::mantel.partial(bdiv_sor[[1]], carb_diff, geo_diff)
vegan::mantel.partial(rowMeans(bdiv_sor,dims = 2), carb_diff, geo_diff)

vegan::mantel.correlog(bdiv_sor[[1]], geo_diff)

vegan::cca(cbind(test[,1,]), AGBlist_birds$AGBdraws[,2])

vegan::cca(cbind(test[,1,]), AGBlist_birds$AGBdraws[,2], geo_diff)

test2 <- cbind(test[,1,])

plot(rowMeans(bdiv_sor,dims = 2) ~ carb_diff)
points(rowMeans(test) ~ bmod_data$carb, col = "red")
abline(lm(y ~ x, data.frame(x = as.numeric(carb_diff), y = as.numeric(rowMeans(bdiv_sor,dims = 2))), col = "blue"))

ggplot(data = data.frame(x = as.numeric(carb_diff), y = as.numeric(as.matrix(bdiv_sor[[1]]))), aes(x=x, y=y)) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "lm", formula = y ~ poly(x,3)) + 
  theme_bw()

ggplot(data = data.frame(x = as.numeric(geo_diff), y = as.numeric(as.matrix(bdiv_sor[[1]]))), aes(x=x, y=y)) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "lm", formula = y ~ poly(x,3)) + 
  theme_bw()

ggplot(data = data.frame(x = as.numeric(geo_diff), y = as.numeric(carb_diff)), aes(x=x, y=y)) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "lm", formula = y ~ poly(x,3)) + 
  theme_bw()




test2 <- lapply(1:ncol(Zdraws[,-c(1:2)]), function(i) cbind(test[,i,]))
test3 <- data.frame(SiteNum = AGBlist_birds$plotdata$SiteNum, test2[[1]])

hei <- betapart::beta.pair(test2[[4]])$beta.sim
test3 <- data.frame(SiteNum = AGBlist_birds$plotdata$SiteNum, as.matrix(hei))

sitepairs <- gdm::formatsitepair(bioData = test3, bioFormat = 3, dist = "clark", siteColumn = "SiteNum", XColumn = "long", YColumn = "lat",
                    predData = data.frame(SiteNum = AGBlist_birds$AGBdraws[,1], long = AGBlist_birds$plotdata$long, lat = AGBlist_birds$plotdata$lat, AGB = AGBlist_birds$AGBdraws[,2]))
plot(gdm::gdm(sitepairs, geo = T, splines = c(3,3)))
test4 <- gdm::gdm(sitepairs, geo = F, splines = 10)
test5 <- gdm::gdm(sitepairs, geo = T, splines = c(3,3))

summary(lm(rowMeans(AGBlist$AGBdraws[,-1]) ~ AGBlist$plotdata$ALOSelev))

plot(rowMeans(AGBlist$AGBdraws[,-1]) ~ AGBlist$plotdata$ALOSelev)
abline(lm(rowMeans(AGBlist$AGBdraws[,-1]) ~ AGBlist$plotdata$ALOSelev))








