#################
################# Carbon and biodiversity relationships
#################
if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}

setwd(dir.path)

library(dplyr)
library(brms)
library(ggplot2)

##########################
## Read and filter data ##
##########################
AGBlist_clus <- readRDS("Output\\CarbDiv\\AGBlist_clus.RDS")
#AGBlist <- readRDS("Output\\CarbonDiversity\\AGBlist.RDS")

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

# Remove point without bird data
AGBlist_clus$plotdata <- AGBlist_clus$plotdata[which(AGBlist_clus$plotdata$SiteCode %in% z_info$point),]
AGBlist_clus$AGBdraws <- AGBlist_clus$AGBdraws[which(AGBlist_clus$AGBdraws$ClusNum %in% unique(AGBlist_clus$plotdata$ClusNum)),]
AGBlist_clus$clusdata <- AGBlist_clus$clusdata[which(AGBlist_clus$clusdata$ClusNum %in% unique(AGBlist_clus$plotdata$ClusNum)),]

# # Add biogeography
# biogeo <- as.data.frame(AGBlist_clus$plotdata %>% group_by(AreaCode) %>% summarise(region = first(region)))
# biogeo <- biogeo %>% mutate(biogeo = case_when(
#           region == "central" ~ "Central Andes",
#           region == "central-eastern" ~ "Central Andes",
#           region == "santa_marta" ~ "Santa Marta",
#           region == "amazonia" ~ "Amazonia",
#           region == "oriental eastern slope" ~ "East Andes",
#           region == "oriental western slope" ~ "East Andes",
#           region == "occidental" ~ "West Andes"))
# biogeo <- biogeo %>% mutate(biogeo = case_when(
#           AreaCode == "EP" ~ "Magdalena valley",
#           AreaCode == "PP" ~ "Magdalena valley",
#           TRUE ~ biogeo
# ))

# AGBlist_clus$plotdata <- left_join(AGBlist_clus$plotdata, biogeo[,c("AreaCode","biogeo")])
# AGBlist_clus$clusdata <- left_join(AGBlist_clus$clusdata, biogeo[,c("AreaCode","biogeo")])

####################
## Draw bird data ##
####################
# Additional bird data
birdlife_extents <- read.csv("C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD\\Diversity_code_data\\birdlife_scraper\\birdlife_extents.csv")
birdlife_extents <- birdlife_extents[which(!is.na(birdlife_extents$extent_breeding_resident)),]
birdlife_extents$invrange <- 1/scale(birdlife_extents$extent_breeding_resident, center = F)
birdlife_extents$X <- gsub("[ ]", "_", birdlife_extents$X)

redlist_status <- read.csv("C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD\\Diversity_code_data\\birdlife_scraper\\redlist_status.csv")
redlist_status <- dplyr::filter(redlist_status, !V1  %in% c("Extinct","Extinct in the Wild","Data Deficient"))
redlist_status$X <- gsub("[ ]", "_", redlist_status$X)

# Remove point without carbon data
z_info <- dplyr::filter(z_info, point %in% AGBlist_clus$plotdata$SiteCode)

# Draw bird data
Zdraws <- get_specZ_draws_expand(10, draws, z_info, pointid = AGBlist_clus$plotdata[,c("SiteCode","SiteNum")])
Zclus <- aggregate(Zdraws, list(Zdraws$species, left_join(Zdraws, AGBlist_clus$plotdata[,c("SiteNum","ClusNum")])$ClusNum), max)[,-c(1,3)]

# Restricted-range species
Zclus_rr <- dplyr::filter(Zclus, species %in% dplyr::filter(birdlife_extents, (!migratory_status  %in% c("full migrant","unknown")) & extent_breeding_resident < 50000)$X)

# Redlist status
Zclus_endangered <- dplyr::filter(Zclus, species %in% filter(redlist_status, V1 %in% c("Critically Endangered","Endangered"))$X)

# Inverse range weighted species richness
Zclus_rangew <- dplyr::filter(Zclus, species %in% dplyr::filter(birdlife_extents, !migratory_status  %in% c("full migrant","unknown"))$X)
invrange <- as.numeric(dplyr::left_join(Zclus_rangew, dplyr::filter(birdlife_extents, !migratory_status  %in% c("full migrant","unknown"))[,c("X","invrange")], by = c("species" = "X"))$invrange)
Zclus_rangew[,3:12] <- Zclus_rangew[,3:12] * invrange

###########################
## Draw alpha model data ##
###########################
AGBlist_clus_temp <- AGBlist_clus_temp
Zclus_temp <- Zclus
AGBlist_clus$clusdata <- AGBlist_clus$clusdata[which(!AGBlist_clus$clusdata$biogeo == "East Andes"),]
AGBlist_clus$AGBdraws <- AGBlist_clus$AGBdraws[which(AGBlist_clus$AGBdraws$ClusNum %in% AGBlist_clus$clusdata$ClusNum),]
Zclus <- Zclus[which(Zclus$Group.2 %in% AGBlist_clus$clusdata$ClusNum),]

# Model data
moddat <- list()
moddat$richness <- list("n_cluster" = length(unique(Zclus[,1])),
                        "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                        "niter" = dim(Zclus[,-c(1,2)])[2],
                        "div" = rowMeans(apply(simplify2array(by(Zclus[,-c(1:2)], Zclus$species, as.matrix)),2,rowSums)), ### Sorted by cluster number
                        "carbon_cluster" = rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]),
                        "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)))
moddat$ranger <- list("n_cluster" = length(unique(Zclus_rr[,1])),
                      "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                      "niter" = dim(Zclus_rr[,-c(1,2)])[2],
                      "div" = rowMeans(apply(simplify2array(by(Zclus_rr[,-c(1:2)], Zclus_rr$species, as.matrix)),2,rowSums)), ### Sorted by cluster number
                      "carbon_cluster" = colMeans(AGBlist_clus$AGBdraws),
                      "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata$biogeo)))
moddat$endangered <- list("n_cluster" = length(unique(Zclus_endangered[,1])),
                          "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                          "niter" = dim(Zclus_endangered[,-c(1,2)])[2],
                          "div" = rowMeans(apply(simplify2array(by(Zclus_endangered[,-c(1:2)], Zclus_endangered$species, as.matrix)),2,rowSums)), ### Sorted by cluster number
                          "carbon_cluster" = colMeans(AGBlist_clus$AGBdraws),
                          "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata$biogeo)))

### MODEL 1
mod <- "data{
              int<lower=1> n_cluster;            // number of clusters
              vector[n_cluster] carbon_cluster;  // carbon at each cluster
              vector[n_cluster] div;             // diversity metric at each cluster
              //int<lower=1> biogeo[n_cluster];  // biogeographic region of each cluster
              //int<lower=1> n_biogeo;           // number of biogeographic regions
            }
            
            parameters{
              real alpha;
              real<lower=0> sigma_cluster;
              real beta_carbon;
              real beta2;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + beta_carbon * carbon_cluster + beta2 * (carbon_cluster .^ 2);
            }
            
            model{
              // priors
              alpha ~ normal(0, 50);
              beta_carbon ~ normal(0, 10);
              beta2 ~ normal(0, 10);
              sigma_cluster ~ normal(0,20);
              // likelihood
              div ~ normal(cluster_mean, sigma_cluster);
            }"

testmod_div <- rstan::stan(model_code = mod, data = moddat$richness, chains = 3, iter = 2000)

pairs(testmod_div, pars=c("alpha","beta_carbon","beta2","sigma_cluster"))

test <- t(matrix(rep(rstan::extract(testmod_div, "alpha")[[1]], length(moddat$richness$carbon_cluster)), ncol = length(moddat$richness$carbon_cluster))) + 
  t(matrix(rep(rstan::extract(testmod_div,"beta_carbon")[[1]], length(moddat$richness$carbon_cluster)), ncol = length(moddat$richness$carbon_cluster))) * moddat$richness$carbon_cluster + 
  t(matrix(rep(rstan::extract(testmod_div,"beta2")[[1]], length(moddat$richness$carbon_cluster)), ncol = length(moddat$richness$carbon_cluster))) * moddat$richness$carbon_cluster ^ 2
plot(moddat$richness$div ~ moddat$richness$carbon_cluster, col = moddat$richness$biogeo, pch = 16) ; points(rowMeans(test) ~ moddat$richness$carbon_cluster)

plot(moddat$richness$div ~ AGBlist_clus$clusdata$meanAGB[order(AGBlist_clus$clusdata$ClusNum)], col = moddat$richness$biogeo, pch = 16) ; points(rowMeans(test) ~ moddat$richness$carbon_cluster)

moddat$richness <- list("n_cluster" = length(unique(Zclus[,1])),
                        "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                        "niter" = dim(Zclus[,-c(1,2)])[2],
                        "div" = rowMeans(apply(simplify2array(by(Zclus[,-c(1:2)], Zclus$species, as.matrix)),2,rowSums)),              ### Sorted by cluster number
                        "carbon_cluster" = colMeans(AGBlist_clus$AGBdraws)[order(AGBlist_clus$clusdata$ClusNum)],
                        "n_knots" = 3)

mod <- "data{
              int<lower=1> n_cluster;            // number of clusters
              vector[n_cluster] carbon_cluster;  // carbon at each cluster
              vector[n_cluster] div;             // diversity metric at each cluster
              int<lower=1> n_knots;
            }
            
            parameters{
              real alpha;
              real<lower=0> sigma_cluster;
              vector[n_knots] beta;
              
              real<lower=0> kprob;
              vector[n_knots] kz;
              
              real<lower=0> add_mu[n_knots];
              real<lower=0> add_sd[n_knots];
              real<lower=0> theta_add[n_knots];
            }
            
            transformed parameters{
              real<lower=0> theta[n_knots];
              
              theta[1] = 0;
              for(i in 2:n_knots){
                theta[i] <- theta[i-1] + theta_add[i];
              }
              
              matrix[n_cluster, n_knots] carbspline;
              
              for(i in 1:n_knots){
                for(k in 1:n_cluster){
                  carbspline[i,k] = carbon_cluster[i] - theta[k] < 0 ? 0 : carbon_cluster[k] - theta[k];
                }
              }
              
              vector[n_cluster] carbfunc = carbspline * (beta .* kz);
              vector[n_cluster] cluster_mean = alpha + carbfunc;
            }
            
            model{
              // priors
              alpha ~ normal(0, 50);
              beta ~ normal(0, 10);
              sigma_cluster ~ normal(0,20);
              
              add_mu ~ normal(0,50);
              add_sd ~ normal(0,20);
              theta_add ~ normal(add_mu,add_sd);
              
              kprob ~ uniform(0,1);
              kz ~ bernoulli(kprob);
              
              // likelihood
              div ~ normal(cluster_mean, sigma_cluster);
            }"

testmod_div <- rstan::stan(model_code = mod, data = moddat$richness, chains = 3, iter = 2000)

mod <- "data{
              int<lower=1> n_knots;
            }
            
            parameters{

              real<lower=0> kprob;
              vector[n_knots] kz;

            }
            
            transformed parameters{
            }
            
            model{
              // priors

              kprob ~ uniform(0,1);
              kz ~ bernoulli(1);

            }"

#######################################
## Alpha diversity vs. carbon models ##
#######################################
#modinits <- function(){list(alpha = )

moddat <- list()
moddat$richness <- list("nsite" = length(unique(Zclus[,1])),
                        "niter" = dim(Zclus[,-c(1,2)])[2],
                        "div" = rowMeans(apply(simplify2array(by(Zclus[,-c(1:2)], Zclus$species, as.matrix)),2,rowSums)),              ### Sorted by cluster number
                        "carb" = colMeans(AGBlist_clus$AGBdraws),
                        "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata$biogeo)))

n.chains <- 4
n.burnin <- 50000
n.iter <- n.burnin + 2000
n.thin <- 1
n.cores <- n.chains

start.time <- Sys.time()
divmod_richness <- jagsUI::jags(model.file = "JAGS//CarbDiv_linear.txt", data = moddat$richness, #inits = modinits,
                              parameters.to.save = c("alpha","b1","kz","nz","kprob","theta","sd","mu","carbspline","carbfunc","div_pred"),
                              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                              parallel = T, n.cores = n.cores, codaOnly = c("loglik"))
t.splinemod <- Sys.time()-start.time

start.time <- Sys.time()
divmod_ranger <- jagsUI::jags(model.file = "JAGS//CarbDiv_spline.txt", data = moddat$ranger, #inits = modinits,
                              parameters.to.save = c("alpha","kz","nz","kprob","theta","b","sd","mu","carbspline","carbfunc","div_pred"),
                              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                              parallel = T, n.cores = n.cores, codaOnly = c("loglik"))
t.splinemod <- Sys.time()-start.time

start.time <- Sys.time()
divmod_endangered <- jagsUI::jags(model.file = "JAGS//CarbDiv_spline.txt", data = moddat$endangered, #inits = modinits,
                              parameters.to.save = c("alpha","kz","nz","kprob","theta","b","sd","mu","carbspline","carbfunc","div_pred"),
                              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                              parallel = T, n.cores = n.cores, codaOnly = c("loglik"))
t.splinemod <- Sys.time()-start.time

#######
plot(rowMeans(moddat$richness$div[,-1]) ~ moddat$richness$carb, col = as.factor(AGBlist_clus$clusdata$biogeo), pch = 18)
points(divmod_richness$mean$div_pred ~ moddat$richness$carb, col = "blue", pch = 18)

plot(rowMeans(moddat$ranger$div[,-1]) ~ moddat$ranger$carb, col = as.factor(AGBlist_clus$clusdata$Habitat), pch = 18)
points(divmod_ranger$mean$div_pred ~ moddat$ranger$carb, col = "blue", pch = 18)

plot(rowMeans(moddat$endangered$div[,-1]) ~ moddat$endangered$carb, col = as.factor(AGBlist_clus$clusdata$RegionNum), pch = 18)
points(divmod_endangered$mean$div_pred ~ moddat$endangered$carb, col = "blue", pch = 18)





























# Choose spatial index
Ztemp <- Zdraws
#spatIndex <- left_join(Ztemp, AGBlist_birds$plotdata[,c("SiteNum","AreaNum")])$AreaNum ; carb <- carbmod_area$mean$carb.area[-which(carbmod_area$mean$carb.area == 0)]
#spatIndex <- left_join(Ztemp, AGBlist_birds$plotdata[,c("SiteNum","ClusNum")])$ClusNum ; carb <- carbmod_clus$mean$carb.clus[-which(carbmod_clus$mean$carb.clus == 0)]
spatIndex <- Ztemp$SiteNum ; carb <- rowMeans(AGBlist_birds$AGBdraws[order(AGBlist_birds$AGBdraws$SiteNum),-1])
spatIndex <- left_join(Ztemp, AGBlist_birds$plotdata[,c("SiteNum","ClusNum")])$ClusNum ; carb <- aggregate(rowMeans(AGBlist_birds$AGBdraws[order(AGBlist_birds$AGBdraws$SiteNum),-1]), list(AGBlist_birds$plotdata$ClusNum), mean)[,2]

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

plot(rowMeans(moddat$div[,-1]) ~ moddat$carb, col = as.factor(AGBlist_birds$plotdata$RegionNum))
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

plot(rowMeans(moddat$div[,-1]) ~ moddat$carb, col = as.factor(AGBlist_birds$plotdata$RegionNum), pch = 18)
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
divmod_polynomial <- jagsUI::jags(model.file = "JAGS//CarbDiv_polynomial.txt", data = moddat$richness, #inits = modinits,
                                  parameters.to.save = c("alpha","b1","b2","sd","mu","div_pred"),
                                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                                  parallel = T, n.cores = n.cores, codaOnly = c("loglik"))
t.polymod <- Sys.time()-start.time

plot(rowMeans(moddat$richness$div[,-1]) ~ moddat$richness$carb, col = as.factor(AGBlist_clus$plotdata$AreaNum), pch = 16)
points(divmod_polynomial$mean$div_pred ~ moddat$richness$carb, col = "blue")



##########################
## Draw beta model data ##
##########################
## Remove upper/lower triangle - don't include self-distances in vector - permute

test <- simplify2array(by(Zdraws[,-c(1:2)], Zdraws$species, as.matrix))

Zclus <- aggregate(Zdraws, list(Zdraws$species, left_join(Zdraws, AGBlist_clus$plotdata[,c("SiteNum","ClusNum")])$ClusNum), max)
test <- simplify2array(by(Zclus[,-c(1:4)], Zclus$species, as.matrix))

bdiv_sor <- sapply(1:ncol(Zdraws[,-c(1:2)]), function(x) as.matrix(vegan::betadiver(cbind(test[,x,]), "sor")), simplify = "array")
carb_diff <- as.matrix(dist(AGBlist_clus$AGBdraws[2,], diag = NA))
geo_diff <- as.matrix(dist(cbind(AGBlist_clus$clusdata$long,AGBlist_clus$clusdata$lat)))



rowMeans(test)

ape::mantel.test(bdiv_sor[,,1],carb_diff)

vegan::mantel.partial(bdiv_sor[[1]], carb_diff, geo_diff)
vegan::mantel.partial(rowMeans(bdiv_sor,dims = 2), carb_diff, geo_diff)

vegan::mantel.correlog(bdiv_sor[,,1], geo_diff)

vegan::cca(cbind(test[,1,]), AGBlist_clus$AGBdraws[2,])

vegan::cca(cbind(test[,1,]), AGBlist_clus$AGBdraws[2,], geo_diff)

test2 <- cbind(test[,1,])

plot(rowMeans(bdiv_sor,dims = 2) ~ carb_diff)
points(rowMeans(test) ~ bmod_data$carb, col = "red")
abline(lm(y ~ x, data.frame(x = as.numeric(carb_diff), y = as.numeric(rowMeans(bdiv_sor,dims = 2))), col = "blue"))

ggplot(data = data.frame(x = as.numeric(carb_diff), y = as.numeric(as.matrix(bdiv_sor[,,1]))), aes(x=x, y=y)) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "lm", formula = y ~ poly(x,3)) + 
  theme_bw()

ggplot(data = data.frame(x = as.numeric(geo_diff), y = as.numeric(as.matrix(bdiv_sor[,,1]))), aes(x=x, y=y)) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "lm", formula = y ~ poly(x,3)) + 
  theme_bw()

ggplot(data = data.frame(x = as.numeric(geo_diff), y = as.numeric(carb_diff)), aes(x=x, y=y)) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "lm", formula = y ~ poly(x,3)) + 
  theme_bw()



test2 <- lapply(1:ncol(Zclus[,-c(1:4)]), function(i) cbind(test[,i,]))
test3 <- data.frame(SiteNum = AGBlist_clus$clusdata$ClusNum, test2[[1]])

hei <- betapart::beta.pair(test2[[4]])$beta.sim
test3 <- data.frame(SiteNum = AGBlist_clus$clusdata$ClusNum, as.matrix(hei))

sitepairs <- gdm::formatsitepair(bioData = test3, bioFormat = 3, dist = "clark", siteColumn = "SiteNum", XColumn = "long", YColumn = "lat",
                    predData = data.frame(SiteNum = AGBlist_clus$clusdata$ClusNum, long = AGBlist_clus$clusdata$long, lat = AGBlist_clus$clusdata$lat, AGB = AGBlist_clus$AGBdraws[2,]))
plot(gdm::gdm(sitepairs, geo = T, splines = c(3,3)))
test4 <- gdm::gdm(sitepairs, geo = F, splines = 10)
test5 <- gdm::gdm(sitepairs, geo = T, splines = c(3,3))

summary(lm(rowMeans(AGBlist_clus$AGBdraws[,-1]) ~ AGBlist_clus$clusdata$ALOSelev))

plot(rowMeans(AGBlist$AGBdraws[,-1]) ~ AGBlist$plotdata$ALOSelev)
abline(lm(rowMeans(AGBlist$AGBdraws[,-1]) ~ AGBlist$plotdata$ALOSelev))








