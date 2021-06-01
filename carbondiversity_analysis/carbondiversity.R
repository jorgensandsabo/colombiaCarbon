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

# Remove points
#AGBlist_clus$plotdata <- dplyr::filter(AGBlist_clus$plotdata, biogeo!= "West Andes")
#AGBlist_clus$plotdata <- dplyr::filter(AGBlist_clus$plotdata, HabitatP != "Paramo" & AreaCode != "TU")

# Remove point without bird data
AGBlist_clus$plotdata <- AGBlist_clus$plotdata[which(AGBlist_clus$plotdata$SiteCode %in% z_info$point),]
AGBlist_clus$AGBdraws <- AGBlist_clus$AGBdraws[which(AGBlist_clus$AGBdraws$ClusNum %in% unique(AGBlist_clus$plotdata$ClusNum)),]
AGBlist_clus$clusdata <- AGBlist_clus$clusdata[which(AGBlist_clus$clusdata$ClusNum %in% unique(AGBlist_clus$plotdata$ClusNum)),]

# Changes to biogeography
AGBlist_clus$clusdata  <- AGBlist_clus$clusdata %>%  mutate(biogeo = case_when(
                                            AreaCode %in% c("AL","SM","SE","JU","CH","MO") ~ "East - east slope",
                                            AreaCode %in% c("CQ","EP","TN","IG","GA","AT","TS","TA","VI","VC","BE","TU","SO","CC","RA","PP","HT") ~ "East - west slope",
                                            # AreaCode == "RA" ~ "Magdalena valley",
                                            # AreaCode == "HT" ~ "Magdalena valley",
                                            # AreaCode == "AL" ~ "East slope",
                                            # AreaCode == "SM" ~ "East slope",
                                            # 
                                            # 
                                            # Cluster %in% c("BE1","BE2","BE3","BE4","BE5","BE6","BE7","CH1","CH2","CH3","CH4","CH5","CH6","CH7","CH8","CH9","IG8",
                                            #                "TU1","TU2","TU3","TU4","TU5","TU6","TU7") ~ "East paramo",
                                            # Cluster %in% c("PN1","PU2") ~ "Central paramo",
                                            TRUE ~ biogeo))
AGBlist_clus$plotdata$biogeo <- NULL
AGBlist_clus$plotdata <- left_join(AGBlist_clus$plotdata, AGBlist_clus$clusdata %>% group_by(AreaCode) %>% summarise(biogeo = first(biogeo)))

# AGBlist_clus$plotdata <- left_join(AGBlist_clus$plotdata, biogeo[,c("AreaCode","biogeo")])
# AGBlist_clus$clusdata <- left_join(AGBlist_clus$clusdata, biogeo[,c("AreaCode","biogeo")])

####################
## Draw bird data ##
####################
# Additional bird data
birdlife_extents <- read.csv("C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD\\Diversity_code_data\\birdlife_scraper\\birdlife_extents.csv")
birdlife_extents <- birdlife_extents[which(!is.na(birdlife_extents$extent_breeding_resident)),]
birdlife_extents$X <- gsub("[ ]", "_", birdlife_extents$X)

redlist_status <- read.csv("C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD\\Diversity_code_data\\birdlife_scraper\\redlist_status.csv")
redlist_status <- dplyr::filter(redlist_status, !V1  %in% c("Extinct","Extinct in the Wild","Data Deficient"))
redlist_status$X <- gsub("[ ]", "_", redlist_status$X)

# Remove point without carbon data
z_info <- dplyr::filter(z_info, point %in% AGBlist_clus$plotdata$SiteCode)
#z_info_pasture <- z_info ; z_info_pasture$pasture <- 1

# Draw bird data - forest and equivalent pasture
Zdraws <- get_specZ_draws_expand(10, draws, z_info, pointid = AGBlist_clus$plotdata[,c("SiteCode","SiteNum")])
#Zdraws_pasture <- get_specZ_draws_expand(10, draws, z_info_pasture, pointid = AGBlist_clus$plotdata[,c("SiteCode","SiteNum")])

# Filter out full migrants
Zdraws <- dplyr::filter(Zdraws, species %in% dplyr::filter(birdlife_extents, (!migratory_status  %in% c("full migrant","unknown")))$X)
#Zdraws_pasture <- dplyr::filter(Zdraws_pasture, species %in% dplyr::filter(birdlife_extents, (!migratory_status  %in% c("full migrant","unknown")))$X)

# Weigh by pasture occurrence probability
#pastureprob <- cbind(Zdraws_pasture[,c(1:2)], pastureprob = rowMeans(Zdraws_pasture[,-c(1:2)]))
#Zdraws <- cbind(Zdraws[,c(1:2)],apply(Zdraws[,-c(1:2)], 2, function(x) x*(1-pastureprob$pastureprob)))

# Point to cluster occurrence
Zclus <- aggregate(Zdraws, list(Zdraws$species, left_join(Zdraws, AGBlist_clus$plotdata[,c("SiteNum","ClusNum")])$ClusNum), max)[,-c(1,3)]

# Restricted-range species
Zclus_rr <- dplyr::filter(Zclus, species %in% dplyr::filter(birdlife_extents, extent_breeding_resident < 50000)$X)

# Redlist status
Zclus_endangered <- dplyr::filter(Zclus, species %in% filter(redlist_status, V1 %in% c("Critically Endangered","Endangered","Vulnerable"))$X)

# Range weighted species richness
extents <- dplyr::filter(birdlife_extents, !migratory_status  %in% c("full migrant","unknown"))[,c("X","extent_breeding_resident")]
extents <- dplyr::filter(extents, X %in% Zclus$species)
extents$extent_breeding_resident[which(extents$extent_breeding_resident > 10^7.5)] <- 10^7.5
Zclus_rangew <- cbind(Zclus[,c(1:2)],Zclus[,-c(1:2)] * (1 - (extents$extent_breeding_resident / max(extents$extent_breeding_resident))))

###########################
## Draw alpha model data ##
###########################
# Model data
moddat <- list()
moddat$richness <- list("n_cluster" = length(unique(Zclus[,1])),
                        "n_area" = length(unique(AGBlist_clus$clusdata$AreaNum)),
                        "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                        "niter" = dim(Zclus[,-c(1,2)])[2],
                        "carbon_region" = aggregate(x = rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]), 
                                                    by = list(as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo))),
                                                    mean)$x,
                        "div" = rowMeans(apply(simplify2array(by(Zclus[,-c(1:2)], Zclus$species, as.matrix)),2,rowSums)), ### Sorted by cluster number
                        "carbon_cluster" = rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]),
                        "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)),
                        "area" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$AreaNum)),
                        "biogeo_area" = as.numeric(as.factor(AGBlist_clus$clusdata %>% group_by(AreaNum) %>% summarise(reg = first(bioNum)) %>% pull(reg))))
moddat$ranger <- list("n_cluster" = length(unique(Zclus_rr[,1])),
                      "n_area" = length(unique(AGBlist_clus$clusdata$AreaNum)),
                      "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                      "niter" = dim(Zclus_rr[,-c(1,2)])[2],
                      "div" = rowMeans(apply(simplify2array(by(Zclus_rr[,-c(1:2)], Zclus_rr$species, as.matrix)),2,rowSums)), ### Sorted by cluster number
                      "carbon_cluster" = rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]),
                      "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)),
                      "area" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$AreaNum)),
                      "biogeo_area" = as.numeric(as.factor(AGBlist_clus$clusdata %>% group_by(AreaNum) %>% summarise(reg = first(bioNum)) %>% pull(reg))))
moddat$endangered <- list("n_cluster" = length(unique(Zclus_endangered[,1])),
                          "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                          "n_area" = length(unique(AGBlist_clus$clusdata$AreaNum)),
                          "niter" = dim(Zclus_endangered[,-c(1,2)])[2],
                          "div" = rowMeans(apply(simplify2array(by(Zclus_endangered[,-c(1:2)], Zclus_endangered$species, as.matrix)),2,rowSums)), ### Sorted by cluster number
                          "carbon_cluster" = rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]),
                          "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)),
                          "area" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$AreaNum)),
                          "biogeo_area" = as.numeric(as.factor(AGBlist_clus$clusdata %>% group_by(AreaNum) %>% summarise(reg = first(bioNum)) %>% pull(reg))))
moddat$rangew <- list("n_cluster" = length(unique(Zclus_rangew[,1])),
                          "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                          "n_area" = length(unique(AGBlist_clus$clusdata$AreaNum)),
                          "niter" = dim(Zclus_rangew[,-c(1,2)])[2],
                          "div" = rowMeans(apply(simplify2array(by(Zclus_rangew[,-c(1:2)], Zclus_rangew$species, as.matrix)),2,rowSums)), ### Sorted by cluster number
                          "carbon_cluster" = rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]),
                          "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)),
                          "area" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$AreaNum)),
                          "biogeo_area" = as.numeric(as.factor(AGBlist_clus$clusdata %>% group_by(AreaNum) %>% summarise(reg = first(bioNum)) %>% pull(reg))))

## Linmod
linmod_div <- rstan::stan(model_code = linmod, data = moddat$richness, chains = 3, iter = 5000)

post_lin <- rstan::extract(linmod_div)

seqlist <- lapply(as.list(sort(unique(moddat$richness$biogeo))), function(x) seq(min(moddat$richness$carbon_cluster[which(moddat$richness$biogeo == x)]),max(moddat$richness$carbon_cluster[which(moddat$richness$biogeo == x)]),0.01))

clusmean <- list()
for(i in 1:length(seqlist)){
  clusmean_reg <- matrix(nrow = length(post_lin$alpha), ncol = length(seqlist[[i]]))
  for(j in 1:length(seqlist[[i]])){
    for(k in 1:length(post_lin$alpha)){
      clusmean_reg[k,j] <- post_lin$alpha_reg[k,i] + post_lin$beta1[k,i] * seqlist[[i]][j]
    }
  }
  clusmean[[i]] <- clusmean_reg
}

globseq <- seq(min(moddat$richness$carbon_cluster),max(moddat$richness$carbon_cluster),1)
clusmean_glob <- matrix(nrow = length(post_lin$alpha), ncol = length(globseq))
b1 <- rowMeans(post_lin$beta_reg)
for(j in 1:length(globseq)){
  for(k in 1:length(post_lin$alpha)){
    clusmean_glob[k,j] <- post_lin$alpha[k] + b1[k] * globseq[j]
  }
}

plot(moddat$richness$div ~ moddat$richness$carbon_cluster, col = as.numeric(as.factor(AGBlist_clus$plotdata %>% group_by(ClusNum) %>% summarise(elev = first(HabitatP)) %>% pull(elev))), pch = 16)
plot(moddat$richness$div ~ moddat$richness$carbon_cluster, col = moddat$richness$biogeo, pch = 16); title("Richness")
plot(moddat$rangew$div ~ moddat$rangew$carbon_cluster, col = moddat$rangew$biogeo, pch = 16); title("Range-weighted")
plot(moddat$endangered$div ~ moddat$endangered$carbon_cluster, col = moddat$endangered$biogeo, pch = 16); title("endangered")
plot(moddat$ranger$div ~ moddat$ranger$carbon_cluster, col = moddat$ranger$biogeo, pch = 16); title("range-restricted")

lines(colMeans(clusmean_glob) ~ globseq)
for(i in 1:length(clusmean)){
  lines(colMeans(clusmean[[i]]) ~ seqlist[[i]], col = i)
}
legend(x = 15, legend = levels(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)), fill = c(1:length(moddat$richness$biogeo)))

text(moddat$richness$carbon_cluster, moddat$richness$div , labels = AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$AreaCode)

### POLYMOD
polymod_div <- rstan::stan(model_code = polymod, data = moddat$richness, chains = 3, iter = 5000)

pairs(polymod_div, pars=c("alpha","beta1","beta2","sigma_cluster","cluster_mean[1]"))

test <- t(matrix(rep(rstan::extract(polymod_div, "alpha")[[1]], length(moddat$richness$carbon_cluster)), ncol = length(moddat$richness$carbon_cluster))) + 
  t(matrix(rep(rstan::extract(polymod_div,"beta1")[[1]], length(moddat$richness$carbon_cluster)), ncol = length(moddat$richness$carbon_cluster))) * moddat$richness$carbon_cluster + 
  t(matrix(rep(rstan::extract(polymod_div,"beta2")[[1]], length(moddat$richness$carbon_cluster)), ncol = length(moddat$richness$carbon_cluster))) * moddat$richness$carbon_cluster ^ 2

test <- t(rstan::extract(polymod_div, "alpha")[[1]][,moddat$richness$biogeo]) + 
  t(rstan::extract(polymod_div, "beta1")[[1]][,moddat$richness$biogeo]) * moddat$richness$carbon_cluster + 
  t(rstan::extract(polymod_div, "beta2")[[1]][,moddat$richness$biogeo]) * (moddat$richness$carbon_cluster ^ 2) 
plot(moddat$richness$div ~ moddat$richness$carbon_cluster, col = moddat$richness$biogeo, pch = 16) ; points(rowMeans(test) ~ moddat$richness$carbon_cluster)

post_poly <- rstan::extract(polymod_div)

seqlist <- lapply(as.list(sort(unique(moddat$richness$biogeo))), function(x) seq(min(moddat$richness$carbon_cluster[which(moddat$richness$biogeo == x)]),max(moddat$richness$carbon_cluster[which(moddat$richness$biogeo == x)]),1))

clusmean <- list()
for(i in 1:length(seqlist)){
  clusmean_reg <- matrix(nrow = length(post_poly$alpha), ncol = length(seqlist[[i]]))
  for(j in 1:length(seqlist[[i]])){
    for(k in 1:length(post_poly$alpha)){
      clusmean_reg[k,j] <- post_poly$alpha_reg[k,i] + post_poly$beta1[k,i] * seqlist[[i]][j] + post_poly$beta2[k,i] * (seqlist[[i]][j] ^ 2)
    }
  }
  clusmean[[i]] <- clusmean_reg
}

globseq <- seq(min(moddat$richness$carbon_cluster),max(moddat$richness$carbon_cluster),1)
clusmean_glob <- matrix(nrow = length(post_poly$alpha), ncol = length(globseq))
b1 <- rowMeans(post_poly$beta1)
b2 <- rowMeans(post_poly$beta2)
for(j in 1:length(globseq)){
  for(k in 1:length(post_poly$alpha)){
    clusmean_glob[k,j] <- post_poly$alpha[k] + b1[k] * globseq[j] + b2[k] * (globseq[j] ^ 2)
  }
}

plot(moddat$richness$div ~ moddat$richness$carbon_cluster, col = moddat$richness$biogeo, pch = 16)
lines(colMeans(clusmean_glob) ~ globseq)
for(i in 1:length(clusmean)){
  lines(colMeans(clusmean[[i]]) ~ seqlist[[i]], col = i)
}
legend(x = 130, legend = levels(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)), fill = c(1:7))

#### BREAKPOINTMOD
bpmod_div <- rstan::stan(model_code = bpmod, data = moddat$richness, chains = 3, iter = 2000)

pairs(bpmod_div, pars=c("alpha","beta1","beta2","sigma_cluster","cluster_mean[1]","breakpoint","test[1]","reg_bp"))

test <- t(matrix(rep(rstan::extract(bpmod_div, "alpha")[[1]], length(moddat$richness$carbon_cluster)), ncol = length(moddat$richness$carbon_cluster))) + 
  t(matrix(rep(rstan::extract(bpmod_div,"beta1")[[1]], length(moddat$richness$carbon_cluster)), ncol = length(moddat$richness$carbon_cluster))) * moddat$richness$carbon_cluster + 
  t(matrix(rep(rstan::extract(bpmod_div,"beta2")[[1]], length(moddat$richness$carbon_cluster)), ncol = length(moddat$richness$carbon_cluster))) * t(rstan::extract(bpmod_div,"test")[[1]])
plot(moddat$richness$div ~ moddat$richness$carbon_cluster, col = moddat$richness$biogeo, pch = 16) ; points(rowMeans(test) ~ moddat$richness$carbon_cluster, col = moddat$richness$biogeo)

test <- t(rstan::extract(bpmod_div, "alpha")[[1]][,moddat$richness$biogeo]) + 
  t(matrix(rep(rstan::extract(bpmod_div,"beta1")[[1]], length(moddat$richness$carbon_cluster)), ncol = length(moddat$richness$carbon_cluster))) * moddat$richness$carbon_cluster + 
  t(matrix(rep(rstan::extract(bpmod_div,"beta2")[[1]], length(moddat$richness$carbon_cluster)), ncol = length(moddat$richness$carbon_cluster))) * t(rstan::extract(bpmod_div,"test")[[1]])
plot(moddat$richness$div ~ moddat$richness$carbon_cluster, col = moddat$richness$biogeo, pch = 16) ; points(rowMeans(test) ~ moddat$richness$carbon_cluster, col = moddat$richness$biogeo)
plot(moddat$richness$div ~ moddat$richness$carbon_cluster, col = moddat$richness$area, pch = moddat$richness$biogeo) ; points(rowMeans(test) ~ moddat$richness$carbon_cluster, col = moddat$richness$biogeo)

plot(moddat$richness$div ~ AGBlist_clus$clusdata$meanAGB[order(AGBlist_clus$clusdata$ClusNum)], col = moddat$richness$biogeo, pch = 16)
points(rowMeans(test) ~ moddat$richness$carbon_cluster, col = moddat$richness$biogeo)

post_bp <- rstan::extract(bpmod_div)

seqlist <- lapply(as.list(sort(unique(moddat$richness$biogeo))), function(x) seq(min(moddat$richness$carbon_cluster[which(moddat$richness$biogeo == x)]),max(moddat$richness$carbon_cluster[which(moddat$richness$biogeo == x)]),1))

clusmean <- list()
for(i in 1:length(seqlist)){
  stepfunc <- matrix(nrow = length(post_bp$alpha), ncol = length(seqlist[[i]]))
  clusmean_reg <- matrix(nrow = length(post_bp$alpha), ncol = length(seqlist[[i]]))
  for(j in 1:length(seqlist[[i]])){
    for(k in 1:length(post_bp$alpha)){
    stepfunc[k,j] <- ifelse(seqlist[[i]][j] > post_bp$breakpoint[k], seqlist[[i]][j] - post_bp$breakpoint[k], 0) + ifelse(seqlist[[i]][j] > post_bp$reg_bp[k,i], seqlist[[i]][j] - post_bp$reg_bp[k,i], 0)
    clusmean_reg[k,j] <- post_bp$alpha_reg[k,i] + post_bp$beta1_reg[k,i] * seqlist[[i]][j] + post_bp$beta2_reg[k,i] * stepfunc[k,j]
    }
  }
  clusmean[[i]] <- clusmean_reg
}

globseq <- seq(min(moddat$richness$carbon_cluster),max(moddat$richness$carbon_cluster),1)
stepfunc <- matrix(nrow = length(post_bp$alpha), ncol = length(globseq))
clusmean_glob <- matrix(nrow = length(post_bp$alpha), ncol = length(globseq))
for(j in 1:length(globseq)){
  for(k in 1:length(post_bp$alpha)){
    stepfunc[k,j] <- ifelse(globseq[j] > post_bp$breakpoint[k], globseq[j] - post_bp$breakpoint[k], 0)
    #clusmean_glob[k,j] <- post_bp$alpha[k] + post_bp$beta1[k] * globseq[j] + post_bp$beta2[k] * stepfunc[k,j]
    clusmean_glob[k,j] <- post_bp$alpha[k] + post_bp$beta1[k] * globseq[j] + post_bp$beta2[k] * stepfunc[k,j]
  }
}

plot(moddat$richness$div ~ moddat$richness$carbon_cluster, col = moddat$richness$biogeo, pch = 16)
lines(colMeans(clusmean_glob) ~ globseq)
for(i in 1:length(clusmean)){
  lines(colMeans(clusmean[[i]]) ~ seqlist[[i]], col = i)
}
legend(x = 130, legend = levels(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)), fill = c(1:7))

## OTHER
plot(moddat$richness$div ~ moddat$richness$carbon_cluster, col = RColorBrewer::brewer.pal(9, "Blues")[round(scales::rescale(AGBlist_clus$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev), c(1,9)))])


plot(moddat$richness$div ~ AGBlist_clus$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))
plot(moddat$richness$carbon_cluster ~ AGBlist_clus$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))





## Economics
moddat <- list()
value <- read.csv("area_dep.csv")

moddat$richness <- list("n_cluster" = length(unique(Zclus[,1])),
                        "n_area" = length(unique(AGBlist_clus$clusdata$AreaNum)),
                        "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                        "niter" = dim(Zclus[,-c(1,2)])[2],
                        "div" = rowMeans(apply(simplify2array(by(Zclus[,-c(1:2)], Zclus$species, as.matrix)),2,rowSums)), ### Sorted by cluster number
                        "carbon_cluster" = left_join(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),], value)$value / rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]),
                        "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)))

plot(rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]), 
     left_join(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),], value)$value, 
     col = moddat$richness$biogeo, pch = 16)

plot(moddat$richness$div ~ log(moddat$richness$carbon_cluster), col = moddat$richness$biogeo, pch = 16, xaxt = 'n')
axis(1, at= log(c(0.5,1,2,5,10,20)), labels=c(0.5,1,2,5,10,20))



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








