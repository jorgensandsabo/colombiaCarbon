#################
################# Carbon and biodiversity relationships
#################
if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}

setwd(dir.path)

library(dplyr)
library(brms)
library(ggplot2)

source("Code\\diversity_code\\get_posterior_z_v6.R")
source("Code\\colombiaCarbon\\carbondiversity_analysis\\carbondiversity_functions.R")

##########################
## Read and filter data ##
##########################
# Read AGB data
AGBlist_clus <- readRDS("Output\\CarbDiv\\AGBlist_clus.RDS")

# Read bird data
bird_data <- readRDS("Data\\diversity\\bird_stan_data6_package.RDS")
birds <- readRDS("Data\\diversity\\birds.RDS")
draws <- posterior::as_draws_df(readRDS("Data\\diversity\\draws_thinned_500.RDS"))

# create z_info object for computing posterior Z (see get_posterior_z.R)
z_info <- data.frame(bird_data$data[8:41])
z_info$point <- birds$point
z_info$species <- birds$species
#z_info$cl_q_real <- cluster_q(z_info, z_info$Q)[z_info$id_spCl]

# Remove points
#AGBlist_clus$plotdata <- dplyr::filter(AGBlist_clus$plotdata, biogeo!= "West Andes")
#AGBlist_clus$plotdata <- dplyr::filter(AGBlist_clus$plotdata, HabitatP != "Paramo" & AreaCode != "TU")

# Remove point without bird data
AGBlist_clus$plotdata <- AGBlist_clus$plotdata[which(AGBlist_clus$plotdata$SiteCode %in% z_info$point),]
AGBlist_clus$AGBdraws <- AGBlist_clus$AGBdraws[which(AGBlist_clus$AGBdraws$ClusNum %in% unique(AGBlist_clus$plotdata$ClusNum)),]
AGBlist_clus$clusdata <- AGBlist_clus$clusdata[which(AGBlist_clus$clusdata$ClusNum %in% unique(AGBlist_clus$plotdata$ClusNum)),]

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

# Draw bird data (for points in carbon plotdata only)
Zdraws <- get_specZ_draws_expand(niter = 10, draws, z_info, pointid = AGBlist_clus$plotdata[,c("SiteCode","SiteNum")], spatial_effect = "include")
psis <- get_pasture_forest_psi(niter = 10, draws, z_info)
psis <- lapply(psis, function(x) dplyr::left_join(x, data.frame(z_info %>% group_by(id_sp) %>% summarise(species = unique(species)))) %>% relocate("species",.before = "id_sp"))

# Remove points without carbon data
Zdraws <- dplyr::filter(Zdraws, SiteNum %in% AGBlist_clus$plotdata$SiteNum)

# Filter out full migrants
Zdraws <- dplyr::filter(Zdraws, species %in% dplyr::filter(birdlife_extents, (!migratory_status  %in% c("full migrant","unknown")))$X)
psis <- lapply(psis, function(x) dplyr::filter(x, species %in% dplyr::filter(birdlife_extents, (!migratory_status  %in% c("full migrant","unknown")))$X))

# Point to cluster occurrence
Zclus <- aggregate(formula = .~ ClusNum + species, data = left_join(Zdraws, AGBlist_clus$plotdata[,c("SiteNum","ClusNum")])[,-1], FUN = max)

# Weight by pasture occurrence
psis_weight <- lapply(psis, function(x) do.call("rbind", lapply(Zclus$species, function(y) x[which(x$species == y),])))
pasture_weight <- log(psis_weight[[1]][,-c(1:2)] / psis_weight[[2]][,-c(1:2)])
pasture_weight <- max(pasture_weight) - pasture_weight
pasture_weight <- pasture_weight / max(pasture_weight)

# Range weighted species richness
extent <- dplyr::filter(birdlife_extents, !migratory_status  %in% c("full migrant","unknown"))[,c("X","extent_breeding_resident")]
extent <- dplyr::filter(extent, X %in% Zclus$species)
extent$inv <- max(log(extent$extent_breeding_resident)) - log(extent$extent_breeding_resident)
extent$weight <- extent$inv / max(extent$inv)
extent_weights <- left_join(Zclus[,c(1:2)], extent, by = c("species" = "X"))

# Weight by pasture occurrence and range
Zclus_weighted <- cbind(Zclus[,c(1:2)], Zclus[,-c(1:2)] * extent_weights$weight * pasture_weight)

# Species of conservation concern
Zclus_concern <- dplyr::filter(Zclus, species %in% c(dplyr::filter(birdlife_extents, extent_breeding_resident < 50000)$X, filter(redlist_status, V1 %in% c("Critically Endangered","Endangered","Vulnerable"))$X))

#############################
## Carbon-diversity models ##
#############################
# Model data
moddat <- list()
moddat$richness <- list("n_cluster" = length(unique(Zclus[,1])),
                        "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                        "niter" = dim(Zclus[,-c(1,2)])[2],
                        "div" = rowMeans(apply(simplify2array(by(Zclus[,-c(1:2)], Zclus$species, as.matrix)),2,rowSums)), ### Sorted by cluster number
                        "predictor_cluster" = scale(rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]))[,1],
                        "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)))
moddat$concern <- list("n_cluster" = length(unique(Zclus_concern[,1])),
                       "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                       "niter" = dim(Zclus_concern[,-c(1,2)])[2],
                       "div" = rowMeans(apply(simplify2array(by(Zclus_concern[,-c(1:2)], Zclus_concern$species, as.matrix)),2,rowSums)), ### Sorted by cluster number
                       "predictor_cluster" = scale(rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]))[,1],
                       "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)))
moddat$weighted <- list("n_cluster" = length(unique(Zclus_weighted[,1])),
                        "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                        "niter" = dim(Zclus_weighted[,-c(1,2)])[2],
                        "div" = (rowMeans(apply(simplify2array(by(Zclus_weighted[,-c(1:2)], Zclus_weighted$species, as.matrix)),2,rowSums)) / max(rowMeans(apply(simplify2array(by(Zclus_weighted[,-c(1:2)], Zclus_weighted$species, as.matrix)),2,rowSums)))) * 100, ### Sorted by cluster number
                        "predictor_cluster" = scale(rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]))[,1],
                        "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)))
moddat$elevation <- list("n_cluster" = length(unique(Zclus[,1])),
                         "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                         "niter" = dim(Zclus[,-c(1,2)])[2],
                         "predictor_cluster" = scale(as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$ALOSelev)))[,1],
                         "div" = (rowMeans(apply(simplify2array(by(Zclus_weighted[,-c(1:2)], Zclus_weighted$species, as.matrix)),2,rowSums)) / max(rowMeans(apply(simplify2array(by(Zclus_weighted[,-c(1:2)], Zclus_weighted$species, as.matrix)),2,rowSums)))) * 100, ### Sorted by cluster number
                         "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)))

models <- list(richness = list(), concern = list(), weighted = list(), elevation = list())

models$richness$poly_intercept <- rstan::stan(file = "STAN\\carbdiv_poly_intercept.stan", data = moddat$richness, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$richness$poly_regintercept <- rstan::stan(file = "STAN\\carbdiv_poly_regintercept.stan", data = moddat$richness, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$richness$poly_regslope <- rstan::stan(file = "STAN\\carbdiv_poly_regslope.stan", data = moddat$richness, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15))
models$richness$lin_intercept <- rstan::stan(file = "STAN\\carbdiv_lin_intercept.stan", data = moddat$richness, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$richness$lin_regintercept <- rstan::stan(file = "STAN\\carbdiv_lin_regintercept.stan", data = moddat$richness, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$richness$lin_regslope <- rstan::stan(file = "STAN\\carbdiv_lin_regslope.stan", data = moddat$richness, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15))

models$concern$poly_intercept <- rstan::stan(file = "STAN\\carbdiv_poly_intercept.stan", data = moddat$concern, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$concern$poly_regintercept <- rstan::stan(file = "STAN\\carbdiv_poly_regintercept.stan", data = moddat$concern, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$concern$poly_regslope <- rstan::stan(file = "STAN\\carbdiv_poly_regslope.stan", data = moddat$concern, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15))
models$concern$lin_intercept <- rstan::stan(file = "STAN\\carbdiv_lin_intercept.stan", data = moddat$concern, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$concern$lin_regintercept <- rstan::stan(file = "STAN\\carbdiv_lin_regintercept.stan", data = moddat$concern, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$concern$lin_regslope <- rstan::stan(file = "STAN\\carbdiv_lin_regslope.stan", data = moddat$concern, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15))

models$weighted$poly_intercept <- rstan::stan(file = "STAN\\carbdiv_poly_intercept.stan", data = moddat$weighted, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$weighted$poly_regintercept <- rstan::stan(file = "STAN\\carbdiv_poly_regintercept.stan", data = moddat$weighted, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$weighted$poly_regslope <- rstan::stan(file = "STAN\\carbdiv_poly_regslope.stan", data = moddat$weighted, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15))
models$weighted$lin_intercept <- rstan::stan(file = "STAN\\carbdiv_lin_intercept.stan", data = moddat$weighted, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$weighted$lin_regintercept <- rstan::stan(file = "STAN\\carbdiv_lin_regintercept.stan", data = moddat$weighted, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$weighted$lin_regslope <- rstan::stan(file = "STAN\\carbdiv_lin_regslope.stan", data = moddat$weighted, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15))

models$elevation$poly_intercept <- rstan::stan(file = "STAN\\carbdiv_poly_intercept.stan", data = moddat$elevation, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$elevation$poly_regintercept <- rstan::stan(file = "STAN\\carbdiv_poly_regintercept.stan", data = moddat$elevation, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$elevation$poly_regslope <- rstan::stan(file = "STAN\\carbdiv_poly_regslope.stan", data = moddat$elevation, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15))
models$elevation$lin_intercept <- rstan::stan(file = "STAN\\carbdiv_lin_intercept.stan", data = moddat$elevation, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$elevation$lin_regintercept <- rstan::stan(file = "STAN\\carbdiv_lin_regintercept.stan", data = moddat$elevation, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
models$elevation$lin_regslope <- rstan::stan(file = "STAN\\carbdiv_lin_regslope.stan", data = moddat$elevation, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15))

#pairs(models$richness$poly_regintercept2, pars = c("alpha","beta1","beta2","sigma_cluster","cluster_mean[1]"))
#pairs(models$weighted$poly_regintercept, pars = c("alpha","beta1","beta2","sigma_cluster","cluster_mean[1]"))

## Elevation + carbon model
carb_elevation_data <- list("n_cluster" = length(unique(Zclus[,1])),
                            "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                            "niter" = dim(Zclus[,-c(1,2)])[2],
                            "elevation_cluster" = scale(as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$ALOSelev)))[,1],
                            "div" = (rowMeans(apply(simplify2array(by(Zclus_weighted[,-c(1:2)], Zclus_weighted$species, as.matrix)),2,rowSums)) / max(rowMeans(apply(simplify2array(by(Zclus_weighted[,-c(1:2)], Zclus_weighted$species, as.matrix)),2,rowSums)))) * 100, ### Sorted by cluster number
                            "carbon_cluster" = scale(rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]))[,1],
                            "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)))

##################
## Model output ##
##################
## Model comparison table
loos <- lapply(models, function(x) lapply(x, rstan::loo, moment_match = T, k_threshold = 0.5))
loo_comp <- lapply(loos, function(x) loo::loo_compare(x))
print(loo_comp, simplify=F)
write.csv(do.call(rbind, loo_comp), "Output\\CarbDiv\\loo_table_divmod.csv")

# Compare elevation and carbon models
carbon_elevation_model <- rstan::stan(file = "STAN\\carb_elevation_div_model.stan", data = carb_elevation_data, chains = 4, iter = 2000, cores = 4, control = list(max_treedepth = 15), init = 0)
finalloos <- loo::loo_compare(lapply(list(carbon_elevation_model = carbon_elevation_model, 
                                          carbonmodel = models$weighted$poly_regslope, 
                                          elevationmodel = models$elevation$poly_regintercept), 
                                     loo, moment_match = T, threshold = 0.5))
write.csv(finalloos, "Output\\CarbDiv\\loo_table_final.csv")

## Prepare plotting data
postsamp <- lapply(models, function(x) lapply(x, rstan::extract))
seqlist <- lapply(moddat, function(y) lapply(as.list(sort(unique(y$biogeo))), function(x) seq(min(y$predictor_cluster[which(y$biogeo == x)]),max(y$predictor_cluster[which(y$biogeo == x)]),0.01)))

clusmeans <- list(richness = list(), concern = list(), weighted = list(), elevation = list())
for(m in 1:length(clusmeans)){ # Number of indices
  clusmean_dat <- list()
  seq_temp <- seqlist[[m]]
  dat_temp <- moddat[[m]]
  for(l in 1:length(postsamp[[m]])){ # Number of models
    clusmean_mod <- list()
    post_temp <- postsamp[[m]][[l]]
    if(!is.null(post_temp$alpha_area)){
      test <- t(apply(post_temp$alpha_area, 1, function(x) aggregate(x, list(dat_temp$biogeo_area), mean)[,2]))
    }
    for(i in 1:length(seq_temp)){ # Number of regions
      clusmean_reg <- matrix(nrow = length(post_temp$lp__), ncol = length(seq_temp[[i]]))
      for(j in 1:length(seq_temp[[i]])){ # Number of clusters
        for(k in 1:length(post_temp$lp__)){ # Number of iterations
          clusmean_reg[k,j] <- 
            ifelse(is.null(post_temp$alpha_biogeo), post_temp$alpha[k], post_temp$alpha_biogeo[k,i]) + 
            ifelse(is.null(post_temp$beta1), 0,
                   ifelse(is.null(post_temp$beta1_biogeo), post_temp$beta1[k], post_temp$beta1_biogeo[k,i]) * seq_temp[[i]][j]) +
            ifelse(is.null(post_temp$beta2), 0, 
                   ifelse(is.null(post_temp$beta2_biogeo), post_temp$beta2[k], post_temp$beta2_biogeo[k,i])) * (seq_temp[[i]][j] ^ 2)+
            ifelse(is.null(post_temp$betaelev1), 0,
                   ifelse(is.null(post_temp$betaelev1_biogeo), post_temp$betaelev1[k], post_temp$betaelev1_biogeo[k,i]) * seq_temp[[i]][j]) + 
            ifelse(is.null(post_temp$betaelev2), 0,
                   ifelse(is.null(post_temp$betaelev2_biogeo), post_temp$betaelev2[k], post_temp$betaelev2_biogeo[k,i]) * seq_temp[[i]][j])
        }
      }
      clusmean_mod[[i]] <- clusmean_reg
    }
    clusmean_dat[[l]] <- clusmean_mod
  }
  clusmeans[[m]] <- clusmean_dat
}

# Plot all
par(mfrow=c(length(seqlist),6))
for(i in 1:length(seqlist)){
  toplot <- clusmeans[[i]]
  dattoplot <- moddat[[i]]
  seqtoplot <- seqlist[[i]]
  for(j in 1:length(toplot)){
    plot(dattoplot$div ~ dattoplot$predictor_cluster, col = dattoplot$biogeo, pch = 16)
    for(k in 1:length(toplot[[j]])){
      lines(colMeans(toplot[[j]][[k]]) ~ seqtoplot[[k]], col = k)
    }
  }
}

## Plot final models
tiff(file = "Output\\CarbDiv\\Divregression_plots.tiff", width = 7000, height = 5000, res = 900)
  layout(matrix(c(1,2,5,3,4,5), 2, 3, byrow = TRUE), widths = c(0.4,0.4,0.2))
  par(mar = c(4, 4, 0.5, 0.5))
  # Diversity plots
  modplot <- list(data = c(3,2,4),
                  models = c(3,3,2),
                  xlabs = c("Landscape biomass","Landscape biomass","Elevation"),
                  ylabs = c("Relative weighted species richness","Species of conservation concern","Relative weighted species richness"))
  for(i in 1:length(modplot[[1]])){
    toplot <- clusmeans[[modplot$data[i]]]
    dattoplot <- moddat[[modplot$data[i]]]
    seqtoplot <- seqlist[[modplot$data[i]]]
    for(j in modplot$models[i]){
      plot(dattoplot$div ~ dattoplot$predictor_cluster, col = dattoplot$biogeo, pch = 16, 
           xlab = modplot$xlabs[i], ylab = modplot$ylabs[i])
      for(k in 1:length(toplot[[j]])){
        lines(colMeans(toplot[[j]][[k]]) ~ seqtoplot[[k]], col = k)
      }
    }
  }
  # Carbon/elevation plot
  plot(carb_elevation_data$carbon_cluster ~ carb_elevation_data$elevation_cluster, col = carb_elevation_data$biogeo, pch = 16,
       xlab = "Elevation", ylab = "Landscape biomass")
  legend("topright", bty = "n",
         paste("Pearson's r =", round(cor(carb_elevation_data$carbon_cluster,carb_elevation_data$elevation_cluster), 2)))
  # Legend
  par(mar = c(0,0,0,0))
  plot(0, xaxt = 'n', yaxt = 'n', pch = '', ylab = '', xlab = '', bty = 'n')
  legend("left", legend = unique(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo), 
         fill = unique(unique(as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)))), box.col = 'white')
dev.off()

## 
test <- lapply(list("Richness~carbon" = models$weighted$poly_regslope,
                    "Elevation~carbon" = models$elevation$poly_regintercept,
                    "Richness~elevation+carbon" = carbon_elevation_model,
                    "Richness~priority" = models$concern$poly_regslope),
                    function(x) rstan::extract(x))

test2 <- lapply(test, function(x) lapply(x[which(names(x) %in% c("alpha0","sigma_cluster","beta1","beta2","betaelev1","betaelev2"))], function(y) if(length(dim(y)) > 1){cbind(colMeans(y),apply(y, 2, function(y) y[order(y)][length(y) * 0.025]),apply(y, 2, function(y) y[order(y)][length(y) * 0.975]))}
                                               else{cbind(mean(y), y[order(y)][length(y)*0.025],  y[order(y)][length(y)*0.975])}))
test3 <- lapply(test2, function(x) data.frame(parameter = rep(names(x), times = unlist(lapply(x, nrow))), do.call(rbind, x)))

test4 <- Reduce(function(x,y) merge(x = x, y = y, by = "parameter", all = TRUE), test3)
names(test4)[-1] <- rep(names(test2), each = 3)
write.csv(test4, "Output\\CarbDiv\\finalmodoutput.csv")

### TO DO:
# 4: Turnover
# 5: Economics
# 6: Check Moran's I

### ELEVATION
moddat$richness2 <- list("n_cluster" = length(unique(Zclus[,1])),
                        "n_area" = length(unique(AGBlist_clus$clusdata$AreaNum)),
                        "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                        "niter" = dim(Zclus[,-c(1,2)])[2],
                        "carbon_region" = aggregate(x = rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]), 
                                                    by = list(as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo))),
                                                    mean)$x,
                        "div" = rowMeans(apply(simplify2array(by(Zclus[,-c(1:2)], Zclus$species, as.matrix)),2,rowSums)), ### Sorted by cluster number
                        "elevation" = scale(AGBlist_clus$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev), center = F)[,1],
                        "carbon_cluster" = rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]),
                        "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)),
                        "area" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$AreaNum)),
                        "biogeo_area" = as.numeric(as.factor(AGBlist_clus$clusdata %>% group_by(AreaNum) %>% summarise(reg = first(bioNum)) %>% pull(reg))))

polymod_elev <- rstan::stan(model_code = polymod_elev, data = moddat$richness2, chains = 3, iter = 5000)
polymod_elev2 <- rstan::stan(model_code = polymod_elev2, data = moddat$richness2, chains = 3, iter = 5000)
polymod_elev3 <- rstan::stan(model_code = polymod_elev3, data = moddat$richness2, chains = 3, iter = 5000)
polymod_elev4 <- rstan::stan(model_code = polymod_elev4, data = moddat$richness2, chains = 3, iter = 5000)

plot(moddat$richness$div ~ moddat$richness2$div)

## Economics
## 1 t carbon = 3.67 t CO2, 1 t biomass = 0.5 t C?
value <- read.csv("area_dep.csv")

moddat$richness <- list("n_cluster" = length(unique(Zclus[,1])),
                        "n_area" = length(unique(AGBlist_clus$clusdata$AreaNum)),
                        "n_biogeo" = length(unique(AGBlist_clus$clusdata$biogeo)),
                        "niter" = dim(Zclus[,-c(1,2)])[2],
                        "div" = (rowMeans(apply(simplify2array(by(Zclus_weighted[,-c(1:2)], Zclus_weighted$species, as.matrix)),2,rowSums)) / max(rowMeans(apply(simplify2array(by(Zclus_weighted[,-c(1:2)], Zclus_weighted$species, as.matrix)),2,rowSums)))) * 100, ### Sorted by cluster number
                        "carbon_cluster" = left_join(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),], value)$value / 
                          (rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]) / 0.47) * 3.67,
                        "elevation_cluster" = AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$ALOSelev,
                        "biogeo" = as.numeric(as.factor(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),]$biogeo)))

plot(rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]), 
     left_join(AGBlist_clus$clusdata[order(AGBlist_clus$clusdata$ClusNum),], value)$value, 
     col = moddat$richness$biogeo, pch = 16)

par(mfrow=c(1,2))
plot(moddat$richness$div ~ log(moddat$richness$carbon_cluster), col = moddat$richness$biogeo, pch = 16, xaxt = 'n',
     xlab = "Break-even CO2 price (log scale)", ylab = "Species richness")
axis(1, at= log(c(0.5,1,2,5,10,20)), labels=c(0.5,1,2,5,10,20))
abline(v = log(3.2))
plot(moddat$richness$elevation_cluster ~ log(moddat$richness$carbon_cluster), col = moddat$richness$biogeo, pch = 16, xaxt = 'n',
     xlab = "Break-even CO2 price (log scale)", ylab = "Elevation")
axis(1, at= log(c(0.5,1,2,5,10,20)), labels=c(0.5,1,2,5,10,20))
abline(v = log(3.2))

####################
## Beta diversity ##
####################
Zspecsite <- reshape2::dcast(Zclus, ClusNum ~ species, value.var = "2")

Zbetadiv <- betapart::beta.pair(Zspecsite[,-1])$beta.sim
Zbetadiv <- data.frame(SiteNum = Zspecsite$ClusNum, as.matrix(Zbetadiv))

predData = data.frame(SiteNum = AGBlist_clus$clusdata$ClusNum[order(AGBlist_clus$clusdata$ClusNum)], 
                      long = AGBlist_clus$clusdata$long[order(AGBlist_clus$clusdata$ClusNum)], 
                      lat = AGBlist_clus$clusdata$lat[order(AGBlist_clus$clusdata$ClusNum)],
                      AGB = rowMeans(AGBlist_clus$AGBdraws[order(AGBlist_clus$AGBdraws$ClusNum),-1]), 
                      elev = AGBlist_clus$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))

sitepairs <- gdm::formatsitepair(bioData = Zbetadiv[order(Zbetadiv$SiteNum),], bioFormat = 3, dist = "bray", siteColumn = "SiteNum", XColumn = "long", YColumn = "lat",
                                 predData = predData[c("SiteNum","elev","AGB","long","lat")])

##
gdmmods <- list()
gdmmods$geo <- gdm::gdm(gdm::formatsitepair(bioData = Zbetadiv[order(Zbetadiv$SiteNum),], bioFormat = 3, dist = "bray", siteColumn = "SiteNum", XColumn = "long", YColumn = "lat",
                                        predData = predData[c("SiteNum","long","lat")]), geo = T, splines = c(3))
gdmmods$AGB <- gdm::gdm(gdm::formatsitepair(bioData = Zbetadiv[order(Zbetadiv$SiteNum),], bioFormat = 3, dist = "bray", siteColumn = "SiteNum", XColumn = "long", YColumn = "lat",
                                       predData = predData[c("SiteNum","AGB","long","lat")]), geo = T, splines = c(3,3))
gdmmods$elev <- gdm::gdm(gdm::formatsitepair(bioData = Zbetadiv[order(Zbetadiv$SiteNum),], bioFormat = 3, dist = "bray", siteColumn = "SiteNum", XColumn = "long", YColumn = "lat",
                                       predData = predData[c("SiteNum","elev","long","lat")]), geo = T, splines = c(3,3))
gdmmods$full <- gdm::gdm(gdm::formatsitepair(bioData = Zbetadiv[order(Zbetadiv$SiteNum),], bioFormat = 3, dist = "bray", siteColumn = "SiteNum", XColumn = "long", YColumn = "lat",
                                       predData = predData[c("SiteNum","elev","AGB","long","lat")]), geo = T, splines = c(3,3,3))

# Plot
test(x = gdmmods, plot.layout = c(3,2), 
              point.cols = c("yellow","blue","red","green"), alpha = 0.1, ppch = 20, pcex = 2,
              line.cols = rep("black", 4), plot.linewidth = 2)

gdmtab <- rbind(data.frame(null = gdmmods$AGB$nulldeviance, do.call(data.frame, lapply(gdmmods, function(x) x$gdmdeviance))),
                data.frame(null = NA, do.call(data.frame, lapply(gdmmods, function(x) x$explained))))
gdmtab <- as.data.frame(t(gdmtab), row.names = names(gdmtab))
names(gdmtab) <- c("deviance","%explained")
write.csv(gdmtab, "Output\\CarbDiv\\gdmtab.csv")



