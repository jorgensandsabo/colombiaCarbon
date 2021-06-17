if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}
setwd(dir.path)
library(dplyr)

##############################
## Upscale to cluster level ##
##############################
AGBlist <- readRDS("Output\\CarbDiv\\With_choco\\AGBlist.RDS")
spatial <- read.csv("Data\\vegetation\\Spatialdata.csv")
AGBlist$plotdata <- dplyr::left_join(AGBlist$plotdata, spatial[,c("SiteCode","ALOSelev","TotPrec","TempVar","PrecVar")])

AGBlist_forest <- lapply(AGBlist, function(x) x[which(x$SiteNum %in% AGBlist$plotdata[AGBlist$plotdata$HabitatP == "Forest",]$SiteNum),])
AGBlist_paramo <- lapply(AGBlist, function(x) x[which(x$SiteNum %in% AGBlist$plotdata[AGBlist$plotdata$HabitatP == "Paramo",]$SiteNum),])

AGBlist_forest$plotdata <- dplyr::filter(AGBlist_forest$plotdata, Cluster != "JU5")

AGBlist_forest$plotdata <- AGBlist_forest$plotdata[-which(AGBlist_forest$plotdata$Cluster == "JU5"),]
AGBlist_forest$AGBdraws <- AGBlist_forest$AGBdraws[-which(AGBlist_forest$plotdata$Cluster == "JU5"),]
forestpredictors <- cbind(scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))[,1],
                          scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(prec = mean(TotPrec)) %>% pull(prec))[,1],
                          scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(tempvar = mean(TempVar)) %>% pull(tempvar))[,1],
                          scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(precvar = mean(PrecVar)) %>% pull(precvar))[,1])
forestpredictors <- matrix(scale(AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))[,1])
paramopredictors <- matrix(scale(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(elev = mean(ALOSelev)) %>% pull(elev))[,1])

forestdata <- list(n_point = nrow(AGBlist_forest$plotdata),                                                                          # Number of points
                   n_cluster = length(unique(AGBlist_forest$plotdata$ClusNum)),                                                      # Number of clusters
                   n_area = length(unique(AGBlist_forest$plotdata$AreaNum)),                                                         # Number of areas
                   n_pred = ncol(forestpredictors),                                                                                  # Number of predictors
                   lcarbon_point = log(as.numeric(rowMeans(AGBlist_forest$AGBdraws[order(AGBlist_forest$AGBdraws$SiteNum),-1]))),    # Carbon at each point
                   cluster = as.numeric(as.factor(AGBlist_forest$plotdata[order(AGBlist_forest$plotdata$SiteNum),]$ClusNum)),        # Cluster at each point
                   plotsize = AGBlist_forest$plotdata[order(AGBlist_forest$plotdata$SiteNum),] 
                                %>% mutate(plotsize = case_when(Size < 100 ~ 1, Size >= 100 ~ 0)) %>% pull(plotsize),                  # Plotsize for each point
                   predictor = forestpredictors,                                                                                     # Predictors at each cluster
                   area_cluster = AGBlist_forest$plotdata %>% group_by(ClusNum) %>% summarise(area = first(AreaNum)) %>% pull(area)) # Area at each cluster

paramodata <- list(n_point = nrow(AGBlist_paramo$plotdata),                                                                          # Number of points
                   n_cluster = length(unique(AGBlist_paramo$plotdata$ClusNum)),                                                      # Number of clusters
                   n_area = length(unique(AGBlist_paramo$plotdata$AreaNum)),                                                         # Number of areas
                   n_pred = ncol(paramopredictors),                                                                                  # Number of predictors
                   lcarbon_point = log(as.numeric(rowMeans(AGBlist_paramo$AGBdraws[order(AGBlist_paramo$AGBdraws$SiteNum),-1]))),    # Carbon at each point
                   predictor = paramopredictors,                                                                                     # Predictors at each cluster
                   area_cluster = as.numeric(as.factor(AGBlist_paramo$plotdata %>% group_by(ClusNum) %>% summarise(area = first(AreaNum)) %>% pull(area))), # Area at each cluster
                   cluster = as.numeric(as.factor(AGBlist_paramo$plotdata[order(AGBlist_paramo$plotdata$SiteNum),]$ClusNum)))        # Cluster at each point

# Intercept only
mod1_forest <- rstan::stan(file = "STAN\\carbdiv_clus_forest_1.stan", data = forestdata, chains = 4, iter = 5000, thin = 5, cores = 4, control = list(adapt_delta = 0.99))
# Elevation covariate
mod2_forest <- rstan::stan(file = "STAN\\carbdiv_clus_forest_2.stan", data = forestdata, chains = 4, iter = 5000, thin = 5, cores = 4, control = list(adapt_delta = 0.99))
# Elevation + elevation ^ 2 
#mod3_forest <- rstan::stan(file = "STAN\\carbdiv_clus_forest_3.stan", data = forestdata, chains = 4, iter = 5000, thin = 5, cores = 4, control = list(adapt_delta = 0.99))
# Elevation + heteroscedasticity (sigma = fraction of the mean)
mod4_forest <- rstan::stan(file = "STAN\\carbdiv_clus_forest_4.stan", data = forestdata, chains = 4, iter = 5000, thin = 5, cores = 4, control = list(adapt_delta = 0.99))
# Elevation + heteroscedasticity (sigma = sigma0 * mean ^ term)
mod5_forest <- rstan::stan(file = "STAN\\carbdiv_clus_forest_5.stan", data = forestdata, chains = 4, iter = 5000, thin = 5, cores = 4, control = list(adapt_delta = 0.99))
# Plot size sigma
#mod6_forest <- rstan::stan(file = "STAN\\carbdiv_clus_forest_6.stan", data = forestdata, chains = 3, iter = 2000, cores = 3, control = list(adapt_delta = 0.99))
# Plot size sigma + heteroscedasticity (sigma = fraction of the mean)
#mod7_forest <- rstan::stan(file = "STAN\\carbdiv_clus_forest_7.stan", data = forestdata, chains = 3, iter = 2000, cores = 3, control = list(adapt_delta = 0.99))

# Pairs plots
pairs(mod1_forest, pars = c("alpha","lsigma_point","sigma_cluster","cluster_mean","cluster_lmean[1]","carbon_cluster[1]"))
pairs(mod2_forest, pars = c("alpha","lsigma_point","sigma_cluster","cluster_mean[1]","cluster_lmean[1]","carbon_cluster[1]","beta"))
pairs(mod3_forest, pars = c("alpha","lsigma_point","sigma_cluster","cluster_mean[1]","cluster_lmean[1]","carbon_cluster[1]","beta","beta2"))
pairs(mod4_forest, pars = c("alpha","lsigma_point","sigma_cluster[1]","meanfraction","cluster_mean[1]","cluster_lmean[1]","carbon_cluster[1]","beta"))
#pairs(mod5_forest, pars = c("alpha","lsigma_point","sigma_cluster[1]","sigma_cluster0","meanfraction","cluster_mean[1]","cluster_lmean[1]","carbon_cluster[1]","beta"))
#pairs(mod6_forest, pars = c("alpha","lsigma_point[1]","lsigma_small","lsigma_offset","sigma_cluster","cluster_mean[1]","cluster_lmean[1]","carbon_cluster[1]","beta"))
#pairs(mod7_forest, pars = c("alpha","lsigma_point[103]","lsigma_small","meanfraction","lsigma_offset","sigma_cluster[103]","cluster_mean[103]","cluster_lmean[103]","carbon_cluster[103]","beta"))

clusmod_forest <- rstan::stan(file = "STAN\\carbdiv_clus_forest_intercept.stan", data = forestdata, init = function(){list(alpha = rnorm(1, 200, 20))}, 
                              chains = 4, cores = 4, iter = 5000, thin = 1, control = list(adapt_delta = 0.99, max_treedepth = 15))
clusmod_paramo <- rstan::stan(file = "STAN\\carbdiv_clus_paramo.stan", data = paramodata, init = function(){list(alpha = rnorm(1, 30, 5))}, 
                              chains = 4, cores = 4, iter = 2000, thin = 1, control = list(adapt_delta = 0.99, max_treedepth = 15))

## LOO comparison
models <- list(mod1_forest, mod2_forest, mod4_forest)#, mod3_forest, mod5_forest, mod6_forest, mod7_forest)
loo::loo_compare(lapply(models, loo::loo))

## Plot differences in cluster carbon estimates
plot(colMeans(rstan::extract(mod4_forest,"carbon_cluster")[[1]]), 
     colMeans(rstan::extract(mod5_forest,"carbon_cluster")[[1]]),
     xlab = "fixed variance", ylab = "heteroscedastic model", main = "model estimate difference")
abline(0,1)

# Errors plot
postvalues <- rstan::extract(mod4_forest)
par(mfrow = c(2,3))
plot(rowMeans(forestdata$lcarbon_point - t(postvalues$cluster_lmean)) ~ rowMeans(t(postvalues$cluster_lmean)),
     xlab = "Estimated expected point carbon", ylab = "Observed minus expected point carbon") ; abline(h = 0)
plot(aggregate(exp(forestdata$lcarbon_point), list(forestdata$cluster), mean)[,2] - colMeans(postvalues$cluster_mean) ~ colMeans(postvalues$cluster_mean), 
     xlab = "Estimated expected cluster carbon", ylab = "Observed minus expected cluster carbon") ; abline(h = 0)
plot(colMeans(postvalues$carbon_cluster) - colMeans(postvalues$cluster_mean) ~ colMeans(postvalues$cluster_mean), 
     xlab = "Estimated expected cluster carbon", ylab = "Estimated minus expected cluster carbon") ; abline(h = 0)
plot(aggregate(exp(forestdata$lcarbon_point), list(forestdata$cluster), mean)[,2] ~ colMeans(postvalues$carbon_cluster), 
     xlab = "Estimated cluster carbon", ylab = "Observed cluster carbon") ; abline(0,1)
plot(aggregate(exp(forestdata$lcarbon_point), list(forestdata$cluster), mean)[,2] ~ forestdata$predictor[,1],
     xlab = "Elevation", ylab = "Observed average cluster value")
        points(mean(postvalues$alpha) + mean(postvalues$beta) * forestdata$predictor[,1] ~ forestdata$predictor[,1], col = "red")
plot(colMeans(postvalues$carbon_cluster) ~ forestdata$predictor[,1],
     xlab = "Elevation", ylab = "Estimated cluster value")
        points(mean(postvalues$alpha) + mean(postvalues$beta) * forestdata$predictor[,1] ~ forestdata$predictor[,1], col = "red")

## Histogram
postvalues <- rstan::extract(mod5_forest)
hist(postvalues$meanfraction, main = "fraction of the mean parameter")

## ppc_dens
test <- rstan::extract(mod4_forest, pars = "cluster_lmean")[[1]]
test2 <- rstan::extract(mod4_forest, pars = "lsigma_point")[[1]]
test3 <- sapply(1:nrow(test), function(x) rnorm(ncol(test), test[x,], test2[x,]))
bayesplot::ppc_dens_overlay(forestdata$lcarbon_point, t(test3))
bayesplot::ppc_dens_overlay(forestdata$lcarbon_point, rstan::extract(mod2_forest, pars = "cluster_lmean")[[1]])


##### MODELS
# Only random cluster effect
write("data{
              int<lower=1> n_point;                 // number of points
              int<lower=1> n_cluster;               // number of clusters
              int<lower=1> n_area;                  // number of areas
              int<lower=1> n_pred;                  // number of predictors
              vector[n_point] lcarbon_point;        // log-carbon at each point
              int<lower=1> cluster[n_point];        // cluster at each point
              vector[n_point] plotsize;             // plot size at each point (1 = 75m2, 0 = 300m2)
              matrix[n_cluster, n_pred] predictor;  // predictors at each cluster
              int<lower=1> area_cluster[n_cluster]; // area at each cluster
            }
            
            parameters{
              real alpha;                                 // global intercept
              real<lower=0> sigma_cluster;                // between-cluster sd
              vector<lower=0> [n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
              real<lower=0> lsigma_point;
            }
            
            transformed parameters{
              real cluster_mean = alpha;
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0,200);
              sigma_cluster ~ normal(0,100);
              lsigma_point ~ normal(0,10);
              // likelihood
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[i], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest_1.stan")

# Elevation covariate
write("data{
              int<lower=1> n_point;                 // number of points
              int<lower=1> n_cluster;               // number of clusters
              int<lower=1> n_area;                  // number of areas
              int<lower=1> n_pred;                  // number of predictors
              vector[n_point] lcarbon_point;        // log-carbon at each point
              int<lower=1> cluster[n_point];        // cluster at each point
              vector[n_point] plotsize;             // plot size at each point (1 = 75m2, 0 = 300m2)
              matrix[n_cluster, n_pred] predictor;  // predictors at each cluster
              int<lower=1> area_cluster[n_cluster]; // area at each cluster
            }
            
            parameters{
              real alpha;                                // global intercept
              vector[n_pred] beta;
              vector<lower=0>[n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
              real<lower=0> lsigma_point;
              real<lower=0> sigma_cluster;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + predictor * beta;
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0,200);
              beta ~ normal(0,100);
              lsigma_point ~ normal(0,10);
              sigma_cluster ~ normal(0,100);
              // likelihood
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[i], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest_2.stan")

# Elevation covariate
write("data{
              int<lower=1> n_point;                 // number of points
              int<lower=1> n_cluster;               // number of clusters
              int<lower=1> n_area;                  // number of areas
              int<lower=1> n_pred;                  // number of predictors
              vector[n_point] lcarbon_point;        // log-carbon at each point
              int<lower=1> cluster[n_point];        // cluster at each point
              vector[n_point] plotsize;             // plot size at each point (1 = 75m2, 0 = 300m2)
              matrix[n_cluster, n_pred] predictor;  // predictors at each cluster
              int<lower=1> area_cluster[n_cluster]; // area at each cluster
            }
            
            parameters{
              real alpha;                                 // global intercept
              vector[n_pred] beta;
              vector[n_pred] beta2;
              vector<lower=0> [n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
              real<lower=0> lsigma_point;
              real<lower=0> sigma_cluster;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + predictor * beta + (predictor ^ 2) * beta2;
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0,200);
              beta ~ normal(0,100);
              beta2 ~ normal(0,100);
              lsigma_point ~ normal(0,10);
              sigma_cluster ~ normal(0,100);
              // likelihood
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[i], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest_3.stan")

# Heteroscedasticity
write("data{
              int<lower=1> n_point;                 // number of points
              int<lower=1> n_cluster;               // number of clusters
              int<lower=1> n_area;                  // number of areas
              int<lower=1> n_pred;                  // number of predictors
              vector[n_point] lcarbon_point;        // log-carbon at each point
              int<lower=1> cluster[n_point];        // cluster at each point
              vector[n_point] plotsize;             // plot size at each point (1 = 75m2, 0 = 300m2)
              matrix[n_cluster, n_pred] predictor;  // predictors at each cluster
              int<lower=1> area_cluster[n_cluster]; // area at each cluster
            }
            
            parameters{
              real alpha;                                 // global intercept
              vector[n_pred] beta;
              vector<lower=0> [n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
              real<lower=0> lsigma_point;
              real<lower=0> meanfraction;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + predictor * beta;
              vector[n_cluster] sigma_cluster = cluster_mean * meanfraction;
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0,200);
              beta ~ normal(0,100);
              lsigma_point ~ normal(0,10);
              meanfraction ~ uniform(0,1);
              // likelihood
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[i], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest_4.stan")

# Heteroscedasticity
write("data{
              int<lower=1> n_point;                 // number of points
              int<lower=1> n_cluster;               // number of clusters
              int<lower=1> n_area;                  // number of areas
              int<lower=1> n_pred;                  // number of predictors
              vector[n_point] lcarbon_point;        // log-carbon at each point
              int<lower=1> cluster[n_point];        // cluster at each point
              vector[n_point] plotsize;             // plot size at each point (1 = 75m2, 0 = 300m2)
              matrix[n_cluster, n_pred] predictor;  // predictors at each cluster
              int<lower=1> area_cluster[n_cluster]; // area at each cluster
            }
            
            parameters{
              real alpha;                                 // global intercept
              vector[n_pred] beta;
              vector<lower=0> [n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
              real<lower=0> lsigma_point;
              //real<lower=0> meanfraction;
              real meanfraction;
              real<lower=0> sigma_cluster0;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + predictor * beta;
              vector[n_cluster] sigma_cluster = ((sigma_cluster0 ^ 2) * (cluster_mean ^ meanfraction)) ^ 0.5;
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0,200);
              beta ~ normal(0,100);
              lsigma_point ~ normal(0,10);
              //meanfraction ~ uniform(0,1);
              meanfraction ~ normal(0,10);
              sigma_cluster0 ~ normal(0,100);
              // likelihood
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[i], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest_5.stan")

# Plot size variance
write("data{
              int<lower=1> n_point;                 // number of points
              int<lower=1> n_cluster;               // number of clusters
              int<lower=1> n_area;                  // number of areas
              int<lower=1> n_pred;                  // number of predictors
              vector[n_point] lcarbon_point;        // log-carbon at each point
              int<lower=1> cluster[n_point];        // cluster at each point
              vector[n_point] plotsize;             // plot size at each point (1 = 75m2, 0 = 300m2)
              matrix[n_cluster, n_pred] predictor;  // predictors at each cluster
              int<lower=1> area_cluster[n_cluster]; // area at each cluster
            }
            
            parameters{
              real alpha;                                 // global intercept
              vector[n_pred] beta;
              vector<lower=0> [n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
              real<lower=0> lsigma_small;
              real<lower=0> lsigma_offset;
              real<lower=0> sigma_cluster;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + predictor * beta;
              vector[n_point] lsigma_point = lsigma_small * plotsize + lsigma_small * lsigma_offset * (1-plotsize);
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point .^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0,200);
              beta ~ normal(0,100);
              lsigma_small ~ normal(0,10);
              lsigma_offset ~ uniform(0,1);
              sigma_cluster ~ normal(0,100);
              // likelihood
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[i], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest_6.stan")

# Plot size variance + heteroscedasticity
write("data{
              int<lower=1> n_point;                 // number of points
              int<lower=1> n_cluster;               // number of clusters
              int<lower=1> n_area;                  // number of areas
              int<lower=1> n_pred;                  // number of predictors
              vector[n_point] lcarbon_point;        // log-carbon at each point
              int<lower=1> cluster[n_point];        // cluster at each point
              vector[n_point] plotsize;             // plot size at each point (1 = 75m2, 0 = 300m2)
              matrix[n_cluster, n_pred] predictor;  // predictors at each cluster
              int<lower=1> area_cluster[n_cluster]; // area at each cluster
            }
            
            parameters{
              real alpha;                                 // global intercept
              vector[n_pred] beta;
              vector<lower=0> [n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
              real<lower=0> lsigma_small;
              real<lower=0> lsigma_offset;
              real<lower=0> meanfraction;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + predictor * beta;
              vector[n_cluster] sigma_cluster = cluster_mean * meanfraction;
              vector[n_point] lsigma_point = lsigma_small * plotsize + lsigma_small * lsigma_offset * (1-plotsize);
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point .^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0,200);
              beta ~ normal(0,100);
              lsigma_small ~ normal(0,10);
              lsigma_offset ~ uniform(0,1);
              meanfraction ~ uniform(0,1);
              // likelihood
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[i], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest_7.stan")


