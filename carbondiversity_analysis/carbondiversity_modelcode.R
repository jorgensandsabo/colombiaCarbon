if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}

setwd(dir.path)

##################################
#### Carbon preparation models ###
##################################
## Model to estimate individual tree WSG
write("data{
              int<lower=1> n_trees;            // number of trees
              int<lower=1> n_cores;            // number of cores
              int<lower=1> n_spec;             // number of species
              int<lower=1> n_area;             // number of areas
              int<lower=1> n_cluster;          // number of clusters      
              int<lower=1> n_site;             // number of sites 
              vector[n_cores] WSG;             // measured WSG of cores
              int core_treenum[n_cores];       // treenumber for each core
              int species[n_trees];            // species number for each tree
              vector[n_trees] speciesid;       // tree identified?
              int site[n_trees];               // site for each tree
              int cluster_site[n_site];        // cluster at each site
              int area_cluster[n_cluster];     // area at each cluster
            }
            
            parameters{
              real<lower=0> alpha;              // intercept
              real<lower=0> sigma1;             // sigma for identified trees
              real<lower=0> sigma2;             // sigma for non-identified trees
              real mu_spec;                     // mean of species effect
              vector [n_area] b_area;           // coefficient for area effect
              real<lower=0> sigma_spec;         // standard deviation of species effect
              real<lower=0> sigma_area;         // standard deviation for area effect
              real<lower=0> sigma_cluster;      // standard deviation for cluster effect
              real<lower=0> sigma_site;         // standard deviation for site effect
              vector[n_site] sigma_site_raw;
              vector[n_cluster] sigma_cluster_raw;
              vector[n_spec] sigma_spec_raw;
            }
            
            transformed parameters{
              vector[n_cluster] b_cluster = b_area[area_cluster] + sigma_cluster * sigma_cluster_raw; // coefficient for cluster effect
              vector[n_site] b_site = b_cluster[cluster_site] + sigma_site * sigma_site_raw;          // coefficient for site effect
              vector[n_spec] b_spec = mu_spec + sigma_spec * sigma_spec_raw;                          // coefficient for species effect
              
              vector[n_trees] mu = alpha + b_spec[species] .* speciesid + b_site[site];               // estimated mean wsg
              vector[n_trees] sigma = sigma1 * speciesid + sigma2 * fabs(1-speciesid);                // estimated standard deviation
            }
            
            model{
              // priors
              alpha ~ normal(0, 30);
              sigma1 ~ normal(0, 0.25);
              sigma2 ~ normal(0, 0.25);
              
              mu_spec ~ normal(0, 0.25);
              sigma_spec ~ normal(0, 0.1);
              
              sigma_area ~ normal(0, 0.25);
              sigma_cluster ~ normal(0, 0.25);
              sigma_site ~ normal(0, 0.25);
              b_area ~ normal(0, sigma_area);
              
              sigma_cluster_raw ~ std_normal();
              sigma_site_raw ~ std_normal();
              sigma_spec_raw ~ std_normal();
              
              // likelihood
              WSG ~ normal(mu[core_treenum], sigma[core_treenum]);
            }
            
            generated quantities{
              real log_lik[n_cores];
              for(i in 1:n_cores){
                log_lik[i] = normal_lpdf(WSG[i] | mu[core_treenum[i]], sigma[core_treenum[i]]);
              }
            }
      ", file = "STAN\\carbdiv_wsgtreemodel.stan")

## Model to estimate grass biomass
write("data{
              int<lower=1> n_site;             // number of sites
              int<lower=1> n_grass;            // number of grass samples
              int<lower=1> n_habitat;          // number of habitats
              int<lower=1> n_cluster;          // number of clusters
              vector[n_grass] grass;           // measured grass AGB
              int habitat[n_site];             // habitat for each grass sample
              int grass_sitenum[n_grass];      // sitenum for each grass sample
              int cluster[n_site];             // cluster for each site
              int hab_cluster[n_cluster];      // habitat for each cluster
            }
            
            parameters{
              vector[n_habitat] alpha_hab;                  // intercept by habitat
              vector<lower=0>[n_habitat] sigma_hab;         // standard deviation by habitat
              vector<lower=0>[n_habitat] sigma_cluster_hab; // standard deviation of cluster effect by habitat
              vector[n_cluster] sigma_cluster_hab_raw;      // non-centred parametrization parameter for cluster random effect
            }
            
            transformed parameters{
              vector[n_cluster] alpha_clus = alpha_hab[hab_cluster] + sigma_cluster_hab[hab_cluster] .* sigma_cluster_hab_raw;   // cluster random effect
              vector[n_cluster] sigma_cluster = sigma_cluster_hab[hab_cluster];                                                  // standard deviation of cluster
              vector[n_site] mu = alpha_clus[cluster];                                                                           // mean of site biomass within cluster
              vector[n_site] sigma = sigma_hab[habitat];                                                                         // standard deviation of site biomass within habitat
            }
            
            model{
              // priors
              alpha_hab ~ normal(0,2);
              sigma_hab ~ normal(0,1);
              sigma_cluster_hab ~ normal(0,1);
              sigma_cluster_hab_raw ~ std_normal();
              
              // likelihood
              grass ~ normal(mu[grass_sitenum], sigma[grass_sitenum]);
            }
            
            generated quantities{
              real log_lik[n_grass];
              for(i in 1:n_grass){
                log_lik[i] = normal_lpdf(grass[i] | mu[grass_sitenum[i]], sigma[grass_sitenum[i]]);
              }
            }
      ", file = "STAN\\carbdiv_grassmodel.stan")

################################
#### Landscape carbon models ###
################################
## Forest cluster carbon model - no covariates
write("data{
              int<lower=1> n_point;                 // number of points
              int<lower=1> n_cluster;               // number of clusters
              vector[n_point] lcarbon_point;        // log-carbon at each point
              int<lower=1> cluster[n_point];        // cluster at each point
            }
            
            parameters{
              real alpha;                          // intercept
              real<lower=0> lsigma_point;          // standard deviation of plots
              real<lower=0> sigma_cluster;         // standard deviation of clusters
              vector<lower=-alpha/sigma_cluster>[n_cluster] sigma_cluster_raw; // centering parameter of cluster effect
            }
            
            transformed parameters{
              real landscape_mean = alpha;                                                           // landscape average carbon
              vector[n_cluster] carbon_cluster = landscape_mean + sigma_cluster * sigma_cluster_raw; // estimated cluster carbon stock
              vector[n_cluster] cluster_lmean = log(carbon_cluster) - (lsigma_point ^ 2)/2;          // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0,200);
              sigma_cluster_raw ~ std_normal();
              lsigma_point ~ normal(0,10);
              sigma_cluster ~ normal(0,100);
              // likelihood
              lcarbon_point ~ normal(cluster_lmean[cluster], lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[cluster[i]], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest_intercept.stan")

## Forest cluster carbon model - covariates
write("data{
              int<lower=1> n_point;                 // number of points
              int<lower=1> n_cluster;               // number of clusters
              int<lower=1> n_pred;                  // number of predictors
              vector[n_point] lcarbon_point;        // log-carbon at each point
              int<lower=1> cluster[n_point];        // cluster at each point
              matrix[n_cluster, n_pred] predictor;  // predictors at each cluster
            }
            
            parameters{
              real alpha;                         
              vector[n_pred] beta;
              real<lower=0> lsigma_point;
              real<lower=0> meanfraction;
              vector<lower=-1/meanfraction>[n_cluster] sigma_cluster_raw;
            }
            
            transformed parameters{
              vector[n_cluster] landscape_mean = alpha + predictor * beta;
              vector[n_cluster] sigma_cluster = landscape_mean * meanfraction;
              vector[n_cluster] carbon_cluster = landscape_mean + sigma_cluster .* sigma_cluster_raw;
              vector[n_cluster] cluster_lmean = log(carbon_cluster) - (lsigma_point ^ 2)/2;
            }
            
            model{
              // priors
              alpha ~ normal(0,200);
              beta ~ normal(0,100);
              sigma_cluster_raw ~ std_normal();
              //sigma_cluster_raw ~ normal(0,1);
              lsigma_point ~ normal(0,10);
              meanfraction ~ uniform(0,1);
              //meanfraction ~ normal(0.5,0.12);
              // likelihood
              lcarbon_point ~ normal(cluster_lmean[cluster], lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[cluster[i]], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest_covariates.stan")

## Paramo cluster biomass model - no covariates
write("data{
            int<lower=1> n_point;                // number of points
            int<lower=1> n_cluster;              // number of clusters
            vector[n_point] lcarbon_point;       // log-carbon at each point
            int<lower=1> cluster[n_point];       // cluster at each point
          }
          
          parameters{
            real<lower=0> alpha;                   // intercept
            real<lower=0> lsigma_cluster;          // between-cluster sd
            real<lower=0> lsigma_point;            // between-point sd
            vector[n_cluster] lsigma_cluster_raw;
          }
          
          transformed parameters{
            real<lower=0> landscape_mean = alpha;                                                     // mean landscape carbon
            vector[n_cluster] lcarbon_cluster = landscape_mean + lsigma_cluster * lsigma_cluster_raw; // carbon stocks at cluster level
          }
          
          model{
            // priors
            alpha ~ normal(0, 50);
            lsigma_cluster ~ normal(0,50);
            lsigma_cluster_raw ~ std_normal();
            lsigma_point ~ normal(0,5);
            // likelihood
            lcarbon_point ~ normal(lcarbon_cluster[cluster], lsigma_point);
          }
          
          generated quantities{
            vector[n_cluster] carbon_cluster = exp(lcarbon_cluster + (lsigma_point ^ 2) / 2);
            real log_lik[n_point];
            for(i in 1:n_point){
              log_lik[i] = normal_lpdf(lcarbon_point[i] | lcarbon_cluster[cluster[i]], lsigma_point);
            }
          }", file = "STAN\\carbdiv_clus_paramo_intercept.stan")

## Paramo cluster biomass model - covariates
write("data{
            int<lower=1> n_point;                // number of points
            int<lower=1> n_cluster;              // number of clusters
            int<lower=1> n_pred;                 // number of predictors
            vector[n_point] lcarbon_point;       // log-carbon at each point
            int<lower=1> cluster[n_point];       // cluster at each point
            matrix[n_cluster, n_pred] predictor; // predictors at each cluster
          }
          
          parameters{
            real<lower=0> alpha;                   // intercept
            vector[n_pred] beta;                 // coefficient for predictors
            real<lower=0> lsigma_cluster;          // between-cluster sd
            real<lower=0> lsigma_point;            // between-point sd
            vector[n_cluster] lsigma_cluster_raw;
          }
          
          transformed parameters{
            vector[n_cluster] landscape_mean = alpha + predictor * beta;                              // mean landscape carbon
            vector[n_cluster] lcarbon_cluster = landscape_mean + lsigma_cluster * lsigma_cluster_raw; // carbon stocks at cluster level
          }
          
          model{
            // priors
            alpha ~ normal(0, 50);
            beta ~ normal(0, 50);
            lsigma_cluster ~ normal(0,50);
            lsigma_cluster_raw ~ std_normal();
            lsigma_point ~ normal(0,5);
            // likelihood
            lcarbon_point ~ normal(lcarbon_cluster[cluster], lsigma_point);
          }
          
          generated quantities{
            vector[n_cluster] carbon_cluster = exp(lcarbon_cluster + (lsigma_point ^ 2) / 2);
            real log_lik[n_point];
            for(i in 1:n_point){
              log_lik[i] = normal_lpdf(lcarbon_point[i] | lcarbon_cluster[cluster[i]], lsigma_point);
            }
          }", file = "STAN\\carbdiv_clus_paramo_covariates.stan")

#########################
#### Diversity models ###
#########################
## Diversity models
write("data{
              int<lower=1> n_cluster;            // number of clusters
              vector[n_cluster] predictor_cluster;  // predictor at each cluster
              vector[n_cluster] div;             // diversity metric at each cluster
              int<lower=1> biogeo[n_cluster];  // biogeographic region of each cluster
              int<lower=1> n_biogeo;           // number of biogeographic regions
            }
            
            parameters{
              real alpha;
              vector [n_biogeo] alpha_biogeo;
              real<lower=0> sigma_cluster;
              real beta1;
              real beta2;
              vector[n_biogeo] beta1_biogeo;
              vector[n_biogeo] beta2_biogeo;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha_biogeo[biogeo] + beta1_biogeo[biogeo] .* predictor_cluster + beta2_biogeo[biogeo] .* (predictor_cluster .^ 2);
            }
            
            model{
              // priors
              alpha ~ normal(0, 100);
              beta1 ~ normal(0, 50);
              beta2 ~ normal(0, 50);
              sigma_cluster ~ normal(0,100);
              // likelihood
              alpha_biogeo ~ normal(alpha, 100);
              beta1_biogeo ~ normal(beta1, 50);
              beta2_biogeo ~ normal(beta2, 50);
              div ~ normal(cluster_mean, sigma_cluster);
            }
            
            generated quantities{
              real log_lik[n_cluster];
              for(i in 1:n_cluster){
                log_lik[i] = normal_lpdf(div[i] | cluster_mean[i], sigma_cluster);
              }
            }", file = "STAN\\carbdiv_poly_regslope.stan")

write("data{
              int<lower=1> n_cluster;            // number of clusters
              vector[n_cluster] predictor_cluster;  // predictor at each cluster
              vector[n_cluster] div;             // diversity metric at each cluster
              int<lower=1> biogeo[n_cluster];  // biogeographic region of each cluster
              int<lower=1> n_biogeo;           // number of biogeographic regions
            }
            
            parameters{
              real alpha;
              vector[n_biogeo] alpha_biogeo;
              real<lower=0> sigma_cluster;
              real beta1;
              real beta2;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha_biogeo[biogeo] + beta1 * predictor_cluster + beta2 * (predictor_cluster .^ 2);
            }
            
            model{
              // priors
              alpha ~ normal(0, 100);
              beta1 ~ normal(0, 50);
              beta2 ~ normal(0, 50);
              sigma_cluster ~ normal(0,100);
              // likelihood
              alpha_biogeo ~ normal(alpha, 100);
              div ~ normal(cluster_mean, sigma_cluster);
            }
            
            generated quantities{
              real log_lik[n_cluster];
              for(i in 1:n_cluster){
                 log_lik[i] = normal_lpdf(div[i] | cluster_mean[i], sigma_cluster);
              }
            }", file = "STAN\\carbdiv_poly_regintercept.stan")

write("data{
            int<lower=1> n_cluster;            // number of clusters
            vector[n_cluster] predictor_cluster;  // predictor at each cluster
            vector[n_cluster] div;             // diversity metric at each cluster
           }
          
           parameters{
            real alpha;
            real<lower=0> sigma_cluster;
            real beta1;
            real beta2;
           }
          
           transformed parameters{
            vector[n_cluster] cluster_mean = alpha + beta1 * predictor_cluster + beta2 * (predictor_cluster .^ 2);
           }
          
           model{
            // priors
            alpha ~ normal(0, 100);
            beta1 ~ normal(0, 50);
            beta2 ~ normal(0, 50);
            sigma_cluster ~ normal(0,100);
            // likelihood
            div ~ normal(cluster_mean, sigma_cluster);
           }
      
          generated quantities{
            real log_lik[n_cluster];
            for(i in 1:n_cluster){
              log_lik[i] = normal_lpdf(div[i] | cluster_mean[i], sigma_cluster);
            }
          }", file = "STAN\\carbdiv_poly_intercept.stan")



write("data{
              int<lower=1> n_cluster;            // number of clusters
              vector[n_cluster] predictor_cluster;  // predictor at each cluster
              vector[n_cluster] div;             // diversity metric at each cluster
              int<lower=1> biogeo[n_cluster];  // biogeographic region of each cluster
              int<lower=1> n_biogeo;           // number of biogeographic regions
            }
            
            parameters{
              real alpha;
              vector [n_biogeo] alpha_biogeo;
              real<lower=0> sigma_cluster;
              real beta1;
              vector[n_biogeo] beta1_biogeo;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha_biogeo[biogeo] + beta1_biogeo[biogeo] .* predictor_cluster;
            }
            
            model{
              // priors
              alpha ~ normal(0, 100);
              beta1 ~ normal(0, 50);
              sigma_cluster ~ normal(0,100);
              // likelihood
              alpha_biogeo ~ normal(alpha, 100);
              beta1_biogeo ~ normal(beta1, 50);
              div ~ normal(cluster_mean, sigma_cluster);
            }
            
            generated quantities{
              real log_lik[n_cluster];
              for(i in 1:n_cluster){
                log_lik[i] = normal_lpdf(div[i] | cluster_mean[i], sigma_cluster);
            }}", file = "STAN\\carbdiv_lin_regslope.stan")


write("data{
              int<lower=1> n_cluster;            // number of clusters
              vector[n_cluster] predictor_cluster;  // predictor at each cluster
              vector[n_cluster] div;             // diversity metric at each cluster
              int<lower=1> biogeo[n_cluster];  // biogeographic region of each cluster
              int<lower=1> n_biogeo;           // number of biogeographic regions
            }
            
            parameters{
              real alpha;
              vector[n_biogeo] alpha_biogeo;
              real<lower=0> sigma_cluster;
              real beta1;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha_biogeo[biogeo] + beta1 * predictor_cluster;
            }
            
            model{
              // priors
              alpha ~ normal(0, 100);
              beta1 ~ normal(0, 50);
              sigma_cluster ~ normal(0,100);
              // likelihood
              alpha_biogeo ~ normal(alpha, 100);
              div ~ normal(cluster_mean, sigma_cluster);
            }
            
            generated quantities{
              real log_lik[n_cluster];
              for(i in 1:n_cluster){
                log_lik[i] = normal_lpdf(div[i] | cluster_mean[i], sigma_cluster);
            }}", file = "STAN\\carbdiv_lin_regintercept.stan")

write("data{
              int<lower=1> n_cluster;            // number of clusters
              vector[n_cluster] predictor_cluster;  // predictor at each cluster
              vector[n_cluster] div;             // diversity metric at each cluster
              int<lower=1> biogeo[n_cluster];  // biogeographic region of each cluster
              int<lower=1> n_biogeo;           // number of biogeographic regions
            }
            
            parameters{
              real alpha;
              real<lower=0> sigma_cluster;
              real beta1;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + beta1 * predictor_cluster;
            }
            
            model{
              // priors
              alpha ~ normal(0, 100);
              beta1 ~ normal(0, 50);
              sigma_cluster ~ normal(0,100);
              // likelihood
              div ~ normal(cluster_mean, sigma_cluster);
            }
            
            generated quantities{
              real log_lik[n_cluster];
              for(i in 1:n_cluster){
                log_lik[i] = normal_lpdf(div[i] | cluster_mean[i], sigma_cluster);
            }}", file = "STAN\\carbdiv_lin_intercept.stan")

#################################
#### Carbon + elevation model ###
#################################
write("data{
              int<lower=1> n_cluster;            // number of clusters
              vector[n_cluster] carbon_cluster;  // carbon at each cluster
              vector[n_cluster] div;             // diversity metric at each cluster
              vector[n_cluster] elevation_cluster;
              int<lower=1> biogeo[n_cluster];  // biogeographic region of each cluster
              int<lower=1> n_biogeo;           // number of biogeographic regions
            }
            
            parameters{
              real alpha;
              vector[n_biogeo] alpha_biogeo;
              real<lower=0> sigma_cluster;
              real beta1;
              real beta2;
              real betaelev1;
              real betaelev2;
              vector[n_biogeo] beta1_biogeo;
              vector[n_biogeo] beta2_biogeo;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha_biogeo[biogeo] + beta1_biogeo[biogeo] .* carbon_cluster + beta2_biogeo[biogeo] .* (carbon_cluster .^ 2) + betaelev1 * elevation_cluster + betaelev2 * (elevation_cluster .^ 2);
            }
            
            model{
              // priors
              alpha ~ normal(0,100);
              beta1 ~ normal(0,50);
              beta2 ~ normal(0,50);
              betaelev1 ~ normal(0,50);
              betaelev2 ~ normal(0,50);
              sigma_cluster ~ normal(0,100);
              // likelihood
              beta1_biogeo ~ normal(beta1,50);
              beta2_biogeo ~ normal(beta2, 50);
              alpha_biogeo ~ normal(alpha, 100);
              div ~ normal(cluster_mean, sigma_cluster);
            }
            
            generated quantities{
              real log_lik[n_cluster];
              for(i in 1:n_cluster){
                 log_lik[i] = normal_lpdf(div[i] | cluster_mean[i], sigma_cluster);
              }
            }", file = "STAN\\carb_elevation_div_model.stan")

