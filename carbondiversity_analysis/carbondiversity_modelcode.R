

if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}

setwd(dir.path)

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
              vector[n_habitat] alpha_hab;         // intercept by habitat
              vector<lower=0>[n_habitat] sigma_hab;         // standard deviation by habitat
              vector<lower=0>[n_habitat] sigma_cluster_hab; // standard deviation of cluster effect by habitat
              vector[n_cluster] alpha_clus;        // intercept by cluster within habitat
            }
            
            transformed parameters{
              vector[n_cluster] sigma_cluster = sigma_cluster_hab[hab_cluster];  // standard deviation of cluster
              vector[n_site] mu = alpha_clus[cluster];                           // mean of site biomass within cluster
              vector[n_site] sigma = sigma_hab[habitat];                         // standard deviation of site biomass within habitat
            }
            
            model{
              // priors
              alpha_hab ~ normal(0,30);
              sigma_hab ~ normal(0, 10);
              sigma_cluster_hab ~ normal(0, 10);
              alpha_clus ~ normal(alpha_hab[hab_cluster], sigma_cluster);
              
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

## Model to estimate cluster biomass - forest - WIGH CHOCO (plot sizes)
write("data{
              int<lower=1> n_point;                 // number of points
              int<lower=1> n_cluster;               // number of clusters
              int<lower=1> n_area;                  // number of areas
              int<lower=1> n_pred;                  // number of predictors
              vector[n_point] lcarbon_point;        // log-carbon at each point
              int<lower=1> cluster[n_point];        // cluster at each point
              //vector[n_point] plotsize;             // plot size at each point (1 = 75m2, 0 = 300m2)
              matrix[n_cluster, n_pred] predictor;  // predictors at each cluster
              int<lower=1> area_cluster[n_cluster]; // area at each cluster
            }
            
            parameters{
              real alpha;                        // global intercept
              vector[n_area] alpha_area;         // intercept by area
              vector[n_pred] beta;                        // coefficient for each predictor
              real<lower=0> sigma_cluster;                // between-cluster sd
              //vector<lower=0>[n_cluster] sigma_cluster;
              real<lower=0> sigma_area;                   // standard deviation of area effect
              real<lower=0> lsigma_small;                 // between-point sd for small plots
              real<lower=0> lsigma_offset;                // large plot proportional sd of small plots
              vector<lower=0> [n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha_area[area_cluster] + predictor * beta;
              vector[n_point] lsigma_point = lsigma_small * plotsize + lsigma_small * lsigma_offset * (1 - plotsize);
              //vector[n_point] lsigma_point = lsigma_small[cluster] .* plotsize + lsigma_small[cluster] .* (lsigma_offset * (1 - plotsize));
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point .^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(200,50);
              sigma_area ~ normal(0,50);
              beta ~ normal(0,10);
              sigma_cluster ~ normal(0,50);
              lsigma_small ~ normal(0,1);
              lsigma_offset ~ normal(0,0.5);
              // likelihood
              alpha_area ~ normal(alpha, sigma_area);
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
      
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[i], lsigma_point[i]);
              }
            }", file = "STAN\\carbdiv_clus_forest.stan")

## Model to estimate cluster biomass - forest - NOCHOCO
write("data{
              int<lower=1> n_point;                 // number of points
              int<lower=1> n_cluster;               // number of clusters
              int<lower=1> n_area;                  // number of areas
              int<lower=1> n_pred;                  // number of predictors
              vector[n_point] lcarbon_point;        // log-carbon at each point
              int<lower=1> cluster[n_point];        // cluster at each point
              matrix[n_cluster, n_pred] predictor;  // predictors at each cluster
              int<lower=1> area_cluster[n_cluster]; // area at each cluster
            }
            
            parameters{
              real alpha;                                 // global intercept
              vector[n_area] alpha_area;                  // intercept by area
              vector[n_pred] beta;                        // coefficient for each predictor
              real<lower=0> sigma_cluster1;                // between-cluster sd
              real <lower=0> b;
              real<lower=0> sigma_area;                   // standard deviation of area effect
              vector<lower=0> [n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
              real<lower=0> lsigma_point;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha_area[area_cluster] + predictor * beta;
              vector[n_cluster] sigma_cluster = sigma_cluster1 + b * cluster_mean;
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point .^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(200,50);
              sigma_area ~ normal(0,50);
              beta ~ normal(0,10);
              sigma_cluster1 ~ normal(0,30);
              b ~ normal(0,30);
              lsigma_point ~ normal(0,10);
              // likelihood
              alpha_area ~ normal(alpha, sigma_area);
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }
      
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[i], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest2.stan")

## Model to estimate cluster biomass - paramo
write("data{
              int<lower=1> n_point;                // number of points
              int<lower=1> n_cluster;              // number of clusters
              int<lower=1> n_area;                 // number of areas
              int<lower=1> n_pred;                 // number of predictors
              vector[n_point] lcarbon_point;       // log-carbon at each point
              int<lower=1> cluster[n_point];       // cluster at each point
              //int<lower=1> area_cluster[n_cluster];// area at each cluster
              matrix[n_cluster, n_pred] predictor; // predictors at each cluster
            }
            
            parameters{
              real<lower=0> alpha;               // intercept
              //vector<lower=0>[n_area] alpha_area;// intercept for areas
              //vector[n_pred] beta;               // coefficient for predictors
              //real<lower=0> sigma_area;          // standard deviation for areas
              real<lower=0> sigma_cluster;       // between-cluster sd by habitat (should this grow with the mean?)
              real<lower=0> lsigma_point;        // between-point sd by habitat for large plots
              //vector[n_cluster] sigma_c_raw;
              vector<lower=0>[n_cluster] carbon_cluster;
              //vector[n_area] sigma_area_raw;
            }
            
            transformed parameters{
              //vector[n_area] alpha_area = alpha + sigma_area * sigma_area_raw;
              //vector[n_cluster] cluster_mean = alpha_area[area_cluster] + predictor * beta;
              //vector[n_cluster] cluster_mean = alpha + predictor * beta;
              real cluster_mean = alpha;
              //vector[n_cluster] carbon_cluster = cluster_mean + sigma_cluster * sigma_c_raw;
              vector[n_cluster] cluster_lmean = log(carbon_cluster) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0, 50);
             // sigma_area ~ normal(0, 50);
            //  sigma_area_raw ~ std_normal();
             // beta ~ normal(0, 10);
              sigma_cluster ~ normal(0,20);
              //sigma_c_raw ~ std_normal();
              lsigma_point ~ normal(0,5);
              // likelihood
              //alpha_area ~ normal(alpha, sigma_area);
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean[cluster], lsigma_point);
            }", file = "STAN\\carbdiv_clus_paramo.stan")










