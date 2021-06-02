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
              real alpha;                        // global intercept
              //vector[n_area] alpha_area;         // intercept by area
              vector[n_pred] beta;                        // coefficient for each predictor
              real<lower=0> sigma_cluster;                // between-cluster sd
              //real<lower=0> sigma_area;                   // standard deviation of area effect
              real<lower=0> lsigma_small;                 // between-point sd for small plots
              real<lower=0> lsigma_offset;                // large plot proportional sd of small plots
              vector<lower=0> [n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
              //real<lower=0> lsigma_point;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + predictor * beta;
              vector[n_point] lsigma_point = lsigma_small * plotsize + lsigma_small * lsigma_offset * (1 - plotsize);
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point .^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0,200);
              //sigma_area ~ normal(0,50);
              beta ~ normal(0,10);
              sigma_cluster ~ normal(0,50);
              lsigma_small ~ normal(0,1);
              lsigma_offset ~ normal(0,0.5);
              //lsigma_point ~ normal(0,1);
              // likelihood
              //alpha_area ~ normal(alpha, sigma_area);
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[i], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest.stan")

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
              real alpha;                           // intercept
              //vector[n_area] alpha_area;          // intercept for areas
              vector[n_cluster] alpha_cluster;
              vector[n_pred] beta;                  // coefficient for predictors
              //real<lower=0> sigma_area;           // standard deviation for areas
              real<lower=0> sigma_cluster;
              //real<lower=0> lsigma_small;         // between-plot sd for small plots
              //real<lower=0> lsigma_offset;        // proportional offset for large plots
              real<lower=0> lsigma_point;
            }
            
            transformed parameters{
              //vector[n_cluster] cluster_lmean = alpha_area[area_cluster] + predictor * beta;
              vector[n_cluster] cluster_lmean = alpha_cluster + predictor * beta;
              //vector[n_point] lsigma_point = lsigma_small * plotsize + lsigma_small * lsigma_offset * (1 - plotsize);
              vector[n_point] carbon_cluster = cluster_lmean[cluster] + (lsigma_point ^ 2)/2;
            }
            
            model{
              // priors
              alpha ~ normal(0,20);
              //sigma_area ~ normal(0,1);
              sigma_cluster ~ normal(0,1);
              beta ~ normal(0,5);
              //lsigma_offset ~ normal(0,0.5);
              //lsigma_small ~ normal(0,1);
              lsigma_point ~ normal(0,1);
              // likelihood
              //alpha_area ~ normal(alpha, sigma_area);
              alpha_cluster ~ normal(alpha, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean[cluster], lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[cluster[i]], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_forest2.stan")
            
write("data{
              int<lower=1> n_point;                // number of points
              int<lower=1> n_cluster;              // number of clusters
              int<lower=1> n_area;                 // number of areas
              int<lower=1> n_pred;                 // number of predictors
              vector[n_point] lcarbon_point;       // log-carbon at each point
              int<lower=1> cluster[n_point];       // cluster at each point
              int<lower=1> area_cluster[n_cluster];// area at each cluster
              matrix[n_cluster, n_pred] predictor; // predictors at each cluster
            }
            
            parameters{
              real alpha;                        // intercept
              //vector[n_area] alpha_area;       // intercept for areas
              vector[n_pred] beta;               // coefficient for predictors
              //real<lower=0> sigma_area;        // standard deviation for areas
              real<lower=0> sigma_cluster;       // between-cluster sd by habitat (should this grow with the mean?)
              real<lower=0> lsigma_point;        // between-point sd by habitat for large plots
              vector[n_cluster] sigma_c_raw;
              //vector<lower=0>[n_cluster] carbon_cluster;
              //vector[n_area] sigma_area_raw;
            }
            
            transformed parameters{
              //vector[n_area] alpha_area = alpha + sigma_area * sigma_area_raw;
              //vector[n_cluster] cluster_mean = alpha_area[area_cluster] + predictor * beta;
              vector[n_cluster] cluster_mean = alpha + predictor * beta;
              vector[n_cluster] carbon_cluster = cluster_mean + sigma_cluster * sigma_c_raw;
              vector[n_cluster] cluster_lmean = log(carbon_cluster) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0, 30);
              //sigma_area ~ normal(0, 20);
              //sigma_area_raw ~ std_normal();
              beta ~ normal(0, 10);
              sigma_cluster ~ normal(0,20);
              sigma_c_raw ~ std_normal();
              lsigma_point ~ normal(0,5);
              // likelihood
              //alpha_area ~ normal(alpha, sigma_area);
              //carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean[cluster], lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[cluster[i]], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_paramo.stan")

write("data{
              int<lower=1> n_point;            // number of points
              int<lower=1> n_cluster;          // number of clusters
              int<lower=1> n_area;             // number of areas
              int<lower=1> n_pred;             // number of predictors
              vector[n_point] lcarbon_point;   // log-carbon at each point
              int<lower=1> cluster[n_point];   // cluster at each point
              int<lower=1> area_cluster[n_cluster];// area at each cluster
              matrix[n_cluster, n_pred] predictor; // predictors at each cluster
            }
            
            parameters{
              real alpha;                 // intercept by habitat
              //vector[n_area] alpha_area;// intercept for areas
              vector[n_cluster] alpha_cluster;
              real<lower=0> sigma_cluster;
              //real<lower=0> sigma_area;          // standard deviation for areas
              vector[n_pred] beta;                 // coefficient for each predictor
              real<lower=0> lsigma_point;          // between-plot sd
              //vector[n_area] sigma_area_raw;
            }
            
            transformed parameters{
              //vector[n_area] alpha_area = alpha + sigma_area * sigma_area_raw;
              vector[n_cluster] cluster_lmean = alpha_cluster + predictor * beta;
              vector[n_cluster] carbon_cluster = cluster_lmean + (lsigma_point ^ 2)/2;
            }
            
            model{
              // priors
              alpha ~ normal(0,5);
              sigma_cluster ~ normal(0,20);
              //sigma_area ~ normal(0,20);
              //sigma_area_raw ~ std_normal();
              lsigma_point ~ normal(0,5);
              beta ~ normal(0,5);
              // likelihood
              //alpha_area ~ normal(alpha, sigma_area);
              alpha_cluster ~ normal(alpha, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean[cluster], lsigma_point);
            }
                  
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[cluster[i]], lsigma_point);
              }
            }", file = "STAN\\carbdiv_clus_paramo2.stan")



