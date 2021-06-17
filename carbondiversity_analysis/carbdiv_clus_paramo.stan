data{
              int<lower=1> n_point;                // number of points
              int<lower=1> n_cluster;              // number of clusters
              int<lower=1> n_pred;                 // number of predictors
              vector[n_point] lcarbon_point;       // log-carbon at each point
              int<lower=1> cluster[n_point];       // cluster at each point
              matrix[n_cluster, n_pred] predictor; // predictors at each cluster
            }
            
            parameters{
              real alpha;                          // intercept
              vector[n_pred] beta;                 // coefficient for predictors
              real<lower=0> sigma_cluster;         // between-cluster sd by habitat (should this grow with the mean?)
              real<lower=0> lsigma_point;          // between-point sd by habitat for large plots
              vector[n_cluster] sigma_cluster_raw;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + predictor * beta;                            // mean landscape carbon
              vector[n_cluster] carbon_cluster = cluster_mean + sigma_cluster * sigma_cluster_raw;  // carbon stocks at cluster level
              vector[n_cluster] cluster_lmean = log(carbon_cluster) - (lsigma_point ^ 2)/2;         // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0, 50);
              beta ~ normal(0, 50);
              sigma_cluster ~ normal(0,50);
              sigma_cluster_raw ~ std_normal();
              lsigma_point ~ normal(0,5);
              // likelihood
              lcarbon_point ~ normal(cluster_lmean[cluster], lsigma_point);
            }
            
            generated quantities{
              real log_lik[n_point];
              for(i in 1:n_point){
                log_lik[i] = normal_lpdf(lcarbon_point[i] | cluster_lmean[cluster[i]], lsigma_point);
              }
            }
