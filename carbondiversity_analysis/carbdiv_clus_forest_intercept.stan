data{
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
            }
