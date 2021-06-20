data{
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
            }
