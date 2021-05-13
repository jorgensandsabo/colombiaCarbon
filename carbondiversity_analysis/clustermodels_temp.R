write("data{
              int<lower=1> n_point;            // number of points
              int<lower=1> n_cluster;          // number of clusters
              vector[n_point] lcarbon_point;   // log-carbon at each point
              int<lower=1> cluster[n_point];   // cluster at each point
              vector[n_point] plotsize;        // plot size at each point (1 = 75m2, 0 = 300m2)
              vector[n_cluster] elevation;     // elevation at each point
              int<lower=1> n_area;
              int<lower=1> area_cluster[n_cluster];
            }
            
            parameters{
              real<lower=0> alpha;                        // global intercept
              vector<lower=0>[n_area] alpha_area;         // intercept by area
              real beta_elev;                             // coefficient for elevation
              real<lower=0> sigma_cluster_zero;           // between-cluster sd at 0 carbon
              real<lower=0> powvar;                       // power mean-variance term
              //real<lower=0> powvar2;
              real<lower=0> sigma_area;                   // standard deviation of area effect
              real<lower=0> lsigma_small;                 // between-point sd for small plots
              real<lower=0> lsigma_offset;                // large plot proportional sd of small plots
              vector<lower=0> [n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha_area[area_cluster] + beta_elev * elevation;
              //vector[n_cluster] sigma_cluster = ((sigma_cluster_zero ^ 2) * (cluster_mean .^ (2 * powvar))) .^ 0.5;  // VarPower
              //vector[n_cluster] sigma_cluster = ((sigma_cluster_zero ^ 2) * exp(2 * powvar * cluster_mean)) .^ 0.5; // VarExp
              //vector[n_cluster] sigma_cluster = ((sigma_cluster_zero ^ 2) * ((powvar2 + (cluster_mean .^ powvar)) ^ 2)) .^ 0.5;
              vector[n_cluster] sigma_cluster = sigma_cluster_zero + powvar * cluster_mean; // SD linear;
              vector[n_point] lsigma_point = lsigma_small * plotsize + lsigma_small * lsigma_offset * (1 - plotsize);
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(200,50);
              sigma_area ~ normal(0,50);
              beta_elev ~ normal(0,10);
              sigma_cluster_zero ~ normal(0,20);
              lsigma_small ~ normal(0,1);
              lsigma_offset ~ normal(0,0.5);
              powvar ~ normal(0,1);
              //powvar2 ~ normal(0,1);
              // likelihood
              alpha_area ~ normal(alpha, sigma_area);
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }", file = "STAN\\carbdiv_clus_forest1.stan")

write("data{
              int<lower=1> n_point;            // number of points
              int<lower=1> n_cluster;          // number of clusters
              vector[n_point] lcarbon_point;   // log-carbon at each point
              int<lower=1> cluster[n_point];   // cluster at each point
              vector[n_point] plotsize;        // plot size at each point (1 = 75m2, 0 = 300m2)
              vector[n_cluster] elevation;     // elevation at each cluster (scaled)
              int<lower=1> n_area;
              int<lower=1> area_cluster[n_cluster];
            }
            
            parameters{
              real<lower=0> alpha;              // intercept
              real beta_elev;                   // coefficient for elevation
              real<lower=0> lsigma_small;       // between-plot sd for small plots
              real<lower=0> lsigma_offset;      // proportional offset for large plots
              real<lower=0> sigma_area;
              vector[n_area] alpha_area;
            }
            
            transformed parameters{
              vector[n_point] sigma_cluster = sigma_cluster_zero + powvar * cluster_mean; // SD linear;
              vector[n_point] cluster_lmean = alpha_area[area_cluster][cluster] + beta_elev * elevation[cluster];
              vector[n_point] lsigma_point = lsigma_small * plotsize + lsigma_small * lsigma_offset * (1 - plotsize);
              vector[n_point] carbon_cluster = cluster_lmean + (lsigma_point[cluster] .^ 2)/2;
            }
            
            model{
              // priors
              alpha ~ normal(0,5);
              sigma_area ~ normal(0, 20);
              lsigma_offset ~ uniform(0,1);
              lsigma_small ~ normal(0,1);
              beta_elev ~ normal(0,5);
              // likelihood
              alpha_area ~ normal(alpha, sigma_area);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point[cluster]);
            }", file = "STAN\\carbdiv_clus_forest2.stan")
            
write("data{
              int<lower=1> n_point;            // number of points
              int<lower=1> n_cluster;          // number of clusters
              vector[n_point] lcarbon_point;   // log-carbon at each point
              int<lower=1> cluster[n_point];   // cluster at each point'
              vector[n_cluster] elevation;     // elevation at each point
            }
            
            parameters{
              real<lower=0> alpha;              // intercept
              real beta_elev;                   // coefficient for elevation
              real<lower=0> sigma_cluster;      // between-cluster sd by habitat (should this grow with the mean?)
              real<lower=0> lsigma_point;       // between-point sd by habitat for large plots
              vector [n_cluster] sigma_c_raw;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + beta_elev * elevation;
              vector[n_cluster] carbon_cluster = cluster_mean + sigma_cluster * sigma_c_raw;
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0, 30);
              beta_elev ~ normal(0, 5);
              sigma_cluster ~ normal(0,20);
              sigma_c_raw ~ std_normal();
              lsigma_point ~ normal(0,5);
              // likelihood
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }", file = "STAN\\carbdiv_clus_paramo1.stan")

write("data{
              int<lower=1> n_point;            // number of points
              int<lower=1> n_cluster;          // number of clusters
              vector[n_point] lcarbon_point;   // log-carbon at each point
              int<lower=1> cluster[n_point];   // cluster at each point
              vector[n_cluster] elevation;     // elevation at each cluster (scaled)
            }
            
            parameters{
              real<lower=0> alpha;              // intercept by habitat
              real beta_elev;                   // coefficient for elevation by habitat
              real<lower=0> lsigma_point;       // between-plot sd
            }
            
            transformed parameters{
              vector[n_point] cluster_lmean = alpha + beta_elev * elevation[cluster];
              vector[n_point] carbon_cluster = cluster_lmean + (lsigma_point ^ 2)/2;
            }
            
            model{
              // priors
              alpha ~ normal(0,5);
              lsigma_point ~ normal(0,5);
              beta_elev ~ normal(0,5);
              // likelihood
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }", file = "STAN\\carbdiv_clus_paramo2.stan")










write("data{
              int<lower=1> n_point;            // number of points
              int<lower=1> n_cluster;          // number of clusters
              vector[n_point] lcarbon_point;    // carbon at each point
              int<lower=1> cluster[n_point];   // cluster at each point'
              vector[n_cluster] elevation;     // elevation at each point
            }
            
            parameters{
              real<lower=0> alpha;              // intercept by habitat
              real beta_elev;                   // coefficient for elevation
              real<lower=0> sigma_cluster;      // between-cluster sd by habitat (should this grow with the mean?)
              real<lower=0> sigma_point;       // between-point sd by habitat for large plots
              vector<lower=0> [n_cluster] carbon_cluster;     // mean carbon stocks at the cluster scale
              //vector [n_cluster] sigma_c_raw;
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha + beta_elev * elevation;
              //vector[n_cluster] carbon_cluster = cluster_mean + sigma_cluster * sigma_c_raw;
              //vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0, 30);
              beta_elev ~ normal(0, 5);
              sigma_cluster ~ normal(0,20);
              //sigma_c_raw ~ std_normal();
              sigma_point ~ normal(0,5);
              // likelihood
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(carbon_cluster[cluster], sigma_point);
            }", file = "STAN\\carbdiv_clus_paramo1.stan")


