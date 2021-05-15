

if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}

setwd(dir.path)

## Model to estimate individual tree WSG
write("data{
              int<lower=1> n_trees;            // number of trees
              int<lower=1> n_cores;            // number of cores
              int<lower=1> n_spec;             // number of species
              int<lower=1> n_area;             // number of areas
              int<lower=1> n_biogeo;           // number of biogeographic regions
              int<lower=1> n_cluster;          // number of clusters      
              int<lower=1> n_site;             // number of sites 
              vector[n_cores] WSG;             // measured WSG of cores
              int core_treenum[n_cores];       // treenumber for each core
              int species[n_trees];            // species number for each tree
              vector [n_trees] speciesid;      // tree identified?
              int site[n_trees];               // site for each tree
              int cluster_site[n_site];        // cluster at each site
              int area_cluster[n_cluster];     // area at each cluster
              int biogeo_area[n_area];         // biogeographic region at each area
            }
            
            parameters{
              real<lower=0> alpha;              // intercept
              real<lower=0> sigma1;             // sigma for identified trees
              real<lower=0> sigma2;             // sigma for non-identified trees
              real mu_spec;                     // mean of species effect
              vector [n_biogeo] b_biogeo;       // coefficient for biogeographic region effect
              real<lower=0> sigma_spec;         // standard deviation of species effect
              real<lower=0> sigma_biogeo;       // standard deviation for region effect
              real<lower=0> sigma_area;         // standard deviation for area effect
              real<lower=0> sigma_cluster;      // standard deviation for cluster effect
              real<lower=0> sigma_site;         // standard deviation for site effect
              vector [n_site] sigma_site_raw;
              vector [n_cluster] sigma_cluster_raw;
              vector [n_area] sigma_area_raw;
              vector [n_spec] sigma_spec_raw;
            }
            
            transformed parameters{
              vector[n_area] b_area = b_biogeo[biogeo_area] + sigma_area * sigma_area_raw;            // coefficient for area effect
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
              
              sigma_biogeo ~ normal(0, 0.25);
              sigma_area ~ normal(0, 0.25);
              sigma_cluster ~ normal(0, 0.25);
              sigma_site ~ normal(0, 0.25);
              b_biogeo ~ normal(0, sigma_biogeo);
              
              sigma_area_raw ~ std_normal();
              sigma_cluster_raw ~ std_normal();
              sigma_site_raw ~ std_normal();
              sigma_spec_raw ~ std_normal();
              
              // likelihood
              WSG ~ normal(mu[core_treenum], sigma[core_treenum]);
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
              vector[n_habitat] sigma_hab;         // standard deviation by habitat
              vector[n_habitat] sigma_cluster_hab; // standard deviation of cluster effect by habitat
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
      ", file = "STAN\\carbdiv_grassmodel.stan")

## Model to estimate cluster biomass - forest
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
              real<lower=0> alpha;                        // global intercept
              vector<lower=0>[n_area] alpha_area;         // intercept by area
              vector[n_pred] beta;                        // coefficient for each predictor
              real<lower=0> sigma_cluster;                // between-cluster sd
              real<lower=0> sigma_area;                   // standard deviation of area effect
              real<lower=0> lsigma_small;                 // between-point sd for small plots
              real<lower=0> lsigma_offset;                // large plot proportional sd of small plots
              vector<lower=0> [n_cluster] carbon_cluster; // mean carbon stocks at the cluster scale
            }
            
            transformed parameters{
              vector[n_cluster] cluster_mean = alpha_area[area_cluster] + predictor * beta;
              vector[n_point] lsigma_point = lsigma_small * plotsize + lsigma_small * lsigma_offset * (1 - plotsize);
              vector[n_point] cluster_lmean = log(carbon_cluster[cluster]) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(200,50);
              sigma_area ~ normal(0,50);
              beta ~ normal(0,10);
              sigma_cluster ~ normal(0,20);
              lsigma_small ~ normal(0,1);
              lsigma_offset ~ normal(0,0.5);
              // likelihood
              alpha_area ~ normal(alpha, sigma_area);
              carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean, lsigma_point);
            }", file = "STAN\\carbdiv_clus_forest.stan")

## Model to estimate cluster biomass - paramo
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
              real<lower=0> alpha;               // intercept
              //vector<lower=0>[n_area] alpha_area;// intercept for areas
              vector[n_pred] beta;               // coefficient for predictors
              real<lower=0> sigma_area;          // standard deviation for areas
              real<lower=0> sigma_cluster;       // between-cluster sd by habitat (should this grow with the mean?)
              real<lower=0> lsigma_point;        // between-point sd by habitat for large plots
              vector[n_cluster] sigma_c_raw;
              //vector<lower=0>[n_cluster] carbon_cluster;
              vector[n_area] sigma_area_raw;
            }
            
            transformed parameters{
              vector[n_area] alpha_area = alpha + sigma_area * sigma_area_raw;
              vector[n_cluster] cluster_mean = alpha_area[area_cluster] + predictor * beta;
              //vector[n_cluster] cluster_mean = alpha + predictor * beta;
              vector[n_cluster] carbon_cluster = cluster_mean + sigma_cluster * sigma_c_raw;
              vector[n_cluster] cluster_lmean = log(carbon_cluster) - (lsigma_point ^ 2)/2; // mean of the logarithm of point-level carbon stocks at the cluster scale
            }
            
            model{
              // priors
              alpha ~ normal(0, 30);
              sigma_area ~ normal(0, 20);
              sigma_area_raw ~ std_normal();
              beta ~ normal(0, 10);
              sigma_cluster ~ normal(0,20);
              sigma_c_raw ~ std_normal();
              lsigma_point ~ normal(0,5);
              // likelihood
              //alpha_area ~ normal(alpha, sigma_area);
              //carbon_cluster ~ normal(cluster_mean, sigma_cluster);
              lcarbon_point ~ normal(cluster_lmean[cluster], lsigma_point);
            }", file = "STAN\\carbdiv_clus_paramo.stan")












## Carbon diversity model - linear
sink("JAGS//CarbDiv_linear.txt")
cat("model {
    
    #Likelihood model
      for(i in 1:n_cluster){
        for(j in 1:niter){
          div[i,j] ~ dnorm(mu[i],sd)
        }
        
      mu[i] <- alpha + b1 * carb[i]
      div_pred[i] ~ dnorm(mu[i], sd)
      }
    
    #Priors
    alpha ~ dnorm(0,0.001)
    b1 ~ dnorm(0,0.001)
    sd ~ dunif(0,1)
    }
",fill=TRUE)
sink()

## Carbon diversity model - spline
sink("JAGS//CarbDiv_spline.txt")
cat("model {
    
    #Likelihood model
      for(i in 1:n_cluster){
        for(j in 1:niter){
          div[i,j] ~ dnorm(mu[i],tau)
        }
        
        mu[i] <- alpha + sum(carbfunc[i,])
        
        for(k in 1:3){
          carbspline[i,k] <- ifelse(carb[i] - theta[k] < 0, 0, carb[i] - theta[k])
          carbfunc[i,k] <- b[k] * carbspline[i,k] * kz[k]
        }
        
        div_pred[i] ~ dnorm(mu[i], tau)
      }
    
    for(k in 2:3){
      kz[k] <- ifelse(k > nz + 1, 0, 1)
      theta_add[k] ~ dnorm(add_mu[k],add_tau[k])
      theta[k] <- theta[k-1] + theta_add[k]
      b[k] ~ dnorm(0, 0.01)
      
      add_mu[k] ~ dnorm(0, 0.0005) T(0, )
      add_sd[k] ~ dgamma(1, 0.1)
      add_tau[k] <- 1/(add_sd[k]^2)
    }
    kz[1] <- 1
    theta[1] <- 0
    b[1] ~ dnorm(0, 0.001)
    
    nz ~ dbin(kprob, 2)
    kprob ~ dbeta(1,1)
    
    alpha ~ dnorm(0,0.001)
    sd ~ dgamma(1,0.1)
    tau <- 1/(sd^2)
    
    # for(l in 1:seqlength){
    #   for(k in 1:4){
    #     predspline[l,k] <- ifelse(carbseq[l] - theta[k] < 0, 0, carbseq[l] - theta[k])
    #     predfunc[l,k] <- b[k] * predspline[l,k] * kz[k]
    #   }
    #   predmu[l] <- alpha + sum(predfunc[l,])
    #   div_pred[l] ~ dnorm(predmu[l], sd)
    # }
  }
",fill=TRUE)
sink()

## Carbon diversity model - polynomial
sink("JAGS//CarbDiv_polynomial.txt")
cat("model {
    
    #Likelihood model
      for(i in 1:n_cluster){
        for(j in 1:niter){
          div[i,j] ~ dnorm(mu[i],sd)
        }
        
      mu[i] <- alpha + b1 * carb[i] + b2 * (carb[i] ^ 2)
      div_pred[i] ~ dnorm(mu[i], sd)
      }
    
    #Priors
    alpha ~ dnorm(0,0.001)
    b1 ~ dnorm(0, 0.001)
    b2 ~ dnorm(0, 0.001)
    sd ~ dunif(0,1)
    }
",fill=TRUE)
sink()


