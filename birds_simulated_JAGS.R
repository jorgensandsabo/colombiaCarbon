# Use machine identifier to automatically set directory with points file
jorgen.desktop <- file.exists('C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD\\machine_identifier_lu847jp1o.txt')
if(jorgen.desktop){
  dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
setwd(dir.path)


# i indexes points, j indexes visits, k indexes species.

# Define data size
#n_cluster <- 300
n_cluster <- 30
ppc <- 3
n_point <- n_cluster*ppc
n_species <- 30
n_visit <- rep(4, n_point)

# Define covariates
clusterID <- rep(c(1:n_cluster), ppc)
sp_cov1 <- runif(n_species) - .5
sp_cov2 <- runif(n_species) - .5
pt_cov1 <- runif(n_point) - .5
cl.cov.1 <- runif(n_cluster) - .5
pt_cov2 <- rep(NA, n_point)
for(i in 1:n_point){
  pt_cov2[i] <- cl.cov.1[clusterID[i]]
}

vis_cov1 <- matrix(data = NA, nrow = n_point, ncol = max(n_visit))
for(i in 1:n_point){
  vis_cov1[i, ] <- runif(n_visit[i]) - .5
}

exists_visit <- matrix(data = 0, nrow = n_point, ncol = max(n_visit))
for(i in 1:n_point){
  for(j in 1:n_visit[i]){
    exists_visit[i,j] <- 1
  }
}

# Define hyperparameters
occ.hyper <- list(b0 = c(0, .5), b1 = c(0, 1), b2 = c(1, .5), b3 = c(-1, 1), b4 = c(0, .2),
                  b5 = c(1,2), b6 = c(2, 2))
b0 <- rnorm(n_species, occ.hyper$b0[1], occ.hyper$b0[2])
b1 <- matrix(data = rnorm(n_cluster*n_species, occ.hyper$b1[1], occ.hyper$b1[2]), nrow = n_cluster)
b2 <- rnorm(n_species, occ.hyper$b2[1], occ.hyper$b2[2])
b3 <- rnorm(n_species, occ.hyper$b3[1], occ.hyper$b3[2])
b4 <- rnorm(n_species, occ.hyper$b4[1], occ.hyper$b4[2])
b5 <- rnorm(n_species, occ.hyper$b5[1], occ.hyper$b5[2])
b6 <- rnorm(n_species, occ.hyper$b6[1], occ.hyper$b6[2])

det.hyper <- list(d0 = c(-2, .5), d1 = c(0, 1), d2 = c(0, .5), d3 = c(1, 1), d4 = c(0, .2),
                  d5 = c(2, .5))
d0 <- rnorm(n_species, det.hyper$d0[1], det.hyper$d0[2])
d1 <- matrix(data = rnorm(n_point*n_species, det.hyper$d1[1], det.hyper$d1[2]), nrow = n_point)
d2 <- rnorm(n_species, det.hyper$d2[1], det.hyper$d2[2])
d3 <- rnorm(n_species, det.hyper$d3[1], det.hyper$d3[2])
d4 <- rnorm(n_species, det.hyper$d4[1], det.hyper$d4[2])
d5 <- rnorm(n_species, det.hyper$d5[1], det.hyper$d5[2])

# Simulate parameters from hyperparameters
logit.occ <- psi <- matrix(NA, nrow = n_point, ncol = n_species)
for(i in 1:n_point){
  for(k in 1:n_species){
    logit.occ[i, k] <- b0[k] + b1[clusterID[i],k] + b2[k]*sp_cov1[k] + b3[k]*sp_cov2[k] + 
      b4[k]*pt_cov1[i] + b5[k]*pt_cov2[i] + b6[k]*sp_cov1[k]*pt_cov2[i]
    psi[i, k] <- boot::inv.logit(logit.occ[i, k])
  }
}

logit.det <- theta <- array(NA, dim = c(n_point, max(n_visit), n_species))
for(i in 1:n_point){
  for(j in 1:n_visit[i]){
    for(k in 1:n_species){
      logit.det[i, j, k] <- d0[k] + d1[i,k] + d2[k]*sp_cov1[k] + d3[k]*sp_cov2[k] + 
        d4[k]*pt_cov1[i] + d5[k]*vis_cov1[i,j]
      theta[i, j, k] <- boot::inv.logit(logit.det[i, j, k])
    }
  }
}


# simulate data from parameters
Z <- matrix(NA, nrow = n_point, ncol = n_species)
for(i in 1:n_point){
  for(k in 1:n_species){
    Z[i, k] <- rbinom(1, 1, psi[i, k])
  }
}

det_data <- array(NA, dim = c(n_point, max(n_visit), n_species))
for(i in 1:n_point){
  for(j in 1:n_visit[i]){
    for(k in 1:n_species){
      det_data[i,j,k] <- Z[i, k] * rbinom(1, 1, theta[i,j,k])
    }
  }
}


# Detection/non-detection site-level
Q <- apply(det_data, c(1,3), function(x){return(as.numeric(sum(x) > 0))})


###################################
# JAGS MODELS
###################################
priors <-  " 
  #Priors
  d2 ~ dnorm(0,0.1)
  d3 ~ dnorm(0,0.1)
  
  b2 ~ dnorm(0,0.1)
  b3 ~ dnorm(0,0.1)
  
  for(k in 1:nspecies){
    d0[k] ~ dnorm(d0mu,d0tau)
    d4[k] ~ dnorm(d4mu,d4tau)
    d5[k] ~ dnorm(d5mu,d5tau)
    
    b0[k] ~ dnorm(b0mu,b0tau)
    b4[k] ~ dnorm(b4mu,b4tau)
    b5[k] ~ dnorm(b5mu,b5tau)
    b6[k] ~ dnorm(b6mu,b6tau)
      
    for(l in 1:ncluster){
      d1[k,l] ~ dnorm(d1mu,d1tau)
      b1[k,l] ~ dnorm(b1mu,b1tau)        
    }
  }
  # Hyper-priors
  d0mu ~ dnorm(0,0.1)
  d0sd ~ dunif(0,100)
  d0tau <- pow(d0sd,-2)
  d1mu ~ dnorm(0,0.1)
  d1sd ~ dunif(0,100)
  d1tau <- pow(d1sd,-2)
  d4mu ~ dnorm(0,0.1)
  d4sd ~ dunif(0,100)
  d4tau <- pow(d4sd,-2)
  d5mu ~ dnorm(0,0.1)
  d5sd ~ dunif(0,100)
  d5tau <- pow(d5sd,-2)
  
  b0mu ~ dnorm(0,0.1)
  b0sd ~ dnorm(0,0.1) T(0,)
  b0tau <- pow(b0sd,-2)
  b1mu ~ dnorm(0,0.1)
  b1sd ~ dnorm(0,0.1) T(0,)
  b1tau <- pow(b1sd,-2)
  b4mu ~ dnorm(0,0.1)
  b4sd ~ dnorm(0,0.1) T(0,)
  b4tau <- pow(b4sd,-2)
  b5mu ~ dnorm(0,0.1)
  b5sd ~ dnorm(0,0.1) T(0,)
  b5tau <- pow(b5sd,-2)
  b6mu ~ dnorm(0,0.1)
  b6sd ~ dnorm(0,0.1) T(0,)
  b6tau <- pow(b6sd,-2)
}"

# Non-marginalized model
model1.string <- paste("
  model {
    for(k in 1:nspecies){
      for(i in 1:nsites){
        for(j in 1:nvisits[i]){
          
          count[i,j,k] ~ dbern(p.occ[k,i,j])
          p.occ[k,i,j] <- z[i,k]*p[k,i,j]
          logit(p[k,i,j]) <- d0[k] + d2 * sp_cov1[k] + d3 * sp_cov2[k] + d4[k] * pt_cov1[i] + d5[k] * vis_cov1[i,j]
        }
    
        z[i,k] ~ dbern(psi[k,i])
        logit(psi[k,i]) <- b0[k] + b1[k,clusterID[i]] + b2 * sp_cov1[k] + b3 * sp_cov2[k] + b4[k] * pt_cov1[i] + b5[k] * pt_cov2[i] + b6[k] * sp_cov1[k] * pt_cov2[i]
      }
    }
    ",priors,sep="")

sink("JAGS\\mod1.txt")
cat(paste(model1.string))
sink()

# Marginalized model if/else in JAGS
model1.string <- paste("
  model {
  
    # Process model
    for(i in 1:nsites){
      for(k in 1:nspecies){

        ones[k,i] ~ dbern(ifZ[k, i, Qadd[i,k]])
        # then
        ifZ[k, i, 2] <- psi[k,i]
        # else
        ifZ[k, i, 1] <- 1 - psi[k,i] * prod(pcomp[k,i,1:nvisits[i]])

        for(j in 1:nvisits[i]){
          pcomp[k,i,j] <- 1-p[k,i,j]
          
          count[i,j,k] ~ dbern(ifCOUNT[k, i, j, Qadd[i,k]])
          
          # then
          ifCOUNT[k,i,j,2] <- p[k,i,j]
          # else
          ifCOUNT[k,i,j,1] <- 0
          
          # outside ifelse:
          logit(p[k,i,j]) <- d0[k] + d2 * sp_cov1[k] + d3 * sp_cov2[k] + d4[k] * pt_cov1[i] + d5[k] * vis_cov1[i,j]
        }

        # outside ifelse:   
        logit(psi[k,i]) <- b0[k] + b1[k,clusterID[i]] + b2 * sp_cov1[k] + b3 * sp_cov2[k] + b4[k] * pt_cov1[i] + b5[k] * pt_cov2[i] + b6[k] * sp_cov1[k] * pt_cov2[i]
      }
    }
",priors,sep="")

sink("JAGS\\mod1marg.txt")
cat(paste(model1.string))
sink()

# Marginalized model if/else in R
model1.string <- paste("
  model {
  
    # Process model
    for(i in 1:nsites){
      for(k in 1:nspecies){

        ones[k,i] ~ dbern(psi[k,i] * Q[i,k] + (1 - psi[k,i] * prod(pcomp[k,i,1:nvisits[i]]))*Qinv[i,k])

        for(j in 1:nvisits[i]){
          pcomp[k,i,j] <- 1-p[k,i,j]
          count[i,j,k] ~ dbern(p[k,i,j]*Q[i,k])
          logit(p[k,i,j]) <- d0[k] + d2 * sp_cov1[k] + d3 * sp_cov2[k] + d4[k] * pt_cov1[i] + d5[k] * vis_cov1[i,j]
        }

        # outside ifelse:   
        logit(psi[k,i]) <- b0[k] + b1[k,clusterID[i]] + b2 * sp_cov1[k] + b3 * sp_cov2[k] + b4[k] * pt_cov1[i] + b5[k] * pt_cov2[i] + b6[k] * sp_cov1[k] * pt_cov2[i]
      }
    }
",priors,sep="")

sink("JAGS\\mod2marg.txt")
cat(paste(model1.string))
sink()

# Marginalized model complex data structure
det_data_split <- det_data
N1 <- colSums(Q)
N <- apply(det_data_split, c(1,3), function(x){sum(x)})
pt_cov1_split <- matrix(rep(pt_cov1,n_species),nrow=n_point,ncol=n_species)
pt_cov2_split <- matrix(rep(pt_cov2,n_species),nrow=n_point,ncol=n_species)
vis_cov1_split <- array(rep(vis_cov1,n_species),dim=c(n_point,max(n_visit),n_species))
nvisit_split <- matrix(rep(n_visit,n_species),nrow=n_point,ncol=n_species)
clusterID_split <- matrix(rep(clusterID,n_species),nrow=n_point,ncol=n_species)
for(i in 1:n_species){
  neworder <- rev(order(N[,i]))
  pt_cov1_split[,i] <- pt_cov1_split[neworder,i]
  pt_cov2_split[,i] <- pt_cov2_split[neworder,i]
  nvisit_split[,i] <- nvisit_split[neworder,i]
  clusterID_split[,i] <- clusterID_split[neworder,i]
  for(j in 1:4){
    vis_cov1_split[,j,i] <- vis_cov1_split[neworder,j,i]
    det_data_split[,j,i] <- det_data_split[neworder,j,i]
  }
}

model1.string <- paste("
  model {
  
    # Process model
    for(k in 1:nspecies){
      for(i in 1:N1[k]){
      
        ones[k,i] ~ dbern(psi[k,i])
        logit(psi[k,i]) <- b0[k] + b1[k,clusterID[i,k]] + b2 * sp_cov1[k] + b3 * sp_cov2[k] + b4[k] * pt_cov1[i,k] + b5[k] * pt_cov2[i,k] + b6[k] * sp_cov1[k] * pt_cov2[i,k]
        
        for(j in 1:nvisits[i,k]){
        
          count[i,j,k] ~ dbern(p[k,i,j])
          logit(p[k,i,j]) <- d0[k] + d2 * sp_cov1[k] + d3 * sp_cov2[k] + d4[k] * pt_cov1[i,k] + d5[k] * vis_cov1[i,j,k]
        }
      }
      for(i in (N1[k]+1):nsites){
      
        ones[k,i] ~ dbern(1 - psi[k,i] * prod(pcomp[k,i,1:nvisits[i,k]]))
        logit(psi[k,i]) <- b0[k] + b1[k,clusterID[i,k]] + b2 * sp_cov1[k] + b3 * sp_cov2[k] + b4[k] * pt_cov1[i,k] + b5[k] * pt_cov2[i,k] + b6[k] * sp_cov1[k] * pt_cov2[i,k]
        
        for(j in 1:nvisits[i,k]){
        
          pcomp[k,i,j] <- 1-p[k,i,j]
          logit(p[k,i,j]) <- d0[k] + d2 * sp_cov1[k] + d3 * sp_cov2[k] + d4[k] * pt_cov1[i,k] + d5[k] * vis_cov1[i,j,k]
        }
      }
    }
",priors,sep="")

sink("JAGS\\mod3marg.txt")
cat(paste(model1.string))
sink()

###################################
# STAN MODEL
###################################
stan.model <- '
data {
  int<lower=1> n_point; //number of sites
  int<lower=1> n_visit; //fixed number of visits
  int<lower=1> n_species; //number of species
  int<lower=1> n_cluster; //number of clusters
  int<lower=0, upper=1> det_data[n_point, n_visit, n_species]; //detection history
  int<lower=0, upper=1> Q[n_point, n_species]; //at least one detection
  int<lower=0, upper=n_cluster> clusterID[n_point]; //cluster identifier (for random effects)
  vector[n_species] sp_cov1; //species covariate 1
  vector[n_species] sp_cov2; //species covariate 2
  vector[n_point] pt_cov1; //point covariate 1
  vector[n_point] pt_cov2; //point covariate 2
  matrix[n_point, n_visit] vis_cov1; //visit covariate 1
}
parameters {
  real mu_b0;
  real<lower=0> sigma_b0;
  vector[n_species] b0_raw;
  
  real<lower=0> sigma_b1;
  matrix[n_species, n_cluster] b1_raw;
  
  real b2;
  
  real b3;
  
  real mu_b4;
  real<lower=0> sigma_b4;
  vector[n_species] b4_raw;
  
  real mu_b5;
  real<lower=0> sigma_b5;
  vector[n_species] b5_raw;
  
  real mu_b6;
  real<lower=0> sigma_b6;
  vector[n_species] b6_raw;
  
  real mu_d0;
  real<lower=0> sigma_d0;
  vector[n_species] d0_raw;
  
  real<lower=0> sigma_d1;
  matrix[n_species, n_point] d1_raw;
  
  real d2;
  real d3;
  real mu_d4;
  real<lower=0> sigma_d4;
  vector[n_species] d4_raw;
  
  real mu_d5;
  real<lower=0> sigma_d5;
  vector[n_species] d5_raw;
}
transformed parameters{
  vector[n_species] b0 = mu_b0 + b0_raw * sigma_b0;
  matrix[n_species, n_cluster] b1 = b1_raw * sigma_b1;
  vector[n_species] b4 = mu_b4 + b4_raw * sigma_b4;
  vector[n_species] b5 = mu_b5 + b5_raw * sigma_b5;
  vector[n_species] b6 = mu_b6 + b6_raw * sigma_b6;
  
  vector[n_species] d0 = mu_d0 + d0_raw * sigma_d0;
  matrix[n_species, n_point] d1 = d1_raw * sigma_d1;
  vector[n_species] d4 = mu_d4 + d4_raw * sigma_d4;
  vector[n_species] d5 = mu_d5 + d5_raw * sigma_d5;
  real logit_psi[n_point, n_species];
  real logit_theta[n_point, n_visit, n_species];
  matrix[n_point, n_species] log_prob_increment;
  for(i in 1:n_point){
    for(k in 1:n_species){
      logit_psi[i,k] = b0[k] + b1[k, clusterID[i]] + b2*sp_cov1[k] + b3*sp_cov2[k] + b4[k]*pt_cov1[i] + b5[k]*pt_cov2[i] + b6[k]*sp_cov1[k]*pt_cov2[i];
    }
  }
  for(i in 1:n_point){
    for(j in 1:n_visit){
      for(k in 1:n_species){
        logit_theta[i,j,k] = d0[k] + d2*sp_cov1[k] + d3*sp_cov2[k] + d4[k]*pt_cov1[i] + d5[k]*vis_cov1[i,j];
      }
    }
  }
  
  for(i in 1:n_point){
    for(k in 1:n_species){
      if(Q[i,k] == 1)
        log_prob_increment[i,k] = log_inv_logit(logit_psi[i,k]) + 
                                    bernoulli_logit_lpmf(det_data[i,1,k] | logit_theta[i,1,k]) + 
                                    bernoulli_logit_lpmf(det_data[i,2,k] | logit_theta[i,2,k]) +
                                    bernoulli_logit_lpmf(det_data[i,3,k] | logit_theta[i,3,k]) +
                                    bernoulli_logit_lpmf(det_data[i,4,k] | logit_theta[i,4,k]);
      else
        log_prob_increment[i,k] = log_sum_exp(log_inv_logit(logit_psi[i,k]) + log1m_inv_logit(logit_theta[i,1,k]) + 
                                                    log1m_inv_logit(logit_theta[i,2,k]) + log1m_inv_logit(logit_theta[i,3,k]) + 
                                                    log1m_inv_logit(logit_theta[i,4,k]), 
                                                    log1m_inv_logit(logit_psi[i,k]));
    }
  }
}
 
model {
  //Hyper-priors:
  mu_b0 ~ normal(0,10);
  b2 ~ normal(0,10);
  b3 ~ normal(0,10);
  mu_b4 ~ normal(0,10);
  mu_b5 ~ normal(0,10);
  mu_b6 ~ normal(0,10);
  
  mu_d0 ~ normal(0,10);
  d2 ~ normal(0,10);
  d3 ~ normal(0,10);
  mu_d4 ~ normal(0,10);
  mu_d5 ~ normal(0,10);
  
  sigma_b0 ~ normal(0,10);
  sigma_b1 ~ normal(0,10);
  sigma_b4 ~ normal(0,10);
  sigma_b5 ~ normal(0,10);
  sigma_b6 ~ normal(0,10);
  
  sigma_d0 ~ normal(0,10);
  sigma_d1 ~ normal(0,10);
  sigma_d4 ~ normal(0,10);
  sigma_d5 ~ normal(0,10);
  
  //Random Effects
  b0_raw ~ normal(0, 1);
  to_vector(b1_raw) ~ normal(0, 1);
  b4_raw ~ normal(0, 1);
  b5_raw ~ normal(0, 1);
  b6_raw ~ normal(0, 1);
  
  d0_raw ~ normal(0, 1);
  to_vector(d1_raw) ~ normal(0, 1);
  d4_raw ~ normal(0, 1);
  d5_raw ~ normal(0, 1);
  
  //Likelihood (data level)
  target += sum(log_prob_increment);
}'

###################################
# RUN JAGS MODELS
###################################
nc <- 4
n.adapt <- 5000
n.burn <- 5000
n.iter <- 1000
thin <- 1

moddata <- list(count = det_data, ones = matrix(1,n_species,n_point),
                Q = Q, Qadd = Q + 1, Qinv = ifelse(Q > 0, 0, 1), clusterID = clusterID,
                nspecies = n_species, nsites = n_point, ncluster = n_cluster, nvisits = n_visit,
                sp_cov1 = sp_cov1, sp_cov2 = sp_cov2, 
                pt_cov1 = pt_cov1, pt_cov2 = pt_cov2,
                vis_cov1 = vis_cov1)
mod1params <- c("psi","p","z",
                "d0","d1","d2","d3","d4","d5",
                "b0","b1","b2","b3","b4","b5","b6",
                "d0mu","d0sd","d1mu","d1sd","d4mu","d4sd","d5mu","d5sd",
                "b0mu","b0sd","b1mu","b1sd","b4mu","b4sd","b5mu","b5sd","b6mu","b6sd")
modmargparams <- c("psi","p",
                   "d0","d1","d2","d3","d4","d5",
                   "b0","b1","b2","b3","b4","b5","b6",
                   "d0mu","d0sd","d1mu","d1sd","d4mu","d4sd","d5mu","d5sd",
                   "b0mu","b0sd","b1mu","b1sd","b4mu","b4sd","b5mu","b5sd","b6mu","b6sd")
mod1inits <- list(z = Q)

moddata2 <- list(count = det_data, ones = matrix(1,n_species,n_point),
                N1 = N1, clusterID = clusterID_split,
                nspecies = n_species, nsites = n_point, ncluster = n_cluster, nvisits = nvisit_split,
                sp_cov1 = sp_cov1, sp_cov2 = sp_cov2, 
                pt_cov1 = pt_cov1_split, pt_cov2 = pt_cov2_split,
                vis_cov1 = vis_cov1_split)

# Non-marginalized model
mod.start.time = proc.time()
cl <- parallel::makePSOCKcluster(nc)
tmp <- parallel::clusterEvalQ(cl, library(dclone))
dclone::parLoadModule(cl, "glm")
dclone::parListModules(cl)

mod.samples <- dclone::jags.parfit(cl, moddata, 
                                   mod1params, "JAGS\\mod1.txt", 
                                   inits=mod1inits, 
                                   n.chains=nc, n.adapt=n.adapt, n.update = n.burn, thin = thin, n.iter = n.iter)

parallel::stopCluster(cl)
mod.end.time = proc.time()
mod.dtime = mod.end.time - mod.start.time

summary.mod <- MCMCvis::MCMCsummary(mod.samples)

# Marginalized model if/else in JAGS
marg1.start.time = proc.time()
cl <- parallel::makePSOCKcluster(nc)
tmp <- parallel::clusterEvalQ(cl, library(dclone))
dclone::parLoadModule(cl, "glm")
dclone::parListModules(cl)

marg1.samples <- dclone::jags.parfit(cl, moddata, 
                                      modmargparams, "JAGS\\mod1marg.txt", 
                                      #inits=inits, 
                                      n.chains=nc, n.adapt=n.adapt, n.update = n.burn, thin = thin, n.iter = n.iter)

parallel::stopCluster(cl)
marg1.end.time = proc.time()
marg1.dtime = marg1.end.time - marg1.start.time

summary.marg1 <- MCMCvis::MCMCsummary(marg1.samples)

# Marginalized model if/else in R
marg2.start.time = proc.time()
cl <- parallel::makePSOCKcluster(nc)
tmp <- parallel::clusterEvalQ(cl, library(dclone))
dclone::parLoadModule(cl, "glm")
dclone::parListModules(cl)

marg2.samples <- dclone::jags.parfit(cl, moddata, 
                                     modmargparams, "JAGS\\mod2marg.txt", 
                                     #inits=inits, 
                                     n.chains=nc, n.adapt=n.adapt, n.update = n.burn, thin = thin, n.iter = n.iter)

parallel::stopCluster(cl)
marg2.end.time = proc.time()
marg2.dtime = marg2.end.time - marg2.start.time

summary.marg2 <- MCMCvis::MCMCsummary(marg2.samples)

# Marginalized model complex data structure
marg3.start.time = proc.time()
cl <- parallel::makePSOCKcluster(nc)
tmp <- parallel::clusterEvalQ(cl, library(dclone))
dclone::parLoadModule(cl, "glm")
dclone::parListModules(cl)

marg3.samples <- dclone::jags.parfit(cl, moddata2, 
                                     modmargparams, "JAGS\\mod3marg.txt", 
                                     #inits=inits, 
                                     n.chains=nc, n.adapt=n.adapt, n.update = n.burn, thin = thin, n.iter = n.iter)

parallel::stopCluster(cl)
marg3.end.time = proc.time()
marg3.dtime = marg3.end.time - marg3.start.time

summary.marg3 <- MCMCvis::MCMCsummary(marg3.samples)

###################################
# RUN STAN MODEL
###################################
library(rstan)
rstan_options(auto_write = T)

stan.data <- list(n_point = n_point, n_species = n_species, 
                  n_cluster = n_cluster, n_visit = 4,
                  det_data = det_data, 
                  Q = Q, 
                  clusterID = clusterID, 
                  sp_cov1 = sp_cov1, sp_cov2 = sp_cov2, 
                  pt_cov1 = pt_cov1, pt_cov2 = pt_cov2,
                  vis_cov1 = vis_cov1)

nc <- nc

stan.start.time <- Sys.time()
stan.samples <- stan(model_code = stan.model, data = stan.data, iter = 2000, chains = nc, cores = nc,
                     pars = c('logit_psi', 'logit_theta', 'log_prob_increment', 'd1_raw', 'd1', 'b1_raw', 'b1'),
                     include = FALSE)
stan.elapsed <- Sys.time() - stan.start.time

summary.samples <- summary(stan.samples)                                




# Extract benchmarks
mod.dtime
min(summary.mod[,7])
summary.mod[which(summary.mod[,7] == min(summary.mod[,7])),]
max(summary.mod[,6],na.rm=TRUE)
summary.mod[which(summary.mod[,6] == max(summary.mod[,6])),]
object.size(mod.samples)

marg1.dtime
min(summary.marg1[,7])
summary.marg1[which(summary.marg1[,7] == min(summary.marg1[,7])),]
max(summary.marg1[,6])
summary.marg1[which(summary.marg1[,6] == max(summary.marg1[,6])),]
object.size(marg1.samples)

marg2.dtime
min(summary.marg2[,7])
summary.marg2[which(summary.marg2[,7] == min(summary.marg2[,7])),]
max(summary.marg2[,6])
summary.marg2[which(summary.marg2[,6] == max(summary.marg2[,6])),]
object.size(marg2.samples)

marg3.dtime
min(summary.marg3[,7])
summary.marg3[which(summary.marg3[,7] == min(summary.marg3[,7])),]
max(summary.marg3[,6])
summary.marg3[which(summary.marg3[,6] == max(summary.marg3[,6])),]
object.size(marg3.samples)

stan.elapsed
which(summary.samples$summary[,9] == min(summary.samples$summary[,9]))  
min(summary.samples$summary[,9])                                        # min n_eff
max(summary.samples$summary[,10])                                       # max R-hat
object.size(stan.samples)

###### benchmarks
# JAGS non-marginalized
# 30 species, 30 clusters, many parameters saved
# 10.000 adaptations, 1000 burn-in, 1000 iterations, 4 chains
#    min neff = 11 (+ 0 for several Zs) = b1mu
#    max R-hat = 15.69 (+ NaN for several Zs)
#    elapsed = 3739.70
#    object.size: 588633976 bytes

# JAGS marginalized - if in JAGS
# 30 species, 30 clusters, many parameters saved
# 10.000 adaptations, 1000 burn-in, 1000 iterations, 4 chains
#    min neff = 46 = b0
#    max R-hat = 1.53 = b1
#    elapsed: 4694.44
#    object.size: 500907736

# JAGS marginalized - if in R
# 30 species, 30 clusters, many parameters saved
#    min neff = 
#    max R-hat = 
#    elapsed = 
#    object.size: 

# JAGS marginalized - complex
# 30 species, 30 clusters, many parameters saved
# 10.000 adaptations, 1000 burn-in, 1000 iterations, 4 chains
#    min neff = 49 = b0
#    max R-hat = 3.14 = b1
#    elapsed = 4272.49
#    object.size: 500907736

# STAN
# 30 species, 30 clusters, many parameters saved
#    min neff (lp__) = 
#    max R-hat = 
#    grad. eval. range: 
#    elapsed range: 
#    object.size: 