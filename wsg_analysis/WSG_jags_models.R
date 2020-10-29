
if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}

setwd(dir.path)

## Individual tree model
sink("JAGS//WSG_trees.txt")
cat("model {
    
    #Likelihood model
    for(n in 1:Nobs){
      WSG[n] ~ dnorm(mu[n], tau1 * speciesid[n] + tau0 * abs(1-speciesid[n]))
      mu[n] <- alpha + b_vol * vol[n] + b_spec[Species[n]] * speciesid[n] + b_site[SiteNum[n]]
      
      res[n] <- WSG[n] - mu[n]
      WSG.pred[n] ~ dnorm(mu[n], tau1 * speciesid[n] + tau0 * abs(1-speciesid[n]))
      res.pred[n] <- WSG.pred[n] - mu[n]
      
      loglik[n] <- logdensity.norm(WSG[n], mu[n], tau1 * speciesid[n] + tau0 * abs(1-speciesid[n]))
    }
    
    #Priors
    tau0 ~ dgamma(1,0.001)
    tau1 ~ dgamma(1,0.001)
    alpha ~ dnorm(0.5,1)
    b_vol ~ dnorm(0,2)
    
    for(i in 1:nspec){
      b_spec[i] ~ dnorm(muspec, tauspec)
    }
    muspec ~ dnorm(0,2)
    tauspec ~ dgamma(2,0.001)
    
    for(i in 1:narea){
      b_area[i] ~ dnorm(0,tauarea)
    }
    for(j in 1:ncluster){
      b_cluster[j] ~ dnorm(b_area[AreaNumC[j]],taucluster)
    }
    for(k in 1:nsites){
      b_site[k] ~ dnorm(b_cluster[ClusNumS[k]],tausite)
    }
    
    tauarea ~ dgamma(1,0.001)
    taucluster ~ dgamma(1,0.001)
    tausite ~ dgamma(1,0.001)
    
    # Derived parameters
    var0 <- 1/tau0
    var1 <- 1/tau1
    varspec <- 1/tauspec
    vararea <- 1/tauarea
    varcluster <- 1/taucluster
    varsite <- 1/tausite
    }
",fill=TRUE)
sink()

#####################
## Plot WSG models ##
#####################

# # Non-spatial models
# sink("JAGS//WSG_plots_null.txt")
# cat("model {
# 
#     #Likelihood model
#     for(n in 1:Nobs){
#       WSG[n] ~ dnorm(mu[n], tau)
#       mu[n] <- alpha
# 
#      res[n] <- WSG[n] - mu[n]
#      WSG.pred[n] ~ dnorm(mu[n], tau)
#      res.pred[n] <- WSG.pred[n] - mu[n]
# 
#      loglik[n] <- logdensity.norm(WSG[n], mu[n], tau)
#    }
# 
#    #Priors
#    var0 <- 1/tau
#    tau ~ dgamma(1,0.001)
#    alpha ~ dnorm(0.5,2)T(0,)
#    }
# ",fill=TRUE)
# sink()

sink("JAGS//WSG_plots_null.txt")
cat("model {
    
    #Likelihood model
    for(n in 1:Nobs){
      mu[n] <- alpha
      
      for(i in 1:ndraws){
        WSG[n,i] ~ dnorm(mu[n], tau)
        
        res_i[n,i] <- WSG[n,i] - mu[n]
        loglik_i[n,i] <- logdensity.norm(WSG[n,i], mu[n], tau)
      }
      
      MSE[n] <- sum(pow(res_i[n,],2)) / ndraws
      loglik[n] <- sum(loglik_i[n,])
      
      WSG.pred[n] ~ dnorm(mu[n], tau)
    }
    
    #Priors
    tau ~ dgamma(1,0.001)
    var0 <- 1/tau
    alpha ~ dnorm(0.5,2)T(0,)
    }
",fill=TRUE)
sink()

sink("JAGS//WSG_plots_cov.txt")
cat("model {
    
    #Likelihood model
    for(n in 1:Nobs){
      mu[n] <- alpha + inprod(b_cov[],cov[n,])
      
      for(i in 1:ndraws){
        WSG[n,i] ~ dnorm(mu[n], tau)
        
        res_i[n,i] <- WSG[n,i] - mu[n]
        loglik_i[n,i] <- logdensity.norm(WSG[n,i], mu[n], tau)
      }
      
      MSE[n] <- sum(pow(res_i[n,],2)) / ndraws
      loglik[n] <- sum(loglik_i[n,])
      
      WSG.pred[n] ~ dnorm(mu[n], tau)
    }
    
    #Priors
    var0 <- 1/tau
    tau ~ dgamma(1,0.001)
    alpha ~ dnorm(0.5,2)T(0,)
    
    for(i in 1:ncov){
      b_cov[i] ~ dnorm(0,2)
    }
    
    varspatial <- 0
  }
",fill=TRUE)
sink()


# Spatial models

sink("JAGS//WSG_plots_null_spatial.txt")
cat("model {
    
    #Likelihood model
    for(n in 1:Nobs){
      mu[n] <- alpha + b_spatial[SpatialNum[n]]
      
      for(i in 1:ndraws){
        WSG[n,i] ~ dnorm(mu[n], tau)
        
        res_i[n,i] <- WSG[n,i] - mu[n]
        loglik_i[n,i] <- logdensity.norm(WSG[n,i], mu[n], tau)
      }
      
      MSE[n] <- sum(pow(res_i[n,],2)) / ndraws
      loglik[n] <- sum(loglik_i[n,])
      
      WSG.pred[n] ~ dnorm(mu[n], tau)
    }
    
    #Priors
    var0 <- 1/tau
    tau ~ dgamma(1,0.001)
    alpha ~ dnorm(0.5,2)T(0,)
    
    for(i in 1:nspatial){
      b_spatial[i] ~ dnorm(0,tauspatial)
    }
    varspatial <- 1/tauspatial
    tauspatial ~ dgamma(1,0.001)
  }
",fill=TRUE)
sink()

sink("JAGS//WSG_plots_spatial_cov.txt")
cat("model {

    #Likelihood model
    for(n in 1:Nobs){
      mu[n] <- alpha + inprod(b_cov[],cov[n,]) + b_spatial[SpatialNum[n]]
      
      for(i in 1:ndraws){
        WSG[n,i] ~ dnorm(mu[n], tau)
        
        res_i[n,i] <- WSG[n,i] - mu[n]
        loglik_i[n,i] <- logdensity.norm(WSG[n,i], mu[n], tau)
      }
      
      MSE[n] <- sum(pow(res_i[n,],2)) / ndraws
      loglik[n] <- sum(loglik_i[n,])
      
      WSG.pred[n] ~ dnorm(mu[n], tau)
    }
    
    #Priors
    sd0 ~ dt(0, 0.75, 1)T(0,)
    var0 <- pow(sd0, 2)
    tau <- 1/var0
    
    #var0 <- 1/tau
    #tau ~ dgamma(1,0.001)
    alpha ~ dnorm(0.5,2)T(0,)
    
    for(i in 1:ncov){
      b_cov[i] ~ dnorm(0,2)
    }
    
    for(i in 1:nspatial){
      b_spatial[i] ~ dnorm(0,tauspatial)
    }
    sdspatial ~ dt(0, 0.75, 1)T(0,)
    varspatial <- pow(sdspatial, 2)
    tauspatial <- 1/varspatial
    
    #varspatial <- 1/tauspatial
    #tauspatial ~ dgamma(1,0.001)
  }
",fill=TRUE)
sink()

## Plot WSG model without spatial covariates for cross-validation
sink("JAGS//WSG_plots_CV.txt")
cat("model {
    
    #Likelihood model
    for(n in 1:Nobst){
      mut[n] <- alpha + inprod(b_cov[],covt[n,])
      
      for(i in 1:ndraws){
        WSGt[n,i] ~ dnorm(mut[n], tau)
      }
    }
    
    for(n in 1:Nobsp){
      mup[n] <- alpha + inprod(b_cov[],covp[n,])
      WSG.pred[n] ~ dnorm(mup[n], tau)
    
      for(i in 1:ndraws){
        res_i[n,i] <- WSGp[n,i] - mup[n]
        loglik_i[n,i] <- logdensity.norm(WSGp[n,i], mup[n], tau)
      }
      
      MSE[n] <- sum(pow(res_i[n,],2)) / ndraws
      elpd[n] <- sum(loglik_i[n,])
    }
    
    #Priors
    var0 <- 1/tau
    tau ~ dgamma(1,0.001)
    alpha ~ dnorm(0.5,2)T(0,)
    
    for(i in 1:ncov){
      b_cov[i] ~ dnorm(0,2)
    }
  }
",fill=TRUE)
sink()

## Plot WSG model without spatial covariates for cross-validation
sink("JAGS//WSG_plots_null_CV.txt")
cat("model {
    
    #Likelihood model
    for(n in 1:Nobst){
      mut[n] <- alpha + inprod(b_cov[],covt[n,])
      
      for(i in 1:ndraws){
        WSGt[n,i] ~ dnorm(mut[n], tau)
      }
    }
    
    for(n in 1:Nobsp){
      mup[n] <- alpha + inprod(b_cov[],covp[n,])
      WSG.pred[n] ~ dnorm(mup[n], tau)
      
      for(i in 1:ndraws){
        res_i[n,i] <- WSGp[n,i] - mup[n]
        loglik_i[n,i] <- logdensity.norm(WSGp[n,i], mup[n], tau)
      }
      
      MSE[n] <- sum(pow(res_i[n,],2)) / ndraws
      elpd[n] <- sum(loglik_i[n,])
    }
    
    #Priors
    var0 <- 1/tau
    tau ~ dgamma(1,0.001)
    alpha ~ dnorm(0.5,2)T(0,)
    
    for(i in 1:ncov){
      b_cov[i] ~ dnorm(0,2)
    }
  }
",fill=TRUE)
sink()

