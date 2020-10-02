
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
    sigma1 ~ dgamma(1,5)
    var1 <- pow(sigma1,2)
    tau1 <- 1/var1
    sigma0 ~ dgamma(1,5)
    var0 <- pow(sigma0,2)
    tau0 <- 1/var0
    alpha ~ dnorm(0.5,2)T(0,)
    b_vol ~ dnorm(0,2)
    
    for(i in 1:nspec){
      b_spec[i] ~ dnorm(muspec, tauspec)
    }
    muspec ~ dnorm(0,2)
    sdspec ~ dgamma(1,5)
    varspec <- pow(sdspec,2)
    tauspec <- 1/varspec
    
    for(i in 1:narea){
      b_area[i] ~ dnorm(0,tauarea)
    }
    for(j in 1:ncluster){
      b_cluster[j] ~ dnorm(b_area[AreaNumC[j]],taucluster)
    } 
    for(k in 1:nsites){
      b_site[k] ~ dnorm(b_cluster[ClusNumS[k]],tausite)
    }
    
    sdarea ~ dgamma(1,5)
    vararea <- pow(sdarea,2)
    tauarea <- 1/vararea
    sdcluster ~ dgamma(1,5)
    varcluster <- pow(sdcluster,2)
    taucluster <- 1/varcluster
    sdsite ~ dgamma(1,5)
    varsite <- pow(sdsite,2)
    tausite <- 1/varsite
    
    # Derived parameters
    fit <- sum(res[])
    fit.pred <- sum(res.pred[])
    }
",fill=TRUE)
sink()

## Individual tree model with variance structure on site intersite variation
sink("JAGS//WSG_trees_varstruct.txt")
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
    sigma1 ~ dgamma(1,5)
    var1 <- sigma1^2
    tau1 <- 1/var1
    sigma0 ~ dgamma(1,5)
    var0 <- sigma0^2
    tau0 <- 1/var0
    alpha ~ dnorm(0.5,2)T(0,)
    b_vol ~ dnorm(0,2)
    
    for(i in 1:nspec){
      b_spec[i] ~ dnorm(muspec, tauspec)
    }
    muspec ~ dnorm(0,2)
    sdspec ~ dgamma(1,5)
    varspec <- sdspec^2
    tauspec <- 1/varspec
    
    for(i in 1:narea){
      b_area[i] ~ dnorm(0, tauarea)
    }
    for(j in 1:ncluster){
      b_cluster[j] ~ dnorm(b_area[AreaNumC[j]], taucluster)
    } 
    for(k in 1:nsites){
      b_site[k] ~ dnorm(b_cluster[ClusNumS[k]], tausite[k])
      tausite[k] <- 1 / (varsite * exp(2 * dsite * sitecores[k]))
    }
    
    sdarea ~ dgamma(1,5)
    vararea <- pow(sdarea,2)
    tauarea <- 1/vararea
    sdcluster ~ dgamma(1,5)
    varcluster <- pow(sdcluster,2)
    taucluster <- 1/varcluster
    sdsite ~ dgamma(1,5)
    varsite <- pow(sdsite,2)
    dsite ~ dnorm(0,0.1)
    
    # Derived parameters
    fit <- sum(res[])
    fit.pred <- sum(res.pred[])
    }
",fill=TRUE)
sink()

## Plot WSG models
sink("JAGS//WSG_plots_null.txt")
cat("model {
    
    #Likelihood model
    for(n in 1:Nobs){
      WSG[n] ~ dnorm(mu[n], tau)
      mu[n] <- alpha + b_area[AreaNum[n]]

      res[n] <- WSG[n] - mu[n]
      WSG.pred[n] ~ dnorm(mu[n], tau)
      res.pred[n] <- WSG.pred[n] - mu[n]

      loglik[n] <- logdensity.norm(WSG[n], mu[n], tau)
    }
    
    #Priors
    sigma ~ dgamma(1,5)
    tau <- pow(sigma,-2)
    alpha ~ dnorm(0.5,2)T(0,)
    
    for(i in 1:narea){
      b_area[i] ~ dnorm(0,tauarea)
    }
    sdarea ~ dgamma(1,5)
    tauarea <- pow(sigma,-2)
    
    # Derived parameters
    fit <- sum(res[])
    fit.pred <- sum(res.pred[])
    }
",fill=TRUE)
sink()

sink("JAGS//WSG_plots_cov.txt")
cat("model {
    
    #Likelihood model
    for(n in 1:Nobs){
      WSG[n] ~ dnorm(mu[n], tau)
      mu[n] <- alpha + inprod(b_cov[],cov[n,]) + b_area[AreaNum[n]]
      
      res[n] <- WSG[n] - mu[n]
      WSG.pred[n] ~ dnorm(mu[n], tau)
      res.pred[n] <- WSG.pred[n] - mu[n]

      loglik[n] <- logdensity.norm(WSG[n], mu[n], tau)
    }
    
    #Priors
    sigma ~ dgamma(1,5)
    tau <- pow(sigma,-2)
    alpha ~ dnorm(0.5,2)T(0,)
    
    for(i in 1:ncov){
      b_cov[i] ~ dnorm(0,2)
    }
    
    for(i in 1:narea){
      b_area[i] ~ dnorm(0,tauarea)
    }
    sdarea ~ dgamma(1,5)
    tauarea <- pow(sigma,-2)
    
    # Derived parameters
    fit <- sum(res[])
    fit.pred <- sum(res.pred[])
    }
",fill=TRUE)
sink()

