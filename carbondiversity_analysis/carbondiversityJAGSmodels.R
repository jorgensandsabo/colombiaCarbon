

if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}

setwd(dir.path)

## Model to estimate individual tree WSG
sink("JAGS//CarbDiv_WSGmod.txt")
cat("model {
    
    #Likelihood model
    for(n in 1:Nobs){
      WSG[n] ~ dnorm(mu[n], tau1 * speciesid[n] + tau0 * abs(1-speciesid[n]))
      mu[n] <- alpha + b_spec[Species[n]] * speciesid[n] + b_site[SiteNum[n]] + b_habitat[habitat[n]]
      
      res[n] <- WSG[n] - mu[n]
      WSG.pred[n] ~ dnorm(mu[n], tau1 * speciesid[n] + tau0 * abs(1-speciesid[n]))
      res.pred[n] <- WSG.pred[n] - mu[n]
      
      loglik[n] <- logdensity.norm(WSG[n], mu[n], tau1 * speciesid[n] + tau0 * abs(1-speciesid[n]))
    }
    
    #Priors
    tau0 ~ dgamma(1,0.001)
    tau1 ~ dgamma(1,0.001)
    alpha ~ dnorm(0.5,300)
    
    for(i in 1:nspec){
      b_spec[i] ~ dnorm(muspec, tauspec)
    }
    muspec ~ dnorm(0,2)
    tauspec ~ dgamma(2,0.001)
    
    for(i in 1:nhabitat){
      b_habitat[i] ~ dnorm(0, 0.01)
    }
    
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

## Model to estimate grass biomass
sink("JAGS//CarbDiv_grassBM.txt")
cat("model {
    
    #Likelihood model
    for(n in 1:Nobs){
      Grass[n] ~ dlnorm(mu[n], tau[habitat[n]])
      mu[n] <- alpha[habitat[n]] + b_area[AreaNum[n]]
      
      res[n] <- Grass[n] - mu[n]
      Grass.pred[n] ~ dlnorm(mu[n], tau[habitat[n]])
      res.pred[n] <- Grass.pred[n] - mu[n]
      
      loglik[n] <- logdensity.norm(Grass[n], mu[n], tau[habitat[n]])
    }
    
    #Priors
    for(i in 1:nhabitat){
      alpha[i] ~ dnorm(0, 1)
      tau[i] ~ dgamma(1,6)
    }
    
    for(i in 1:narea){
      b_area[i] ~ dnorm(muarea, tauarea)
    }
    muarea ~ dnorm(0, 0.2)
    tauarea ~ dgamma(5,5)
    
    # Derived parameters
    var1 <- 1/tau
    vararea <- 1/tauarea
    }
",fill=TRUE)
sink()

### Carb - plot to cluster
sink("JAGS//Carbdiv_plot2clus.txt")
cat("model {
    
    #Likelihood model
    for(s in 1:nsite){
      carb[s] ~ dlnorm(mu.clus[clus[s]], tau.clus[clus[s]])
    }
    
    for(i in 1:nclus){
      #carb.clus[i] <- exp(mu.clus[i] + ((1/sqrt(tau.clus[i]))/2))
      carb.clus[i] <- exp(mu.clus[i])
      
      #tau.clus[i] ~ dunif(0.1,100)
      mu.clus[i] ~ dnorm(log(mediancarb[i]), 10)
      
      tau.clus[i] <- 1 / pow((sd.clus * exp(2 * d * size[i])), 2)
    }
    
    sd.clus ~ dunif(0,1)
    d ~ dnorm(0,0.0001)
  }
",fill=TRUE)
sink()

### Carb - plot to cluster - new
sink("JAGS//Carbdiv_plot2clus.txt")
cat("model {
    
    #Likelihood model
    for(s in 1:nsite){
      carb[s] ~ dlnorm(mu.clus[clus[s]], tau.clus[clus[s]])
    }
    
    for(i in 1:nclus){
      #carb.clus[i] <- exp(mu.clus[i] + (0.5 * sqrt(1/tau.clus[i])))
      #carb.clus[i] <- exp(mu.clus[i] + (0.5 * (1/tau.clus[i])))
      
      tau.clus[i] ~ dgamma(10,5)
      mu.clus[i] ~ dnorm(muc[i], sdclus)
      
      muc[i] <- alpha + bhab * habitat[i] + belev * elev[i] + barea[area[i]]
    }
    
    sdclus ~ dunif(0,10)
    alpha ~ dnorm(5.5, 1)
    belev ~ dnorm(0, 3)
    bhab ~ dnorm(0, 3)
    
    for(k in 1:narea){
      barea[k] ~ dnorm(0,3)
    }
    
    #alpha ~ dnorm(4, 1/(pow(0.6, 2)))
    #belev ~ dnorm(0, 1/(pow(1.5,2)))
  }
",fill=TRUE)
sink()


# ### Carb - plot to cluster - spatial decay
# sink("JAGS//Carbdiv_plot2clus.txt")
# cat("model {
#     
#     #Likelihood model
#     for(s in 1:nsite){
#       carb[s] ~ dlnorm(mu.clus[clus[s]], tau.clus[clus[s]])
#     }
#     
#     for(i in 1:nclus){
#       carb.clus[i] <- exp(mu.clus[i] + (0.5 * sqrt(1/tau.clus[i])))
#       
#       tau.clus[i] ~ dgamma(10,5)
#       mu.clus[i] <- alpha + bhab * habitat[i] + belev * elev[i] + W[i]
#       
#       muW[i] <- 0
#       
#       for(j in 1:nclus){
#         H[i,j] <- (1/tauw) * exp(-phi * pow(d[i,j], 2))
#       }
#     }
#     
#     W[1:nclus] ??? dmnorm(muW[], Omega[,])
#     tauw ??? dgamma(0.001,0.001)T(1,)
#     phi ??? dgamma(0.01,0.01)T(1,)
#     
#     Omega[1:nclus,1:nclus] <- inverse(H[1:nclus,1:nclus])
#     
#     sdclus ~ dunif(0,10)
#     alpha ~ dnorm(5.5, 1)
#     belev ~ dnorm(0, 3)
#     bhab ~ dnorm(0, 3)
#   }
# ",fill=TRUE)
# sink()

# ### Carb - plot to cluster - spatial decay
# sink("JAGS//test.txt")
# cat("model {
#     
#     #Likelihood model
#     for(i in 1:nclus){
#       tau.clus[i] ~ dgamma(10,5)
#       #mu.clus[i] <- alpha + W[i]
#     }
#     
#     #W[1:nclus] ??? dmnorm(muW[], Omega[,])
#     #tauw ??? dgamma(0.001,0.001)
#     #phi ??? dgamma(0.01,0.01)
#     phi ~ dgamma(1, 0.1)
#     tauw ~ dgamma(2,1)
#     W[1:nclus] ??? dmnorm.vcov(muW[], H[,])
#     
#     #alpha ~ dnorm(5.5, 1)
#     for(i in 1:nclus){
#      muW[i] <- 0
#      H[i,i] <- 1/tauw
#      for(j in 1:(i-1)){
#         H[i,j] <- (1/tauw) * exp(-phi * pow(d[i,j], 2))
#         H[j,i] <- H[i,j]
#       }}
#     #Omega[1:nclus,1:nclus] <- inverse(H[1:nclus, 1:nclus])
#     #Omega <- inverse(H)
#   }
# ",fill=TRUE)
# sink()
# 
# carbmod_clus <- jagsUI::jags(model.file = "JAGS//test.txt", data = carbdat, #inits = modinits,
#                              parameters.to.save = c("H","Omega","muW","tauw","phi","alpha","bhab","belev","barea","carb.clus","tau.clus","mu.clus","d","sd.clus"),
#                              n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
#                              parallel = T, n.cores = n.cores, codaOnly = c("loglik"))

### Carb - plot to area
sink("JAGS//Carbdiv_plot2area.txt")
cat("model {
    
    #Likelihood model
    for(s in 1:nsite){
      carb[s] ~ dlnorm(mu.area[area[s]], tau.area[area[s]])
    }
    
    for(i in 1:narea){
      #carb.area[i] <- exp(mu.area[i] + ((1/sqrt(tau.area[i]))/2))
      carb.area[i] <- exp(mu.area[i])
      
      #tau.area[i] ~ dunif(0.1,100)
      mu.area[i] ~ dnorm(log(mediancarb[i]), 10)
      
      tau.area[i] <- 1 / pow((sd.area * exp(2 * d * size[i])), 2)
    }
    
    sd.area ~ dunif(0,1)
    d ~ dnorm(0,0.0001)
  }
",fill=TRUE)
sink()

## Carbon diversity model - linear
sink("JAGS//CarbDiv_linear.txt")
cat("model {
    
    #Likelihood model
      for(i in 1:nsite){
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
      for(i in 1:nsite){
        for(j in 1:niter){
          div[i,j] ~ dnorm(mu[i],sd)
        }
        
        mu[i] <- alpha + sum(carbfunc[i,])
      
        for(k in 1:4){
          carbspline[i,k] <- ifelse(carb[i] - theta[k] < 0, 0, carb[i] - theta[k])
          carbfunc[i,k] <- b[k] * carbspline[i,k] * kz[k]
        }
        
        div_pred[i] ~ dnorm(mu[i], sd)
      }
    
    for(k in 2:4){
      kz[k] <- ifelse(k > nz + 1, 0, 1)
      theta_add[k] ~ dunif(0,500)
      theta[k] <- theta[k-1] + theta_add[k]
      b[k] ~ dnorm(0, 0.001)
    }
    kz[1] <- 1
    theta[1] <- 0
    b[1] ~ dnorm(0, 0.001)
    
    nz ~ dbin(kprob, 3)
    kprob ~ dunif(0,1)
    
    alpha ~ dnorm(0,0.001)
    sd ~ dunif(0,1)
    
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
      for(i in 1:nsite){
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
