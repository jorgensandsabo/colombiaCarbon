# Temporary script to do first checks of dung beetle data

library(dplyr)
library(lme4)

# Use machine identifier to automatically set directory with points file
jorgen.desktop <- file.exists('C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD\\machine_identifier_lu847jp1o.txt')

if(jorgen.desktop){
  dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}

setwd(dir.path)

# Read file with cluster names
jorg <- read.csv(paste(dir.path,"\\Data\\raw\\Plots_BM.csv",sep=""))
names(jorg)[which(names(jorg)=="SiteCode")] <- "Cod_Loc"

# Read beetle data - summarise abudance for repeated species/day combinations
dung_beetles <- read.csv(paste(dir.path,"\\Data_raw\\Additional_data\\Diego_dung_beetles\\DBAndes_9-01-20_V1.3-Diego.csv",sep=""),stringsAsFactors = F)
dung_beetles$Cod_Loc <- as.character(dung_beetles$Cod_Loc)
dung_beetles[which(substr(dung_beetles$Cod_Loc,1,2) %in% c("BF","BP")),]$Cod_Loc <- paste("BR",substring(dung_beetles[which(substr(dung_beetles$Cod_Loc,1,2) %in% c("BF","BP")),]$Cod_Loc,2,3),sep="")
dung_beetles <- data.frame(dung_beetles %>% group_by(.dots=names(dung_beetles)[-which(names(dung_beetles) %in% c("Abundance","Preparation"))]) %>% summarise(Abundance=sum(Abundance)))

dung_beetles <- left_join(dung_beetles,jorg[,c("Cod_Loc","Cluster")], by="Cod_Loc")

# Add missing clusters
dung_beetles$Cluster <- as.character(dung_beetles$Cluster)
dung_beetles[which(is.na(dung_beetles$Cluster)),][which(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),]$Cod_Loc,4,5) %in% c(1:3)),]$Cluster <- paste(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),][which(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),]$Cod_Loc,4,5) %in% c(1:3)),]$Cod_Loc,1,3),1,sep="")
dung_beetles[which(is.na(dung_beetles$Cluster)),][which(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),]$Cod_Loc,4,5) %in% c(4:6)),]$Cluster <- paste(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),][which(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),]$Cod_Loc,4,5) %in% c(4:6)),]$Cod_Loc,1,3),2,sep="")
dung_beetles[which(is.na(dung_beetles$Cluster)),][which(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),]$Cod_Loc,4,5) %in% c(7:9)),]$Cluster <- paste(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),][which(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),]$Cod_Loc,4,5) %in% c(7:9)),]$Cod_Loc,1,3),3,sep="")
dung_beetles[which(is.na(dung_beetles$Cluster)),][which(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),]$Cod_Loc,4,5) %in% c(10:12)),]$Cluster <- paste(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),][which(substr(dung_beetles[which(is.na(dung_beetles$Cluster)),]$Cod_Loc,4,5) %in% c(10:12)),]$Cod_Loc,1,3),4,sep="")

#### Filter as wanted
db <- dung_beetles
db$Day <- as.character(db$Day)
db$Cientific.Name <- as.character(db$Cientific.Name)

db <- filter(db,Habitat=="Pasture")
db <- filter(db,Habitat=="Forest")
db <- filter(db,County=="Legizamo")

#### Temporary: remove rows with no species/day
db <- filter(db,!Day == "")
db <- filter(db,!Cientific.Name == "")


#### Convert to 3d table
dbtab <- table(db$Cod_Loc,db$Day,db$Cientific.Name)
#dbtab <- dbtab[,,-1]

    #### Total numbers
    nspecies <- length(dbtab[1,1,])
    nsites <- length(dbtab[,1,1])
    ntrapdays <- length(dbtab[1,,1])
    hab <- as.numeric(as.factor(as.data.frame(db %>% group_by(db$Cod_Loc) %>% summarise(hab=first(Habitat)))$hab))
    specname <- sort(unique(db$Cientific.Name))#;specname <- specname[-1]
    
for(i in 1:nspecies){
  for(j in 1:nsites){
    for(k in 1:ntrapdays){
      if(dbtab[j,k,i] == 1){
        dbtab[j,k,i] <- db[which(db$Cientific.Name == specname[i]
                                  & db$Cod_Loc == sort(unique(db$Cod_Loc))[j]
                                  & db$Day == sort(unique(db$Day))[k])
                            ,]$Abundance
      }}}}    


# Proportion of captures per trap-day for each species and trap
proptrap <- apply(dbtab, c(1,3), function(i) i/sum(i,na.rm=TRUE))
propindp <- apply(proptrap,c(2,3),mean)

# Relative change from day 1
test <- apply(dbtab, c(1,3), function(i) i - i[1])

  # Check:
  dbtab[58:68,1:4,"Uroxys microcularis"]
  test[1:4,58:68,"Uroxys microcularis"]
  
test2 <- reshape2::melt(test)
test2 <- filter(test2,value != 0)
par(mfrow=c(2,2))
hist(filter(test2,Var1=="R2")$value,breaks=1000)
hist(filter(test2,Var1=="R3")$value,breaks=1000)
hist(filter(test2,Var1=="R4")$value,breaks=1000)

nrow(filter(filter(test2,Var1 == "R2"),value>0))

## Decline across trapdays / baiting day?
proptab <- reshape2::melt(proptrap)
plot(proptab$value~proptab$Var1)
par(mfrow=c(2,2));hist(filter(proptab,Var1=="R1")$value);hist(filter(proptab,Var1=="R2")$value);hist(filter(proptab,Var1=="R3")$value);hist(filter(proptab,Var1=="R4")$value)

## Number of days where the species has been captured
trapdays <- proptrap ; trapdays[which(trapdays > 0)] <- 1
trapdays <- apply(trapdays,c(2,3),sum)
traptab <- reshape2::melt(trapdays)
table(traptab$value)

# Number of traps where the species has been captures in clusters where it is captured
names(traptab) <- c("Cod_Loc","Cientific.Name","captures")
clustab <- dplyr::left_join(traptab,unique(db[,c('Cod_Loc','Cluster')]))
clustab[which(clustab$captures > 1),]$captures <- 1
clustab <- as.data.frame(clustab %>% group_by(Cluster,Cientific.Name) %>% summarise(captures = sum(captures,na.rm=TRUE)))
table(clustab$captures)

# Distribution of counts
hist(reshape2::melt(dbtab)$value)
# Distribution of counts in clusters where the species is observed
Ntrap <- left_join(reshape2::melt(dbtab),unique(db[,c('Cod_Loc','Cluster')]),by=c("Var1"="Cod_Loc"))
Ntrap <- left_join(Ntrap,clustab,by=c("Cluster", "Var3"="Cientific.Name"))
Ntrap <- Ntrap[which(Ntrap$captures>0),]
hist(Ntrap$value,breaks=1000)

###### GLM tests
Ntrap2 <- dplyr::left_join(Ntrap,unique(db[c("Cod_Loc","Habitat")]),by=c("Var1"="Cod_Loc"))
Ntrap2$index <- c(1:nrow(Ntrap2))
Ntrap2$rebait <- ifelse(Ntrap2$Var2 %in% c("R1","R3"),1,2)

# Intercept abundance models with species/habitat random effects
poismod <- lme4::glmer(value ~ 1 + (1|Var3/Habitat), data=Ntrap2, family=poisson)                  # Poisson
poisovrdisp <- lme4::glmer(value ~ 1 + (1|Var3/Habitat) + (1|index),data=Ntrap2, family=poisson)   # Overdispersed poisson
nbmod <- lme4::glmer.nb(value ~ 1 + (1|Var3/Habitat),data=Ntrap2)                                  # Negative binomial

# Model diagnostics
#devtools::install_github(repo = "DHARMa", username = "florianhartig", subdir = "DHARMa")
modfit_res_overdisp <- DHARMa::simulateResiduals(poisovrdisp, refit=T)
DHARMa::testDispersion(modfit_res_overdisp)
DHARMa::testZeroInflation(modfit_res_overdisp)
DHARMa::plotSimulatedResiduals(modfit_res_overdisp)


###### Abundance model, no detection
model1.string <- "
  model {

  # Process model
    for(i in 1:nspecies){
      for(j in 1:nsites){
        for(k in 1:ntrapdays){
          
          #count[j,k,i] ~ dpois(lambda[j,k,i]*ORLE[j,k,i])
          log(lambda[j,k,i]) <- alpha[i] + beta1[i,hab[j]] + beta2[i,hab[j]] * rebait[k]

          # Zero-inflation:
          count[j,k,i] ~ dpois(mu[j,k,i])
          mu[j,k,i] <- lambda[j,k,i] *  z[j,k,i] + 0.000001
          
          z[j,k,i] ~ dbern(psi[i,hab[j]])
        }
      }
    }

  # Priors
  mua ~ dnorm(0,1)
  sda ~ dunif(0,10)
  taua <- pow(sda,-2)
  
  mub1 ~ dnorm(0,1)
  sdb1 ~ dunif(0,10)
  taub1 <- pow(sdb1,-2)
  
  mub2 ~ dnorm(0,1)
  sdb2 ~ dunif(0,10)
  taub2 <- pow(sdb2,-2)
  
  for(i in 1:nspecies){
    alpha[i] ~ dnorm(mua,taua)
    
    for(h in 1:2){
      beta1[i,h] ~ dnorm(mub1,taub1)
      beta2[i,h] ~ dnorm(mub2,taub2)
    }
  }
  
  # Priors zero-inflation term
  for(i in 1:nspecies){
    for(m in 1:2){
        psi[i,m] ~ dunif(0,1)
    }
  }
  
  # Priors for overdispersion term
 # tauORLE ~ dgamma(1,1)

  #for(i in 1:nspecies){
   # for(j in 1:nsites){
    #  for(k in 1:ntrapdays){
     #     ORLE[j,k,i] ~ dgamma(1,1)
      #}
    #}
  #}
}"

sink("JAGS\\abunmod.txt")
cat(paste(model1.string))
sink()

DBmod_abunmod <- rjags::jags.model(file = "JAGS\\abunmod.txt",
                                   data=list(count = dbtab, hab=hab, rebait=c(0,1,0,1), nspecies= nspecies, nsites = nsites, ntrapdays = ntrapdays),
                                   #inits = list(abun = apply(dbtab,c(1,3),sum)),
                                   n.chains=3,n.adapt=10000)
update(object = DBmod_abunmod, n.iter = 10000)
DBmod_abunmod <- rjags::coda.samples(model = DBmod_abunmod,
                                     variable.names = c("lambda","psi","mua","sda","mub1","sdb1","mub2","sdb2","alpha","beta1","beta2"),
                                     n.iter=10000, thin = 1)

### Model checks
MCMCvis::MCMCplot(DBmod_abunmod)
MCMCvis::MCMCtrace(DBmod_abunmod,params=c("mua","sda","mub1","sdb1","mub2","sdb2","alpha","beta1","beta2"),wd="C:\\Users\\Jorgen\\Desktop",filename="DBmod_abunmod.pdf")
MCMCvis::MCMCtrace(DBmod_abunmod,params=c("lambda"),wd="C:\\Users\\Jorgen\\Desktop",filename="DBmod_abunmod_lambda.pdf")
MCMCvis::MCMCsummary(DBmod_abunmod)


###### Removal model
model1.string <- "
  model {

  # Process model
    for(i in 1:nspecies){
      for(j in 1:nsites){
        for(k in 1:ntrapdays){
        
          abun[j,k+1,i] <- abun[j,k,i] - count[j,k,i]
          count[j,k,i] ~ dbin(q[i,hab[j],rebait[k]],abun[j,k,i])
          #logit(q[i,j,k]) <- alpha[i,hab[j]] + beta[i,hab[j]] * rebait[k]
        }
        abun[j,1,i] ~  dpois(u[j,i])
      }
    }

  # Priors
  for(i in 1:nspecies){

    for(h in 1:2){
      #alpha[i,h] ~ dnorm(0,1)
      #beta[i,h] ~ dnorm(0,1)
      
      for(r in 1:2){
        q[i,h,r] ~ dbeta(3,3)
      }
    }

    for(j in 1:nsites){
      u[j,i] ~ dgamma(1,0.001)
    }
  }
}"

sink("JAGS\\removalmod.txt")
cat(paste(model1.string))
sink()

DBmod_removal <- rjags::jags.model(file = "JAGS\\removalmod.txt",
                                 data=list(count = dbtab, hab=hab, rebait=c(1,2,1,2), nspecies= nspecies, nsites = nsites, ntrapdays = ntrapdays),
                                 #inits = list(abun = apply(dbtab,c(1,3),sum)),
                                 n.chains=3,n.adapt=1000)
update(object = DBmod_removal, n.iter = 1000)
DBmod_removal <- rjags::coda.samples(model = DBmod_removal,
                                   variable.names = c("u","beta","alpha","q","abun"),
                                   n.iter=1000, thin = 1)

### Model checks
MCMCvis::MCMCplot(DBmod_removal)
MCMCvis::MCMCtrace(DBmod_removal,wd="C:\\Users\\Jorgen\\Desktop",filename="DBmod_removal.pdf")



###### Poisson abundance model - no shared species parameter distribution
model1.string <- "
  model {

  # Process model

    # State process
    for(i in 1:nspecies){
      for(j in 1:nsites){
        abun[j,i] ~ dpois(lambda[i,hab[j]])

        # Observation process
        for(k in 1:ntrapdays){
          count[j,k,i] ~ dbin(d[i,j,k],abun[j,i])
          logit(d[i,j,k]) <- alpha[i,hab[j]] + betad[i,k] * rebait[k]
        }
      }
    }

  # Priors
  for(i in 1:nspecies){
    for(k in 1:ntrapdays){
      betad[i,k] ~ dnorm(0,1)
    }
    for(m in 1:2){
        #d[i,m] ~ dbeta(1,1)
        alpha[i,m] ~ dnorm(0,1)
        lambda[i,m] ~ dgamma(1,0.001)
    }
  }
}"

sink("JAGS\\Nmod_poisson1.txt")
cat(paste(model1.string))
sink()

DBmod_Pois1 <- rjags::jags.model(file = "JAGS\\Nmod_poisson1.txt",
                                data=list(count = dbtab, hab=hab, rebait=c(1,2,1,2), nspecies= nspecies, nsites = nsites, ntrapdays = ntrapdays),
                                inits =list(abun = apply(dbtab,c(1,3),sum)),
                                n.chains=3,n.adapt=1000)
update(object = DBmod_Pois1, n.iter = 1000)
DBmod_Pois1 <- rjags::coda.samples(model = DBmod_Pois1,
                                  variable.names = c("d","lambda","abun"),
                                  n.iter=1000, thin = 1)

### Model checks
MCMCvis::MCMCplot(DBmod_Pois1)
MCMCvis::MCMCtrace(DBmod_Pois1,wd="C:\\Users\\Jorgen\\Desktop",filename="DBmod_Pois1.pdf")

###### Poisson abundance model - shared species parameter distribution
model1.string <- "
  model {

  # Process model
  
    # State process
    for(i in 1:nspecies){
      for(j in 1:nsites){
        abun[j,i] ~ dpois(lambda[i,hab[j]])
        
        # Observation process
        for(k in 1:ntrapdays){
          count[j,k,i] ~ dbin(d[i,j,k],abun[j,i])
          logit(d[i,j,k]) <- alpha[i,hab[j]] + betad[i,k] * rebait[k]
        }
      }
    }

  # Priors
  meand ~ dnorm(0,1)
  sdd ~ dunif(0,10)
  taud <- pow(sdd,-2)
  meanl ~ dnorm(0,1)
  sdl ~ dunif(0,10)
  taul <- pow(sdl,-2)

  for(i in 1:nspecies){
    for(k in 1:ntrapdays){
      betad[i,k] ~ dnorm(0,1)
    }
    for(m in 1:2){
        alpha[i,m] ~ dnorm(meand,taud)
        #a[i,m] ~ dnorm(meand,taud)
        #logit(d[i,m]) <- a[i,m]
        b[i,m] ~ dnorm(meanl,taul)
        log(lambda[i, m]) <- b[i,m]
    }
  }
}"

sink("JAGS\\Nmod_poisson.txt")
cat(paste(model1.string))
sink()

DBmod_Pois <- rjags::jags.model(file = "JAGS\\Nmod_poisson.txt",
                              data=list(count = dbtab, hab=hab,rebait=c(1,2,1,2), nspecies= nspecies, nsites = nsites, ntrapdays = ntrapdays),
                              inits =list(abun = apply(dbtab,c(1,3),sum)),
                              n.chains=3,n.adapt=10000)
update(object = DBmod_Pois, n.iter = 10000)
DBmod_Pois <- rjags::coda.samples(model = DBmod_Pois,
                                variable.names = c("d","lambda","meand","sdd","meanl","sdl","abun"),
                                n.iter=10000, thin = 1)

### Model checks
MCMCvis::MCMCplot(DBmod_Pois)
MCMCvis::MCMCtrace(DBmod_Pois,wd="C:\\Users\\Jorgen\\Desktop",filename="DBmod_Pois1.pdf")




###### Overdispersed poisson abundance model
model1.string <- "
  model {

  # Process model
  
    # State process
    for(i in 1:nspecies){
      for(j in 1:nsites){
        # Observation process
        for(k in 1:ntrapdays){
          abun[j,k,i] ~ dpois(lambda[i,hab[j]] * ORLE[j,k,i])
          count[j,k,i] ~ dbin(d[i,j,k],abun[j,k,i])
          logit(d[i,j,k]) <- alpha[i,hab[j]] + betad[i,k] * rebait[k]
        }
      }
    }

  # Priors
  meand ~ dnorm(0,1)
  sdd ~ dgamma(1,0.001)
  taud <- pow(sdd,-2)
  meanl ~ dnorm(0,0.001)
  sdl ~ dunif(0,10)
  taul <- pow(sdl,-2)
  
  for(i in 1:nspecies){
    for(k in 1:ntrapdays){
      betad[i,k] ~ dnorm(0,1)
    }
    for(m in 1:2){
        #a[i,m] ~ dnorm(meand,1/sdd)
        #logit(d[i,m]) <- a[i,m]
        alpha[i,m] ~ dnorm(meand,taud)
        b[i,m] ~ dnorm(meanl,taul)
        log(lambda[i, m]) <- b[i,m]
    }
  }

  # Priors for overdispersion term
  tauORLE ~ dgamma(0.001, 0.001)

  for(i in 1:nspecies){
    for(j in 1:nsites){
      for(k in 1:ntrapdays){
          ORLE[j,k,i] ~ dnorm(0,tauORLE)
      }
    }
  }
}"

sink("JAGS\\Nmod_overdisp.txt")
cat(paste(model1.string))
sink()

DBmod_OD <- rjags::jags.model(file = "JAGS\\Nmod_overdisp.txt",
                                data=list(count = dbtab, hab=hab,rebait=c(1,2,1,2), nspecies= nspecies, nsites = nsites, ntrapdays = ntrapdays),
                                inits =list(abun=dbtab),
                                n.chains=3,n.adapt=10000)
update(object = DBmod_OD, n.iter = 10000)
DBmod_OD <- rjags::coda.samples(model = DBmod_OD,
                                 variable.names = c("d","lambda","meand","sdd"),
                                 n.iter=10000, thin = 1)

MCMCvis::MCMCtrace(DBmod_OD,wd="C:\\Users\\Jorgen\\Desktop",filename="DBmod_OD_trace.pdf")

###### Zero-inflated poisson abundance model
model1.string <- "
  model {

  # Process model
  
    # State process
    for(i in 1:nspecies){
      for(j in 1:nsites){
        # Observation process
        for(k in 1:ntrapdays){
          abun[j,k,i] ~ dpois(mu[j,k,i])
          count[j,k,i] ~ dbin(d[i,hab[j]],abun[j,k,i])
          
          mu[j,k,i] <- lambda[i,hab[j]]*z[j,k,i] + 0.00001
          z[j,k,i] ~ dbern(psi[i,hab[j]])
        }
      }
    }

  # Priors
  meand ~ dnorm(0,1)
  sdd ~ dgamma(1,0.001)
  meanl ~ dnorm(0,0.001)
  sdl ~ dunif(0,10)
  taul <- pow(sdl,-2)
  
  for(i in 1:nspecies){
    for(m in 1:2){
        a[i,m] ~ dnorm(meand,1/sdd)
        logit(d[i,m]) <- a[i,m]
        b[i,m] ~ dnorm(meanl,taul)
        log(lambda[i, m]) <- b[i,m]
        psi[i,m] ~ dunif(0,1)
    }
  }
}"

sink("JAGS\\Nmod_zeroinf.txt")
cat(paste(model1.string))
sink()

DBmod_ZI <- rjags::jags.model(file = "JAGS\\Nmod_zeroinf.txt",
                           data=list(count = dbtab, hab=hab, nspecies= nspecies, nsites = nsites, ntrapdays = ntrapdays),
                           inits =list(abun=dbtab),
                           n.chains=3,n.adapt=1000)
update(object = DBmod_ZI, n.iter = 1000)
DBmod_ZI <- rjags::coda.samples(model = DBmod_ZI,
                             variable.names = c("d","lambda","meand","sdd","meanl","sdl","psi"),
                             n.iter=1000, thin = 1)



###### Zero-inflated overdispered poisson abundance model
model1.string <- "
  model {

  # Process model
  
    # State process
    for(i in 1:nspecies){
      for(j in 1:nsites){
        # Observation process
        for(k in 1:ntrapdays){
          abun[j,k,i] ~ dpois(mu[j,k,i])
          count[j,k,i] ~ dbin(d[i,hab[j]],abun[j,k,i])

          mu[j,k,i] <- lambda[i,hab[j]]*ORLE[j,k,i]*z[j,k,i] + 0.00001
          z[j,k,i] ~ dbern(psi[i,hab[j]])
        }
      }
    }

  # Priors
  meand ~ dnorm(0,1)
  sdd ~ dgamma(1,0.001)
  meanl ~ dnorm(0,0.001)
  sdl ~ dunif(0,10)
  taul <- pow(sdl,-2)

  for(i in 1:nspecies){
    for(m in 1:2){
        a[i,m] ~ dnorm(meand,1/sdd)
        logit(d[i,m]) <- a[i,m]
        b[i,m] ~ dnorm(meanl,taul)
        log(lambda[i, m]) <- b[i,m]
        psi[i,m] ~ dunif(0,1)
    }
  }

  # Priors for overdispersion term
  tauORLE ~ dgamma(0.001, 0.001)

  for(i in 1:nspecies){
    for(j in 1:nsites){
      for(k in 1:ntrapdays){
          ORLE[j,k,i] ~ dnorm(0,tauORLE)
      }
    }
  }
}"

sink("JAGS\\Nmod_overdisp_zeroinf.txt")
cat(paste(model1.string))
sink()

DBmod_OD_ZI <- rjags::jags.model(file = "JAGS\\Nmod_overdisp_zeroinf.txt",
                              data=list(count = dbtab, hab=hab, nspecies= nspecies, nsites = nsites, ntrapdays = ntrapdays),
                              inits =list(abun=dbtab),
                              n.chains=3,n.adapt=1000)
update(object = DBmod_OD_ZI, n.iter = 1000)
DBmod_OD_ZI <- rjags::coda.samples(model = DBmod_OD_ZI,
                                variable.names = c("d","lambda","meand","sdd","meanl","sdl","psi","tauORLE"),
                                n.iter=1000, thin = 1)



### Model checks
MCMCvis::MCMCplot(DBmod_OD_ZI)
MCMCvis::MCMCtrace(DBmod_OD_ZI,wd="C:\\Users\\Jorgen\\Desktop",filename="DBmod_Pois1.pdf")
MCMCvis::MCMCsummary(DBmod_OD_ZI)










######## FROM OTHER SCRIPT, MAY USE
# Output
table <- cbind(summary(DBmod_OD)[[1]][,1:2],summary(DBmod_OD)[[2]][,c(1,5)])
round(table,3)

# Diagnostics
coda::gelman.plot(DBmod_OD)
coda::gelman.diag(DBmod_OD)
coda::autocorr.plot(DBmod_OD)
par(mfrow=c(3,4));coda::traceplot(DBmod_OD);par(mfrow=c(1,1))
coda::geweke.diag(DBmod_OD)
coda::geweke.plot(DBmod_OD)

# Model selection
mcall <- rbind(linmodWSG[[1]],linmodWSG[[2]],linmodWSG[[3]])
mc_ll <- mcall[,paste0("loglik[",1:nrow(plotF_wsg),"]")]
waicWSG <- loo::waic(mc_ll)
looWSG <- loo::loo(mc_ll,r_eff = loo::relative_eff(mc_ll,c(rep(1,1000),rep(2,1000),rep(3,1000))))

# Residuals vs. fitted
mc_st_resid <- mcall[,paste0("st_resid[",1:nrow(plotF_wsg),"]")]
mc_fit <- mcall[,paste0("mu[",1:nrow(plotF_wsg),"]")]
plot(colMeans(mc_fit),colMeans(mc_st_resid));abline(h=0,col="red")

car::qqPlot(colMeans(mc_st_resid),"norm")
