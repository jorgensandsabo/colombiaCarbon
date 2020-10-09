########################################
## DATA PREPARATION OF RAW DATAFRAMES ##
########################################
# Code to prepare dataset
# Uses the clean versions of the raw datasets (errors fixed in the field or fixes traceable in 'Data_raw' folder)

# Use machine identifier to automatically set correct working path
if(file.exists('C:\\Users\\Jorgen\\Documents\\machine_identifier_lu847jp1o.txt')){dir.path <- "C:\\Users\\Jorgen\\OneDrive - Norwegian University of Life Sciences\\PhD"}
if(file.exists('C:\\Users\\jorgesan\\Documents\\machine_identifier_lksj7842.txt')){dir.path <- "C:\\Users\\jorgesan\\OneDrive - Norwegian University of Life Sciences\\PhD"}

setwd(dir.path)

library(dplyr)

Deaddata <- read.csv("Data\\vegetation\\Dead_raw.csv")
Grass <- read.csv("Data\\vegetation\\Grass_raw.csv")
Paramo <- read.csv("Data\\vegetation\\Paramo_raw.csv")
PTrees <- read.csv("Data\\vegetation\\PTrees_raw.csv")
Trees <- read.csv("Data\\vegetation\\Trees_raw.csv")
Plotdata <- read.csv("Data\\vegetation\\Plots_raw.csv", stringsAsFactors = FALSE)
metadata <- read.csv("Data\\vegetation\\RAWdata_samplingpointsCO_MASTER.csv", stringsAsFactors = FALSE)

######################
### MISSING VALUES ###
######################
# Code to check for and replace missing values
# Deadwood data
levels(Deaddata$SF)[which(levels(Deaddata$SF)=="")] <- NA
levels(Deaddata$Shape)[which(levels(Deaddata$Shape)=="")] <- NA

Deaddata[which(is.na(Deaddata$SiteCode)),]
Deaddata[which(is.na(Deaddata$N)),]
Deaddata[which(is.na(Deaddata$D1)),]

Deaddata[which(is.na(Deaddata$SF)),]
Deaddata[which(is.na(Deaddata$Shape)),]
Deaddata[which(is.na(Deaddata$Length)),]
Deaddata[which(is.na(Deaddata$DecayClass)),]

## SF: Doesn't matter

## Shape: The only standing tree is short --> assume cylinder if the tree only has one diameter measurement
Deaddata[which(is.na(Deaddata$Shape) & is.na(Deaddata$D2)),]$Shape <- "Cyl"
Deaddata[which(is.na(Deaddata$Shape) & !is.na(Deaddata$D2)),]$ Shape <- "Cone"

## Length: The only standing tree is decayclass 5 --> Don't expect it to be a full tree. Only three NAs. Use average length, expect it to be better than removing the pieces.
Deaddata[which(is.na(Deaddata$Length)),]$Length <- mean(filter(Deaddata,Length != 0)$Length)

## DecayClass: Plots often had several pieces of similar class (e.g. treefall event)
## Use the rounded plot median (note that median and mean give same values except for two small pieces)
for(i in 1:nrow(Deaddata)){
  if(is.na(Deaddata$DecayClass[i])){
    Deaddata$DecayClass[i] <- round(mean(filter(Deaddata,SiteCode==Deaddata$SiteCode[i])$DecayClass,na.rm=TRUE))
  }
}

# Grass data
Grass[which(is.na(Grass$SiteCode)),]
Grass[which(is.na(Grass$kgm2)),]

# Paramo data
Paramo[which(is.na(Paramo$SiteCode)),]
Paramo[which(is.na(Paramo$PlantN)),]
Paramo[which(is.na(Paramo$Type)),]

Paramo[which(is.na(Paramo$Diameter)),]
Paramo[which(is.na(Paramo$Height)),]

# Pasture tree data
PTrees[which(is.na(PTrees$SiteCode)),]
PTrees[which(is.na(PTrees$TreeN)),]
PTrees[which(is.na(PTrees$DBH)),]
PTrees[which(is.na(PTrees$Dist)),]

# Tree data
Trees[which(is.na(Trees$SiteCode)),]
Trees[which(is.na(Trees$TreeN)),]
Trees[which(is.na(Trees$DBH)),]

## TreeN: Add as the last tree number
for(i in 1:nrow(Trees)){
  if(is.na(Trees$TreeN[i])){
    Trees$TreeN[i] <- max(filter(Trees,SiteCode==Trees$SiteCode[i])$TreeN,na.rm=TRUE)
  }
}

#############################
### WOOD SPECIFIC GRAVITY ### REMOVED
#############################
# Code to fill missing wood specific gravity values
# Trees$Species <- as.character(Trees$Species)
# Trees <- left_join(Trees,Plotdata[c("SiteCode","AreaCode","Cluster")],by="SiteCode")

# Change WSG values (remove outliers for two species - remove all below 1.0 - remove palms/ferns
# Trees[which(Trees$Species=="Quercus Humboldtii" & Trees$WSG==min(Trees[which(Trees$Species=="Quercus Humboldtii"),]$WSG,na.rm=TRUE)),]$WSG <- NA
# Trees[which(Trees$Species=="TAF12sp1" & Trees$WSG==max(Trees[which(Trees$Species=="TAF12sp1"),]$WSG,na.rm=TRUE)),]$WSG <- NA
# Trees[which(Trees$WSG==1.0),]$WSG <- NA
# Trees[which(Trees$Species=="Palm"),]$WSG <- NA
# Trees[which(Trees$Species=="Fern"),]$WSG <- NA

# Averages
# AreaAvg <- as.data.frame(Trees %>% group_by(AreaCode) %>% summarise(n=n(),cores=length(which(!is.na(WSG))),perc=length(which(!is.na(WSG)))/n(),Areaavg=mean(WSG,na.rm=TRUE)))
# SiteAvg <- as.data.frame(Trees %>% group_by(SiteCode) %>% summarise(n=n(),cores=length(which(!is.na(WSG))),perc=length(which(!is.na(WSG)))/n(),Siteavg=mean(WSG,na.rm=TRUE)))
# ClusAvg <- as.data.frame(Trees %>% group_by(Cluster) %>% summarise(n=n(),cores=length(which(!is.na(WSG))),perc=length(which(!is.na(WSG)))/n(),Clusavg=mean(WSG,na.rm=TRUE)))
# SpecAvg <- as.data.frame(Trees %>% group_by(Species) %>% summarise(n=n(),cores=length(which(!is.na(WSG))),perc=length(which(!is.na(WSG)))/n(),Specavg=mean(WSG,na.rm=TRUE)))

    # This codeblock plots the differences between averages against number of cores. 
    # Choose minimum core number by visual assessment: ~ stable variance of the difference = sufficient cores to represent the plot/cluster
#    ClusAvg$AreaCode <- NA
#    for(i in 1:nrow(ClusAvg)){ClusAvg$AreaCode[i] <- as.character(Plotdata[which(Plotdata$Cluster == ClusAvg$Cluster[i]),]$AreaCode[1])}
#    ClusAvg <- left_join(ClusAvg,AreaAvg[c("AreaCode","Areaavg")],by="AreaCode")
#    ClusAvg$diff <- ClusAvg$Areaavg-ClusAvg$Clusavg
#    plot(ClusAvg$cores,ClusAvg$diff);abline(h=0,v=30)
    
#    SiteAvg$Cluster <- NA
#    for(i in 1:nrow(SiteAvg)){SiteAvg$Cluster[i] <- as.character(Plotdata[which(Plotdata$SiteCode == SiteAvg$SiteCode[i]),]$Cluster[1])}
#    SiteAvg <- left_join(SiteAvg,ClusAvg[c("Cluster","Clusavg")],by="Cluster")
#    SiteAvg$diff <- SiteAvg$Clusavg-SiteAvg$Siteavg
#    plot(SiteAvg$cores,SiteAvg$diff):abline(h=0,v=10)
    
#    SiteAvg$AreaCode <- NA
#    for(i in 1:nrow(SiteAvg)){SiteAvg$AreaCode[i] <- as.character(Plotdata[which(Plotdata$SiteCode == SiteAvg$SiteCode[i]),]$AreaCode[1])}
#    SiteAvg <- left_join(SiteAvg,AreaAvg[c("AreaCode","Areaavg")],by="AreaCode")
#    SiteAvg$diffA <- SiteAvg$Areaavg-SiteAvg$Siteavg
#    plot(SiteAvg$cores,SiteAvg$diffA);abline(h=0,v=10)
    # Conclusion: > 10 cores for SiteCode, > 30 cores for Cluster

# Shorten by rules (e.g. what plot WSGs should override area WSGs)
#ClusAvg <- ClusAvg[-which(ClusAvg$cores <= 30 | ClusAvg$cores == ClusAvg$n),]
#SiteAvg <- SiteAvg[-which(SiteAvg$cores <= 10 | SiteAvg$cores == SiteAvg$n),]
#SpecAvg <- SpecAvg[-which(SpecAvg$cores < 2 | SpecAvg$cores == SpecAvg$n),]
#SpecAvg <- filter(SpecAvg,!is.na(Species))

# Fill NA values
#Project average
#Trees$WSG_used <- mean(Trees$WSG,na.rm=TRUE)
# Area average
#for(i in 1:nrow(Trees)){if(Trees$AreaCode[i] %in% AreaAvg$AreaCode){Trees$WSG_used[i] <- AreaAvg[which(AreaAvg$AreaCode == Trees$AreaCode[i]),]$Areaavg}}
# Cluster average
#for(i in 1:nrow(Trees)){if(Trees$Cluster[i] %in% ClusAvg$Cluster){Trees$WSG_used[i] <- ClusAvg[which(ClusAvg$Cluster == Trees$Cluster[i]),]$Clusavg}}
# Site average
#for(i in 1:nrow(Trees)){if(Trees$SiteCode[i] %in% SiteAvg$SiteCode){Trees$WSG_used[i] <- SiteAvg[which(SiteAvg$SiteCode == Trees$SiteCode[i]),]$Siteavg}}
# Species average
#for(i in 1:nrow(Trees)){if(Trees$Species[i] %in% SpecAvg$Species){Trees$WSG_used[i] <- SpecAvg[which(SpecAvg$Species == Trees$Species[i]),]$Specavg}}
# Measured WSG
#for(i in 1:nrow(Trees)){if(!is.na(Trees$WSG[i])){Trees$WSG_used[i] <- Trees$WSG[i]}}

# Remove for Palms and Ferns
#Trees[which(Trees$Species %in% c("Palm","Fern")),]$WSG_used <- NA


#### REPEAT FOR PTREES
#PTrees$Species <- as.character(PTrees$Species)
#PTrees <- left_join(PTrees,Plotdata[c("SiteCode","AreaCode","Cluster")],by="SiteCode")

# Change WSG values
#PTrees[which(PTrees$WSG==1.0),]$WSG <- NA
#PTrees[which(PTrees$Species=="Palm"),]$WSG <- NA
#PTrees[which(PTrees$Species=="Fern"),]$WSG <- NA

# Averages
#AreaAvg <- as.data.frame(PTrees %>% group_by(AreaCode) %>% summarise(n=n(),cores=length(which(!is.na(WSG))),perc=length(which(!is.na(WSG)))/n(),Areaavg=mean(WSG,na.rm=TRUE)))
#SiteAvg <- as.data.frame(PTrees %>% group_by(SiteCode) %>% summarise(n=n(),cores=length(which(!is.na(WSG))),perc=length(which(!is.na(WSG)))/n(),Clusavg=mean(WSG,na.rm=TRUE)))
#ClusAvg <- as.data.frame(PTrees %>% group_by(Cluster) %>% summarise(n=n(),cores=length(which(!is.na(WSG))),perc=length(which(!is.na(WSG)))/n(),Siteavg=mean(WSG,na.rm=TRUE)))
#SpecAvg <- as.data.frame(PTrees %>% group_by(Species) %>% summarise(n=n(),cores=length(which(!is.na(WSG))),perc=length(which(!is.na(WSG)))/n(),Specavg=mean(WSG,na.rm=TRUE)))

# Shorten by rules
#AreaAvg <- AreaAvg[-which(AreaAvg$cores < 2),]
#ClusAvg <- ClusAvg[-which(ClusAvg$cores <= 30 | ClusAvg$cores == ClusAvg$n),]
#SiteAvg <- SiteAvg[-which(SiteAvg$cores <= 10 | SiteAvg$cores == SiteAvg$n),]
#SpecAvg <- SpecAvg[-which(SpecAvg$cores == 0 | SpecAvg$cores == SpecAvg$n),]
#SpecAvg <- filter(SpecAvg,!is.na(Species))

# Fill NA values
#Project average
#PTrees$WSG_used <- mean(PTrees$WSG,na.rm=TRUE)
# Area average
#for(i in 1:nrow(PTrees)){if(PTrees$AreaCode[i] %in% AreaAvg$AreaCode){PTrees$WSG_used[i] <- AreaAvg[which(AreaAvg$AreaCode == PTrees$AreaCode[i]),]$Areaavg}}
# Cluster average
#for(i in 1:nrow(PTrees)){if(PTrees$Cluster[i] %in% ClusAvg$Cluster){PTrees$WSG_used[i] <- ClusAvg[which(ClusAvg$Cluster == PTrees$Cluster[i]),]$Clusavg}}
# Site average
#for(i in 1:nrow(PTrees)){if(PTrees$SiteCode[i] %in% SiteAvg$SiteCode){PTrees$WSG_used[i] <- SiteAvg[which(SiteAvg$SiteCode == PTrees$SiteCode[i]),]$Siteavg}}
# Species average
#for(i in 1:nrow(PTrees)){if(PTrees$Species[i] %in% SpecAvg$Species){PTrees$WSG_used[i] <- SpecAvg[which(SpecAvg$Species == PTrees$Species[i]),]$Specavg}}
# Measured WSG
#for(i in 1:nrow(PTrees)){if(!is.na(PTrees$WSG[i])){PTrees$WSG_used[i] <- PTrees$WSG[i]}}

# Remove for Palms and Ferns
#PTrees[which(PTrees$Species %in% c("Palm","Fern")),]$WSG_used <- NA

## Remove temporary objects from the environment
#rm(AreaAvg,ClusAvg,SiteAvg,SpecAvg,i)

#############################
### DBH / HEIGHT / VOLUME ###
#############################
# Code to change DBH values of trees measured at other heights (POM = point of measurement)
# estimate height of trees and standing deadwood
# and estimate volume of deadwood

# Treeheight from DBH:H equation (Feldpausch 2011)
#Trees <- dplyr::mutate(Trees, Height = exp(1.2229+0.5320*log(DBH))*100)
#PTrees <- dplyr::mutate(PTrees, Height = exp(1.2229+0.5320*log(DBH))*100)

# DBH (for trees with POM=/=1.3) from taper function (Chambers 2000)
Trees <- dplyr::mutate(Trees, DBH_used = dplyr::case_when(
  is.na(POM) ~ DBH,
  !is.na(POM) ~ DBH / (1.59 * ((POM * 100)^-0.091))
))

PTrees <- dplyr::mutate(PTrees, DBH_used = dplyr::case_when(
  is.na(POM) ~ DBH,
  !is.na(POM) ~ DBH / (1.59 * ((POM * 100)^-0.091))
))

# Remove for fern/palms
Trees[which(Trees$Species %in% c("Palm","Fern")),]$DBH_used <- Trees[which(Trees$Species %in% c("Palm","Fern")),]$DBH
#Trees[which(Trees$Species %in% c("Palm","Fern")),]$Height <- NA

PTrees[which(PTrees$Species %in% c("Palm","Fern")),]$DBH_used <- PTrees[which(PTrees$Species %in% c("Palm","Fern")),]$DBH
#PTrees[which(PTrees$Species %in% c("Palm","Fern")),]$Height <- NA

#### Deadwood
## Add common length column - estimate for full trees (Feldpausch 2011)
Deaddata <- dplyr::mutate(Deaddata,L_used = dplyr::case_when(
  Length==0 ~ exp(1.2229+0.5320*log(D1))*100,
  Length!=0 ~ as.numeric(as.character(Length))))

## Add volume - taper function from Chambers et al. 2000
Deaddata <- dplyr::mutate(Deaddata,Volume = dplyr::case_when(
  Shape=="Cyl" ~ pi*L_used*(D1/2)^2,
  Shape=="Cone" ~ (pi*L_used*((D1/2)^2+(D1/2)*(D2/2)+(D2/2)^2))/3,
  Shape=="Tap" ~ (pi*L_used*(((1.59*D1*0.05^-0.091)/2)+((1.59*D1*0.05^-0.091)/2)*((1.59*(D1)*(L_used^-0.091))/2)+((1.59*(D1)*(L_used^-0.091))/2)^2))/3))

##########################
### BIOMASS ESTIMATION ###
##########################
# Trees
#Trees <- dplyr::mutate(Trees,
#                          BMAlvII2 = exp(3.103-1.794*log(DBH_used)+1.290*(log(DBH_used))^2-0.128*(log(DBH_used))^3+0.819*log(WSG_used)),
#                          BMAlvDry = exp(3.652-1.697*log(DBH_used)+1.169*(log(DBH_used))^2-0.122*(log(DBH_used))^3+1.285*log(WSG_used)),  
#                          ### Palm equation (Goodman et al. 2013)
#                          BMpalm = dplyr::case_when(Species=="Palm" ~ exp(-3.3488+2.7483*log(DBH_used))),
#                          ### Polylepis equation (Espinoza & Quispe 2005)
#                          BMpoly = dplyr::case_when(Species=="Polylepis" ~ 0.0694*(DBH_used^2.35996)))

#Trees <- dplyr::mutate(Trees,BMfinal = case_when(
#  !is.na(BMpalm) ~ BMpalm,
#  !is.na(BMpoly) ~ BMpoly,
#  AreaCode == "CC" ~ BMAlvDry,
#  !is.na(BMAlvII2) ~ BMAlvII2
#))

# PTrees
#PTrees <- dplyr::mutate(PTrees,
#                       BMAlvII2 = exp(3.103-1.794*log(DBH_used)+1.290*(log(DBH_used))^2-0.128*(log(DBH_used))^3+0.819*log(WSG_used)),
#                       BMAlvDry = exp(3.652-1.697*log(DBH_used)+1.169*(log(DBH_used))^2-0.122*(log(DBH_used))^3+1.285*log(WSG_used)),
#                       ### Palm equation (Goodman et al. 2013)
#                       BMpalm = dplyr::case_when(Species=="Palm" ~ exp(-3.3488+2.7483*log(DBH_used))),
#                       ### Polylepis equation (Espinoza & Quispe 2005)
#                       BMpoly = dplyr::case_when(Species=="Polylepis" ~ 0.0694*(DBH_used^2.35996)))

#PTrees <- dplyr::mutate(PTrees,BMfinal = case_when(
#  !is.na(BMpalm) ~ BMpalm,
#  !is.na(BMpoly) ~ BMpoly,
#  AreaCode == "CC" ~ BMAlvDry,
#  !is.na(BMAlvII2) ~ BMAlvII2
#))

# Deadwood (values Pfeifer et. al. 2015)
Deaddata <- dplyr::mutate(Deaddata,Biomass = dplyr::case_when(
  DecayClass==1 ~ (Volume*0.40)/1000,
  DecayClass==2 ~ (Volume*0.58)/1000,
  DecayClass==3 ~ (Volume*0.37)/1000,
  DecayClass==4 ~ (Volume*0.26)/1000,
  DecayClass==5 ~ (Volume*0.16)/1000))

# Paramo
Paramo <- dplyr::mutate(Paramo,AGB = dplyr::case_when(
  Type=="Espeletia" ~  exp(0.23)*(Diameter*Height)^0.38,
  Type=="Shrub" ~ exp(0.94)*((pi*(Diameter)^2/4)*Height)^0.21),
  BGB = dplyr::case_when(
    Type=="Espeletia" ~ exp(0.24)*(Diameter*Height)^0.54,
    Type=="Shrub" ~ exp(2.23)*(Diameter*Height)^0.50),
  TBM = dplyr::case_when(
    Type=="Espeletia" ~ exp(1.55)*((pi*(Diameter)^2/4)*Height)^0.32,
    Type=="Shrub" ~ exp(1.12)*((pi*(Diameter)^2/4)*Height)^0.22))

# Grass
Grass$GrassAGB <- (Grass$kgm2/1000)*10000


####################
### PLOT BIOMASS ###
####################
### Plot biomass dataframe
#Trees
#PlotBM <- dplyr::left_join(Plotdata,dplyr::summarise(dplyr::group_by(Trees,SiteCode),TreeAGB = sum(BMfinal,na.rm=TRUE)))
#PlotBM[which(is.na(PlotBM$TreeAGB) & PlotBM$Plot == "X"),]$TreeAGB <- 0
#PlotBM$TreeAGB <- ((PlotBM$TreeAGB/PlotBM$Size)*10000)/1000

#PTrees
#PlotBM <- dplyr::left_join(PlotBM, dplyr::mutate(dplyr::left_join(dplyr::left_join(PTrees %>% dplyr::group_by(SiteCode) %>% dplyr::summarise(PtBM = sum(BMfinal,na.rm=TRUE)),
#                                                                                   PTrees %>% dplyr::group_by(SiteCode) %>% dplyr::summarise(Dist = max(Dist))),
#                                                                  PTrees %>% dplyr::group_by(SiteCode) %>% dplyr::summarise(NTrees = length(unique(TreeN)))),
#                                                 PTha = dplyr::case_when(NTrees < 10 ~ (PtBM/(pi*(50^2)))*10000,
#                                                                         NTrees == 10 ~ (PtBM/(pi*(Dist)^2))*10000)
#)[,c("SiteCode","PTha")],by="SiteCode")
#PlotBM[which(is.na(PlotBM$PTha) & PlotBM$Ptrees == "X"),]$PTha <- 0
#PlotBM$PTha <- PlotBM$PTha/1000

#Paramo
#PlotBM <- dplyr::left_join(PlotBM,dplyr::summarise(dplyr::group_by(Paramo,SiteCode),ParamoAGB = sum(AGB,na.rm=TRUE)))
#PlotBM[which(is.na(PlotBM$ParamoAGB) & PlotBM$Plot == "X"),]$ParamoAGB <- 0
#PlotBM$ParamoAGB <- ((PlotBM$ParamoAGB/PlotBM$Size)*10000)/1000

#Dead
#PlotBM <- dplyr::left_join(PlotBM,dplyr::summarise(dplyr::group_by(Deaddata,SiteCode),DeadAGB = sum(Biomass,na.rm=TRUE)))
#PlotBM[which(is.na(PlotBM$DeadAGB) & PlotBM$Plot == "X"),]$DeadAGB <- 0
#PlotBM$DeadAGB <- ((PlotBM$DeadAGB/PlotBM$Size)*10000)/1000

#Grass
#PlotBM <- dplyr::left_join(PlotBM,Grass[c("SiteCode","GrassAGB")])
#PlotBM[which(is.na(PlotBM$GrassAGB) & PlotBM$Grass == "X"),]$GrassAGB <- 0

#PlotBM <- PlotBM[,c("AreaCode","Cluster","SiteCode","HabitatS","HabitatP","Ptrees","Grass","Plot","Size","Dataset",
#                    "FAge","TreeAGB","PTha","ParamoAGB","DeadAGB","GrassAGB")]



#####################################
### COMBINE PLOTDATA AND METADATA ###
#####################################
### Name fixes
metadata[which(metadata$point_id == "CHA10d"),]$point_id <- "CHA10D"
metadata[which(metadata$point_id == "CHA11d"),]$point_id <- "CHA11D"
metadata[which(metadata$point_id == "CHA12d"),]$point_id <- "CHA12D"

metadata[which(metadata$site == "Mesenia"),]$site <- "ChA"
metadata[which(metadata$site == "LasTangaras"),]$site <- "ChC" 
metadata[which(metadata$site == "Montezuma"),]$site <- "ChCore"

metadata <- mutate(metadata, region = case_when(
  region == "oriental" ~ paste(region,mountain_slope,sep=" "),
  TRUE ~ region
))
metadata[which(metadata$region == "oriental na"),]$region <- "magdalenavalley"

### Select columns and rename
plotsdf <- metadata[,c("point_id","site","cluster","lat","long","region")]
names(plotsdf)[c(1,2,3)] <- c("SiteCode","AreaCode","Cluster")

### Add missing cluster numbers
Plotdata$clusterN <- as.numeric(gsub("[^0-9.]", "",  Plotdata$Cluster))
plotsdf <- dplyr::left_join(plotsdf, Plotdata[,c("SiteCode","clusterN")])

plotsdf[which(plotsdf$Cluster != plotsdf$clusterN),]$Cluster <- plotsdf[which(plotsdf$Cluster != plotsdf$clusterN),]$clusterN
plotsdf[which(is.na(plotsdf$Cluster)),]$Cluster <- plotsdf[which(is.na(plotsdf$Cluster)),]$clusterN
plotsdf[plotsdf$SiteCode %in% c("FVF1","FVF2","FVF3","ORF1","ORF2","ORF3","RCF1","RCF2","RCF3"),]$Cluster <- 1
plotsdf[plotsdf$SiteCode %in% c("FVP1","FVP2","FVP3","ORP1","ORP2","ORP3","RCF4","RCF5","RCF6"),]$Cluster <- 2
plotsdf[plotsdf$SiteCode %in% c("RCP1","RCP2","RCP3"),]$Cluster <- 3
plotsdf[plotsdf$SiteCode %in% c("RCP4","RCP5","RCP6"),]$Cluster <- 4

plotsdf[which(is.na(plotsdf$Cluster)),]
plotsdf$clusterN <- NULL

### Combine spatial data and plot data
plotsdf <- dplyr::left_join(plotsdf, Plotdata[,c("SiteCode","Size","HabitatS","HabitatP","HabitatChoco","Dataset","FAge","Ptrees","Grass","Plot")])
rm(Plotdata)
plotsdf$Cluster <- paste(plotsdf$AreaCode,plotsdf$Cluster,sep="")

####################
## WRITE DATASETS ##
####################
write.csv(Deaddata,"Data\\vegetation\\Dead_prep.csv",row.names=F)
write.csv(Grass,"Data\\vegetation\\Grass_prep.csv",row.names=F)
write.csv(Paramo,"Data\\vegetation\\Paramo_prep.csv",row.names=F)
write.csv(plotsdf,"Data\\vegetation\\Plots_prep.csv",row.names=F)
write.csv(PTrees,"Data\\vegetation\\PTrees_prep.csv",row.names=F)
write.csv(Trees,"Data\\vegetation\\Trees_prep.csv",row.names=F)






