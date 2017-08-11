
# This code reads in the hourly occupancy table for dolphin detections and models
# daily occupancy using GEEGLM's

# 1)  Load packages and declare local functions      ##############

# This code investigates the various models that look at the different occupancy distributions
# incorporated by the data 
rm(list=ls())
library(boot)            # for inv.logit
library(mgcv)
library(ggplot2)
library(lme4)
library(dplyr)           # for distinct function 
library(geepack)         # To make the GEE's
#install_version("geepack", version = "1.0-7", repos = "http://cran.us.r-project.org") for anova.gee, which is gone- Email Soren 5/24/2017
library(splines)
library(RColorBrewer)
library(MuMIn)           # for QIC
library(MASS)            # for mvrnorm in boostrapping intervals 
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(mvtnorm)         # for rmvnorm used in predictions/plotting
library(geosphere)       # for distance between lat long points, estimating distance between river and locs

river_locs=data.frame(Rivername=c("Esk", "Dee", "South Esk", "Spey", "Tay", "Thurso", "Tweed", "Cromarty"),
                      Lat=c("555638", "570723", "564259", "574014", "562243", "583536", "554512", "5741454"),
                      Lon=c("-045656", "-035257"," -032659", "-045432", "-043518", "-033056", "-035341", "-03592429"))

river_locs$LatDeg=as.numeric(substr(river_locs$Lat, 1,2)) +
  as.numeric(substr(river_locs$Lat, 3,4))/60 + 
  as.numeric(substr(river_locs$Lat, 5,6))/60/60

river_locs$lonDeg=as.numeric(substr(river_locs$Lon, 1,3)) -
  as.numeric(substr(river_locs$Lon, 4,5))/60 - 
  as.numeric(substr(river_locs$Lon, 6,7))/60/60




# Calculate AUC, %0's id'ed and %1's ided, from Pirotta sperm whale paper- 
# modified 
CalcAUC<-function(mod, data_sub){
  
  pr <- predict(mod, data_sub, type="response")                          # the final model is used to predict the data on the response scale (i.e. a value between 0 and 1)
  pred <- prediction(pr, data_sub$BBOcc)                                    # to specify the vector of predictions (pr) and the vector of labels (i.e. the observed values "Pres")
  perf <- performance(pred, measure="tpr", x.measure="fpr")          # to assess model performance in the form of the true positive rate and the false positive rate
  plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5)) # to plot the ROC curve
  
  
  # Choice of the best cut-off probability
  
  y<-as.data.frame(perf@y.values)
  x<-as.data.frame(perf@x.values)
  fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45° line and the line joining the origin with the point (x;y) on the ROC curve
  L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
  d <- L*sin(fi)                                                     # to calculate the distance between the 45° line and the ROC curve
  # write.table(d,"C:\\distances.txt")                                 # to write a table with the computed distances
  
  # The table should then be opened in Microsoft Excel to find the maximum distance with the command "Sort", and the relative position (i.e. the number of the corresponding record)
  # MAX d= 0.1127967 --> position 39
  
  alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
  Best_cutoff=alpha[which.max(unlist(d)),]                                                       # to identify the alpha value (i.e. the cut-off) that corresponds to the maximum distance between the 45° line and the curve
  
  # Best cutoff:   0.3464173
  # This value can now be used to build the confusion matrix:
  
  DATA<-matrix(0,nrow(data_sub),3)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 919 - the number of rows can be checked with dim(dat)) 
  DATA<-as.data.frame(DATA)
  names(DATA)<-c("plotID","Observed","Predicted")
  DATA$plotID<-1:nrow(data_sub)                                                # the first column is filled with an ID value that is unique for each row
  DATA$Observed<-data_sub[,names(attr(modlist[[1]]$terms, 'dataClasses')[1])]                                            # the second column reports the observed response (0s and 1s)
  DATA$Predicted<-predict(mod,data_sub,type="response")                 # the third column reports the predictions
  cmx(DATA, threshold = Best_cutoff)                                   # the identified cut-off must be used here
  
  # Area under the Curve 
  auc <- unlist(performance(pred, measure="auc")@y.values)
  
  # Proportion of the presences correctly identified 
  pres=prop.table(cmx(DATA, threshold = Best_cutoff))[1,1]
  
  # Proportion of the absences correctly idenified
  abs=prop.table(cmx(DATA, threshold = Best_cutoff))[2,2]
  
  
  return(c(auc, pres, abs))
}




# colorblind palette with black for plotting
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




# 2)  Data Prep #################################################################################

# This section processes the hourly occupancy table to correct level orders (GroupID and UnitLoc)
# Hourly Bottlenose dolphin offsets are delcared where offset represnts probability that a click detection
# (OccAll) was produced by a broadband species. If only frequency banded clicks, P(BND)=0.06 if Broadband
# clicks P(BND)=0.77 and if unknown clicks P(BND)=0.50


OccTable= read.csv('W:/KJP PHD/4-Bayesian Habitat Use/R Code/OccupancyTable_ThreePdets.csv')
level_names=c( "Lat_05", "Lat_10", "Lat_15",
               "Hel_05", "Hel_10", "Hel_15",
               "Cro_05", "Cro_10", "Cro_15",
               "SpB_05", "SpB_10", "SpB_15",
               "Fra_05", "Fra_10", "Fra_15",
               "Cru_05", "Cru_10", "Cru_15",
               "Sto_05", "Sto_10", "Sto_15",
               "Abr_05", "Abr_10", "Abr_15",
               "StA_05", "StA_10", "StA_15",
               "Stb_05", "Stb_10", "Stb_15")

OccTable$UnitLoc=factor(OccTable$UnitLoc, levels=level_names)

OccTable$UnitLoc=(droplevels(OccTable)$UnitLoc)


OccTable$ShoreDist=as.character(OccTable$ShoreDist)
OccTable$ShoreDist[OccTable$ShoreDist=="5"]='05'

OccTable$ShoreDist=as.factor(OccTable$ShoreDist)

OccTable$Year=as.factor(OccTable$Year)
meta=read.csv('W:/KJP PHD/Deployment Information/CPODs for Kaitlin.csv')
meta$UnitLoc=factor(meta$UnitLoc, levels=level_names)



meta2=read.csv('W:\\KJP PHD\\Deployment Information\\SlopeAndAspect.csv')
meta2$UnitLoc=factor(meta2$UnitLoc, levels=level_names)
meta2$DistToSalmonRun=0
meta2$RiverName='blarg'

# calculate distance to nearest salmon river
for(ii in 1:30){
  temp=rep(0, nrow(river_locs))
  
  for(jj in 1:nrow(river_locs)){
    temp[jj]=distm (c(meta2$Lon[ii], meta2$Lat[ii]), c(river_locs$lonDeg[jj], river_locs$LatDeg[jj]), fun = distHaversine)
  }
  
  meta2$DistToSalmonRun[ii]=min(temp)
  meta2$RiverName[ii]=as.character(river_locs$Rivername[which.min(temp)])
  
}


meta=merge(meta, meta2, by='UnitLoc', all.x = TRUE)


meta_sub=subset(meta, select=c('UnitLoc', 'Slope2', 'DistToSalmonRun','RiverName', 'Depth_m' ))

OccTable=merge(OccTable, meta_sub, all.x = TRUE)
rm(meta, meta2, meta_sub)

OccTable$GroupId=unlist(strsplit(as.character(OccTable$UnitLoc), split = "_"))[seq(1,(nrow(OccTable)*2)-1,2)]

level_names=c( "Lat", "Hel", "Cro",
               "SpB", "Fra", "Cru",
               "Sto", "Abr", "StA",
               "Stb")

OccTable$GroupId=factor(OccTable$GroupId, levels=level_names)
OccTable$JD_scale=scale(OccTable$JulienDay)

OccTable$FBOcc[is.na(OccTable$FBOcc)]=0
OccTable$BBOcc[is.na(OccTable$BBOcc)]=0
OccTable$UNKOcc[is.na(OccTable$UNKOcc)]=0


# Species offset for Occupancy
# If two species same hour and same unit then uncertain
OccTable$SpeciesOffset=OccTable$BBOcc+OccTable$FBOcc+OccTable$UNKOcc
OccTable$SpeciesOffset[OccTable$SpeciesOffset>1]=0.5
OccTable$SpeciesOffset[OccTable$SpeciesOffset==1 & OccTable$BBOcc==1]=0.77
OccTable$SpeciesOffset[OccTable$SpeciesOffset==1 & OccTable$FBOcc==1]=0.06
OccTable$SpeciesOffset[OccTable$SpeciesOffset==1 & OccTable$UNKOcc==1]=0.5
#OccTable$SpeciesOffset[OccTable$SpeciesOffset==0] = 1

OccTable$BNDTotOffset=(OccTable$BBOcc*.77+OccTable$FBOcc*.06+OccTable$UNKOcc*.5)/(OccTable$BBOcc+OccTable$FBOcc+OccTable$UNKOcc)


OccTable$BNDTotOffset[is.nan(OccTable$BNDTotOffset)]=0


# 2a) Aggregate hourly Detections into Daily Occupancy table ######################################################

# This section of code aggregates the hourly occupancy into daily occupancy and
# Assigns BND offset scores to each day indicating the probability that on any given day,
# There was a bottlenose dolphin in the detections

# Offset scores where multiple species present calculated using a weighted average of the hourly 
# detection scores.(Next iteration migh consider likelihood ratios instead, but enough stuff to deal with at present)


# Add total number of detections 
mm.bbtot=as.data.frame(aggregate(BBOcc~UnitLoc+Date, FUN=sum, data = OccTable))
colnames(mm.bbtot)[3]='BBTot'
mm.fbtot=as.data.frame(aggregate(FBOcc~UnitLoc+Date, FUN=sum, data = OccTable))
colnames(mm.fbtot)[3]='FBTot'
mm.unktot=as.data.frame(aggregate(UNKOcc~UnitLoc+Date, FUN=sum, data = OccTable))
colnames(mm.unktot)[3]='UNKTot'


mm=distinct(OccTable, Date, UnitLoc, JulienDay, GroupId, ShoreDist, Slope2, Year, Month, RiverName, DistToSalmonRun, Depth_m)
mm.bb=distinct(subset(OccTable, BBOcc>0), Date, BBOcc, UnitLoc)
mm.fb=distinct(subset(OccTable, FBOcc>0), Date, FBOcc, UnitLoc)
mm.unk=distinct(subset(OccTable, UNKOcc>0), Date, UNKOcc, UnitLoc)


OccTable_daily=merge(mm, mm.bb, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.fb, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.unk, by = c('Date', 'UnitLoc'), all.x = TRUE)

OccTable_daily=merge(OccTable_daily, mm.bbtot, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.fbtot, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.unktot, by = c('Date', 'UnitLoc'), all.x = TRUE)


# Set NA scores to 0
OccTable_daily[is.na(OccTable_daily)] <- 0


OccTable_daily$SpeciesOffset=OccTable_daily$BBOcc+OccTable_daily$FBOcc+OccTable_daily$UNKOcc
OccTable_daily$OccAll=ifelse(OccTable_daily$SpeciesOffset>=1,1,0)

# Species offset for Daily Occupancy
# If two species same day and same unit then uncertain
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset>1]=0.5
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==1 & OccTable_daily$BBOcc==1]=0.77
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==1 & OccTable_daily$FBOcc==1]=0.06
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==1 & OccTable_daily$UNKOcc==1]=0.5
#OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==0] = .01

# Total offset
OccTable_daily$BNDTotOffset=(OccTable_daily$BBTot*.77+OccTable_daily$FBTot*.06+OccTable_daily$UNKTot*.5)/
  (OccTable_daily$BBTot+OccTable_daily$FBTot+OccTable_daily$UNKTot)

OccTable_daily$TotDet=(OccTable_daily$BBTot+OccTable_daily$FBTot+OccTable_daily$UNKTot)

OccTable_daily$BNDTotOffset[is.na(OccTable_daily$BNDTotOffset)]=0

rm(mm, mm.bb, mm.fb, mm.unk, mm.bbtot, mm.fbtot, mm.unktot)


# 3)  Data Exploration/Viz ################################################################################

# Much of this was done elsewhere, cursory exploration of colinearity between deployment distance from
# shore factor(ShoreDist) and the slope of the bathymetry. 

meta$ShoreDist=substr(meta$UnitLoc,5,6)

boxplot(meta$Slope2~meta$ShoreDist)

library(heplots) # for eta
model.aov <- aov(Slope2 ~ ShoreDist, data = meta)
summary(model.aov)
etasq(model.aov, partial = FALSE)

# Proportion of variance explained is .65
# https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable

rm(meta_sub, meta)

# Bin the data for visualisation #
OccTable_daily$DayBin=cut(OccTable_daily$JulienDay, breaks=20)


mm=data.frame(aggregate(data=OccTable_daily, BBOcc~DayBin+GroupId+ShoreDist, FUN=mean))
mm=cbind(mm, aggregate(data=OccTable_daily, JulienDay~DayBin+GroupId+ShoreDist, FUN=median)[,4])
mm=cbind(mm, aggregate(data=OccTable_daily, BBOcc~DayBin+GroupId+ShoreDist, FUN=sum)[,4])
mm=cbind(mm, aggregate(data=OccTable_daily, BBOcc~DayBin+GroupId+ShoreDist, FUN=length)[,4])

colnames(mm)[5:7]=c('med', 'sum', 'n')

library('Hmisc')
mm=cbind(mm, binconf(x=mm$sum, n=mm$n)[,2:3])

library(ggplot2)
ggplot(data=mm, aes(x=med, y=BBOcc, color=ShoreDist)) +
  theme_bw() +
  facet_wrap(~GroupId) +
  geom_point() +
  xlab('Julien Day') +
  ylab('Detection Probability') +
  ggtitle('BND Occupancy')

rm(mm)


# 5)  Visualise the daily autocorrelation for all units ################################################################################




acf_val=data.frame(UnitLoc=factor(),
                   lag=numeric(),
                   acf_score=numeric,
                   GroupId=factor(),
                   Shoredist=numeric())

for(ii in 1:30){
  #idx=which(data_detections$UnitLoc==unique(data_detections$UnitLoc)[ii])
  #data_detections$DiffDays[idx]=c(0,data_detections$JulienDay[idx[2:length(idx)]]-data_detections$JulienDay[idx[1:(length(idx)-1)]])
  
  mm=acf(OccTable_daily$BNDTotOffset[OccTable_daily$UnitLoc== levels(OccTable$UnitLoc)[ii]], plot = F)
  
  acf_val=rbind(acf_val, data.frame(acf_score=mm$acf,
                                    lag=mm$lag,
                                    GroupId= strsplit(levels(OccTable$UnitLoc)[ii], split = '_')[[1]][1],
                                    ShoreDist= strsplit(levels(OccTable$UnitLoc)[ii], split = '_')[[1]][2],
                                    UnitLoc=levels(OccTable$UnitLoc)[ii]))
  
  
  rm(mm)
}


# Interleaved histograms

acf_val$jlag=jitter(acf_val$lag)

ggplot(acf_val, aes(x=lag, y=acf_score)) +
  facet_wrap(~GroupId) +
  geom_segment(aes(x = jlag, y = 0, xend = jlag, yend = acf_score, color=ShoreDist)) +
  scale_color_manual(name="Shore Distance",
                     breaks=levels(OccTable$ShoreDist),
                     labels=c("Near", "Mid", "Off"),
                     values=cbbPalette) +
  geom_hline(yintercept = c(0.1, -.1), lty=2, color='red') +
  theme_bw() +
  xlab('Number of Days')+
  ylab('ACF')








# 2) Spatial model #####
OccTable_daily$dateunit=as.factor(paste(OccTable_daily$Date, OccTable_daily$UnitLoc))

mod=gamm(BNDTotOffset~ scale(Slope2) + scale(DistToSalmonRun) + scale(Depth_m) ,
         correlation=corAR1(form = ~1|dateunit),
         family=binomial,
         data=OccTable_daily,
         random=list(UnitLoc=~1))

mod1=gamm(BNDTotOffset~ s(scale(Slope2), k = 3) + s(scale(DistToSalmonRun), k = 3) + s(scale(Depth_m), k = 3) ,
         correlation=corAR1(form = ~1|dateunit),
         family=binomial,
         data=OccTable_daily,
         random=list(UnitLoc=~1))

mod2=gamm(BNDTotOffset~ scale(Slope2) + s(scale(DistToSalmonRun), k = 3) + s(scale(Depth_m), k = 3) ,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))

mod3=gamm(BNDTotOffset~ s(scale(Slope2), k = 3) + scale(DistToSalmonRun) + s(scale(Depth_m), k = 3) ,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))

mod4=gamm(BNDTotOffset~ s(scale(Slope2), k=3) + s(scale(DistToSalmonRun), k=3) + scale(Depth_m) ,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))


mod5=gamm(BNDTotOffset~ s(scale(Slope2), k = 3) + scale(DistToSalmonRun) + scale(Depth_m) ,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))
mod6=gamm(BNDTotOffset~ scale(Slope2) + s(scale(DistToSalmonRun), k = 3) + scale(Depth_m) ,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))
mod7=gamm(BNDTotOffset~ scale(Slope2) + scale(DistToSalmonRun) + s(scale(Depth_m), k = 3) ,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))



modlist_daily=list()
modlist_daily[[1]]=mod
modlist_daily[[2]]=mod1
modlist_daily[[3]]=mod2
modlist_daily[[4]]=mod3
modlist_daily[[5]]=mod4
modlist_daily[[6]]=mod5
modlist_daily[[7]]=mod6
modlist_daily[[8]]=mod7

mm=unlist(lapply(modlist_daily, AIC))

saveRDS(modlist_daily, "Daily_gamm.rds")
modlist=readRDS('C:\\Users\\charlotte\\Documents\\GitHub\\PhD\\Chapter 4\\R Code\\ModelLists.rds')





           
           