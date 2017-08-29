
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

river_locs=data.frame(Rivername=factor(c("Esk",
                                  "Dee", 
                                  "South Esk", 
                                  "Spey",
                                  "Tay Firth", 
                                  "Tweed",
                                  "Cromarty Firth",
                                  "Ness")),
                      Lat=c("564216",
                            "570841",
                            "564216",
                            "574038",
                            "562702",
                            "554550",
                            "5741454",
                             "573439"),
                      Lon=c("-022641",
                            "-020337",
                            "-022641",
                            "-030555",
                            "-024623",
                            "-015901",
                            "-035924",
                            "-040504"))


river_locs$LatDeg=as.numeric(substr(river_locs$Lat, 1,2)) +
  as.numeric(substr(river_locs$Lat, 3,4))/60 + 
  as.numeric(substr(river_locs$Lat, 5,6))/60/60

river_locs$lonDeg=as.numeric(substr(river_locs$Lon, 1,3)) -
  as.numeric(substr(river_locs$Lon, 4,5))/60 - 
  as.numeric(substr(river_locs$Lon, 6,7))/60/60




# Calculate AUC, %0's id'ed and %1's ided, from Pirotta sperm whale paper- 
# modified for new model
CalcAUC<-function(mod, data_sub, BinaryResponse_var){
  
  pr_dat=as.data.frame(as.numeric(predict(mod, data_sub, type="response")))
  pr_dat$labels=data_sub[, BinaryResponse_var]
  
  colnames(pr_dat)[1]='predictions'
  pred_dat=prediction(pr_dat$predictions, pr_dat$labels)
  perf <- performance(pred_dat, measure="tpr", x.measure="fpr") 
  plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5)) # to plot the ROC curve
  
  y<-as.data.frame(perf@y.values)
  x<-as.data.frame(perf@x.values)
  fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45° line and the line joining the origin with the point (x;y) on the ROC curve
  L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
  d <- L*sin(fi)        
  
  
  alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
  Best_cutoff=alpha[which.max(unlist(d)),] 
  
  
  DATA<-matrix(0,nrow(data_sub),3)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 919 - the number of rows can be checked with dim(dat)) 
  DATA<-as.data.frame(DATA)
  names(DATA)<-c("plotID","Observed","Predicted")
  DATA$plotID<-1:nrow(data_sub)                                                # the first column is filled with an ID value that is unique for each row
  DATA$Observed<-data_sub$BBOcc                                            # the second column reports the observed response (0s and 1s)
  DATA$Predicted<-predict(mod,data_sub,type="response")                 # the third column reports the predictions
  cmx(DATA, threshold = Best_cutoff)   
  
  
  # Area under the Curve 
  auc <- unlist(performance(pred_dat, measure="auc")@y.values)
  
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



# Load meta data and merge slope, depth and distance to nearest river
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


meta=merge(meta[,c('UnitLoc', 'Depth_m')], meta2, by='UnitLoc', all.x = FALSE)
meta_sub=subset(meta, select=c('UnitLoc', 'Slope2', 'DistToSalmonRun','RiverName', 'Depth_m' ))
meta_sub=meta_sub[!duplicated(meta_sub[,1:ncol(meta_sub)]),]
meta_sub$RiverDistFactor='blarg'
# Create factor for distance to each river

for(ii in 1:nrow(river_locs)){
  
  meta_sub$RiverDistFactor[meta_sub$RiverName==unique(river_locs$Rivername[ii])]=
    order(meta_sub$DistToSalmonRun[meta_sub$RiverName==unique(river_locs$Rivername[ii])])
  
}


OccTable=merge(OccTable, meta_sub, by='UnitLoc', all.x = TRUE)
rm(meta, meta2, meta_sub)

OccTable$RiverName=as.factor(OccTable$RiverName)




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

OccTable_daily$dateunit=as.factor(paste(OccTable_daily$Date, OccTable_daily$UnitLoc))
OccTable_daily$ScaleJd=scale(OccTable_daily$JulienDay)
# Add season
OccTable_daily$Season=NA
OccTable_daily$Season[OccTable_daily$Month < 6]='Spring'
OccTable_daily$Season[OccTable_daily$Month > 8]='Autum'
OccTable_daily$Season[which(is.na(OccTable_daily$Season))]='Summer'
OccTable_daily$Season=as.factor(OccTable_daily$Season)

OccTable_daily$CentredMonth=OccTable_daily$Month-median(median(unique(OccTable_daily$Month)))
OccTable_daily$yearseason=as.factor(paste(OccTable_daily$Year, OccTable_daily$Season))
OccTable_daily$SeasonRiver=as.factor(paste(OccTable_daily$Season, OccTable_daily$RiverName))

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


# 5)  Visualise the daily autocorrelation for all units and look into correlation ################################################################################




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

acf_val$jlag=ifelse(acf_val$ShoreDist=='05', acf_val$lag, acf_val$lag+.2)
acf_val$jlag[acf_val$ShoreDist=='15']=acf_val$lag[acf_val$ShoreDist=='15']+.4

  jitter(acf_val$lag)


png(filename = paste('ACF Values.png'),
      units="in", 
       width=6, 
      height=6, 
      pointsize=12,res = 400)
  

ggplot(acf_val, aes(x=lag, y=acf_score)) +
  facet_wrap(~GroupId) +
  geom_segment(aes(x = jlag, y = 0, xend = jlag, yend = acf_score, color=ShoreDist)) +
  scale_color_manual(name="Shore Distance",
                     breaks=levels(OccTable$ShoreDist),
                     labels=c("Near", "Mid", "Off"),
                     values=cbbPalette) +
  geom_hline(yintercept = c(0.1, -.1), lty=2, color='red') +
  theme_bw() +
    xlab('Lag (Days)')+
    ylab('Autocorrelation Score')

dev.off()




cor(OccTable_daily$Depth_m, OccTable_daily$DistToSalmonRun, method='pearson') # correlation -0.4
cor(OccTable_daily$Depth_m, OccTable_daily$Slope2, method='pearson') 
cor(OccTable_daily$DistToSalmonRun, OccTable_daily$Slope2, method='pearson') 




# 2) Spatial model #####


mod=gamm(BNDTotOffset~ scale(Slope2) + scale(DistToSalmonRun) + scale(Depth_m) + RiverName,
         correlation=corAR1(form = ~1|dateunit),
         family=binomial,
         data=OccTable_daily,
         random=list(UnitLoc=~1))

mod1=gamm(BNDTotOffset~ s(scale(Slope2), k = 3) +
            s(scale(DistToSalmonRun), k = 3) + s(scale(Depth_m), k = 3) + RiverName,
         correlation=corAR1(form = ~1|dateunit),
         family=binomial,
         data=OccTable_daily,
         random=list(UnitLoc=~1))

mod2=gamm(BNDTotOffset~ scale(Slope2) + s(scale(DistToSalmonRun), k = 3) + s(scale(Depth_m), k = 3) + RiverName,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))

mod3=gamm(BNDTotOffset~ s(scale(Slope2), k = 3) + scale(DistToSalmonRun) + s(scale(Depth_m), k = 3) + RiverName ,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))

mod4=gamm(BNDTotOffset~ s(scale(Slope2), k=3) + s(scale(DistToSalmonRun), k=3) + scale(Depth_m) + RiverName,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))


mod5=gamm(BNDTotOffset~ s(scale(Slope2), k = 3) + scale(DistToSalmonRun) + scale(Depth_m) + RiverName,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))

mod6=gamm(BNDTotOffset~ scale(Slope2) + s(scale(DistToSalmonRun), k = 3) + scale(Depth_m) + RiverName ,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))

mod7=gamm(BNDTotOffset~ scale(Slope2) + scale(DistToSalmonRun) + s(scale(Depth_m), k = 3) + RiverName ,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))

mod8=gamm(BNDTotOffset~ scale(Slope2) + scale(DistToSalmonRun):RiverName + s(scale(Depth_m), k = 3),
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))

mod9=gamm(BNDTotOffset~ scale(Slope2) + scale(DistToSalmonRun):RiverName,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))

mod10=gamm(BNDTotOffset~ scale(Slope2) + scale(Depth_m):RiverName,
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))




modlist_daily=list()
modlist_daily[[1]]=mod
modlist_daily[[2]]=mod1
modlist_daily[[3]]=mod2
modlist_daily[[4]]=mod2
modlist_daily[[5]]=mod4
modlist_daily[[6]]=mod5
modlist_daily[[7]]=mod6
modlist_daily[[8]]=mod7
modlist_daily[[9]]=mod8
modlist_daily[[10]]=mod9
modlist_daily[[11]]=mod10


saveRDS(modlist_daily, "Daily_gamm.rds")
modlist_daily=readRDS('C:\\Users\\charlotte\\Documents\\GitHub\\PhD\\Chapter 4\\R Code\\Daily_gamm.rds')


# Make Model Selection table
modelselection=data.frame(Formula= unlist(lapply(modlist_daily, FUN = function(x){Reduce(paste, deparse(formula(x)))})),
                          AIC=unlist(lapply(modlist_daily, AIC)))

modelselection$dAIC=modelselection$AIC-min(modelselection$AIC)


# Use AUC to assess model fit
vals=CalcAUC(mod6, OccTable_daily, 'BBOcc')


anova(mod9$lme)
anova(modlist_daily[[10]]$gam)



# Temporal Models ################################################################


smod=gamm(BNDTotOffset~ s(ScaleJd, by=UnitLoc, k=5),
         correlation=corAR1(form = ~1|dateunit),
         family=binomial,
         data=OccTable_daily,
         random=list(UnitLoc=~1))

smod1=gamm(BNDTotOffset~ s(ScaleJd, by=GroupId, k=10),
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1), niterPQL = 10)

smod2=gamm(BNDTotOffset~ s(ScaleJd, by=ShoreDist, k=5),
          correlation=corAR1(form = ~1|dateunit),
          family=binomial,
          data=OccTable_daily,
          random=list(UnitLoc=~1))


modlist_daily_time=list()
modlist_daily_time[[1]]=smod
modlist_daily_time[[2]]=smod1
modlist_daily_time[[3]]=smod2


AIC(smod,  smod2)


saveRDS(modlist_daily_time, "Daily_gamm_time.rds")
modlist=readRDS('C:\\Users\\charlotte\\Documents\\GitHub\\PhD\\Chapter 4\\R Code\\ModelLists.rds')



#############################################################################





# Fit full model Combination of the spatial and temporal models



fullmodel0= gamm(BNDTotOffset ~  scale(DistToSalmonRun):RiverName + RiverName+# spatial
                                     s(ScaleJd, k=3, bs='ts') + Year,                         # temporal
                                 correlation=corAR1(form = ~1|dateunit),
                                 family=binomial,
                                 data=OccTable_daily)


fullmodel1= gamm(BNDTotOffset ~ scale(Slope2) + scale(DistToSalmonRun):RiverName + RiverName+# spatial
                  s(ScaleJd, by=RiverName, bs='ts') + Year,                                # temporal
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily)

fullmodel2= gamm(BNDTotOffset ~ scale(Slope2) + scale(DistToSalmonRun)+RiverName + # spatial
                   s(ScaleJd, by=RiverName, bs='ts') + Year,                                # temporal
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily)

fullmodel3= gamm(BNDTotOffset ~ scale(Slope2) + scale(DistToSalmonRun):RiverName + RiverName + # spatial
                   s(ScaleJd, by=RiverName, bs='ts') ,                                # temporal
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily, method = 'REML')


fullmodel4= gamm(BNDTotOffset ~ scale(Slope2) + scale(DistToSalmonRun)+Year:RiverName + RiverName+# spatial
                  s(ScaleJd, k=3, bs='ts') + Year,                         # temporal
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily)
# 
# fullmodel5= gamm(BNDTotOffset ~ scale(Slope2) + s(scale(DistToSalmonRun), k=3,bs='ts')+ RiverName+# spatial
#                   s(ScaleJd, by=RiverName, k=3),                         # temporal
#                 correlation=corAR1(form = ~1|dateunit),
#                 family=binomial,
#                 data=OccTable_daily)  
# 
# fullmodel6= gamm(BNDTotOffset ~ scale(Slope2) +  RiverName+# spatial
#                   s(ScaleJd, by=RiverName, k=3, bs='ts') ,                         # temporal
#                 correlation=corAR1(form = ~1|dateunit),
#                 family=binomial,
#                 data=OccTable_daily)  

# fullmodel5= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), k=3, , bs='ts')+ RiverName+# spatial
#                   s(ScaleJd, by=RiverName, k=3, bs='ts'),                         # temporal
#                 correlation=corAR1(form = ~1|dateunit),
#                 family=binomial,
#                 data=OccTable_daily)  

fullmodel7= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), k=3, bs='ts')+ ShoreDist+# spatial
                  s(ScaleJd, by=ShoreDist, k=3, bs='ts'),                         # temporal
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily)  


fullmodel15= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), by=Season, k=3, bs='ts')  +  Season,                         # temporal
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily)  

fullmodel21= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), by=Season, k=3, bs='ts')  +  Season + Year,                         # temporal
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily)  

fullmodel17= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun),  by=Season, k=3, bs='ts')  +  Season + scale(Depth_m),                         # temporal
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily)  


fullmodel18= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), by=Season, k=3, bs='ts')  +  Season + s(scale(Depth_m), k=3, bs='ts'),                         # temporal
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily)  

fullmodel20= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), by=Season, k=3, bs='ts')  +  Season + s(scale(Slope2), k=3, bs='ts'),                         # temporal
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily)  

fullmodel19= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), k=3, bs='ts')  +  Season + s(scale(Slope2), k=3, bs='ts'),                         # temporal
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily)  


fullmodel16= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), k=3, bs='ts')  +  Season + scale(Slope2),                         # temporal
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily)  


fullmodel14= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), by=Season, k=3, bs='ts') + scale(Slope2) +  Season,                         # temporal
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily)  

fullmodel8= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), by=Season, k=3, bs='ts') + Year+ scale(Slope2) + Season,                         # temporal
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily)  

fullmodel12= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), by=Season, k=3, bs='ts') + RiverName + scale(Slope2) + Season,                         # temporal
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily)  

fullmodel13= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), by=Season, k=3, bs='ts') + RiverName +  Season,                         # temporal
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily) 



fullmodel9= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), k=3, bs='ts') + SeasonRiver,                         # temporal
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily)  

fullmodel10= gamm(BNDTotOffset ~ s(scale(DistToSalmonRun), k=3, bs='ts')+ RiverName+# spatial
                  s(ScaleJd, by=RiverName, k=3, bs='ts'),                         # temporal
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily)  

fullmodel11= gamm(BNDTotOffset ~ scale(Slope2)  + s(scale(DistToSalmonRun), k=3, bs='ts')+ # spatial
                SeasonRiver,                         # temporal
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily) 

fullmodel22 =gamm(BNDTotOffset ~ 
                te(ScaleJd, scale(DistToSalmonRun), k=3, by=RiverName, bs = 'ts') +
                  RiverName, method="REML",   select=TRUE,                     # temporal
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily)





fullmodel24 =gamm(BNDTotOffset ~ 
                    te(ScaleJd, scale(DistToSalmonRun), k=3, by=RiverName, bs = 'ts') +
                    RiverName, , method="REML",   select=TRUE,                     # temporal
                  correlation=corAR1(form = ~1|dateunit),
                  family=binomial,
                  data=OccTable_daily, 
                  random=list(UnitLoc=~1))

fullmodel23 =gamm(BNDTotOffset ~ 
                  s(ScaleJd, k=3, bs='ts', by=SeasonRiver) +
                 SeasonRiver,                         # temporal
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily)


fullmodel24 =gamm(BNDTotOffset ~ 
                    s(ScaleJd, k=5, bs='ts', by=RiverName) +
                    RiverName + Slope2,                         # temporal
                  correlation=corAR1(form = ~1|dateunit),
                  family=binomial,
                  data=OccTable_daily)

# 
# fullmodel24 =gamm(BNDTotOffset ~ 
#                   ti(ScaleJd, scale(DistToSalmonRun), k=3, bs='ts', by=RiverName) + #can't converge
#                   RiverName,                         # temporal
#                 correlation=corAR1(form = ~1|dateunit),
#                 family=binomial,
#                 data=OccTable_daily)





# 
# fullmode9= gamm(BNDTotOffset ~ UnitLoc +# spatial
#                   s(ScaleJd, by=UnitLoc, k=5),                         # Cant Converge
#                 correlation=corAR1(form = ~1|dateunit),
#                 family=binomial,
#                 data=OccTable_daily)  


modlist=list(fullmodel0, fullmodel1, fullmodel2, fullmodel3, fullmodel4,  
                           fullmodel8, fullmodel9, fullmodel10, fullmodel11, fullmodel12, fullmodel13, fullmodel14,fullmodel15,
                          fullmodel16, fullmodel17, fullmodel18, fullmodel19, fullmodel20, fullmodel21,
             fullmodel22,fullmodel23)

mm=data.frame(AIC(fullmodel0, fullmodel1, fullmodel2, fullmodel3, fullmodel4,  
                  fullmodel8, fullmodel9, fullmodel10, fullmodel11, fullmodel12, fullmodel13, fullmodel14,fullmodel15,
                  fullmodel16, fullmodel17, fullmodel18, fullmodel19, fullmodel20, fullmodel21,
                  fullmodel22,fullmodel23))

mm[order(mm$AIC),]


ModTable_AIC=data.frame(Formula=unlist(lapply(modlist, function(x){Reduce(paste,deparse(formula(x)))})))
ModTable_AIC$AIC=as.numeric(lapply(modlist, function(x){round(AIC(x),digits = 0)}))
ModTable_AIC$DeltaAIC=ModTable_AIC$AIC-min(ModTable_AIC$AIC)
ModTable_AIC$ModNames=rownames(mm)
ModTable_AIC$AUC=0 # Area under the Curve 
ModTable_AIC$Pres=0 # Proportion of the presences correctly identified 
ModTable_AIC$Abs=0 # Proportion of the absences correctly idenified




for(ii in 1:length(modlist)){
  AUC_vals=CalcAUC(modlist[[ii]], OccTable_daily, 'BBOcc')
  ModTable_AIC$AUC[ii]=round(AUC_vals[1], 2)
  ModTable_AIC$Pres[ii]=round(AUC_vals[2], 2)
  ModTable_AIC$Abs[ii]=round(AUC_vals[3], 2)
  
}




# Get STrath... coefficients 
coefs <- data.frame(coef(summary(fullmode15$lme)))
coefs$df.Satt <- coef(summary(fullmode15$lme))[, 5]
coefs$df.Satt <- coef(summary(fullmode15))[, 3]


## Look at Model Performance #################################################################

# Table to store model performance #

ModelTable=data.frame(GroupId=levels(OccTable_daily$Unit))
ModelTable$ModelFormula='none'
ModelTable$AIC=0
ModelTable$AUC=0 # Area under the Curve 
ModelTable$Pres=0 # Proportion of the presences correctly identified 
ModelTable$Abs=0 # Proportion of the absences correctly idenified



# Fill in the Model Table
for(ii in 1:130){
  
  data_sub=subset(OccTable_daily, UnitLoc==unique(OccTable_daily$UnitLoc)[ii])
  data_sub$ShoreDist=factor(data_sub$ShoreDist,
                            levels=c('05', '10', '15'))
  data_sub <- droplevels(data_sub)
  
  
  mod=fullmode15
  ModelTable$ModelFormula[ii]=Reduce(paste, deparse(formula(mod)))
  ModelTable$AIC[ii]= AIC(mod)
  
  AUC_vals= CalcAUC(mod, data_sub, "BBOcc")
  ModelTable$AUC[ii]=AUC_vals[1]
  ModelTable$Pres[ii]=AUC_vals[2]
  ModelTable$Abs[ii]=AUC_vals[3]
  rm(mod)
}






plot(fullmodel22$gam, scheme=2, asp=1, select=1, trans=inv.logit)




# generate data for plotting

RiverName=factor()
dist=numeric()
JulienDay=numeric()

# Add distance ranges 
for (ii in 1:length(unique(OccTable_daily$RiverName))){
  
  
  data_sub=OccTable_daily[(OccTable_daily$RiverName==unique(OccTable_daily$RiverName)[ii]),]
  
  
  # River Names
  RiverName=c(RiverName, 
              rep(as.character(data_sub$RiverName[1]), 200))
  
  # Distance Ranges
  range_vals=range(data_sub$DistToSalmonRun)
  
  # Monitoring Ranges
  Jrange=range(data_sub$JulienDay)
  
  dist=c(dist,
         seq(from=range_vals[1], to=range_vals[2], length.out = 100))
  
  JulienDay=c(JulienDay,
              seq(from=Jrange[1], to=Jrange[2], length.out = 100))
}



fitdata=data.frame(JulienDay=JulienDay,
                  DistToSalmonRun=dist,
                  RiverName=RiverName)

fitdata$ScaledDist=scale(fitdata$DistToSalmonRun)
fitdata$ScaleJd=scale(fitdata$JulienDay)

  # # predictions for simulated data 
  # preds=predict(fullmodel22, fitdata , type="response", se=TRUE)
  # fitdata_preds=cbind(fitdata, preds)



library(akima)
library(colorspace)
my.heat.colors <- function(x) { (diverge_hcl(x)) }

# Create the matricies and plot 
for(ii in 1:length(unique(OccTable_daily$RiverName))){
  

  data_sub=subset(OccTable_daily, RiverName==unique(OccTable_daily$RiverName)[ii])
  data_sub=droplevels(data_sub)
  
  myjd <- matrix(data=seq(from=min(data_sub$ScaleJd), 
                          to=max(data_sub$ScaleJd), length=100),
                 nrow=100, ncol=100)    
  
  mydist <- t(matrix(data=seq(from=min(data_sub$DistToSalmonRun),
                              to=max(data_sub$DistToSalmonRun), length=100),
                     nrow=100, ncol=100)) 
  
  fitdata <- data.frame(DistToSalmonRun=as.vector(mydist),
                        ScaleJd=as.vector(myjd), 
                        RiverName=as.character(data_sub$RiverName[1]) )
  
  mypredict.exit <- predict(fullmodel22, fitdata, type="response")  
  
  mypredict.exit <- matrix(mypredict.exit, nrow=c(100,100))

  mypredict.exit=mypredict.exit[!duplicated(mypredict.exit)]
  
  filled.contour(
                    z=mypredict.exit)
  ,  nlevels=15, color=my.heat.colors)
  ,
                 main=fitdata_temp$RiverName[1], 
                 cex.main=2, las=2,
                 ylab="Distance from Feature (km)", 
                 xlab="Julian Day of Year",
                 plot.axes={
                   axis(1, at=c(-2,-1,0,1,2), pos=0, 
                        labels=as.numeric(round(quantile(fitdata_temp$JulienDay,probs = seq(0,1,length.out = 5)))),
                        las=0, col="black")
                   axis(2, at=c(0,1,2,3,4,5), pos=-2, labels=c(0,1,2,3,4,5), las=0, col="black")
                   rect(-0.83, 1.75, 0.83, 3.45, border="black", lty="dashed", lwd=2)
                 },)
  
  # rug(x=unique(OccTable_daily$DistToSalmonRun[OccTable_daily$RiverName==fitdata_preds$RiverName[1]])/1000,
  #     side = 2)
  # 
  
}




  
mypredict.exit <- matrix(mypredict.exit, nrow=c(100,100))








filled.contour()
fitdata_temp=fitdata_temp[order(fitdata_temp$x),]
contour(x = fitdata_temp$x, y = fitdata_temp$y, fitdata_temp$fit)


library(RColorBrewer)
#png(file="AstrosExitVelo.png", width=600, height=675)






filled.contour(z=matrix(mypredict.exit, nrow=c(100,100)))



,
               
               zlim=c(85,92), nlevels=15,
               color=colorRampPalette(rev(brewer.pal(11, "RdYlBu"))), 
               main="Houston Astros 2015 Exit Velocity", cex.main=2, xlab="Horizontal Location (ft., Umpire's View)", ylab="Vertical Location (ft.)",
               
               plot.axes={
                 axis(1, at=c(-2,-1,0,1,2), pos=0, labels=c(-2,-1,0,1,2), las=0, col="black")
                 axis(2, at=c(0,1,2,3,4,5), pos=-2, labels=c(0,1,2,3,4,5), las=0, col="black")
                 rect(-0.83, 1.75, 0.83, 3.45, border="black", lty="dashed", lwd=2)
               },
               key.axes={
                 ylim=c(0,1.0)
                 axis(4, at=c(85,86,87,88,89,90,91,92), labels=c(85,86,87,88,89,90,91,92), pos=1, las=0, col="black")
               })
text(1.4, 2.5, "Exit Velocity", cex=1.1, srt=90)
#dev.off()





