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
                            "573957",
                            "562702",
                            "554550",
                            "5741454",
                            "573439"),
                      Lon=c("-022641",
                            "-020337",
                            "-022641",
                            "-033819",
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




OccTable= read.csv('W:/KJP PHD/4-Bayesian Habitat Use/R Code/OccupancyTable_ThreePdets.csv')
level_names_unit=c( "Lat_05", "Lat_10", "Lat_15",
               "Hel_05", "Hel_10", "Hel_15",
               "Cro_05", "Cro_10", "Cro_15",
               "SpB_05", "SpB_10", "SpB_15",
               "Fra_05", "Fra_10", "Fra_15",
               "Cru_05", "Cru_10", "Cru_15",
               "Sto_05", "Sto_10", "Sto_15",
               "Abr_05", "Abr_10", "Abr_15",
               "StA_05", "StA_10", "StA_15",
               "Stb_05", "Stb_10", "Stb_15")

OccTable$UnitLoc=factor(OccTable$UnitLoc, levels=level_names_unit)
OccTable$DateUnitloc=as.factor(paste(OccTable$Date, OccTable$UnitLoc))

# Distance from shore factor
OccTable$ShoreDist=as.character(OccTable$ShoreDist)
OccTable$ShoreDist[OccTable$ShoreDist=="5"]='05'
OccTable$ShoreDist=as.factor(OccTable$ShoreDist)

# Group ID factor
level_names=c( "Lat", "Hel", "Cro",
               "SpB", "Fra", "Cru",
               "Sto", "Abr", "StA",
               "Stb")
OccTable$GroupId=unlist(strsplit(as.character(OccTable$UnitLoc), split = "_"))[seq(1,(nrow(OccTable)*2)-1,2)]
OccTable$GroupId=factor(OccTable$GroupId, levels=level_names)


# YEAR
OccTable$Year=as.factor(OccTable$Year)

# Load meta data and merge slope, depth and distance to nearest river
meta=read.csv('W:/KJP PHD/Deployment Information/CPODs for Kaitlin.csv')
meta$UnitLoc=factor(meta$UnitLoc, levels=level_names_unit)


meta2=read.csv('W:\\KJP PHD\\Deployment Information\\SlopeAndAspect.csv')
meta2$UnitLoc=factor(meta2$UnitLoc, levels=level_names_unit)
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


OccTable=merge(OccTable, meta_sub, by='UnitLoc', all.x = TRUE)
rm(meta, meta2, meta_sub)

# Create Grouping Factors 
id1=as.numeric(as.factor(OccTable$Date)) ## used to make id4 but not used in correlation
id2=as.numeric(as.factor(OccTable$UnitLoc)) ## used to make id4 but not used in correlation
id3=as.numeric(as.factor(OccTable$HourAfterPeakSolEle))
id4=as.numeric(paste(id1, id2, sep=""))  

OccTable$id3=id3
OccTable$id4=id4



# Calculate hourly  detection rate

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




# Add a dummy variable and remove all days with no detections
mm=aggregate(data=OccTable, OccAll~Date+Year+UnitLoc, FUN = sum)
colnames(mm)[4]='SumHrlyDet'



OccTable=merge(OccTable, mm, all.x=TRUE)
OccTable_DPD=subset(OccTable, SumHrlyDet>0 )



OccTable_DPD$IsCro=as.factor(ifelse(OccTable_DPD$UnitLoc=='Cro_05', 'Cro_05', 'Other'))
OccTable_DPD_nocro=OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',]
OccTable_DPD_nocro=droplevels(OccTable_DPD_nocro)


OccTable_Cro=OccTable[OccTable$UnitLoc == 'Cro_05',]
OccTable_Cro=droplevels(OccTable_Cro)





################################################################################################################



OccTable$IsCro=as.factor(ifelse(OccTable$UnitLoc=='Cro_05', 'Cro_05', 'Other'))
OccTable_NoCro=OccTable[OccTable$IsCro=='Other',]

FullModel_NoCro=gamm(BNDTotOffset ~ShoreDist+GroupId+s(HourAfterHigh, bs='cc', k=10, by=GroupId) +
                 te(HourAfterPeakSolEle, scale(DistToSalmonRun), bs='cc', k=10, by=GroupId),
               correlation=corAR1(form= ~id3|id4),
               data=OccTable_NoCro,
               family=binomial,
               random=list(UnitLoc=~1))








FullModel=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc', k=10)+ s(HourAfterPeakSolEle, bs='cc', k=10, by=IsCro)+IsCro,
           correlation=corAR1(form= ~id3|id4),
             data=OccTable_DPD_nocro,
             family=binomial,
             random=list(UnitLoc=~1))

FullModel_1=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ s(HourAfterPeakSolEle, bs='cc', k=10, by=GroupId)+GroupId,
               correlation=corAR1(form=~id3|id4),
               data=OccTable_DPD_nocro,
               family=binomial,
               random=list(UnitLoc=~1))

FullModel_2=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ s(HourAfterPeakSolEle, bs='cc', k=10, by=ShoreDist)+ShoreDist,
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))

FullModel_3=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ te(HourAfterPeakSolEle, DistToSalmonRun, bs='cc', k=10),
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))

FullModel_4=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ te(HourAfterPeakSolEle, Depth_m, bs='cc', k=10),
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))
FullModel_5=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ te(HourAfterPeakSolEle, scale(DistToSalmonRun), bs='cc', k=10),
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))

FullModel_6=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ s(HourAfterPeakSolEle, bs='cc', k=10, by=GroupId) 
                 + scale(Depth_m) + GroupId,
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))


FullModel_8=gamm(BNDTotOffset ~te(Z,Depth_m, k=10)+ s(HourAfterPeakSolEle, bs='cc', k=10, by=GroupId),+GroupId
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))


FullModel_6=gamm(BNDTotOffset ~scale(Depth_m) + te(HourAfterPeakSolEle, HourAfterHigh, bs='cc', k=10),
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))

FullModel_7=gamm(BNDTotOffset ~ HourAfterHigh + s(HourAfterPeakSolEle, bs='cc', k=10, by=GroupId) + Depth_m + DistToSalmonRun,
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))

FullModel_9=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ s(HourAfterPeakSolEle, bs='cc', k=10, by=GroupId) + Depth_m +DistToSalmonRun,
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))

FullModel_10=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ s(HourAfterPeakSolEle, bs='cc', k=10, by=GroupId),
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))



mm=data.frame(AIC(FullModel, FullModel_1, FullModel_2, FullModel_3, FullModel_4, 
                  FullModel_5, FullModel_6, FullModel_7, FullModel_8, FullModel_9, 
                  FullModel_10))
mm=mm[order(mm$AIC),]

modlist=list(FullModel, FullModel_1, FullModel_2, FullModel_3, FullModel_4, 
              FullModel_5, FullModel_6, FullModel_7, FullModel_8, FullModel_9, 
              FullModel_10)

# Make Model Selection table
modelselection=data.frame(Formula= unlist(lapply(modlist, FUN = function(x){Reduce(paste, deparse(formula(x)))})),
                          AIC=unlist(lapply(modlist, AIC)))
  
modelselection$dAIC=modelselection$AIC-min(modelselection$AIC)
  




# Make model performance table 
ModelTable=data.frame(GroupId=levels(OccTable_DPD_nocro$GroupId))
ModelTable$ModelFormula='none'
ModelTable$AIC=0
ModelTable$AUC=0 # Area under the Curve 
ModelTable$Pres=0 # Proportion of the presences correctly identified 
ModelTable$Abs=0 # Proportion of the absences correctly idenified



# Fill in the Model Table
for(ii in 1:10){
  
  data_sub=subset(OccTable_DPD_nocro, GroupId==unique(OccTable_DPD_nocro$GroupId)[ii])
  data_sub$ShoreDist=factor(data_sub$ShoreDist,
                            levels=c('05', '10', '15'))
  data_sub <- droplevels(data_sub)
  
  
  mod=FullModel_10
  ModelTable$ModelFormula[ii]=Reduce(paste, deparse(formula((mod))))
  ModelTable$AIC[ii]= AIC(FullModel_10)
  
  AUC_vals=CalcAUC(FullModel_10, OccTable_DPD_nocro, 'BBOcc')
  ModelTable$AUC[ii]=AUC_vals[1]
  ModelTable$Pres[ii]=AUC_vals[2]
  ModelTable$Abs[ii]=AUC_vals[3]
  rm(mod)
}






