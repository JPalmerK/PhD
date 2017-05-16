#############
# Setup     #
#############
# This code investigates the various models that look at the different occupancy distributions
# incorporated by the data 
rm(list=ls())
library(boot) # for inv.logit
library(mgcv)
library(ggplot2)
library(lme4)
library(dplyr) # for distinct function 
library(geepack)
library(splines)
library(RColorBrewer)
library(MuMIn) # for QIC

setwd("W:/KJP PHD/4-Bayesian Habitat Use/R Code")


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

meta=read.csv('W:/KJP PHD/CPOD Processing/2013 to 2016 SM deployments.csv')
meta$UnitLoc=factor(meta$UnitLoc, levels=level_names)

meta_sub=subset(meta, select=c('UnitLoc', 'Slope'))
OccTable=merge(OccTable, meta_sub, all.x = TRUE)
rm(meta_sub)


OccTable$IsCroFactor=ifelse(OccTable$UnitLoc=='Cro_05', 'Cro05','Other')

################################################################################
# General Data Prep #
################################################################################

OccTable$GroupId=unlist(strsplit(as.character(OccTable$UnitLoc), split = "_"))[seq(1,(nrow(OccTable)*2)-1,2)]

level_names=c( "Lat", "Hel", "Cro",
               "SpB", "Fra", "Cru",
               "Sto", "Abr", "StA",
               "Stb")

OccTable$GroupId=factor(OccTable$GroupId, levels=level_names)



OccTable$ShoreDist=unlist(strsplit(as.character(OccTable$UnitLoc), split = "_"))[seq(2,(nrow(OccTable)*2),2)]


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
OccTable$SpeciesOffset[OccTable$SpeciesOffset==0] = 1

OccTable$BNDTotOffset=(OccTable$BBOcc*.77+OccTable$FBOcc*.06+OccTable$UNKOcc*.5)/(OccTable$BBOcc+OccTable$FBOcc+OccTable$UNKOcc)


OccTable$BNDTotOffset[is.nan(OccTable$BNDTotOffset)]=1



#####################################################
# Daily Occupancy #
#####################################################



# Add total number of detections 
mm.bbtot=as.data.frame(aggregate(BBOcc~UnitLoc+Date, FUN=sum, data = OccTable))
colnames(mm.bbtot)[3]='BBTot'
mm.fbtot=as.data.frame(aggregate(FBOcc~UnitLoc+Date, FUN=sum, data = OccTable))
colnames(mm.fbtot)[3]='FBTot'
mm.unktot=as.data.frame(aggregate(UNKOcc~UnitLoc+Date, FUN=sum, data = OccTable))
colnames(mm.unktot)[3]='UNKTot'


mm=distinct(OccTable, Date, UnitLoc, JulienDay, GroupId, ShoreDist, Slope, Year,
            Month, IsCroFactor, Hr)
mm.bb=distinct(subset(OccTable, BBOcc>0), Date, BBOcc, UnitLoc)
mm.fb=distinct(subset(OccTable, FBOcc>0), Date, FBOcc, UnitLoc)
mm.unk=distinct(subset(OccTable, UNKOcc>0), Date, UNKOcc, UnitLoc)


OccTable_daily=merge(mm, mm.bb, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.fb, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.unk, by = c('Date', 'UnitLoc'), all.x = TRUE)

OccTable_daily=merge(OccTable_daily, mm.bbtot, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.fbtot, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.unktot, by = c('Date', 'UnitLoc'), all.x = TRUE)



OccTable_daily[is.na(OccTable_daily)] <- 0



OccTable_daily$SpeciesOffset=OccTable_daily$BBOcc+OccTable_daily$FBOcc+OccTable_daily$UNKOcc
OccTable_daily$OccAll=ifelse(OccTable_daily$SpeciesOffset>=1,1,0)




# Species offset for Daily Occupancy
# If two species same day and same unit then uncertain
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset>1]=0.5
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==1 & OccTable_daily$BBOcc==1]=0.77
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==1 & OccTable_daily$FBOcc==1]=0.06
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==1 & OccTable_daily$UNKOcc==1]=0.5
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==0] = 1
# Total offset
OccTable_daily$BNDTotOffset=(OccTable_daily$BBTot*.77+OccTable_daily$FBTot*.06+OccTable_daily$UNKTot*.5)/
  (OccTable_daily$BBTot+OccTable_daily$FBTot+OccTable_daily$UNKTot)

OccTable_daily$TotDet=(OccTable_daily$BBTot+OccTable_daily$FBTot+OccTable_daily$UNKTot)

OccTable_daily$BNDTotOffset[is.na(OccTable_daily$BNDTotOffset)]=1

rm(mm, mm.bb, mm.fb, mm.unk, mm.bbtot, mm.fbtot, mm.unktot)



# Add a dummy variable and remove all days with no detections
mm=aggregate(data=OccTable, OccAll~Date+Year+UnitLoc, FUN = sum)
colnames(mm)[4]='SumHrlyDet'



OccTable=merge(OccTable, mm, all.x=TRUE)
OccTable_DPD=subset(OccTable, SumHrlyDet>0 )

################################################################################
# Function to calculate AUC #
################################################################################

# This function calculates AUC for the model (to access model fit) taken from  
# Pirotta E, Matthiopoulos J, MacKenzie M, Scott-Hayward L, Rendell L Modelling sperm whale habitat preference: a novel approach combining transect and follow data

CalcAUC<-function(mod, data_sub){
  
  pr <- predict(mod,data_sub, type="response")                          # the final model is used to predict the data on the response scale (i.e. a value between 0 and 1)
  pred <- prediction(pr,data_sub$BBOcc)                                    # to specify the vector of predictions (pr) and the vector of labels (i.e. the observed values "Pres")
  perf <- performance(pred, measure="tpr", x.measure="fpr")          # to assess model performance in the form of the true positive rate and the false positive rate
  plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5)) # to plot the ROC curve
  
  
  # Choice of the best cut-off probability
  
  y<-as.data.frame(perf@y.values)
  x<-as.data.frame(perf@x.values)
  fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45° line and the line joining the origin with the point (x;y) on the ROC curve
  L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
  d <- L*sin(fi)                                                     # to calculate the distance between the 45° line and the ROC curve
  # write.table(d,"C:\\distances.txt")                                # to write a table with the computed distances
  
  # The table should then be opened in Microsoft Excel to find the maximum distance with the command "Sort", and the relative position (i.e. the number of the corresponding record)
  # MAX d= 0.1127967 --> position 39
  
  alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
  Best_cutoff=alpha[which.max(unlist(d)),]                           # to identify the alpha value (i.e. the cut-off) that corresponds to the maximum distance between the 45° line and the curve
  
  # Best cutoff:   0.3464173
  # This value can now be used to build the confusion matrix:
  
  DATA<-matrix(0,nrow(data_sub),3)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 919 - the number of rows can be checked with dim(dat)) 
  DATA<-as.data.frame(DATA)
  names(DATA)<-c("plotID","Observed","Predicted")
  DATA$plotID<-1:nrow(data_sub)                                                # the first column is filled with an ID value that is unique for each row
  DATA$Observed<-data_sub$BBOcc                                                # the second column reports the observed response (0s and 1s)
  DATA$Predicted<-predict(modlist[[ii]],data_sub,type="response")              # the third column reports the predictions
  cmx(DATA, threshold = Best_cutoff)                                           # the identified cut-off must be used here
  
  # Area under the Curve 
  auc <- unlist(performance(pred, measure="auc")@y.values)
  
  # Proportion of the presences correctly identified 
  pres=prop.table(cmx(DATA, threshold = Best_cutoff))[1,1]
  
  # Proportion of the absences correctly idenified
  abs=prop.table(cmx(DATA, threshold = Best_cutoff))[2,2]
  
  
  return(c(auc, pres, abs))
}

###########################################################################################################################
# Model Fitting
###########################################################################################################################

# 1 Determine what form hour of day should take (linear offset, interaction or smooth)
# 2 Do the same for tide heigh
# 3 Use backwards selection to get the model order


#############################################################################################################################
# 1) Determine what form hour of day should take (linear offset, interaction or smooth)
#############################################################################################################################

# An empty model is fitted: the binary response "OccAll" is modelled as a function of year, shore dist and GroupID only. These are expressed as B-splines with 
# one knot positioned at the average value. The autoregressive  correlation  is used and the block is defined on the basis of the "GroupID and the Date" 
# such that within each day the probability of detecting a click depends on whether a click was detected in the previous hour 

empty=geeglm(OccAll ~ Year+ ShoreDist + GroupId, 
       corstr = 'ar1', 
       family = binomial, # leave out constrains
       id=GroupId:Date, 
       offset = BNDTotOffset, 
       data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])




## Test linear or smooth for hour of day and with or without an interaction with GroupID 
HoDl=geeglm(OccAll ~Year+ ShoreDist + GroupId + HourAfterPeakSolEle,
       corstr = 'ar1', 
       family = binomial, # leave out constrains
       id=GroupId:Date, 
       offset = BNDTotOffset, 
       data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])

HoDs=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
            corstr = 'ar1', 
            family = binomial, # leave out constrains
            id=GroupId:Date, 
            offset = BNDTotOffset, 
            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])

HoDIntS=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))*GroupId,
            corstr = 'ar1', 
            family = binomial, # leave out constrains
            id=GroupId:Date, 
            offset = BNDTotOffset, 
            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])

HoDIntl=geeglm(OccAll ~Year+ ShoreDist + GroupId*HourAfterPeakSolEle,
               corstr = 'ar1', 
               family = binomial, # leave out constrains
               id=GroupId:Date, 
               offset = BNDTotOffset, 
               data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])


# Try group ID as a random effect
HoDsREGroup=geeglm(OccAll ~Year+ ShoreDist + GroupId:bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
            corstr = 'ar1', 
            family = binomial, # leave out constrains
            id=GroupId:Date, 
            offset = BNDTotOffset, 
            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])






QIC(empty, HoDl, HoDs, HoDIntS, HoDIntl, HoDsREGroup)


# QIC
# empty   19746.04
# HoDl    19756.12
# HoDs    19215.70 # Winner 
# HoDIntS 19266.08
# HoDIntl 19787.06 

#############################################################################################################################
# 2) Determine what form tidal phase should take (linear offset, interaction or smooth)
#############################################################################################################################

## Test linear or smooth for hour of day and with or without an interaction with GroupID 
TideDl=geeglm(OccAll ~Year+ ShoreDist + GroupId + HourAfterHigh,
              corstr = 'ar1', 
              family = binomial, # leave out constrains
              id=GroupId:Date, 
              offset = BNDTotOffset, 
              data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])

Tides=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterHigh, knots = mean(HourAfterHigh)),
             corstr = 'ar1', 
             family = binomial, # leave out constrains
             id=GroupId:Date, 
             offset = BNDTotOffset, 
             data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])

TideIntS=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterHigh, knots = mean(HourAfterHigh))*HourAfterHigh,
                corstr = 'ar1', 
                family = binomial, # leave out constrains
                id=GroupId:Date, 
                offset = BNDTotOffset, 
                data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',]) #Fail

TideIntl=geeglm(OccAll ~Year+ ShoreDist + GroupId*HourAfterHigh,
                corstr = 'ar1', 
                family = binomial, # leave out constrains
                id=GroupId:Date, 
                offset = BNDTotOffset, 
                data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',]) 

# Try group ID as a random effect
TideREGroup=geeglm(OccAll ~Year+ ShoreDist + GroupId:bs(HourAfterHigh, knots = mean(HourAfterHigh)),
                   corstr = 'ar1', 
                   family = binomial, # leave out constrains
                   id=GroupId:Date, 
                   offset = BNDTotOffset, 
                   data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])


QIC(empty, TideDl, Tides, TideIntl, TideREGroup) #TideIntS

# QIC
# empty       19746.04
# TideDl      19708.60 #Winner
# Tides       19717.72
# TideIntl    19728.36
# TideREGroup 19754.79

#########################################################################################################################
## 3 Use backwards selection to get the model order ##
#########################################################################################################################

QIC(geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+ShoreDist+Year+HourAfterHigh,
           corstr = 'ar1', 
           family = binomial, # leave out constrains
           id=GroupId:Date, 
           offset = BNDTotOffset, 
           data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])) #19171.21 

# knock out hour after high
QIC(geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+ShoreDist+Year,
           corstr = 'ar1', 
           family = binomial, # leave out constrains
           id=GroupId:Date, 
           offset = BNDTotOffset, 
           data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])) #19215.7 (2)

# knock out Year
QIC(geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+ShoreDist+HourAfterHigh,
           corstr = 'ar1', 
           family = binomial, # leave out constrains
           id=GroupId:Date, 
           offset = BNDTotOffset, 
           data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])) #19171.14 (5) 

# Knock out ShoreDist
QIC(geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+Year+HourAfterHigh,
           corstr = 'ar1', 
           family = binomial, # leave out constrains
           id=GroupId:Date, 
           offset = BNDTotOffset, 
           data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])) #19174.96 (4)

# KNock out GroupId
QIC(geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+ShoreDist+Year+HourAfterHigh,
           corstr = 'ar1', 
           family = binomial, # leave out constrains
           id=GroupId:Date, 
           offset = BNDTotOffset, 
           data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])) ] #19211.52 (3)

# Knock out hour of day
QIC(geeglm(OccAll ~GroupId+ShoreDist+Year+HourAfterHigh,
           corstr = 'ar1', 
           family = binomial, # leave out constrains
           id=GroupId:Date, 
           offset = BNDTotOffset, 
           data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])) # 19708.6 (1) 


mod=geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+HourAfterHigh+GroupId+ShoreDist+Year,
           corstr = 'ar1', 
           family = binomial, # leave out constrains
           id=GroupId:Date, 
           offset = BNDTotOffset, 
           data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# Since yar does't give us much lets compare the models
mod1=geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+HourAfterHigh+GroupId+ShoreDist,
            corstr = 'ar1', 
            family = binomial, # leave out constrains
            id=GroupId:Date, 
            offset = BNDTotOffset, 
            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])

QIC(mod, mod1)

# Year doesn't do too much one way or the other (delta QIC=.07)
summary(mod)
#####################################################################################################################

