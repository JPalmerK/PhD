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
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(mvtnorm)

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

# Set model variables to factors
OccTable$BNDTotOffset[is.nan(OccTable$BNDTotOffset)]=1
OccTable$Year=as.factor(OccTable$Year)
OccTable$ShoreDist=as.factor(OccTable$ShoreDist)

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
  DATA$Observed<-data_sub$OccAll                                                # the second column reports the observed response (0s and 1s)
  DATA$Predicted<-predict(mod,data_sub,type="response")              # the third column reports the predictions
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
# 2 Do the same for tide height- also check factor or continuous
# 3 Use backwards selection to get the model order


#############################################################################################################################
# 1) Determine what form hour of day should take (linear offset, interaction or smooth)
#############################################################################################################################

# # An empty model is fitted: the binary response "OccAll" is modelled as a function of year, shore dist and GroupID only. These are expressed as B-splines with 
# # one knot positioned at the average value. The autoregressive  correlation  is used and the block is defined on the basis of the "GroupID and the Date" 
# # such that within each day the probability of detecting a click depends on whether a click was detected in the previous hour 
# 
# empty=geeglm(OccAll ~ Year+ ShoreDist + GroupId,
#        corstr = 'ar1',
#        family = binomial, # leave out constrains
#        id=GroupId:Date,
#        offset = BNDTotOffset,
#        data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 

# 
# 
# ## Test linear or smooth for hour of day and with or without an interaction with GroupID 
# HoDl=geeglm(OccAll ~Year+ ShoreDist + GroupId + HourAfterPeakSolEle,
#        corstr = 'ar1', 
#        family = binomial, # leave out constrains
#        id=GroupId:Date, 
#        offset = BNDTotOffset, 
#        data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# HoDs=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
#             corstr = 'ar1', 
#             family = binomial, # leave out constrains
#             id=GroupId:Date, 
#             offset = BNDTotOffset, 
#             data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# HoDIntS=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))*GroupId,
#             corstr = 'ar1', 
#             family = binomial, # leave out constrains
#             id=GroupId:Date, 
#             offset = BNDTotOffset, 
#             data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# HoDIntl=geeglm(OccAll ~Year+ ShoreDist + GroupId*HourAfterPeakSolEle,
#                corstr = 'ar1', 
#                family = binomial, # leave out constrains
#                id=GroupId:Date, 
#                offset = BNDTotOffset, 
#                data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# 
# # Try group ID as a random effect
# HoDsREGroup=geeglm(OccAll ~Year+ ShoreDist + GroupId:bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
#             corstr = 'ar1', 
#             family = binomial, # leave out constrains
#             id=GroupId:Date, 
#             offset = BNDTotOffset, 
#             data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# 
# 
# 
# 
# 
# QIC(empty, HoDl, HoDs, HoDIntS, HoDIntl, HoDsREGroup)
# 
# 
# # QIC
# # empty   19746.04
# # HoDl    19756.12
# # HoDs    19215.70 # Winner 
# # HoDIntS 19266.08
# # HoDIntl 19787.06 
# 
# #############################################################################################################################
# # 2) Determine what form tidal phase should take (linear offset, interaction or smooth)
# #############################################################################################################################
# 
## Test linear or smooth for hour of day and with or without an interaction with GroupID
# TideDl=geeglm(OccAll ~Year+ ShoreDist + GroupId + HourAfterHigh,
#               corstr = 'ar1',
#               family = binomial, # leave out constrains
#               id=GroupId:Date,
#               offset = BNDTotOffset,
#               data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# Tides=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterHigh, knots = mean(HourAfterHigh)),
#              corstr = 'ar1',
#              family = binomial, # leave out constrains
#              id=GroupId:Date,
#              offset = BNDTotOffset,
#              data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# TidesPh=geeglm(OccAll ~Year+ ShoreDist + GroupId + Phase,
#               corstr = 'ar1',
#               family = binomial, # leave out constrains
#               id=GroupId:Date,
#               offset = BNDTotOffset,
#               data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# TidesPhInt=geeglm(OccAll ~Year+ ShoreDist*Phase + GroupId ,
#                corstr = 'ar1',
#                family = binomial, # leave out constrains
#                id=GroupId:Date,
#                offset = BNDTotOffset,
#                data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# TideIntS=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterHigh, knots = mean(HourAfterHigh))*HourAfterHigh,
#                 corstr = 'ar1',
#                 family = binomial, # leave out constrains
#                 id=GroupId:Date,
#                 offset = BNDTotOffset,
#                 data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',]) #Fail
# 
# TideIntl=geeglm(OccAll ~Year+ ShoreDist + GroupId*HourAfterHigh,
#                 corstr = 'ar1',
#                 family = binomial, # leave out constrains
#                 id=GroupId:Date,
#                 offset = BNDTotOffset,
#                 data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# # Try group ID as a random effect
# TideREGroup=geeglm(OccAll ~Year+ ShoreDist + GroupId:bs(HourAfterHigh, knots = mean(HourAfterHigh)),
#                    corstr = 'ar1',
#                    family = binomial, # leave out constrains
#                    id=GroupId:Date,
#                    offset = BNDTotOffset,
#                    data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
#  
#  QIC(empty, TideDl, Tides, TideIntl, TideREGroup, TidesPh, TidesPhInt) #TideIntS
# 
# # QIC
# # empty       19746.04
# # TideDl      19708.60 #Winner
# # Tides       19717.72
# # TideIntl    19728.36
# # TideREGroup 19754.79


 
 
# #########################################################################################################################
# ## 3 Use backwards selection to get the model order ##
# #########################################################################################################################
# # 
# QIC(geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+ShoreDist+Year+HourAfterHigh,
#            corstr = 'ar1',
#            family = binomial, # leave out constrains
#            id=GroupId:Date,
#            offset = BNDTotOffset,
#            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])) #19172
# 
# # knock out hour after high
# QIC(geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+ShoreDist+Year,
#            corstr = 'ar1',
#            family = binomial, # leave out constrains
#            id=GroupId:Date,
#            offset = BNDTotOffset,
#            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])) #19214 (3)
# 
# # knock out Year
# QIC(geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+ShoreDist+HourAfterHigh,
#            corstr = 'ar1',
#            family = binomial, # leave out constrains
#            id=GroupId:Date,
#            offset = BNDTotOffset,
#            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])) #19171 (5)
# 
# # Knock out ShoreDist
# QIC(geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+Year+HourAfterHigh,
#            corstr = 'ar1',
#            family = binomial, # leave out constrains
#            id=GroupId:Date,
#            offset = BNDTotOffset,
#            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])) #19176 (4)
# 
# # KNock out GroupId
# QIC(geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+ShoreDist+Year+HourAfterHigh,
#            corstr = 'ar1',
#            family = binomial, # leave out constrains
#            id=GroupId:Date,
#            offset = BNDTotOffset,
#            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',]))  #19215 (2)
# 
# # Knock out hour of day
# QIC(geeglm(OccAll ~GroupId+ShoreDist+Year+HourAfterHigh,
#            corstr = 'ar1',
#            family = binomial, # leave out constrains
#            id=GroupId:Date,
#            offset = BNDTotOffset,
#            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])) # 19708.6 (1)
# 
# 
mod=geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+HourAfterHigh+ShoreDist+Year,
           corstr = 'ar1',
           family = binomial, # leave out constrains
           id=GroupId:Date,
           offset = BNDTotOffset,
           data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# mod1b=geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+Phase+ShoreDist,
#            corstr = 'ar1',
#            family = binomial, # leave out constrains
#            id=GroupId:Date,
#            offset = BNDTotOffset,
#            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# 
# 
# Since yar does't give us much lets compare the models
mod1=geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+HourAfterHigh+ShoreDist,
            corstr = 'ar1',
            family = binomial, # leave out constrains
            id=GroupId:Date,
            offset = BNDTotOffset,
            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])

mod9=geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+bs(HourAfterHigh, degree = 2)+GroupId+ShoreDist+Year,
            corstr = 'ar1',
            family = binomial, # leave out constrains
            id=GroupId:Date,
            offset = BNDTotOffset,
            data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])

 QIC(mod, mod1)

# Year doesn't do too much one way or the other (delta QIC=.07) so knock it out (Also the ROC plot does much better)

summary(mod)
mod=mod1
#####################################################################################################################
# Explore AUC values #
#####################################################################################################################

data=OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',]

CalcAUC(mod, data_sub=data)


##############################################################################################################
# Create the Partial Plots #
##############################################################################################################



#######################################################
# Julien Date Smoothes #
#######################################################
# Get the smoothed index terms using MRSea data

partialDF=function(mod, data, Variable){
  
  coefpos <- c(1, grep(Variable, colnames(model.matrix(mod))))
  xvals <- data[, which(names(data) == Variable)]
  newX <- seq(min(xvals), max(xvals), length = 500)
  eval(parse(text = paste(Variable, "<- newX", 
                          sep = "")))
  response <- rep(1, 500)
  newBasis <- eval(parse(text = labels(terms(mod))[grep(Variable, 
                                                        labels(terms(mod)))]))
  partialfit <- cbind(rep(1, 500), newBasis) %*% coef(mod)[coefpos]
  rcoefs <- NULL
  try(rcoefs <- rmvnorm(1000, coef(mod), summary(mod)$cov.scaled), 
      silent = T)
  if (is.null(rcoefs) || length(which(is.na(rcoefs) == 
                                      T)) > 0) {
    rcoefs <- rmvnorm(1000, coef(mod), as.matrix(nearPD(summary(mod)$cov.scaled)$mat))
  }
  rpreds <- cbind(rep(1, 500), newBasis) %*% t(rcoefs[, 
                                                      coefpos])
  quant.func <- function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }
  cis <- t(apply(rpreds, 1, quant.func))
  
  partialfit <- mod$family$linkinv(partialfit)
  cis <- mod$family$linkinv(cis)
  
  fitdf=data.frame(x=newX, y=partialfit, LCI=cis[,1], UCI=cis[,2]) 
  colnames(fitdf)[1]=Variable
  return(fitdf)
}

partialdf_factor=function(mod, data, variable){
  coeffac <- c(1,grep(variable, colnames(model.matrix(mod))))
  coefradial <- c(grep("LocalRadialFunction", colnames(model.matrix(mod))))
  coefpos <- coeffac[which(is.na(match(coeffac, coefradial)))]
  xvals <- data[, which(names(data) == variable)]
  newX <- sort(unique(xvals))
  newX <- newX[2:length(newX)]
  partialfit <- coef(mod)[c(coefpos)]
  rcoefs <- NULL
  try(rcoefs <- rmvnorm(1000, coef(mod), summary(mod)$cov.scaled), 
      silent = T)
  if (is.null(rcoefs) || length(which(is.na(rcoefs) == T)) > 0) {
    rcoefs <- rmvnorm(1000, coef(mod), as.matrix(nearPD(summary(mod)$cov.scaled)$mat))
  }
  
  if((length(coefpos))>1){
    
    rpreds <- as.data.frame(rcoefs[, c(coefpos)])
    BootstrapCoefs3=data.frame(vals=rpreds[,1])
    BootstrapCoefs3$FactorVariable=as.factor(paste(variable, levels(xvals)[1], sep = ''))
    
    ############################
    # Recompile for plotting #
    #############################
    
    
    for(jj in 2:ncol(rpreds)){
      temp=data.frame(vals=rpreds[,jj])
      temp$FactorVariable=as.factor(colnames(rpreds)[jj])
      BootstrapCoefs3=rbind(BootstrapCoefs3, temp)
      rm(temp)
    }
    
  }else{
    
    rpreds <- rcoefs[,coefpos]
    BootstrapCoefs3=data.frame(vals=rpreds)
    BootstrapCoefs3$FactorVariable=colnames(model.matrix(mod))[coefpos[2:length(coefpos)]]
    
  }
  
  fitdf=BootstrapCoefs3
  colnames(fitdf)[2]=variable
  
  return(fitdf)
  
  
}



fitdf_Hour=partialDF(mod, data, 'HourAfterPeakSolEle')
fitdf_tide=partialDF(mod, data, 'HourAfterHigh')
fitdf_GroupId=partialdf_factor(mod, data, 'GroupId')
fitdf_ShoreDist=partialdf_factor(mod, data, 'ShoreDist')

# Year knocked out
# fitdf_Year=partialdf_factor(mod, data, 'Year')


# Aggregate the data for plotting
aggdata_hour=data.frame(aggregate(data=data, BBOcc~HourAfterPeakSolEle, FUN=mean))
aggdata_tide=data.frame(aggregate(data=data, BBOcc~HourAfterHigh, FUN=mean))
aggdata_GroupId=data.frame(aggregate(data=data, BBOcc~GroupId, FUN=mean))


# 
data_sub=OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',]
data_sub$DecHour=data_sub$MatlabDate[1]-floor(data_sub$MatlabDate[1])*24

# Partial plot for hour after solar noon
ggplot(data=fitdf_Hour) +
  theme_bw()+
  scale_colour_manual(values=cbbPalette) +
  geom_line(aes(HourAfterPeakSolEle, y), size=1) +  
  geom_ribbon(aes(x=HourAfterPeakSolEle, ymin=LCI, ymax=UCI),alpha=.2,linetype= 'blank') +
  geom_point(data=aggdata_hour, aes(x=HourAfterPeakSolEle, y=BBOcc)) +
  geom_rug(data=subset(data_sub, BBOcc==1), 
           aes(x=jitter(HourAfterPeakSolEle,2), y=OccAll*.1),
           sides='t', alpha=.8) +
  geom_rug(data=subset(data_sub, BBOcc==0), 
           aes(x=jitter(HourAfterPeakSolEle,2), y=OccAll),
           sides='b', alpha=.8) +
  xlab('Hour Relative to Solar Noon') +
  ylab('Occupancy Probability') 

# Partial plot for hour after high tide
ggplot(data=fitdf_tide) +
  theme_bw()+
  scale_colour_manual(values=cbbPalette) +
  geom_line(aes(HourAfterHigh, y), size=1) +  
  geom_ribbon(aes(x=HourAfterHigh, ymin=LCI, ymax=UCI),alpha=.2,linetype= 'blank') +
  geom_point(data=aggdata_tide, aes(x=HourAfterHigh, y=BBOcc)) +
  xlab('Hour After High Tide') +
  ylab('Hour') 

# Partial plot for shore dist
ggplot(data=fitdf_ShoreDist) +
  theme_bw() +
  geom_boxplot(aes(x=ShoreDist, y=inv.logit(vals))) +
  scale_x_discrete(breaks=unique(fitdf_ShoreDist$ShoreDist),
                   labels=c("Near", "Mid", "Off")) +
  xlab("") +
  
  ylab("")

# # Partial plot for year
# Year knocked out
# ggplot(data=fitdf_Year) +
#   theme_bw() +
#   geom_boxplot(aes(x=Year, y=inv.logit(vals)))+
#   scale_x_discrete(breaks=unique(fitdf_Year$Year),
#                    labels=c("2013", "2014", "2015")) +
#   xlab("") +
#   ylab("")

ggplot(data=fitdf_GroupId) +
  theme_bw() +
  geom_boxplot(aes(x=GroupId, y=inv.logit(vals)))+
  scale_x_discrete(breaks=unique(fitdf_GroupId$GroupId),
                   labels=unique(data$GroupId)) +
  xlab("") +
  ylab("")





######################################################################################################################
# Get and plot the fit  #
######################################################################################################################
# data=OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',]
# 
# OneYearAggs=data.frame(aggregate(data=subset(data, Year==2013), OccAll~HourAfterPeakSolEle+GroupId+ShoreDist, FUN=mean))
# AggData=cbind(OneYearAggs, predictvcv(mod, newdata = OneYearAggs))
# newdata=subset(OccTable_DPD, Year==2013)
# fit=cbind(newdata,  predictvcv(mod = mod, newdata = newdata))
# 
# 
# 
# # Color Bline Paette
# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# 
# # Plot all the data
# ggplot(data=fit) +
#   theme_bw() +
#   facet_wrap(~GroupId) +
#   scale_colour_manual(values=cbbPalette) +
#   geom_line(aes(Hr, inv.logit(fit), colour=ShoreDist), size=1) +
#   geom_ribbon(aes(x=HourAfterHigh, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
#               alpha=.2,linetype= 'blank') +
# 
#   xlab("") +
#   ylab("")
# 



  
  
 