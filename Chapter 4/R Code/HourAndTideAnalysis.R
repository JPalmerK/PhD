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


OccTable= read.csv('W:/KJP PHD/4-Bayesian Habitat Use/R Code/OccupancyTable_ThreePdets1.csv')
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

OccTable$TimeofDay='Day'
OccTable$TimeofDay[OccTable$elevation > -5 & OccTable$elevation < 5]='Crepuscular'
OccTable$TimeofDay[OccTable$elevation < -5]='Night'

###############################################
# Function for backwards stepwise selection #
###############################################

SelectModel=function(ModelFull){
  
  # Calculate the QIC of the full model
  fullmodQ=QIC(ModelFull)
  newQIC=0
  terms=attr(ModelFull$terms,"term.labels")
  
  
  while(newQIC != fullmodQ & length(terms)>1){
    
    
    # get all the terms for the full model
    terms <- attr(ModelFull$terms,"term.labels")
    n=length(terms)
    
    newmodel=list()
    newQIC=list()
    
    newmodel[[1]]=ModelFull
    newQIC[[1]]=fullmodQ
    
    # Make n models with selection
    for (ii in 1:n){
      dropvar=terms[ii]
      newTerms <- terms[-match(dropvar,terms)]
      newform <- as.formula(paste(".~.-",dropvar))
      newmodel[[ii+1]] <- update(ModelFull,newform)
      newQIC[[ii+1]] =QIC(newmodel[[ii]])
      
    }
    
    # Get the model with the lowest QIC
    LowestMod=which.min(unlist(newQIC))
    
    if (LowestMod != 1){
      ModelFull=newmodel[[LowestMod]]
      newQIC=min(unlist(newQIC))
    } else {
      ModelFull=ModelFull
      newQIC=min(unlist(newQIC))
    }
    
    
    #end the model selection
    
    
  }
  return(ModelFull)
  
}

######################################################
# Function for walds signficance #
######################################################

DropVarsWalds=function(ModelFull){
  
  # If no terms included return 
  if (length(attr(ModelFull$terms,"term.labels"))<2){
    NewModel='No Covariates to select from'
    
  }else{
    
    
    OldModel=ModelFull
    # Get the anova values
    temp=anova(ModelFull)
    
    # Make n models with selection
    while(length(which(temp$`P(>|Chi|)`>.05))>0 & is.data.frame(temp)){
      
      
      # get the maximum value
      dropvar=rownames(temp)[which.max(temp$`P(>|Chi|)`)]
      
      # new formula for the full model
      newform <- as.formula(paste(".~.-",dropvar))
      
      # new full model
      ModelFull= update(ModelFull,newform) 
      
      # Get the model covariate names
      terms <- attr(ModelFull$terms,"term.labels")
      
      # # Get the anova values
      # temp=anova(ModelFull)
      
      temp=tryCatch({anova(ModelFull)}, error=function(e){e})
      
      
    }
    
    NewModel=ModelFull
  }
  
  return(NewModel)
}


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
  DATA$Observed<-data_sub$BBOcc                                            # the second column reports the observed response (0s and 1s)
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





########################################################
# Fit Models and Predictions #
#######################################################




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

# An empty model is fitted: the binary response "OccAll" is modelled as a function of year, shore dist and GroupID only. These are expressed as B-splines with
# one knot positioned at the average value. The autoregressive  correlation  is used and the block is defined on the basis of the "GroupID and the Date"
# such that within each day the probability of detecting a click depends on whether a click was detected in the previous hour

OccTable_DPD_nocro=OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',]
OccTable_DPD_nocro=droplevels(OccTable_DPD_nocro)
OccTable_DPD_nocro$BNDTotOffset[OccTable_DPD_nocro$BNDTotOffset==1]=0


empty=geeglm(OccAll ~ Year+ ShoreDist + GroupId,
       corstr = 'ar1',
       family = binomial, # leave out constrains
       id=UnitLoc:Date,
       offset = BNDTotOffset,
       data = OccTable_DPD_nocro)




## Test linear or smooth for hour of day and with or without an interaction with GroupID
HoDl=geeglm(OccAll ~Year+ ShoreDist + GroupId + HourAfterPeakSolEle,
            corstr = 'ar1',
            family = binomial, # leave out constrains
            id=UnitLoc:Date,
            offset = BNDTotOffset,
            data = OccTable_DPD_nocro)

HoDs=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
            corstr = 'ar1',
            family = binomial, # leave out constrains
            id=UnitLoc:Date,
            offset = BNDTotOffset,
            data = OccTable_DPD_nocro)

HoDIntS=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))*GroupId,
               corstr = 'ar1',
               family = binomial, # leave out constrains
               id=UnitLoc:Date,
               offset = BNDTotOffset,
               data = OccTable_DPD_nocro)

HoDIntl=geeglm(OccAll ~Year+ ShoreDist + GroupId*HourAfterPeakSolEle,
               corstr = 'ar1',
               family = binomial, # leave out constrains
               id=UnitLoc:Date,
               offset = BNDTotOffset,
               data = OccTable_DPD_nocro)


# Try group ID as a random effect
HoDsREGroup=geeglm(OccAll ~Year+ ShoreDist + GroupId:bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
                   corstr = 'ar1',
                   family = binomial, # leave out constrains
                   id=UnitLoc:Date,
                   offset = BNDTotOffset,
                   data = OccTable_DPD_nocro)






QIC(empty, HoDl, HoDs, HoDIntS, HoDIntl , HoDsREGroup)


# QIC
# empty       14129.63
# HoDl        14140.88
# HoDs        13736.24 # winner
# HoDIntS     13824.43 
# HoDIntl     14169.58
# HoDsREGroup 13800.58



# #############################################################################################################################
# # 2) Determine what form tidal phase should take (linear offset, interaction or smooth)
# #############################################################################################################################

# Test linear or smooth for hour of day and with or without an interaction with GroupID
TideDl=geeglm(OccAll ~Year+ ShoreDist + GroupId + HourAfterHigh,
              corstr = 'ar1',
              family = binomial, # leave out constrains
              id=UnitLoc:Date,
              offset = BNDTotOffset,
              data = OccTable_DPD_nocro)

Tides=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterHigh, knots = mean(HourAfterHigh)),
             corstr = 'ar1',
             family = binomial, # leave out constrains
             id=UnitLoc:Date,
             offset = BNDTotOffset,
             data = OccTable_DPD_nocro)

TidesPh=geeglm(OccAll ~Year+ ShoreDist + GroupId + Phase,
               corstr = 'ar1',
               family = binomial, # leave out constrains
               id=UnitLoc:Date,
               offset = BNDTotOffset,
               data = OccTable_DPD_nocro)

TidesPhInt=geeglm(OccAll ~Year+ ShoreDist*Phase + GroupId ,
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=UnitLoc:Date,
                  offset = BNDTotOffset,
                  data = OccTable_DPD_nocro)

TideIntS=geeglm(OccAll ~Year+ ShoreDist + GroupId + bs(HourAfterHigh, knots = mean(HourAfterHigh))*ShoreDist,
                corstr = 'ar1',
                family = binomial, # leave out constrains
                id=UnitLoc:Date,
                offset = BNDTotOffset,
                data = OccTable_DPD_nocro) #Fail

TideIntl=geeglm(OccAll ~Year+ ShoreDist + GroupId*HourAfterHigh,
                corstr = 'ar1',
                family = binomial, # leave out constrains
                id=UnitLoc:Date,
                offset = BNDTotOffset,
                data = OccTable_DPD_nocro)

# Try group ID as a random effect
TideREGroup=geeglm(OccAll ~Year+ ShoreDist + GroupId:bs(HourAfterHigh, knots = mean(HourAfterHigh)),
                   corstr = 'ar1',
                   family = binomial, # leave out constrains
                   id=UnitLoc:Date,
                   offset = BNDTotOffset,
                   data = OccTable_DPD_nocro)


 QIC(empty, TideDl, Tides, TideIntl, TideREGroup, TideIntS, TidesPh, TidesPhInt) #TideIntS


 # QIC
 # empty       14129.63
 # TideDl      14109.25 #winner
 # Tides       14121.73
 # TideIntl    14120.54
 # TideREGroup 14222.72
 # TideIntS    14130.29
 # TidesPh     14120.46
 # TidesPhInt  14127.63
 
 
# #########################################################################################################################
# ## 3 Use backwards selection to get the model order ##
# #########################################################################################################################

 
 
 ModelFull=geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+GroupId+ShoreDist+Year+HourAfterHigh,
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=UnitLoc:Date,
                  offset = BNDTotOffset,
                  data = OccTable_DPD_nocro)
 
##########################################################################
# Check group id or slope #
########################################################################

 # In sufficient DOF to look at unit loc as a fixed effect 
 ModelFull_slope=geeglm(OccAll ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+Slope2+ShoreDist+Year+HourAfterHigh,
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=UnitLoc:Date,
                  offset = BNDTotOffset,
                  data = OccTable_DPD_nocro)
 
 QIC(ModelFull, ModelFull_slope)
 
# Model with group ID wins
 
 
 
 
 Model_out=SelectModel(ModelFull)
 
 ##############################################################################################
 ## Use repeated Walds tests to select drivers  ####
 #############################################################################################
 
 DropVarsWalds=function(ModelFull){
   
   # If no terms included return 
   if (length(attr(ModelFull$terms,"term.labels"))<2){
     NewModel='No Covariates to select from'
     
   }else{
     
     
     OldModel=ModelFull
     # Get the anova values
     temp=anova(ModelFull)
     
     # Make n models with selection
     while(length(which(temp$`P(>|Chi|)`>.05))>0 & is.data.frame(temp)){
       
       
       # get the maximum value
       dropvar=rownames(temp)[which.max(temp$`P(>|Chi|)`)]
       
       # new formula for the full model
       newform <- as.formula(paste(".~.-",dropvar))
       
       # new full model
       ModelFull= update(ModelFull,newform) 
       
       # Get the model covariate names
       terms <- attr(ModelFull$terms,"term.labels")
       
       # # Get the anova values
       # temp=anova(ModelFull)
       
       temp=tryCatch({anova(ModelFull)}, error=function(e){e})
       
       
     }
     
     NewModel=ModelFull
   }
   
   return(NewModel)
 }
 
 Model_dropped=DropVarsWalds(Model_out)
 

#####################################################################################################################
# Explore AUC values #
#####################################################################################################################



CalcAUC(Model_dropped, data_sub=OccTable_DPD_nocro)


##############################################################################################################
# Create the Partial Plots #
##############################################################################################################



#######################################################
# Julien Date Smoothes #
#######################################################

 # Functions for calculating partial contributions for
 # factors and continuous variables 
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
 

fitdf_Hour=partialDF(Model_dropped, OccTable_DPD_nocro, 'HourAfterPeakSolEle')
fitdf_tide=partialDF(Model_dropped, OccTable_DPD_nocro, 'HourAfterHigh')
fitdf_GroupId=partialdf_factor(Model_dropped, OccTable_DPD_nocro, 'GroupId')
fitdf_ShoreDist=partialdf_factor(Model_dropped, OccTable_DPD_nocro, 'ShoreDist')


# Aggregate the data for plotting
aggdata_hour=data.frame(aggregate(data=OccTable_DPD_nocro, BBOcc~HourAfterPeakSolEle, FUN=mean))
aggdata_tide=data.frame(aggregate(data=OccTable_DPD_nocro, BBOcc~HourAfterHigh, FUN=mean))
aggdata_GroupId=data.frame(aggregate(data=OccTable_DPD_nocro, BBOcc~GroupId, FUN=mean))


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
  # geom_rug(data=subset(data, BBOcc==1), 
  #          aes(x=jitter(HourAfterPeakSolEle, 2.5), y=BBOcc*.1),
  #          sides='t', alpha=.8) +
  # geom_rug(data=subset(data, BBOcc==0), 
  #          aes(x=jitter(HourAfterPeakSolEle, 2.5), y=BBOcc),
  #          sides='b', alpha=.8) +
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


ggplot(data=fitdf_ShoreDist) +
  theme_bw() +
  geom_boxplot(aes(x=ShoreDist, y=inv.logit(vals))) +
  scale_x_discrete(breaks=unique(fitdf_ShoreDist$ShoreDist),
                   labels=c("Near", "Mid", "Off")) +
  xlab("") +
  ylab("")

ggplot(data=fitdf_GroupId) +
  theme_bw() +
  geom_boxplot(aes(x=GroupId, y=inv.logit(vals)))+
  scale_x_discrete(breaks=unique(fitdf_GroupId$GroupId),
                   labels=unique(data$GroupId)) +
  xlab("") +
  ylab("")



#################################################################################################
# Repeat the Analysis for the Cromarty 5 location #
##################################################################################################

# 1 Whether hour of the day, Elevation or Julien Day given elevation is the most important
# 2 Do the same for tide height- also check factor or continuous
# 3 Use backwards selection to get the model order




##################################################################################################
# Exploratory Analysis #
##################################################################################################

Cro_data=OccTable[OccTable$UnitLoc=='Cro_05',]
Cro_data=droplevels(Cro_data)
Cro_data$elevationBin=cut(Cro_data$elevation, breaks = 12, labels = round(seq(min(Cro_data$elevation), max(Cro_data$elevation), length=12)))

aggdata_hour=data.frame(aggregate(data=Cro_data, BBOcc~HourAfterPeakSolEle+Year, FUN=mean))
aggdata_tide=data.frame(aggregate(data=Cro_data, BBOcc~HourAfterHigh+Year, FUN=mean))
aggdata_elevation=data.frame(aggregate(data=Cro_data, BBOcc~elevationBin+Year, FUN=mean))

ggplot(aggdata_hour, aes(HourAfterPeakSolEle, BBOcc, color=Year))+geom_point()
ggplot(aggdata_tide, aes(HourAfterHigh, BBOcc, color=Year))+geom_point()
ggplot(aggdata_elevation, aes(elevationBin, BBOcc, color=Year))+geom_point()




#############################################################################################################################
# 1) Determine what form hour of day should take (linear offset, interaction or smooth)
#############################################################################################################################

# An empty model is fitted: the binary response "OccAll" is modelled as a function of year, shore dist and GroupID only. These are expressed as B-splines with
# one knot positioned at the average value. The autoregressive  correlation  is used and the block is defined on the basis of the "GroupID and the Date"
# such that within each day the probability of detecting a click depends on whether a click was detected in the previous hour

###########################################
# Also check solar elevation/Julien Day 
##########################################

empty=geeglm(OccAll ~ Year,
             corstr = 'ar1',
             family = binomial, # leave out constrains
             id=Date,
             offset = BNDTotOffset,
             data = Cro_data)


## Test linear or smooth for hour of day and with or without an interaction with GroupID
HoDl=geeglm(OccAll ~Year + HourAfterPeakSolEle,
            corstr = 'ar1',
            family = binomial, # leave out constrains
            id=Date,
            offset = BNDTotOffset,
       data = Cro_data)

HoDs=geeglm(OccAll ~Year + bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
            corstr = 'ar1',
            family = binomial, # leave out constrains
            id=Date,
            offset = BNDTotOffset,
            data = Cro_data)


HoDIntL=geeglm(OccAll ~Year*HourAfterPeakSolEle ,
               corstr = 'ar1',
               family = binomial, # leave out constrains
               id=Date,
               offset = BNDTotOffset,
               data = Cro_data)


HoDsREYear=geeglm(OccAll ~Year*bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=Date,
                  offset = BNDTotOffset,
            data = Cro_data)


## Test linear or smooth for hour of day and with or without an interaction with GroupID
Elevationl=geeglm(OccAll ~ Year+ elevation,
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=Date,
                  offset = BNDTotOffset,
            data = Cro_data)

elevations=geeglm(OccAll ~Year + bs(elevation, knots = mean(elevation)),
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=Date,
                  offset = BNDTotOffset,
            data = Cro_data)


elevationIntL=geeglm(OccAll ~Year*elevation ,
                     corstr = 'ar1',
                     family = binomial, # leave out constrains
                     id=Date,
                     offset = BNDTotOffset,
               data = Cro_data)


elevationsREYear=geeglm(OccAll ~Year*bs(elevation, knots = mean(elevation)),
                        corstr = 'ar1',
                        family = binomial, # leave out constrains
                        id=Date,
                        offset = BNDTotOffset,
                  data = Cro_data)


## Test linear for julien day 
JulienDayl=geeglm(OccAll ~ Year+ JulienDay,
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=Date,
                  offset = BNDTotOffset,
                  data = Cro_data)

JulienDays=geeglm(OccAll ~Year + bs(JulienDay, knots = mean(JulienDay)),
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=Date,
                  offset = BNDTotOffset,
                  data = Cro_data)


JulienDayIntL=geeglm(OccAll ~Year*JulienDay ,
                     corstr = 'ar1',
                     family = binomial, # leave out constrains
                     id=Date,
                     offset = BNDTotOffset,
                     data = Cro_data)


JulienDaysREYear=geeglm(OccAll ~Year*bs(JulienDay, knots = mean(JulienDay)) + bs(elevation, knots = mean(elevation)),
                        corstr = 'ar1',
                        family = binomial, # leave out constrains
                        id=Date,
                        offset = BNDTotOffset,
                        data = Cro_data)

# Interactions between julien day and year 

JulienDayIntYear=geeglm(OccAll ~Year:bs(JulienDay, knots = mean(JulienDay)),
                        corstr = 'ar1',
                        family = binomial, # leave out constrains
                        id=Date,
                        offset = BNDTotOffset,
       data = Cro_data)

JulienDayIntElevation=geeglm(OccAll ~elevation*bs(JulienDay, knots = mean(JulienDay)),
                             corstr = 'ar1',
                             family = binomial, # leave out constrains
                             id=Date,
                             offset = BNDTotOffset,
       data = Cro_data)

ElevationIntJulienDay=geeglm(OccAll ~JulienDay*bs(elevation, knots = mean(elevation)),
                             corstr = 'ar1',
                             family = binomial, # leave out constrains
                             id=Date,
                             offset = BNDTotOffset,
       data = Cro_data)

# Combine Elevation and Julien Day of Year
YearIntJdateLelevation= geeglm(OccAll ~Year*bs(JulienDay, knots = mean(JulienDay))+ bs(elevation, knots = mean(elevation)),
                               corstr = 'ar1',
                               family = binomial, # leave out constrains
                               id=Date,
                               offset = BNDTotOffset,
           data = Cro_data)

YearJdateLelevation= geeglm(OccAll ~Year+bs(JulienDay, knots = mean(JulienDay))+ bs(elevation, knots = mean(elevation)),
     corstr = 'ar1',
      family = binomial, # leave out constrains
      id=Date,
      offset = BNDTotOffset,
      data = Cro_data)

YearIntJdateSHour=geeglm(OccAll ~Year*bs(JulienDay, knots = mean(JulienDay))+ bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
        corstr = 'ar1',
        family = binomial, # leave out constrains
        id=Date,
        offset = BNDTotOffset,
        data = Cro_data)

YearJdateSHourS=geeglm(OccAll ~Year+bs(JulienDay, knots = mean(JulienDay)) +
                         bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
           corstr = 'ar1',
           family = binomial, # leave out constrains
           id=Date,
           offset = BNDTotOffset,
           data = Cro_data)

QIC(empty, HoDl, HoDs, HoDIntL, HoDsREYear,
    Elevationl, elevations, elevationIntL, elevationsREYear,
    JulienDayl, JulienDays, JulienDayIntL,JulienDaysREYear,
    JulienDayIntYear, JulienDayIntElevation, ElevationIntJulienDay,
    YearIntJdateLelevation, YearJdateLelevation, YearIntJdateSHour, YearJdateSHourS)



# QIC
# empty                  5899.534
# HoDl                   5896.878
# HoDs                   5889.971
# HoDIntL                5898.250
# HoDsREYear             5896.166
# Elevationl             5877.255
# elevations             5871.729
# elevationIntL          5878.752
# elevationsREYear       5872.244
# JulienDayl             5886.435
# JulienDays             5827.869
# JulienDayIntL          5887.680
# JulienDaysREYear       5795.541
# JulienDayIntYear       5803.119
# JulienDayIntElevation  5828.210
# ElevationIntJulienDay  5865.519
# YearIntJdateLelevation 5795.541
# YearJdateLelevation    5824.088
# YearIntJdateSHour      5789.177
# YearJdateSHourS        5817.988 #winner

# #############################################################################################################################
# # 2) Determine what form hour/solar elevation should take (linear offset, interaction or smooth)
# #############################################################################################################################

# Test linear or smooth for hour of day and with or without an interaction with GroupID
HodDl=geeglm(OccAll ~HourAfterPeakSolEle,
              corstr = 'ar1',
              family = binomial, # leave out constrains
              id=Date,
              offset = BNDTotOffset,
              data = Cro_data) #5906 


HodS=geeglm(OccAll ~ bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
             corstr = 'ar1',
             family = binomial, # leave out constrains
             id=Date,
             offset = BNDTotOffset,
             data =Cro_data) #5903

ElevationL=geeglm(OccAll ~ elevation,
               corstr = 'ar1',
               family = binomial, # leave out constrains
               id=Date,
               offset = BNDTotOffset,
               data = Cro_data) #5883 

ElevationLS=geeglm(OccAll ~ bs(elevation, knots=median(elevation)),
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=Date,
                  offset = BNDTotOffset,
                  data = Cro_data) # 5878

# Solar elevation or julien day?- Julien Day

JdayS=geeglm(OccAll ~ bs(JulienDay, knots=median(JulienDay)),
                   corstr = 'ar1',
                   family = binomial, # leave out constrains
                   id=Date,
                   offset = BNDTotOffset,
                   data = Cro_data) # 5829 

# Julien day AND elevation?

JdayAndElevation=geeglm(OccAll ~ bs(JulienDay, knots=median(JulienDay))+bs(elevation, knots=median(elevation)),
                        corstr = 'ar1',
                        family = binomial, # leave out constrains
                        id=Date,
                        offset = BNDTotOffset,
                        data = Cro_data) #5826


QIC(JdayAndElevation, JdayS, ElevationLS, ElevationL, HodS, HodDl)

# 
# QIC
# JdayAndElevation 5825.630
# JdayS            5829.169
# ElevationLS      5877.808
# ElevationL       5882.547
# HodS             5902.570
# HodDl            5906.209

# #############################################################################################################################
# # 3) Determine what form tidal phase should take (linear offset, interaction or smooth)
# #############################################################################################################################

# Test linear or smooth for hour of day and with or without an interaction with GroupID
TideDl=geeglm(OccAll ~Year+ HourAfterHigh,
              corstr = 'ar1',
              family = binomial, # leave out constrains
              id=Date,
              offset = BNDTotOffset,
              data = Cro_data)

Tides=geeglm(OccAll ~Year+ bs(HourAfterHigh, knots = mean(HourAfterHigh)),
             corstr = 'ar1',
             family = binomial, # leave out constrains
             id=Date,
             offset = BNDTotOffset,
             data =Cro_data)

TidesPh=geeglm(OccAll ~Year+ Phase,
              corstr = 'ar1',
              family = binomial, # leave out constrains
              id=Date,
              offset = BNDTotOffset,
              data = Cro_data)

TidesPhiNT=geeglm(OccAll ~Year*Phase,
               corstr = 'ar1',
               family = binomial, # leave out constrains
               id=Date,
               offset = BNDTotOffset,
               data = Cro_data)




 QIC(empty, TideDl, Tides, TidesPh, TidesPhiNT) 

 # empty      5899.534
 # TideDl     5897.996
 # Tides      5894.389 # Winner
 # TidesPh    5895.536
 # TidesPhiNT 5894.984 
 
#########################################################################################################################
## 3 Use backwards selection to Model Covariates ##
#########################################################################################################################
 
 # Full Model
 fullMod=geeglm(OccAll ~ 
                  bs(elevation, knots=median(elevation))+
                  bs(HourAfterHigh, knots = mean(HourAfterHigh))+
                  bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+
                  Year,
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=Date,
                  offset = BNDTotOffset,
                  data = Cro_data) 
 
 modelout=SelectModel(ModelFull = fullMod)
 
 modelsig=DropVarsWalds(modelout)
 
 CalcAUC(data=Cro_data, mod = modelsig)

 ##############################################################################################################
 # Create the Partial Plots #
 ##############################################################################################################
 
 
 # Aggregate the data for plotting
 aggdata_hour=data.frame(aggregate(data=Cro_data, BBOcc~HourAfterPeakSolEle, FUN=function(x){mean(x)/2}))
 aggdata_tide=data.frame(aggregate(data=Cro_data, BBOcc~HourAfterHigh, FUN=function(x){mean(x)/6}))
 aggdata_tide$Nobs=aggregate(data=Cro_data, UnitLoc~HourAfterHigh, FUN=length)[2]
 aggdata_elevation=data.frame(aggregate(data=Cro_data, BBOcc~elevationBin, FUN=function(x){mean(x)}))

 
 aggdata_GroupId=data.frame(aggregate(data=Cro_data, BBOcc~GroupId, FUN=mean))
 
 
 #######################################################
 # Plot data #
 #######################################################
 

 fitdf_tide=partialDF(modelsig, Cro_data, 'HourAfterHigh')
 fitdf_Year=partialdf_factor(modelsig, Cro_data, 'Year')
# fitdf_jdate=partialDF(modelsig, Cro_data, 'JulienDay')
 fitdf_ele=partialDF(modelsig, Cro_data,'elevation')
 
 ggplot(data=fitdf_tide) +
   theme_bw()+
   scale_colour_manual(values=cbbPalette) +
   geom_line(aes(x=HourAfterHigh, logit(y)), size=1) +  
   geom_ribbon(aes(x=HourAfterHigh, ymin=logit(LCI), ymax=logit(UCI)),alpha=.2,linetype= 'blank') +
   geom_point(data=aggdata_tide, aes(x=HourAfterHigh, y=logit(BBOcc))) +
   xlab('Hour Relative High Tide') +
   ylab('Hour')
 
 # ggplot(data=fitdf_jdate) +
 #   theme_bw()+
 #   scale_colour_manual(values=cbbPalette) +
 #   geom_line(aes(JulienDay, logit(y)), size=1) +  
 #   geom_ribbon(aes(x=JulienDay, ymin=logit(LCI), ymax=logit(UCI)),alpha=.2,linetype= 'blank') +
 #   #geom_point(data=aggdata_tide, aes(x=JulienDay, y=BBOcc)) +
 #   xlab('Hour Relative High Tide') +
 #   ylab('Hour')
 
 ggplot(data=fitdf_ele) +
   theme_bw()+
   scale_colour_manual(values=cbbPalette) +
   geom_line(aes(elevation, y), size=1) +  
   geom_ribbon(aes(x=elevation, ymin=LCI, ymax=UCI),alpha=.2,linetype= 'blank') +
    geom_point(data=aggdata_elevation, aes(x=as.numeric(as.character(aggdata_elevation$elevationBin)), y=BBOcc)) +
   xlab('Solar Elevation Angle Above Horizon') +
   ylab('Detection Probability')
 
 
 
 

  
  
 