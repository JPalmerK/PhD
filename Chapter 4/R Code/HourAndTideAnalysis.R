# This code undertakes model selection for the sub-daily occupancy analysis

# 1) Load packages and functions ################################################################################
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




# Backwards stepwise selection
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

# Walds signficance #
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

# Function to calculate AUC #
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



# 2) Load and prep the data ################################################################

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

# meta=read.csv('W:/KJP PHD/CPOD Processing/2013 to 2016 SM deployments.csv')
# meta$UnitLoc=factor(meta$UnitLoc, levels=level_names)
# 
# meta_sub=subset(meta, select=c('UnitLoc', 'Slope'))
# OccTable=merge(OccTable, meta_sub, all.x = TRUE)
# rm(meta_sub)

# Add the new slope values and update

meta2=read.csv('W:\\KJP PHD\\Deployment Information\\SlopeAndAspect.csv')
meta2$UnitLoc=factor(meta2$UnitLoc, levels=level_names)



meta_sub=subset(meta2, select=c('UnitLoc', 'Slope2'))

OccTable=merge(OccTable, meta_sub, all.x = TRUE)
rm( meta2, meta_sub)


OccTable$IsCroFactor=ifelse(OccTable$UnitLoc=='Cro_05', 'Cro05','Other')

OccTable$TimeofDay='Day'
OccTable$TimeofDay[OccTable$elevation > -5 & OccTable$elevation < 5]='Crepuscular'
OccTable$TimeofDay[OccTable$elevation < -5]='Night'


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


OccTable$BNDTotOffset=(OccTable$BBOcc*.77+OccTable$FBOcc*.06+OccTable$UNKOcc*.5)/(OccTable$BBOcc+OccTable$FBOcc+OccTable$UNKOcc)

# Set model variables to factors
OccTable$BNDTotOffset[is.nan(OccTable$BNDTotOffset)]=0
OccTable$Year=as.factor(OccTable$Year)
OccTable$ShoreDist=as.factor(OccTable$ShoreDist)


rm(mm, mm.bb, mm.fb, mm.unk, mm.bbtot, mm.fbtot, mm.unktot)



# Add a dummy variable and remove all days with no detections
mm=aggregate(data=OccTable, OccAll~Date+Year+UnitLoc, FUN = sum)
colnames(mm)[4]='SumHrlyDet'



OccTable=merge(OccTable, mm, all.x=TRUE)
OccTable_DPD=subset(OccTable, SumHrlyDet>0 )
# OccTable_DPD$BNDTotOffset=(ifelse(OccTable_DPD$OccAll==0,
#                                               0,
#                                               (1-OccTable_DPD$BNDTotOffset)))


# 3) Build the models for non-Cromarty ################################################################

# An empty model is fitted: the binary response "OccAll" is modelled as a function of year, shore dist and GroupID only. These are expressed as B-splines with
# one knot positioned at the average value. The autoregressive  correlation  is used and the block is defined on the basis of the "GroupID and the Date"
# such that within each day the probability of detecting a click depends on whether a click was detected in the previous hour





OccTable_DPD_nocro=OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',]
OccTable_DPD_nocro=droplevels(OccTable_DPD_nocro)

Null=geeglm(BNDTotOffset ~ Year ,
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=UnitLoc:Date,
                  #offset = BNDTotOffset,
                  data = OccTable_DPD_nocro)

empty_Unit=geeglm(BNDTotOffset ~ Year + UnitLoc,
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=UnitLoc:Date,
                  #offset = BNDTotOffset,
                  data = OccTable_DPD_nocro)

empty_Slope=geeglm(BNDTotOffset ~ Year+ Slope2,
                   corstr = 'ar1',
                   family = binomial, # leave out constrains
                   id=UnitLoc:Date,
                   #offset = BNDTotOffset,
                   data = OccTable_DPD_nocro)

empty_SlopeBS=geeglm(BNDTotOffset ~ Year+ bs(Slope2, knots=mean(Slope2)),
                      corstr = 'ar1',
                      family = binomial, # leave out constrains
                      id=UnitLoc:Date,
                      #offset = BNDTotOffset,
                      data = OccTable_DPD_nocro)

empty_GrupIdShoreDist=geeglm(BNDTotOffset ~ Year+ GroupId + ShoreDist,
                    corstr = 'ar1',
                    family = binomial, # leave out constrains
                   id=UnitLoc:Date,
               #offset = BNDTotOffset,
               data = OccTable_DPD_nocro)



QIC(Null, empty_Unit, empty_SlopeBS, empty_Slope, empty_GrupIdShoreDist)

# # 
# QIC
# empty_Unit            8030.293
# empty_SlopeBS         8159.221
# empty_Slope           8166.644
# empty_GrupIdShoreDist 8024.853 #winner


New_Null= empty_GrupIdShoreDist

## Test linear or smooth for hour of day 

HoDl=geeglm(BNDTotOffset ~ Year+  HourAfterPeakSolEle,
            corstr = 'ar1',
            family = binomial, # leave out constrains
            id=UnitLoc:Date,
            #offset = BNDTotOffset,
            data = OccTable_DPD_nocro)

HoDs=geeglm(BNDTotOffset ~ Year+ bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
            corstr = 'ar1',
            family = binomial, # leave out constrains
            id=UnitLoc:Date,
            #offset = BNDTotOffset,
            data = OccTable_DPD_nocro)

SolarEles=geeglm(BNDTotOffset ~ Year+ bs(elevation, knots = mean(elevation)),
            corstr = 'ar1',
            family = binomial, # leave out constrains
            id=UnitLoc:Date,
            #offset = BNDTotOffset,
            data = OccTable_DPD_nocro)

SolarElel=geeglm(BNDTotOffset ~ Year+ elevation,
                 corstr = 'ar1',
                 family = binomial, # leave out constrains
                 id=UnitLoc:Date,
                 #offset = BNDTotOffset,
                 data = OccTable_DPD_nocro)




QIC( HoDl, HoDs, SolarElel, SolarEles)

# QIC
# empty 13504.77
# HoDl  13514.98
# HoDs  13108.20 # winner

# Determine what form tidal phase should take (linear offset, interaction or smooth)
# Test linear or smooth for hour of day and with or without an interaction with GroupID
Tidel=geeglm(BNDTotOffset ~Year + HourAfterHigh,
             corstr = 'ar1',
             family = binomial, # leave out constrains
             id=Date,
             #offset = BNDTotOffset,
             data = OccTable_DPD_nocro)

Tides=geeglm(BNDTotOffset ~Year + bs(HourAfterHigh, knots = mean(HourAfterHigh)),
             corstr = 'ar1',
             family = binomial, # leave out constrains
             id=Date,
             #offset = BNDTotOffset,
             data = OccTable_DPD_nocro)

TideHeithl=geeglm(BNDTotOffset ~Year + Z,
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=Date,
                  #offset = BNDTotOffset,
                  data = OccTable_DPD_nocro)

TideHeights=geeglm(BNDTotOffset ~Year+  bs(Z, knots = mean(Z)),
                   corstr = 'ar1',
                   family = binomial, # leave out constrains
                   id=Date,
                   #offset = BNDTotOffset,
                   data = OccTable_DPD_nocro)

TidesPh=geeglm(BNDTotOffset ~Year+  Phase,
               corstr = 'ar1',
               family = binomial, # leave out constrains
               id=Date,
               #offset = BNDTotOffset,
               data = OccTable_DPD_nocro)


QIC(Null, TidesPh, TideHeights, TideHeithl, Tides, Tidel) #TideIntS


# QIC
# Null       8185.373
# TidesPh     8197.
# TideHeights 8211.708
# TideHeithl  8193.622
# Tides       8200.669
# Tidel       8186.807 # winner

## 3 Use backwards QIC selection to get the model fit  
ModelFull=geeglm(BNDTotOffset ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+
                   UnitLoc + 
                   Year,
                 corstr = 'ar1',
                 family = binomial, # leave out constrains
                 id=Date,
                 #offset = BNDTotOffset,
                 data = OccTable_DPD_nocro)


# Select best model 
Model_out=SelectModel(ModelFull)

# Use repeated Walds tests to select drivers  
Model_dropped=DropVarsWalds(Model_out)

# Get AUC score
CalcAUC(Model_dropped, data_sub=OccTable_DPD_nocro)


# 4) Plot the non-Cromarty Models ################################################################

fitdf_Hour=partialDF(Model_dropped, OccTable_DPD_nocro, 'HourAfterPeakSolEle')
fitdf_tide=partialDF(Model_dropped, OccTable_DPD_nocro, 'HourAfterHigh')
fitdf_UnitLoc=partialdf_factor(Model_dropped, OccTable_DPD_nocro, 'UnitLoc')

fitdf_UnitLoc$ShoreDist=substr(fitdf_UnitLoc$UnitLoc, 12,13)
fitdf_UnitLoc$GroupId=substr(fitdf_UnitLoc$UnitLoc, 8,10)


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


ggplot(data=fitdf_UnitLoc) +
  theme_bw() +
  geom_boxplot(aes(x=UnitLoc, y=inv.logit(vals))) +
  scale_x_discrete(breaks=unique(fitdf_UnitLoc$ShoreDist),
                   labels=c("Near", "Mid", "Off")) +
  xlab("") +
  ylab("")



# 5) Build a model from Cromarty 05 ################################################################

# 1 Whether hour of the day, Elevation or Julien Day given elevation is the most important
# 2 Do the same for tide height- also check factor or continuous
# 3 Use backwards selection to get the model order




# Exploratory Analysis #

Cro_data=OccTable_DPD[OccTable_DPD$UnitLoc=='Cro_05',]
Cro_data=droplevels(Cro_data)
Cro_data$elevationBin=cut(Cro_data$elevation, breaks = 12, labels = round(seq(min(Cro_data$elevation), max(Cro_data$elevation), length=12)))

aggdata_hour=data.frame(aggregate(data=Cro_data, BNDTotOffset~HourAfterPeakSolEle+Year, FUN=mean))
aggdata_tide=data.frame(aggregate(data=Cro_data, BNDTotOffset~HourAfterHigh+Year, FUN=mean))
aggdata_elevation=data.frame(aggregate(data=Cro_data, BNDTotOffset~elevationBin+Year, FUN=mean))

ggplot(aggdata_hour, aes(HourAfterPeakSolEle, BNDTotOffset, color=Year))+geom_point()
ggplot(aggdata_tide, aes(HourAfterHigh, BNDTotOffset, color=Year))+geom_point()
ggplot(aggdata_elevation, aes(elevationBin, BNDTotOffset, color=Year))+geom_point()



# 1) Determine what forms the covariates should take

# An empty model is fitted: the binary response "OccAll" is modelled as a function of year, shore dist and GroupID only. These are expressed as B-splines with
# one knot positioned at the average value. The autoregressive  correlation  is used and the block is defined on the basis of the "GroupID and the Date"
# such that within each day the probability of detecting a click depends on whether a click was detected in the previous hour



## Test linear or smooth for hour of day 
empty=geeglm(BNDTotOffset ~Year,
             corstr = 'ar1',
             family = binomial, # leave out constrains
             id=Date,
             #offset = BNDTotOffset,
             data = Cro_data)

HoDl=geeglm(BNDTotOffset ~Year + HourAfterPeakSolEle,
            corstr = 'ar1',
            family = binomial, # leave out constrains
            id=Date,
            #offset = BNDTotOffset,
            data = Cro_data)

HoDs=geeglm(BNDTotOffset ~Year + bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle)),
            corstr = 'ar1',
            family = binomial, # leave out constrains
            id=Date,
            #offset = BNDTotOffset,
            data = Cro_data)

SolarEles=geeglm(BNDTotOffset ~ Year+ bs(elevation, knots = mean(elevation)),
                 corstr = 'ar1',
                 family = binomial, # leave out constrains
                 id=Date,
                 #offset = BNDTotOffset,
                 data = Cro_data)

SolarElel=geeglm(BNDTotOffset ~ Year+ elevation,
                 corstr = 'ar1',
                 family = binomial, # leave out constrains
                 id=Date,
                 #offset = BNDTotOffset,
                 data = Cro_data)




QIC(empty, HoDl, HoDs, SolarElel, SolarEles)

# QIC
# empty 13504.766
# HoDl   4606.496
# HoDs   4602.937 #winner


# 2) Determine what form tidal phase should take (linear offset, interaction or smooth)

# Test linear or smooth for hour of day and with or without an interaction with GroupID
Tidel=geeglm(BNDTotOffset ~Year+ HourAfterHigh,
             corstr = 'ar1',
             family = binomial, # leave out constrains
             id=Date,
             #offset = BNDTotOffset,
             data = Cro_data)

Tides=geeglm(BNDTotOffset ~Year + bs(HourAfterHigh, knots = mean(HourAfterHigh)),
             corstr = 'ar1',
             family = binomial, # leave out constrains
             id=Date,
             #offset = BNDTotOffset,
             data = Cro_data)

TideHeithl=geeglm(BNDTotOffset ~Year+ Z ,
                  corstr = 'ar1',
                  family = binomial, # leave out constrains
                  id=Date,
                  #offset = BNDTotOffset,
                  data = Cro_data)

TideHeights=geeglm(BNDTotOffset ~Year + bs(Z, knots = mean(Z)),
                   corstr = 'ar1',
                   family = binomial, # leave out constrains
                   id=Date,
                   #offset = BNDTotOffset,
                   data = Cro_data)

TidesPh=geeglm(BNDTotOffset ~Year+ Phase,
               corstr = 'ar1',
               family = binomial, # leave out constrains
               id=Date,
               #offset = BNDTotOffset,
               data = Cro_data)


QIC(empty, TidesPh, TideHeights, TideHeithl, Tides, Tidel) #TideIntS

# QIC
# empty        3770.003
# TidesPh      3770.883
# TideHeights  3775.503
# TideHeithl   3770.282 
# Tides        3771.453
# Tidel        3769.510 # winner (barely)




## 3 Use backwards QIC selection to get the model fit  ##

ModelFull=geeglm(BNDTotOffset ~bs(HourAfterPeakSolEle, knots = mean(HourAfterPeakSolEle))+
                   Z+
                   Year,
                 corstr = 'ar1',
                 family = binomial, # leave out constrains
                 id=Date,
                 #offset = BNDTotOffset,
                 data = Cro_data)

Cro_Model_out=SelectModel(ModelFull)

Cro_Model_sig=DropVarsWalds(Cro_Model_out)

CalcAUC(Cro_Model_sig, data_sub =Cro_data )

# 6) Plot the Cromarty Model ################################################################

# Aggregate the data for plotting
aggdata_hour=data.frame(aggregate(data=Cro_data, BNDTotOffset~HourAfterPeakSolEle, 
                                  FUN=function(x){mean(x)/(length(formula(Cro_Model_sig)))}))


aggdata_elevation=data.frame(aggregate(data=Cro_data, BNDTotOffset~HourAfterPeakSolEle, 
                                       FUN=function(x){mean(x)/(length(formula(Cro_Model_sig)))}))

aggdata_GroupId=data.frame(aggregate(data=Cro_data, BNDTotOffset~GroupId,
                                     FUN=function(x){mean(x)/(length(formula(Cro_Model_sig)))}))

Cro_data$TideCut=cut(x = Cro_data$Z, breaks =  quantile(Cro_data$Z, seq(0,1, by = .1)))
aggdata_tide=data.frame(aggregate(data=Cro_data, BNDTotOffset~TideCut, 
                                  FUN=function(x){mean(x)/(length(formula(Cro_Model_sig)))}))
aggdata_tide$Z=aggregate(data=Cro_data, Z~TideCut, FUN=median)[,2]




fitdf_tide=partialDF(Cro_Model_sig, Cro_data, 'Z')
fitdf_Year=partialdf_factor(Cro_Model_sig, Cro_data, 'Year')
fitdf_ele=partialDF(Cro_Model_sig, Cro_data,'HourAfterPeakSolEle')
fitdf_year=partialdf_factor(Cro_Model_sig, Cro_data, 'Year')

ggplot(data=fitdf_ele) +
  theme_bw()+
  scale_colour_manual(values=cbbPalette) +
  geom_line(aes(x=HourAfterPeakSolEle, y=y), size=1) +  
  geom_ribbon(aes(x=HourAfterPeakSolEle, ymin=LCI, ymax=UCI),alpha=.2,linetype= 'blank') +
  geom_point(data=aggdata_hour,
             aes(x=HourAfterPeakSolEle, y=BNDTotOffset)) +
  xlab('Hour Relative to Solar Noon') +
  ylab('Detection Probability')


ggplot(data=fitdf_tide) +
  theme_bw()+
  scale_colour_manual(values=cbbPalette) +
  geom_line(aes(x=Z, (y)), size=1) +  
  geom_ribbon(aes(x=Z, ymin=(LCI), ymax=(UCI)),alpha=.2,linetype= 'blank') +
  geom_point(data=aggdata_tide, aes(x=Z, y=(BNDTotOffset))) +
  xlab('Hour Relative High Tide') +
  ylab('Hour')
# 
# ggplot(data=fitdf_year) +
#   theme_bw()+
#   scale_colour_manual(values=cbbPalette) +
#   geom_boxplot(aes(x=Year, y=inv.logit(vals)), size=1)
# 




