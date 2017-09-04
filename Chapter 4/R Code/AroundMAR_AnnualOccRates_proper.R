#########################################################################
# Descriptive Statistics For Occupancy #
#########################################################################

# This code reads in the hourly occupancy table for dolphin detections and models
# daily occupancy using GEEGLM's

# 1)  Load packages and declare local functions      ##############

# This code investigates the various models that look at the different occupancy distributions
# incorporated by the data 
# Outdated !!! Use Mapping
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


river_locs=data.frame(Rivername=c("Esk", "Dee", "South Esk", "Spey", "Tay", "Thurso", "Tweed"),
                 Lat=c("555638", "570723", "564259", "574014", "562243", "583536", "554512"),
                 Lon=c("-045656", "-035257"," -032659", "-045432", "-043518", "-033056", "-035341"))

river_locs$LatDeg=as.numeric(substr(river_locs$Lat, 1,2)) +
  as.numeric(substr(river_locs$Lat, 3,4))/60 + 
  as.numeric(substr(river_locs$Lat, 5,6))/60/60

river_locs$lonDeg=as.numeric(substr(river_locs$Lon, 1,3)) -
  as.numeric(substr(river_locs$Lon, 4,5))/60 - 
  as.numeric(substr(river_locs$Lon, 6,7))/60/60




### Functions  

# These function calculates do backwards stepwise QIC to select the best model,
# then repeated Walds tests to retain meaningful variables and CalcAUC to caluclate
# the area under the ROC curve for each fitted model. All code based on
# Pirotta E, Matthiopoulos J, MacKenzie M, Scott-Hayward L, Rendell L Modelling sperm whale habitat preference: a novel approach combining transect and follow data

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

# Function to catch the warnings, used with make.model to nix
# crappy models
withWarnings <-function (expr) { 
  warnings <- character() 
  retval <- withCallingHandlers(expr, warning = function(ex) { 
    warnings <<- c(warnings, conditionMessage(ex)) 
    invokeRestart("muffleWarning") 
  }) 
  list(Value = retval, Warnings = warnings) 
} 

# Function to do QIC selection for models with interactions only
Make.model <- function (data_sub) {
  modformula=list()
  modformula[[1]]=as.formula(BNDTotOffset~Year+ShoreDist*bs(JulienDay, knots = mean(JulienDay)))
  modformula[[2]]=as.formula(BNDTotOffset~ShoreDist+Year*bs(JulienDay, knots = mean(JulienDay)))
  modformula[[3]]=as.formula(BNDTotOffset~ShoreDist+Year+bs(JulienDay, knots = mean(JulienDay)))
  
  QIC_val=c(Inf,Inf,Inf)
  
  for(ii in 1:3){
    
    mod <- try(geeglm(modformula[[ii]], 
                      family = binomial, 
                      data   = data_sub,
                      weights = rep(1, nrow(data_sub)),
                      id     = UnitLoc,
                      corstr = 'ar1'), silent = T)
    
    
    if (class(mod) == "try-error") {
      cat("Insufficient data for full model")
      
      
    }else{
      
      mod_warns=withWarnings(geeglm(modformula[[ii]], 
                                    family = binomial, 
                                    data   = data_sub,
                                    weights = rep(1, nrow(data_sub)),
                                    id     = UnitLoc,
                                    corstr = 'ar1'))['Warnings']
      
      # If there are two or more model warnings then it's likely the fitted probabilities didn't 
      # converge " fitted probabilities numerically 0 or 1 occurred ", in whcih case nix the model
      if(length(mod_warns[[1]])<2){ QIC_val[ii]=QIC(mod)}
      
      
      
    }
  }  
  
  # Get the model with the lowest QIC
  modformula[[which.min(QIC_val)]]
  
  out.mod=geeglm(modformula[[which.min(QIC_val)]], 
                 family = binomial, 
                 data   = data_sub,
                 weights = rep(1, nrow(data_sub)),
                 id     = UnitLoc,
                 corstr = 'ar1')
  warnings()
  return(out.mod)
}


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



# functions to produce partial plots (largely based on MRSea, Scott-Hawyard)
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
    fitdf_year=data.frame(vals=rpreds[,1])
    fitdf_year$FactorVariable=as.factor(paste(variable, levels(xvals)[1], sep = ''))
    
    ############################
    # Recompile for plotting #
    #############################
    
    
    for(jj in 2:ncol(rpreds)){
      temp=data.frame(vals=rpreds[,jj])
      temp$FactorVariable=as.factor(colnames(rpreds)[jj])
      fitdf_year=rbind(fitdf_year, temp)
      rm(temp)
    }
    
  }else{
    
    rpreds <- rcoefs[,coefpos]
    fitdf_year=data.frame(vals=rpreds)
    fitdf_year$FactorVariable=colnames(model.matrix(mod))[coefpos[2:length(coefpos)]]
    
  }
  
  fitdf=fitdf_year
  colnames(fitdf)[2]=variable
  
  return(fitdf)
  
  
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
meta=read.csv('W:/KJP PHD/CPOD Processing/2013 to 2016 SM deployments.csv')
meta$UnitLoc=factor(meta$UnitLoc, levels=level_names)



meta2=read.csv('W:\\KJP PHD\\Deployment Information\\SlopeAndAspect.csv')
meta2$UnitLoc=factor(meta2$UnitLoc, levels=level_names)

meta=merge(meta, meta2, by='UnitLoc', all.x = TRUE)
colnames(meta)[26]='Slope'

meta_sub=subset(meta, select=c('UnitLoc', 'Slope2'))

OccTable=merge(OccTable, meta_sub, all.x = TRUE)


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


mm=distinct(OccTable, Date, UnitLoc, JulienDay, GroupId, ShoreDist, Slope2, Year, Month)
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


# 4)  Make Chapter tables for Overall Occupancy Rates ##################################################################################

# Table for proportion of days occupied
Basic_table=aggregate(data=subset(OccTable_daily, Year==2013), OccAll~UnitLoc, FUN=sum)
temp1=aggregate(data=subset(OccTable_daily, Year==2014), OccAll~UnitLoc, FUN=sum)
temp2=aggregate(data=subset(OccTable_daily, Year==2015), OccAll~UnitLoc, FUN=sum)

temp3=aggregate(data=subset(OccTable_daily, Year==2013), Date~UnitLoc, FUN=length)
temp4=aggregate(data=subset(OccTable_daily, Year==2014), Date~UnitLoc, FUN=length)
temp5=aggregate(data=subset(OccTable_daily, Year==2015), Date~UnitLoc, FUN=length)


temp6=aggregate(data=subset(OccTable_daily, Year==2013), OccAll~UnitLoc, FUN=mean)
temp7=aggregate(data=subset(OccTable_daily, Year==2014), OccAll~UnitLoc, FUN=mean)
temp8=aggregate(data=subset(OccTable_daily, Year==2015), OccAll~UnitLoc, FUN=mean)


Basic_table=merge(Basic_table, temp1, all = TRUE, by='UnitLoc')
Basic_table=merge(Basic_table, temp2, all = TRUE, by='UnitLoc')

Basic_table=merge(Basic_table, temp3, all = TRUE, by='UnitLoc')
Basic_table=merge(Basic_table, temp4, all = TRUE, by='UnitLoc')
Basic_table=merge(Basic_table, temp5, all = TRUE, by='UnitLoc')

Basic_table=merge(Basic_table, temp6, all = TRUE, by='UnitLoc')
Basic_table=merge(Basic_table, temp7, all = TRUE, by='UnitLoc')
Basic_table=merge(Basic_table, temp8, all = TRUE, by='UnitLoc')

# Reorder

Basic_table=Basic_table[, c(1,2,5,8, 3,6,9, 4,7,10)]

rm(temp1, temp2, temp3, temp4,temp5,temp6, temp7, temp8)

colnames(Basic_table)=c('UnitLoc', 'NDetections2013', 'NDays2013', 'PropOccupied2013'
                        , 'NDetections2014', 'NDays2014', 'PropOccupied2014'
                        , 'NDetections2015', 'NDays2015', 'PropOccupied2015')
Basic_table[is.na(Basic_table)]=0

library(binom)
Basic_table$BinConf2013L=binconf(x=Basic_table$NDetections2013, n = Basic_table$NDays2013)[,2]
Basic_table$BinConf2014L=binconf(x=Basic_table$NDetections2014, n = Basic_table$NDays2014)[,2]
Basic_table$BinConf2015L=binconf(x=Basic_table$NDetections2015, n = Basic_table$NDays2015)[,2]

Basic_table$BinConf2013U=binconf(x=Basic_table$NDetections2013, n = Basic_table$NDays2013)[,3]
Basic_table$BinConf2014U=binconf(x=Basic_table$NDetections2014, n = Basic_table$NDays2014)[,3]
Basic_table$BinConf2015U=binconf(x=Basic_table$NDetections2015, n = Basic_table$NDays2015)[,3]



# write.csv(Basic_table, 'W:/KJP PHD/4-Bayesian Habitat Use/Figures/Basic Ocuppancy 2013-2015.csv')

# Table for BND occupancy  #

Basic_table_bb=aggregate(data=subset(OccTable_daily, Year==2013), BBOcc~UnitLoc, FUN=sum)
temp1=aggregate(data=subset(OccTable_daily, Year==2014), BBOcc~UnitLoc, FUN=sum)
temp2=aggregate(data=subset(OccTable_daily, Year==2015), BBOcc~UnitLoc, FUN=sum)

temp3=aggregate(data=subset(OccTable_daily, Year==2013), Date~UnitLoc, FUN=length)
temp4=aggregate(data=subset(OccTable_daily, Year==2014), Date~UnitLoc, FUN=length)
temp5=aggregate(data=subset(OccTable_daily, Year==2015), Date~UnitLoc, FUN=length)


temp6=aggregate(data=subset(OccTable_daily, Year==2013), BBOcc~UnitLoc, FUN=mean)
temp7=aggregate(data=subset(OccTable_daily, Year==2014), BBOcc~UnitLoc, FUN=mean)
temp8=aggregate(data=subset(OccTable_daily, Year==2015), BBOcc~UnitLoc, FUN=mean)


Basic_table_bb=merge(Basic_table_bb, temp1, all = TRUE, by='UnitLoc')
Basic_table_bb=merge(Basic_table_bb, temp2, all = TRUE, by='UnitLoc')

Basic_table_bb=merge(Basic_table_bb, temp3, all = TRUE, by='UnitLoc')
Basic_table_bb=merge(Basic_table_bb, temp4, all = TRUE, by='UnitLoc')
Basic_table_bb=merge(Basic_table_bb, temp5, all = TRUE, by='UnitLoc')

Basic_table_bb=merge(Basic_table_bb, temp6, all = TRUE, by='UnitLoc')
Basic_table_bb=merge(Basic_table_bb, temp7, all = TRUE, by='UnitLoc')
Basic_table_bb=merge(Basic_table_bb, temp8, all = TRUE, by='UnitLoc')

# Reorder (this could be done programatically but there ar only 10 groups, so it's fine)
Basic_table_bb=Basic_table_bb[, c(1,2,5,8, 3,6,9, 4,7,10)]

rm(temp1, temp2, temp3, temp4,temp5,temp6, temp7, temp8)

colnames(Basic_table_bb)=c('UnitLoc', 'NDetections2013', 'NDays2013', 'PropOccupied2013'
                           , 'NDetections2014', 'NDays2014', 'PropOccupied2014'
                           , 'NDetections2015', 'NDays2015', 'PropOccupied2015')

Basic_table_bb[is.na(Basic_table_bb)]=0

library(binom)
Basic_table_bb$BinConf2013L=binconf(x=Basic_table_bb$NDetections2013, n = Basic_table_bb$NDays2013)[,2]
Basic_table_bb$BinConf2014L=binconf(x=Basic_table_bb$NDetections2014, n = Basic_table_bb$NDays2014)[,2]
Basic_table_bb$BinConf2015L=binconf(x=Basic_table_bb$NDetections2015, n = Basic_table_bb$NDays2015)[,2]

Basic_table_bb$BinConf2013U=binconf(x=Basic_table_bb$NDetections2013, n = Basic_table_bb$NDays2013)[,3]
Basic_table_bb$BinConf2014U=binconf(x=Basic_table_bb$NDetections2014, n = Basic_table_bb$NDays2014)[,3]
Basic_table_bb$BinConf2015U=binconf(x=Basic_table_bb$NDetections2015, n = Basic_table_bb$NDays2015)[,3]







# 5) Final pre-processing-only include >2 units with some detections in modelling process ###

## Finaslf ##
OccTable_daily$YrUnitLoc=paste(OccTable_daily$UnitLoc, OccTable_daily$Year, sep='')

OccTable_daily_wDetections= OccTable_daily[OccTable_daily$YrUnitLoc %in%
                                             names(which(tapply(OccTable_daily$OccAll,OccTable_daily$YrUnitLoc,sum)>2)),]
OccTable_daily_wDetections=droplevels(OccTable_daily_wDetections)




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

  





# X) Fit models to each deployment gorups, as before (make sure to load predictvcv first!) #######################################

  # list to store the models #
  modlist=list()
  
  
  ModelTable=data.frame(GroupId=levels(OccTable_daily_wDetections$GroupId))
  ModelTable$Data2013=0
  ModelTable$Data2014=0
  ModelTable$Data2015=0
  ModelTable$Nunits2013=0
  ModelTable$Nunits2014=0
  ModelTable$Nunits2015=0
  ModelTable$DeltaQIC=0
  
  modlist=list()
  par(ask=F)
  # Build Models for 10 Variables
  for(ii in 1:10){
  
    
    data_sub=subset(OccTable_daily_wDetections, GroupId==unique(OccTable$GroupId)[ii])
    data_sub$ShoreDist=factor(data_sub$ShoreDist, levels=c('05', '10', '15'))
    data_sub <- droplevels(data_sub)
    
    
    ModelTable$Nunits2013[ii]=length(unique(data_sub$UnitLoc[data_sub$Year==2013]))
    ModelTable$Nunits2014[ii]=length(unique(data_sub$UnitLoc[data_sub$Year==2014]))
    ModelTable$Nunits2015[ii]=length(unique(data_sub$UnitLoc[data_sub$Year==2015]))
    
    newdat=subset(data_sub, select=c('JulienDay', 'ShoreDist', 'GroupId', 'UnitLoc', 'OccAll', 'Year','BNDTotOffset'))
    
    tryCatch({
      ModelTable$Data2013[ii]<-as.character(Reduce(paste, unique(data_sub$UnitLoc[data_sub$Year==2013])))}, error=function(e){})
    
    tryCatch({
      ModelTable$Data2014[ii]<-as.character(Reduce(paste, unique(data_sub$UnitLoc[data_sub$Year==2014])))}, error=function(e){})
    
    
    tryCatch({
      ModelTable$Data2015[ii]<-as.character(Reduce(paste, unique(data_sub$UnitLoc[data_sub$Year==2015])))}, error=function(e){})
    
    
    # Determine whether linear model or spline for Julien Day #
  
    modformula=list()
    modformula[[1]]=as.formula(BNDTotOffset~Year+ShoreDist*bs(JulienDay, knots = mean(JulienDay)))
    modformula[[2]]=as.formula(BNDTotOffset~ShoreDist+Year*bs(JulienDay, knots = mean(JulienDay)))
    modformula[[3]]=as.formula(BNDTotOffset~ShoreDist+Year+bs(JulienDay, knots = mean(JulienDay)))
    modformula[[4]]=as.formula(BNDTotOffset~ShoreDist+Year+JulienDay)
    
    QIC_val=c(Inf,Inf,Inf)
    print(ii)
    for(jj in 1:4){
      mod=try(geeglm(modformula[[jj]], 
                     family = binomial, 
                     data   = data_sub,
                     weights = rep(1, nrow(data_sub)),
                     id     = UnitLoc,
                     corstr = 'ar1'), silent = T)
      
      
      
      if (class(mod) == "try-error"){QIC_val[jj]=Inf} else {
        temp=withWarnings(geeglm(modformula[[jj]], 
                                family = binomial, 
                                data   = data_sub,
                                weights = rep(1, nrow(data_sub)),
                                id     = UnitLoc,
                                corstr = 'ar1'))
        
        if(length(temp$Warnings)<2){QIC_val[jj]=QIC(temp$Value)}
        
      }
      print(jj)
    }
    
    finalmodel=geeglm(modformula[[which.min(QIC_val)]], 
                      family = binomial, 
                      data   = data_sub,
                      weights = rep(1, nrow(data_sub)),
                      id     = UnitLoc,
                      corstr = 'ar1')
    
    
    ModelTable$DeltaQIC[ii]=sort(QIC_val)[2]-sort(QIC_val)[1]
  
    
    
    # Document the formula
    ModelTable$ModelFormula[ii]=Reduce(paste, deparse(formula(finalmodel)))
    
  
      
      terms_orig=attr(finalmodel$terms,"term.labels")
      terms_idx= grep('JulienDay', terms_orig)
      terms=terms_orig
      
      
      # Expand terms
      
      if(length(grep('JulienDay', terms_orig))>0){jday=seq(min(newdat$JulienDay), max(newdat$JulienDay))
      }else{
        jday=median(newdat$JulienDay)}
      
      if(length(grep('Year', terms_orig))>0){Year=data_sub$Year[1]}else{
        Year=unique(data_sub$Year)}
      
      if(length(grep('ShoreDist', terms_orig))>0){ShoreDist=data_sub$ShoreDist[1]}else{
        ShoreDist=unique(data_sub$ShoreDist)}
      
      
      # Keep the final model 
      modlist[[ii]]=finalmodel
      
      
      # Report AUC values
      AUCvals=CalcAUC(finalmodel, data_sub = data_sub)
      ModelTable$AUC[ii]=round(AUCvals[1],2)
      ModelTable$Pres[ii]=round(AUCvals[2],2)
      ModelTable$Abs[ii]=round(AUCvals[3],2)  
      

      # Make Predictions
      if(ii==1){
        fit=cbind(newdat,  predictvcv(finalmodel))
        
      }else {
        fit=rbind(fit, cbind(newdat, predictvcv(finalmodel)) )
            }
    }  
    

  # Look at model performance at all 30 locations
  Modperformance=data.frame(GroupId=character(length = length(unique(OccTable_daily_wDetections$UnitLoc))),
                            UnitLoc=character(length = length(unique(OccTable_daily_wDetections$UnitLoc))),
                            AUC_unit=numeric(length = length(unique(OccTable_daily_wDetections$UnitLoc))),
                            Pres_unit=numeric(length = length(unique(OccTable_daily_wDetections$UnitLoc))),
                            Abs_unit=numeric(length = length(unique(OccTable_daily_wDetections$UnitLoc))))
  
  levels(Modperformance$GroupId)= levels(OccTable_daily_wDetections$GroupId)
  levels(Modperformance$UnitLoc)= levels(OccTable_daily_wDetections$UnitLoc)
  
  count=1
  for(ii in 1:10){
    
    data_sub=subset(OccTable_daily_wDetections, GroupId==unique(OccTable$GroupId)[ii])
    
    n=length(unique(data_sub$UnitLoc))
    for(jj in 1:n){
      
      
      subsub=subset(data_sub, UnitLoc==unique(data_sub$UnitLoc)[jj])
      mm=CalcAUC(modlist[[ii]], data_sub = subsub)
      Modperformance$AUC_unit[count]=mm[1]
      Modperformance$Pres_unit[count]=mm[2]
      Modperformance$Abs_unit[count]=mm[3]
      
      Modperformance$GroupId[count]=subsub$GroupId[1]
      Modperformance$UnitLoc[count]=subsub$UnitLoc[1]
      count=count+1
    }
    rm(subsub, mm)
  }    


  MM=merge(ModelTable, Modperformance, by='GroupId', all.x = TRUE)
  
  

  
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
  
  
  
  MM$UnitLoc=factor(MM$UnitLoc, levels=level_names)
  MM=with(MM, MM[order(UnitLoc),])
  MM=MM[! duplicated(MM),]
  
  
# X)  Create coefficents table ###
  

  
  mod=modlist[[1]]
  coeffvals=t(data.frame(coefficients(mod)))
  
  for(ii in 2:10){
    mod=modlist[[ii]]
    NewVals=data.frame(coefficients(mod))
    
    coeffvals=merge(coeffvals, t(NewVals), all=TRUE)

    rm(NewVals)
  }
  
  # Make the order something sensable 
  coeffvals=coeffvals[,c(1,5,3,4,2,seq(6,ncol(coeffvals)))]
  
  # Add group id
  coeffvals$GroupId=levels(OccTable$GroupId)
  coeffvals=coeffvals[, c(ncol(coeffvals),seq(2,(ncol(coeffvals)-1)))]
  

  
# Xx) Plot models and fits for each location ####
  
  
  
  
  
  # Add column for year and shoredist
  OccTable_daily_wDetections$YearShoreDist=paste(OccTable_daily_wDetections$ShoreDist, OccTable_daily_wDetections$Year)

  # Three years one plot
  ggplot(data=OccTable_daily_wDetections) +
    theme_bw() +
    facet_wrap(~GroupId) +
    scale_colour_manual(values=cbbPalette) +
    geom_point(aes(x=JulienDay, y=BNDTotOffset, color=ShoreDist), size=.0005) +
    geom_line(data=subset(fit, Year==2013), aes(JulienDay, inv.logit(fit), colour=ShoreDist)) +
    geom_line(data=subset(fit, Year==2014), aes(JulienDay, inv.logit(fit), colour=ShoreDist)) +
    geom_line(data=subset(fit, Year==2015), aes(JulienDay, inv.logit(fit), colour=ShoreDist)) +
    geom_ribbon(data=fit,aes(x=JulienDay, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
                alpha=.3,linetype= 'blank') +
    ggtitle('Daily Occupancy ' ) +
    ylab('P(BND)')+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  
  # Only 5km Units
  ggplot(data=OccTable_daily_wDetections) +
    theme_bw() +
    facet_wrap(~GroupId, drop = FALSE) +
    scale_colour_manual(values=cbbPalette) +
    geom_point(data=subset(OccTable_daily_wDetections, ShoreDist=='05'),
               aes(x=JulienDay, y=BNDTotOffset, color=Year), size=.0005) +
    geom_line(data=fit[fit$ShoreDist=='05',],
              aes(JulienDay, inv.logit(fit), color=Year)) +
    geom_ribbon(data=fit[fit$ShoreDist=='05',],
                aes(x=JulienDay, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=Year),
                alpha=.3,linetype= 'blank') +
    scale_x_continuous(breaks = pretty(x = fit$JulienDay, n = 3)) +
    ggtitle('Daily Occupancy 05 Sites' )+
    ylab('P(BND)')+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  # Only 10km Units
  ggplot(OccTable_daily_wDetections) +
    theme_bw() +
    facet_wrap(~GroupId) +
    scale_colour_manual(values=cbbPalette) +
    geom_point(data=subset(OccTable_daily_wDetections, ShoreDist=='10'),
               aes(x=JulienDay, y=BNDTotOffset, color=Year), size=.0005) +
    geom_line(data=fit[fit$ShoreDist=='10',],
              aes(JulienDay, inv.logit(fit), color=Year)) +
    geom_ribbon(data=fit[fit$ShoreDist=='10',],
                aes(x=JulienDay, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=Year),
                alpha=.3,linetype= 'blank') +
     scale_x_continuous(breaks = pretty(x = fit$JulienDay, n = 3)) +
    ggtitle('Daily Occupancy 10 Sites' ) +
   
    ylab('P(BND)')+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  # Only 15km Units
  ggplot(OccTable_daily_wDetections) +
    theme_bw() +
    facet_wrap(~GroupId) +
    scale_colour_manual(values=cbbPalette) +
    geom_point(data=subset(OccTable_daily_wDetections, ShoreDist=='15'),
               aes(x=JulienDay, y=BNDTotOffset, color=Year), size=.0005) +
    geom_line(data=fit[fit$ShoreDist=='15',],
              aes(JulienDay, inv.logit(fit), color=Year)) +
    geom_ribbon(data=fit[fit$ShoreDist=='15',],
                aes(x=JulienDay, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=Year),
                alpha=.3,linetype= 'blank') +
    scale_x_continuous(breaks = pretty(x = fit$JulienDay, n = 3)) +
    ggtitle('Daily Occupancy 15 Sites' ) +
    ylab('P(BND)')+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  
  
  
  
  # Only 2013
  png(filename = paste('2013 Occupancy.png'),
      units="in", 
      width=7, 
      height=9, 
      pointsize=12,res = 400)
  
  ggplot(data=subset(OccTable_daily_wDetections, Year==2013)) +
    theme_bw() +
    facet_wrap(~GroupId) +
    scale_colour_manual(values=cbbPalette) +
    geom_point(aes(x=JulienDay, y=BNDTotOffset, color=ShoreDist), size=.0005) +
    geom_line(data=subset(fit, Year==2013), aes(JulienDay, inv.logit(fit), colour=ShoreDist)) +
    geom_ribbon(data=subset(fit, Year==2013),aes(x=JulienDay, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
                alpha=.3,linetype= 'blank') +
    ggtitle('Daily Occupancy Rates 2013' )+
    ylab('P(BND)')+
    theme(plot.title = element_text(hjust = 0.5))
  
  dev.off()
    
  
  png(filename = paste('2014 Occupancy.png'),
      units="in", 
      width=7, 
      height=9, 
      pointsize=12,res = 400)
  
  ggplot(data=subset(OccTable_daily_wDetections, Year==2014)) +
    theme_bw() +
    facet_wrap(~GroupId) +
    scale_colour_manual(values=cbbPalette) +
    geom_point(aes(x=JulienDay, y=BNDTotOffset, color=ShoreDist), size=.0005) +
    geom_line(data=subset(fit, Year==2014), aes(JulienDay, inv.logit(fit), colour=ShoreDist)) +
    geom_ribbon(data=subset(fit, Year==2014),aes(x=JulienDay, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
                alpha=.3,linetype= 'blank') +
    ggtitle('Daily Occupancy Rates 2014')+
    ylab('P(BND)')+
    theme(plot.title = element_text(hjust = 0.5))
  dev.off()
  
  
  png(filename = paste('2015 Occupancy.png'),
      units="in", 
      width=7, 
      height=9, 
      pointsize=12,res = 400)
  ggplot(data=subset(OccTable_daily_wDetections, Year==2015)) +
    theme_bw() +
    facet_wrap(~GroupId) +
    scale_colour_manual(values=cbbPalette) +
    geom_point(aes(x=JulienDay, y=BNDTotOffset, color=ShoreDist), size=.0005) +
    geom_line(data=subset(fit, Year==2015), aes(JulienDay, inv.logit(fit), colour=ShoreDist)) +
    geom_ribbon(data=subset(fit, Year==2015),aes(x=JulienDay, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
                alpha=.3,linetype= 'blank') +
    ggtitle('Daily Occupancy Rates 2015')+
    ylab('P(BND)') +
    theme(plot.title = element_text(hjust = 0.5))
    
  dev.off()

# Xi) Make partial plots for all groups as before ###########
  
  
  
  # Plot Storage
  p=list()
  Sd_P=list()
  Yr_P=list()
  
  # Make partial plots for all deployment groups
  for(ii in 1:10){
    
    # Again subset the data
    data_sub=subset(OccTable_daily_wDetections, GroupId==unique(OccTable$GroupId)[ii])
    data_sub$ShoreDist=factor(data_sub$ShoreDist, levels=c('05', '10', '15'))
    data_sub <- droplevels(data_sub)
    
    # Add a dummy date for plotting
    data_sub$DummyDate=as.Date(data_sub$JulienDay, origin=as.Date("2013-01-01"))
    
    # Extract model (just for clarity)
    mod=modlist[[ii]]
    
    if(!is.character(mod)){
      
      #######################################################
      # Julien Date Smoothes #
      #######################################################
      
      if(length(grep(x = as.character( Reduce(paste, deparse(formula(mod))) ), pattern = 'JulienDay'))>0){
        fitdf_jdate=partialDF(mod = mod, data = data_sub, Variable = 'JulienDay')
        fitdf_jdate$DummyDate=as.Date(fitdf_jdate$JulienDay, origin=as.Date("2013-01-01"))
        
      }else{
        fitdf_jdate=data.frame(JulienDay=data_sub$JulienDay, 
                               DummyDate=as.Date(data_sub$JulienDay, origin=as.Date("2013-01-01")),
                               y=NaN,
                               LCI=NaN,
                               UCI=NaN)
      }
      
      
      
      p[[ii]]=ggplot(data=fitdf_jdate) + 
        geom_line(aes(x=DummyDate, y=y)) +
        geom_ribbon(aes(x=DummyDate, ymin=LCI, ymax=UCI), alpha=.3) +
        geom_rug(data=subset(data_sub, OccAll==1), 
                 aes(x=DummyDate, y=OccAll*.3),
                 sides='t', alpha=.8) +
        geom_rug(data=subset(data_sub, OccAll==0), 
                 aes(x=DummyDate, y=OccAll),
                 sides='b', alpha=.8) +
        xlab('Julien Day') +
        ylab('Detection Probability') +
        theme_minimal() +
        ggtitle(paste('Partial Plot of Julien Day', as.character(unique(data_sub$GroupId))))
      
      
      #######################################################
      # Shore Dist Factors #
      #######################################################
      
      if(length(grep(x = as.character( Reduce(paste, deparse(formula(mod))) ), pattern = 'ShoreDist'))>0){
        
        fitdf_shoredist=partialdf_factor(mod = mod, data = data_sub, variable = 'ShoreDist')
        fitdf_shoredist$GroupId=data_sub$GroupId[1]
        
        
      }else{
        fitdf_shoredist=data.frame(ShoreDist=factor(c('ShoreDist05','ShoreDist10','ShoreDist15')), 
                                   vals=NaN,
                                   GroupId=data_sub$GroupId[1])
      }
      
      
      
      ggplot(fitdf_shoredist, aes(x=ShoreDist, y=inv.logit(vals)))+geom_violin()
      
      #######################################################
      # Year as Factors #
      #######################################################
      
      if(length(grep(x = as.character( Reduce(paste, deparse(formula(mod))) ), pattern = 'Year'))>0){
        fitdf_year = partialdf_factor(mod = mod, data = data_sub, variable = 'Year')
        fitdf_year$GroupId=data_sub$GroupId[1]
        
      }else{
        fitdf_year=data.frame(Year=factor(c('Year2013','Year2014','Year2015')), 
                              vals=NaN,
                              GroupId=data_sub$GroupId[1])}
      ggplot(fitdf_year, aes(x=Year, y=inv.logit(vals)))+geom_violin()
      
      
      
      
      
      
      Yr_P[[ii]]=ggplot(fitdf_year, aes(x=Year, y=vals)) +
        geom_boxplot() +
        theme_minimal()+
        ggtitle(paste('Partial Plot of Years', as.character(unique(data_sub$GroupId))))
      #############################
      # Add all data for plotting #
      ############################
      
      fitdf_jdate$GroupId=data_sub$GroupId[1]
      
      OneYearAggs=data.frame(aggregate(data=data_sub, BBOcc~DayBin + GroupId, FUN=mean))
      OneYearAggs=cbind(OneYearAggs,aggregate(data=data_sub, JulienDay~DayBin+GroupId, FUN=median)[,3])
      
      colnames(OneYearAggs)[4]=c('med')
      
      if(ii==1){
        fitdf_jdate_out=fitdf_jdate
        fitdf_shoredist_out=fitdf_shoredist
        fitdf_Year_out=fitdf_year
        
        AggData=OneYearAggs
        
      }else {
        fitdf_jdate_out=rbind(fitdf_jdate_out, fitdf_jdate)
        fitdf_shoredist_out=rbind(fitdf_shoredist_out, fitdf_shoredist)
        fitdf_Year_out=rbind(fitdf_Year_out, fitdf_year)
        
        AggData=rbind(AggData, OneYearAggs)
      }
      
      
      rm(data_sub, mod, JdateForPlotting, x1, test, BootstrapCoefs, Basis, BootstrapFits, BootstrapParameters)
    }
  }
  
  
  
  
  
  # Bin the data for visualisation #
  OccTable_daily_wDetections$DayBin=cut(OccTable_daily_wDetections$JulienDay, breaks=20)
  
  
  AggData=data.frame(aggregate(data=OccTable_daily_wDetections, OccAll~DayBin+GroupId+ShoreDist, FUN=mean))
  AggData=cbind(AggData, aggregate(data=OccTable_daily_wDetections, JulienDay~DayBin+GroupId+ShoreDist,
                                   FUN=median)[,4])
  
  colnames(AggData)[5]=c('med')
  
  AggData$DummyDate=as.Date(AggData$med, origin=as.Date("2013-01-01"))
  OccTable_daily_wDetections$DummyDate=as.Date(OccTable_daily_wDetections$JulienDay, origin=as.Date("2013-01-01"))
  
  
  ggplot(data=fitdf_jdate_out) +
    theme_bw() +
    facet_wrap(~GroupId) +
    scale_colour_manual(values=cbbPalette) +
    #geom_point(data=AggData, aes(x=DummyDate, y=OccAll), size=.9) +
    geom_line(aes(DummyDate, y), size=1) +
    geom_ribbon(aes(x=DummyDate, ymin=LCI, ymax=UCI),
                alpha=.2,linetype= 'blank') +
    geom_rug(data=subset(OccTable_daily_wDetections, OccAll==1), 
             aes(x=DummyDate, y=BBOcc*.3),
             sides='t', alpha=.8) +
    geom_rug(data=subset(OccTable_daily_wDetections, OccAll==0), 
             aes(x=DummyDate, y=BBOcc),
             sides='b', alpha=.8) +
    xlab("") +
    ylab("") 
  
  
  
  ggplot(data=fitdf_shoredist_out) +
    theme_bw() +
    facet_wrap(~GroupId) +
    geom_boxplot(aes(x=ShoreDist, y=inv.logit(vals))) +
    scale_x_discrete(breaks=unique(fitdf_shoredist_out$ShoreDist),
                     labels=c("Near", "Mid", "Off")) +
    ylab("Occupancy Probability") +
    xlab("")
  
  
  
  ggplot(data=fitdf_Year_out) +
    theme_bw() +
    facet_wrap(~GroupId) +
    geom_boxplot(aes(x=Year, y=inv.logit(vals)))+
    scale_x_discrete(breaks=unique(fitdf_Year_out$Year),
                     labels=c("2013", "2014", "2015")) +
    ylab("Occupancy Probability") +
    xlab("")
  
  
  
  
  
  
  

# 6)  Fit Models and Predictions (make sure to load predictvcv first!) #########################################################



# Add Occupancy threshold
OccTable_daily_wDetections$OccThresh=ifelse(OccTable_daily_wDetections$BNDTotOffset>=0.5, 1,0)



# Address weights for deployment locations where data for 2+ survey days?


mod1=geeglm(BNDTotOffset~GroupId+ShoreDist*bs(JulienDay, knots = mean(JulienDay)), 
         family = binomial, 
         data   = OccTable_daily_wDetections,
         weights = rep(1, nrow(OccTable_daily_wDetections)),
         id     = UnitLoc,
         corstr = 'ar1')

mod2=geeglm(BNDTotOffset~ShoreDist+GroupId*bs(JulienDay, knots = mean(JulienDay)), 
            family = binomial, 
            data   = OccTable_daily_wDetections,
            weights = rep(1, nrow(OccTable_daily_wDetections)),
            id     = UnitLoc,
            corstr = 'ar1')

mod3=geeglm(BNDTotOffset~UnitLoc+bs(JulienDay, knots = mean(JulienDay)), 
            family = binomial, 
            data   = OccTable_daily_wDetections,
            weights = rep(1, nrow(OccTable_daily_wDetections)),
            id     = UnitLoc,
            corstr = 'ar1')


QIC(mod1, mod2, mod3)

# Version 2- select only BB days 
mod4=geeglm(BBOcc~GroupId+ShoreDist*bs(JulienDay, knots = mean(JulienDay)), 
            family = binomial, 
            data   = OccTable_daily_wDetections,
            #weights = rep(1, nrow(OccTable_daily_wDetections)),
            id     = UnitLoc,
            corstr = 'ar1')

mod5=geeglm(BBOcc~ShoreDist+GroupId*bs(JulienDay, knots = mean(JulienDay)), 
            family = binomial, 
            data   = OccTable_daily_wDetections,
            #weights = rep(1, nrow(OccTable_daily_wDetections)),
            id     = UnitLoc,
            corstr = 'ar1')

mod6=geeglm(BBOcc~UnitLoc+bs(JulienDay, knots = mean(JulienDay)), 
            family = binomial, 
            data   = OccTable_daily_wDetections,
            #weights = rep(1, nrow(OccTable_daily_wDetections)),
            id     = UnitLoc,
            corstr = 'ar1')
            
QIC(mod4, mod5, mod6)

# Version 3- select days where probability is greater than .4

mod7=geeglm(OccThresh~GroupId+ShoreDist*bs(JulienDay, knots = mean(JulienDay)), 
            family = binomial, 
            data   = OccTable_daily_wDetections,
            #weights = rep(1, nrow(OccTable_daily_wDetections)),
            id     = UnitLoc,
            corstr = 'ar1')

mod8=geeglm(OccThresh~ShoreDist+GroupId*bs(JulienDay, knots = mean(JulienDay)), 
            family = binomial, 
            data   = OccTable_daily_wDetections,
            #weights = rep(1, nrow(OccTable_daily_wDetections)),
            id     = UnitLoc,
            corstr = 'ar1')

mod9=geeglm(OccThresh~UnitLoc+bs(JulienDay, knots = mean(JulienDay)), 
            family = binomial, 
            data   = OccTable_daily_wDetections,
            #weights = rep(1, nrow(OccTable_daily_wDetections)),
            id     = UnitLoc,
            corstr = 'ar1')

QIC(mod7, mod8, mod9)



# Fill in the Model Table

# Table to store model performance #
ModelTable=data.frame(Model=ls(pattern = 'mod'))
ModelTable$ModelFormula='none'
ModelTable$QIC=0
ModelTable$AUC=0 # Area under the Curve 
ModelTable$Pres=0 # Proportion of the presences correctly identified 
ModelTable$Abs=0 # Proportion of the absences correctly idenified



# Fill in the Model Table
for(ii in 1:9){
  mod=ls(pattern = 'mod')[ii]
  ModelTable$ModelFormula[ii]=Reduce(paste, deparse(formula(c(mod))))
  ModelTable$QIC[ii]= QIC(get(eval(mod)))
  
  AUC_vals= CalcAUC(get(eval(mod)), OccTable_daily_wDetections)
  ModelTable$AUC[ii]=AUC_vals[1]
  ModelTable$Pres[ii]=AUC_vals[2]
  ModelTable$Abs[ii]=AUC_vals[3]
  rm(mod)
}



# 3) Partial Plots for the top three GEEGLM's ######################################################



# Third model retained in each group

mod1=mod3
mod2=mod6
mod3=mod9

# Plot model 1
ncol_df= ncol(OccTable_daily_wDetections)
OccTable_daily_wDetections[,ncol_df:(ncol_df+3)]=predictvcv(mod1, newdata = OccTable_daily_wDetections)

ggplot(data=OccTable_daily_wDetections) +
  theme_bw() +
  facet_wrap(~GroupId) +
  scale_colour_manual(values=cbbPalette) +
  geom_point(aes(x=JulienDay, y=BNDTotOffset, color=ShoreDist), size=.0005) +
  geom_line(aes(JulienDay, inv.logit(fit), colour=ShoreDist), size=1) +
  geom_ribbon(aes(x=JulienDay, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
              alpha=.3,linetype= 'blank') +
  ggtitle('Daily Classification Rates')



# Plot model 2
OccTable_daily_wDetections[,ncol_df:(ncol_df+3)]=predictvcv(mod2, newdata = OccTable_daily_wDetections)

ggplot(data=OccTable_daily_wDetections) +
  theme_bw() +
  facet_wrap(~GroupId) +
  scale_colour_manual(values=cbbPalette) +
  geom_point(aes(x=JulienDay, y=BBOcc, color=ShoreDist), size=.0005) +
  geom_line(aes(JulienDay, inv.logit(fit), colour=ShoreDist), size=1) +
  geom_ribbon(aes(x=JulienDay, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
              alpha=.3,linetype= 'blank') +
  ggtitle('Only BND Days')


# Plot model 3
OccTable_daily_wDetections[,ncol_df:(ncol_df+3)]=predictvcv(mod3, newdata = OccTable_daily_wDetections)

ggplot(data=OccTable_daily_wDetections) +
  theme_bw() +
  facet_wrap(~GroupId) +
  scale_colour_manual(values=cbbPalette) +
  geom_point(aes(x=JulienDay, y=OccThresh, color=ShoreDist), size=.0005) +
  geom_line(aes(JulienDay, inv.logit(fit), colour=ShoreDist), size=1) +
  geom_ribbon(aes(x=JulienDay, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
              alpha=.3,linetype= 'blank') +
  ggtitle('BND days with prob>.4')


# Partial plots for factors
yeardat=partialdf_factor(mod1, OccTable_daily_wDetections, 'Year')
shoredat=partialdf_factor(mod1, OccTable_daily_wDetections, 'ShoreDist')
DateDat=partialDF(mod1, OccTable_daily_wDetections, 'JulienDay')

ggplot(yeardat, aes(x=Year, y=vals)) +
  geom_boxplot() +
  theme_minimal()+
  ggtitle(paste('Partial Plot of Years', as.character(unique(OccTable_daily_wDetections$GroupId))))

ggplot(shoredat, aes(x=ShoreDist, y=vals)) +
  geom_boxplot() +
  theme_minimal()+
  ggtitle(paste('Partial Plot of Years', as.character(unique(OccTable_daily_wDetections$GroupId))))


ggplot(data=DateDat) + 
  geom_line(aes(x=JulienDay, y=y)) +
  geom_ribbon(aes(x=JulienDay, ymin=LCI, ymax=UCI), alpha=.3) +
  geom_rug(data=subset(OccTable_daily_wDetections, OccAll==1), 
           aes(x=JulienDay, y=OccAll*.01),
           sides='t', alpha=.8) +
  geom_rug(data=subset(OccTable_daily_wDetections, OccAll==0), 
           aes(x=JulienDay, y=OccAll),
           sides='b', alpha=.8) +
  xlab('Julien Day') +
  ylab('Detection Probability') +
  theme_minimal() +
  ggtitle(paste('Partial Plot of Julien Day'))




# 7)  Make Unique Plots for all of the Group ID's ###################################
# 7)  Make Partial plots for QIC selected Models ##################################################################


# Plot Storage
p=list()
Sd_P=list()
Yr_P=list()

# Make partial plots for all deployment groups
for(ii in 1:10){
  
  # Again subset the data
  data_sub=subset(OccTable_daily_wDetections, GroupId==unique(OccTable$GroupId)[ii])
  data_sub$ShoreDist=factor(data_sub$ShoreDist, levels=c('05', '10', '15'))
  data_sub <- droplevels(data_sub)
  
  # Add a dummy date for plotting
  data_sub$DummyDate=as.Date(data_sub$JulienDay, origin=as.Date("2013-01-01"))
  
  # Extract model (just for clarity)
  mod=modsig[[ii]]
  
  if(!is.character(mod)){
  
  #######################################################
  # Julien Date Smoothes #
  #######################################################
  
  if(length(grep(x = as.character( Reduce(paste, deparse(formula(mod))) ), pattern = 'JulienDay'))>0){
  fitdf_jdate=partialDF(mod = mod, data = data_sub, Variable = 'JulienDay')
  fitdf_jdate$DummyDate=as.Date(fitdf_jdate$JulienDay, origin=as.Date("2013-01-01"))
  
  }else{
    fitdf_jdate=data.frame(JulienDay=data_sub$JulienDay, 
                           DummyDate=as.Date(data_sub$JulienDay, origin=as.Date("2013-01-01")),
                           y=NaN,
                           LCI=NaN,
                           UCI=NaN)
  }
  
  
  
  p[[ii]]=ggplot(data=fitdf_jdate) + 
    geom_line(aes(x=DummyDate, y=y)) +
    geom_ribbon(aes(x=DummyDate, ymin=LCI, ymax=UCI), alpha=.3) +
    geom_rug(data=subset(data_sub, OccAll==1), 
             aes(x=DummyDate, y=OccAll*.3),
             sides='t', alpha=.8) +
    geom_rug(data=subset(data_sub, OccAll==0), 
             aes(x=DummyDate, y=OccAll),
             sides='b', alpha=.8) +
    xlab('Julien Day') +
    ylab('Detection Probability') +
    theme_minimal() +
    ggtitle(paste('Partial Plot of Julien Day', as.character(unique(data_sub$GroupId))))
  
  
  #######################################################
  # Shore Dist Factors #
  #######################################################
  
  if(length(grep(x = as.character( Reduce(paste, deparse(formula(mod))) ), pattern = 'ShoreDist'))>0){
    
    fitdf_shoredist=partialdf_factor(mod = mod, data = data_sub, variable = 'ShoreDist')
    fitdf_shoredist$GroupId=data_sub$GroupId[1]
    
  
    }else{
    fitdf_shoredist=data.frame(ShoreDist=factor(c('ShoreDist05','ShoreDist10','ShoreDist15')), 
                           vals=NaN,
                           GroupId=data_sub$GroupId[1])
  }
  

  
  ggplot(fitdf_shoredist, aes(x=ShoreDist, y=inv.logit(vals)))+geom_violin()
  
  #######################################################
  # Year as Factors #
  #######################################################

 if(length(grep(x = as.character( Reduce(paste, deparse(formula(mod))) ), pattern = 'Year'))>0){
     fitdf_year = partialdf_factor(mod = mod, data = data_sub, variable = 'Year')
     fitdf_year$GroupId=data_sub$GroupId[1]
     
   }else{
     fitdf_year=data.frame(Year=factor(c('Year2013','Year2014','Year2015')), 
                                vals=NaN,
                                GroupId=data_sub$GroupId[1])}
  ggplot(fitdf_year, aes(x=Year, y=inv.logit(vals)))+geom_violin()
  
  



  
   Yr_P[[ii]]=ggplot(fitdf_year, aes(x=Year, y=vals)) +
     geom_boxplot() +
     theme_minimal()+
     ggtitle(paste('Partial Plot of Years', as.character(unique(data_sub$GroupId))))
  #############################
  # Add all data for plotting #
  ############################
   
  fitdf_jdate$GroupId=data_sub$GroupId[1]
   
  OneYearAggs=data.frame(aggregate(data=data_sub, BBOcc~DayBin + GroupId, FUN=mean))
  OneYearAggs=cbind(OneYearAggs,aggregate(data=data_sub, JulienDay~DayBin+GroupId, FUN=median)[,3])
   
   colnames(OneYearAggs)[4]=c('med')
  
  if(ii==1){
    fitdf_jdate_out=fitdf_jdate
    fitdf_shoredist_out=fitdf_shoredist
    fitdf_Year_out=fitdf_year
    
    AggData=OneYearAggs
    
  }else {
    fitdf_jdate_out=rbind(fitdf_jdate_out, fitdf_jdate)
    fitdf_shoredist_out=rbind(fitdf_shoredist_out, fitdf_shoredist)
    fitdf_Year_out=rbind(fitdf_Year_out, fitdf_year)
    
    AggData=rbind(AggData, OneYearAggs)
  }
   
   
  rm(data_sub, mod, JdateForPlotting, x1, test, BootstrapCoefs, Basis, BootstrapFits, BootstrapParameters)
  }
}





# Bin the data for visualisation #
OccTable_daily_wDetections$DayBin=cut(OccTable_daily_wDetections$JulienDay, breaks=20)


AggData=data.frame(aggregate(data=OccTable_daily_wDetections, OccAll~DayBin+GroupId+ShoreDist, FUN=mean))
AggData=cbind(AggData, aggregate(data=OccTable_daily_wDetections, JulienDay~DayBin+GroupId+ShoreDist,
                                 FUN=median)[,4])

colnames(AggData)[5]=c('med')

AggData$DummyDate=as.Date(AggData$med, origin=as.Date("2013-01-01"))
OccTable_daily_wDetections$DummyDate=as.Date(OccTable_daily_wDetections$JulienDay, origin=as.Date("2013-01-01"))


ggplot(data=fitdf_jdate_out) +
  theme_bw() +
  facet_wrap(~GroupId) +
  scale_colour_manual(values=cbbPalette) +
  #geom_point(data=AggData, aes(x=DummyDate, y=OccAll), size=.9) +
  geom_line(aes(DummyDate, y), size=1) +
  geom_ribbon(aes(x=DummyDate, ymin=LCI, ymax=UCI),
              alpha=.2,linetype= 'blank') +
  geom_rug(data=subset(OccTable_daily_wDetections, OccAll==1), 
           aes(x=DummyDate, y=BBOcc*.3),
           sides='t', alpha=.8) +
  geom_rug(data=subset(OccTable_daily_wDetections, OccAll==0), 
           aes(x=DummyDate, y=BBOcc),
           sides='b', alpha=.8) +
  xlab("") +
  ylab("") 



ggplot(data=fitdf_shoredist_out) +
  theme_bw() +
  facet_wrap(~GroupId) +
  geom_boxplot(aes(x=ShoreDist, y=inv.logit(vals))) +
  scale_x_discrete(breaks=unique(fitdf_shoredist_out$ShoreDist),
                     labels=c("Near", "Mid", "Off")) +
  ylab("Occupancy Probability") +
  xlab("")

 

ggplot(data=fitdf_Year_out) +
  theme_bw() +
  facet_wrap(~GroupId) +
  geom_boxplot(aes(x=Year, y=inv.logit(vals)))+
  scale_x_discrete(breaks=unique(fitdf_Year_out$Year),
                   labels=c("2013", "2014", "2015")) +
  ylab("Occupancy Probability") +
  xlab("")



# Ass ton of model selection code that Luke Nixed ##################





# list to store the models #
modlist=list()

# List to store the signficant models
modsig=list()



# Model selection for ten variables
for(ii in 1:10){
  
  
  data_sub=subset(OccTable_daily_wDetections, GroupId==unique(OccTable$GroupId)[ii])
  data_sub$ShoreDist=factor(data_sub$ShoreDist, levels=c('05', '10', '15'))
  data_sub <- droplevels(data_sub)
  
  
  ModelTable$Nunits2013[ii]=length(unique(data_sub$UnitLoc[data_sub$Year==2013]))
  ModelTable$Nunits2014[ii]=length(unique(data_sub$UnitLoc[data_sub$Year==2014]))
  ModelTable$Nunits2015[ii]=length(unique(data_sub$UnitLoc[data_sub$Year==2015]))
  
  newdat=subset(data_sub, select=c('JulienDay', 'ShoreDist', 'GroupId', 'UnitLoc', 'OccAll', 'Year','BNDTotOffset'))
  
  tryCatch({
    ModelTable$Data2013[ii]<-as.character(Reduce(paste, unique(data_sub$UnitLoc[data_sub$Year==2013])))}, error=function(e){})
  
  tryCatch({
    ModelTable$Data2014[ii]<-as.character(Reduce(paste, unique(data_sub$UnitLoc[data_sub$Year==2014])))}, error=function(e){})
  
  
  tryCatch({
    ModelTable$Data2015[ii]<-as.character(Reduce(paste, unique(data_sub$UnitLoc[data_sub$Year==2015])))}, error=function(e){})
  
  
  # Determine whether linear model or spline for Julien Day #
  
  
  
  
  mod1=geeglm(BNDTotOffset~ShoreDist+bs(JulienDay, knots = mean(JulienDay)), 
              corstr = 'ar1', 
              family = binomial, 
              id     = UnitLoc, 
              data   = data_sub,
              weights = rep(1, nrow(data_sub))) 
  
  
  
  
  
  
  null=geeglm(OccAll~ShoreDist, 
              corstr = 'ar1', 
              offset = tempoffset, 
              family = binomial, 
              id     = UnitLoc, 
              data   = data_sub) 
  
  mod1=geeglm(OccAll~ShoreDist+bs(JulienDay, knots = mean(JulienDay)), 
              corstr = 'ar1', 
              offset  = tempoffset, 
              family = binomial, 
              id     = UnitLoc, 
              data   = data_sub) 
  
  mod2=geeglm(OccAll~ShoreDist+JulienDay, 
              corstr = 'ar1', 
              offset = tempoffset, 
              family = binomial, 
              id     = UnitLoc, 
              data   = data_sub) 
  
  QICvals=QIC(null, mod1, mod2)
  
  if (sum(QICvals>10e+15)>2){
    
    null=geeglm(OccAll~Year, 
                corstr = 'ar1', 
                offset = tempoffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub) 
    
    mod1=geeglm(OccAll~Year+bs(JulienDay, knots = mean(JulienDay)), 
                corstr = 'ar1', 
                offset = tempoffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub)
    
    mod2=geeglm(OccAll~Year+JulienDay, 
                corstr = 'ar1', 
                offset = tempoffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub)
    
    QICvals=QIC(null, mod1, mod2)}
  
  JdateForm=which.min(QICvals[2:3,1]-QICvals[1,1])
  if(JdateForm==1){
    JdateForm='bs(JulienDay, knots = mean(JulienDay))'
  }else{
    JdateForm='JulienDay'
  }
  
  
  ModelFull=geeglm(as.formula(paste('OccAll~', JdateForm, '+ ShoreDist + Year')), 
                   corstr = 'ar1', 
                   offset = tempoffset, 
                   family = binomial, 
                   id     = UnitLoc, 
                   data   = data_sub)
  
  
  # Bakcwards selection for QIC
  modlist[[ii]]=SelectModel(ModelFull)
  
  # Document the formula
  ModelTable$ModelFormula[ii]=Reduce(paste, deparse(formula(modlist[[ii]])))
  
  # The Data are too sparse to support removal of unsignifcant terms
  # However, significant terms can be displayed
  
  # Walds test function to remove unsignificant terms
  tempmod=DropVarsWalds(modlist[[ii]])
  
  # store the mode
  modsig[[ii]]=tempmod
  
  
  if(is.character(tempmod)){
    ModelTable$WaldsSigVars[ii]=tempmod
    
    # If covariates were retained then report them, otherwise fill in model table with N values
    ModelTable$AUC[ii]='NA'
    ModelTable$Pres[ii]='NA'
    ModelTable$Abs[ii]='NA'
    
  }else{
    
    ModelTable$WaldsSigVars[ii]=Reduce(paste, deparse(formula(tempmod)))
    
    terms_orig=attr(tempmod$terms,"term.labels")
    terms_idx= grep('JulienDay', terms_orig)
    terms=terms_orig
    
    # If Julien day present in the model replace it for aggregating
    if(length(grep('JulienDay', terms_orig))>0){
      terms[terms_idx]='DayBin'
    }
    rm(terms_idx)
    
    
    # Expand terms
    
    if(length(grep('JulienDay', terms_orig))>0){jday=seq(min(newdat$JulienDay), max(newdat$JulienDay))
    }else{
      jday=median(newdat$JulienDay)}
    
    if(length(grep('Year', terms_orig))>0){Year=data_sub$Year[1]}else{
      Year=unique(data_sub$Year)}
    
    if(length(grep('ShoreDist', terms_orig))>0){ShoreDist=data_sub$ShoreDist[1]}else{
      ShoreDist=unique(data_sub$ShoreDist)}
    
    
    # Make new prediction grid using the terms from the model
    newdat_perdOnly=expand.grid(JulienDay=jday, Year=Year, ShoreDist=ShoreDist, GroupId=unique(newdat$GroupId),OccAll=0)
    
    
    # Aggregated data for plotting
    form=as.formula(paste('OccAll ~', paste(terms, collapse=" + ")))
    form1=as.formula(paste('JulienDay ~', paste(terms, collapse=" + ")))
    form2=as.formula(paste('BNDTotOffset ~', paste(terms, collapse=" + ")))
    
    OneYearAggs=data.frame(aggregate(data=data_sub, form, FUN=mean))
    OneYearAggs=cbind(OneYearAggs, aggregate(data=data_sub, form1, FUN=median)[,(length(terms)+1)])
    OneYearAggs=cbind(OneYearAggs, aggregate(data=data_sub, form2, FUN=mean)[,(length(terms)+1)])
    
    
    colnames(OneYearAggs)[(length(terms)+2)]=c('med')
    colnames(OneYearAggs)[(length(terms)+3)]='BNDTotOffset'
    
    if(length(grep('ShoreDist', colnames(OneYearAggs)))==0){OneYearAggs$ShoreDist=data_sub$ShoreDist[1]}
    if(length(grep('Year', colnames(OneYearAggs)))==0){OneYearAggs$Year=data_sub$Year[1]}
    if(length(grep('JulienDay', colnames(OneYearAggs)))==0){OneYearAggs$JulienDay=data_sub$JulienDay[1]
    OneYearAggs$DayBin= data_sub$JulienDay[1]}
    
    
    
    OneYearAggs$JulienDay=OneYearAggs$med
    OneYearAggs$GroupId=unique(newdat$GroupId)
    
    # Report AUC values
    AUCvals=CalcAUC(tempmod, data_sub = data_sub)
    ModelTable$AUC[ii]=round(AUCvals[1],2)
    ModelTable$Pres[ii]=round(AUCvals[2],2)
    ModelTable$Abs[ii]=round(AUCvals[3],2)  
    
    # Make Predictions
    if(ii==1){
      fit=cbind(newdat,  predictvcv(tempmod))
      dummyfit=cbind(newdat_perdOnly, predictvcv(tempmod, newdata = newdat_perdOnly))
      AggData=cbind(OneYearAggs, predictvcv(tempmod, newdata = OneYearAggs))
      
    }else {
      fit=rbind(fit, cbind(newdat, predictvcv(tempmod)) )
      dummyfit=rbind(dummyfit, cbind(newdat_perdOnly, predictvcv(tempmod, newdata = newdat_perdOnly)))
      AggData=rbind(AggData, cbind(OneYearAggs, predictvcv(tempmod, newdata = OneYearAggs)))
    }
  }  
  
  rm(null, mod1, mod2, ModelFull, tempmod, newdat_perdOnly, OneYearAggs, data_sub, AUCvals, form, form1, form2)
  
  
}


# Add dummy dates for plotting (to fix the X axis)
AggData$DummyDate=as.Date(AggData$med, origin=as.Date("2013-01-01"))
AggData$UnitLoc=paste(AggData$GroupId, AggData$ShoreDist, sep="_")
AggData$BBEst=AggData$BNDTotOffset*AggData$OccAll

fit$DummyDate=as.Date(fit$JulienDay, origin=as.Date("2013-01-01"))
dummyfit$DummyDate=as.Date(dummyfit$JulienDay, origin=as.Date("2013-01-01"))




temp=dummyfit[!duplicated(dummyfit$GroupId),]

# Plot all the data
ggplot(data=dummyfit) +
  theme_bw() +
  facet_wrap(~GroupId) +
  scale_colour_manual(values=cbbPalette) +
  geom_line(aes(DummyDate, inv.logit(fit), colour=ShoreDist), size=1) +
  annotate("text", x=as.Date("2013-08-15"), y=1, label= as.character(temp$Year)) +
  geom_ribbon(aes(x=DummyDate, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
              alpha=.2,linetype= 'blank') +
  geom_point(data=AggData, aes(x=DummyDate, y=(BBEst),
                               color=ShoreDist), size=.9) +
  xlab("") +
  ylab("")

# Re-plot but limit the x axis to points where data were collected
ggplot(data=AggData, aes(x=DummyDate, y=(BBEst),
                         color=ShoreDist), size=.9) +
  theme_bw() +
  facet_wrap(~GroupId) +
  geom_point() +
  scale_colour_manual(values=cbbPalette) +
  annotate("text", x=as.Date("2013-08-15"), y=1, label= as.character(temp$Year)) +
  geom_ribbon(aes(x=DummyDate, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
              alpha=.2,linetype= 'blank') +
  geom_line(aes(x=DummyDate, y=inv.logit(fit)), size=.5) +
  xlab("") +
  ylab("")
















