#########################################################################
# Descriptive Statistics For Occupancy#
#########################################################################



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
#library("mvtnorm", lib.loc="~/R/win-library/3.3") # For mtvnorm- partial plots
library(MASS) # for mvrnorm in boostrapping intervals 
library(mvtnorm)
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix

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

meta_sub=subset(meta, select=c('UnitLoc', 'Slope'))
OccTable=merge(OccTable, meta_sub, all.x = TRUE)
rm(meta_sub)


################################################################################
# Function to calculate AUC #
################################################################################

# This function calculates AUC for the model (to access model fit) taken from  
# Pirotta E, Matthiopoulos J, MacKenzie M, Scott-Hayward L, Rendell L Modelling sperm whale habitat preference: a novel approach combining transect and follow data

CalcAUC<-function(mod, data_sub){
  
  pr <- predict(mod,data_sub, type="response")                          # the final model is used to predict the data on the response scale (i.e. a value between 0 and 1)
  pred <- prediction(pr,data_sub$OccAll)                                    # to specify the vector of predictions (pr) and the vector of labels (i.e. the observed values "Pres")
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
  DATA$Predicted<-predict(modlist[[ii]],data_sub,type="response")                 # the third column reports the predictions
  cmx(DATA, threshold = Best_cutoff)                                   # the identified cut-off must be used here
  
  # Area under the Curve 
  auc <- unlist(performance(pred, measure="auc")@y.values)
  
  # Proportion of the presences correctly identified 
  pres=prop.table(cmx(DATA, threshold = Best_cutoff))[1,1]
  
  # Proportion of the absences correctly idenified
  abs=prop.table(cmx(DATA, threshold = Best_cutoff))[2,2]
  
  
  return(c(auc, pres, abs))
}


################################################################################
# General Data Prep #
################################################################################

OccTable$GroupId=unlist(strsplit(as.character(OccTable$UnitLoc), split = "_"))[seq(1,(nrow(OccTable)*2)-1,2)]

level_names=c( "Lat", "Hel", "Cro",
               "SpB", "Fra", "Cru",
               "Sto", "Abr", "StA",
               "Stb")

OccTable$GroupId=factor(OccTable$GroupId, levels=level_names)



#OccTable$ShoreDist=unlist(strsplit(as.character(OccTable$UnitLoc), split = "_"))[seq(2,(nrow(OccTable)*2),2)]


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


mm=distinct(OccTable, Date, UnitLoc, JulienDay, GroupId, ShoreDist, Slope, Year, Month)
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


###############################################################################
# Data vis #
##############################################################################


# Bin the data for visualisation #
OccTable_daily$DayBin=cut(OccTable_daily$JulienDay, breaks=20)


mm=data.frame(aggregate(data=OccTable_daily, BBOcc~DayBin+GroupId+ShoreDist, FUN=mean))
mm=cbind(mm, aggregate(data=OccTable_daily, JulienDay~DayBin+GroupId+ShoreDist, FUN=median)[,4])
mm=cbind(mm, aggregate(data=OccTable_daily, BBOcc~DayBin+GroupId+ShoreDist, FUN=sum)[,4])
mm=cbind(mm, aggregate(data=OccTable_daily, BBOcc~DayBin+GroupId+ShoreDist, FUN=length)[,4])

colnames(mm)[5:7]=c('med', 'sum', 'n')

library('Hmisc')
mm=cbind(mm, binconf(x=mm$sum, n=mm$n, )[,2:3])

library(ggplot2)
ggplot(data=mm, aes(x=med, y=BBOcc, color=ShoreDist)) +
  theme_bw() +
  facet_wrap(~GroupId) +
  geom_point() +
  xlab('Julien Day') +
  ylab('Detection Probability') +
  ggtitle('BND Occupancy')


#################################################################################
# Table For Number of Detections at Each Unit Loc #
#################################################################################

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

##################################################################################
# Basic Table for BND Trains #
##################################################################################

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

# Reorder

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






#################################################################################
# Different Spline for all Locs #
#################################################################################


OccTable_daily$YrUnitLoc=paste(OccTable_daily$UnitLoc, OccTable_daily$Year, sep='')

OccTable_daily_wDetections= OccTable_daily[OccTable_daily$YrUnitLoc %in%
                                             names(which(tapply(OccTable_daily$OccAll,OccTable_daily$YrUnitLoc,sum)>1)),]


OccTable_daily_wDetections=droplevels(OccTable_daily_wDetections)


ModelTable=data.frame(DeplotymentLoc=unique(OccTable$GroupId))
ModelTable$ModelFormula='dunno'
ModelTable$Data2013='none'
ModelTable$Data2014='none'
ModelTable$Data2015='none'
ModelTable$Nunits2013= 0
ModelTable$Nunits2014= 0
ModelTable$Nunits2015= 0
ModelTable$CorrStuct='none'
ModelTable$RsquaredAdj=0
ModelTable$AUC=0
ModelTable$Pres=0
ModelTable$Abs=0


# list to store the models
# list to store the models
modlist=list()

for(ii in 1:10){
    data_sub=subset(OccTable_daily_wDetections, GroupId==unique(OccTable$GroupId)[ii])
    data_sub$ShoreDist=factor(data_sub$ShoreDist, levels=c('05', '10', '15'))
    data_sub <- droplevels(data_sub)
    
    
    
    
    ModelTable$Nunits2013[ii]=length(unique(data_sub$UnitLoc[data_sub$Year==2013]))
    ModelTable$Nunits2014[ii]=length(unique(data_sub$UnitLoc[data_sub$Year==2014]))
    ModelTable$Nunits2014[ii]=length(unique(data_sub$UnitLoc[data_sub$Year==2015]))
  
    newdat=subset(data_sub, select=c('JulienDay', 'ShoreDist', 'GroupId', 'UnitLoc', 'OccAll', 'Year'))
    newdat_perdOnly=expand.grid(JulienDay=seq(min(newdat$JulienDay), max(newdat$JulienDay)),
                                OccAll=0,
                                ShoreDist=unique(newdat$ShoreDist),
                                GroupId=unique(newdat$GroupId),
                                Year=aggregate(data=data_sub, BBOcc~Year, FUN=mean)[ which.max(aggregate(data=data_sub, BBOcc~Year, FUN=mean)[,2]),1])
    
    
    tryCatch({
      ModelTable$Data2013[ii]<-as.character(Reduce(paste, unique(data_sub$UnitLoc[data_sub$Year==2013])))}, error=function(e){})
    
    tryCatch({
      ModelTable$Data2014[ii]<-as.character(Reduce(paste, unique(data_sub$UnitLoc[data_sub$Year==2014])))}, error=function(e){})
    
    
    tryCatch({
      ModelTable$Data2015[ii]<-as.character(Reduce(paste, unique(data_sub$UnitLoc[data_sub$Year==2015])))}, error=function(e){})
    
    
    
    # At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
    # (the ones that, if removed, determine the biggest increase in QIC enter the model first). # Pilfered from Priotta Sperm Whales 
    mod1=geeglm(OccAll~bs(JulienDay)+ShoreDist+Year, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub) 
    
    mod2=geeglm(OccAll~bs(JulienDay)+ShoreDist, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub)   
    mod3=geeglm(OccAll~bs(JulienDay)+Year, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub)
    mod4=geeglm(OccAll~ShoreDist+Year, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub)
    
    Qicdf=data.frame(QIC(mod1, mod2, mod3, mod4))
    Qicdf$deltaQIC=abs(Qicdf$QIC-Qicdf$QIC[1])
    Qicdf$Varnames=c('All',  'Year', 'ShoreDist' ,'bs(JulienDay)')
    Qicdf=Qicdf[order(Qicdf$deltaQIC,decreasing = TRUE),]
    
    gee_form=as.formula(paste('OccAll~', Qicdf$Varnames[1], '+',
                              Qicdf$Varnames[2], '+',
                              Qicdf$Varnames[3]))
    if(ii==6 |ii==5|ii==7){
      modlist[[ii]]=geeglm(formula = gee_form, 
                           corstr = 'independence', 
                           offset = BNDTotOffset, 
                           family = binomial, 
                           id     = UnitLoc, 
                           data   = data_sub) 
    }else{
      modlist[[ii]]=geeglm(formula = gee_form, 
                           corstr = 'ar1', 
                           offset = BNDTotOffset, 
                           family = binomial, 
                           id     = UnitLoc, 
                           data   = data_sub)}
    
  
  
  
  if(ii==1){
    fit=cbind(newdat,  predictvcv(modlist[[ii]]))
    dummyfit=cbind(newdat_perdOnly, predictvcv(modlist[[ii]], newdata = newdat_perdOnly))
  }else {
    fit=rbind(fit, cbind(newdat, predictvcv(modlist[[ii]])) )
    dummyfit=rbind(dummyfit,cbind(newdat_perdOnly, predictvcv(modlist[[ii]], newdata = newdat_perdOnly)))
  }
  
    
  ModelTable$ModelFormula[ii]=Reduce(paste, deparse(formula(modlist[[ii]])))  
  ModelTable$CorrStuct[ii]=modlist[[ii]]$corstr
  
  #Calculate conditional and marginal coefficient of determination for Generalized mixed-effect models (R_GLMM²).
  ModelTable$RsquaredAdj[ii]=r.squaredGLMM(modlist[[ii]])
  
  AUCvals=CalcAUC(modlist[[ii]], data_sub = data_sub)
  ModelTable$AUC[ii]=AUCvals[1]
  ModelTable$Pres[ii]=AUCvals[2]
  ModelTable$Abs[ii]=AUCvals[3]

rm(mod1,mod2, mod3, mod4, Qicdf)
}


mm$DummyDate=as.Date(mm$med, origin=as.Date("2013-01-01"))
mm$UnitLoc=paste(mm$GroupId, mm$ShoreDist, sep="_")
mm_nocro=mm[mm$UnitLoc != 'Cro_05',]

fit$DummyDate=as.Date(fit$JulienDay, origin=as.Date("2013-01-01"))

dummyfit$DummyDate=as.Date(dummyfit$JulienDay, origin=as.Date("2013-01-01"))



# colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
  geom_point(data=mm, aes(x=DummyDate, y=BBOcc,
                          color=ShoreDist), size=.9) 




# Plot the data without Cromarty 05
dummyfit$UnitLoc=paste(dummyfit$GroupId, dummyfit$ShoreDist, sep = '_')
ggplot(data=dummyfit[dummyfit$UnitLoc != 'Cro_05',]) +
  theme_bw() +
  facet_wrap(~GroupId) +
  scale_colour_manual(values=cbbPalette) +
  geom_line(aes(DummyDate, inv.logit(fit), colour=ShoreDist), size=1) +
  geom_ribbon(aes(x=DummyDate, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
              alpha=.2,linetype= 'blank') +
  geom_point(data=mm[mm$UnitLoc !='Cro_05',], aes(x=DummyDate, y=BBOcc,
                                                  color=ShoreDist), size=.9) 





# Get the Partial Residuals for each plot Group ID

# Plot Storage
p=list()
Sd_P=list()
Yr_P=list()

# in the model coeficients there is not 2013 or 05 shroe distance
for(ii in 1:10){
  
  data_sub=subset(OccTable_daily_wDetections, GroupId==unique(OccTable$GroupId)[ii])
  data_sub$ShoreDist=factor(data_sub$ShoreDist, levels=c('05', '10', '15'))
  data_sub <- droplevels(data_sub)
  
  
  data_sub$DummyDate=as.Date(data_sub$JulienDay, origin=as.Date("2013-01-01"))
  mod=modlist[[ii]]
  
  JdateForPlotting<- seq(min(data_sub$JulienDay), max(data_sub$JulienDay))
  YearsForPlotting=unique(data_sub$Year)
  ShoreDistForPlotting=unique(data_sub$ShoreDist)
  
  # Trim Julien Date to be divisible by 10 
  if (length(JdateForPlotting) %% 10>0){
    JdateForPlotting=JdateForPlotting[-runif(length(JdateForPlotting) %% 10, min=1, max=length(JdateForPlotting))]
  }
  
  BootstrapParameters<-mvrnorm(10000, coef(mod), summary(mod)$cov.unscaled)
  #######################################################
  # Julien Date Smoothes #
  #######################################################
  # Get the smoothed index terms 
  BS_idx=which(grepl("bs", colnames(BootstrapParameters)))
  
  
  
  test<- glm(formula(mod),family=binomial, data=data_sub)
  x1<-model.matrix(test)[,BS_idx]%*%coef(mod)[BS_idx]
  
  
  BootstrapCoefs<- BootstrapParameters[,c(1, BS_idx)]
  Basis<- cbind(rep(1,10), bs(JdateForPlotting))
  
  RealFit<- Basis%*%coef(mod)[c(1, BS_idx)]
  RealFitCenter1<- RealFit-mean(x1)-coef(mod)[1]
  BootstrapFits<- Basis%*%t(BootstrapCoefs)
  quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
  
  cis<-apply(BootstrapFits, 1, quant.func)
  MinimumYlim1<- min(cis-mean(x1)-coef(mod)[1])
  MaximumYlim1<- max(cis-mean(x1)-coef(mod)[1])
  cil1<-cis[1,]-mean(x1)-coef(mod)[1]
  ciu1<-cis[2,]-mean(x1)-coef(mod)[1]
  
  
  fitdf_jdate=data.frame(x=JdateForPlotting, y=inv.logit(RealFitCenter1), LCI=inv.logit(cil1), UCI=inv.logit(ciu1))  
  fitdf_jdate$DummyDate=as.Date(JdateForPlotting, origin=as.Date("2013-01-01"))
  
  
  
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
  
  SD_idx=c(1, which(grepl("Shore", colnames(BootstrapParameters))))
  
  coeffac <- c(1,grep("ShoreDist", colnames(model.matrix(mod))))
  coefradial <- c(grep("LocalRadialFunction", colnames(model.matrix(mod))))
  coefpos <- coeffac[which(is.na(match(coeffac, coefradial)))]
  xvals <- data_sub[, which(names(data_sub) == "ShoreDist")]
  newX <- sort(unique(xvals))
  newX <- newX[2:length(newX)]
  partialfit <- coef(mod)[c(coefpos)]
  rcoefs <- NULL
  try(rcoefs <- rmvnorm(1000, coef(mod), summary(mod)$cov.scaled), 
      silent = T)
  if (is.null(rcoefs) || length(which(is.na(rcoefs) == T)) > 0) {
    rcoefs <- rmvnorm(1000, coef(mod), as.matrix(nearPD(summary(mod)$cov.scaled)$mat))
  }
  
  
  if(length(SD_idx)>1){
    
    rpreds <- as.data.frame(rcoefs[, c(coefpos)])
    BootstrapCoefs2=data.frame(vals=rpreds[,1])
    BootstrapCoefs2$ShoreDist=colnames(rpreds)[1]
    
    # Recompile for plotting
    for(jj in 2:length(SD_idx)){
      temp=data.frame(vals=rpreds[,jj])
      temp$ShoreDist=colnames(rpreds)[jj]
      BootstrapCoefs2=rbind(BootstrapCoefs2, temp)
      rm(temp)
    }
    
  }else{
    
    rpreds <- rcoefs[,coefpos]
    BootstrapCoefs2=data.frame(vals=rpreds)
    BootstrapCoefs2$ShoreDist=colnames(BootstrapParameters)[SD_idx]
  }
  BootstrapCoefs2$GroupId=unique(data_sub$GroupId)
  
  ggplot(BootstrapCoefs2, aes(x=ShoreDist, y=inv.logit(vals)))+geom_violin()
  
  #######################################################
  # Year as Factors #
  #######################################################
  Yr_idx=c(1,which(grepl("Year", colnames(BootstrapParameters))))
  
  
  coeffac <- c(1,grep("Year", colnames(model.matrix(mod))))
  coefradial <- c(grep("LocalRadialFunction", colnames(model.matrix(mod))))
  coefpos <- coeffac[which(is.na(match(coeffac, coefradial)))]
  xvals <- data_sub[, which(names(data_sub) == "Year")]
  newX <- sort(unique(xvals))
  newX <- newX[2:length(newX)]
  partialfit <- coef(mod)[c(coefpos)]
  rcoefs <- NULL
  try(rcoefs <- rmvnorm(1000, coef(mod), summary(mod)$cov.scaled), 
      silent = T)
  if (is.null(rcoefs) || length(which(is.na(rcoefs) == T)) > 0) {
    rcoefs <- rmvnorm(1000, coef(mod), as.matrix(nearPD(summary(mod)$cov.scaled)$mat))
  }
  
  
  if(length(Yr_idx)>1){
    
    rpreds <- as.data.frame(rcoefs[, c(coefpos)])
    BootstrapCoefs3=data.frame(vals=rpreds[,1])
    BootstrapCoefs3$Year=colnames(rpreds)[1]
    
    ############################
    # Recompile for plotting #
    #############################
    
    
    for(jj in 2:ncol(rpreds)){
      temp=data.frame(vals=rpreds[,jj])
      temp$Year=colnames(rpreds)[jj]
      BootstrapCoefs3=rbind(BootstrapCoefs3, temp)
      rm(temp)
    }
    
  }else{
    
    rpreds <- rcoefs[,coefpos]
    BootstrapCoefs3=data.frame(vals=rpreds)
    BootstrapCoefs3$Year=colnames(model.matrix(mod))[Yr_idx]
  }
  
  BootstrapCoefs3$GroupId=unique(data_sub$GroupId)
  
   Yr_P[[ii]]=ggplot(BootstrapCoefs3, aes(x=Year, y=vals)) +
     geom_boxplot() +
     theme_minimal()+
     ggtitle(paste('Partial Plot of Years', as.character(unique(data_sub$GroupId))))


  
  # Add all data for plotting
  fitdf_jdate$GroupId=data_sub$GroupId[1]
  
  
  if(ii==1){
    fitdf_jdate_out=fitdf_jdate
    fitdf_shoredist_out=BootstrapCoefs2
    fitdf_Year_out=BootstrapCoefs3
    
  }else {
    fitdf_jdate_out=rbind(fitdf_jdate_out, fitdf_jdate)
    fitdf_shoredist_out=rbind(fitdf_shoredist_out, BootstrapCoefs2)
    fitdf_Year_out=rbind(fitdf_Year_out, BootstrapCoefs3)
  }
  rm(data_sub, mod, JdateForPlotting, x1, test, BootstrapCoefs, Basis, BootstrapFits, BootstrapParameters)
  
}


# Visualise partial plots with data


mm=data.frame(aggregate(data=OccTable_daily, BBOcc~DayBin+GroupId, FUN=mean))
mm$med=aggregate(data=OccTable_daily, JulienDay~DayBin+GroupId, FUN=median)[,3]
mm$DummyDate=as.Date(mm$med, origin=as.Date("2013-01-01"))
OccTable_daily$DummyDate=as.Date(OccTable_daily$JulienDay, origin=as.Date("2013-01-01"))


ggplot(data=fitdf_jdate_out) +
  theme_bw() +
  facet_wrap(~GroupId) +
  scale_colour_manual(values=cbbPalette) +
  geom_line(aes(DummyDate, y), size=1) +
  geom_ribbon(aes(x=DummyDate, ymin=LCI, ymax=UCI),
              alpha=.2,linetype= 'blank') +
  geom_rug(data=subset(OccTable_daily, OccAll==1), 
           aes(x=DummyDate, y=OccAll*.3),
           sides='t', alpha=.8) +
  geom_rug(data=subset(OccTable_daily, OccAll==0), 
           aes(x=DummyDate, y=OccAll),
           sides='b', alpha=.8) 

ggplot(data=fitdf_shoredist_out) +
  theme_bw() +
  facet_wrap(~GroupId) +
  geom_boxplot(aes(x=ShoreDist, y=(vals)))

ggplot(data=fitdf_Year_out) +
  theme_bw() +
  facet_wrap(~GroupId) +
  geom_boxplot(aes(x=Year, y=(vals)))

