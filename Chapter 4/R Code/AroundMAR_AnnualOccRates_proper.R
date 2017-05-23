#########################################################################
# Descriptive Statistics For Occupancy#
#########################################################################



#############
# Setup     #
#############
# This code investigates the various models that look at the different occupancy distributions
# incorporated by the data 
rm(list=ls())
library(boot)            # for inv.logit
library(mgcv)
library(ggplot2)
library(lme4)
library(dplyr)           # for distinct function 
library(geepack)
library(splines)
library(RColorBrewer)
library(MuMIn)           # for QIC
library(MASS)            # for mvrnorm in boostrapping intervals 
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(mvtnorm)         # for rmvnorm used in predictions/plotting

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
  
  pr <- predict(mod, data_sub, type="response")                          # the final model is used to predict the data on the response scale (i.e. a value between 0 and 1)
  pred <- prediction(pr, data_sub$BBOcc)                                    # to specify the vector of predictions (pr) and the vector of labels (i.e. the observed values "Pres")
  perf <- performance(pred, measure="tpr", x.measure="fpr")          # to assess model performance in the form of the true positive rate and the false positive rate
  plot(perf, colorize=TRUE), print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5)) # to plot the ROC curve
  
  
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
# Data Exploration/Viz #
##############################################################################


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






#################################################################################
# Different Spline for all Locs #
#################################################################################


OccTable_daily$YrUnitLoc=paste(OccTable_daily$UnitLoc, OccTable_daily$Year, sep='')

OccTable_daily_wDetections= OccTable_daily[OccTable_daily$YrUnitLoc %in%
                                             names(which(tapply(OccTable_daily$OccAll,OccTable_daily$YrUnitLoc,sum)>1)),]


OccTable_daily_wDetections=droplevels(OccTable_daily_wDetections)

# Table to store model performance #
ModelTable=data.frame(DeplotymentLoc=unique(OccTable$GroupId))
ModelTable$ModelFormula='none'
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


# list to store the models #
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
    
    ##########################################################################################
    # First determine whether linear model or spline for Julien Day #
    #########################################################################################
    
    null=geeglm(OccAll~ShoreDist+Year, 
                         corstr = 'ar1', 
                         offset = BNDTotOffset, 
                         family = binomial, 
                         id     = UnitLoc, 
                         data   = data_sub) 
    
    mod1=geeglm(OccAll~ShoreDist+Year+bs(JulienDay, knots = mean(JulienDay)), 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub) 
    
    mod2=geeglm(OccAll~ShoreDist+Year+JulienDay, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub) 
    
    QICvals=QIC(null, mod1, mod2)
    
    if (sum(QICvals>10e+15)>2){
   
    null=geeglm(OccAll~ShoreDist, 
                  corstr = 'ar1', 
                  offset = BNDTotOffset, 
                  family = binomial, 
                  id     = UnitLoc, 
                  data   = data_sub) 
    
    mod1=geeglm(OccAll~ShoreDist+bs(JulienDay, knots = mean(JulienDay)), 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub) 
    
    mod2=geeglm(OccAll~ShoreDist+JulienDay, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
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
    
    # At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
    # (the ones that, if removed, determine the biggest increase in QIC enter the model first). # Pilfered from Priotta Sperm Whales 
    mod1=geeglm(as.formula(paste('OccAll~', JdateForm, '+ ShoreDist + Year')), 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub) 
    
    mod2=geeglm(as.formula(paste('OccAll~', JdateForm, '+ ShoreDist')), 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, 
                id     = UnitLoc, 
                data   = data_sub)   
    
    mod3=geeglm(as.formula(paste('OccAll~', JdateForm,'+ Year')), 
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

    
    # Also compare model with linear smooth
    
    Qicdf=data.frame(QIC(mod1, mod2, mod3, mod4))
    Qicdf$deltaQIC=abs(Qicdf$QIC-Qicdf$QIC[1])
    Qicdf$Varnames=c('All',  'Year', 'ShoreDist' , JdateForm)
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
                           data   = data_sub)
      }
    
    # Aggregate data for plotting
    OneYearAggs=data.frame(aggregate(data=subset(data_sub, Year==unique(newdat_perdOnly$Year)),
                                                      BBOcc~DayBin+GroupId+ShoreDist, FUN=mean))
    OneYearAggs=cbind(OneYearAggs,
                      aggregate(data=subset(data_sub, Year==unique(newdat_perdOnly$Year)), 
                                JulienDay~DayBin+GroupId+ShoreDist, FUN=median)[,4])
   
    OneYearAggs=cbind(OneYearAggs,
                      aggregate(data=subset(data_sub, Year==unique(newdat_perdOnly$Year)),
                                     BNDTotOffset~DayBin+GroupId+ShoreDist, FUN=mean))
    
    
    colnames(OneYearAggs)[5]=c('med')
    OneYearAggs$JulienDay=OneYearAggs$med
    OneYearAggs$Year=unique(newdat_perdOnly$Year)
    
  if(ii==1){
    fit=cbind(newdat,  predictvcv(modlist[[ii]]))
    dummyfit=cbind(newdat_perdOnly, predictvcv(modlist[[ii]], newdata = newdat_perdOnly))
    AggData=cbind(OneYearAggs, predictvcv(modlist[[ii]], newdata = OneYearAggs))
    
    
  }else {
    fit=rbind(fit, cbind(newdat, predictvcv(modlist[[ii]])) )
    dummyfit=rbind(dummyfit,cbind(newdat_perdOnly, predictvcv(modlist[[ii]], newdata = newdat_perdOnly)))
    AggData=rbind(AggData, cbind(OneYearAggs, predictvcv(modlist[[ii]], newdata = OneYearAggs)))
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





# Add dummy dates for plotting (to fix the X axis)
AggData$DummyDate=as.Date(AggData$med, origin=as.Date("2013-01-01"))
AggData$UnitLoc=paste(AggData$GroupId, AggData$ShoreDist, sep="_")
AggData$BBEst=AggData$BNDTotOffset*AggData$OccAll

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
  


##################################################################
# Get the Partial Plots for each plot Group ID               #
##################################################################

# Code based on MUMIN run.partials Scott-Hayward 2015
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


# Plot Storage
p=list()
Sd_P=list()
Yr_P=list()


for(ii in 1:10){
  
  # Again subset the data
  data_sub=subset(OccTable_daily_wDetections, GroupId==unique(OccTable$GroupId)[ii])
  data_sub$ShoreDist=factor(data_sub$ShoreDist, levels=c('05', '10', '15'))
  data_sub <- droplevels(data_sub)
  
  # Add a dummy date for plotting
  data_sub$DummyDate=as.Date(data_sub$JulienDay, origin=as.Date("2013-01-01"))
  
  # Extract model (just for clarity)
  mod=modlist[[ii]]
  
  #######################################################
  # Julien Date Smoothes #
  #######################################################
  fitdf_jdate=partialDF(mod = mod, data = data_sub, Variable = 'JulienDay')
  
  fitdf_jdate$DummyDate=as.Date(fitdf_jdate$JulienDay, origin=as.Date("2013-01-01"))
  
  
  
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
  
  BootstrapCoefs2=partialdf_factor(mod = mod, data = data_sub, variable = 'ShoreDist')
  BootstrapCoefs2$GroupId=data_sub$GroupId[1]
  
  ggplot(BootstrapCoefs2, aes(x=ShoreDist, y=inv.logit(vals)))+geom_violin()
  
  #######################################################
  # Year as Factors #
  #######################################################
  BootstrapCoefs3 = partialdf_factor(mod = mod, data = data_sub, variable = 'Year')
  BootstrapCoefs3$GroupId=data_sub$GroupId[1]
  
  ggplot(BootstrapCoefs3, aes(x=Year, y=inv.logit(vals)))+geom_violin()
  
  
   Yr_P[[ii]]=ggplot(BootstrapCoefs3, aes(x=Year, y=vals)) +
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
    fitdf_shoredist_out=BootstrapCoefs2
    fitdf_Year_out=BootstrapCoefs3
    
    AggData=OneYearAggs
    
  }else {
    fitdf_jdate_out=rbind(fitdf_jdate_out, fitdf_jdate)
    fitdf_shoredist_out=rbind(fitdf_shoredist_out, BootstrapCoefs2)
    fitdf_Year_out=rbind(fitdf_Year_out, BootstrapCoefs3)
    
    AggData=rbind(AggData, OneYearAggs)
  }
  rm(data_sub, mod, JdateForPlotting, x1, test, BootstrapCoefs, Basis, BootstrapFits, BootstrapParameters)
  
}


# Visualise partial plots with data


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
  geom_rug(data=subset(OccTable_daily_wDetections, BBOcc==1), 
           aes(x=DummyDate, y=BBOcc*.3),
           sides='t', alpha=.8) +
  geom_rug(data=subset(OccTable_daily_wDetections, BBOcc==0), 
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



################################################################################
# Occupancy rate seems to be correlated with AUC scores indicating that the models need
# more detections in order to fit
# Run Mann-Whiteney U test to determine if that's the case
###############################################################################

Group_OccRates=aggregate(data=OccTable_daily_wDetections, BBOcc~GroupId, FUN=mean)
Group_OccRates$OccOrder=order(Group_OccRates$BBOcc)
Group_OccRates$AucOrder=order(ModelTable$AUC)

x=c(ModelTable$AUC, Group_OccRates$BBOcc)
g=as.factor(c(rep('AUC',10),rep('OccRate',10)))

wilcox.test(x ~ g) 

