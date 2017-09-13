# 1)  Load packages and declare local functions      ##############

# This code investigates the various models that look at the different occupancy distributions
# incorporated by the data 
# rm(list=ls())
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
library(car)             # for big A anova (re Dr. Nora V Carlson )
library(viridis)




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



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
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


# Create grouping factors for group + is cro
OccTable$GroupIsCro=as.factor(paste(OccTable$GroupId, OccTable$IsCro))

id1_cro=as.numeric(as.factor(OccTable$Date)) ## used to make id4 but not used in correlation
id2_cro=as.numeric(as.factor(OccTable$UnitLoc)) ## used to make id4 but not used in correlation
id3_cro=as.numeric(as.factor(OccTable$HourAfterPeakSolEle))
id4_cro=as.numeric(paste(id1_cro, id2_cro, sep=""))  


OccTable$id3_cro=id3_cro
OccTable$id4_cro=id4_cro
id4_cro[OccTable$UnitLoc=='Cro_05']=max(id4_cro)+as.numeric(id1_cro[OccTable$UnitLoc=='Cro_05'])  


OccTable$id3_cro=id3_cro
OccTable$id4_cro=id4_cro

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

temp=gamm(BNDTotOffset~s(HourAfterPeakSolEle, bs='cc', k=10, by=IsCro)+IsCro, 
                     correlation=corAR1(form= ~id3|id4),
                     data=OccTable_DPD,
                     family=binomial,
                     random=list(UnitLoc=~1))







FullModel=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc', k=10)+ 
                 s(HourAfterPeakSolEle, bs='cc', k=10, by=IsCro)+IsCro,
           correlation=corAR1(form= ~id3|id4),
             data=OccTable_DPD_nocro,
             family=binomial,
             random=list(UnitLoc=~1))

FullModel_1=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ 
                   s(HourAfterPeakSolEle, bs='cc', k=10, by=GroupId)+GroupId,
               correlation=corAR1(form=~id3|id4),
               data=OccTable_DPD_nocro,
               family=binomial,
               random=list(UnitLoc=~1))

FullModel_2=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10) +
                   s(HourAfterPeakSolEle, bs='cc', k=10, by=ShoreDist)+ShoreDist,
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

FullModel_12=gamm(BNDTotOffset ~ s(HourAfterPeakSolEle, bs='cc', k=10, by=GroupId)+GroupId,
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))

FullModel_13=gamm(BNDTotOffset ~ s(HourAfterPeakSolEle, bs='cc', k=10),
                  correlation=corAR1(form=~id3|id4),
                  data=OccTable_DPD_nocro,
                  family=binomial,
                  random=list(UnitLoc=~1))

FullModel_14=gamm(BNDTotOffset ~ s(HourAfterPeakSolEle, bs='cc', k=10, by=GroupId)+GroupId,
                  correlation=corAR1(form=~id3|id4),
                  data=OccTable_DPD_nocro,
                  family=binomial,
                  random=list(UnitLoc=~1))

FullModel_15=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ 
                   s(HourAfterPeakSolEle, bs='cc', k=10),
                 correlation=corAR1(form=~id3|id4),
                 data=OccTable_DPD_nocro,
                 family=binomial,
                 random=list(UnitLoc=~1))

FullModel_15=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ s(Depth_m)
                    s(HourAfterPeakSolEle, bs='cc', k=10),
                  correlation=corAR1(form=~id3|id4),
                  data=OccTable_DPD_nocro,
                  family=binomial,
                  random=list(UnitLoc=~1))


FullModel_11=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc' , k=10)+ 
                    s(HourAfterPeakSolEle, bs='cc', k=10, by=ShoreDist) + Depth_m +DistToSalmonRun,
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



FullModel_12=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc',by=IsCro, k=10)+
                    s(HourAfterPeakSolEle, bs='cc', k=10, by=IsCro),
                  correlation=corAR1(form=~id3|id4),
                  data=OccTable_DPD,
                  family=binomial,
                  random=list(UnitLoc=~1))


FullModel_13=gamm(BNDTotOffset ~s(HourAfterHigh, bs='cc',by=GroupIsCro, k=3)+
                    s(HourAfterPeakSolEle, bs='cc',by=GroupIsCro, k=3),
                  correlation=corAR1(form=~id3_cro|id4_cro),
                  data=OccTable_DPD,
                  family=binomial,
                  random=list(UnitLoc=~1))




mm=data.frame(AIC( FullModel_1, FullModel_2, FullModel_3, FullModel_4, 
                  FullModel_5, FullModel_6, FullModel_7, FullModel_9, 
                  FullModel_10, FullModel_11, FullModel_12, FullModel_13,
                  FullModel_14, FullModel_15))

mm=mm[order(mm$AIC),]

modlist=list( FullModel_1, FullModel_2, FullModel_3, FullModel_4, 
              FullModel_5, FullModel_6, FullModel_7, FullModel_9, 
              FullModel_10, )

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







mod1=gamm(BNDTotOffset ~ s(HourAfterPeakSolEle, bs='cc') + s(HourAfterHigh,  bs='cc'),
          correlation=corAR1(form=~id3|id4),
          data=OccTable_DPDtem_nocro,
          family=binomial,
          random=list(UnitLoc=~1))


mod2=gamm(BNDTotOffset ~ ShoreDist + s(HourAfterPeakSolEle, bs='cc') + s(HourAfterHigh,  bs='cc', by = ShoreDist) 
          correlation=corAR1(form=~id3|id4),
          data=OccTable_DPD_nocro,
          family=binomial,
          random=list(UnitLoc=~1))

mod3=gamm(BNDTotOffset ~ GroupId + s(HourAfterPeakSolEle,  bs = 'cc',  by = GroupId) + s(HourAfterHigh,  bs = 'cc'),
          correlation=corAR1(form=~id3|id4),
          data=OccTable_DPD_nocro,
          family=binomial,
          random=list(UnitLoc=~1))

mod4=gamm(BNDTotOffset ~ te(HourAfterPeakSolEle,  HourAfterHigh,  by=ShoreDist,  bs='cc'),
          correlation=corAR1(form=~id3|id4),
          data=OccTable_DPD_nocro,
          family=binomial,
          random=list(UnitLoc=~1))

mod5=gamm(BNDTotOffset ~ ShoreDist + te(HourAfterPeakSolEle,  HourAfterHigh,  by=ShoreDist,  bs='cc'),
          correlation=corAR1(form=~id3|id4),
          data=OccTable_DPD_nocro,
          family=binomial,
          random=list(UnitLoc=~1))




# Cromarty Subet Modelling#####################################################################################################

OccTable_DPD_cro=OccTable_DPD[OccTable_DPD$UnitLoc =='Cro_05',]
OccTable_DPD_cro=droplevels(OccTable_DPD_cro)

# Correlation structure

id1_cro=as.numeric(as.factor(OccTable_DPD_cro$Date)) ## used to make id4 but not used in correlation
id3_cro=as.numeric(as.factor(OccTable_DPD_cro$HourAfterPeakSolEle))

OccTable_DPD_cro$id1_cro=id1_cro
OccTable_DPD_cro$id3_cro=id3_cro


Cro_mod1=gamm(BNDTotOffset~s(HourAfterPeakSolEle, bs='cc') + s(HourAfterHigh,  bs='cc'),
              correlation=corAR1(form=~id3_cro|id1_cro),
              data=OccTable_DPD_cro,
              family=binomial)

Cro_mod2=gamm(BNDTotOffset ~  te(HourAfterPeakSolEle,  HourAfterHigh,  bs='cc'),
              correlation=corAR1(form=~id3_cro|id1_cro),
              data=OccTable_DPD_cro,
              family=binomial)





# 

# Do some plotting ###############

# Set up the predictions dataframe
PredDat=expand.grid(GroupId=unique(OccTable$GroupId),
                    ShoreDist=unique(OccTable$ShoreDist),
                    HourAfterPeakSolEle=unique(OccTable$HourAfterPeakSolEle))

PredDat$HourAfterHigh=median(OccTable$HourAfterHigh)
PredDat$UnitLoc=paste(PredDat$GroupId, PredDat$ShoreDist, sep = '_')

PredDat$UnitLoc=factor(PredDat$UnitLoc, levels=level_names_unit)

# Add actual data
Datapoints=aggregate(data=OccTable_DPD, BNDTotOffset ~ UnitLoc+HourAfterPeakSolEle, FUN=mean)
# Merge with prediction data
PredDat=merge(Datapoints, PredDat, all.y=TRUE, by=c('UnitLoc','HourAfterPeakSolEle'))


PredDat_nocro=subset(PredDat, UnitLoc != 'Cro_05')
PredDat_cro=subset(PredDat, UnitLoc == 'Cro_05')


# Non-cromarty predictions
NoCroFit=as.data.frame(predict(mod3, PredDat_nocro, se.fit=TRUE))
PredDat_nocro$fit=inv.logit(NoCroFit$fit)
PredDat_nocro$UCI=inv.logit(NoCroFit$fit+(1.96*NoCroFit$se.fit))
PredDat_nocro$LCI=inv.logit(NoCroFit$fit-(1.96*NoCroFit$se.fit))

# Cromarty predictions
CroFit=as.data.frame(predict(Cro_mod2, PredDat_cro, se.fit=TRUE))
PredDat_cro$fit=inv.logit(CroFit$fit)
PredDat_cro$UCI=inv.logit(CroFit$fit+(1.96*CroFit$se.fit))
PredDat_cro$LCI=inv.logit(CroFit$fit-(1.96*CroFit$se.fit))

# FUUUUUUUUUUUUUUZE them. Yes, that's an Invader Zim reference.




# Create figure of predictions and data
png(filename = paste('HourlyDet.png'),
    units="in",
    width=7,
    height=9,
    pointsize=12,res = 400)

ggplot(PredDat_nocro) +
  facet_wrap(~GroupId) +
  geom_point(aes(HourAfterPeakSolEle, BNDTotOffset, col=ShoreDist), size=.5)+
  geom_line(aes(HourAfterPeakSolEle, fit)) +
  geom_ribbon(aes(x=HourAfterPeakSolEle, ymin=LCI, ymax=UCI),alpha=.3,linetype= 'blank') +
  geom_point(data=PredDat_cro,aes(HourAfterPeakSolEle, BNDTotOffset, col=ShoreDist), size=.5)+
  geom_line(data=PredDat_cro,aes(HourAfterPeakSolEle, fit)) +
  geom_ribbon(data=PredDat_cro,aes(x=HourAfterPeakSolEle, ymin=LCI, ymax=UCI),alpha=.3,linetype= 'blank')+
  theme_bw() +
  ylim(c(0,.2)) +
  scale_colour_manual(values=cbbPalette) +
  ggtitle('Hourly Detection Rate') + 
  xlab('Hour Relative to Solar Noon') +
  ylab('Broadband Clcik Train Encounter Rate') +
  theme(plot.title = element_text(hjust = 0.5))
  
dev.off()

# Figure for Cromarty Model time and Tide #####

Cro_Model_data=expand.grid(HourAfterPeakSolEle=unique(OccTable$HourAfterPeakSolEle),
                      HourAfterHigh=unique(OccTable$HourAfterHigh))
Cro_Model_data$UnitLoc='Cro_05'
Cro_Model_data$BNDTotOffset=0

Cro_Model_fit=predict(Cro_mod2, Cro_Model_data, se.fit=TRUE)
Cro_Model_data$fit=inv.logit(Cro_Model_fit$fit)
Cro_Model_data$LCI=inv.logit(Cro_Model_fit$fit-(1.96*Cro_Model_fit$se.fit))
Cro_Model_data$UCI=inv.logit(Cro_Model_fit$fit+(1.96*Cro_Model_fit$se.fit))




# Create figure of predictions and data
png(filename = paste('CromartyModel.png'),
    units="in",
    width=10,
    height=9,
    pointsize=12,res = 400)

p1<-ggplot(data=Cro_Model_data, aes(HourAfterPeakSolEle, HourAfterHigh))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))  +
  ggtitle('Model Fit') +
  geom_tile(aes(fill = fit)) +
  #scale_fill_distiller(palette = "Spectral") +
  scale_fill_viridis(option ='inferno') +
  ylab('Hour Relative to Hight Tide') +
  xlab('') +
  geom_point(data=subset(OccTable_DPD_cro,BNDTotOffset>0) ,
             aes(HourAfterPeakSolEle, HourAfterHigh, color=BNDTotOffset)) +
   #theme(legend.position='none') +
  scale_color_continuous(low = 'black', high = 'white') +
  guides(color = FALSE, fill=FALSE, size = FALSE)

p3<-ggplot(data=Cro_Model_data, aes(HourAfterPeakSolEle, HourAfterHigh))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))  +
  ggtitle('Model Upper 95% Confidence Interval') +
  geom_tile(aes(fill = UCI)) +
  #scale_fill_distiller(palette = "Spectral") +
  scale_fill_viridis(option ='inferno') +
  ylab('') +
  xlab('') +
  geom_point(data=subset(OccTable_DPD_cro,BNDTotOffset>0) ,
             aes(HourAfterPeakSolEle, HourAfterHigh, color=BNDTotOffset)) +
  #theme(legend.position='none') +
  scale_color_continuous(low = 'black', high = 'white') +
  guides(color = FALSE, fill=FALSE, size = FALSE)


p2<-ggplot(data=Cro_Model_data, aes(HourAfterPeakSolEle, HourAfterHigh))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))  +
  ggtitle('Model Lower 95%  Confidence Interval') +
  geom_tile(aes(fill = UCI)) +
  #scale_fill_distiller(palette = "Spectral") +
  scale_fill_viridis(option ='inferno') +
  ylab('Hour Relative to Hight Tide') +
  xlab('Hour Relative to Solar Noon') +
  geom_point(data=subset(OccTable_DPD_cro,BNDTotOffset>0) ,
             aes(HourAfterPeakSolEle, HourAfterHigh, color=BNDTotOffset)) +
  #theme(legend.position='none') +
  scale_color_continuous(low = 'black', high = 'white') +
  guides(color = FALSE, size = FALSE)

  
 multiplot(p1, p2, p3,cols = 2)

dev.off()





