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



OccTable_daily_ModelVars=OccTable[,c('BBOcc', 'GroupId', 'IsCroFactor', 'Hr', 'elevation', 'ShoreDist', 'Z',
                                     'BNDTotOffset', 'Date')]


# fgeeMod=geeglm(BBOcc~GroupId+IsCroFactor+ShoreDist+
#                  IsCroFactor*bs(Hr)+ShoreDist*bs(Hr)+GroupId*bs(Hr)+
#                  IsCroFactor*bs(elevation)+ShoreDist*bs(elevation)+GroupId*bs(elevation)+
#                  IsCroFactor*bs(Z)+ShoreDist*bs(Z)+GroupId*bs(Z),
#                corstr = 'ar1', 
#                family = binomial, # leave out constrains
#                id=GroupId:Date, 
#                offset = BNDTotOffset, 
#                data = OccTable_daily_ModelVars)

# fullmod=geeglm(BBOcc~GroupId+IsCroFactor+ShoreDist+Z+elevation,
#                corstr = 'ar1', 
#                family = binomial, # leave out constrains
#                id=GroupId:Date, 
#                offset = BNDTotOffset, 
#                data = OccTable_daily_ModelVars)

# # Test for autocorrelated residuals
# library(lawstat)
# library(MRSea)
# runs.test(residuals(fgeeMod, type = "pearson"),alternative = c("two.sided"))

# The small p-value (p << 0:05) indicates that there is an issue with cor-
# relation in the residuals

# OccTable_daily_ModelVars$blockid <- paste(OccTable_daily_ModelVars$GroupId, OccTable_daily_ModelVars$Date, sep = "")
# 
# runACF(OccTable_daily_ModelVars$blockid, fgeeMod, store = F)
#  
# # Correlation in all blocks declines to approximately zero, which we want to
# # see if our blocking is appropriate
# 
# splineParams<-makesplineParams(data=OccTable_daily_ModelVars,
#                                varlist=c('Hr', 'elevation', 'Z'),
#                                predictionData=OccTable_daily_ModelVars)
# 
# 
# # Plotting Cumulative Residauals
# plotCumRes(fullmod, varlist= c("Hr"), splineParams=splineParams)
# plotCumRes(fullmod, varlist= c("elevation"), splineParams=splineParams)
# plotCumRes(fullmod, varlist= c("Z"), splineParams=splineParams)
# 
# ##


mod1<- geeglm(BBOcc ~Year+ ShoreDist + GroupId + 
                GroupId*bs(elevation) + ShoreDist*bs(elevation),# + IsCroFactor*bs(elevation),
              corstr = 'ar1', 
              family = binomial, # leave out constrains
              id=GroupId:Date, 
              offset = BNDTotOffset, 
              data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# tie
mod2<- geeglm(BBOcc ~Year+ ShoreDist + GroupId + 
                GroupId*bs(HourAfterPeakSolEle) + ShoreDist*bs(HourAfterPeakSolEle),# + IsCroFactor*bs(HourAfterPeakSolEle),
              corstr = 'ar1', 
              family = binomial, # leave out constrains
              id=GroupId:Date, 
              offset = BNDTotOffset, 
              data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])
# tie win
mod3<- geeglm(BBOcc ~Year+ ShoreDist + GroupId + 
                GroupId*bs(Hr) + ShoreDist*bs(Hr),# + IsCroFactor*bs(Hr),
              corstr = 'ar1', 
              family = binomial, # leave out constrains
              id=GroupId:Date, 
              offset = BNDTotOffset, 
              data = OccTable_DPD[OccTable_DPD$UnitLoc!='Cro_05',])


# Tide Effects

# win
mod4<- geeglm(BBOcc ~Year+ ShoreDist + GroupId + 
                GroupId*bs(HourAfterHigh) + ShoreDist*bs(HourAfterHigh),# + IsCroFactor*bs(HourAfterHigh),
              corstr = 'ar1', 
              family = binomial, # leave out constrains
              id=GroupId:Date, 
              offset = BNDTotOffset, 
              data = OccTable_DPD[OccTable_DPD$UnitLoc !='Cro_05',])

#nope
mod5<- geeglm(BBOcc ~Year+ ShoreDist + GroupId + 
                GroupId*Phase + ShoreDist*Phase,# + IsCroFactor*Phase,
              corstr = 'ar1', 
              family = binomial, # leave out constrains
              id=GroupId:Date, 
              offset = BNDTotOffset, 
              data = OccTable_DPD[OccTable_DPD$UnitLoc !='Cro_05',])

mod6<- geeglm(BBOcc ~Year+ ShoreDist + GroupId + 
                GroupId*bs(Z) + ShoreDist*bs(Z), #+ IsCroFactor*bs(Z),
              corstr = 'ar1', 
              family = binomial, # leave out constrains
              id=GroupId:Date, 
              offset = BNDTotOffset, 
              data = OccTable_DPD[OccTable_DPD$UnitLoc !='Cro_05',])

###############################################################

mod7<- geeglm(BBOcc ~Year+ ShoreDist + GroupId + 
                GroupId*bs(HourAfterPeakSolEle) + ShoreDist*bs(HourAfterPeakSolEle) +
                GroupId*bs(Z) + ShoreDist*bs(Z),
              corstr = 'ar1', 
              family = binomial, # leave out constrains
              id=GroupId:Date, 
              offset = BNDTotOffset, 
              data = OccTable_DPD[OccTable_DPD$UnitLoc !='Cro_05',])


# test model seven for autocorrelation
runs.test(residuals(mod7, type = "pearson"),alternative = c("two.sided"))


# Model names

QIC_scores=as.data.frame(QIC(mod1, mod2,
                             mod3,
                             mod4,
                             mod6,
                             mod7))

QIC_scores$Model=rownames(QIC_scores)
rownames(QIC_scores)=seq(1, nrow(QIC_scores))
QIC_scores=rbind(QIC_scores, c(NaN, 'mod5'), c(NaN, 'mod10'))
QIC_scores$Formula='Unk'


# Export Pattern
for (ii in 1:nrow(QIC_scores)){QIC_scores$Formula[QIC_scores$Model==ls(pattern = 'mod')[ii]] = 
  Reduce(paste, deparse(formula(ls(pattern = 'mod')[ii])))}

QIC_scores=QIC_scores[order(QIC_scores$QIC),]
QIC_scores$DeltaQIC=c(as.numeric(QIC_scores$QIC[1:sum(!is.nan(as.numeric(QIC_scores$QIC)))])-as.numeric(QIC_scores$QIC[1]),
                                                               rep(NaN,sum(is.nan(as.numeric(QIC_scores$QIC)))))

# Export rsquared (actually not useful)
# for (ii in 1:nrow(QIC_scores)){QIC_scores$RsquaredGLMM[ii]=r.squaredGLMM(get(ls(pattern = 'mod')[ii]))[1]}

# Write the final file for model Selection
write.csv(x = QIC_scores, file = "W:/KJP PHD/4-Bayesian Habitat Use/Figures/HourlyOccupancyModelSelection.csv")  




##################################################################################
# Data Visualization #
##################################################################################


newdat<- expand.grid(
  HourAfterPeakSolEle = seq(-12,11),
  Year=2013,
  GroupId=unique(OccTable_DPD$GroupId),
  #IsCroNumeric=c(1, 2),
  ShoreDist=unique(OccTable_DPD$ShoreDist),
  Z=seq(min(OccTable_DPD$Z), max(OccTable_DPD$Z), length.out = 12),
  BBOcc=0)

newdat$IsCroFactor=ifelse(newdat$GroupId=='Cro' & newdat$ShoreDist=='05', 'Cro05', 'Other')

newdat=cbind(newdat, predictvcv(mod7, newdata = newdat))


# Add Counts as function of Hour and Is Cromarty to Model Selection
Agg_counts=as.data.frame(aggregate(data=subset(OccTable_DPD, UnitLoc !='Cro_05'), 
                     BBOcc~GroupId+Hr, FUN=mean))
Agg_counts=rbind(Agg_counts, as.data.frame(aggregate(data=subset(OccTable_DPD, UnitLoc =='Cro_05'), 
                                   BBOcc~GroupId+Hr, FUN=mean)))

Agg_counts$IsCroFactor=c(rep('Other',nrow(Agg_counts)-24),  rep('Cro05',24))
Agg_counts$BBOcc=Agg_counts$BBOc/3


# Add Counts as function of tidal Height and Is Cromarty to Model Selection

OccTable_DPD$Zcat<-cut(OccTable_DPD$Z, quantile(OccTable_DPD$Z, probs = seq(0,1, by = .15)), right=FALSE)


Agg_counts_tide=as.data.frame(aggregate(data=OccTable_DPD[!is.na(OccTable_DPD$HourAfterHigh),], 
                                        BBOcc~GroupId+HourAfterHigh, FUN=mean))

Agg_counts_tide=as.data.frame(aggregate(data=OccTable_DPD[!is.na(OccTable_DPD$HourAfterHigh),], 
                                        BBOcc~GroupId+HourAfterHigh, FUN=mean))

Agg_counts_noon=as.data.frame(aggregate(data=OccTable_DPD[!is.na(OccTable_DPD$HourAfterPeakSolEle),], 
                                        BBOcc~GroupId+HourAfterPeakSolEle+IsCroFactor, FUN=mean))

ggplot()+
  theme_bw() +
  facet_wrap(~GroupId) +
  scale_color_brewer(type = "div", palette = "Set1", direction = 1, color=IsCro05) +  
  geom_point(data = Agg_counts_tide, aes(x=HourAfterHigh, y=BBOcc)) 



ggplot()+
  theme_bw() +
  facet_wrap(~GroupId) +
 # scale_color_brewer(type = "div", palette = "Set1", direction = 1) +  
  geom_point(data = Agg_counts_noon, aes(x=HourAfterPeakSolEle, y=BBOcc, color=IsCro05)) 






ggplot(subset(OccTable_DPD, BBOcc==1), aes(x=elevation, fill=GroupId)) +
  facet_wrap(~GroupId) +
  geom_histogram(binwidth=5, alpha=.5, position="identity")

ggplot(subset(OccTable_DPD, BBOcc==1), aes(x=HourAfterPeakSolEle, fill=GroupId)) +
  facet_wrap(~GroupId) +
  #geom_histogram(binwidth=1, alpha=.5, position="identity")+
  geom_density()


library(ggplot2)
ggplot()+
  theme_bw() +
  facet_wrap(~GroupId) +
  scale_color_brewer(type = "div", palette = "Set1", direction = 1) +  
  geom_point(data = Agg_counts, aes(x=Hr, y=BBOcc, color=IsCroFactor)) +
  # geom_point(data = subset(Agg_counts, IsCroFactor=='Cro05'),
  #            aes(x=Hr, y=BBOcc)) +
  geom_line(data=newdat, aes(Hr, inv.logit(fit), colour=IsCroFactor), size=1) +
  geom_ribbon(data=newdat, aes(x=Hr, ymin=inv.logit(lwr), ymax=inv.logit(upr),
                               colour=as.factor(IsCroFactor)), alpha=.4,linetype= 'blank') +
  
  xlab('Hour of the Day') +
  ylab('Detection Probability')
  

###################################################################################################################






scores=numeric(length = length(ls(pattern = 'mod')) )
Formulas=character(length = length(ls(pattern = 'mod')) )
Names=ls(pattern = 'mod')
for (ii in 1:length(scores)){scores[ii]=QIC(get(ls(pattern = 'mod')[ii]))}
for (ii in 1:length(scores)){Formulas[ii]= Reduce(paste, deparse(formula(get(ls(pattern = 'mod')[ii]))))}

QIC_scores=data.frame(Names=Names, Scores=scores, Formulas=Formulas)
QIC_scores=QIC_scores[order(QIC_scores$Scores),]



model.sel(get(as.character(QIC_scores$Names[1])), 
          get(as.character(QIC_scores$Names[2])),
          get(as.character(QIC_scores$Names[3])),
          get(as.character(QIC_scores$Names[4])),
          get(as.character(QIC_scores$Names[5])), rank = QIC)

