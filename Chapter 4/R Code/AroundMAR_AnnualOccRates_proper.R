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

# list to store the models
# list to store the models
modlist=list()

for(ii in 1:10){
  data_sub=subset(OccTable_daily_wDetections, GroupId==unique(OccTable$GroupId)[ii])
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
  geom_point(data=mm, aes(x=DummyDate, y=BBOcc,
                          color=ShoreDist), size=.9) +
  geom_ribbon(aes(x=DummyDate, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
              alpha=.2,linetype= 'blank') 
  


# Plot the data without Cromarty 05
dummyfit$UnitLoc=paste(dummyfit$GroupId, dummyfit$ShoreDist, sep = '_')
ggplot(data=dummyfit[dummyfit$UnitLoc != 'Cro_05',]) +
  theme_bw() +
  facet_wrap(~GroupId) +
  scale_colour_manual(values=cbbPalette) +
  geom_line(aes(DummyDate, inv.logit(fit), colour=ShoreDist), size=1) +
  geom_point(data=mm[mm$UnitLoc !='Cro_05',], aes(x=DummyDate, y=BBOcc,
                                                  color=ShoreDist), size=.9) +
  geom_ribbon(aes(x=DummyDate, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
              alpha=.2,linetype= 'blank') 





# Get the Partial Residuals for each plot Group ID

# Plot Storage
p=list()
Sd_P=list()
Yr_P=list()

for(ii in 1:10){
  
  data_sub=subset(OccTable_daily_wDetections, GroupId==unique(OccTable$GroupId)[ii])
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

  
  fitdf=data.frame(x=JdateForPlotting, y=inv.logit(RealFitCenter1), LCI=inv.logit(cil1), UCI=inv.logit(ciu1))  
  fitdf$DummyDate=as.Date(JdateForPlotting, origin=as.Date("2013-01-01"))
  

  
  p[[ii]]=ggplot(data=fitdf) + 
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
  SD_idx=which(grepl("Shore", colnames(BootstrapParameters)))
  if(length(SD_idx)>1){
  x2<-model.matrix(test)[,SD_idx]%*%coef(mod)[SD_idx]
  
  BootstrapCoefs2<- BootstrapParameters[,c(1, SD_idx)]
  RealFit<- coef(mod)[c(1, SD_idx)]
  RealFitCenter1<- RealFit-mean(x2)-coef(mod)[1]
  SDs=apply(BootstrapCoefs2, 2, sd)
  sd1=RealFitCenter1-SDs
  sd2=RealFitCenter1+SDs
  
  cis=apply(BootstrapCoefs2, 2, quant.func)

  MinimumYlim1<- min(cis-mean(x2)-coef(mod)[1])
  MaximumYlim1<- max(cis-mean(x2)-coef(mod)[1])
  cil1<-cis[1,]-mean(x2)-coef(mod)[1]
  ciu1<-cis[2,]-mean(x2)-coef(mod)[1]

  
  fitdf=data.frame(x=ShoreDistForPlotting, y=inv.logit(RealFitCenter1), 
                   LCI=inv.logit(cil1), UCI=inv.logit(ciu1),
                   lsd=inv.logit(sd1), usd=inv.logit(sd2))  
  
    BootstrapCoefs2_invlogit=data.frame(vals=inv.logit(c(BootstrapCoefs2[,1],
                                               BootstrapCoefs2[,2],
                                               BootstrapCoefs2[,3])))
    
    BootstrapCoefs2_invlogit$ShoreDist=as.character(sort(rep(c(5,10,15), length.out=nrow(BootstrapCoefs2_invlogit))))
    
    
    Sd_P[[ii]]=ggplot(BootstrapCoefs2_invlogit, aes(x=ShoreDist, y=vals)) + 
      geom_boxplot() + 
      theme_minimal() +
      ggtitle(paste('Partial Plot of ShoreDist', as.character(unique(data_sub$GroupId))))}

 
  #######################################################
  # Year as Factors #
  #######################################################
  Yr_idx=which(grepl("Year", colnames(BootstrapParameters)))
  
  if(length(Yr_idx)>1){
    x2<-model.matrix(test)[,Yr_idx]%*%coef(mod)[Yr_idx]
    BootstrapCoefs2<- BootstrapParameters[,c(1, Yr_idx)]
    RealFit<- coef(mod)[c(1, Yr_idx)]
    RealFitCenter1<- RealFit-mean(x2)-coef(mod)[1]
    SDs=apply(BootstrapCoefs2, 2, sd)
    sd1=RealFitCenter1-SDs
    sd2=RealFitCenter1+SDs
    
    cis=apply(BootstrapCoefs2, 2, quant.func)
    
    MinimumYlim1<- min(cis-mean(x2)-coef(mod)[1])
    MaximumYlim1<- max(cis-mean(x2)-coef(mod)[1])
    cil1<-cis[1,]-mean(x2)-coef(mod)[1]
    ciu1<-cis[2,]-mean(x2)-coef(mod)[1]
    
    
    fitdf=data.frame(x=YearsForPlotting, y=inv.logit(RealFitCenter1), 
                     LCI=inv.logit(cil1), UCI=inv.logit(ciu1),
                     lsd=inv.logit(sd1), usd=inv.logit(sd2))  
    
    BootstrapCoefs2_invlogit=data.frame(vals=inv.logit(c(BootstrapCoefs2[,1],
                                                         BootstrapCoefs2[,2],
                                                         BootstrapCoefs2[,3])))
    
    BootstrapCoefs2_invlogit$Year=as.character(sort(rep(YearsForPlotting, length.out=nrow(BootstrapCoefs2_invlogit))))
    
    
    Yr_P[[ii]]=ggplot(BootstrapCoefs2_invlogit, aes(x=Year, y=vals)) + 
      geom_boxplot() +
      theme_minimal()+
      ggtitle(paste('Partial Plot of Years', as.character(unique(data_sub$GroupId))))}
  

  
    rm(data_sub, mod, JdateForPlotting, x1, test, BootstrapCoefs, Basis, BootstrapFits, BootstrapParameters)
    
}




####################################
# Modelling #
####################################
OccTable_daily$IsCroFactor=ifelse(OccTable_daily$GroupId=='Cro' & OccTable_daily$ShoreDist=='05', 'Cro05', 'Other')
OccTable_daily$IsCroorStOFactor=IsCroFactor
OccTable_daily$IsCroorStOFactor[OccTable_daily$UnitLoc=='Cro_05']='Sto_05'

OccTable_daily_nocro=OccTable_daily[OccTable_daily$UnitLoc != 'Cro_05',]



mod1 <- geeglm(OccAll~bs(JulienDay), 
                      corstr = 'ar1', 
                      offset = BNDTotOffset, 
                      family = binomial, # leave out constrains
                      id=UnitLoc, 
                      data = OccTable_daily)

mod2 <- geeglm(OccAll~GroupId*bs(JulienDay), 
                      corstr = 'ar1', 
                      offset = BNDTotOffset, 
                      family = binomial, # leave out constrains
                      id=UnitLoc, 
                      data = OccTable_daily)

mod3 <- geeglm(OccAll~GroupId*bs(JulienDay)+Year, 
               corstr = 'ar1', 
               offset = BNDTotOffset, 
               family = binomial, # leave out constrains
               id=UnitLoc, 
               data = OccTable_daily)

mod4 <- geeglm(OccAll~GroupId*bs(JulienDay)+ShoreDist, 
               corstr = 'ar1', 
               offset = BNDTotOffset, 
               family = binomial, # leave out constrains
               id=UnitLoc, 
               data = OccTable_daily)

mod5 <- geeglm(OccAll~GroupId*bs(JulienDay)+ShoreDist+Year, 
               corstr = 'ar1', 
               offset = BNDTotOffset, 
               family = binomial, # leave out constrains
               id=UnitLoc, 
               data = OccTable_daily)




# nope
# mod4 <- geeglm(OccAll~UnitLoc*bs(JulienDay)+Year, 
#                corstr = 'ar1', 
#                offset = BNDTotOffset, 
#                family = binomial, # leave out constrains
#                id=UnitLoc, 
#                data = OccTable_daily)

# nope
# mod6 <- geeglm(OccAll~UnitLoc*bs(JulienDay), 
#                corstr = 'ar1', 
#                offset = BNDTotOffset, 
#                family = binomial, # leave out constrains
#                id=UnitLoc, 
#                data = OccTable_daily)

mod6 <- geeglm(OccAll~ShoreDist*bs(JulienDay), 
               corstr = 'ar1', 
               offset = BNDTotOffset, 
               family = binomial, # leave out constrains
               id=UnitLoc, 
               data = OccTable_daily)

mod7 <- geeglm(OccAll~ShoreDist*bs(JulienDay)+GroupId, 
               corstr = 'ar1', 
               offset = BNDTotOffset, 
               family = binomial, # leave out constrains
               id=UnitLoc, 
               data = OccTable_daily)

mod8 <- geeglm(OccAll~ShoreDist*bs(JulienDay)+Year, 
               corstr = 'ar1', 
               offset = BNDTotOffset, 
               family = binomial, # leave out constrains
               id=UnitLoc, 
               data = OccTable_daily)

mod9 <-geeglm(OccAll~ShoreDist*bs(JulienDay)+Year+GroupId, 
               corstr = 'ar1', 
               offset = BNDTotOffset, 
               family = binomial, # leave out constrains
               id=UnitLoc, 
               data = OccTable_daily)

##########################################################


mod10 <- geeglm(OccAll~IsCroFactor*bs(JulienDay), 
               corstr = 'ar1', 
               offset = BNDTotOffset, 
               family = binomial, # leave out constrains
               id=UnitLoc, 
               data = OccTable_daily)

mod11 <- geeglm(OccAll~IsCroFactor*bs(JulienDay)+ShoreDist*bs(JulienDay), 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, # leave out constrains
                id=UnitLoc, 
                data = OccTable_daily)

mod12 <- geeglm(OccAll~IsCroFactor*bs(JulienDay)+ShoreDist*bs(JulienDay)+Year, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, # leave out constrains
                id=UnitLoc, 
                data = OccTable_daily)

mod13 <- geeglm(OccAll~IsCroFactor*bs(JulienDay)+ShoreDist*bs(JulienDay)+GroupId, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, # leave out constrains
                id=UnitLoc, 
                data = OccTable_daily)

mod14 <- geeglm(OccAll~IsCroFactor*bs(JulienDay)+ShoreDist*bs(JulienDay)+GroupId+Year, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, # leave out constrains
                id=UnitLoc, 
                data = OccTable_daily)


mod14a <- geeglm(OccAll~ShoreDist*bs(JulienDay)+GroupId+Year+ShoreDist, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, # leave out constrains
                id=UnitLoc, 
                data = OccTable_daily_nocro)

mod14b <- geeglm(OccAll~bs(JulienDay, knots=mean(JulienDay)), 
                 random = ~1 | UnitLoc,
                 corstr = 'ar1', 
                 offset = BNDTotOffset, 
                 family = binomial, # leave out constrains
                 id=UnitLoc, 
                 data = OccTable_daily_nocro)



mod15 <- geeglm(OccAll~IsCroFactor*bs(JulienDay)+GroupId, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, # leave out constrains
                id=UnitLoc, 
                data = OccTable_daily)

mod16 <- geeglm(OccAll~IsCroFactor*bs(JulienDay)+Year, 
                corstr = 'ar1', 
                offset = BNDTotOffset, 
                family = binomial, # leave out constrains
                id=UnitLoc, 
                data = OccTable_daily)


QIC_scores=data.frame(rbind(QIC(mod1, mod2, mod3, mod4, mod5,
                          mod6, mod7, mod8, mod9, mod10,
                          mod11, mod12, mod13, mod14, mod15,
                          mod16)))
# Model 14 Wins.
QIC_scores$Model=rownames(QIC_scores)
rownames(QIC_scores)=seq(1, nrow(QIC_scores))
QIC_scores$Formula='Unk'


# Get Forumla
for (ii in 1:nrow(QIC_scores)){QIC_scores$Formula[QIC_scores$Model==ls(pattern = 'mod')[ii]] = 
  Reduce(paste, deparse(formula(ls(pattern = 'mod')[ii])))}

QIC_scores=QIC_scores[order(QIC_scores$QIC),]
QIC_scores$DeltaQIC=QIC_scores$QIC-QIC_scores$QIC[1]

# Export rsquared (actually not useful)
# QIC_scores$RsquaredGLMM=NA
# for (ii in 1:nrow(QIC_scores)){QIC_scores$RsquaredGLMM[ii]=r.squaredGLMM(get(ls(pattern = 'mod')[ii]))[1]}

# Write the final file for model Selection
write.csv(x = QIC_scores, file = "W:/KJP PHD/4-Bayesian Habitat Use/Figures/DailyOccupancyModelSelection.csv")  
  



#  
#   
#   
#   
#   newdat=expand.grid()
#   
#   newdat=expand.grid(JulienDay=seq(min(data_sub$JulienDay),
#                                    max(data_sub$JulienDay)),
#                      Year=aggregate(data=data_sub, BBOcc~Year, FUN=mean)[ which.max(aggregate(data=data_sub, BBOcc~Year, FUN=mean)[,2]),1],
#                      ShoreDist=unique(data_sub$ShoreDist),
#                      GroupId=unique(data_sub$GroupId))
#   newdat$UnitLoc=paste(newdat$GroupId, newdat$ShoreDist, sep='_')
#   newdat$OccAll=0
#   
#   ModelTable$DeplotymentLoc[ii]=(unique(data_sub$GroupId))
#   
#   # More than one year and more than one shore dist
#   if(length(unique(data_sub$ShoreDist))>1 & length(unique(data_sub$Year))>1){
#     mod=geeglm(OccAll~bs(JulienDay)+ShoreDist+Year, 
#                corstr = 'ar1', 
#                offset = BNDTotOffset, 
#                family = binomial, 
#                id     = Year, 
#                data   = data_sub)
#     
#     ModelTable$ModelFormula[ii]=Reduce(paste, deparse(formula(mod)))
#     # Else if only one year
#   }else if (length(unique(data_sub$ShoreDist))>1 & length(unique(data_sub$Year))==1){
#     mod=geeglm(OccAll~bs(JulienDay)+ShoreDist, 
#                corstr = 'ar1', 
#                offset = BNDTotOffset, 
#                family = binomial, 
#                id=Year, 
#                data = data_sub)
#     
#     ModelTable$ModelFormula[ii]=Reduce(paste, deparse(formula(mod)))
#     
#     # Else if only one Shore Dist drop shore dist but keep year
#   }else if(length(unique(data_sub$ShoreDist))==1 & length(unique(data_sub$Year))>1){
#     
#     mod=geeglm(OccAll~bs(JulienDay)+Year, 
#                corstr = 'ar1', 
#                offset = BNDTotOffset, 
#                family = binomial, 
#                id=Year, 
#                data = data_sub)
#     
#     ModelTable$ModelFormula[ii]=Reduce(paste, deparse(formula(mod)))
#     # else if only one year and one shore dist level
#   } else if (length(unique(data_sub$ShoreDist))==1 & length(unique(data_sub$Year))==1){
#     
#     mod=geeglm(OccAll~bs(JulienDay), 
#                corstr = 'ar1', 
#                offset = BNDTotOffset, 
#                family = binomial, 
#                id=Year, 
#                data = data_sub)
#     ModelTable$ModelFormula[ii]=Reduce(paste, deparse(formula(mod)))
#     
#   }
#   
# 
#   if(ii==1){
#     fit=cbind(newdat,  predictvcv(mod, newdata = newdat))
#   }else {
#     fit=rbind(fit,cbind(newdat,  predictvcv(mod, newdata = newdat)) )
#   }
#   
# }


fit$DummyDate=as.Date(fit$JulienDay, origin=as.Date("2013-01-01"))
mm$DummyDate=as.Date(mm$med, origin=as.Date("2013-01-01"))
mm$UnitLoc=paste(mm$GroupId, mm$ShoreDist, sep="_")
mm_nocro=mm[mm$UnitLoc != 'Cro_05',]


ggplot(data=fit) +
  theme_bw() +
  facet_wrap(~GroupId) +
  geom_line(aes(DummyDate, inv.logit(fit), colour=ShoreDist), size=1) +
  geom_point(data=subset(mm), aes(x=DummyDate, y=BBOcc,
                                        color=ShoreDist)) +

  geom_ribbon(aes(x=DummyDate, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
              alpha=.2,linetype= 'blank') +
  geom_point(data=subset(mm_nocro), aes(x=DummyDate, y=BBOcc,
                                        color=ShoreDist)) 






##################################################################################
# Data Visualization #
##################################################################################

newdat=OccTable_daily
newdat$IsCroFactor=ifelse(newdat$GroupId=='Cro' & newdat$ShoreDist=='05', 'Cro05', 'Other')
newdat=cbind(OccTable_daily_nocro, predictvcv(mod14a, newdata = OccTable_daily_nocro))



newdat$DummyDate=as.Date(newdat$JulienDay, origin=as.Date("2013-01-01"))
OccTable_daily$DummyDate=as.Date(OccTable_daily$JulienDay,
                                 origin=as.Date("2013-01-01"))

mm$DummyDate=as.Date(mm$med, origin=as.Date("2013-01-01"))
mm$UnitLoc=paste(mm$GroupId, mm$ShoreDist, sep="_")
mm_nocro=mm[mm$UnitLoc != 'Cro_05',]


library(ggplot2)
ggplot(data=subset(newdat,IsCroFactor=='Other'))+
  theme_bw() +
  facet_wrap(~GroupId) +
  scale_color_brewer(type = "div", palette = "Set1", direction = 1) +  
  geom_line(aes(DummyDate, inv.logit(fit), colour=ShoreDist), size=1) +
  geom_ribbon(aes(x=DummyDate, ymin=inv.logit(lwr), ymax=inv.logit(upr), color=ShoreDist),
              alpha=.3,linetype= 'blank') +
  geom_point(data=subset(mm_nocro), aes(x=DummyDate, y=BBOcc,
                                        color=ShoreDist)) 
  
  library(ggplot2)
ggplot(data=mm, aes(x=med, y=BBOcc, color=ShoreDist)) +
  theme_bw() +
  facet_wrap(~GroupId) +
  
  xlab('Julien Day') +
  ylab('Detection Probability') +
  ggtitle('BND Occupancy')



  
  
  
  
  
  
  
  
  geom_rug(data=subset(OccTable_daily, OccAll==1), 
           aes(x=DummyDate, y=OccAll*.3, colour= ShoreDist),
           sides='t', alpha=.8) +
  geom_rug(data=subset(OccTable_daily, OccAll==0), 
           aes(x=DummyDate, y=OccAll,  colour= ShoreDist),
           sides='b', alpha=.5) +
  xlab('Julien Day') +
  ylab('Detection Probability') +
  ggtitle('BND Occupancy')

  


ggplot(data=subset(newdat, ShoreDist=='05'))+
  theme_bw() +
  facet_wrap(~GroupId, nrow=2) +
  scale_color_brewer(type = "div", palette = "Set1", direction = 1, name='Year') +  
  geom_line(aes(DummyDate, inv.logit(fit), colour=as.factor(Year)), size=1) +
  geom_ribbon(aes(x=DummyDate, ymin=inv.logit(lwr), ymax=inv.logit(upr)),
              alpha=.5,linetype= 'blank') +
  geom_rug(data=subset(OccTable_daily, ShoreDist=='05' & OccAll==1), 
           aes(x=DummyDate, y=OccAll*.3, colour= as.factor(Year)),
           sides='t', alpha=.8) +
  geom_rug(data=subset(OccTable_daily, ShoreDist=='05' & OccAll==0), 
           aes(x=DummyDate, y=OccAll, colour= as.factor(Year)),
           sides='b', alpha=.5) +
  xlab('Julien Day') +
  ylab('Detection Probability') +
  ggtitle('Inshore (e.g. 5km)')

ggplot(data=subset(newdat, ShoreDist=='10'))+
  theme_bw() +
  facet_wrap(~GroupId, nrow=2) +
  scale_color_brewer(type = "div", palette = "Set1", direction = 1, name='Year') +  
  geom_line(aes(DummyDate, inv.logit(fit), colour=as.factor(Year)), size=1) +
  geom_ribbon(aes(x=DummyDate, ymin=inv.logit(lwr), ymax=inv.logit(upr)),
              alpha=.6,linetype= 'blank') +
  geom_rug(data=subset(OccTable_daily, ShoreDist=='10' & OccAll==1), 
           aes(x=DummyDate, y=OccAll*.3, colour= as.factor(Year)),
           sides='t', alpha=.8) +
  geom_rug(data=subset(OccTable_daily, ShoreDist=='10' & OccAll==0), 
           aes(x=DummyDate, y=OccAll, colour= as.factor(Year)),
           sides='b', alpha=.5) +
  xlab('Julien Day') +
  ylab('Detection Probability') +
  ggtitle('Midshore (e.g. 10km)')

ggplot(data=subset(newdat, ShoreDist=='15'))+
  theme_bw() +
  facet_wrap(~GroupId, nrow=2) +
  scale_color_brewer(type = "div", palette = "Set1", direction = 1, name='Year') +  
  geom_line(aes(DummyDate, inv.logit(fit), colour=as.factor(Year)), size=1) +
  geom_ribbon(aes(x=DummyDate, ymin=inv.logit(lwr), ymax=inv.logit(upr)),
              alpha=.6,linetype= 'blank') +
  geom_rug(data=subset(OccTable_daily, ShoreDist=='15' & OccAll==1), 
           aes(x=DummyDate, y=OccAll*.3, colour= as.factor(Year)),
           sides='t', alpha=.8) +
  geom_rug(data=subset(OccTable_daily, ShoreDist=='15' & OccAll==0), 
           aes(x=DummyDate, y=OccAll, colour= as.factor(Year)),
           sides='b', alpha=.5) +
  xlab('Julien Day') +
  ylab('Detection Probability') +
  ggtitle('Offshore (e.g. 15km)')









