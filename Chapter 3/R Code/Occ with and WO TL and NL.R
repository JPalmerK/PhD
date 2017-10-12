
#########################################################################
# Detection Probability Sensitivity Analysis #
#########################################################################

# 1) Investigate relative occupancy rates
#   a) Without conisering relative area monitored
#   b) With area monitored as a function of time and place
#   c) With transmission loss only 


#############
# Setup     #
#############
# This code investigates the various models that look at the different occupancy distributions
# incorporated by the data 
rm(list=ls())
library(rstan)
library(coda)
library(solidearthtide) #For datenum, convert date to matlabdate
library(boot) # for inv.logit
library(rethinking)
set.seed(20151204)


setwd("W:/KJP PHD/3-Detection Function/R Code")
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



meta=read.csv('W:/KJP PHD/CPOD Processing/2013 to 2016 SM deployments.csv')
meta$UnitLoc=factor(meta$UnitLoc, levels=level_names)


################################################################################
# General Data Prep #
################################################################################



############################################################################
# Trim the data such that only C-PODs with an adjacen SM are incorporated  #
############################################################################

meta_sm=subset(meta, Year != 2015)
meta_sm$MatlabEnddate=Datenum(as.Date(meta_sm$Usable.from.date, format= '%d/%m/%Y'))

# Pull out time periods with SM coverage

OccTable_SM=subset(OccTable, SMCoverage==1)
OccTable_SM$UnitLoc=(droplevels(OccTable_SM)$UnitLoc)

OccTable_SM$CV=NA
OccTable_SM$CV[OccTable_SM$UnitLoc=="Lat_05"]=0.02
OccTable_SM$CV[OccTable_SM$UnitLoc=="Hel_15"]=0.01
OccTable_SM$CV[OccTable_SM$UnitLoc=="Cro_15"]=0.03
OccTable_SM$CV[OccTable_SM$UnitLoc=="SpB_10"]=0.04
OccTable_SM$CV[OccTable_SM$UnitLoc=="Fra_05"]=0.15
OccTable_SM$CV[OccTable_SM$UnitLoc=="Cru_05"]=0.01
OccTable_SM$CV[OccTable_SM$UnitLoc=="Sto_05"]=0.06
OccTable_SM$CV[OccTable_SM$UnitLoc=="Abr_10"]=0.02
OccTable_SM$CV[OccTable_SM$UnitLoc=="StA_10"]=0.0001
OccTable_SM$CV[OccTable_SM$UnitLoc=="Stb_05"]=0.16

OccTable_SM$MedianArea1_sd=OccTable_SM$CV*OccTable_SM$MedianArea1
OccTable_SM$MedianArea2_sd=OccTable_SM$CV*OccTable_SM$MedianArea2
OccTable_SM$MedianArea3_sd=OccTable_SM$CV*OccTable_SM$MedianArea3




################################################################################



mod1=glm(OccAll~UnitLoc, family = 'binomial', data=OccTable_SM)
mod2=glm(OccAll~UnitLoc, family = 'binomial', data=OccTable_SM, 
         offset = inv.logit(scale(MedianPdet1)))
mod3=glm(OccAll~UnitLoc, family = 'binomial', data=OccTable_SM, 
         offset =  inv.logit(scale(MedianPdet2)))
mod4=glm(OccAll~UnitLoc, family = 'binomial', data=OccTable_SM, 
         offset =  inv.logit(scale(MedianPdet3)))

mod1.coef=summary(mod1)
mod1.coef=mod1.coef$coefficients
mod2.coef=summary(mod2)
mod2.coef=mod2.coef$coefficients
mod3.coef=summary(mod3)
mod3.coef=mod3.coef$coefficients
mod4.coef=summary(mod4)
mod4.coef=mod4.coef$coefficients

NewDf=data.frame(aggregate(MedianPdet1~UnitLoc, data=OccTable_SM, FUN=median))


NewDf=cbind(NewDf, aggregate(MedianPdet2~UnitLoc, data=OccTable_SM, FUN=median)[2])
NewDf=cbind(NewDf, aggregate(MedianPdet3~UnitLoc, data=OccTable_SM, FUN=median)[2])


Predvalues=OccTable_SM[, c('UnitLoc', 'MedianPdet1', 'MedianPdet2', 'MedianPdet3', 'OccAll')]
Predvalues$M1.pred=as.numeric(inv.logit(predict(object = mod1, newdata = Predvalues)))
Predvalues$M2.pred=as.numeric(inv.logit(predict(object = mod2, newdata = Predvalues)))
Predvalues$M3.pred=as.numeric(inv.logit(predict(object = mod3, newdata = Predvalues)))
Predvalues$M4.pred=as.numeric(inv.logit(predict(object = mod4, newdata = Predvalues)))

Predvalues$M1.occ=rbinom(nrow(Predvalues), 1, Predvalues$M1.pred)
Predvalues$M2.occ=rbinom(nrow(Predvalues), 1, Predvalues$M2.pred)
Predvalues$M3.occ=rbinom(nrow(Predvalues), 1, Predvalues$M3.pred)
Predvalues$M4.occ=rbinom(nrow(Predvalues), 1, Predvalues$M4.pred)

as.numeric(inv.logit(predict(object = mod1, newdata = Predvalues,se.fit = T)))

################################################################################
# Aggrage Table #
###############################################################################

AggTable=data.frame(UnitLoc= as.matrix(aggregate(Predvalues$OccAll~Predvalues$UnitLoc, FUN=mean))[,1],
                    Obs=  as.numeric(unlist((aggregate(Predvalues$OccAll~Predvalues$UnitLoc, FUN=mean))[,2])))

AggTable$M1= as.numeric(unlist((aggregate(Predvalues$M1.occ~Predvalues$UnitLoc, FUN=mean))[,2]))
AggTable$M2= as.numeric(unlist((aggregate(Predvalues$M2.occ~Predvalues$UnitLoc, FUN=mean))[,2]))
AggTable$M3= as.numeric(unlist((aggregate(Predvalues$M3.occ~Predvalues$UnitLoc, FUN=mean))[,2]))
AggTable$M4= as.numeric(unlist((aggregate(Predvalues$M4.occ~Predvalues$UnitLoc, FUN=mean))[,2]))
AggTable$Nhrs= as.numeric(unlist((aggregate(Predvalues$M4.occ~Predvalues$UnitLoc, FUN=length))[,2]))



###############################################################################

ScaledValues=data.frame(UnitLoc= AggTable$UnitLoc)
ScaledValues=cbind(ScaledValues, apply(AggTable[,2:5], FUN=scale, MARGIN = 2))
###############################################################################

# Days occupied
PredictedDays=data.frame(UnitLoc= NewDf$UnitLoc)
PredictedDays=cbind(PredictedDays,aggregate(OccTable_SM$OccAll~OccTable_SM$UnitLoc, FUN=length)[2])
PredictedDays=cbind(PredictedDays, aggregate(OccTable_SM$OccAll~OccTable_SM$UnitLoc, FUN=sum)[2])

colnames(PredictedDays)<-c('UnitLoc', 'Nhrs', 'NDet')

PredictedDays$M1=round(PredictedDays$Nhrs*Predvalues$M1.pred, 2)
PredictedDays$M2=round(PredictedDays$Nhrs*Predvalues$M2.pred, 2)
PredictedDays$M3=round(PredictedDays$Nhrs*Predvalues$M3.pred, 2)
PredictedDays$M4=round(PredictedDays$Nhrs*Predvalues$M4.pred, 2)
##########################################################################################################################

##############################
# Data prep for sim 1 #
##############################


library(R2jags)
library(rjags)
library(runjags)

Site=as.numeric(OccTable_SM$UnitLoc)
N=length(Site)
NSites=length(unique(Site))
y=OccTable_SM$OccAll
AreaM1Scaled=OccTable_SM$MedianArea1
sd=OccTable_SM$MedianArea1_sd
###################
#  initial values #
###################
intercept <- list(-20, 20)
alpha.occ <- list(chain1=c(NA, rep(1,NSites-1)), chain2=c(NA, rep(-1,NSites-1)))

########################################
# Model 1 Pretend Nothing is Important #
########################################

Results.M1 <- autorun.jags(model='M1_jags.txt', n.chains=2)
Results.M1.samps=as.data.frame(do.call(rbind, Results.M1[["mcmc"]]))
colnames(Results.M1.samps)[1:NSites]=levels(OccTable_SM$UnitLoc)

posterior_name=paste('Results_M1_SMonly_samps.csv')
write.csv(Results.M1.samps, posterior_name)

##############################
# Data prep for sim 2 #
##############################

MedianRange1_scaled=as.numeric(scale(OccTable_SM$MedianArea1))
MedianArea=OccTable_SM$MedianArea1
sd=OccTable_SM$MedianArea1_sd

Results.M2  <- autorun.jags(model='M2_jags_b.txt', n.chains=2)
Results.M2.samps=as.data.frame(do.call(rbind, Results.M2[["mcmc"]]))
colnames(Results.M2.samps)[1:NSites]=levels(OccTable_SM$UnitLoc)




posterior_name=paste('Results_M2_SMonly_samps_b.csv') # model b has standard deviation accounted
write.csv(Results.M2.samps, posterior_name)

##############################
# Data prep for sim 3 #
##############################

MedianRange1_scaled=as.numeric(scale(OccTable_SM$MedianArea2))
MedianArea=OccTable_SM$MedianArea2
sd=OccTable_SM$MedianArea2_sd



Results.M3 <- autorun.jags(model='M2_jags_b.txt', n.chains=2)
Results.M3.samps=as.data.frame(do.call(rbind, Results.M3[["mcmc"]]))
colnames(Results.M3.samps)[1:NSites]=levels(OccTable_SM$UnitLoc)

posterior_name=paste('Results_M3_SMonly_samps_b.csv')
write.csv(Results.M3.samps, posterior_name)



##############################
# Data prep for sim 4 #
##############################

MedianRange1_scaled=as.numeric(scale(OccTable_SM$MedianArea3))
MedianArea=OccTable_SM$MedianArea3
sd=OccTable_SM$MedianArea3_sd

Results.M4 <- run.jags(model='M2_jags_b.txt', n.chains=2)
Results.M4.samps=as.data.frame(do.call(rbind, Results.M4[["mcmc"]]))
colnames(Results.M4.samps)[1:NSites]=levels(OccTable_SM$UnitLoc)

posterior_name=paste('Results_M4_SMonly_samps_b.csv')
write.csv(Results.M4.samps, posterior_name)



#######################################################################
# TL only #
#######################################################################

MinArea=data.frame(aggregate(MedianArea1~UnitLoc, data=OccTable_SM, FUN=min))
MinArea=cbind(MinArea, aggregate(MedianArea2~UnitLoc, data=OccTable_SM, FUN=min)[2])
MinArea=cbind(MinArea, aggregate(MedianArea3~UnitLoc, data=OccTable_SM, FUN=min)[2])
colnames(MinArea)[2:4]=c('MinAreaM1','MinAreaM2','MinAreaM3')

OccTable_SM=merge(OccTable_SM, MinArea)

##############################
# Data prep for sim 5 TL Only #
##############################

MedianRange1_scaled=as.numeric(scale(OccTable_SM$MinAreaM1))
MedianArea=OccTable_SM$MinAreaM1
sd=OccTable_SM$MedianArea1_sd

Results.M5  <- autorun.jags(model='M2_jags_b.txt', n.chains=2)
Results.M5.samps=as.data.frame(do.call(rbind, Results.M5[["mcmc"]]))
colnames(Results.M5.samps)[1:NSites]=levels(OccTable_SM$UnitLoc)

posterior_name=paste('Results_M5_TLonly_samps_b.csv')
write.csv(Results.M5.samps, posterior_name)




##############################
# Data prep for sim 6 TL Only #
##############################

MedianRange1_scaled=as.numeric(scale(OccTable_SM$MinAreaM2))
MedianArea=OccTable_SM$MinAreaM2
sd=OccTable_SM$MedianArea2_sd

Results.M6  <- autorun.jags(model='M2_jags_b.txt', n.chains=2)
Results.M6.samps=as.data.frame(do.call(rbind, Results.M6[["mcmc"]]))
colnames(Results.M6.samps)[1:NSites]=levels(OccTable_SM$UnitLoc)


posterior_name=paste('Results_M6_TLonly_samps_b.csv')
write.csv(Results.M6.samps, posterior_name)


###############################
# Data prep for sim 7 TL Only #
###############################

MedianRange1_scaled=as.numeric(scale(OccTable_SM$MinAreaM3))
MedianArea=OccTable_SM$MinAreaM3
sd=OccTable_SM$MedianArea3_sd

Results.M7  <- autorun.jags(model='M2_jags_b.txt', n.chains=2)
Results.M7.samps=as.data.frame(do.call(rbind, Results.M7[["mcmc"]]))
colnames(Results.M7.samps)[1:NSites]=levels(OccTable_SM$UnitLoc)


posterior_name=paste('Results_M7_TLonly_samps_b.csv')
write.csv(Results.M7.samps, posterior_name)

####################################################################################################
# Expand the analysis to include all data with noise levels #
####################################################################################################

OccTable_NL=subset(OccTable, !is.na(MedianNoiseLevel))
OccTable_NL$UnitLoc=(droplevels(OccTable_NL)$UnitLoc)

Site=as.numeric(OccTable_NL$UnitLoc)
N=length(Site)
NSites=length(unique(Site))
y=OccTable_NL$OccAll
AreaM1Scaled=scale(OccTable_NL$MedianArea1)

###################
#  initial values #
###################
intercept <- list(-20, 20)
alpha.occ <- list(chain1=c(NA, rep(3,NSites-1)), chain2=c(NA, rep(-3,NSites-1)))



##############################
# Data prep for NL sim 1 #
##############################

Results.M1 <- autorun.jags(model='M1_jags.txt', n.chains=2)
Results.M1.NL=as.data.frame(do.call(rbind, Results.M1[["mcmc"]]))
colnames(Results.M1.NL)[1:NSites]=levels(OccTable_NL$UnitLoc)

posterior_name=paste('Results_M1_NLAll_samps.csv')
write.csv(Results.M1.NL, posterior_name)






##############################
# Data prep for NL sim 2 #
##############################

MedianRange1_scaled=as.numeric(scale(OccTable_NL$MedianArea1))
MedianArea=OccTable_NL$MedianArea1
sd=OccTable_NL$MedianArea1_sd


Results.M2.NL  <- autorun.jags(model='M2_jags_b.txt', n.chains=2)
Results_M2_NLsamps=as.data.frame(do.call(rbind, Results.M2.NL[["mcmc"]]))
colnames(Results_M2_NLsamps)[1:NSites]=levels(OccTable_NL$UnitLoc)

posterior_name=paste('Results_M2_NLAll_samps_b.csv')
write.csv(Results_M2_NLsamps, posterior_name)

##############################
# Data prep for NL sim 3 #
##############################

MedianRange1_scaled=as.numeric(scale(OccTable_NL$MedianArea2))
MedianArea=OccTable_NL$MedianArea2
sd=OccTable_NL$MedianArea2_sd


Results.M3.NL  <- autorun.jags(model='M2_jags_b.txt', n.chains=2)
Results_M3_NLsamps=as.data.frame(do.call(rbind, Results.M3.NL[["mcmc"]]))
colnames(Results_M3_NLsamps)[1:NSites]=levels(OccTable_NL$UnitLoc)

posterior_name=paste('Results_M3_NLAll_samps_b.csv')
write.csv(Results_M3_NLsamps, posterior_name)

##############################
# Data prep for NL sim 4 #
##############################

MedianRange1_scaled=as.numeric(scale(OccTable_NL$MedianArea3))
MedianArea=OccTable_NL$MedianArea3
sd=OccTable_NL$MedianArea3_sd


Results.M4.NL  <- autorun.jags(model='M2_jags_b.txt', n.chains=2)
Results_M4_NLsamps=as.data.frame(do.call(rbind, Results.M4.NL[["mcmc"]]))
colnames(Results_M4_NLsamps)[1:NSites]=levels(OccTable_NL$UnitLoc)

posterior_name=paste('Results_M4_NLAll_samps_b.csv')
write.csv(Results_M4_NLsamps, posterior_name)





  
  