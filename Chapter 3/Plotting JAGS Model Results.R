
# Load and plot the Jags output for the multiple Models #

rm(list=ls())
library(ggplot2)
library(boot) # for inverse logit
library(ggjoy) # density plots
library(grid) # arranging plots
library(viridis) # sexy colors
library(rethinking) #HDPI 

setwd("W:/KJP PHD/3-Detection Function/R Code")

# Read in Occupancy Table

OccTable= read.csv('W:/KJP PHD/4-Bayesian Habitat Use/R Code/OccupancyTable_ThreePdets.csv')

# Pull out time with SM coverage
OccTable_SM=subset(OccTable, SMCoverage==1)
OccTable_SM$UnitLoc=(droplevels(OccTable_SM)$UnitLoc)

# Pull out time with SM coverage
OccTable_NL=subset(OccTable, !is.na(MedianNoiseLevel))
OccTable_NL$UnitLoc=(droplevels(OccTable_NL)$UnitLoc)


########################
# Load the model files #
########################
# Tau values 1500, 5000, 500
# SNR threshold values 1, 10,15

# No Pdet consideration at all #
M1.SMonly=read.csv("Results_M1_SMonly_samps.csv" )


# TL and NL values area monitored defined by M1 constraints
M2.SMonly=read.csv("Results_M2_SMonly_samps_b.csv")
M3.SMonly=read.csv("Results_M3_SMonly_samps_b.csv")
M4.SMonly=read.csv("Results_M4_SMonly_samps_b.csv")

# TL values area monitored defined by M1 constraints
M2.TLonly=read.csv("Results_M5_TLonly_samps.csv")
M3.TLonly=read.csv("Results_M6_TLonly_samps.csv")
M4.TLonly=read.csv("Results_M7_TLonly_samps.csv")

# Analysis with all noise levels
M1.NLAll=read.csv("Results_M2_NLAll_samps.csv")
M2.NLAll=read.csv("Results_M2_NLAll_samps.csv")
M3.NLAll=read.csv("Results_M3_NLAll_samps.csv")
M4.NLAll=read.csv("Results_M4_NLAll_samps.csv")


#################################################################################

# Function for creating ordered occupancy values for plotting #
Ordered_Occ<-function(M1.SMonly){
  
  tempdf=(M1.SMonly[,2:ncol(M1.SMonly)])+(M1.SMonly$intercept)
  newdf=data.frame(UnitLoc=(rep(colnames(tempdf)[1], nrow(tempdf))), 
                   alpha=tempdf[,1])
  
  
  # Reshape the Datafrme So it can be plotted
  for (ii in 2:ncol(tempdf)-1){
    tempdf1=data.frame(UnitLoc=(rep(colnames(tempdf)[ii], nrow(tempdf))), 
                        alpha=tempdf[,(ii)])
    newdf<-rbind(tempdf1, newdf)
    print(ii)
    
  }
  
  newdf$P.occ=inv.logit(newdf$alpha)
  newdf$Scaled=scale(inv.logit(newdf$P.occ))
  newdf$Scaled_mean=inv.logit(newdf$P.occ)-mean(inv.logit(newdf$P.occ))
  newdf$Model=substr(deparse(substitute(M1.SMonly)),start = 1, stop = 2)

  return(newdf)
  
}



####################################################################
# Model 1 Unprocessed C-POD Data #
####################################################################
# Plotting 
M1.output<-Ordered_Occ(M1.SMonly)
# ggplot(data=M1.output, aes(UnitLoc, Scaled)) + 
#   geom_abline(intercept = 0, slope = 0, colour='red')+
#   geom_boxplot() + coord_flip()+
#   theme_bw()+
#   scale_y_continuous(limits = c(-3, 3)) +
#   #theme(legend.position="none") +
#   labs(y='Relative Hourly Occupancy', x="Deployment Site") +
#   ggtitle('Model 1 Unprocessed C-POD Detections') +
#   theme(plot.title = element_text(hjust = 0.5)) 

##################################################################
# Model 2 Processed CPOD data with tau 1 #
####################################################
M2.output<-Ordered_Occ(M2.SMonly)

# Plotting 
# ggplot(data=M2.output, aes(UnitLoc, Scaled)) + 
#   geom_abline(intercept = 0, slope = 0, colour='red')+
#   geom_boxplot() + coord_flip()+
#   theme_bw()+
#   scale_y_continuous(limits = c(-3, 3)) +
#   #theme(legend.position="none") +
#   labs(y='Relative Hourly Occupancy', x="Deployment Site") +
#   ggtitle('Model 2 SM Coverage Only') +
#   theme(plot.title = element_text(hjust = 0.5)) 

M3.output<-Ordered_Occ(M3.SMonly)

# # Plotting 
# ggplot(data=M3.output, aes(UnitLoc, Scaled)) + 
#   geom_abline(intercept = 0, slope = 0, colour='red')+
#   geom_boxplot() + coord_flip()+
#   theme_bw()+
#   scale_y_continuous(limits = c(-3, 3)) +
#   #theme(legend.position="none") +
#   labs(y='Relative Hourly Occupancy', x="Deployment Site") +
#   ggtitle('Model 3 SM Coverage Only') +
#   theme(plot.title = element_text(hjust = 0.5)) 

M4.output<-Ordered_Occ(M4.SMonly)
# # Plotting 
# ggplot(data=M4.output, aes(UnitLoc, Scaled)) + 
#   geom_abline(intercept = 0, slope = 0, colour='red')+
#   geom_boxplot() + coord_flip()+
#   theme_bw()+
#   scale_y_continuous(limits = c(-3, 3)) +
#   #theme(legend.position="none") +
#   labs(y='Relative Hourly Occupancy', x="Deployment Site") +
#   ggtitle('Model 4 SM Coverage Only') +
#   theme(plot.title = element_text(hjust = 0.5)) 

################################################################################
# Facet grid plot #
################################################################################
M1.output$Model='Null'
M2.output$Model='T1'
M3.output$Model='T2'
M4.output$Model='T3'

SM_samps=rbind(M1.output,M2.output,M3.output, M4.output)
SM_samps$UnitLoc_order = factor(SM_samps$UnitLoc, levels=rev(levels(SM_samps$UnitLoc)))
SM_samps$Model=factor(SM_samps$Model, levels=c('Null', 'T1', 'T2', 'T3'))

SM_samps$isMod1=as.factor((SM_samps$Model=='M1')*1)


# Get Bayesian Credability intervals
mm=data.frame(aggregate(data=SM_samps, alpha~UnitLoc+Model, FUN = median))
mm=cbind(mm, as.numeric(aggregate(data=SM_samps, 
                                  alpha~UnitLoc+Model,
                                  FUN = function(x){HPDI(x, 0.90)})[[3]])[1:40])
mm=cbind(mm, as.numeric(aggregate(data=SM_samps, 
                                  alpha~UnitLoc+Model,
                                  FUN = function(x){HPDI(x, 0.90)})[[3]])[41:80])

colnames(mm)<-c('UnitLoc','Model','Median', 'LCI', 'UCI' )

# make overlapping matrix for each model
OverlapList=list()

for(ii in 1:4){
  datasub=subset(mm, Model==unique(mm$Model)[ii])
  
  OverlapMat=matrix(data = 0, nrow=10, ncol=10)
  for(jj in 1:10){
    
    OverlapMat[jj,
    which(datasub$Median[jj]>datasub$LCI & datasub$Median[jj]<datasub$UCI)]=1
    
  }
  
  OverlapList[[ii]]=OverlapMat
  
  
  
}







# True Occupancy 
 ggplot(data=SM_samps, aes(Model, P.occ)) +
   facet_wrap(~UnitLoc_order, nrow = 2) + 
   geom_violin(aes(fill=Model), draw_quantiles = c(.05, .5, 0.95))+
   #geom_boxplot(aes(fill=Model),  outlier.size = .1) +
   scale_fill_manual(values=c('steelblue3',rep("white", length(unique(SM_samps$Model))-1))) +
   theme_bw()+
   theme(legend.position="none")+
   labs(y='Occupancy Rate', x='') 
 
# Relative Occupancy  
  # mm=SM_samps[! SM_samps$UnitLoc %in% c('Hel_15', 'SpB_10', 'Abr_10', 'StA_10'), ]
  # ggplot(data=mm, aes(Model, Scaled)) +

  ggplot(data=SM_samps, aes(Model, Scaled)) +
     facet_wrap(~UnitLoc_order, nrow = 2) + 
     geom_violin(aes(fill=Model), draw_quantiles = c(.05, .5, 0.95))+
     #geom_boxplot(aes(fill=Model),  outlier.size = .1) +
     #scale_color_brewer(type='div')
     scale_fill_manual(values=c('steelblue3',rep("white", length(unique(SM_samps$Model))-1))) +
     theme_bw()+
     theme(legend.position="none")+
     labs(y='Relative Occupancy', x='') 



   ggplot(data=SM_samps, aes(UnitLoc, P.occ)) +
     facet_wrap(~Model, nrow = 1, scale='free_x') +
     guides(fill=FALSE)+
     #geom_violin(aes(UnitLoc, Scaled,fill=Model)) +
     scale_fill_manual(values=c('steelblue3',rep("white", length(unique(SM_samps$Model))-1))) +
     geom_boxplot(aes(fill=Model)) + 
     coord_flip() +
     theme_bw() +
     theme(plot.title = element_text(hjust = 0.5)) +
     ggtitle('Posterior Occupancy Rates')+
     labs(x='', y='Occupancy Rate ')  
   
     
    

    
##################################################################################
# Models considering Transmission loss differences only
##################################################################################

 M2.TLonly.output<-Ordered_Occ(M2.TLonly)

# Plotting 
# ggplot(data=M2.SMonly.output, aes(UnitLoc, Scaled)) + 
#   geom_abline(intercept = 0, slope = 0, colour='red')+
#   geom_boxplot() + coord_flip()+
#   theme_bw()+
#   scale_y_continuous(limits = c(-3, 3)) +
#   #theme(legend.position="none") +
#   labs(y='Relative Hourly Occupancy', x="Deployment Site") +
#   ggtitle('Model 5 Transmissin Loss Only') +
#   theme(plot.title = element_text(hjust = 0.5)) 

M3.TLonly.output<-Ordered_Occ(M3.TLonly)

# # Plotting 
# ggplot(data=M3.SMonly.output, aes(UnitLoc, Scaled)) + 
#   geom_abline(intercept = 0, slope = 0, colour='red')+
#   geom_boxplot() + coord_flip()+
#   theme_bw()+
#   scale_y_continuous(limits = c(-3, 3)) +
#   #theme(legend.position="none") +
#   labs(y='Relative Hourly Occupancy', x="Deployment Site") +
#   ggtitle('Model 6 Transmissin Loss Only') +
#   theme(plot.title = element_text(hjust = 0.5)) 

M4.TLonly.output<-Ordered_Occ(M4.TLonly)

# # Plotting 
# ggplot(data=M7.output, aes(UnitLoc, Scaled)) + 
#   geom_abline(intercept = 0, slope = 0, colour='red')+
#   geom_boxplot() + coord_flip()+
#   theme_bw()+
#   scale_y_continuous(limits = c(-3, 3)) +
#   #theme(legend.position="none") +
#   labs(y='Relative Hourly Occupancy', x="Deployment Site") +
#   ggtitle('Model 7 Transmissin Loss Only') +
#   theme(plot.title = element_text(hjust = 0.5)) 

########################################################################
# Facit plots TL only #
########################################################################

SM_samps_TL=rbind(M1.output, M2.TLonly.output, M3.TLonly.output, M4.TLonly.output)
SM_samps_TL$UnitLoc_order = factor(SM_samps_TL$UnitLoc, levels=rev(levels(SM_samps_TL$UnitLoc)))
SM_samps_TL$isMod1=as.factor((SM_samps_TL$Model=='M1')*1)

#nn=SM_samps_TL[! SM_samps_TL$UnitLoc %in% c('Hel_15', 'SpB_10', 'Abr_10', 'StA_10'), ]


# True Occupancy 
ggplot(data=SM_samps_TL, aes(Model, P.occ)) +
  facet_wrap(~UnitLoc_order, nrow = 2) + 
  #geom_abline(intercept = vlindata, slope = 0, colour='red') +
  #geom_violin(aes(fill=Model))+
  geom_boxplot(aes(fill=Model),  outlier.size = .1) +
  #scale_color_brewer(type='div')
  scale_fill_manual(values=c('steelblue3',rep("white", length(unique(SM_samps$Model))-1))) +
  theme_bw()+
  theme(legend.position="none")+
  labs(y='Occupancy Rate', x='') 

# Relative Occupancy   
ggplot(data=SM_samps_TL, aes(Model, Scaled)) +
  facet_wrap(~UnitLoc_order, nrow = 2) + 
  #geom_abline(intercept = vlindata, slope = 0, colour='red') +
  #geom_violin(aes(fill=Model))+
  geom_boxplot(aes(fill=Model),  outlier.size = .1) +
  #scale_color_brewer(type='div')
  scale_fill_manual(values=c('steelblue3',rep("white", length(unique(SM_samps_TL$Model))-1))) +
  theme_bw()+
  theme(legend.position="none")+
  labs(y='Relative Occupancy', x='') 




ggplot(data=SM_samps_TL, aes(UnitLoc, P.occ)) +
  facet_wrap(~Model, nrow = 1, scale='free_x') +
  guides(fill=FALSE)+
  #geom_violin(aes(UnitLoc, Scaled,fill=Model)) +
  scale_fill_manual(values=c('steelblue3',rep("white", length(unique(SM_samps$Model))-1))) +
  geom_boxplot(aes(fill=Model)) + 
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Run 2- C-POD Detections and Transmission Loss Only')+
  labs(x='', y='Posterior Occupancy Rate Distribution')  


#########################################################################
# Look at all groups with Noise Levels #
#########################################################################

M1NLAll.output<-Ordered_Occ(M1.NLAll)
M2NLAll.output<-Ordered_Occ(M2.NLAll)
M3NLAll.output<-Ordered_Occ(M3.NLAll)
M4NLAll.output<-Ordered_Occ(M4.NLAll)

Samps_NLAll=rbind(M1NLAll.output, M2NLAll.output, M3NLAll.output, M4NLAll.output)
Samps_NLAll$UnitLoc_order = factor(Samps_NLAll$UnitLoc, levels=rev(levels(Samps_NLAll$UnitLoc)))
Samps_NLAll$isMod1=as.factor((Samps_NLAll$Model=='M1')*1)

Samps_NLAll_sub=subset(Samps_NLAll, UnitLoc != 'Cro_05')

# True Occupancy 
ggplot(data=Samps_NLAll_sub, aes(Model, P.occ)) +
  facet_wrap(~UnitLoc_order, nrow = 3) + 
  #geom_abline(intercept = vlindata, slope = 0, colour='red') +
  #geom_violin(aes(fill=Model))+
  geom_boxplot(aes(fill=Model),  outlier.size = .1) +
  #scale_color_brewer(type='div')
  scale_fill_manual(values=c('steelblue3',rep("white", length(unique(Samps_NLAll_sub$Model))-1))) +
  theme_bw()+
  theme(legend.position="none")+
  labs(y='Occupancy Rate', x='') 

# Relative Occupancy   
ggplot(data=Samps_NLAll_sub, aes(Model, Scaled)) +
  facet_wrap(~UnitLoc_order, nrow = 3) + 
  #geom_abline(intercept = vlindata, slope = 0, colour='red') +
  #geom_violin(aes(fill=Model))+
  geom_boxplot(aes(fill=Model),  outlier.size = .1) +
  #scale_color_brewer(type='div')
  scale_fill_manual(values=c('steelblue3',rep("white", length(unique(Samps_NLAll$Model))-1))) +
  theme_bw()+
  theme(legend.position="none")+
  labs(y='Relative Occupancy', x='') 





ggplot(data=Samps_NLAll_sub, aes(UnitLoc_order, P.occ)) +
  facet_wrap(~Model, nrow = 1, scale='free_x') +
  guides(fill=FALSE)+
  #geom_violin(aes(UnitLoc, Scaled,fill=Model)) +
  scale_fill_manual(values=c('steelblue3',rep("white", length(unique(SM_samps$Model))-1))) +
  geom_boxplot(aes(fill=Model)) + 
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Run 3- C-POD Detections and NL from nearest SM2M')+
  labs(x='', y='Posterior Occupancy Rate Distribution')  





ggplot(data=Samps_NLAll_sub, aes(UnitLoc, Scaled)) +
  facet_wrap(~Model, nrow = 1) +
  geom_boxplot(aes(fill=Model)) + 
  coord_flip() +
  theme_bw() +
  labs(y='', x='Relative Occupancy') +
  scale_fill_manual(values=c('steelblue3',rep("white", length(unique(Samps_NLAll_sub$Model))-1))) 


#################################################################################
# Plot the Area Monitored #
################################################################################

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

OccTable_SM=subset(OccTable, SMCoverage==1)
OccTable_SM$UnitLoc=(droplevels(OccTable_SM)$UnitLoc)


newdf=data.frame(UnitLoc=OccTable_SM$UnitLoc, Area=OccTable_SM$MedianArea1, Model='T1')
newdf=rbind(newdf, data.frame(UnitLoc=OccTable_SM$UnitLoc, Area=OccTable_SM$MedianArea2, Model='T2'))
newdf=rbind(newdf, data.frame(UnitLoc=OccTable_SM$UnitLoc, Area=OccTable_SM$MedianArea3, Model='T3'))

newdf$UnitLoc_order = factor(newdf$UnitLoc, levels=rev(levels(newdf$UnitLoc)))

ggplot(data=newdf, aes(UnitLoc_order, Area, Model)) +
  facet_wrap(~Model, nrow = 1, scale='free') +
  guides(fill=FALSE) +
  #geom_boxplot()+
  geom_violin(fill='black')+
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Area Monitored Under 3 Threshold Scenarios')+
  labs(x='', y=expression(Area~Monitored~" "~(km^{2})))


# This is ok but not super nice

ggplot(data=newdf, aes(Area, UnitLoc_order, Model)) +
  facet_wrap(~Model, nrow = 1, scale='free_x') +
  geom_joy(rel_min_height = 0.001, scale = 2) +
  # geom_joy(aes(fill=UnitLoc),rel_min_height = 0.005, scale = 2) + 
  # scale_fill_viridis(discrete = T, direction = -1,
  #                    begin = .1, end = .9) +
  theme_bw() + 
  theme(legend.position="none") +
  ggtitle('Area Monitored Under 3 Threshold Scenarios')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y='', x=expression(Area~Monitored~" "~(km^{2})))


# Make Individual Plots

p1=ggplot(OccTable_SM, aes(x = MedianArea1, y=UnitLoc, fill=UnitLoc)) + 
  geom_joy(aes(fill=UnitLoc),rel_min_height = 0.005, scale = 3) + 
  scale_fill_viridis(discrete = T, direction = -1, 
                     begin = .1, end = .9) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  xlim(c(15,41))+
  labs(y='', x=expression(Area~Monitored~""~(km^{2}))) + 
  ylab("") +
  ggtitle('M1')


p2=ggplot(OccTable_SM, aes(x = MedianArea2, y=UnitLoc)) + 
  geom_joy(aes(fill=UnitLoc),rel_min_height = 0.005, scale = 3) + 
  scale_fill_viridis(discrete = T, direction = -1, 
                     begin = .1, end = .9) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  xlim(c(25,47))+
  labs(y='', x=expression(Area~Monitored~""~(km^{2}))) + 
  ylab("") +
  ggtitle('M2')


p3=ggplot(OccTable_SM, aes(x = MedianArea3, y=UnitLoc)) + 
  geom_joy(aes(fill=UnitLoc),rel_min_height = 0.005, scale = 3) + 
  scale_fill_viridis(discrete = T, direction = -1, 
                     begin = .1, end = .9) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  xlim(c(12,30))+
  labs(y='', x=expression(Area~Monitored~""~(km^{2}))) + 
  ylab("") +
  ggtitle('M3')


grid.arrange(p1, p2, p3, nrow=1)




