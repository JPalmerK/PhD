
# This file reads in processed CPOD data then does some stuff


rm(list=ls())
library(ggplot2)
library(RColorBrewer)
#setwd("W:/KJP PHD/4-Bayesian Habitat Use/R Code")

# Deployment information
meta=read.csv("/home/kpalmer/PHD/PhD/Deployment Information/CPODs for Kaitlin.csv")

meta$Usable.from.date=as.Date(meta$Usable.from.date, "%d/%m/%Y", tz = "GMT")
meta$Usable.until.date=as.Date(meta$Usable.until.date,"%d/%m/%Y", tz = "GMT")
meta$Year=2014
meta$Year[meta$Deployment.number<31]=2013
meta$Year[meta$Deployment.number>=61]=2015


Trains2013=read.csv("/home/kpalmer/PHD/PhD/Chapter 4/DolCPODFiles/2013_ClickTrainAnalysis.csv")
Trains2014=read.csv("/home/kpalmer/PHD/PhD/Chapter 4/DolCPODFiles/2014_ClickTrainAnalysis.csv")
Trains2015=read.csv("/home/kpalmer/PHD/PhD/Chapter 4/DolCPODFiles/2015_ClickTrainAnalysis.csv")

Trains2013$Year=2013
Trains2014$Year=2014
Trains2015$Year=2015
Trains2014$EncounterID=Trains2014$EncounterID+max(Trains2013$EncounterID)
Trains2015$EncounterID=Trains2015$EncounterID+max(Trains2014$EncounterID)

Trains2015$VerfiedSpp=NA



# Fix level names
level_names=c( "Lat_05", "Lat_10", "Lat_15",
               "Hel_05", "Hel_10", "Hel_15",
               "Cro_05", "Cro_10" ,"Cro_15",
               "SpB_05", "SpB_10", "SpB_15",
               "Fra_05", "Fra_10", "Fra_15",      
               "Cru_05", "Cru_10", "Cru_15",
               "Sto_05", "Sto_10", "Sto_15",
               "Abr_05", "Abr_10", "Abr_15",
               "StA_05", "StA_10", "StA_15",
               "Stb_05", "Stb_10", "Stb_15")

# Set proper level names
# meta$UnitLoc<-factor(meta$UnitLoc, levels=rev(level_names))
# meta=meta[order(meta$UnitLoc), ]



# 
# Trains2015$UnitLoc=factor(Trains2015$UnitLoc, levels=rev(level_names))
# Trains2014$UnitLoc=factor(Trains2014$UnitLoc, levels=rev(level_names))
# Trains2013$UnitLoc=factor(Trains2013$UnitLoc, levels=rev(level_names))

####################################
# For plotting horizontally
meta$UnitLoc<-factor(meta$UnitLoc, levels=(level_names))
meta=meta[order(meta$UnitLoc), ]


Trains2015$UnitLoc=factor(Trains2015$UnitLoc, levels=(level_names))
Trains2014$UnitLoc=factor(Trains2014$UnitLoc, levels=(level_names))
Trains2013$UnitLoc=factor(Trains2013$UnitLoc, levels=(level_names))

##################################

Trains2015$TrainDateNoTime= as.Date(Trains2015$TrainDateNoTime, "%d/%m/%Y")
Trains2015=Trains2015[Trains2015$TrainDateNoTime> as.Date("01/01/2012", format="%d/%m/%Y"),]
Trains2014$TrainDateNoTime= as.Date(Trains2014$TrainDateNoTime, "%d/%m/%Y")
Trains2013$TrainDateNoTime= as.Date(Trains2013$TrainDateNoTime, "%d/%m/%Y")


Trains2015$TrainDay = NULL

Trains=rbind(Trains2013, Trains2014, Trains2015)
Trains=Trains[sample(nrow(Trains), replace = FALSE),]
###########################################################################
# Add Tidal Phase, Time of Day and a Dummy for Day and Month#
###########################################################################
TideData=read.csv('TidalHeight.csv')
TideData$Phase=as.factor(TideData$Phase)
levels(TideData$Phase)=c('Low', 'Flood', 'High', 'Ebb')

temp=TideData[,c('Z','Date','Phase', 'UnitLoc', 'Time')]
temp$Date=as.Date(temp$Date, "%Y-%m-%d", tz = "GMT")

Trains=merge(Trains, temp, by.x= c('TrainDateNoTime', 'UnitLoc', 'TrainHour'), by.y= c('Date', 'UnitLoc', 'Time'))
rm(temp)

# Add decimal time
Trains$DecimalHour=(Trains$MatlabDecDate-floor(Trains$MatlabDecDate))*24

# morning noon evening
Trains$ToD=0
Trains$ToD[Trains$TrainHour>=4& Trains$TrainHour<10]='Dawn'
Trains$ToD[Trains$TrainHour>=10& Trains$TrainHour<16]='Day'
Trains$ToD[Trains$TrainHour>=16& Trains$TrainHour<22]='Dusk'
Trains$ToD[which(Trains$ToD==0)]='Night'

# Trains Dummy Date (for combining across years)
Trains$DummyDate=as.Date(paste(as.character(format(Trains$TrainDateNoTime, "%d-%m")),
                       '-1900', sep = ""), "%d-%m-%Y") 

############################################################################
# Add day of Year
Trains$JDOY=as.numeric(strftime(Trains$TrainDateNoTime, format = "%j"))
Trains$DOY=strftime(Trains$TrainDateNoTime, format = "%d/%m")

# Add Encounter Duration #
Trains$EncStart=0
Trains$EncEnd=0
Trains$EncDur=0
Trains$EncN=0
Trains$EncMinICI=0
Trains$EncLrICIQuartile=0
Trains$EncPropLowICI=0


# Add encouter quartile
for (ii in 1: max(Trains$EncounterID)){
  idx=which(Trains$EncounterID==ii)
  Trains$EncStart[idx]=min(Trains$MatlabDecDate[idx])
  Trains$EncEnd[idx]=max(Trains$MatlabDecDate[idx])
  Trains$EncDur[idx]=(max(Trains$MatlabDecDate[idx])-min(Trains$MatlabDecDate[idx]))*86400
  Trains$EncN[idx]=length(idx)
  Trains$EncMinICI[idx]=min(Trains$MedICI[idx])
  Trains$EncLrICIQuartile[idx]=as.numeric(quantile(Trains$MedICI[idx], 0.25))
  Trains$EncPropLowICI[idx]=sum(log10(Trains$MedICI[idx])<=-2.5)/length(idx)
  }

###################################################################################
# Do Some Exploring
###################################################################################









###############################################################################################
# Plot the usable data for 2013 Along with Detections#
###############################################################################################
Trains2015=subset(Trains, Year==2015)
Trains2014=subset(Trains, Year==2014)
Trains2013=subset(Trains, Year==2013)

# Trains2015=subset(Trains, Year==2015 & EncounterSpp=='COD/BND' )
# Trains2014=subset(Trains, Year==2014 & EncounterSpp=='COD/BND')
# Trains2013=subset(Trains, Year==2013 & EncounterSpp=='COD/BND')


meta2013=subset(meta, Year==2013)
meta2014=subset(meta, Year==2014)
meta2015=subset(meta, Year==2015)
meta2015a=subset(meta2015, Deployment.number<91)
meta2015b=subset(meta2015, Deployment.number>=91)

# ggplot()+geom_rect(data=meta2013, aes(x=UnitLoc,
#                                 xmin= seq(1:30)+.45,
#                                 xmax= seq(1:30)-.45,
#                          ymin=Usable.from.date,
#                          ymax=Usable.until.date),
#                     alpha=0.4) +
#   #scale_x_date(limits = c(as.Date("20/04/2013", format="%d/%m/%Y"),
#                           as.Date("30/11/2013", format="%d/%m/%Y"))) +
#   #scale_x_discrete(limits = rev(levels(meta2013$UnitLoc))) +
#   coord_flip()+
#   geom_jitter(data=Trains2013, aes(y=UnitLoc, x=TrainDateNoTime,
#                                    color=EncounterSpp)) 
#   
#  
#   theme_bw()+
#   theme(legend.position="none", text = element_text(size=15), 
#         axis.text.x  = element_text(angle=90, vjust=.75, hjust=.75,size=15)) +
#   theme(axis.title = element_text(face="bold", size=16)) +
#  
#  

# # Horizontal Plot
# meta2013$UnitLoc <- factor(meta2013$UnitLoc, levels = rev(level_names))
# meta2013$Order[16:18]=12:14
# meta2013$OrderSwap=rev(meta2013$Order)
# Trains2013$UnitLoc=factor(Trains2013$UnitLoc, levels=rev(level_names))
# 
# 
# levels(Trains2013$UnitLoc)=rev(level_names)
# flevels <- levels(meta2013$UnitLoc)
# flevels <- rev(flevels)

meta2 = meta2013
trains2 = Trains2013
meta2$UnitLoc = with(meta2, factor(UnitLoc, levels = rev(levels(UnitLoc))))
trains2$Unitloc =  with(trains2, factor(UnitLoc, levels = rev(levels(UnitLoc))))


trains3 = trains2
trains3$UnitLoc <- factor(trains3$UnitLoc,
                          levels = rev(as.character(levels(trains3$UnitLoc))))



# Verticle #

ggplot()+
  theme_bw()+
  geom_rect(data=meta2, aes(y= UnitLoc,
                            ymin= seq(1:30)+.45,
                            ymax= seq(1:30)-.45,
                            xmin=rev(Usable.from.date),
                            xmax=rev(Usable.until.date),
                            alpha=0.4)) +
  guides(alpha = FALSE, size = FALSE) +
  ggtitle('2013')+
  scale_x_date(limits = c(as.Date("20/04/2013", format="%d/%m/%Y"),
                        as.Date("30/11/2013", format="%d/%m/%Y"))) +
  geom_jitter(data=trains2, 
              width = 0.01,
              aes(y=UnitLoc, x=TrainDateNoTime, color=EncounterSpp)) +  
  scale_color_discrete(
                      name="Encounter\nClassification",
                      breaks=c("COD/BND", "UNK", "WBD/RSD"),
                      labels=c("Broadband", "Unknown", "Frequency\nBanded")) +
  theme(#panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab('Month of Year') 


  ggsave("/media/kpalmer/My Passport/Figs For Aquatic Conservation/2013 All Occupancy .tiff", 
         units="in", width=18, height=8, dpi=800, compression = 'lzw')


  
  
  ggplot()+
    theme_bw() +
    geom_rect(data=meta2, aes(y= UnitLoc,
                              ymin= seq(1:30)+.45,
                              ymax= seq(1:30)-.45,
                              xmin=rev(Usable.from.date),
                              xmax=rev(Usable.until.date)),
              alpha=0.4) +
    scale_x_date(limits = c(as.Date("20/04/2013", format="%d/%m/%Y"),
                            as.Date("30/11/2013", format="%d/%m/%Y"))) +
    ggtitle('2013')+
    #guides(alpha = FALSE, size = FALSE)   + 
    geom_rect(data=trains3, aes(xmin = TrainDateNoTime,
                                xmax = TrainDateNoTime + 1, 
                                ymin = as.numeric(trains3$UnitLoc)-.45,
                                ymax = as.numeric(trains3$UnitLoc)+.45, 
                                fill=EncounterSpp), alpha=.7) +
    scale_fill_manual(values = c("#2B80C2","#a3a3c2", "#FFA42A"),
                      labels=c("Broadband","Unknown","Frequency Banded"),
                      name="Encounter Category") +
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      plot.title = element_text(hjust = 0.5)) 
  
  ggsave("/media/kpalmer/My Passport/Figs For Aquatic Conservation/2013 All Occupancy rect.tiff",
         units="in", width=18, height=8, dpi=800, compression = 'lzw')
  
  
  

###$#########
# 2014 data
############
meta2 = meta2014
trains2 = Trains2014
trains2 = trains2[which(as.Date(trains2$TrainDateNoTime) > min(as.Date(meta2$Usable.from.date), na.rm = TRUE)),]

trains3 = trains2
trains3$UnitLoc <- factor(trains3$UnitLoc,
                          levels = rev(as.character(levels(trains3$UnitLoc))))


meta2$UnitLoc = with(meta2, factor(UnitLoc, levels = rev(levels(UnitLoc))))
trains2$Unitloc =  with(trains2, factor(UnitLoc, levels = rev(levels(UnitLoc))))


ggplot()+
  theme_bw()+
  geom_rect(data=meta2, aes(y= UnitLoc,
                            ymin= seq(1:30)+.45,
                            ymax= seq(1:30)-.45,
                            xmin=rev(Usable.from.date),
                            xmax=rev(Usable.until.date),
                            alpha=0.4)) +
  guides(alpha = FALSE, size = FALSE)  +
  
  ggtitle('2014')+
  scale_x_date(limits = c(as.Date("20/04/2014", format="%d/%m/%Y"),
                          as.Date("30/11/2014", format="%d/%m/%Y"))) +
  geom_jitter(data=trains2, 
              width = 0.01,
              aes(y=UnitLoc, x=TrainDateNoTime, color=EncounterSpp)) +  
  scale_color_discrete(
    name="Encounter\nClassification",
    breaks=c("COD/BND", "UNK", "WBD/RSD"),
    labels=c("Broadband", "Unknown", "Frequency\nBanded")) +
  
  theme(#panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      plot.title = element_text(hjust = 0.5)) 

  ggsave("/media/kpalmer/My Passport/Figs For Aquatic Conservation/2014 All Occupancy .tiff",
       units="in", width=18, height=8, dpi=800, compression = 'lzw')

  
  
  
  
  ggplot()+
    theme_bw() +
    geom_rect(data=meta2, aes(y= UnitLoc,
                              ymin= seq(1:30)+.45,
                              ymax= seq(1:30)-.45,
                              xmin=rev(Usable.from.date),
                              xmax=rev(Usable.until.date)),
                              alpha=0.4) +
    ggtitle('2014')+
    scale_x_date(limits = c(as.Date("20/04/2014", format="%d/%m/%Y"),
                            as.Date("30/11/2014", format="%d/%m/%Y"))) +
    #guides(alpha = FALSE, size = FALSE)   + 
    geom_rect(data=trains3, aes(xmin = TrainDateNoTime,
                                xmax = TrainDateNoTime + 1, 
                                ymin = as.numeric(trains3$UnitLoc)-.45,
                                ymax = as.numeric(trains3$UnitLoc)+.45, 
                                fill=EncounterSpp), alpha=.7) +
    scale_fill_manual(values = c("#2B80C2","#a3a3c2", "#FFA42A"),
                      labels=c("Broadband","Unknown","Frequency Banded"),
                      name="Encounter Category") +
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      plot.title = element_text(hjust = 0.5)) 
  
  ggsave("/media/kpalmer/My Passport/Figs For Aquatic Conservation/2014 All Occupancy rect.tiff",
         units="in", width=18, height=8, dpi=800, compression = 'lzw')
  
  
  
  
  
  
  
  
  
############
# 2015 data
############

# Verticle
  
meta2 = meta2015a
trains2 = Trains2015
trains2 = trains2[-which(trains2$UnitLoc == 'StA_15' & trains2$TrainMonth>7),]
trains3 = trains2
trains3$UnitLoc <- factor(trains3$UnitLoc, levels = rev(as.character(levels(trains3$UnitLoc))))


meta2$UnitLoc = with(meta2, factor(UnitLoc, levels = rev(levels(UnitLoc))))
trains2$Unitloc =  with(trains2, factor(UnitLoc, levels = rev(levels(UnitLoc))))

meta2b = meta2015b
meta2b$UnitLoc = with(meta2b,
                      factor(UnitLoc, levels = rev(levels(UnitLoc))))
meta2b$Usable.until.date[meta2b$UnitLoc == 'Cru_10'] = meta2b$Usable.until.date[1]-15


ggplot()+
  geom_rect(data=meta2, aes(
              y=UnitLoc,
              ymin= seq(1:30)+.45,
              ymax= seq(1:30)-.45,
              xmin=rev(Usable.from.date),
              xmax=rev(Usable.until.date),
              alpha=0.4))+
  ggtitle('2015')+
  theme_bw()+
  guides(alpha = FALSE)  +
  geom_rect(data=meta2b, aes(
    y=UnitLoc,
    ymin= seq(1:30)+.45,
    ymax= seq(1:30)-.45,
    xmin=rev(Usable.from.date),
    xmax=rev(Usable.until.date),
    alpha=0.4)) + 
  geom_jitter(data=trains2, 
              width = 0.0001,
              aes(y=UnitLoc, x=TrainDateNoTime, color=EncounterSpp)) +  
  geom_jitter(data=trains2, 
              aes(y=UnitLoc, x=TrainDateNoTime, color=EncounterSpp)) +
  scale_x_date(limits = c(as.Date("01/04/2015", format="%d/%m/%Y"),
                          as.Date("30/11/2015", format="%d/%m/%Y")))+
  scale_color_discrete(
    name="Encounter\nClassification",
    breaks=c("COD/BND", "UNK", "WBD/RSD"),
    labels=c("Broadband", "Unknown", "Frequency\nBanded")) +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    plot.title = element_text(hjust = 0.5)) 


ggsave("/media/kpalmer/My Passport/Figs For Aquatic Conservation/2015 All Occupancy .tiff",
       units="in", width=18, height=8, dpi=800, compression = 'lzw')




ggplot()+
  geom_rect(data=meta2, aes(
    y=UnitLoc,
    ymin= seq(1:30)+.45,
    ymax= seq(1:30)-.45,
    xmin=rev(Usable.from.date),
    xmax=rev(Usable.until.date)),
    alpha=0.4)+
  ggtitle('2015')+
  theme_bw()+
  geom_rect(data=meta2b, aes(
    y=UnitLoc,
    ymin= seq(1:30)+.45,
    ymax= seq(1:30)-.45,
    xmin=rev(Usable.from.date),
    xmax=rev(Usable.until.date)),
    alpha=0.4) + 
  scale_x_date(limits = c(as.Date("01/04/2015", format="%d/%m/%Y"),
                          as.Date("30/11/2015", format="%d/%m/%Y")))+
  geom_rect(data=trains3, aes(xmin = TrainDateNoTime,
                              xmax = TrainDateNoTime + 1, 
                              ymin = as.numeric(trains3$UnitLoc)-.45,
                              ymax = as.numeric(trains3$UnitLoc)+.45, 
                         fill=EncounterSpp), alpha=.7) +
  scale_fill_manual(values = c("#2B80C2","#a3a3c2", "#FFA42A"),
                    labels=c("Broadband","Unknown","Frequency Banded"),
                    name="Encounter Category") +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    plot.title = element_text(hjust = 0.5)) 
  
ggsave("/media/kpalmer/My Passport/Figs For Aquatic Conservation/2015 All Occupancy rect.tiff",
       units="in", width=18, height=8, dpi=800, compression = 'lzw')














  
  
  geom_jitter(data=trains2, 
              width = 0.0001,
              aes(y=UnitLoc, x=TrainDateNoTime, color=EncounterSpp)) +  
  geom_jitter(data=trains2, 
              aes(y=UnitLoc, x=TrainDateNoTime, color=EncounterSpp)) +
  scale_x_date(limits = c(as.Date("01/04/2015", format="%d/%m/%Y"),
                          as.Date("30/11/2015", format="%d/%m/%Y")))+
  scale_color_discrete(
    name="Encounter\nClassification",
    breaks=c("COD/BND", "UNK", "WBD/RSD"),
    labels=c("Broadband", "Unknown", "Frequency\nBanded")) +
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    plot.title = element_text(hjust = 0.5)) 

ggplot() + 
  geom_rect(data=meta_cpod, aes(
    y=UnitLoc,
    ymin= seq(1:30)+.45, ymax= seq(1:30)-.45,
    xmin=Usable.from.date, xmax=Usable.until.date), alpha=0.3) +
  
  geom_rect(data=mm, aes(xmin = TrainDateNoTime, xmax = TrainDateNoTime + 1, 
                         ymin = as.numeric(mm$UnitLoc)-.45, ymax = as.numeric(mm$UnitLoc)+.45, 
                         fill=EncounterSpp), alpha=.9) +
  geom_rect(data=meta_sm, aes(
    y=UnitLoc,
    ymin= c(3,5,8,12,15,18,20,22,25,30)+.45, ymax= c(3,5,8,12,15,18,20,22,25,30)-.45,
    xmin=Usable.from.date, xmax=Usable.until.date),
    size=.6, color="black", fill=NA) +
  
  scale_fill_manual(values = c("#2B80C2","#a3a3c2", "#FFA42A"),
                    labels=c("Broadband","Unknown","Frequency Banded"),
                    name="Encounter Category")+
  labs(y="", x="") +
  theme_bw()+
  theme(legend.position="none")









############################################################################################################
# 2015 data
# Color by tidal phase



ggplot()+geom_rect(data=meta2015a, aes(x=UnitLoc,
                                       xmin= seq(1:30)+.45,
                                       xmax= seq(1:30)-.45,
                                       ymin=Usable.from.date,
                                       ymax=Usable.until.date),
                   alpha=0.4) +
  geom_rect(data=meta2015b, aes(x=UnitLoc,
                                xmin= seq(1:30)+.45,
                                xmax= seq(1:30)-.45,
                                ymin=Usable.from.date,
                                ymax=Usable.until.date),
            alpha=0.4) +
  geom_jitter(data=Trains2015a, aes(x=UnitLoc, TrainDateNoTime, color=Phase)) 




Trains$DecimalHour=(Trains$MatlabDecDate-floor(Trains$MatlabDecDate))*24
Trains=Trains[sample(seq(1,nrow(Trains)), nrow(Trains), replace = FALSE),]

# Plot BND trains for all years as a function of phase
BNDTrains=subset(Trains, EncounterSpp=='COD/BND')
BNDTrains$Year=as.factor(BNDTrains$Year)
RSDTrains=subset(Trains, EncounterSpp=='WBD/RSD')

# All data, look at phase
ggplot()+geom_jitter(data=BNDTrains, aes(x=UnitLoc, DummyDate, color=Phase)) 

# All BND look at time of day
ggplot()+geom_jitter(data=BNDTrains, aes(x=UnitLoc, DummyDate, color=ToD)) 

ggplot()+geom_jitter(data=BNDTrains,
                     aes(y=UnitLoc, x=DecimalHour)) 

ggplot()+geom_jitter(data=RSDTrains,
                     aes(y=UnitLoc, x=DecimalHour)) 

#ggplot()+geom_density(data=BNDTrains, aes(y=UnitLoc, x=DecimalHour)) + scale_y_continuous('DecimalHour')


BNDTrains_nocro=subset(BNDTrains, UnitLoc !='Cro_05')
qplot(DecimalHour, data = BNDTrains_nocro, geom = "density", 
      fill=Lat, alpha=I(0.2), ylim = c(0,.25))+ ggtitle("BND/CPD")


RSDTrains_nocro=subset(RSDTrains, UnitLoc !='Cro_05')
qplot(DecimalHour, data = RSDTrains, geom = "density", 
      fill=Lat, alpha=I(0.2), ylim = c(0,.25))+ ggtitle("RSD")



ggplot()+geom_jitter(data=BNDTrains,
                     aes(x=UnitLoc, DecimalHour, color=Year))

# Look at cromarty 5
crodata=subset(Trains, UnitLoc=='Cro_05')

ggplot()+geom_jitter(data=crodata,
                     aes(x=Phase, y=log10(MedICI)))






RSDTrains=subset(Trains, EncounterSpp== "WBD/RSD")
ggplot()+geom_jitter(data=RSDTrains, aes(x=UnitLoc, DummyDate, color=ToD)) 
ggplot()+geom_jitter(data=RSDTrains, aes(x=UnitLoc, DummyDate, color=Phase))
ggplot()+geom_jitter(data=RSDTrains, aes(x=UnitLoc, DecimalHour)) 

hist(RSDTrains$DecimalHour)

# Not cromarty
hist(BNDTrains$DecimalHour[BNDTrains$UnitLoc != 'Cro_05'])
xx=seq(0,23, .01)
yy=1200+1000*cos((2*pi)/24*xx)
lines(xx, yy)

# Cromarty
hist(crodata$DecimalHour)
yy1=1200-700*cos((2*pi)/24*xx)
lines(xx, yy1)
######################################################################################################






boxplot(Trains$EncMinICI~Trains$UnitLoc, axt='n', xlab="", las=2)


BNDTrains=subset(Trains, EncounterSpp=="COD/BND")
WBDTrains=subset(Trains, EncounterSpp=="WBD/RSD")

boxplot(log10(BNDTrains$EncMinICI)~BNDTrains$UnitLoc, axt='n', xlab="", las=2)
boxplot(log10(BNDTrains$Enc)~BNDTrains$UnitLoc, axt='n', xlab="", las=2)
boxplot((BNDTrains$EncDur/60)~BNDTrains$UnitLoc, axt='n', xlab="", las=2)



#Broadband Dolphins
ggplot(data=BNDTrains, aes(UnitLoc, log10(MedICI),  color= log10(MedICI))) +
  geom_jitter()

# Broadband Banded Dolphins- color by month
ggplot(data=BNDTrains, aes(UnitLoc, log10(MedICI))) +
  geom_jitter(aes(color=as.factor(TrainMonth)))  + 
  scale_color_brewer(palette = 'greens')


# Broadband Banded Dolphins- color by month
ggplot(data=BNDTrains, aes(UnitLoc, DOY)) +
  geom_jitter(aes(color=as.factor(Year)))  + 
  scale_color_brewer(palette = 'greens')

  
# Frequency Banded Banded Dolphins- color by month
ggplot(data=WBDTrains, aes(UnitLoc, log10(MedICI))) +
  geom_jitter(aes(color=as.factor(TrainMonth)))  + 
  scale_color_brewer(palette = 'reds')


ggplot(data=subset(BNDTrains, Year==2013), aes(UnitLoc, TrainDateNoTime)) +
  geom_jitter(aes(color=log10(MedICI))) 
ggplot(data=subset(BNDTrains, Year==2014), aes(UnitLoc, TrainDateNoTime)) +
  geom_jitter(aes(color=log10(MedICI))) 




ggplot(data=subset(BNDTrains, Year==2014), aes(UnitLoc, MatlabDecDate)) +
  geom_jitter(aes(color=log10(MedICI)))

# Look at atime of year
ggplot(data= subset(BNDTrains, Year==2013), aes(UnitLoc, TrainDateNoTime)) +
  geom_jitter(aes(color=log10(MedICI))) 
scale_y_date(limits = c(as.Date('01/05/2013', "%d/%m/%Y"), as.Date('01/12/2013', "%d/%m/%Y")))

ggplot(data= subset(BNDTrains, Year==2014), aes(UnitLoc, TrainDateNoTime)) +
  geom_jitter(aes(color=log10(MedICI)))  +
  scale_y_date(limits = c(as.Date('01/05/2014', "%d/%m/%Y"), as.Date('01/12/2014', "%d/%m/%Y")))

# Look at time of year
ggplot(data= BNDTrains, aes(UnitLoc, DOY)) +
  geom_jitter(aes(color=log10(MedICI)))

ggplot(data= WBDTrains, aes(UnitLoc, DOY)) +
  geom_jitter(aes(color=log10(MedICI)))


# for WBD
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
ggplot(data= subset(WBDTrains, Year==2013), aes(UnitLoc, TrainDateNoTime)) +
  geom_jitter(aes(color=log10(MedICI))) 
 scale_y_date(limits = c(as.Date('01/05/2013', "%d/%m/%Y"), as.Date('01/12/2013', "%d/%m/%Y")))

ggplot(data= subset(WBDTrains, Year==2014), aes(UnitLoc, TrainDateNoTime)) +
  geom_jitter(aes(color=log10(MedICI)))  +
  scale_y_date(limits = c(as.Date('01/05/2014', "%d/%m/%Y"), as.Date('01/12/2014', "%d/%m/%Y")))



ggplot(data=subset(BNDTrains, Year==2013), aes(TrainDateNoTime, log10(MedICI))) +
  geom_jitter(aes(color=log10(MedICI))) 


# Look at distance from shore

ggplot(data= subset(BNDTrains, Year==2014), aes(ShoreDist, TrainDateNoTime)) +
  geom_jitter(aes(color=log10(MedICI)))
ggplot(data= subset(WBDTrains, Year==2013), aes(ShoreDist, TrainDateNoTime)) +
  geom_jitter(aes(color=log10(MedICI))) 

#############################################
# Fucking around #
Trains1=transform(Trains, gz=(as.numeric(EncounterSpp)-2)*2*(DaysSinceMay1*.5))
qplot(UnitLoc, log10(MedICI), data=Trains1, colour=gz) + scale_colour_gradient2()

ggplot(data=Trains1, aes(UnitLoc, log10(MedICI))) +
  geom_jitter(aes(color=gz))  + 
  scale_colour_gradient2()



###############################################





ggplot(Trains, aes(UnitLoc, log10(MedICI))) +
  geom_jitter(data = subset(Trains, EncounterSpp="COD/BND"))


# Frequency Banded Dolphins
ggplot(data=WBDTrains, aes(UnitLoc, log10(MedICI),  color= log10(MedICI))) +geom_jitter() + scale_fill_brewer(palette = "Greens")
# Frequency Banded Dolphins- color by month
ggplot(data=WBDTrains, aes(UnitLoc, log10(MedICI))) +
 geom_jitter(aes(color=TrainMonth)) +
 scale_fill_gradientn(colors=(blues9[1:length(blues9)]))



# All Dolphins
ggplot(data=Trains, aes(UnitLoc, log10(MedICI),  color= log10(MedICI))) +geom_jitter() + scale_fill_brewer(palette = "Greens")

# ICI values between dolphin spp
ggplot(Trains, aes(UnitLoc, log10(MedICI))) + 
  geom_jitter(aes(color=EncounterSpp), alpha=.3, stroke=0, size=2) + 
  scale_color_manual(values = c("orange", "gray","purple")) + 
  theme(legend.position=c(1,1),legend.justification=c(1,1)) 









# Pull out the unique encounters
EncTrains=Trains[!duplicated(Trains$EncounterID),]
###########################################################################################
# Plot the encounter duration, Proportion of low ICI trains for BND and WBD groups#
###########################################################################################


BNDEnc=subset(EncTrains, EncounterSpp=="COD/BND")
WBDEnc=subset(EncTrains, EncounterSpp=="WBD/RSD")

#Broadband Dolphins
ggplot(data=BNDEnc, aes(UnitLoc, log10(MedICI), color= log10(MedICI))) +geom_jitter() + scale_fill_brewer(palette = "Greens")
# Frequency Banded Dolphins
ggplot(data=WBDTrains, aes(UnitLoc, log10(MedICI),  color= log10(MedICI))) +geom_jitter() + scale_fill_brewer(palette = "Greens")

# All Dolphins
ggplot(data=Trains, aes(UnitLoc, log10(MedICI),  color= log10(MedICI))) +geom_jitter() + scale_fill_brewer(palette = "Greens")



##########################################################################################
# Look at the timing of the Encoutners #
##########################################################################################


# ICI values between dolphin spp
ggplot(subset(Trains, Year==2013), aes(UnitLoc, MatlabDecDate)) + 
  geom_jitter(aes(color=EncounterSpp)) + 
  scale_color_manual(values = c("orange", "gray","purple")) + 
  theme(legend.position=c(1,1),legend.justification=c(1,1)) 

ggplot(subset(Trains, Year==2014), aes(UnitLoc, MatlabDecDate+TrainTime/(24*60*60*1000))) + 
  geom_jitter(aes(color=EncounterSpp)) + 
  scale_color_manual(values = c("orange", "gray","purple")) + 
  theme(legend.position=c(1,1),legend.justification=c(1,1)) 


ggplot(Trains, aes(UnitLoc, MatlabDecDate+TrainTime/(24*60*60*1000))) + 
  geom_jitter(aes(color=EncounterSpp)) + 
  scale_color_manual(values = c("orange", "gray","purple")) + 
  theme(legend.position=c(1,1),legend.justification=c(1,1)) 




########################################################################################################
# Look at the ratio of low ICI click trains to total number of click trains
###################################################################################

# Number of low ICI trains
  low_iciBND=subset(BNDTrains, log10(MedICI)<=(-2.25))
  low_iciWBD=subset(WBDTrains, log10(MedICI)<=(-2.25))

# Number of trains total
  n_BND=table(BNDTrains$UnitLoc)
  n_WBD=table(WBDTrains$UnitLoc)

#
  
# Ratio of low ICI trains to the Rest of the Trains
  Prop_BND=table(low_iciBND$UnitLoc)/table(BNDTrains$UnitLoc)
  Prop_WBD=table(low_iciWBD$UnitLoc)/table(WBDTrains$UnitLoc)
  
 # Plot the ratio of low ICI trains to the rest of the trains (regular and log)
   plot(as.matrix(n_BND), as.matrix(table(low_iciBND$UnitLoc)), 
        xlab='(n) Broadband trains', ylab='Number of Buzzes')
   
   plot(as.matrix(n_WBD), as.matrix(table(low_iciWBD$UnitLoc)), 
        xlab='(n) Broadband trains', ylab='Number of Buzzes')
   
   
   
   
 # Plot the number of buzzes at each location
   plot(as.matrix(table(low_iciBND$UnitLoc)), xaxt = "n", xlab="", pch=19, col='Maroon', ylab='Number of Low ICI')
   points(as.matrix(table(low_iciWBD$UnitLoc)), xaxt = "n", xlab="", pch=19, col='CornflowerBlue')
   axis(1, at=1:length(level_names), labels=level_names, las=2)
   legend("topright", legend=c("Broadband", "Frequency Banded"), pch=19, col=c('Maroon', 'CornflowerBlue'))
   
   
   
   
   
   
   plot(as.matrix(n_BND), as.matrix(Prop_BND), 
        xlab='(n) Broadband trains', ylab='Prop <10e-2')
   
#  plot(as.matrix(n_BND), as.matrix(Prop_BND), log='x', 
#      xlab='log (n) Broadband trains', ylab='Prop <10e-2')
#   plot(as.matrix(n_BND), (as.matrix(Ratio_BND)), log="xy", 
#        xlab='log(n) Broadband trains', ylab='log Prop <10e-2')

kk=lm(as.matrix(Prop_BND)~log(as.matrix(n_BND))) # R^2 is shit, no relationship!
ll=lm(as.matrix(Prop_WBD)~log(as.matrix(n_WBD))) # R^2 is ALSO shit, no relationship!

plot(as.matrix(Prop_BND), xaxt = "n", xlab="", pch=19, col='Maroon', ylab='Proportion Low ICI')
points(as.matrix(Prop_WBD), xaxt = "n", xlab="", pch=19, col='CornflowerBlue')
axis(1, at=1:length(level_names), labels=level_names, las=2)
legend("topright", legend=c("Broadband", "Frequency Banded"), pch=19, col=c('Maroon', 'CornflowerBlue'))


ggplot(data=BNDEnc, aes(TrainMonth, log10(EncLrICIQuartile)))+ geom_jitter()

