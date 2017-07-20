


rm(list=ls())
library(ggplot2)
library(rethinking)

setwd("W:/KJP PHD/5-Pp Tt interactions/R Code")


# colorblind palette with black for plotting
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



# Load Dolphin Clicks ################################################################################################



Trains2013=read.csv("W:/KJP PHD/4-Bayesian Habitat Use/DolCPODFiles/2013_ClickTrainAnalysis.csv")
Trains2014=read.csv("W:/KJP PHD/4-Bayesian Habitat Use/DolCPODFiles/2014_ClickTrainAnalysis.csv")
Trains2015=read.csv("W:/KJP PHD/4-Bayesian Habitat Use/DolCPODFiles/2015_ClickTrainAnalysis.csv")

Trains2013$TrainDay= as.numeric(substr((Trains2013$TrainDateNoTime), 1,2))
Trains2014$TrainDay= as.numeric(substr((Trains2014$TrainDateNoTime), 1,2))

Trains2013$Year=2013
Trains2014$Year=2014
Trains2015$Year=2015


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

Trains2015$UnitLoc=factor(Trains2015$UnitLoc, levels=level_names)
Trains2014$UnitLoc=factor(Trains2014$UnitLoc, levels=level_names)
Trains2013$UnitLoc=factor(Trains2013$UnitLoc, levels=level_names)


Trains2015$TrainDateNoTime= as.Date(Trains2015$TrainDateNoTime, "%d/%m/%Y", tz = "GMT")
Trains2014$TrainDateNoTime= as.Date(Trains2014$TrainDateNoTime, "%d/%m/%Y", tz = "GMT")
Trains2015$TrainDateNoTime= as.Date(Trains2015$TrainDateNoTime, "%d/%m/%Y", tz = "GMT")


Trains2014$EncounterID=Trains2014$EncounterID+max(Trains2013$EncounterID)
Trains2015$EncounterID=Trains2015$EncounterID+max(Trains2014$EncounterID)

Trains=rbind(Trains2013, Trains2014, Trains2015 )

rm(Trains2013, Trains2014, Trains2015)




# Load Porposie Clicks # #######################################################################################
 # This works! Export Train details for NBHF
 PP2013_trains=read.table("W:/KJP PHD/5-Pp Tt interactions/R Code/2013Por_Train_details.txt", 
                          sep="\t",  row.names=NULL, header=TRUE)
 PP2013_trains$Year=2013
 PP2014_trains=read.table("W:/KJP PHD/5-Pp Tt interactions/R Code/2014Por_Train_details.txt", 
                          sep="\t",  row.names=NULL, header=TRUE)
 PP2014_trains$Year=2014
 
 PP2015A_trains=read.table("W:/KJP PHD/5-Pp Tt interactions/R Code/2015APor_Train_details.txt", 
                          sep="\t",  row.names=NULL, header=TRUE)
 PP2015A_trains$Year=2015
 
 PP2015B_trains=read.table("W:/KJP PHD/5-Pp Tt interactions/R Code/2015BPor_Train_details.txt", 
                           sep="\t",  row.names=NULL, header=TRUE)
 PP2015B_trains$Year=2015
 
 
 

HPTrains=rbind(PP2013_trains, PP2014_trains, PP2015A_trains)
HPTrains$UnitLoc=as.factor(substr(as.character(HPTrains$row.names), 1,6))
HPTrains$TrainDateNoTime= as.Date(HPTrains$TrnID, "%d/%m/%Y")
HPTrains$MatlabDecDate= Datenum(HPTrains$TrainDateNoTime) +
  as.numeric(substr((HPTrains$TrnID[1]), 11,12))/24+
  as.numeric(substr((HPTrains$TrnID[1]), 14,15))/24/60
HPTrains$EncounterSpp='HP'



# Create HP Encounters ###

HPTrains$EncounterID=array(NA, dim=nrow(HPTrains))
EncounterVal=10
MaxETime=5/60/24
MaxETime_min=.05


HPTrains$EncounterID[1]=1
HPTrains$EncounterDuration=0

  
EncounterID_out <- vector(mode="numeric", length=0)




for (ii in 1:length(unique(HPTrains$UnitLoc))){
  
  diff_vals=diff(HPTrains$Time[HPTrains$UnitLoc==unique(HPTrains$UnitLoc)[ii]])
  encounterIDX=c(1, which(diff_vals>=MaxETime_min))
  EncounterIDs=c(EncounterVal, diff_vals*0)
  
  
  for (jj in 1:(length(encounterIDX))){
    
    
    
    if (jj==length(encounterIDX)){
      
      EncounterIDs[encounterIDX[jj]:length(EncounterIDs)]=EncounterVal
                   
      }else{
      
      EncounterIDs[encounterIDX[jj]:encounterIDX[jj+1]-1]=EncounterVal
    }
    
    EncounterVal=EncounterVal+1
  }
  
  
  EncounterID_out=c(EncounterID_out, EncounterIDs)
  
}

HPTrains$EncounterID=EncounterID_out



rm(EncounterID_out, EncounterVal, encounterIDX, diff_vals, EncounterIDs, ii, jj, 
   PP2013_trains, PP2014_trains,PP2015A_trains,PP2015B_trains)

# Add Dates in Date Format


## Pre-processing ##########################################################################################################

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

HPTrains$UnitLoc=factor(HPTrains$UnitLoc, levels=level_names)
HPTrains$Hour=as.numeric(substr(HPTrains$TrnID, 11,12))

# Min (here time) is number of minutes since 1900
# Start time within minute (here Trclass).
# In MATLAB The datenum function creates a numeric array that represents 
# each point in time as the number of days from January 0, 0000. 

# Convert C-POD datetime to matlab date time

# 1 CPOD time in days
CPODTime=(HPTrains$Time+HPTrains$TrClass/6e7)/1440


# Matlabdatenum for 1900 is  693962 
HPTrains$MatlabDecDate=CPODTime+693962 

rm(CPODTime)


HPTrains$EncounterSpp

meta=read.csv('W:/KJP PHD/Deployment Information/CPODs for Kaitlin.csv')
meta$Usable.from.date=as.Date(meta$Usable.from.date, "%d/%m/%Y", tz = "GMT")
meta$Usable.until.date=as.Date(meta$Usable.until.date, "%d/%m/%Y", tz = "GMT")
meta$Usable.Days=difftime(meta$Usable.until.date, meta$Usable.from.date,units='days')
meta$Year=substr(meta$Usable.from.date,1,4)
meta=meta[!is.na(meta$Usable.Days),]

# Get meta start and end to matlab dates, go through C-POD times because 000000 isn't really a thing
meta$MatlabStart=difftime(meta$Usable.from.date,  as.Date('1900/01/01'))+693962
meta$MatlabEnd=difftime(meta$Usable.until.date,  as.Date('1900/01/01'))+693962



# DO NOT MERGE THE TABLE, INSUFFICIENT RAM
# AllTrains=merge(Trains, HPTrains, by = c('UnitLoc', 'Year'), all.x = TRUE)


## Select Units and look at random distributions ######################

# Add UnitLoc year for indexing
HPTrains$UnitLocYear=paste(HPTrains$UnitLoc, HPTrains$Year)
Trains$UnitLocYear=paste(Trains$UnitLoc, Trains$Year)
meta$UnitLocYear=paste(meta$UnitLoc, meta$Year)


# For dolphin trains make a new table for the encounters and determine the 
# mean encounter time

OthCeteEncounters=aggregate(data=Trains, MatlabDecDate~EncounterID+UnitLoc+Year+EncounterSpp+UnitLocYear,
                            FUN=function(x){min(x)+((max(x)-min(x))/2)})
colnames(OthCeteEncounters)[ncol(OthCeteEncounters)]='EncounterMid'

# Encounter duration in minutes
OthCeteEncounters$EncounterDur=aggregate(data=Trains, MatlabDecDate~EncounterID+UnitLoc+Year+EncounterSpp+UnitLocYear,
                                        FUN=function(x){(max(x)-min(x))*24*60})[6]

# Select only encounters greater than 5 min
OthCeteEncounters=OthCeteEncounters[OthCeteEncounters$EncounterDur>5,]

CetaceanWaitingTimes=data.frame(Dur=numeric(length=0),
                             UnitLoc=character(length=0),
                             Year=numeric(length=0),
                             EncounterSpp=character(length=0),
                             id=numeric(length=0))

id=1


# Calculate waiting times when dolphin is present and waiting times in the three days
# prior to a dolphin detection
for(ii in 1:length(unique(Trains$UnitLocYear))){

  
  Unit=unique(Trains$UnitLocYear)[ii]
  Cetsub=OthCeteEncounters[OthCeteEncounters$UnitLocYear==Unit,]
  
  
  
  if(isTRUE(nrow(Cetsub)>0)){
  
  PPsub=HPTrains[HPTrains$UnitLocYear==Unit,]
  metasub=meta[meta$UnitLocYear==Unit,]

  
  # Get differences from the midpoint of the encounter until the next PP encounter
  DolWait=numeric(length = nrow(Cetsub))
 
  if(metasub$Year[1]==2015){
    UsableDates=seq(min(metasub$MatlabStart), max(metasub$MatlabEnd))}else{
    UsableDates=seq(metasub$MatlabStart, metasub$MatlabEnd)
  }
  
  
  daily_occ=data.frame(UsableDates= UsableDates)
  DolDaily=aggregate(EncounterSpp~round(EncounterMid), data=Cetsub, FUN=length)
  colnames(DolDaily)[1]='UsableDates'
  daily_occ=merge(daily_occ, DolDaily, all.x = TRUE, by='UsableDates')
  daily_occ[is.na(daily_occ)]=0
  
  # 'Baseline occurrence at each site was characterized by randomly selecting 100 control
  #  points from the week prior to the seismic survey and calculating the waiting times 
  #  from these points to the next porpoise detection'
  
  
    for(jj in 1:nrow(Cetsub)){
      
      DolDetDate=round(Cetsub$EncounterMid[jj])
      
      # Get waiting time from dolpin to next encounter
      aa=PPsub$MatlabDecDate-Cetsub$EncounterMid[jj]
      DolWait[jj]=min(aa[aa>0])
      
      # Get random waiting times from three days before the encounter
      dateidx=which(daily_occ$UsableDates == DolDetDate)
      
      if(isTRUE(dateidx>3) & sum(aa>0)>0){
       
        # If sufficient gap, then grab 30 random points from the prvious days
          if(all(daily_occ$EncounterSpp[(dateidx-3):(dateidx-1)]==0)){
      
            # Dates from which control periods are drawn
            ControlDates=daily_occ$UsableDates[(dateidx-3):(dateidx-1)]
            
            # Control Times
            ControlTimes=runif(n = 30, min = min(ControlDates), max = max(ControlDates))
            
            ControlWait=numeric(length=30)
            
            # Time difference in days
            for(kk in 1:30){
              aa=PPsub$MatlabDecDate-ControlTimes[kk]
              ControlWait[kk]=min(aa[aa>0])
              rm(aa)
            }
            
            controldf= data.frame(Dur=ControlWait,
                                  UnitLoc=substr(Unit, 1,6),
                                  Year=PPsub$Year[1],
                                  EncounterSpp='Control',
                                  id=id)
            
            doldf=data.frame(Dur=DolWait[jj],
                       UnitLoc=substr(Unit, 1,6),
                       Year=PPsub$Year[1],
                       EncounterSpp=Cetsub$EncounterSpp,
                       id=id)
            
            CetaceanWaitingTimes=rbind(CetaceanWaitingTimes,controldf, doldf)
            id=id+1
          }
      }
      
    }
    
  }
 print(paste(ii, 'of', length(unique(Trains$UnitLocYear))))
}
  


# Fix level orders
CetaceanWaitingTimes$EncounterSpp=factor(CetaceanWaitingTimes$EncounterSpp,
                                         levels=c('Control', 'COD/BND', 'WBD/RSD', 'UNK'))

CetaceanWaitingTimes$AllDol=factor(CetaceanWaitingTimes$AllDol,
                                         levels=c('Control', 'Dolphin'))
CetaceanWaitingTimes$GroupId=substr(CetaceanWaitingTimes$UnitLoc, 1,3)
CetaceanWaitingTimes$ShoreDist=substr(CetaceanWaitingTimes$UnitLoc, 5,6)
CetaceanWaitingTimes$AllDol= ifelse(CetaceanWaitingTimes$EncounterSpp=='Control',
                                    'Control', 'Dolphin')


# Remove infinite and NA cases
CetaceanWaitingTimes=CetaceanWaitingTimes[complete.cases(CetaceanWaitingTimes),]
CetaceanWaitingTimes=CetaceanWaitingTimes[!is.infinite(CetaceanWaitingTimes$Dur),]

# Make plots for 10 deployment locations and histogram to match Tompson et al. #########################################


# Calculate the duration in minutes and round
CetaceanWaitingTimes$Dur_min=round(CetaceanWaitingTimes$Dur*24*60)


p_all=list()
p_spp=list()

for(ii in 1:10){
  data_sub=subset(CetaceanWaitingTimes, GroupId==unique(CetaceanWaitingTimes$GroupId)[ii])
  data_sub=data_sub[complete.cases(data_sub),]
  data_sub=data_sub[!is.infinite(data_sub$Dur),]
  
  p_spp[[ii]]=
    ggplot(data_sub, aes(Dur, ..density.., colour = EncounterSpp)) +
    geom_freqpoly(binwidth = .1) +
    facet_wrap(~ShoreDist)+
    scale_color_manual(values=cbbPalette) +
    ggtitle(data_sub$GroupId[1])+
    theme_light()
  
  p_all[[ii]]=
  ggplot(data_sub, aes(Dur, ..density.., colour = AllDol)) +
    geom_freqpoly(binwidth = .1) +
    facet_wrap(~ShoreDist) +
    scale_color_manual(values=cbbPalette) +
    ggtitle(data_sub$GroupId[1])+
    theme_light()
    
  
}

ggplot(data=CetaceanWaitingTimes, aes(x=Dur_min, group=AllDol, fill=AllDol)) +
  geom_histogram(position="dodge", bins = 50) + 
  scale_x_continuous(limits = c(0, 24*60*3)) +
  theme_bw()

ggplot(data=subset(CetaceanWaitingTimes, EncounterSpp != 'UNK'),
        aes(x=Dur_min, group=EncounterSpp, fill=EncounterSpp)) +
  geom_histogram(position="dodge", bins = 50) + 
  scale_x_continuous(limits = c(0, 24*60*3)) +
  theme_bw()

ggplot(data=CetaceanWaitingTimes, aes(x=Dur_min, group=AllDol, fill=AllDol)) +
  geom_density(position="dodge", alpha=.4) + 
  scale_x_continuous(limits = c(0, 24*60*3)) +
  theme_bw()

ggplot(data=subset(CetaceanWaitingTimes, EncounterSpp != 'UNK'),
       aes(x=Dur_min, group=EncounterSpp, fill=EncounterSpp)) +
  geom_density(position="dodge") + 
  scale_x_continuous(limits = c(0, 24*60*3)) +
  theme_bw()

# Build linear models ####################################################################################





mod=glm(Dur~EncounterSpp, 
         family.glmm = negative.binomial(2,link = "log"),
         data=CetaceanWaitingTimes)

m1 <- glm.nb(Dur ~ EncounterSpp + ShoreDist + Year, data = CetaceanWaitingTimes)
m3 <- glm(Dur ~ EncounterSpp + ShoreDist + Year+UnitLoc, family = "poisson", data = CetaceanWaitingTimes)

m3 <- glmer(Dur ~ EncounterSpp+(1|ShoreDist)+(1|GroupId), 
            family =poisson, data = CetaceanWaitingTimes)







data$obs_effect<-1:nrow(data)
overdisp.fit<-lmer(y~1+obs_effect+x+(1|obs_effect)+(1+x|subject_id),data=data,family=poisson)



m2 <- update(m1, . ~ . - prog)
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)



# in progress ###########################################################################








############################################################################################
ggplot(mtlong, aes(value)) + facet_wrap(~variable, scales = 'free_x') +
  geom_histogram(binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)))




DolEncounters_min=aggregate(data=DolTrains, MatlabDecDate~EncounterID+UnitLoc+EncounterSpp+Year, FUN=min)
DolEncounters_max=cbind(DolEncounters_min, aggregate(data=DolTrains, MatlabDecDate~EncounterID+UnitLoc+EncounterSpp+Year, FUN=max)[,5])
PorEncounters_min=aggregate(data=HPTrains,  MatlabDecDate~EncounterID+UnitLoc+EncounterSpp+Year, FUN=min)
PorEncounters_max=cbind(PorEncounters_min, aggregate(data=HPTrains,  MatlabDecDate~EncounterID+UnitLoc+EncounterSpp+Year, FUN=max)[,5])

colnames(DolEncounters_max)[5:6]=c('Encounter_Start','Encounter_End')
colnames(PorEncounters_max)[5:6]=c('Encounter_Start','Encounter_End')

NewTable=rbind(DolEncounters_max, PorEncounters_max)

NewTable$EncounterDuration_sec= (NewTable$Encounter_End-NewTable$Encounter_Start)*24*60*60
NewTable$EncounterDateStart=as.POSIXct((NewTable$Encounter_Start - 719529)*86400, origin = "1970-01-01", tz = "UTC")
NewTable$EncounterDateEnd  =as.POSIXct((NewTable$Encounter_End - 719529)*86400, origin = "1970-01-01", tz = "UTC")
NewTable$EncounterDateStart=as.Date(NewTable$EncounterDateStart, '%Y-%m-%d', tz = 'UTC')
NewTable$EncounterDateEnd=as.Date(NewTable$EncounterDateEnd, '%Y-%m-%d', tz = 'UTC')

# Sort according to date and unit Loc# 
NewTable=NewTable[order(NewTable$UnitLoc, NewTable$Encounter_Start),]


# Trim the dataset #
meta=read.csv('W:/KJP PHD/Deployment Information/CPODs for Kaitlin.csv')
meta[meta[,]=='']=NA
meta=meta[complete.cases(meta$Usable.from.date),]
meta$Usable.from.date= as.Date(meta$Usable.from.date, "%d/%m/%Y")
meta$Usable.until.date= as.Date(meta$Usable.until.date, "%d/%m/%Y")
meta$Year=as.numeric(substr(meta$Usable.from.date, 1,4))

Nix=numeric(length = 0)
for(ii in 1:nrow(meta)){
  
  # get the approperiate unit and year
  idx=which(NewTable$UnitLoc==meta$UnitLoc[ii] & NewTable$Year==meta$Year[ii])
  
  # Find any trains that are earlier than the start date or later than the end date
  Nix1=idx[which(NewTable$EncounterDateStart[idx]<meta$Usable.from.date[ii])]
  Nix=c(Nix, idx[which(NewTable$EncounterDateEnd[idx]>meta$Usable.until.date[ii])])
  
  
}

NewTable1=NewTable[-Nix,]
NewTable1$IsDolphin=ifelse(NewTable1$EncounterSpp=='HP', 'HP', 'DOL')
NewTable1$GroupId=as.factor(unlist(strsplit(as.character(NewTable1$UnitLoc), split = "_"))
                            [seq(1,(nrow(NewTable1)*2)-1,2)])



########################################################################################

# Data Vis # 

library(ggplot2)

mm_sub=subset(NewTable1, Year==2015 & GroupID==unique(NewTable1$GroupId)[1])
mm_sub=NewTable1[NewTable1$GroupId==unique(NewTable1$GroupId)[4] & 
                   NewTable1$Year==2015,]

# 
# png("W:\\KJP PHD\\5-Pp Tt interactions\\Figures and Tables\\PPTTAll.png", 
#     width = 8, height = 24, units = 'in', res = 300)
ggplot(mm_sub, aes(color=IsDolphin)) +
  geom_segment(aes(x=Encounter_Start, xend=Encounter_End, 
                   y=IsDolphin, yend=IsDolphin), size=15) +
  facet_wrap(~UnitLoc, ncol = 1)



dev.off()


ggplot(mm_sub, aes(color=IsDolphin)) +
  geom_point(aes(x=Encounter_Start, y=IsDolphin))
  geom_jitter()

  
  
  
               
               
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








