

# This code creates hourly occupancy and poisson data for full 
# C-POD detections as well as combining covariate data including
# Tidal data and sun elevation



rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(EMD)
library(solidearthtide) #For datenum, convert date to matlabdate
setwd("W:/KJP PHD/4-Bayesian Habitat Use/R Code")

### Sun Elevation and Azmouth Function ###
sunPosition <- function(year, month, day, hour=12, min=0, sec=0,
                        lat=46.5, long=6.5) {
  
  
  twopi <- 2 * pi
  deg2rad <- pi / 180
  
  # Get day of the year, e.g. Feb 1 = 32, Mar 1 = 61 on leap years
  month.days <- c(0,31,28,31,30,31,30,31,31,30,31,30)
  day <- day + cumsum(month.days)[month]
  leapdays <- year %% 4 == 0 & (year %% 400 == 0 | year %% 100 != 0) & day >= 60
  day[leapdays] <- day[leapdays] + 1
  
  # Get Julian date - 2400000
  hour <- hour + min / 60 + sec / 3600 # hour plus fraction
  delta <- year - 1949
  leap <- trunc(delta / 4) # former leapyears
  jd <- 32916.5 + delta * 365 + leap + day + hour / 24
  
  # The input to the Atronomer's almanach is the difference between
  # the Julian date and JD 2451545.0 (noon, 1 January 2000)
  time <- jd - 51545.
  
  # Ecliptic coordinates
  
  # Mean longitude
  mnlong <- 280.460 + .9856474 * time
  mnlong <- mnlong %% 360
  mnlong[mnlong < 0] <- mnlong[mnlong < 0] + 360
  
  # Mean anomaly
  mnanom <- 357.528 + .9856003 * time
  mnanom <- mnanom %% 360
  mnanom[mnanom < 0] <- mnanom[mnanom < 0] + 360
  mnanom <- mnanom * deg2rad
  
  # Ecliptic longitude and obliquity of ecliptic
  eclong <- mnlong + 1.915 * sin(mnanom) + 0.020 * sin(2 * mnanom)
  eclong <- eclong %% 360
  eclong[eclong < 0] <- eclong[eclong < 0] + 360
  oblqec <- 23.429 - 0.0000004 * time
  eclong <- eclong * deg2rad
  oblqec <- oblqec * deg2rad
  
  # Celestial coordinates
  # Right ascension and declination
  num <- cos(oblqec) * sin(eclong)
  den <- cos(eclong)
  ra <- atan(num / den)
  ra[den < 0] <- ra[den < 0] + pi
  ra[den >= 0 & num < 0] <- ra[den >= 0 & num < 0] + twopi
  dec <- asin(sin(oblqec) * sin(eclong))
  
  # Local coordinates
  # Greenwich mean sidereal time
  gmst <- 6.697375 + .0657098242 * time + hour
  gmst <- gmst %% 24
  gmst[gmst < 0] <- gmst[gmst < 0] + 24.
  
  # Local mean sidereal time
  lmst <- gmst + long / 15.
  lmst <- lmst %% 24.
  lmst[lmst < 0] <- lmst[lmst < 0] + 24.
  lmst <- lmst * 15. * deg2rad
  
  # Hour angle
  ha <- lmst - ra
  ha[ha < -pi] <- ha[ha < -pi] + twopi
  ha[ha > pi] <- ha[ha > pi] - twopi
  
  # Latitude to radians
  lat <- lat * deg2rad
  
  # Azimuth and elevation
  el <- asin(sin(dec) * sin(lat) + cos(dec) * cos(lat) * cos(ha))
  az <- asin(-cos(dec) * sin(ha) / cos(el))
  elc <- asin(sin(dec) / sin(lat))
  az[el >= elc] <- pi - az[el >= elc]
  az[el <= elc & ha > 0] <- az[el <= elc & ha > 0] + twopi
  
  el <- el / deg2rad
  az <- az / deg2rad
  lat <- lat / deg2rad
  
  return(list(elevation=el, azimuth=az, JulienDay=day))
}


### Merge with Order Function ####
merge.with.order <- function(x,y, ..., sort = T, keep_order)
{
  # this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
  add.id.column.to.data <- function(DATA)
  {
    data.frame(DATA, id... = seq_len(nrow(DATA)))
  }
  # add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
  order.by.id...and.remove.it <- function(DATA)
  {
    # gets in a data.frame with the "id..." column.  Orders by it and returns it
    if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")
    
    ss_r <- order(DATA$id...)
    ss_c <- colnames(DATA) != "id..."
    DATA[ss_r, ss_c]
  }
  
  # tmp <- function(x) x==1; 1	# why we must check what to do if it is missing or not...
  # tmp()
  
  if(!missing(keep_order))
  {
    if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
    if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
    # if you didn't get "return" by now - issue a warning.
    warning("The function merge.with.order only accepts NULL/1/2 values for the keep_order variable")
  } else {return(merge(x=x,y=y,..., sort = sort, all.x=T))}
}






#########################################################################################################
# Load Deployment Info and Processed Train Data #
#########################################################################################################

#Only CPOD info
DepInfo=read.csv('W:/KJP PHD/Deployment Information/Copy of Marine Scotland 2013-14 CPOD summary data (2).csv')

PdetInfo=read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/MaxPdet.csv')


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
Trains2013$TrainDateNoTime= as.Date(Trains2013$TrainDateNoTime, "%d/%m/%Y", tz = "GMT")


Trains2014$EncounterID=Trains2014$EncounterID+max(Trains2013$EncounterID)

Trains=rbind(Trains2013, Trains2014, Trains2015 )

##################################################################################
# Load the C-POD deployment Information #
##################################################################################



meta=read.csv('W:/KJP PHD/Deployment Information/CPODs for Kaitlin.csv')
meta$Usable.from.date=as.Date(meta$Usable.from.date, "%d/%m/%Y", tz = "GMT")
meta$Usable.until.date=as.Date(meta$Usable.until.date, "%d/%m/%Y", tz = "GMT")
meta$Usable.Days=difftime(meta$Usable.until.date, meta$Usable.from.date,units='days')
meta$Year=substr(meta$Usable.from.date,1,4)
meta=meta[!is.na(meta$Usable.Days),]

# Extra information for lat/lon
#meta=read.csv('W:/KJP PHD/CPOD Processing/2013 to 2016 CPOD deployments.csv')

##########################################################################################
# Load detection probabiliies sheets #
#############################################################################################

# Load the Pdet Files
Pdet1 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2013/Lat_05_NLandPdetVals.csv', header = TRUE)
Pdet1$UnitLoc="Lat_05"
Pdet2 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2013/Hel_15_NLandPdetVals.csv', header = TRUE)
Pdet2$UnitLoc="Hel_15"
Pdet3 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2013/Cro_15_NLandPdetVals.csv', header = TRUE)
Pdet3$UnitLoc="Cro_15"
Pdet4 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2013/SpB_10_NLandPdetVals.csv', header = TRUE)
Pdet4$UnitLoc="SpB_10"
Pdet5 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2013/Fra_05_NLandPdetVals.csv', header = TRUE)
Pdet5$UnitLoc="Fra_05"
Pdet6 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2013/Cru_05_NLandPdetVals.csv', header = TRUE)
Pdet6$UnitLoc="Cru_05"
Pdet7 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2013/Sto_05_NLandPdetVals.csv', header = TRUE)
Pdet7$UnitLoc="Sto_05"
Pdet8 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2013/Arb_10_NLandPdetVals.csv', header = TRUE)
Pdet8$UnitLoc="Arb_10"
Pdet9 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2013/StA_10_NLandPdetVals.csv', header = TRUE)
Pdet9$UnitLoc="StA_10"
Pdet10=read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2013/Stb_05_NLandPdetVals.csv', header = TRUE)
Pdet10$UnitLoc="Stb_05"

Pdetdf=rbind(Pdet1,Pdet2,Pdet3,Pdet4,Pdet5,Pdet6,Pdet7,Pdet8,Pdet9,Pdet10)
Pdetdf$Year=2013

# Fix the Arbroath name
Pdetdf$UnitLoc[Pdetdf$UnitLoc=="Arb_10"]="Abr_10"

# Add date
Pdetdf$Date=as.Date(Pdetdf$MatlabDate, origin = '0000-01-01')
Pdetdf$Hr=round((Pdetdf$MatlabDate-floor(Pdetdf$MatlabDate))*24)

# Add Beta Prior for standard deviation of the detection probability based on the suspected 
# model error


rm(Pdet1,Pdet2,Pdet3,Pdet4,Pdet5,Pdet6,Pdet7,Pdet8,Pdet9,Pdet10)

##################################
# Do the same thing for 2014 data
##################################


# Load the Pdet Files
Pdet1 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2014/Lat_05_NLandPdetVals.csv', header = TRUE)
Pdet1$UnitLoc="Lat_05"
Pdet2 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2014/Hel_15_NLandPdetVals.csv', header = TRUE)
Pdet2$UnitLoc="Hel_15"
Pdet3 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2014/Cro_15_NLandPdetVals.csv', header = TRUE)
Pdet3$UnitLoc="Cro_15"
Pdet4 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2014/SpB_10_NLandPdetVals.csv', header = TRUE)
Pdet4$UnitLoc="SpB_10"
Pdet5 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2014/Fra_05_NLandPdetVals.csv', header = TRUE)
Pdet5$UnitLoc="Fra_05"
Pdet6 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2014/Cru_05_NLandPdetVals.csv', header = TRUE)
Pdet6$UnitLoc="Cru_05"
# No Stonehaven in 2014
#Pdet7 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2014/Sto_05_NLandPdetVals.csv', header = TRUE)
#Pdet7$UnitLoc="Sto_05"
Pdet8 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2014/Arb_10_NLandPdetVals.csv', header = TRUE)
Pdet8$UnitLoc="Arb_10"
Pdet9 =read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2014/StA_10_NLandPdetVals.csv', header = TRUE)
Pdet9$UnitLoc="StA_10"
Pdet10=read.csv('W:/KJP PHD/4-Bayesian Habitat Use/Pdet at Time/2014/Stb_05_NLandPdetVals.csv', header = TRUE)
Pdet10$UnitLoc="Stb_05"

Pdetdf1=rbind(Pdet1,Pdet2,Pdet3,Pdet4,Pdet5,Pdet6,Pdet8,Pdet9,Pdet10)
Pdetdf1$Year=2014

# Fix the Arbroath name
Pdetdf1$UnitLoc[Pdetdf1$UnitLoc=="Arb_10"]="Abr_10"
Pdetdf1$Date=as.Date(Pdetdf1$MatlabDate, origin = '0000-01-01')
Pdetdf1$Hr=round((Pdetdf1$MatlabDate-floor(Pdetdf1$MatlabDate))*24)



Pdetdf=rbind(Pdetdf, Pdetdf1)


Pdetdf$UnitLoc=factor(Pdetdf$UnitLoc, levels=level_names)

rm(Pdet1,Pdet2,Pdet3,Pdet4,Pdet5,Pdet6,Pdet8,Pdet9,Pdet10, Pdetdf1)



##########################################################################################
# Make giant hourly occupancy table #
##########################################################################################

OccTable=data.frame(Date=factor(), UnitLoc=factor(), Hr=numeric(), 
                    Slope=numeric(), ShoreDist=numeric(), Lat=numeric(), Lon=numeric())
meta=subset(meta, Usable.Days>0)

for(ii in 1:nrow(meta)){
  
  ## Create the Hourly Occupancy Table ##
  tempdf=data.frame(Date=sort(rep(seq.Date(from=meta$Usable.from.date[ii], 
                          to= meta$Usable.until.date[ii], 
                          by = 'days'), 24)))
  tempdf$UnitLoc=meta$UnitLoc[ii]
  tempdf$Hr=as.numeric(rep(0:23, nrow(tempdf)/24))
  tempdf$MatlabDate=Datenum(as.Date(tempdf$Date[1], format= '%Y-%m/%d'))
  tempdf$Year=as.numeric(meta$Year[ii])
  tempdf$Slope=meta$Slope[ii]
  tempdf$ShoreDist=as.numeric(substr(meta$UnitLoc[ii],5,7))
  tempdf$Lat=meta$Lat..dec.deg.[ii]
  tempdf$Lon=meta$Long..dec.deg.[ii]

  OccTable=rbind(tempdf, OccTable)
  rm('tempdf')
}


OccTable$MatlabDate=Datenum(OccTable$Date)+OccTable$Hr/24

OccTable$GroupId=substr(OccTable$UnitLoc, start = 1, stop = 3)
OccTable$ShoreDist=unlist(strsplit(as.character(OccTable$UnitLoc), "_"))[2]

OccTable=OccTable[!duplicated(OccTable[,1:3]),]

####################################################################################
# Load and Merge the detection probability sheets with the full dataframe #
################################################################################

  # 2013 #
  # load and merge each deploymet NL values #
  temp=list.files(pattern="*.csv", path='W:\\KJP PHD\\4-Bayesian Habitat Use\\Pdet at Time\\2013')
  
  for (ii in 1:length(temp)){

    if(ii==1){
      NLsheet=read.csv( 
        paste('W:\\KJP PHD\\4-Bayesian Habitat Use\\Pdet at Time\\2013\\', temp[ii],sep=''))
        colnames(NLsheet)[3]='UnitLoc'
        if(NLsheet$UnitLoc[1]=='Arb_10'){
          NLsheet$UnitLoc='Abr_10'
        }
      }else{
        
        NLnew=read.csv( 
          paste('W:\\KJP PHD\\4-Bayesian Habitat Use\\Pdet at Time\\2013\\', temp[ii],sep=''))
        colnames(NLnew)[3]='UnitLoc'
        if(NLsheet$UnitLoc[1]=='Arb_10'){
          NLsheet$UnitLoc='Abr_10'
        }
        NLsheet=rbind(NLsheet, NLnew)
      }
    rm(NLnew)
  }
  

 NLsheet1=NLsheet
 
  # 2014 #
  # load and merge each deploymet NL values #
  temp=list.files(pattern="*.csv", path='W:\\KJP PHD\\4-Bayesian Habitat Use\\Pdet at Time\\2014')
  
  for (ii in 1:length(temp)){
    
    if(ii==1){
      NLsheet=read.csv( 
        paste('W:\\KJP PHD\\4-Bayesian Habitat Use\\Pdet at Time\\2014\\', temp[ii],sep=''))
      colnames(NLsheet)[3]='UnitLoc'
      if(NLsheet$UnitLoc[1]=='Arb_10'){
        NLsheet$UnitLoc='Abr_10'
      }
    }else{
      
      NLnew=read.csv( 
        paste('W:\\KJP PHD\\4-Bayesian Habitat Use\\Pdet at Time\\2014\\', temp[ii],sep=''))
      colnames(NLnew)[3]='UnitLoc'
      if(NLsheet$UnitLoc[1]=='Arb_10'){
        NLsheet$UnitLoc='Abr_10'
      }
      NLsheet=rbind(NLsheet, NLnew)
    }
    rm(NLnew)
  }
  
  NLsheet=rbind(NLsheet, NLsheet1)
  NLsheet$SMCoverage=1
  OccTable=merge(OccTable, NLsheet, by=c('UnitLoc', 'MatlabDate', 'Year'), all.x=T) 
  
#####################################################################################  
  temp=list.files(path = 'W:\\KJP PHD\\3-Detection Function\\Propagation Model\\Area_vs_NL', 
                  pattern = 'csv')
  
  A1=data.matrix(read.csv( paste('W:\\KJP PHD\\3-Detection Function\\Propagation Model\\Area_vs_NL\\', temp[1],sep='')))
  A2=data.matrix(read.csv( paste('W:\\KJP PHD\\3-Detection Function\\Propagation Model\\Area_vs_NL\\', temp[2],sep='')))
  A3=data.matrix(read.csv( paste('W:\\KJP PHD\\3-Detection Function\\Propagation Model\\Area_vs_NL\\', temp[3],sep='')))
  R1=data.matrix(read.csv( paste('W:\\KJP PHD\\3-Detection Function\\Propagation Model\\Area_vs_NL\\', temp[4],sep='')))
  R2=data.matrix(read.csv( paste('W:\\KJP PHD\\3-Detection Function\\Propagation Model\\Area_vs_NL\\', temp[5],sep='')))
  R3=data.matrix(read.csv( paste('W:\\KJP PHD\\3-Detection Function\\Propagation Model\\Area_vs_NL\\', temp[6],sep='')))
  P1=data.matrix(read.csv( paste('W:\\KJP PHD\\3-Detection Function\\Propagation Model\\Area_vs_NL\\', temp[7],sep='')))
  P2=data.matrix(read.csv( paste('W:\\KJP PHD\\3-Detection Function\\Propagation Model\\Area_vs_NL\\', temp[8],sep='')))
  P3=data.matrix(read.csv( paste('W:\\KJP PHD\\3-Detection Function\\Propagation Model\\Area_vs_NL\\', temp[9],sep='')))
  
  
  NLsheet$GroupId=substr(NLsheet$UnitLoc, start = 1, stop = 3)
  NLsheet_NLonly=subset(NLsheet, select=c('MatlabDate', 'MedianNoiseLevel', 'GroupId', 'Year'))
  
  Occ_SMcoverage=subset(OccTable, SMCoverage==1)
  
  # Separate occ table where direct SM coverage and where no SM coverage is available
  Occ_nocoverage=subset(OccTable, is.na(SMCoverage)==TRUE)
  Occ_nocoverage=Occ_nocoverage[, names(Occ_nocoverage) != "MedianNoiseLevel"]
  
  # Merge the noise levels using group id
  Occ_nocoverage=merge(Occ_nocoverage, NLsheet_NLonly, by=c('GroupId', 'MatlabDate', 'Year'), all.x=T)
  
#####################################################################################

  # NLvals=seq(70,150, by=.1)

  Occ_nocoverage$UniLocFactor=as.numeric(factor(Occ_nocoverage$UnitLoc))
  
  idx=as.matrix(which(!is.na(Occ_nocoverage$MedianNoiseLevel)))
  
  row=Occ_nocoverage$UniLocFactor[idx]
  col=(round(Occ_nocoverage$MedianNoiseLevel[idx],1)*10)-700+2
  
  Occ_nocoverage$MedianArea1[idx]=A1[cbind(Occ_nocoverage$UniLocFactor[idx], (round(Occ_nocoverage$MedianNoiseLevel[idx],1)*10)-700+2)]
  Occ_nocoverage$MedianArea2[idx]=A2[cbind(Occ_nocoverage$UniLocFactor[idx], (round(Occ_nocoverage$MedianNoiseLevel[idx],1)*10)-700+2)]
  Occ_nocoverage$MedianArea3[idx]=A3[cbind(Occ_nocoverage$UniLocFactor[idx], (round(Occ_nocoverage$MedianNoiseLevel[idx],1)*10)-700+2)]
  Occ_nocoverage$MedianRange1[idx]=R1[cbind(Occ_nocoverage$UniLocFactor[idx], (round(Occ_nocoverage$MedianNoiseLevel[idx],1)*10)-700+2)]
  Occ_nocoverage$MedianRange2[idx]=R2[cbind(Occ_nocoverage$UniLocFactor[idx], (round(Occ_nocoverage$MedianNoiseLevel[idx],1)*10)-700+2)]
  Occ_nocoverage$MedianRange3[idx]=R3[cbind(Occ_nocoverage$UniLocFactor[idx], (round(Occ_nocoverage$MedianNoiseLevel[idx],1)*10)-700+2)]
  Occ_nocoverage$MedianPdet1[idx]=P1[cbind(Occ_nocoverage$UniLocFactor[idx], (round(Occ_nocoverage$MedianNoiseLevel[idx],1)*10)-700+2)]
  Occ_nocoverage$MedianPdet2[idx]=P2[cbind(Occ_nocoverage$UniLocFactor[idx], (round(Occ_nocoverage$MedianNoiseLevel[idx],1)*10)-700+2)]
  Occ_nocoverage$MedianPdet3[idx]=P3[cbind(Occ_nocoverage$UniLocFactor[idx], (round(Occ_nocoverage$MedianNoiseLevel[idx],1)*10)-700+2)]
  
 OccTable=rbind(Occ_SMcoverage, Occ_nocoverage[,names(Occ_nocoverage) != "UniLocFactor"])

 rm(Occ_nocoverage, Occ_SMcoverage, A1, A2, A3, P1, P2, P3, R1, R2, R3)
 
#####################################################################################
  OccTable$RoundedMatlabDate=floor(OccTable$MatlabDate)
  Trains$RoundedMatlabDate=floor(Trains$MatlabDecDate)
  colnames(Trains)[19]='Hr'

  Bndcod=aggregate(EncounterSpp~UnitLoc+RoundedMatlabDate+Hr, data=subset(Trains, EncounterSpp=='COD/BND'), FUN=length)
  Bndcod$BBOcc=1
  OccTable=merge(OccTable, subset(Bndcod, select=-c(EncounterSpp)), by=c('UnitLoc', 'RoundedMatlabDate', 'Hr'), all.x = T)
  OccTable$BBOcc[is.na(OccTable$BBOcc)]=0
  
  WBDRsd=aggregate(EncounterSpp~UnitLoc+RoundedMatlabDate+Hr, data=subset(Trains, EncounterSpp=='WBD/RSD'), FUN=length)
  WBDRsd$FBOcc=1
  OccTable=merge(OccTable, subset(WBDRsd, select=-c(EncounterSpp)), by=c('UnitLoc', 'RoundedMatlabDate', 'Hr'), all.x = T)
  OccTable$FBOcc[is.na(OccTable$FBOcc)]=0
  
  
  Unk=aggregate(EncounterSpp~UnitLoc+RoundedMatlabDate+Hr, data=subset(Trains, EncounterSpp=='UNK'), FUN=length)
  Unk$UNKOcc=1
  OccTable=merge(OccTable, subset(Unk, select=-c(EncounterSpp)), by=c('UnitLoc', 'RoundedMatlabDate', 'Hr'), all.x = T)
  OccTable$UNKOcc[is.na(OccTable$UNKOcc)]=0
  
  OccTable$OccAll=(OccTable$BBOcc+OccTable$FBOcc+OccTable$UNKOcc)
  OccTable$OccAll=ifelse(OccTable$OccAll>0, OccTable$OccAll/OccTable$OccAll, OccTable$OccAll)
  

  rm(WBDRsd, Unk, Bndcod, Trains2015, Trains2014, Trains2013, Trains)
#######################################################################################     mm= merge.with.order(tempdf, Pdetdf, by=c('UnitLoc','Date',  'Hr', 'Year'),

# Other Covariates #
  OccTable$GroupId=unlist(strsplit(as.character(OccTable$UnitLoc), split = "_"))[seq(1,(nrow(OccTable)*2)-1,2)]
  OccTable$ShoreDist=unlist(strsplit(as.character(OccTable$UnitLoc), split = "_"))[seq(2,(nrow(OccTable)*2),2)]
  
    # Month
   OccTable$Month=as.numeric(do.call(rbind, strsplit(as.character(OccTable$Date), "-"))[,2])
   
   
##########################################################################################
  # Tide Height
   
   TideData=read.csv('W:/KJP PHD/4-Bayesian Habitat Use/R Code/TidalHeight.csv')
   TideData$Phase=as.factor(TideData$Phase)
   TideData$DateFormat=as.Date(TideData$Date,format="%Y-%m-%d")
   TideData$Phase=as.factor(TideData$Phase)
   TideData$DateFormat=as.Date(TideData$Date,format="%Y-%m-%d")
   
   levels(TideData$Phase)=c('Low', 'Flood', 'High', 'Ebb')
   
   temp=TideData[,c('Z','DateFormat','Phase', 'UnitLoc', 'Time', 'Year','HourAfterHigh')]
   colnames(temp)[c(2,5)]=c('Date','Hr')
   temp=temp[!duplicated(temp[,c('Date', 'UnitLoc', 'Hr')]),]
   
   
   OccTable= merge.with.order(OccTable, temp, by=c('UnitLoc','Date',  'Hr', 'Year'),
                              sort = T, keep_order=1)
   
   # # Fix Hour after high tide this was fixed in the PreProcessingTidalData.R file 
   # 19/05/2017 (files on GIT hub)
   # 
   # 
   # 
   # OccTable$yearunitloc=paste(OccTable$Year, OccTable$UnitLoc)
   # OccTable$HourAfterHigh=NA
   # 
   # for(jj in 1:length(unique(OccTable$yearunitloc))){
   #   
   #    OccTable_sub=OccTable[OccTable$yearunitloc==unique(OccTable$yearunitloc)[jj],]
   #    mm=extrema(OccTable_sub$Z)
   #    
   #    OccTable_sub$HourAfterHigh=NA
   #    #OccTable_sub$HourAfterHigh[mm$minindex[,1]]=mm$minindex[,1]
   #    OccTable_sub$HourAfterHigh[mm$maxindex[,1]]=0
   #    OccTable_sub[1:length(seq(-mm$maxindex[1,1],-1))]=seq(-mm$maxindex[1,1],-1)
   #    
   #    kk=1
   #    idx=which(OccTable_sub$HourAfterHigh==0)
   #    
   #    
   #    while(sum(is.na(OccTable_sub$HourAfterHigh))>0 & kk<20){
   #    
   # 
   #      mtry = try(idx[which(is.na(OccTable_sub$HourAfterHigh[idx-kk]))]-kk)
   #      if (!inherits(mtry, "try-error")){
   #        less1_idx=idx[which(is.na(OccTable_sub$HourAfterHigh[idx-kk]))]-kk
   #        OccTable_sub$HourAfterHigh[less1_idx[less1_idx>0 & less1_idx<=nrow(OccTable_sub)]]= -kk
   #      }else{
   #        print('blarg')
   #      }
   #      
   #      mtry = try(idx[which(is.na(OccTable_sub$HourAfterHigh[idx+kk]))]+kk)
   #      if (!inherits(mtry, "try-error")){
   #        plus1_idx=idx[which(is.na(OccTable_sub$HourAfterHigh[idx+kk]))]+kk
   #        OccTable_sub$HourAfterHigh[plus1_idx[plus1_idx>0 & plus1_idx<=nrow(OccTable_sub)]]=  kk
   #      }else{
   #        print('blarg')
   #      }
   #       kk=kk+1
   #       print(kk)
   # 
   #      # less1_idx=idx[which(is.na(OccTable_sub$HourAfterHigh[idx-kk]))]-kk
   #      # plus1_idx=idx[which(is.na(OccTable_sub$HourAfterHigh[idx+kk]))]+kk
   #      # 
   #      # OccTable_sub$HourAfterHigh[less1_idx[less1_idx>0 & less1_idx<nrow(OccTable_sub)]]= -kk
   #      # OccTable_sub$HourAfterHigh[plus1_idx[plus1_idx>0 & plus1_idx<nrow(OccTable_sub)]]=  kk
   #      # 
   #   
   #      
   #    }
   #    
   #    OccTable$HourAfterHigh[OccTable$yearunitloc==unique(OccTable$yearunitloc)[jj]]=OccTable_sub$HourAfterHigh
   #    
   #    rm(OccTable_sub)
   # }

   
   OccTable$yearunitloc=NULL
   
###########################################################################################

   # Add solar elevation
   
   aa=as.data.frame(sunPosition(year = OccTable$Year, month = OccTable$Month,
                                day = as.POSIXlt(OccTable$Date, format="%Y-%m-%d")$mday,
                                hour = OccTable$Hr, min = 0, 
                                sec = 0, lat = OccTable$Lat, long = OccTable$Lon))
   
   OccTable=cbind(OccTable, aa)
   OccTable$HourAfterPeakSolEle=NA
   
   OccTable$dateunitloc=paste(OccTable$Date, OccTable$UnitLoc)
   
   for (ii in 1:length(unique(OccTable$dateunitloc))){
     
     idx=which(OccTable$dateunitloc==unique(OccTable$dateunitloc)[ii])
     
     peaksolar=which.max(OccTable$elevation[idx])
     
     
     OccTable$HourAfterPeakSolEle[idx[peaksolar]]=0
     OccTable$HourAfterPeakSolEle[idx[1]:idx[peaksolar-1]]= seq(-(peaksolar-1),-1)
     OccTable$HourAfterPeakSolEle[idx[peaksolar+1]:max(idx)]= seq(1, length(idx)-peaksolar)
     
   }
   
   OccTable$dateunitloc=NULL
   rm(idx, aa, peaksolar)
   
###########################################################################################

   
####################################################################################
# Export the CSV #
####################################################################################

   write.csv(OccTable,'W:/KJP PHD/4-Bayesian Habitat Use/R Code/OccupancyTable_ThreePdets.csv')   
   
   
   
   