#############################################################################
# Moved to Git Hub 19/05/2017 #
#############################################################################
# This code pre-processes the Tidal height data to trim to C-POD values and 
# add tidal stage

# If/when 2016 data are returned this will need to be re-run with the updated
# meta data

# Update 22/08/2016- Add hour after high tide

rm(list=ls())
library(EMD)
setwd("W:/KJP PHD/4-Bayesian Habitat Use/R Code")


#meta=read.csv('W:/KJP PHD/CPOD Processing/2013 to 2016 CPOD deployments.csv')

meta=read.csv('W:/KJP PHD/Deployment Information/CPODs for Kaitlin.csv')
meta$Usable.from.date=as.Date(meta$Usable.from.date, "%d/%m/%Y", tz = "GMT")
meta$Usable.from.date
meta$Usable.until.date=as.Date(meta$Usable.until.date, "%d/%m/%Y", tz = "GMT")
meta$Usable.Days=difftime(meta$Usable.until.date, meta$Usable.from.date,units='days')
meta$Year=substr(meta$Usable.from.date,1,4)
meta=meta[!is.na(meta$Usable.Days),]

# Remove data for non-recovered units
meta=meta[meta$Usable.from.date !="",]


TidalHeightOut=data.frame()


for(ii in 1:nrow(meta)){
  fname=paste('W:\\KJP PHD\\4-Bayesian Habitat Use\\Covariate Data\\TidalHeight\\', 
              meta$UnitLoc[ii], '.csv', sep="")
  
  TidalDat=read.csv(fname, header = TRUE)

  
  
  TidalDat$Date=as.Date(TidalDat$DateTime_UTC) 
  TidalDat$Time <- as.numeric(substr(as.character(TidalDat$DateTime_UTC),12,13))
  
 
  NewDat=TidalDat
    
  NewDat$Year=as.numeric(format(NewDat$Date, '%Y'))
  NewDat$UnitLoc=meta$UnitLoc[ii]
  
  
  # Calculate the tidal phase
  TideExtrema=extrema(NewDat$Z)
  
  
  #######################
  # Hour After High Tide#
  #######################
  
  
  # Set extrem phases to 0
  NewDat$HourAfterHigh=NA
  NewDat$HourAfterHigh[TideExtrema$maxindex[,1]]=0
  idx=1
  
  while (length(which(is.na(NewDat$HourAfterHigh)))>0 & idx<7){
    
    # Get the indexes of the NA values plus and minus the tidal height hour
    PosTideVals=TideExtrema$maxindex[,1]+idx
    NegTideVals=TideExtrema$maxindex[,1]-idx
    
    PosTideVals=PosTideVals[PosTideVals<nrow(NewDat)]
    NegTideVals=NegTideVals[NegTideVals>0]
    
    
    na_idxsPlus=TideExtrema$maxindex[which(is.na(
      NewDat$HourAfterHigh[PosTideVals]))]+idx
    
    
    na_idxsMinus=TideExtrema$maxindex[which(is.na(
      NewDat$HourAfterHigh[NegTideVals]))]-idx
    
    
    NewDat$HourAfterHigh[na_idxsPlus]=idx
    NewDat$HourAfterHigh[na_idxsMinus[na_idxsMinus>0]]=-idx
    idx=idx+1
    
    #print(idx)
    
    
  }
  
  
 
  ############
  # Phase Analysis
  ############
  
  # Add column for tidal Phase
  NewDat$Phase=0
  
  # Low Tide
  Low=TideExtrema$minindex[,1]
  NewDat$Phase[Low]=1
  textidx=c(Low+1, Low-1)
  NewDat$Phase[textidx[which(NewDat$Z[textidx]<(NewDat$Z[Low+1]-NewDat$Z[Low+1]*.97))]]=1
  NewDat$Phase[textidx[which(NewDat$Z[textidx]<(NewDat$Z[Low+1]-NewDat$Z[Low+1]*.97))]]=1
  
  
  # High Tide
  High=TideExtrema$maxindex[,1]
  NewDat$Phase[High]=3
  textidx=c(High+1, High-1)
  NewDat$Phase[textidx[which(NewDat$Z[textidx]>(NewDat$Z[High+1]-NewDat$Z[High+1]*.97))]]=3
  NewDat$Phase[textidx[which(NewDat$Z[textidx]>(NewDat$Z[High-1]-NewDat$Z[High+1]*.97))]]=3
  
  
  # Flood Tide
  first_low=which(NewDat$Phase==1)[1]
  Flows=seq(from=which(TideExtrema$cross[,1]>first_low)[1], to=length(TideExtrema$cross[,1]), by=2)
  NewDat$Phase[TideExtrema$cross[Flows,]]=2
  
  
  
  NewDat$Phase[TideExtrema$cross[Flows,1]-1]=ifelse( NewDat$Phase[TideExtrema$cross[Flows,1]-1]==0,2,1)
  
  temp=TideExtrema$cross[Flows,2]+1
  NewDat$Phase[temp[temp<nrow(NewDat)]]=ifelse( NewDat$Phase[temp[temp<nrow(NewDat)]]==0,2,3)
  
  # Ebb Tide gets the rest of the values
  NewDat$Phase[NewDat$Phase==0]=4
 
  
 # plot(NewDat$Z, col=NewDat$Phase+1)

  # Trim the data for each year
  IDX=which(NewDat$Date>=as.Date(meta$Usable.from.date[ii], format="%d/%m/%Y") &
              NewDat$Date<=as.Date(meta$Usable.until.date[ii], format="%d/%m/%Y"))

  NewDat=NewDat[IDX,]
  
  # Combine data frames
  TidalHeightOut=rbind(NewDat, TidalHeightOut)
  
  rm(TidalDat)
  
  print((ii))
  }


#########################
# Sort the NA values
#########################

NAIdx=which(is.na(TidalHeightOut$HourAfterHigh))


for (ii in 1:length(NAIdx)){
  
  # Values on either side of the NA idx
  ADJValues=TidalHeightOut$HourAfterHigh[c(NAIdx[ii]-1,NAIdx[ii]+1)]
  
  if (sum(is.na(ADJValues))==0){TidalHeightOut$HourAfterHigh[NAIdx[ii]]=mean(abs(ADJValues))}
  
  
}


##################################################
# Export the Data #
##################################################


write.csv(TidalHeightOut, 'TidalHeight.csv')
