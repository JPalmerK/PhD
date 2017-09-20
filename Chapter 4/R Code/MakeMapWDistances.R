
# 1) Load data and functions ####
rm(list=ls())

library(geosphere)       # for distance between lat long points, estimating distance between river and locs
library(sp)
library(rgeos)
library(mapdata)
library(maptools)
library(ggplot2)
library(raster)
library(marmap)# for bathyemetry
library(landsat) # for calculating the slope
library(elevatr)



# Function to create a circle
SpatialCircle <- function(sp, r=1, npt=100){
  #	Function to generate a circle (SpatialPolygon) with a given radius (m)
  #	sp: SpatialPoint object or matrix of coordinates (in UtM)
  #	r: radius of circle in metres
  #	npt: number of points to draw circle
  if ("Spatial"%in%is(sp)){
    proj <- proj4string(sp)
    coord <- coordinates(sp)
  } else {
    coord <- sp
  }
  x <- coord[,1] + ( r * sin(seq(0,2*pi, length=npt)) )
  y <- coord[,2] + ( r * cos(seq(0,2*pi, length=npt)) )
  poly <- Polygons(srl=list(Polygon(cbind(x,y))), ID=paste("radius_",r,sep=""))
  if ("Spatial"%in%is(sp)){
    return(SpatialPolygons(Srl=list(poly), proj4string=CRS(proj)))
  } else {
    return(poly)
  }
}

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
                            "574038",
                            "562702",
                            "554550",
                            "5741454",
                            "573439"),
                      Lon=c("-022641",
                            "-020337",
                            "-022641",
                            "-030555",
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

# WGS geographic references
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
proj_UTM <- CRS("+proj=utm +zone=30 ellps=WGS84")


# Try higher resolution Scotland Map
gadm <- readRDS("W:/KJP PHD/4-Bayesian Habitat Use/R Code/GBR_adm1.rds") 
gadm=spTransform(gadm, crs.geo)
ScotMap_highRes <- gadm[gadm$NAME_1 == "Scotland", ] #Possibly a ploygone already?



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



# Read in meta data
meta2=read.csv('W:\\KJP PHD\\Deployment Information\\SlopeAndAspect.csv')
meta2$UnitLoc=factor(meta2$UnitLoc, levels=level_names)
meta2$DistToSalmonRun=0
meta2$RiverName='blarg'

# Make fake meta2 for recomendations
meta2_temp=data.frame(Lat=numeric(length = 8),
                      Lon=numeric(length = 8))
meta2_temp$UnitLoc=as.factor(c('New_05', 'Bor_05','Mac_05',
                                                    'Pet_05', 'Bla_05', 'Ivb_05',
                                                    'Leu_05', 'Man_05'))

meta2_temp$UnitLoc[1]='New_05' #Newport
meta2_temp$Lat[1]=58.19
meta2_temp$Lon[1]=-3.448

meta2_temp$UnitLoc[2]='Bor_05'#Brora
meta2_temp$Lat[2]=  57.933
meta2_temp$Lon[2]= -3.66

meta2_temp$UnitLoc[3]='Mac_05'#Macduff
meta2_temp$Lat[3]=  57.70
meta2_temp$Lon[3]= -2.51

meta2_temp$UnitLoc[4]='Pet_05'#Peterhead
meta2_temp$Lat[4]= 57.59
meta2_temp$Lon[4]= -1.68

meta2_temp$UnitLoc[5]='Bla_05'#BlackDog
meta2_temp$Lat[5]=  57.234
meta2_temp$Lon[5]=  -1.962

meta2_temp$UnitLoc[6]='Ivb_05'#Inverbervie
meta2_temp$Lat[6]=   56.764
meta2_temp$Lon[6]=   -2.28


meta2_temp$UnitLoc[7]='Leu_05'#Leuchars
meta2_temp$Lat[7]=    56.3933
meta2_temp$Lon[7]=    -2.73774


meta2_temp$UnitLoc[8]='Man_05'#Mann
meta2_temp$Lat[8]=    56.086
meta2_temp$Lon[8]=     -2.44920


#river_coords=river_locs[, c('lonDeg','LatDeg' )]
coordinates(river_locs)=c('lonDeg','LatDeg' )


proj4string(river_locs) <- crs.geo  # define projection system of our data
summary(river_locs)

coordinates(meta2)=meta2[,c('Lon','Lat' )]
proj4string(meta2) <- crs.geo  # define projection system of our data

coordinates(meta2_temp)=meta2_temp[,c('Lon', 'Lat')]
proj4string(meta2_temp) <- crs.geo  # define projection system of our data


# Load bathymetry (MarMap)
getNOAA.bathy(lon1 = min(meta2$Lon)-1, lon2 =  max(meta2$Lon)+.2,
              lat1 =min(meta2$Lat)-.5, lat2 = max(meta2$Lat)+.5,
              resolution = 1)->NorthSea

# Plot the bathymetry
blues <- c("lightsteelblue4", "lightsteelblue3",
           "lightsteelblue2", "lightsteelblue1")

greys <- c(grey(0.6), grey(0.93), grey(0.99))

plot(NorthSea, n = 0, lwd = 0.5, image=TRUE, 
     bpal = list(c(0, 10, grey(.7), grey(.9), grey(.95)),
                 c(min(NorthSea), 1, "darkblue", "lightblue")))

points(river_locs, pch = 21, col = "black",
       bg = "yellow", cex = 1.3)
points(meta2, pch = 21, col = "black",
       bg = "red", cex = 1.3)

points(meta2_temp, pch = 21, fill = "black",
       bg = "green", cex = 1.3)



# Create a raster version
NorthSea_raster=marmap::as.raster(NorthSea)
NorthSea_raster=projectRaster(NorthSea_raster, crs = crs.geo)
NorSea_spatialGrid=as.SpatialGridDataFrame(bathy = NorthSea)


# Calculate slope values
Slope=terrain(NorthSea_raster, opt = 'Slope', unit = 'Radians', neighbors = 4)


# Distance to Shore
meta2$DistToShore=dist2isobath(NorthSea, coordinates(meta2), isobath = 0)$distance
meta2_temp$DistToShore=dist2isobath(NorthSea, coordinates(meta2_temp), isobath = 0)$distance




# calculate distance to nearest salmon river
for(ii in 1:30){
  temp=rep(0, length(river_locs))
  
  for(jj in 1:length(river_locs)){
    temp[jj]=distm (coordinates(meta2)[ii,], 
                    coordinates(river_locs)[jj,], 
                    fun = distHaversine)
  }
  
  
  meta2$DistToSalmonRun[ii]=min(temp)
  meta2$RiverName[ii]=as.character(river_locs$Rivername[which.min(temp)])
  
  
  if(ii<9){
    for(jj in 1:length(river_locs)){
      temp[jj]=distm (coordinates(meta2_temp)[ii,], 
                      coordinates(river_locs)[jj,], 
                      fun = distHaversine)
    }
    
    
    meta2_temp$DistToSalmonRun[ii]=min(temp)
    meta2_temp$RiverName[ii]=as.character(river_locs$Rivername[which.min(temp)])
    
  }
  
  
}



# 3) Create Data Grid (WGS) and calculate distance to rivers #######


# Create datapoint grids
dat=expand.grid(Lon=seq(min(meta2$Lon)-.2,
                        max(meta2$Lon)+.15, length.out = 250),
                Lat=seq(min(meta2$Lat-.3),
                        max(meta2$Lat+.25), length.out = 250))

coordinates(dat)=c( 'Lon','Lat')
proj4string(dat) <-crs.geo

## Plot checking
# dev.off()
plot(NorthSea, n = 0, lwd = 0.5, image=TRUE,
     bpal = list(c(0, 10, grey(.7), grey(.9), grey(.95)),
                 c(min(NorthSea), 1, "darkblue", "lightblue")))
# points(dat, cex=.5, bg='red')

# Get depth
dat$Depth=get.depth(NorthSea, cbind(dat$Lon, dat$Lat), locator = FALSE)
meta2_temp$Depth_m=get.depth(NorthSea,x = cbind(meta2_temp$Lon, meta2_temp$Lat), locator=FALSE)[,3]

# Filter datapoints where depth is less than 1 meter of water
dat <- dat[dat@data$Depth.depth < 1,]

# Get distane to nearest shore
dat$DistToShore= dist2isobath(NorthSea, dat$Depth.lon , dat$Depth.lat, isobath=0, locator=FALSE)
  
# Get slope values
dat$Slope=raster::extract(Slope, dat)
meta2$SlopeMap=raster::extract(Slope, meta2)
meta2_temp$SlopeMap=raster::extract(Slope, meta2_temp)



# Calculate distance between features and all data
for(ii in 1:length(unique(meta2$RiverName))){
  
  rname=as.character(unique(meta2$RiverName)[ii])
  
  dat@data[, rname]= 
    distRhumb(river_locs[river_locs$Rivername==as.character(unique(meta2$RiverName)[ii]),], dat)
}

# Calculate minimum distances 
dat$DistToSalmonRun=apply(as.data.frame(dat@data[,6:11]), 1, FUN=function(x){min(x, na.rm = TRUE)})
dat$RiverName=apply(as.data.frame(dat@data[,6:11]), 1, 
                    FUN=function(x){colnames(dat@data)[6:11][which.min(x)]})
  
  
  # Create River Database for localisation
  RiverRanges=data.frame(aggregate(data=meta2, DistToSalmonRun~RiverName, FUN=range))
  RiverRanges=merge(RiverRanges, river_locs, 
                    all.x=TRUE, by.x=c('RiverName'), by.y='Rivername')
  
  # Convert to UTM 
  coordinates(RiverRanges)=c('lonDeg','LatDeg' )
  proj4string(RiverRanges) <- crs.geo  # define projection system of our data
  summary(RiverRanges)

  
  # Determine datapoints within feature ranges of the monitoring area
  featureRangeIDX=numeric()

for(ii in 1:length(unique(meta2$RiverName))){
  

  
  ## Create two circles, one for the minimum range one for the max, will interpolate
  r1= RiverRanges$DistToSalmonRun[ii,1]-1000 #2km
  r2= RiverRanges$DistToSalmonRun[ii,2]+1000 #2km
  
  
  river_locs=RiverRanges[ii,]
  river_locs_UTM=spTransform(river_locs, proj_UTM)
  
  
  ## Circles
  circle1 <- SpatialCircle(river_locs_UTM, r1)
  circle2 <- SpatialCircle(river_locs_UTM, r2)
  # convert to wgs
  circle1_wgs=spTransform(circle1, crs.geo)
  circle2_wgs=spTransform(circle2, crs.geo)
  
  
  plot(circle2_wgs, add=TRUE)
  plot(circle1_wgs, add=TRUE)
  
  
  inside.feature= which(!is.na(over(dat, as(circle2_wgs, "SpatialPolygons"))) &
                        is.na(over(dat, as(circle1_wgs, "SpatialPolygons"))))
  
  featureRangeIDX=c(featureRangeIDX, inside.feature)
  
  
  # 
  # # Plot checking
  # points(dat[inside.river,], pch=16, col="red") #good
  # points(meta2, pch=16, col="green") #good
  # points(RiverRanges, pch=16, col="blue") #good
  # 
  #dat[-inside.feature, RiverRanges$RiverName[ii]]=NA

  rm(circle1_wgs, circle2_wgs, circle1, circle2, r1, r2 )
}




# Remove Duplicates
featureRangeIDX=unique(featureRangeIDX)

points(dat[!is.na(dat$`Cromarty Firth`),], pch = 21, col = "Green",
       bg = "green", cex = 1.3)


# Make Spatial Models ######################################################################################

# Clear out some crap 

rm(list=setdiff(ls(), c("dat", 'NorthSea', 'NorthSea_raster',
                        'proj_UTM', 'crs.geo', 'meta2', 'meta2_temp','level_names')))


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
library(geosphere)  
library(fields)
library(marmap)


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




# colorblind palette with black for plotting
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




# 2)  Data Prep #################################################################################

# This section processes the hourly occupancy table to correct level orders (GroupID and UnitLoc)
# Hourly Bottlenose dolphin offsets are delcared where offset represnts probability that a click detection
# (OccAll) was produced by a broadband species. If only frequency banded clicks, P(BND)=0.06 if Broadband
# clicks P(BND)=0.77 and if unknown clicks P(BND)=0.50


OccTable= read.csv('W:/KJP PHD/4-Bayesian Habitat Use/R Code/OccupancyTable_ThreePdets.csv')
meta=read.csv('W:/KJP PHD/Deployment Information/CPODs for Kaitlin.csv')

meta2=merge(meta2, distinct(meta[, c('UnitLoc', 'Depth_m')]), by='UnitLoc',all.x=TRUE)
meta2$UnitLoc=factor(meta2$UnitLoc, levels=level_names)



OccTable$UnitLoc=factor(OccTable$UnitLoc, levels=level_names)
OccTable$DateUnitloc=as.factor(paste(OccTable$Date, OccTable$UnitLoc))

# Distance from shore factor
OccTable$ShoreDist=as.character(OccTable$ShoreDist)
OccTable$ShoreDist[OccTable$ShoreDist=="5"]='05'
OccTable$ShoreDist=as.factor(OccTable$ShoreDist)

# Group ID factor
level_names_group=c( "Lat", "Hel", "Cro",
               "SpB", "Fra", "Cru",
               "Sto", "Abr", "StA",
               "Stb")
OccTable$GroupId=unlist(strsplit(as.character(OccTable$UnitLoc), split = "_"))[seq(1,(nrow(OccTable)*2)-1,2)]
OccTable$GroupId=factor(OccTable$GroupId, levels=level_names_group)


# YEAR
OccTable$Year=as.factor(OccTable$Year)

# Merge meta2 and Occupancy Table
OccTable=merge(OccTable, meta2[,c('UnitLoc', 'DistToSalmonRun', 'RiverName', 'DistToShore', 'SlopeMap', 'Depth_m')],
               by='UnitLoc', all.x = TRUE)

OccTable$RiverName=as.factor(OccTable$RiverName)


OccTable$JD_scale=scale(OccTable$JulienDay)

OccTable$FBOcc[is.na(OccTable$FBOcc)]=0
OccTable$BBOcc[is.na(OccTable$BBOcc)]=0
OccTable$UNKOcc[is.na(OccTable$UNKOcc)]=0


# 2a) Aggregate hourly Detections into Daily Occupancy table ######################################################

# This section of code aggregates the hourly occupancy into daily occupancy and
# Assigns BND offset scores to each day indicating the probability that on any given day,
# There was a bottlenose dolphin in the detections

# Offset scores where multiple species present calculated using a weighted average of the hourly 
# detection scores.(Next iteration migh consider likelihood ratios instead, but enough stuff to deal with at present)


# Add total number of detections 
mm.bbtot=as.data.frame(aggregate(BBOcc~UnitLoc+Date, FUN=sum, data = OccTable))
colnames(mm.bbtot)[3]='BBTot'
mm.fbtot=as.data.frame(aggregate(FBOcc~UnitLoc+Date, FUN=sum, data = OccTable))
colnames(mm.fbtot)[3]='FBTot'
mm.unktot=as.data.frame(aggregate(UNKOcc~UnitLoc+Date, FUN=sum, data = OccTable))
colnames(mm.unktot)[3]='UNKTot'


mm=distinct(OccTable, Date, UnitLoc, JulienDay, GroupId, ShoreDist, SlopeMap, 
            Year, Month, RiverName, DistToSalmonRun, Depth_m, DistToShore)
mm.bb=distinct(subset(OccTable, BBOcc>0), Date, BBOcc, UnitLoc)
mm.fb=distinct(subset(OccTable, FBOcc>0), Date, FBOcc, UnitLoc)
mm.unk=distinct(subset(OccTable, UNKOcc>0), Date, UNKOcc, UnitLoc)


OccTable_daily=merge(mm, mm.bb, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.fb, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.unk, by = c('Date', 'UnitLoc'), all.x = TRUE)

OccTable_daily=merge(OccTable_daily, mm.bbtot, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.fbtot, by = c('Date', 'UnitLoc'), all.x = TRUE)
OccTable_daily=merge(OccTable_daily, mm.unktot, by = c('Date', 'UnitLoc'), all.x = TRUE)


# Set NA scores to 0
OccTable_daily[is.na(OccTable_daily)] <- 0


OccTable_daily$SpeciesOffset=OccTable_daily$BBOcc+OccTable_daily$FBOcc+OccTable_daily$UNKOcc
OccTable_daily$OccAll=ifelse(OccTable_daily$SpeciesOffset>=1,1,0)

# Species offset for Daily Occupancy
# If two species same day and same unit then uncertain
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset>1]=0.5
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==1 & OccTable_daily$BBOcc==1]=0.77
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==1 & OccTable_daily$FBOcc==1]=0.06
OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==1 & OccTable_daily$UNKOcc==1]=0.5
#OccTable_daily$SpeciesOffset[OccTable_daily$SpeciesOffset==0] = .01



# Total offset
OccTable_daily$BNDTotOffset=(OccTable_daily$BBTot*.77+OccTable_daily$FBTot*.06+OccTable_daily$UNKTot*.5)/
  (OccTable_daily$BBTot+OccTable_daily$FBTot+OccTable_daily$UNKTot)
OccTable_daily$TotDet=(OccTable_daily$BBTot+OccTable_daily$FBTot+OccTable_daily$UNKTot)
OccTable_daily$BNDTotOffset[is.na(OccTable_daily$BNDTotOffset)]=0
OccTable_daily$dateunit=as.factor(paste(OccTable_daily$Date, OccTable_daily$UnitLoc))
OccTable_daily$ScaleJd=scale(OccTable_daily$JulienDay)


# Add season
OccTable_daily$Season=NA
OccTable_daily$Season[OccTable_daily$Month < 6]='Spring'
OccTable_daily$Season[OccTable_daily$Month > 8]='Autum'
OccTable_daily$Season[which(is.na(OccTable_daily$Season))]='Summer'
OccTable_daily$Season=as.factor(OccTable_daily$Season)

OccTable_daily$CentredMonth=OccTable_daily$Month-median(median(unique(OccTable_daily$Month)))
OccTable_daily$yearseason=as.factor(paste(OccTable_daily$Year, OccTable_daily$Season))
OccTable_daily$SeasonRiver=as.factor(paste(OccTable_daily$Season, OccTable_daily$RiverName))

rm(mm, mm.bb, mm.fb, mm.unk, mm.bbtot, mm.fbtot, mm.unktot)

##### Spatial Models##########################################################################

# Investagate correlation between covariates

cor(meta2$Depth_m, meta2$DistToSalmonRun, method='pearson') # correlation -0.4
cor(meta2$Depth_m, meta2$SlopeMap, method='pearson') 
cor(meta2$Depth_m, meta2$DistToShore, method='pearson') 

cor(meta2$DistToSalmonRun, meta2$SlopeMap, method='pearson') 
cor(meta2$DistToSalmonRun, meta2$DistToShore, method='pearson') 

cor(meta2$SlopeMap, meta2$DistToShore, method='pearson') 

# All highly fucking correlated! 


# Use PCA make new composite variable
pcadata=as.data.frame(meta2[,c('DistToSalmonRun',  'SlopeMap', 'DistToShore', 'Depth_m')])[,1:4]
# scale PCA data
pcadata=apply(pcadata, 2, scale)
PCAmod=princomp(pcadata)
summary(PCAmod)
#93% of the variation accounted for by the first component of the PCA


# Create a new column for the dat representing the conglomeration of the variables
newdata=data.frame(dat@data)
colnames(newdata)[5]='SlopeMap'
colnames(newdata)[3]='Depth_m'

# Create new variable for the map and for the Occupancy table
newdaa=newdata[,c('DistToSalmonRun',  'SlopeMap', 'DistToShore', 'Depth_m')]
dat$PCAcomp=predict(PCAmod, newdata=newdaa)[,3]

newdaa=OccTable_daily[,c('DistToSalmonRun',  'SlopeMap', 'DistToShore', 'Depth_m')]
OccTable_daily$PCAcomp=predict(PCAmod, newdata=newdaa)[,1]



# New variable for the model
meta2$PCA_predict=predict(PCAmod)[,1]



# Spatial Model
modlist_spatial= gamm(BNDTotOffset~s(PCAcomp, bs='ts'),
                           correlation=corAR1(form = ~1|dateunit),
                           family=binomial,
                           data=OccTable_daily, fx=TRUE ) #54655.14

modlist_spatial_season= gamm(BNDTotOffset~s(PCAcomp, k = 3, bs='ts', by=Season)+Season,
                      correlation=corAR1(form = ~1|dateunit),
                      family=binomial,
                      data=OccTable_daily, fx=TRUE ) #55375.22

modlist_spatial_jdate= gamm(BNDTotOffset~s(PCAcomp, bs='ts', k=3)+s(JulienDay, k=2),
                             correlation=corAR1(form = ~1|dateunit),
                             family=binomial,
                             data=OccTable_daily, fx=TRUE ) # 

modlist_spatial_jdate_te= gamm(BNDTotOffset~te(PCAcomp,JulienDay, bs='ts', k=4),
                            correlation=corAR1(form = ~1|dateunit),
                            family=binomial,
                            data=OccTable_daily, fx=TRUE ) #

modlist_spatial_jdate_te= gamm(BNDTotOffset~te(PCAcomp,JulienDay, bs='ts', k=4),
                               correlation=corAR1(form = ~1|dateunit),
                               family=binomial,
                               data=OccTable_daily, fx=TRUE ) # 


modlist_spatial_jdate_te_wrandom= gamm(BNDTotOffset~te(PCAcomp,JulienDay, bs='ts', k=4),
                               correlation=corAR1(form = ~1|dateunit),
                               family=binomial,
                               random=list(UnitLoc=~1),
                               data=OccTable_daily, fx=TRUE ) #



modlist_spatial_jdate_te_nocro= gamm(BNDTotOffset~te(PCAcomp,JulienDay, bs='ts', k=4),
                               correlation=corAR1(form = ~1|dateunit),
                               family=binomial,
                               data=subset(OccTable_daily, UnitLoc !='Cro_05'), fx=TRUE ) #54012.59 


modlist_spatial_jdate_rivername= gamm(BNDTotOffset~s(PCAcomp, k=3, by=RiverName) +RiverName ,
                                      correlation=corAR1(form = ~1|dateunit),
                                      family=binomial,
                                      data=OccTable_daily, fx=TRUE ) # 55302.02

OccTable_daily_nocro=subset(OccTable_daily, UnitLoc != 'Cro_05')= gamm(BNDTotOffset~s(PCAcomp, k=3, by=RiverName) +RiverName ,
                               correlation=corAR1(form = ~1|dateunit),
                               family=binomial,
                               data=OccTable_daily, fx=TRUE ) # 55302.02

OccTable_daily_nocro=subset(OccTable_daily, UnitLoc != 'Cro_05')


# Ideal Model
ModelFull= gamm(BNDTotOffset~s(JulienDay, k=3, by=RiverName) +RiverName 
                + s(DistToSalmonRun, bs='ts', by=RiverName)+s(Depth_m, bs='ts') +
                  s(DistToShore, bs='ts')+
                  s(SlopeMap, bs='ts')+Year,
                random=list(UnitLoc=~1),
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily, fx=TRUE ) 

# Remove one interaction
ModelFull= gamm(BNDTotOffset~s(JulienDay, k=3, by=RiverName) +RiverName 
                + s(DistToSalmonRun, bs='ts')+s(Depth_m, bs='ts') +
                  s(DistToShore, bs='ts')+
                  s(SlopeMap, bs='ts')+Year,
                random=list(UnitLoc=~1),
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily, fx=TRUE ) 

# Remove other interaction
ModelFull= gamm(BNDTotOffset~ #s(JulienDay, k=3, bs='ts') +
                  s(SlopeMap, bs='ts')+
                  s(DistToSalmonRun, bs='ts', k=3) +
                  s(Depth_m, bs='ts') +
                  s(DistToShore, bs='ts') +
                  RiverName + Season,
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily, Select=TRUE ) 

FinalModel=gamm(BNDTotOffset~ #s(JulienDay, k=3, bs='ts') +
                  #s(SlopeMap, bs='ts')+
                  s(DistToSalmonRun, bs='ts', k=3) +
                  s(Depth_m, bs='ts') +
                  s(DistToShore, bs='ts') +
                  RiverName + Season,
                correlation=corAR1(form = ~1|dateunit),
                family=binomial,
                data=OccTable_daily, Select=TRUE ) 



# Check for co

# Run VIF score
vifmodel=glm(BNDTotOffset~SlopeMap + Depth_m + DistToSalmonRun +DistToShore,
             family=binomial,
             data=OccTable_daily)


library(car)

vif(vifmodel)

# vif suggests colinearity not as big of an issue as thought



# 
# modlist_spatial=list()
# 
# 
# modlist_spatial[[1]]= gamm(BNDTotOffset~s(Depth_m, bs='ts', k=3),
#                            correlation=corAR1(form = ~1|dateunit),
#                            family=binomial,
#                            data=OccTable_daily)
# 
# modlist_spatial[[2]]= gamm(BNDTotOffset~s(SlopeMap, bs='ts', k=3),
#                            correlation=corAR1(form = ~1|dateunit),
#                            family=binomial,
#                            data=OccTable_daily)
# 
# modlist_spatial[[3]]= gamm(BNDTotOffset~s(DistToShore, bs='ts', k=3),
#                            correlation=corAR1(form = ~1|dateunit),
#                            family=binomial,
#                            data=OccTable_daily)
# 
# modlist_spatial[[4]]= gamm(BNDTotOffset~s(DistToSalmonRun, bs='ts', k=3),
#                            correlation=corAR1(form = ~1|dateunit),
#                            family=binomial,
#                            data=OccTable_daily)
# 
# 
# modlist_spatial[[5]]= gamm(BNDTotOffset~s(DistToSalmonRun, bs='ts', k=3)+s(Depth_m, bs='ts', k=3),
#                            correlation=corAR1(form = ~1|dateunit),
#                            family=binomial,
#                            data=OccTable_daily)
# 
# 
# modlist_spatial[[6]]= gamm(BNDTotOffset~s(Depth_m, bs='ts', k=3)+s(SlopeMap, bs='ts', k=3),
#                            correlation=corAR1(form = ~1|dateunit),
#                            family=binomial,
#                            data=OccTable_daily)
# 
# 
# modlist_spatial[[7]]= gamm(BNDTotOffset~s(SlopeMap, bs='ts', k=3, by=RiverName),
#                            correlation=corAR1(form = ~1|dateunit),
#                            family=binomial,
#                            data=OccTable_daily)
# 
# 
# modlist_spatial[[8]]= gamm(BNDTotOffset~s(DistToSalmonRun, bs='ts', k=3)+s(Depth_m, bs='ts', k=3) + Season,
#                            correlation=corAR1(form = ~1|dateunit),
#                            family=binomial,
#                            data=OccTable_daily)
# 
# modlist_spatial[[9]]= gamm(BNDTotOffset~s(DistToSalmonRun, bs='ts', k=3) + 
#                              s(SlopeMap, bs='ts', k=3) + Season,
#                            correlation=corAR1(form = ~1|dateunit),
#                            family=binomial,
#                            data=OccTable_daily,
#                            random=list(UnitLoc=~1))
# 
# 
# modlist_spatial[[10]]= gamm(BNDTotOffset~s(DistToSalmonRun, bs='ts', k=3, by=Season)
#                             +s(SlopeMap, bs='ts', k=3, by=Season) + Season,
#                            correlation=corAR1(form = ~1|dateunit),
#                            family=binomial,
#                            data=OccTable_daily,
#                            random=list(UnitLoc=~1))
# 
# modlist_spatial[[11]]= gamm(BNDTotOffset~s(DistToSalmonRun, bs='ts', k=3)+
#                               s(SlopeMap, bs='ts', k=3, by=Season) + Season,
#                             correlation=corAR1(form = ~1|dateunit),
#                             family=binomial,
#                             data=OccTable_daily,
#                             random=list(UnitLoc=~1))
# 
# modlist_spatial[[12]]= gamm(BNDTotOffset~s(SlopeMap, bs='ts', k=3, by=Season) + Season,
#                             correlation=corAR1(form = ~1|dateunit),
#                             family=binomial,
#                             data=OccTable_daily,
#                             random=list(UnitLoc=~1))
# 
# modlist_spatial[[13]]= gamm(BNDTotOffset~s(DistToSalmonRun, bs='ts', k=3, by=Season)
#                             +s(DistToShore, bs='ts', k=3, by=Season) + Season,
#                             correlation=corAR1(form = ~1|dateunit),
#                             family=binomial,
#                             data=OccTable_daily,
#                             random=list(UnitLoc=~1))
# 
# modlist_spatial[[14]]= gamm(BNDTotOffset~s(SlopeMap, bs='ts', k=3)
#                             +s(DistToShore, bs='ts', k=3, by=Season) + Season,
#                             correlation=corAR1(form = ~1|dateunit),
#                             family=binomial,
#                             data=OccTable_daily,
#                             random=list(UnitLoc=~1))
# 
# 
# modlist_spatial[[15]]= gamm(BNDTotOffset~s(SlopeMap, bs='ts', k=3, by=Season)
#                             + s(DistToShore, bs='ts', k=3) + Season,
#                             correlation=corAR1(form = ~1|dateunit),
#                             family=binomial,
#                             data=OccTable_daily,
#                             random=list(UnitLoc=~1))
# 
# modlist_spatial[[16]]= gamm(BNDTotOffset~s(SlopeMap, bs='ts', k=3, by=Season)
#                             + s(DistToShore, bs='ts', k=3) + Season,
#                             correlation=corAR1(form = ~1|UnitLoc),
#                             family=binomial,
#                             data=OccTable_daily,
#                             random=list(UnitLoc=~1))


lapply(modlist_spatial, AIC)

jdayval=c(150, 212, 275)
#just use model 10

op<-par(no.readonly=TRUE)
for(ii in 1:3){
  
  # Filter by distance
  dat1=dat[dat$DistToShore<(max(meta2$DistToShore)+2000),]
  Preddat=data.frame(SlopeMap=dat1$Slope, 
                     Depth_m=dat1$Depth.depth,
                     DistToSalmonRun=dat1$DistToSalmonRun,
                     BNDTotOffset=rep(0, nrow(dat1)),
                     DistToShore=dat1$DistToShore,
                     Season= unique(OccTable_daily$Season)[ii],
                     JulienDay=jdayval[ii],
                     RiverName=dat1$RiverName,
                     Year=2013)

  # Pca combined model
  Preddat$PCAcomp=predict(PCAmod, Preddat[,c('DistToSalmonRun',  'SlopeMap', 'DistToShore', 'Depth_m')])[,1]
  
  
  
  preds=predict(FinalModel, Preddat, se.fit=TRUE)
  #preds=predict(modlist_spatial[[10]], Preddat, se.fit=TRUE)
  preds$UCI=preds$fit+(1.96*preds$se.fit)
  preds$LCI=preds$fit-(1.96*preds$se.fit)
  
  # Put on the binary scale
  preds[]<-lapply(preds, inv.logit)

  dat1=cbind(dat1, preds)

  rbPal <- colorRampPalette(c('red','blue'))
  #dat1$Col <- rbPal(10)[as.numeric(cut(dat1$fit,breaks = 20))]

  dat1$Col <- heat.colors(20)[as.numeric(cut(dat1$fit,
                                             breaks = as.numeric(quantile(unlist(preds), seq(0, 1, length.out = 20)))))]

  dat1$ColLCI <- heat.colors(20)[as.numeric(cut(dat1$LCI,
                                             breaks = as.numeric(quantile(unlist(preds), seq(0, 1, length.out = 20)))))]
  
  dat1$ColUCI <- heat.colors(20)[as.numeric(cut(dat1$UCI,
                                                breaks = as.numeric(quantile(unlist(preds), seq(0, 1, length.out = 20)))))]
  
  dat1$Season=as.character(Preddat$Season)

  # Plot the bathymetry
  blues <- c("lightsteelblue4", "lightsteelblue3",
           "lightsteelblue2", "lightsteelblue1")

  greys <- c(grey(0.6), grey(0.93), grey(0.99))

   png(filename = paste(Preddat$Season[1], '.png'),
       units="in", 
       width=7, 
       height=9, 
       pointsize=12,res = 400)
  

  # par(oma=c(0,0,0,0),mar=c(3,3,5.2,0),mfrow=c(2,2),pch=16)
   par(op)
   par(oma=c(1,2,4,1),mar=c(3,3,7,0),mfrow=c(2,2),pch=16)

   #  fit
  plot(NorthSea, n = 0, lwd = 0.5, image=TRUE, 
       bpal = list(c(0, 10, grey(.7), grey(.9), grey(.95)),
                   c(min(NorthSea), 1, "lightsteelblue3", "lightsteelblue1")),
       main= 'Fit', 
       xlim=range(coordinates(NorthSea_raster)[,1]),
       ylim=c(56 ,58),  xaxs="i", yaxs="i", frame.plot = FALSE)
  
  scaleBathy(NorthSea, deg=1, x="bottomleft", y=NULL, inset=10, angle=90)
  
  points(x = dat1$Depth.lon, y=dat1$Depth.lat,
       pch = 15, cex=.5,
       col = dat1$Col, main=Preddat$Season[1])
  points(river_locs, pch=19, col='blue')
  
  #image.plot(legend.only = TRUE,zlim=range(unlist(preds)), col = heat.colors(20),)

  points(meta2,
         pch = 20, cex=.5,
         col = 'black')
  #text(river_locs, river_locs$Rivername)

  #  LCI
  plot(NorthSea, n = 0, lwd = 0.5, image=TRUE, 
       bpal = list(c(0, 10, grey(.7), grey(.9), grey(.95)),
                   c(min(NorthSea), 1, "lightsteelblue3", "lightsteelblue1")),
       main= 'Lower 95% CI', 
       xlim=range(coordinates(NorthSea_raster)[,1]),
       ylim=c(56 ,58),  xaxs="i", yaxs="i", frame.plot = FALSE)
  
  scaleBathy(NorthSea, deg=1, x="bottomleft", y=NULL, inset=10, angle=90)
  
  points(x = dat1$Depth.lon, y=dat1$Depth.lat,
         pch = 15, cex=.18,
         col = dat1$ColLCI, main=Preddat$Season[1])
  
  points(meta2,
         pch = 20, cex=.5,
         col = 'black')
  points(river_locs, pch=19, col='blue')
  
  # UCI
  plot(NorthSea, n = 0, lwd = 0.5, image=TRUE, 
       bpal = list(c(0, 10, grey(.7), grey(.9), grey(.95)),
                   c(min(NorthSea), 1, "lightsteelblue3", "lightsteelblue1")),
       main= 'Upper 95% CI',
       xlim=range(coordinates(NorthSea_raster)[,1]),
       ylim=c(56 ,58),  xaxs="i", yaxs="i", frame.plot = FALSE)
  
  scaleBathy(NorthSea, deg=1, x="bottomleft", y=NULL, inset=10, angle=90)
  
  points(x = dat1$Depth.lon, y=dat1$Depth.lat,
         pch = 15, cex=.18,
         col = dat1$ColUCI, main=Preddat$Season[1])
  
  points(meta2,
         pch = 20, cex=.5,
         col = 'black')
  points(river_locs, pch=19, col='blue')
  
  # Plot titel if using Season
  # mtext(text=Preddat$Season[1],side=3,line=0,outer=TRUE, cex=2)
  
  # Plto title otherwise
  mtext(text= 'Daily Occupancy Probability',side=3,line=0,outer=TRUE, cex=2)
  mtext(text="Longitude",side=1,line=0,outer=TRUE)
  mtext(text="Latitude",side=2,line=0,outer=TRUE)
  
  image.plot(legend.only = TRUE,zlim=range(unlist(preds)), 
             col = heat.colors(20), legend.mar = 0, legend.shrink = .5, 
             legend.width = 1)
 
  
  
  
   dev.off()
  
  
  }



# Make Partial models
PmodelFull= gamm(BNDTotOffset~s(SlopeMap, k=3, bs='ts') + s(Depth_m, bs='ts', k=3) +
                   s(DistToSalmonRun,bs='ts', k=3) + s(DistToShore, bs='ts', k=3),
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily, niterPQL = 400) #

PmodelSlope= gamm(BNDTotOffset~s(SlopeMap, k=3, bs='ts') ,
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily, niterPQL = 400) #

PmodelDepth= gamm(BNDTotOffset~ s(Depth_m, bs='ts', k=3),
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily, niterPQL = 400) #

PmodelDisttoSalmon= gamm(BNDTotOffset~s(DistToSalmonRun,bs='ts', k=3),
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily, niterPQL = 400) #

PmodelDisttoShore= gamm(BNDTotOffset~ s(DistToShore, bs='ts', k=3),
                 correlation=corAR1(form = ~1|dateunit),
                 family=binomial,
                 data=OccTable_daily, niterPQL = 400) #






par(mfrow=c(2,3))
# Model output
png(filename = paste('Spatial Model1.png'),
    units="in", 
    width=9, 
    height=7, 
    pointsize=12,res = 400)


plot(modlist_spatial[[10]]$gam, trans=inv.logit)
dev.off()



# Summer UCI
plot(NorthSea, n = 0, lwd = 0.5, image=TRUE, 
     bpal = list(c(0, 10, grey(.7), grey(.9), grey(.95)),
                 c(min(NorthSea), 1, "darkblue", "lightblue")),
     main='Summer')


points(x = dat1$Depth.lon[dat1$Season=='Summer'], y=dat1$Depth.lat[dat1$Season=='Summer'],
       pch = 20,
       col = dat1$Col, main='Summer')

points(meta2,
       pch = 18,
       col = 'black')



# Summer
plot(NorthSea, n = 0, lwd = 0.5, image=TRUE, 
     bpal = list(c(0, 10, grey(.7), grey(.9), grey(.95)),
                 c(min(NorthSea), 1, "darkblue", "lightblue")),
     main='Autum')


points(x = dat1$Depth.lon[dat1$Season=='Autum'], y=dat1$Depth.lat[dat1$Season=='Autum'],
       pch = 20,
       col = dat1$Col)

points(meta2,
       pch = 18,
       col = 'black')



# Compare observed values to modeled values
# using a fucking Man-whitney U

CPOD_preds=data.frame(SlopeMap=meta2$SlopeMap, 
                        Depth_m=meta2$Depth_m,
                        DistToSalmonRun=meta2$DistToSalmonRun,
                        BNDTotOffset=rep(0, nrow(meta2)),
                        DistToShore=meta2$DistToShore,
                        Season='Autum')

preds=data.frame(UnitLoc=meta2$UnitLoc,
  fit=predict(modlist_spatial[[11]], CPOD_preds, type='response'))

CPOD_obs=aggregate(data=subset(OccTable_daily, Season=='Summer'), BNDTotOffset~UnitLoc, FUN=mean)
CPOD_obs=merge(CPOD_obs, preds, by='UnitLoc')
CPOD_obs$fitorder=order(CPOD_obs$fit)
CPOD_obs$obsorder=order(CPOD_obs$BNDTotOffset)

wilcox.test(CPOD_obs$BNDTotOffset, CPOD_obs$fit, paired = TRUE)


# Add actual area monitored by the sensors
plot(NorthSea, 
     deep=-200, 
     shallow=0, step=50, lwd=0.5, drawlabel=TRUE, add=TRUE)


for(ii in 1:30){plot(spTransform(SpatialCircle(sp = spTransform(meta2[ii,], proj_UTM), r=2000), crs.geo), add=TRUE)}





# Try adding in the hypothetical locations
# Create a new dataframe to test autocorrelation of when 8 more deployments added
library(plyr)
meta2_newarray=rbind.fill(as.data.frame(meta2), as.data.frame(meta2_temp))
detach("package:plyr", unload=TRUE)

cor(meta2_newarray$Depth_m, meta2_newarray$DistToSalmonRun, method='pearson') # correlation -0.4
cor(meta2_newarray$Depth_m, meta2_newarray$SlopeMap, method='pearson') 
cor(meta2_newarray$Depth_m, meta2_newarray$DistToShore, method='pearson') 

cor(meta2_newarray$DistToSalmonRun, meta2_newarray$SlopeMap, method='pearson') 
cor(meta2_newarray$DistToSalmonRun, meta2_newarray$DistToShore, method='pearson') 

cor(meta2_newarray$SlopeMap, meta2_newarray$DistToShore, method='pearson') 




