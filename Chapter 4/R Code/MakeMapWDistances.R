
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




map.corunty <- map("world", c("UK"), plot = FALSE, fill = TRUE, res = 0)

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
meta2=read.csv('W:\\KJP PHD\\Deployment Information\\SlopeAndAspect.csv')
meta2$UnitLoc=factor(meta2$UnitLoc, levels=level_names)
meta2$DistToSalmonRun=0
meta2$RiverName='blarg'

# 2) Data prep, make map convert shit ##############


# calculate distance to nearest salmon river
for(ii in 1:30){
  temp=rep(0, nrow(river_locs))
  
  for(jj in 1:nrow(river_locs)){
    temp[jj]=distm (c(meta2$Lon[ii], meta2$Lat[ii]), c(river_locs$lonDeg[jj], river_locs$LatDeg[jj]), fun = distHaversine)
  }
  
  meta2$DistToSalmonRun[ii]=min(temp)
  meta2$RiverName[ii]=as.character(river_locs$Rivername[which.min(temp)])
  
}



# 3) Create Data Grid (WGS) and calculate distance to rivers #######

river_coords=river_locs[, c('lonDeg','LatDeg' )]
coordinates(river_coords)=c('lonDeg','LatDeg' )


crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
proj4string(river_coords) <- crs.geo  # define projection system of our data
summary(river_coords)

coordinates(meta2)=meta2[,c('Lon','Lat' )]
proj4string(meta2) <- crs.geo  # define projection system of our data




dat=expand.grid(Lon=seq(min(meta2$Lon)-1,
                        max(meta2$Lon)+.5, length.out = 500),
                Lat=seq(min(meta2$Lat)-.5,
                        max(meta2$Lat)+.5, length.out = 500))

coordinates(dat)=c( 'Lon','Lat')
proj4string(dat) <-crs.geo

# Calculated distances between feature locations and grid locations
River_dists=matrix(nrow = length(dat), ncol=nrow(river_locs))
for(ii in 1:nrow(river_locs)){River_dists[,ii]=distRhumb(river_coords[ii], dat)}
River_dists=data.frame(River_dists)
colnames(River_dists)=river_locs$Rivername

# Calculate teh distance between CPOD Locations and Grid Data
Cpod_dists=matrix(nrow = length(dat), ncol=nrow(meta2))
for(ii in 1:nrow(meta2)){ Cpod_dists[,ii]=distRhumb(coordinates(meta2)[ii,], dat)}
Cpod_dists=data.frame(Cpod_dists)
colnames(Cpod_dists)=meta2$UnitLoc


# Filter data points that are within 2km of a C-POD
Cpod_detRange=which(Cpod_dists<=2000)



# Projection to use for spatial objects
proj_UTM <- CRS("+proj=utm +zone=30 ellps=WGS84")


# Create Spatial Polygons
# Indeces for points in detection range of CPODs
inside.circle_idx=numeric()

for (ii in 1:nrow(meta2)){
  
  
  # Radius 1 and 2 of the annulus centred on each MARU. To increase legibility of confirmation plots
  # increase 0.5m to 50m
  r=2000 #2km
  
  CPOD_loc=meta2[ii,]
  CPOD_loc_UTM=spTransform(CPOD_loc, proj_UTM)
  
  
  ## Circle
  circle <- SpatialCircle(CPOD_loc_UTM, r)
  
  ## Transform circle to WGS
  circle_wgs=spTransform(circle, crs.geo)
  
  
  
  inside.circle <- !is.na(over(dat, as(circle_wgs, "SpatialPolygons")))
  
  inside.circle_idx=c(inside.circle_idx,
                      which(!is.na(over(dat, as(circle_wgs, "SpatialPolygons")))))
  
}



points(dat[inside.circle_idx, ], pch=16, col="red")



# Extract data points within range of each river system ###
RiverRanges=data.frame(aggregate(data=meta2, DistToSalmonRun~RiverName, FUN=range))
RiverRanges=merge(RiverRanges, river_locs, 
                  all.x=TRUE, by.x=c('RiverName'), by.y='Rivername')

# Convert to spatial 

coordinates(RiverRanges)=c('lonDeg','LatDeg' )
proj4string(RiverRanges) <- crs.geo  # define projection system of our data
summary(river_coords)



# Create Spatial Polygons
# 
# # Ploygone for the UK
# UK_map= map("world", region="UK", fill=TRUE, plot=FALSE)
# IDs <- sapply(strsplit(UK_map$names, ":"), function(x) x[1])
# UkPoly <- map2SpatialPolygons(UK_map, IDs=IDs, 
#                               proj4string=crs.geo)



# Try higher resolution Scotland Map
gadm <- readRDS("W:/KJP PHD/4-Bayesian Habitat Use/R Code/GBR_adm1.rds") 
gadm=spTransform(gadm, crs.geo)
ScotMap <- gadm[gadm$NAME_1 == "Scotland", ] #Possibly a ploygone already?



# Indeces for points in detection range of CPODs
inside.river_idx=numeric()

for (ii in 1:nrow(RiverRanges)){
  
  
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
  # 
  # ## Plot checking
  # plot(circle2_wgs)
  # plot(circle1_wgs, add=TRUE)
  # points(subset(meta2, RiverName==RiverRanges$RiverName[ii]), pch=18)
  
  # Points must be both within range and outside of the UK
  inside.river= which(!is.na(over(dat, as(circle2_wgs, "SpatialPolygons"))) &
                        is.na(over(dat, as(circle1_wgs, "SpatialPolygons"))))
  
  # Make sure points are also in the water
  Swimming=which(!is.na(over(dat[inside.river], as(ScotMap, "SpatialPolygons"))))
  
  
  inside.river_idx=c(inside.river_idx, inside.river[-Swimming])
  
  
  
  River_dists[inside.river[-Swimming],RiverRanges$RiverName[ii]]=NA
  
  # # More plot checking, use 2000 and 5000 for plot checking
  # points(dat[inside.river[-Swimming]], pch=16, col="red") #good
  
  print(ii)
  
}

River_dists=River_dists[, c('Esk', 'Dee', 'Spey', 'Tay Firth', 'Tweed', 'Cromarty Firth')]

plot(ScotMap)
points(dat[inside.river_idx], pch=16, col="red")

plot(ScotMap)
points(dat[is.na(River_dists$Dee)], col='blue')


# Load and extract Bathyemtry marmap #################################################

swimmingdat=dat[inside.river_idx]


# Load the Bathymetry Grid using marmap
getNOAA.bathy(lon1 = min(coordinates(dat)[,1]), lon2 = max(coordinates(dat)[,1]),
              lat1 = min(coordinates(dat)[,2]), lat2 = max(coordinates(dat)[,2]),
              resolution = 1)->NorthSea

# Create a raster version
NorthSea_raster=as.raster(NorthSea)
NorthSea_raster=projectRaster(NorthSea_raster, crs = crs.geo)
NorSea_spatialGrid=as.SpatialGridDataFrame(bathy = NorthSea)

plot(NorthSea, image = TRUE, lwd = 0.3)
plot(NorthSea_raster, image = TRUE, lwd = 0.3,
     xlab = "", ylab = "", axes = FALSE)

# Extract approximate depths
SurveyDepths=get.depth(NorthSea, coordinates(swimmingdat), locator = FALSE)

# Create slope values
Slope=terrain(NorthSea_raster, opt = 'Slope', unit = 'Radians', neighbors = 4)


SurveySlope=raster::extract(Slope, swimmingdat)
CPOD_slope=raster::extract(Slope, meta2)
CPOD_Depth=get.depth(NorthSea,  coordinates(meta2), locator = FALSE)

# Make the dataframe
Slope2=numeric()
RiverName=character()
DistToSalmonRun=numeric()

# Extract river dists for plotting
for(ii in 1:ncol(River_dists)){

  tempDist=  River_dists[!is.na(River_dists), ii] # need to check the river names
  RiverName= c( RiverName, rep(RiverRanges$RiverName[ii], length(tempDist)))
  DistToSalmonRun=c(DistToSalmonRun, tempDist)
  
  
}
tempdf=data.frame(x=tempDist, y=RiverName, z=DistToSalmonRun)

tempdf=tempdf[!duplicated(tempdf),]

# Form dataframe
NewData=data.frame(Slope2=SurveySlope,
                   depth_m=SurveyDepths,
                   RiverName=RiverName)




# Check depths agree both values
meta=read.csv('W:/KJP PHD/Deployment Information/CPODs for Kaitlin.csv')
meta=meta[!duplicated(meta$UnitLoc),]

temp=cbind(meta$Depth_m, CPOD_Depth$depth)
