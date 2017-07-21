
rm(list=ls())
library(corrgram) # for matrix image
library(vegan) # for mantel test
#########################################################################################
# Load the 2013 Data #
########################################################################################


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


rm(Pdet1,Pdet2,Pdet3,Pdet4,Pdet5,Pdet6,Pdet7,Pdet8,Pdet9,Pdet10, cv, sd)

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

# Add Date
Pdetdf1$Date=as.Date(Pdetdf1$MatlabDate, origin = '0000-01-01')
Pdetdf1$Hr=round((Pdetdf1$MatlabDate-floor(Pdetdf1$MatlabDate))*24)

######################################################################################
# Load 2015 Data #

temp = list.files(path='Z:/2015SM_output', pattern="*.txt")

pdetdf2015=Pdetdf[1,]
pdetdf2015=pdetdf2015[-1,]


for (ii in 1:length(temp)){
  
  # stop gap until new Pet files are made =
  df=(t(read.csv(paste('Z:/2015SM_output/', temp[ii], sep = ""), sep=',', header=F)))
  m=as.data.frame(matrix(nrow = nrow(df)-1, ncol=ncol(Pdetdf)))
  colnames(m)=colnames(Pdetdf)
  
  m[,c('MatlabDate', 'MedianNoiseLevel')]=df[2:nrow(df),]
  m$UnitLoc=substr(temp[ii],15,20)
 
  pdetdf2015=rbind(pdetdf2015, m)
  rm(list=c('df', 'm'))
  
  print(ii)

}

pdetdf2015$Date=as.POSIXct((pdetdf2015$MatlabDate - 719529)*86400, origin = "1970-01-01", tz = "UTC")
pdetdf2015$Year=2015

NL=pdetdf2015$MedianNoiseLevel
pdetdf2015$MedianNoiseLevel=NL-10*log10(9287.75)+10*log10(100000)-12



#####################################################################################

Pdetdf=rbind(Pdetdf, Pdetdf1)

rm(Pdet1,Pdet2,Pdet3,Pdet4,Pdet5,Pdet6,Pdet8,Pdet9,Pdet10, cv, sd, Pdetdf1)


# Create large matrix for the correlation test
tempdf=subset(Pdetdf, UnitLoc==unique(Pdetdf$UnitLoc)[1],select=c(MatlabDate, MedianNoiseLevel))
colnames(tempdf)[2]=unique(Pdetdf$UnitLoc)[1]


# Reshape the dataframe... undoubtedy a better way to do this
for (ii in 2:length(unique(Pdetdf$UnitLoc))){
  x1=subset(Pdetdf, UnitLoc==unique(Pdetdf$UnitLoc)[ii],
            select=c(MatlabDate, MedianNoiseLevel))
  colnames(x1)[2]=unique(Pdetdf$UnitLoc)[ii]
  tempdf=merge(tempdf, x1, by='MatlabDate', all = TRUE)
  rm(list=c('x1'))  
}

# Create a new df that accounts for trends using first differences
# y'(t) = y(t) – y(t-1)

df.notrend=tempdf[,1:2]
df.notrend$Lat_05=c(0, df.notrend$Lat_05[2:nrow(df.notrend)]-
                       df.notrend$Lat_05[1:nrow(df.notrend)-1])

for (ii in 2:10){

    df.notrend=cbind(df.notrend, 
                     c(0, tempdf[2:nrow(df.notrend),ii+1]-
                         tempdf[1:nrow(df.notrend)-1,ii+1]))
    colnames(df.notrend)[ii+1]=colnames(tempdf)[ii+1]
 
}


# Correlation matrix
Pairwise_Correlation=cor(tempdf[,2:11], use="pairwise.complete.obs")
summary(Pairwise_Correlation[Pairwise_Correlation<1])

corrgram(Pairwise_Correlation, 
         lower.panel=panel.pie, upper.panel=NULL, 
         type = 'corr',
         main="2013 Noise Levels Correlations ") 

# library(GGally)
# ggpairs(tempdf[1:1095,colnames(tempdf)[2:11]],)
# # 
# pm <- ggpairs(subset(Pdetdf, MedianNoiseLevel>80), columns = 2:3, ggplot2::aes(colour=UnitLoc))
# p_(pm)
# 
# ggpairs(subset(Pdetdf, MedianNoiseLevel>80), columns = 2:3)
# pm <- ggpairs(subset(Pdetdf, MedianNoiseLevel>80), columns = 2:3,
#               ggplot2::aes(colour=UnitLoc))

############################################
# Try mixing up the noise levels completly #
############################################
# First remove the NA values #

df2=tempdf[complete.cases(tempdf),]
df2 <- cbind(tempdf$MatlabDate, apply(tempdf[,2:11], 2, sample))


R_PW_corr=cor(df2[,2:11], use="pairwise.complete.obs")
summary(R_PW_corr[R_PW_corr<1])

corrgram(R_PW_corr, upper.panel=panel.pie, type = 'corr') 

###################################################################
# Mantel test with just the correlation scores #
###################################################################

# 1) create a matrix of distances between units
DepInfo=read.csv('W:/KJP PHD/Deployment Information/Copy of Marine Scotland 2013-14 CPOD summary data (2).csv')
Pdetdf_full=merge(Pdetdf, DepInfo[!duplicated(DepInfo$UnitLoc), c('Lat', 'Long', 'UnitLoc')])


# Dist matrix
temp=DepInfo[!duplicated(DepInfo$UnitLoc), c('Lat', 'Long', 'UnitLoc')]

selected=unique(unique(Pdetdf$UnitLoc))
temp1=temp[temp$UnitLoc %in% selected,]



# Distance matrix
distmat_betweenUnits=scale(as.matrix(dist(temp1)))
distmat_betweenUnits=(as.matrix(dist(temp1)))

# Do mantel test
# null hypothesis: the two matricies are unrelated
mantel(distmat_betweenUnits, Pairwise_Correlation)
plot(mantel.correlog(D.geo = distmat_betweenUnits, D.eco = Pairwise_Correlation, n.class = 7))

# Do mantel test on random correlations
mantel(distmat_betweenUnits, R_PW_corr)

####################################################################
