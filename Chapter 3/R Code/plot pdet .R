



library(plyr)
library(ggplot2)
rm(list=ls())

setwd("W:/KJP PHD/3-Detection Function/R Code")
PdetPath='W:/KJP PHD/6-Bayesian Habitat Use/Pdet at Time/2013'

temp = list.files(path = PdetPath, pattern="*.csv")
idx=c(6,5,2,7,4,3,10,1,8,9)


path=paste(PdetPath,'/',temp[idx[1]], sep = '')
Pdet=read.csv(path)
Pdet$DepLoc=substr(temp[idx[1]], 1, 6)


for (ii in 2:length(temp)){
  path=paste(PdetPath,'/',temp[idx[ii]], sep = '')
  pdettemp=read.csv(path)
  pdettemp$DepLoc=substr(temp[idx[ii]], 1, 6)
  Pdet=rbind(Pdet,pdettemp)
  rm(pdettemp)
}

p <- ggplot(Pdet, aes(DepLoc, MedianPdet))
p + geom_boxplot()





