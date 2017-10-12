library(plyr)
library(ggplot2)
rm(list=ls())
setwd("W:/KJP PHD/3-Detection Function/Propagation Model/Hourly NL 41khz thrdOctBand")
temp = list.files(pattern="*.csv")
loc=c("Stb_05","Fra_05","Arb_10","Cru_05","Hel_15",
      "StA_10", "Sto_05","SpB_10","Lat_05","Cro_15")

idx=c(9,5,10,8,2,4,7,3,6,1)

df=read.csv(temp[idx[1]])
df$DepLoc=loc[idx[1]]


for (ii in 2:length(temp)){
  dftemp=read.csv(temp[idx[ii]])
  dftemp$DepLoc=loc[idx[ii]]
  #dftemp$NLasl=dftemp$NLasl-(min(dftemp$NLasl)-40.44589)
  df=rbind(df,dftemp)
}


df <- within(df, DepLoc <- factor(DepLoc, levels=loc[idx[10:1]]))

min_nl=min(df$NLasl)

df1=subset(df, df$NLasl<41)
dfreal=subset(df, df$NLasl<60)






# 
# # Interleaved histograms
# ggplot(dat, aes(x=rating, fill=cond)) +
#   geom_histogram(binwidth=.5, position="dodge")
# 
# ggplot(df, aes(x=NLasl, fill=DepLoc)) +
#   geom_histogram(binwidth=.5, position="dodge")

cdat <- ddply(df1, "DepLoc", summarise, NLAsl.min=min(NLasl))


ggplot(df1, aes(x=NLasl, colour=DepLoc)) + geom_density(alpha=.3) +
geom_vline(data=cdat, aes(xintercept=NLAsl.min,  colour=DepLoc),
           linetype="dashed", size=1)


ggplot(df1, aes(x=NLasl, fill=DepLoc)) + geom_density(alpha=.3) +
 geom_vline(data=cdat, aes(xintercept=NLAsl.min,  colour=DepLoc),
           linetype="dashed", size=1)


ggplot(data=df1, aes(NLasl, fill=DepLoc)) + geom_bar(position="dodge")





ggplot(data=dfreal, aes(NLasl, DepLoc, color=NLasl)) +geom_jitter()




ggplot(dfreal, aes(NLasl)) + geom_bar() +
  facet_wrap(~ DepLoc)


theme_set(theme_grey(base_size = 18))
qplot(data=dfreal, NLasl, DepLoc, geom = "jitter", ylab='', 
      color=NLasl, main = 'Median Hourly Noise Level in 40kHz Band') + 
  labs(x=expression(NL[ASL]))

ggplot(dfreal, aes(x = name, y = val)) + theme_bw() + geom_bar(stat = "identity")

p=ggplot(dfreal, aes(NLas, DepLoc, color=NLasl))
p=p+layer(geom ="point")

# Plot lines of noise leve
df2=subset(df, MatlabDate<735540)
ggplot(data=df2, aes(x=MatlabDate, y=NLasl, group=DepLoc, colour=DepLoc)) +
  geom_line() +
  geom_point()

# Boxplot of noise levels
p <- ggplot(df2, aes(DepLoc, NLasl))
p + geom_boxplot()


# look at noise relationship between locations
df3=subset(df, NLasl<44 & MatlabDate<735540)
# Reshape dataframe

nldf=data.frame(matrix(0, ncol = 10, nrow = (max(df3$MatlabDate)*24
                                             -min(df3$MatlabDate)*24)))
colnames(nldf)=loc[idx]

outliers=data.frame(matrix(ncol =5, nrow=10), row.names =loc[idx] )
colnames(outliers)=c('Median', 'IQR', 'N', 'Ntime', 'PropTime')

for(ii in 1:10){

  nl=df3$NLasl[df$DepLoc==loc[idx[ii]]]
  nl=nl[!is.na(nl)]
  outliers$Median[ii]=median(nl)
  outliers$IQR[ii]=IQR(nl)
  outliers$N[ii]=sum(nl>quantile(nl, 0.75)+1.5*IQR(nl))
  outliers$Ntime[ii]=length(nl)-sum(is.na(nl))
  outliers$PropTime[ii]=outliers$N[ii]/outliers$Ntime[ii]
  
  
  nldf[1:length(nl),ii]=nl
}

nldf=nldf[1:1095,]
newdata <- nldf[order(nldf$Lat_05),] 
pairs(newdata)

# look at between group and within group variance
medianvar=var(sapply(nldf, median, na.rm=TRUE))
groupvar=sapply(nldf, var, na.rm=TRUE))

# Coefficient of variation for each site
sapply(nldf, sd, na.rm=TRUE)/ sapply(nldf, mean, na.rm=TRUE)



#### Nooo, that's not really considering it
# Look at the IQR divided by the median
sapply(nldf, IQR, na.rm=TRUE)/sapply(nldf, median, na.rm=TRUE)

IQR(sapply(nldf, IQR, na.rm=TRUE))/median(sapply(nldf, IQR, na.rm=TRUE))


