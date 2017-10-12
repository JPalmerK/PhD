

hist_AreaMonitored <- read.csv("C:/Users/charlotte/Desktop/Geek Stick Material/hist_AreaMonitored.csv", 
                               header=FALSE)

hist_AreaMonitored$V1=as.numeric(hist_AreaMonitored$V1)
br=seq(4,24)
h1=hist(subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==1), breaks=br, plot = FALSE)
h2=hist(subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==2), breaks=br, plot = FALSE)
h3=hist(subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==3), breaks=br, plot = FALSE)
h4=hist(subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==4), breaks=br, plot = FALSE)
h5=hist(subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==5), breaks=br, plot = FALSE)
h6=hist(subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==6), breaks=br, plot = FALSE)
h7=hist(subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==7), breaks=br, plot = FALSE)
h8=hist(subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==8), breaks=br, plot = FALSE)
h9=hist(subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==9), breaks=br, plot = FALSE)
h10=hist(subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==10), breaks=br, plot = FALSE)

barplot(height = rbind(h1$counts, h2$counts, h3$counts,
                       h4$counts, h5$counts, h6$counts), 
        col=c('blue', 'gray','black', 'green', 
              'red', 'yellow'), beside=T,
        names.arg=as.character(round(br[2:length(br)],0)),
        las=1, legend = c('StB_05','Fra_05','Cru_05', 'Hel_15'),
        ylab = 'Number of Observations', xlab = 'Area Monitored',
        main = 'Noise Levels at Different Sites')




aa=subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==1)
plot(aa, type='l', main='Area Monitored By SM2Ms',
     ylim=c(4,24), xlab='Time', ylab='Area Monitored (km2)')


for(ii in 2:10){
  aa=subset(hist_AreaMonitored$V1, hist_AreaMonitored$V2==ii)
  lines(aa, type = 'l', col=ii)
}







h2=hist(NL2, breaks=nl.breaks, plot=FALSE)
h4=hist(NL4, breaks=nl.breaks, plot=FALSE)
h7=hist(NL7, breaks=nl.breaks, plot=FALSE)