
# Load files and make dummy date
rm(list=ls())
library(ggplot2)
meta=read.csv('W:\\KJP PHD\\Deployment Information\\CPODs for Kaitlin.csv')

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

meta$UnitLoc=factor(meta$UnitLoc, levels=rev(level_names))


meta$year=substr(as.character(meta$Usable.from.date), 7,10)
meta$dummystart=as.Date(paste(substr(meta$Usable.from.date, 1, 6), '2013', sep=''),
                        "%d/%m/%Y", tz = "GMT")
meta$dummyend=as.Date(paste(substr(meta$Usable.until.date, 1, 6), '2013', sep=''), 
                      "%d/%m/%Y", tz = "GMT")

meta$Usable.from.date=as.Date(meta$Usable.from.date, "%d/%m/%Y", tz = "GMT")
meta$Usable.until.date=as.Date(meta$Usable.until.date, "%d/%m/%Y", tz = "GMT")
meta$dummyend[meta$year==""]=NA
meta$dummystart[meta$year==""]=NA

ggplot()+geom_rect(data=meta[!is.na(meta$dummystart),],aes(
  y=UnitLoc,
  ymin= as.numeric(UnitLoc)+.45,
  ymax= as.numeric(UnitLoc)-.45,
  xmin=dummystart,
  xmax=dummyend, fill=year), alpha=0.4) +
  theme_minimal() +
  scale_x_date(limits = c(as.Date("15/03/2013", format="%d/%m/%Y"),
                          as.Date("28/12/2013", format="%d/%m/%Y"))) +
  labs(y='Deployment Location', x='Date') +
  ggtitle('ECoMMAS Data Coverage')



ggplot(meta, aes(UnitLoc, us)) +
  geom_segment(aes(xend = x, yend = 0), size = 10, lineend = "butt")
