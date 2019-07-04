require(e1071)
require(doParallel)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Get HAB
HAB = read.csv("C:/Users/xswang/HAB Research/West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1

ndbc = read.csv("NDBC_station_data.csv")
ndbc = ndbc[,!is.na(colSums(ndbc))]
ndbc$X = NULL
usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL
#####################################################################################
##Log scale to adjust discharge (add to avoid negatives)
usgs[,1:(ncol(usgs)-1)] = log10(usgs[,1:(ncol(usgs)-1)]+4)
#####################################################################################
par(mfrow=c(2,1))
##Show's seasonlity of wind measurements
mon = as.numeric(format(HAB$Date, "%m"))
nstates = rep(0, 12)
for(i in 1:nrow(HAB)){
  if(states[i]==1) {
    nstates[mon[i]] = nstates[mon[i]]+1
  }
}

par(mfrow=c(2,1))
for(i in 1:ncol(ndbc)){
  novals = rep(NA, 12)
  HABvals = rep(NA, 12)
  for(j in 1:12){
    novals[j] = mean(ndbc[states==0 & mon==j, i])
    HABvals[j] = mean(ndbc[states==1 & mon==j, i])
  }
  plot(1:12,nstates, main=colnames(ndbc)[i],type="l")
  plot(1:12, HABvals, type="l", col="red")
}
#####################################################################################
plot(as.numeric(format(HAB$Date, "%m")), ndbc$venf1.EW)
##Compare wind speeds
tiff(file = "NS v EW.tiff", height = 7.5,width = 7.5,  units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(3,2))
pc = rep("blue", nrow(ndbc))
pc[states==1]="red"
plot(ndbc$smkf1.NS, ndbc$smkf1.EW, col=pc)
plot(ndbc$venf1.NS, ndbc$venf1.EW, col=pc)
plot(ndbc$cdrf1.NS, ndbc$cdrf1.EW, col=pc)
plot(ndbc$X42039.NS, ndbc$X42039.EW, col=pc)
plot(ndbc$X42003.NS, ndbc$X42003.EW, col=pc)
dev.off()
