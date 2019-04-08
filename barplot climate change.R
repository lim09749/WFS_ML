rm(list=ls())
require(e1071)
####################################################################################################
setwd("C:/Users/xswang/HAB Research")

tempmin = seq(from=0, to=1, length.out=83)
tempavg = seq(from=0, to=2.5, length.out=83)
tempmax = seq(from=0, to=6, length.out=83)

prmin   = seq(from=0, to=5, length.out=83)
pravg   = seq(from=0, to=10, length.out=83)
prmax   = seq(from=0, to=20, length.out=83)
####################################################################################################
##Get HAB

##Form alldata
HAB = read.csv("C:/Users/xswang/HAB Research/West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Abundance_cells[is.na(HAB$Abundance_cells)] = 0
HAB$Date = as.Date(HAB$Date)

load(file="radialSVM.RData")
ndbc = read.csv("NDBC_station_data.csv")
ndbc = ndbc[,!is.na(colSums(ndbc))]
ndbc$X = NULL
ndbcold=ndbc

usgsold   = read.csv("USGS_discharge_total.csv")
usgsold$X = NULL
usgs = usgsold

usgs[,1:(ncol(usgs)-1)] = log10(usgs[,1:(ncol(usgs)-1)]+20)
alldata = cbind(ndbc,usgs)
scalefactors = scale(alldata[,1:(ncol(alldata)-1)])

months = as.numeric(format(HAB$Date,"%m"))
inds   = !is.na(match(months,c(1:2,9:12))) & HAB$Date>as.Date("2008-01-01")
#####################################################################################
tiff(file = "comparison.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
statesmin =rep(NA, 83)
statesavg =rep(NA, 83)
statesmax =rep(NA, 83)
statesmin[1] = length(which(HAB$Abundance_cells[inds]>=10^5))
statesavg[1] = length(which(HAB$Abundance_cells[inds]>=10^5))
statesmax[1] = length(which(HAB$Abundance_cells[inds]>=10^5))

for(i in 2:83){
  usgs=usgsold
  ndbc=ndbcold
  ndbc$venf1.ATMP = ndbc$venf1.ATMP+tempmin[i]
  usgs = (prmin[i]+1)*usgs
  usgs[,1:(ncol(usgs)-1)] = log10(usgs[,1:(ncol(usgs)-1)]+20)
  usgs = usgs[inds,]
  ndbc   = ndbc[inds,]
  alldata = cbind(ndbc,usgs)
  alldata[,1:(ncol(alldata)-1)] = scale(alldata[,1:(ncol(alldata)-1)],
                                        center=attr(scalefactors,"scaled:center"),
                                        scale=attr(scalefactors,"scaled:scale"))
  pred = predict(HABsvm,alldata)
  statesmin[i] = length(which(pred=="1"))
}

for(i in 2:83){
  usgs=usgsold
  ndbc=ndbcold
  ndbc$venf1.ATMP = ndbc$venf1.ATMP+tempavg[i]
  usgs = (pravg[i]+1)*usgs
  usgs[,1:(ncol(usgs)-1)] = log10(usgs[,1:(ncol(usgs)-1)]+4)
  usgs = usgs[inds,]
  ndbc   = ndbc[inds,]
  alldata = cbind(ndbc,usgs)
  alldata[,1:(ncol(alldata)-1)] = scale(alldata[,1:(ncol(alldata)-1)],
                                        center=attr(scalefactors,"scaled:center"),
                                        scale=attr(scalefactors,"scaled:scale"))
  pred = predict(HABsvm,alldata)
  statesavg[i] = length(which(pred=="1"))
}

for(i in 2:83){
  usgs=usgsold
  ndbc=ndbcold
  ndbc$venf1.ATMP = ndbc$venf1.ATMP+tempmax[i]
  usgs = (prmax[i]+1)*usgs
  usgs[,1:(ncol(usgs)-1)] = log10(usgs[,1:(ncol(usgs)-1)]+4)
  usgs = usgs[inds,]
  ndbc   = ndbc[inds,]
  alldata = cbind(ndbc,usgs)
  alldata[,1:(ncol(alldata)-1)] = scale(alldata[,1:(ncol(alldata)-1)],
                                        center=attr(scalefactors,"scaled:center"),
                                        scale=attr(scalefactors,"scaled:scale"))
  pred = predict(HABsvm,alldata)
  statesmax[i] = length(which(pred=="1"))
}
####################################################################################
dev.off()
