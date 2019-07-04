require(e1071)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Read in data
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1

alldata = read.csv("totaldata.csv")
alldata$X = NULL
for(i in 1:ncol(alldata)) alldata[is.na(alldata[,i]),i]=mean(alldata[,i],na.rm=T)
alldata[,10:23] = log10(alldata[,10:23]+4)
alldata$State = as.factor(alldata$State)
scalefactors = scale(alldata[,1:(ncol(alldata)-1)])
alldata[,1:(ncol(alldata)-1)] = scale(alldata[,1:(ncol(alldata)-1)])
#####################################################################################
ndbccolumns = 17
meanvals = rep(NA, ncol(alldata)-1)
for(i in 1:(ncol(alldata)-1)){
  meanvals[i] = mean(alldata[alldata$State=="1",i])
}
#####################################################################################
relFreq <- function(counts){
  return(counts/sum(counts))
}
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}
getmode = function(v){
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v,uniqv)))]
}
##Plot histograms
tiff(file = "NDBC histograms.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(5,4))
for(i in 1:9) {
  d1 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i],i],scalefactors)
  d2 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i] & 
                             alldata[,ncol(alldata)]=="1",i],scalefactors)
  h1 = hist(d1, plot = F)
  h2 = hist(d2, plot =F, breaks=h1$breaks)
  h2$counts <- relFreq(h2$counts)
  #h2$counts = h2$counts/h1$counts*relFreq(h1$counts)
  h1$counts <- relFreq(h1$counts)
  plot(h1, main=colnames(alldata)[i], xlab=colnames(alldata)[i])
  plot(h2, main=paste0(colnames(alldata)[i], " during HABs"),col="red", xlab=colnames(alldata)[i])
}
dev.off()

tiff(file = "NDBC NS ATMP histograms.tiff", width =15, height = 9, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(3,4))
for(i in c(1,2,4,6,8)) {
  d1 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i],i],scalefactors)
  d2 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i] & 
                             alldata[,ncol(alldata)]=="1",i],scalefactors)
  h1 = hist(d1, plot = F)
  h2 = hist(d2, plot =F, breaks=h1$breaks)
  h2$counts <- relFreq(h2$counts)
  #h2$counts = h2$counts/h1$counts*relFreq(h1$counts)
  h1$counts <- relFreq(h1$counts)
  plot(h1, main=colnames(alldata)[i], xlab=colnames(alldata)[i])#,ylim=c(0,0.5))
  plot(h2, main=paste0(colnames(alldata)[i], " during HABs"),col="red")
}
dev.off()

tiff(file = "NDBC EW histograms.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(2,4))
for(i in c(2,5,7,9)) {
  d1 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i],i],scalefactors)
  d2 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i] & 
                             alldata[,ncol(alldata)]=="1",i],scalefactors)
  h1 = hist(d1, plot = F)
  h2 = hist(d2, plot =F, breaks=h1$breaks)
  h2$counts <- relFreq(h2$counts)
  #h2$counts = h2$counts/h1$counts*relFreq(h1$counts)
  h1$counts <- relFreq(h1$counts)
  plot(h1, main=colnames(alldata)[i], xlab=colnames(alldata)[i])#,ylim=c(0,0.5))
  plot(h2, main=paste0(colnames(alldata)[i], " during HABs"),col="red")
}
dev.off()

tiff(file = "USGS histograms 1.tiff", width =15, height = 15, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(4,4))
for(i in 10:17) {
  d1 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i],i],scalefactors)
  d2 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i] & 
                             alldata[,ncol(alldata)]=="1",i],scalefactors)
  h1 = hist(d1, plot = F)
  h2 = hist(d2, plot =F, breaks=h1$breaks)
  h2$counts <- relFreq(h2$counts)
  #h2$counts = h2$counts/h1$counts*relFreq(h1$counts)
  h1$counts <- relFreq(h1$counts)
  plot(h1, main=colnames(alldata)[i], xlab=colnames(alldata)[i])#,ylim=c(0,0.5))
  plot(h2, main=paste0(colnames(alldata)[i], " during HABs"),col="red")
}
dev.off()
tiff(file = "USGS histograms 2.tiff", width =15, height = 15, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(3,4))
for(i in 18:23) {
  d1 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i],i],scalefactors)
  d2 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i] & 
                             alldata[,ncol(alldata)]=="1",i],scalefactors)
  h1 = hist(d1, plot = F)
  h2 = hist(d2, plot =F, breaks=h1$breaks)
  h2$counts <- relFreq(h2$counts)
  #h2$counts = h2$counts/h1$counts*relFreq(h1$counts)
  h1$counts <- relFreq(h1$counts)
  plot(h1, main=colnames(alldata)[i], xlab=colnames(alldata)[i])#,ylim=c(0,0.5))
  plot(h2, main=paste0(colnames(alldata)[i], " during HABs"),col="red")
}
dev.off()
############################################################################################
##Histograms during (non-) HAB events
tiff(file = "SJWP NDBC histograms.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(3,2))
for(i in 1:3) {
  d1 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i],i],scalefactors)
  d2 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i] & 
                             alldata[,ncol(alldata)]=="1",i],scalefactors)
  h1 = hist(d1, plot = F)
  h2 = hist(d2, plot =F, breaks=h1$breaks)
  h2$counts <- relFreq(h2$counts)
  #h2$counts = h2$counts/h1$counts*relFreq(h1$counts)
  h1$counts <- relFreq(h1$counts)
  plot(h1, main=colnames(alldata)[i], xlab=colnames(alldata)[i])
  plot(h2, main=paste0(colnames(alldata)[i], " during HABs"),col="red", xlab=colnames(alldata)[i])
}
dev.off()

tiff(file = "SJWP USGS histograms.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(3,2))
for(i in c(11,12,15)) {
  d1 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i],i],scalefactors)
  d2 = undoscale(i,alldata[getmode(alldata[,i])!=alldata[,i] & 
                             alldata[,ncol(alldata)]=="1",i],scalefactors)
  h1 = hist(d1, plot = F)
  h2 = hist(d2, plot =F, breaks=h1$breaks)
  h2$counts <- relFreq(h2$counts)
  #h2$counts = h2$counts/h1$counts*relFreq(h1$counts)
  h1$counts <- relFreq(h1$counts)
  plot(h1, main=colnames(alldata)[i], xlab=colnames(alldata)[i])#,ylim=c(0,0.5))
  plot(h2, main=paste0(colnames(alldata)[i], " during HABs"),col="red")
}
dev.off()
