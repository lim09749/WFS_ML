require(e1071)
require(graphics)
require(grDevices)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Read in data

HAB = read.csv("C:/Users/xswang/HAB Research/West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Abundance_cells[is.na(HAB$Abundance_cells)] = 0
HAB$Date = as.Date(HAB$Date)


load(file="radialSVM.RData")
ndbc = read.csv("NDBC_station_data.csv")
ndbc = ndbc[,!is.na(colSums(ndbc))]
ndbc$X = NULL
usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL
HAB = read.csv("C:/Users/xswang/HAB Research/West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Abundance_cells[is.na(HAB$Abundance_cells)] = 0
HAB$Date = as.Date(HAB$Date)
months = as.numeric(format(HAB$Date,"%m"))

#####################################################################################
##Log scale to adjust discharge (add to avoid negatives)
usgs[,1:(ncol(usgs)-1)] = log10(usgs[,1:(ncol(usgs)-1)]+4)
#####################################################################################
alldata = cbind(ndbc,usgs)
for(i in 1:ncol(alldata)) alldata[,i] = as.numeric(as.character(alldata[,i]))
alldata$State = as.factor(alldata$State)
scalefactors = scale(alldata[,1:(ncol(alldata)-1)])
alldata[,1:(ncol(alldata)-1)] = scale(alldata[,1:(ncol(alldata)-1)])

meanvals = rep(NA, ncol(alldata)-1)
for(i in 1:(ncol(alldata)-1)){
  meanvals[i] = mean(alldata[alldata$State=="1",i])
}
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}

#####################################################################################
tiff(file = "contour.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))

par(mfrow=c(2,2))
numrows = 100
addto = seq(-2, 2, length.out=numrows)
meanvals = colMeans(alldata[!is.na(match(months,8:10)), 1:(ncol(alldata)-1)])
print(meanvals)

#####################################################################################
##Help function
getcontour = function(col1, col2,meanvals, numrows, alldata, addto, HABsvm){
  cmat = matrix(NA, nrow=numrows, ncol=numrows)
  for(i in 1:nrow(cmat)){
    for(j in 1:ncol(cmat)){
      curr = as.data.frame(rbind(meanvals,meanvals))
      colnames(curr) = colnames(alldata)[1:(ncol(alldata)-1)]
      curr[,col1] = curr[,col1]+addto[i]
      curr[,col2] = curr[,col2]+addto[j]
      cmat[i,j] = attr(predict(HABsvm,curr,
                               probability=T),"probabilities")[1,2]
    }
  }
  return(cmat)
}
#####################################################################################
##Collect alldata
cmat = array(NA, dim=c(numrows,numrows,4))

##Hillsborough and Peace
cmat[,,1] = getcontour(15, 11, meanvals, numrows, alldata, addto, HABsvm)

contour(x=undoscale(15, addto+meanvals[15], scalefactors),
        y=undoscale(11, addto+meanvals[11], scalefactors),
        z=cmat[,,1],
        xlab= "Hillsborough River (cfs)",
        ylab= "Peace River (cfs)",
        axes=F,
        cex.main=2,cex.lab=1.5,
        labcex=1.25)
axis(side = 1, at = axTicks(1), labels=round(10^axTicks(1)),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=round(10^axTicks(2)),cex.axis=1.5)
#####################################################################################
##Hillsborough and venf1.NS
cmat[,,2] = getcontour(15, 1, meanvals, numrows, alldata, addto, HABsvm)

contour(x=undoscale(15, addto+meanvals[15], scalefactors),
        y=undoscale(1, addto+meanvals[1], scalefactors),
        z=cmat[,,2],
        xlab= "Hillsborough River (cfs)",
        ylab= "NS wind at venf1 (m/s)",
        axes=F,
        cex.main=2,cex.lab=1.5,
        labcex=1.25)
axis(side = 1, at = axTicks(1), labels=round(10^axTicks(1)),cex.axis=1.5)
axis(side = 2,cex.axis=1.5)
#####################################################################################
##Hillsborough and cdrf1.EW
cmat[,,3] = getcontour(15, 5, meanvals, numrows, alldata, addto, HABsvm)

contour(x=undoscale(15, addto+meanvals[15], scalefactors),
        y=undoscale(5, addto+meanvals[5], scalefactors),
        z=cmat[,,3],
        xlab= "Hillsborough River (cfs)",
        ylab= "EW wind at cdrf1 (m/s)",
        axes=F,
        cex.main=2,cex.lab=1.5,
        labcex=1.25)
axis(side = 1, at = axTicks(1), labels=round(10^axTicks(1)),cex.axis=1.5)
axis(side = 2,cex.axis=1.5)
#####################################################################################
##Hillsborough and ATMP
cmat[,,4] = getcontour(15, 3, meanvals, numrows, alldata, addto, HABsvm)

contour(x=undoscale(15, addto+meanvals[15], scalefactors),
        y=undoscale(3, addto+meanvals[3], scalefactors),
        z=cmat[,,4],
        xlab= "Hillsborough River (cfs)",
        ylab= "ATMP (celsius)",
        axes=F,
        cex.main=2,cex.lab=1.5,
        labcex=1.25)
axis(side = 1, at = axTicks(1), labels=round(10^axTicks(1)),cex.axis=1.5)
axis(side = 2,cex.axis=1.5)
dev.off()
