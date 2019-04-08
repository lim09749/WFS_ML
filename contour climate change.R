rm(list=ls())
require(e1071)
####################################################################################################
setwd("C:/Users/xswang/HAB Research")

tempmin = seq(from=0, to=1, length.out=9)
tempavg = seq(from=0, to=2.5, length.out=9)
tempmax = seq(from=0, to=6, length.out=9)

prmin   = seq(from=0, to=5, length.out=9)
pravg   = seq(from=0, to=10, length.out=9)
prmax   = seq(from=0, to=20, length.out=9)
####################################################################################################
##Get HAB

##Form alldata
load(file="radialSVM.RData")
ndbc = read.csv("NDBC_station_data.csv")
ndbc = ndbc[,!is.na(colSums(ndbc))]
ndbc$X = NULL
usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL
usgs[,1:(ncol(usgs)-1)] = log10(usgs[,1:(ncol(usgs)-1)]+4)

HAB = read.csv("C:/Users/xswang/HAB Research/West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Abundance_cells[is.na(HAB$Abundance_cells)] = 0
HAB$Date = as.Date(HAB$Date)
months = as.numeric(format(HAB$Date,"%m"))

#####################################################################################
alldata = cbind(ndbc,usgs)
for(i in 1:ncol(alldata)) alldata[,i] = as.numeric(as.character(alldata[,i]))
alldata$State = as.factor(alldata$State)
scalefactors = scale(alldata[,1:(ncol(alldata)-1)])
alldata[,1:(ncol(alldata)-1)] = scale(alldata[,1:(ncol(alldata)-1)])

undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}
#####################################################################################
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
tiff(file = "future contour.tiff", width =15, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(1,2),mar=c(6, 5.5, 5, 3) + 0.1)
numrows = 100
addto = seq(-2, 2, length.out=numrows)
meanvals = colMeans(alldata[!is.na(match(months,8:10)), 1:(ncol(alldata)-1)])
print(meanvals)
#####################################################################################
##Hillsborough and Peace
cmat = getcontour(15, 11, meanvals, numrows, alldata, addto, HABsvm)
contour(x=undoscale(15, addto+meanvals[15], scalefactors)+log10(0.0283168),
        y=undoscale(11, addto+meanvals[11], scalefactors)+log10(0.0283168),
        z=cmat,
        xlab= expression(paste('Hillsborough River (m'^3,"s"^-1,")")),
        ylab= expression(paste('Peace River (m'^3,"s"^-1,")")),
        axes=F,
        cex.main=3,cex.lab=2,
        labcex=1.25,
        levels=c(0.2,0.4,0.6,0.8,0.1))
axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)
points(undoscale(15, 0.5+log10(1+pravg), scalefactors)[c(1,5,9)]+log10(0.0283168), 
       undoscale(11, 0.5+log10(1+pravg), scalefactors)[c(1,5,9)]+log10(0.0283168),
       col="blue",pch=19)
text(undoscale(15, 0.5+log10(1+pravg), scalefactors)[c(1,5,9)]+log10(0.0283168), 
     undoscale(11, 0.5+log10(1+pravg), scalefactors)[c(1,5,9)]+log10(0.0283168),
     labels = c("2020","2060","2100"),cex=1.75,
     col="red")
legend(x=log10(3.2), y= log10(5), 
       legend=c("RCP4.5"), 
       col=c("blue"), 
       lty=0, lwd=3, pch=c(19),cex=2)
#####################################################################################

##Peace and ATMP
cmat = getcontour(11, 3, meanvals, numrows, alldata, addto, HABsvm)

contour(x=undoscale(11, addto+meanvals[11], scalefactors)+log10(0.0283168),
        y=undoscale(3, addto+meanvals[3], scalefactors),
        xlim=c(log10(20),2.75),
        z=cmat,
        xlab= expression(paste('Peace River (m'^3,"s"^-1,")")),
        ylab= "Temperature (celsius)",
        axes=F,
        labcex=1.25,
        cex.main=3,cex.lab=2,
        levels=c(0.1,0.2,0.3,0.4))
axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=3),cex.axis=1.5)
axis(side = 2, at = axTicks(2), cex.axis=1.5)
points(undoscale(11,0.5+log10(1+pravg), scalefactors)[c(1,5,9)]+log10(0.0283168), 
       (undoscale(3,meanvals[3],scalefactors)+tempavg)[c(1,5,9)],
       col="blue",pch=19)
text(undoscale(11, 0.5+log10(1+pravg), scalefactors)[c(1,5,9)]+log10(0.0283168), 
       (undoscale(3,meanvals[3],scalefactors)+tempavg)[c(1,5,9)],
       labels = c("2020","2060","2100"),cex=1.75,
       col="red")
dev.off()
