require(e1071)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Get Data
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
months = as.numeric(format(HAB$Date,format="%m"))
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
#####################################################################################

##Read in entire dataset
##Read in totaldata
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(3,11:17,31)]##discharge and temp data
prevdata = alldata
scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)

# load("discharge_temp_svm.RData")
require(beepr)
HABsvm = tune.svm(as.factor(State)~.,
                  data=alldata,
                  type="C-classification",
                  cost = 2^(-5:10),
                  kernel="radial",
                  probability=T)
HABsvm = HABsvm$best.model
#save(HABsvm,file="discharge_temp_svm.RData")
beep(sound=3)

#####################################################################################
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}
getcontour = function(col1, col2,meanvals, numrows, alldata, addto1, addto2, HABsvm){
  cmat = matrix(NA, nrow=numrows, ncol=numrows)
  for(i in 1:nrow(cmat)){
    for(j in 1:ncol(cmat)){
      curr = as.data.frame(rbind(meanvals,meanvals))
      colnames(curr) = colnames(alldata)[1:(ncol(alldata)-1)]
      curr[,col1] = curr[,col1]+addto1[i]
      curr[,col2] = curr[,col2]+addto2[j]
      cmat[i,j] = attr(predict(HABsvm,curr,
                               probability=T),"probabilities")[1,2]
    }
  }
  return(cmat)
}
#####################################################################################
numrows=100
par(mfrow=c(2,2))
xcol = 6
ycol = 1
addto = seq(-1, 1.75, length.out=numrows)
meanvals = colMeans(alldata[,1:(ncol(alldata)-1)])#
cmat = getcontour(xcol, ycol, meanvals, numrows, alldata, addto, addto, HABsvm)
contour(x=log10(10^undoscale(xcol, addto+meanvals[xcol], scalefactors)/(60*60*24*7)*0.0283168),
        y=undoscale(ycol, addto+meanvals[ycol], scalefactors),
        z=cmat,
        xlab= colnames(alldata)[xcol],
        ylab= colnames(alldata)[ycol],
        axes=F,
        lwd=3,
        cex.main=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 2, at = axTicks(2), labels=signif(axTicks(2),digits=2),cex.axis=1.5)
axis(side = 1, at = log10(c(2.5,5,10,50,100)), 
     labels=signif(c(2.5,5,10,50,100),digits=2),cex.axis=1.5)
