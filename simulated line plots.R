require(e1071)
require(neuralnet)
require(kernlab)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Read in data
load(file="radialSVM.RData")
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1

alldata = read.csv("alldata.csv")
alldata$X = NULL
for(i in 1:ncol(alldata)) alldata[is.na(alldata[,i]),i]=mean(alldata[,i],na.rm=T)
alldata[,10:23] = log10(alldata[,10:23]+4)
alldata$State = as.factor(alldata$State)
##########################################################################################################3
scalefactors = scale(alldata[,1:(ncol(alldata)-1)])
alldata[,1:(ncol(alldata)-1)] = scale(alldata[,1:(ncol(alldata)-1)])
months = as.numeric(format(HAB$Date,"%m"))
#####################################################################################################
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}

colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)

desiredmon = c(8:10)
meanvals = colMeans(alldata[!is.na(match(months,desiredmon)), 1:(ncol(alldata)-1)])
maxvals  = colMax(alldata[!is.na(match(months,desiredmon)), 1:(ncol(alldata)-1)])
minvals  = colMin(alldata[!is.na(match(months,desiredmon)), 1:(ncol(alldata)-1)])
####################################################################################
numrows = 1000
names = c("NS wind (venf1)", "EW wind (venf1)", 
          "ATMP (venf1)", 
          "NS wind (cdrf1)", "EW wind (cdrf1)",
          "NS wind (42039)", "EW wind (42039)",
          "NS wind (42003)", "EW wind (42003)",
          "Tampa Bay","Peace River","Lake Okeechobee","Suwanee River",
          "Withlacoochee River","Hillsborough River","Manatee River",
          "Myakka River")
tiff(file = "simulated_usgs.tiff", height = 9, width = 9,  units = "in",
     pointsize=10, res = 600, compression = c("lzw"))
par(mfrow=c(2,2),mai=c(0.75,0.75,0.75,0.75))
for(i in 18:23){#c(11, 13, 14, 15)
  curr = alldata[1,]
  curr[1:length(meanvals)] = meanvals

  preds = rep(NA, numrows)
  sequence = seq(from=minvals[i],
                 to=maxvals[i],
                 length.out=numrows)
  
  for(j in 1:numrows){
    #addto = rep(0, length(curr)-1)
    addto = curr[,1:(length(curr)-1)]
    addto[i] = sequence[j]
    preds[j] = attr(predict(HABsvm,curr[,1:(length(curr)-1)]+addto,
                            probability=T),"probabilities")[1,2]
  }
  plot(undoscale(i, curr[1,i]+sequence, scalefactors),
       preds, type="l",main=colnames(alldata)[i],#main=names[i],
       xlab = expression(paste('Discharge (m'^3,"s"^-1,")")),
       ylab = "HAB probability[%]",axes=F,
       ylim = c(0, 1),
       cex.main=3,cex.lab=2,
       lwd=4)
  axis(1,cex.axis=2,at=axTicks(side=1),labels=signif(10^axTicks(side=1)*0.0283168,digits=1))
  axis(2,cex.axis=2)

}
dev.off()

####################################################################################
bounds = c(-2,2)
meanvals = colMeans(alldata[, 1:(ncol(alldata)-1)])
print(meanvals)

tiff(file = "simulated_ndbc.tiff", width =9, height = 9, units = "in",
     pointsize=10, res = 600, compression = c("lzw"))
par(mfrow=c(2,2),mai=c(0.75,0.75,0.75,0.75))
for(i in c(1,2,4,5)){
  curr = alldata[1,]
  curr[1:length(meanvals)] = meanvals

  preds = rep(NA, numrows)
  sequence = seq(from=bounds[1],to=bounds[2],length.out=numrows)
  for(j in 1:numrows){
    addto = rep(0, length(curr)-1)
    addto[i] = sequence[j]
    preds[j] = attr(predict(HABsvm,curr[,1:(length(curr))]+addto,
                            probability=T),"probabilities")[1,2]
  }
  if(i!= 3){
    plot(undoscale(i, sequence+curr[1,i], scalefactors),preds, main=names[i],type="l",
         xlab = "Wind speed (m/s)",
         ylab = "HAB probability [%]", axes=F,
         cex.main=3,cex.lab=2,
         lwd=4)
    axis(1,cex.axis=2)
    axis(2,cex.axis=2)
  }else{
    plot(undoscale(i, sequence+curr[1,i], scalefactors),preds, main=names[i],type="l",
         xlab = "ATMP (Celsius)",
         ylab = "HAB event")

  }
}
dev.off()
#####################################################################################
