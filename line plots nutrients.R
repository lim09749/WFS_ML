rm(list=ls())
require(e1071)
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Functions
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)
#####################################################################################
##Initialization of variables
numrows = 1000
desiredmon = c(1:12)
x = seq(0.5,1.5,length.out=numrows)
##Read in data
load(file="radialSVM.RData")
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
months = as.numeric(format(HAB$Date,"%m"))
whichright = !is.na(match(months, desiredmon))
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1

##Read in usgs and nutrient data
usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL

##Read in totaldata
alldata = read.csv("totaldata.csv")
alldata$X = NULL
for(i in 1:ncol(alldata)) alldata[is.na(alldata[,i]),i]=mean(alldata[,i],na.rm=T)
prevdata = alldata

alldata[,10:23] = log10(alldata[,10:23]+4)

meanvals = colMeans(alldata[whichright,1:23])
scalefactors = scale(alldata[,1:(ncol(alldata)-1)])
alldata[,1:23] = scale(alldata[,1:23])
stdmeanvals = as.matrix(scale(rbind(meanvals,meanvals),center=attr(scalefactors,"scaled:center"), scale=attr(scalefactors,"scaled:scale")))[1,]
###############################################################################################
##Adjust TN only
##Keep discharge constant and modify 
##Use this function to plot line
simulatedline = function(x, nsd, usgs, whichright, usgs_cols, rivers, suffix,
                         desiredmon, meanvals, numrows, TN_columns, scalefactors, 
                         HABsvm, alldata){
  meandischarge = (mean(usgs[whichright,usgs_cols[i]]*28.318)
                   +nsd*sd(usgs[whichright,usgs_cols[i]]*28.318))
  TN = read.csv(paste0(rivers[i],suffix))
  TN$SampleDate = as.Date(unlist(strsplit(as.character(TN$SampleDate),split = " "))[seq(1,nrow(TN)*2-1,by=2)],
                          format="%m/%d/%Y")
  print(length(which(TN$SampleDate>"1998-01-01" 
                     & !is.na(match(format(TN$SampleDate,"%m"),desiredmon)))))
  if(length(which(TN$SampleDate>"1998-01-01" 
                  & !is.na(match(format(TN$SampleDate,"%m"),desiredmon))))>0){
    TNvals = as.numeric(as.character(TN$Result_Value[TN$SampleDate>"1998-01-01"
                                                     & !is.na(match(format(TN$SampleDate,"%m"),desiredmon))]))*0.001
    simulated = mean(TNvals)*meandischarge*x
    dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
    dat[,TN_columns[i]] = log10(simulated+4)
    dat = scale(dat, center = attr(scalefactors,"scaled:center"),
                scale = attr(scalefactors,"scaled:scale"))
    dat = as.matrix(dat)
    pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
    return(pred)
  }
}
rivers = c("Tampa",
           "Myakka",
           "Peace",
           "Withlacoochee",
           "Manatee",
           "Hillsborough")
TN_columns = 12:17
usgs_cols  = c(1,8,2,5,7,6)
##Get percentages in units
par(mfrow=c(4,3))
for(i in 1:length(rivers)){
  plot(0,0,xlim=c(min(x)-0.1,max(x)+0.1),ylim=c(-0.1,1.1),main=colnames(alldata)[TN_columns[i]],
       xlab=desiredmon)
  
  pred1 = simulatedline(x, 0, usgs, whichright, usgs_cols, rivers, "_TN.csv",
                desiredmon, meanvals, numrows, TN_columns, scalefactors, 
                HABsvm, alldata)
  lines(x, pred1, col="green")
  abline(h=0.5)
}

##Adjust TP
TP_columns = 18:23
for(i in 1:length(rivers)){
  plot(0,0,xlim=c(min(x)-0.1,max(x)+0.1),ylim=c(-0.1,1.1),main=colnames(alldata)[TP_columns[i]])
  pred1 = simulatedline(x, 0, usgs, whichright, usgs_cols, rivers, "_TP.csv",
                desiredmon, meanvals, numrows, TP_columns, scalefactors, 
                HABsvm, alldata)
  lines(x, pred1, col="green")
  abline(h=0.5)
}

#####################################################################################################
##Adjust rivers discharge
for(i in 10:11){
  simulated = sapply(x,
                     FUN = function(x) mean(prevdata[whichright,i])*x)
  
  dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
  dat[,i] = log10(simulated+4)
  dat = scale(dat, center = attr(scalefactors,"scaled:center"),
              scale = attr(scalefactors,"scaled:scale"))
  dat = as.matrix(dat)
  pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
  plot(x,
       pred,main=colnames(alldata)[i],type="l",
       ylim=c(0,1),xlab="")
}

correspondingcols = cbind(12:17, 18:23)
for(i in 1:nrow(correspondingcols)){
  meandischarge = (mean(usgs[whichright,usgs_cols[i]]*28.318)
                   +sd(usgs[whichright,usgs_cols[i]]*28.318))
  TN = read.csv(paste0(rivers[i],"_TN.csv"))
  TN$SampleDate = as.Date(unlist(strsplit(as.character(TN$SampleDate),split = " "))[seq(1,nrow(TN)*2-1,by=2)],
                          format="%m/%d/%Y")
  TNvals = as.numeric(as.character(TN$Result_Value[TN$SampleDate>"1998-01-01"
                                                   & !is.na(match(format(TN$SampleDate,"%m"),desiredmon))]))*0.001
  TP = read.csv(paste0(rivers[i],"_TP.csv"))
  TP$SampleDate = as.Date(unlist(strsplit(as.character(TP$SampleDate),split = " "))[seq(1,nrow(TP)*2-1,by=2)],
                          format="%m/%d/%Y")
  TPvals = as.numeric(as.character(TP$Result_Value[TP$SampleDate>"1998-01-01"
                                                   & !is.na(match(format(TP$SampleDate,"%m"),desiredmon))]))*0.001
  
  simulatedTN = mean(TNvals)*meandischarge*x
  simulatedTP = mean(TPvals)*meandischarge*x
  dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
  dat[,correspondingcols[i,1]] = log10(simulatedTN+4)
  dat[,correspondingcols[i,2]] = log10(simulatedTP+4)
  dat = scale(dat, center = attr(scalefactors,"scaled:center"),
              scale = attr(scalefactors,"scaled:scale"))
  dat = as.matrix(dat)
  pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
  
  plot(x,
       pred,main=rivers[i],type="l",
       ylim=c(0,1),xlab="")
}

#####################################################################################################
##Adjust Wind speeds and temperature
bounds = c(-2,2)
meanvals = colMeans(alldata[,1:(ncol(alldata)-1)])

for(i in 1:9){
  curr = alldata[1,]
  curr[1:length(meanvals)] = meanvals
  
  preds = rep(NA, numrows)
  sequence = seq(from=bounds[1],to=bounds[2],length.out=numrows)
  for(j in 1:numrows){
    addto = rep(0, length(curr)-1)
    addto[i] = sequence[j]
    preds[j] = attr(predict(HABsvm,curr[,1:(length(curr)-1)]+addto,
                            probability=T),"probabilities")[1,2]
  }
  plot(undoscale(i, sequence+curr[1,i], scalefactors),preds,type="l",
       xlab = desiredmon,
       ylab = "HAB probability [%]", main=colnames(alldata)[i])
}
