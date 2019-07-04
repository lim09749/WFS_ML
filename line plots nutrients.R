rm(list=ls())
require(e1071)
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Functions
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}
colMax <- function(data) apply(data, 2, max, na.rm = TRUE)
colMin <- function(data) apply(data, 2, min, na.rm = TRUE)
#####################################################################################
##Initialization of variables
numrows = 1000
desiredmon = c(1:2,9:12)

##Span of percentages
nmat = rbind(seq(0.5,3,length.out=numrows),
             seq(0.5,2,length.out=numrows),
             seq(0.5,2,length.out=numrows),
             seq(0.5,1.5,length.out=numrows),
             seq(0.5,2,length.out=numrows),
             seq(0.5,2,length.out=numrows),
             seq(0.5,1.5,length.out=numrows))

pmat = rbind(seq(0.5,2.5,length.out=numrows),
             seq(0.5,2.5,length.out=numrows),
             seq(0.5,2,length.out=numrows),
             seq(0.5,2,length.out=numrows),
             seq(0.5,4,length.out=numrows),
             seq(0.5,2,length.out=numrows),
             seq(0.5,4,length.out=numrows))

##Read in data
caloosahatchee = read.csv("caloosahatchee_data.csv")
caloosahatchee$X = NULL

load(file="radialSVM.RData")
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1

months = as.numeric(format(HAB$Date,"%m"))
whichright = !is.na(match(months, desiredmon))

##Read in usgs and nutrient data
usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL
usgs = cbind(usgs, caloosahatchee$discharge)

##Read in totaldata
alldata = read.csv("alldata.csv")
alldata$X = NULL
prevdata = alldata

meanvals = rep(NA, ncol(alldata))
for(i in 1:ncol(alldata)){
  meanvals[i]=median(alldata[whichright,i])
}
minvals = colMin(alldata[whichright,])
maxvals = colMax(alldata[whichright,])
scalefactors = scale(alldata)
alldata = scale(alldata)
stdmeanvals = colMeans(alldata[whichright,])
###############################################################################################
##Find percentage span
perspan = function(suffix, i){
  TN = read.csv(paste0(rivers[i],suffix))
  TN$SampleDate = as.Date(unlist(strsplit(as.character(TN$SampleDate),split = " "))[seq(1,nrow(TN)*2-1,by=2)],
                          format="%m/%d/%Y")
  TN = TN[TN$SampleDate>"1998-01-01",]
  return(c(min(TN$Result_Value[TN$Result_Value!=0],na.rm=T),max(TN$Result_Value,na.rm=T))/mean(TN$Result_Value,na.rm=T))
}

##Use this function to plot line
simulatedline = function(x, i, nsd, usgs, whichright, usgs_cols, rivers, suffix,
                         desiredmon, meanvals, numrows, TN_columns, scalefactors, 
                         HABsvm, alldata){
  meandischarge = (median(usgs[whichright,usgs_cols[i]]*28.318)
                   +nsd*sd(usgs[whichright,usgs_cols[i]]*28.318))
  TN = read.csv(paste0(rivers[i],suffix))
  TN = TN[,!is.na(match(colnames(TN),c("SampleDate","Result_Value")))]
  TN = na.omit(TN)
  TN$SampleDate = as.Date(unlist(strsplit(as.character(TN$SampleDate),split = " "))[seq(1,nrow(TN)*2-1,by=2)],
                          format="%m/%d/%Y")
  print(length(which(TN$SampleDate>"1998-01-01"
                     & !is.na(match(format(TN$SampleDate,"%m"),desiredmon)))))
  if(length(which(TN$SampleDate>"1998-01-01" 
                  & !is.na(match(format(TN$SampleDate,"%m"),desiredmon))))>0){
    TNvals = as.numeric(as.character(TN$Result_Value[TN$SampleDate>"1998-01-01"
                                                     & !is.na(match(format(TN$SampleDate,"%m"),desiredmon))]))*0.001
    simulated = mean(TNvals)*meandischarge*x
    print(mean(TNvals))
    print(meandischarge)
    dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
    dat[,TN_columns[i]] = log10(simulated+4)
    dat = scale(dat, center = attr(scalefactors,"scaled:center"),
                scale = attr(scalefactors,"scaled:scale"))
    dat = as.data.frame(dat)
    colnames(dat) = colnames(scalefactors)
    pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
    plot(log10(simulated+4),pred,main=colnames(alldata)[TN_columns[i]])
    return(pred)
  }
}
rivers = c("Tampa",
           "Myakka",
           "Peace",
           "Withlacoochee",
           "Manatee",
           "Hillsborough",
           "Caloosahatchee")
###############################################################################################

##Get percentages in units
##Adjust TN
TN_columns = c(18:23,32)
usgs_cols  = c(1,8,2,5,7,6,10)
TP_columns = c(24:29,33)
# 
# par(mfrow=c(4,3))
# for(i in 1:length(rivers)){
#   x=nmat[i,]
#   #plot(0,0,xlim=c(min(x)-0.1,max(x)+0.1),ylim=c(-0.1,1.1),main=colnames(alldata)[TN_columns[i]])
#   print(perspan("_TN.csv",i))
#   print(colnames(alldata)[TN_columns[i]])
#   pred1 = simulatedline(seq(from=perspan("_TN.csv",i)[1],
#                             to=perspan("_TN.csv",i)[2],
#                             length.out=numrows), 
#                         i, 2, usgs, whichright, usgs_cols, rivers, "_TN.csv",
#                         desiredmon, meanvals, numrows, TN_columns, scalefactors, 
#                         HABsvm, prevdata)
#   x=log10(seq(from=perspan("_TN.csv",i)[1],
#               to=perspan("_TN.csv",i)[2],
#               length.out=numrows)*median(usgs[whichright,usgs_cols[i]])+4)
#   
#   multisimulatedlines(i,usgs,whichright,usgs_cols,rivers,suffix="_TN.csv",desiredmon,
#                       meanvals,numrows,TN_columns,scalefactors,HABsvm,alldata)
#   
#   # plot(0,0,xlim=c(min(x)-0.1,max(x)+0.1),ylim=c(-0.1,1.1),main=colnames(alldata)[TN_columns[i]])
#   # 
#   # lines(log10(seq(from=perspan("_TN.csv",i)[1],
#   #                 to=perspan("_TN.csv",i)[2],
#   #                 length.out=numrows)*median(usgs[whichright,usgs_cols[i]])+4), 
#   #       pred1, col="green")
#   
#   abline(h=0.5)
# }
# 
# ##Adjust TP
# par(mfrow=c(4,3))
# for(i in 1:length(rivers)){
#   x=pmat[i,]
#   #plot(0,0,xlim=c(min(x)-0.1,max(x)+0.1),ylim=c(-0.1,1.1),main=colnames(alldata)[TP_columns[i]])
#   print(perspan("_TP.csv",i))
# 
#   pred1 = simulatedline(seq(from=perspan("_TP.csv",i)[1],
#                             to=perspan("_TP.csv",i)[2],
#                             length.out=numrows), 
#                         i, 0, usgs, whichright, usgs_cols, rivers, "_TP.csv",
#                 desiredmon, meanvals, numrows, TP_columns, scalefactors, 
#                 HABsvm, alldata)
#   multisimulatedlines(i,usgs,whichright,usgs_cols,rivers,suffix="_TP.csv",desiredmon,
#                       meanvals,numrows,TP_columns,scalefactors,HABsvm,alldata)
#   
#   x=log10(seq(from=perspan("_TP.csv",i)[1],
#               to=perspan("_TP.csv",i)[2],
#               length.out=numrows)*median(usgs[whichright,usgs_cols[i]])+4)
#   
#   # plot(0,0,xlim=c(min(x)-0.1,max(x)+0.1),ylim=c(-0.1,1.1),main=colnames(alldata)[TP_columns[i]])
#   # 
#   # lines(x, pred1, col="green")
#   abline(h=0.5)
# }
# 
# #####################################################################################################
# ##Adjust rivers discharge
# par(mfrow=c(4,3))
# for(i in 13:13){
#   x=seq(quantile(usgs[whichright,4])[1]/median(usgs[whichright,4]),
#         max(usgs[whichright,4])/median(usgs[whichright,4]),length.out=numrows)
#   simulated = sapply(x,
#                      FUN = function(x) median(usgs[whichright,4])*x)
#   dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
#   dat[,i] = log10(simulated+4)
#   dat = scale(dat, center = attr(scalefactors,"scaled:center"),
#               scale = attr(scalefactors,"scaled:scale"))
#   dat = as.matrix(dat)
#   pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
#   plot(x,
#        pred,main=colnames(alldata)[i],type="l",
#        ylim=c(0,1),xlab="")
# }
# 
# correspondingcols = cbind(c(10,17,11,14,16,15,31),c(18:23,32), c(24:29,33))
# par(mfrow=c(4,3))
# for(i in 1:nrow(correspondingcols)){
#   x=seq(quantile(usgs[whichright,usgs_cols[i]])[1]/median(usgs[whichright,usgs_cols[i]]),
#         max(usgs[whichright,usgs_cols[i]])/median(usgs[whichright,usgs_cols[i]]),length.out=numrows)
#   meandischarge = median(usgs[whichright,usgs_cols[i]]*28.318)
#   TN = read.csv(paste0(rivers[i],"_TN.csv"))
#   TN = TN[,!is.na(match(colnames(TN),c("SampleDate","Result_Value")))]
#   TN = na.omit(TN)
#   TN$SampleDate = as.Date(unlist(strsplit(as.character(TN$SampleDate),split = " "))[seq(1,nrow(TN)*2-1,by=2)],
#                           format="%m/%d/%Y")
#   TNvals = as.numeric(as.character(TN$Result_Value[TN$SampleDate>"1998-01-01"
#                                                    & !is.na(match(format(TN$SampleDate,"%m"),desiredmon))]))*0.001
#   TP = read.csv(paste0(rivers[i],"_TP.csv"))
#   TP = TP[,!is.na(match(colnames(TP),c("SampleDate","Result_Value")))]
#   TP = na.omit(TP)
#   TP$SampleDate = as.Date(unlist(strsplit(as.character(TP$SampleDate),split = " "))[seq(1,nrow(TP)*2-1,by=2)],
#                           format="%m/%d/%Y")
#   TPvals = as.numeric(as.character(TP$Result_Value[TP$SampleDate>"1998-01-01"
#                                                    & !is.na(match(format(TP$SampleDate,"%m"),desiredmon))]))*0.001
#   
#   simulated_d = median(usgs[whichright,usgs_cols[i]])*x
#   print("##########################")
#   print(summary(median(usgs[whichright,usgs_cols[i]])*x))
#   print(summary(usgs[whichright,usgs_cols[i]]))
#   simulatedTN = median(TNvals,na.rm=T)*meandischarge*x
#   simulatedTP = median(TPvals,na.rm=T)*meandischarge*x
#   dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
#   dat[,correspondingcols[i,1]] = log10(simulated_d+4)
#   dat[,correspondingcols[i,2]] = log10(simulatedTN+4)
#   dat[,correspondingcols[i,3]] = log10(simulatedTP+4)
#   dat = scale(dat, center = attr(scalefactors,"scaled:center"),
#               scale = attr(scalefactors,"scaled:scale"))
#   dat = as.matrix(dat)
#   
#   pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
#   
#   plot(x,
#        pred,main=rivers[i],type="l",
#        ylim=c(0,1),xlab="")
# }
# 
# #####################################################################################################
# ##Adjust Wind speeds and temperature
# bounds = c(-2,2)
# for(i in c(1:9,30)){
#   curr = stdmeanvals
#   preds = rep(NA, numrows)
#   sequence = seq(from=bounds[1],to=bounds[2],length.out=numrows)
#   for(j in 1:numrows){
#     addto = rep(0, length(curr))
#     addto[i] = sequence[j]
#     preds[j] = attr(predict(HABsvm,rbind(curr+addto),
#                             probability=T),"probabilities")[1,2]
#   }
#   plot(undoscale(i, sequence, scalefactors),preds,type="l",
#        xlab = desiredmon,
#        ylab = "HAB probability [%]", main=colnames(alldata)[i])
# }
# ##############################################################################3
par(mfrow=c(3,4))
for(i in 1:(ncol(alldata))){
  xseq = seq(from = min(prevdata[,i]),
             to = max(prevdata[,i]),length.out=numrows)
  dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
  dat[,i] = xseq
  dat = scale(dat, center = attr(scalefactors,"scaled:center"),
              scale = attr(scalefactors,"scaled:scale"))
  dat = as.data.frame(dat)
  colnames(dat) = colnames(scalefactors)
  pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
  plot(xseq,pred,main=c(colnames(alldata)[i],":",i))
}
