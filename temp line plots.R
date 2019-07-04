require(e1071)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}

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
temps = alldata$venf1.ATMP
alldata = alldata[,1:9]
#alldata = alldata[,c(3,11:17,31)]##discharge and temp data
prevdata = alldata
scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)
alldata$venf1.ATMP[alldata$venf1.ATMP<=10]

# load("temp_svm.RData")
require(beepr)
HABsvm = tune.svm(as.factor(State)~.,
                  data=alldata,
                  type="C-classification",
                  cost = 2^(-5:10),
                  kernel="radial",
                  probability=T)
HABsvm = HABsvm$best.model
#save(HABsvm,file="temp_svm.RData")
beep(sound=3)

#####################################################################################
plot_col = c(3)
numrows=1000
rightrows = rep(TRUE, nrow(prevdata))
meanvals = colMeans(prevdata)
colors = c("red","blue","green","cyan","purple")
for(j in 1:length(plot_col)){
  i = plot_col[j]
  xseq = seq(from=quantile(prevdata[,i])[2],
             to=mean(c(quantile(prevdata[,i])[4],quantile(prevdata[,i])[5])),
             length.out=numrows)
  dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
  dat[,i] = xseq
  dat = scale(dat, center = attr(scalefactors,"scaled:center"),
              scale = attr(scalefactors,"scaled:scale"))
  dat = as.data.frame(dat)
  colnames(dat) = colnames(scalefactors)
  pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
  # if(j==1){
  plot(xseq,pred,main=colnames(alldata)[i],
       xlab=colnames(alldata)[plot_col[j]], ylab="HAB probability [%]",
       col=colors[j],cex.main=3,cex.lab=1.5)#,xlim=c(max(c(-1,min(log10(10^xseq/(60*60*24*7)*0.0283168)))),
                                           #        max(log10(10^xseq/(60*60*24*7)*0.0283168))),axes=F)
}
