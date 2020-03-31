rm(list=ls())
require(e1071)
require(kernlab)
require(neuralnet)
require(beepr)
require(DMwR)
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Create functions
fixbounds = function(x){
  y=x
  y[x<=0]=0
  y[x>=1]=1
  y
}
ann_prob = function(HABnet, data){
  fixbounds(predict(HABnet, data))[,1]
}
bay_prob = function(HABnbayes, data){
  fixbounds(predict(HABnbayes, data, "raw")[,2])
}
rvm_prob = function(HABrvm, data){
  fixbounds(predict(HABrvm, data))[,1]
}
svm_prob = function(HABsvm, data){
  fixbounds(attr(predict(HABsvm, data, probability=T),"probabilities")[,2]) 
}
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}
##############################################################################################################

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
alldata = alldata[,c(3,11:17,31)]##wind data
prevdata = alldata
scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)



##Remove NA values
x=is.na(HAB$Abundance_cells)
alldata = alldata[!x,]
HAB = HAB[!x,]

numrows=1000
rightrows = rep(TRUE, nrow(prevdata))
meanvals = colMeans(alldata[,-ncol(alldata)])


########################################################################################
##Train data
# 
# HABsvm = tune.svm(as.factor(State)~.,
#                   data=rbind(alldata,alldata[alldata$State=="1",],alldata[alldata$State=="1",]),
#                              #SMOTE(State~., alldata, perc.over = 100,k=5,perc.under=0)),
#                   type="C-classification",
#                   cost = 2^(-5:10),
#                   kernel="radial",
#                   probability=T)
# svm_wind = HABsvm$best.model
# save(svm_wind,file="new_svm_wind.RData")
# 
# ##Naive Bayes
# bay_wind = naiveBayes(State~.,
#                        data=rbind(alldata,
#                                   alldata[alldata$State=="1",],
#                                   alldata[alldata$State=="1",]))#rbind(alldata,SMOTE(State~., alldata, perc.over = 200,k=5,perc.under=0)))
# save(bay_wind, file = "new_bay_wind.RData")
# 
# ##RVM
# rvm_wind = kernlab::rvm(as.numeric(as.character(State))~.,
#                       data= rbind(alldata,
#                                   alldata[alldata$State=="1",],
#                                   alldata[alldata$State=="1",]),
#                       #SMOTE(State~., alldata, perc.over = 300,k=5,perc.under=0)
#                                   # alldata[alldata$State=="1",]),#rbind(alldata,alldata[alldata$State=="1",]),
#                       kernel = "rbfdot")
# save(rvm_wind, file = "new_rvm_wind.RData")
# beep()
######################################################################################################
load("new_svm_wind.RData")
load("new_bay_wind.RData")
load("new_rvm_wind.RData")


tiff(file = "wind_line_fig_4.tiff", height = 4.5, width = 9,  units = "in",
     pointsize=10, res = 600, compression = c("lzw"))
par(mfrow=c(1,2))
for(i in 1:2){
  xseq = seq(from = -2,
             to   = 2,
             length.out=numrows)
  dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
  dat[,i] = xseq
  # dat = scale(dat, center = attr(scalefactors,"scaled:center"),
  #             scale = attr(scalefactors,"scaled:scale"))
  dat = as.data.frame(dat)
  colnames(dat) = colnames(scalefactors)
  
  x = undoscale(i,xseq,scalefactors)
  ##Print summary
  ##Linear regression adjusts probability value, 
    ## done with lm(seq(0,1, length.out = nrow(alldata))~sort(rvm_prob(rvm_wind, alldata)))
  print(summary(rvm_prob(rvm_wind, dat)*0.8077-0.145))
  plot(xseq, rvm_prob(rvm_wind, dat)*0.8077-0.1459,type="l",
       xlab="Normalized Wind",main="", ylim=c(0.2,0.6),
       ylab="",cex.lab=1.5,
       col="black", lwd=2.5, lty="solid",
       cex.axis=1.5)
  # plot(x, svm_prob(svm_wind, dat),type="l",
  #      xlab="Normalized Wind",axes=F,main="", ylim=c(0.1,1),
  #      ylab="HAB probability",cex.lab=1.5,
  #      col="red", lwd=2.5, lty="dashed")
  # lines(x, rvm_prob(rvm_wind, dat)*1.16-0.1391,
  #       lty="solid",lwd=4,col="green")
  # lines(x, bay_prob(bay_wind, dat),
  #       lty="dotted",lwd=2.5,col="blue")
  
  ##lines(x, ann_prob(ann_wind, dat))
  abline(v=0,col="black", lty="dashed",lwd=2)
  if(i==1){
    mtext(side=2,"HAB probability", line =2.5,cex=1.5)
    text(-0.25,0.21,"N",cex=2.5)
    text(0.25,0.21,"S",cex=2.5)
    
  }else{
    text(-0.25,0.21,"E",cex=2.5)
    text(0.25,0.21,"W",cex=2.5)
  }
}

dev.off()

