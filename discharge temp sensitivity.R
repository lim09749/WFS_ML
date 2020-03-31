require(beepr)
require(kernlab)
require(DMwR)
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
fixbounds = function(x){
  y=x
  y[x<=0]=0
  y[x>=1]=1
  y
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
getcontour = function(col1, col2,meanvals, numrows, alldata, addto1, addto2, model, type){
  cmat = matrix(NA, nrow=numrows, ncol=numrows)
  for(i in 1:nrow(cmat)){
    for(j in 1:ncol(cmat)){
      curr = as.data.frame(rbind(meanvals,meanvals))
      colnames(curr) = colnames(alldata)[1:(ncol(alldata)-1)]
      curr[,col1] = curr[,col1]+addto1[i]
      curr[,col2] = curr[,col2]+addto2[j]
      if(type == "RVM"){
        cmat[i,j] = rvm_prob(model,curr)[1]
      }else if (type=="SVM"){
        cmat[i,j] = svm_prob(model,curr)[1]

      }
    }
  }
  return(cmat)
}
#####################################################################################

##Read in entire dataset
##Read in totaldata
alldata = read.csv("alldata.csv")
alldata$X = NULL
x=is.na(HAB$Abundance_cells)
alldata = alldata[!x,]
HAB = HAB[!x,]

alldata = alldata[,c(3,11:17,31)]##discharge and temp data
prevdata = alldata
scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states[!x])

# HABsvm = tune.svm(as.factor(State)~.,
#                   data=rbind(alldata,
#                              alldata[alldata$State=="1",],
#                              alldata[alldata$State=="1",]),
#                   type="C-classification",
#                   cost = 2^(-5:10),
#                   kernel="radial",
#                   probability=T)
# svm_temp = HABsvm$best.model

rvm_temp = kernlab::rvm(as.numeric(as.character(State))~.,
                      data= rbind(alldata,
                                  alldata[alldata$State=="1",],
                                  alldata[alldata$State=="1",]
                                  #SMOTE(State~., alldata,perc.over=200)
                                  ),
                      kernel = "rbfdot")
# ,
# alldata[alldata$State=="1",],
# alldata[alldata$State=="1",]
#####################################################################################


numrows=100
xcol = which(colnames(alldata)=="Peace")
ycol = which(colnames(alldata)=="venf1.ATMP")
addto1 = seq(-1.25, 1.75, length.out=numrows)
addto2 = seq(-1.25, 1.75, length.out=numrows)
meanvals = colMeans(alldata[alldata$State==1,1:(ncol(alldata)-1)])#
cmat = getcontour(xcol, ycol, meanvals, numrows, alldata, addto1, addto2, rvm_temp, "RVM")
contour(x=log10(10^undoscale(xcol, addto1+meanvals[xcol], scalefactors)/(60*60*24*7)*0.0283168),
        y=undoscale(ycol, addto2+meanvals[ycol], scalefactors),
        z=cmat,
        xlab= colnames(alldata)[xcol],
        ylab= colnames(alldata)[ycol],
        axes=F,
        lwd=3,
        cex.main=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 1, at = log10(c(2.5,5,10,25,50,100,200)),
     labels=signif(c(2.5,5,10,25,50,100,200),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(axTicks(2),digits=2),cex.axis=1.5)
beep()
