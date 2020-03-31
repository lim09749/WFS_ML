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

##Read in data

HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
x = is.na(HAB$Abundance_cells)
#####################################################################################
##Read in entire dataset
alldata = read.csv("alldata.csv")
alldata$X = NULL

scalefactors = scale(alldata)
alldata = as.data.frame(scale(alldata))

alldata$State = as.factor(states)
prevdata$State = as.factor(states)
alldata$State = as.numeric(as.character(alldata$State))

##Create previous values dataset
prevdata = read.csv("alldata.csv")
prevdata$X = NULL
prevdata$State = as.factor(states)

##Remove NA values
alldata = alldata[!x,]
HAB = HAB[!x,]
prevdata = prevdata[!x, ]

#####################################################################################

##SVM
svmdata = alldata[,c(1:9,10:17,31)]
svmdata$State = as.factor(alldata$State)
# svmdata = rbind(svmdata,
#                 svmdata[newdata$States==1,],
#                 svmdata[newdata$States==1,])
                #SMOTE(State~., svmdata, perc.over = 100,k=5,perc.under=0))
HABsvm = tune.svm(as.factor(State)~.,
                  data=svmdata,
                  type="C-classification",
                  cost = 2^(-5:10),
                  kernel="radial",
                  probability=T)
HABsvm = HABsvm$best.model
# save(HABsvm,file="dis_svm.RData")


##Naive Bayes
nbdata = alldata[,c(1:9,10:17,31)]
nbdata$State = as.factor(alldata$State)
HABnbayes = naiveBayes(State~.,data=nbdata)

##RVM
require(kernlab)
rvmdata = alldata[,c(1:9,10:17,31)]
rvmdata$State = as.numeric(alldata$State)
HABrvm = rvm(State~., data= rvmdata,
             kernel = "rbfdot")

##ANN
# require(neuralnet)
# HABnet = neuralnet(State~.,
#                    data=newdata,
#                    linear.output = F,
#                    hidden = c(20, 10))
#####################################################################################################
reliability.plot(as.numeric(alldata$State), svm_prob(HABsvm, alldata))
reliability.plot(as.numeric(alldata$State), rvm_prob(HABrvm, alldata))
reliability.plot(as.numeric(alldata$State), bay_prob(HABnbayes, alldata))

##################################################################################################################
##Columns of discharge
plot_col = c(11,13,15,17,31)

## Number of rows for values
numrows=1000
rightrows = rep(TRUE, nrow(prevdata))

##Find mean values
prevdata$State = NULL
meanvals = colMeans(prevdata)

##Get colors
colors = c("red","blue","green","cyan","purple")

##Values for axes ticks
axes_ticks = rbind(log10(c(0.5,1,2.5,5,25,50,100)),
                   log10(c(15,25,50,75,100,150,200,250,300,350,400,450)),
                   log10(c(0.5,1,2.5,5,10,15,20,25,30)),
                   log10(c(0.25,0.5,1,2.5,5,10,15,20,25,30)),
                   log10(c(2.5,5,10,25,50,100,150)))
par(mfrow=c(3,2))
for(j in 1:length(plot_col)){
  i = plot_col[j]
  
  ##Create range of values explored
  xseq = seq(from=quantile(prevdata[,i])[2]*0.9,
             to=mean(c(quantile(prevdata[,i])[4],quantile(prevdata[,i])[5])),
             length.out=numrows)
  
  ##Create simulated data
  dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
  dat[,i] = xseq
  
  ##Scale data
  dat = scale(dat, center = attr(scalefactors,"scaled:center"),
              scale = attr(scalefactors,"scaled:scale"))
  dat = as.data.frame(dat)
  colnames(dat) = colnames(scalefactors)
  
  ##Plot data
  ##SVM-red
  ##Lin. reg needed to adjust probabilities to same scale
  plot(log10(10^xseq/(60*60*24*7)*0.0283168), -0.001085+1.159881*svm_prob(HABsvm,dat), 
       main=colnames(alldata)[i],
       xlab=expression("Discharge ("~m^{3}/s~")"), 
       ylab="HAB probability [%]",ylim=c(0.1,1),
       col=colors[1],
       cex.main=3,
       cex.lab=1.5,
       xlim=c(max(c(-1,min(log10(10^xseq/(60*60*24*7)*0.0283168)))),max(log10(10^xseq/(60*60*24*7)*0.0283168))),
       axes=F)
  
  ##NB-blue
  ##Lin. reg needed to adjust probabilities to same scale
  points(log10(10^xseq/(60*60*24*7)*0.0283168), 0.1540+0.3415*bay_prob(HABnbayes, dat), col=colors[2], cex.main=3)
  
  ##RVM-green
  ##Lin. reg needed to adjust probabilities to same scale
  points(log10(10^xseq/(60*60*24*7)*0.0283168), -0.04728+1.13998*rvm_prob(HABrvm, dat), col=colors[3], cex.main=3)
  
  # ##ANN-cyan
  # points(log10(10^xseq/(60*60*24*7)*0.0283168), ann_prob(HABnet, dat), col=colors[4], cex.main=3)
  
  ##Create axes
  side1 = axes_ticks[j,]
  # axis(at=axTicks(side=1),labels=signif(10^axTicks(side=1),digits=2),side=1,cex.axis=1.75)
  axis(at=side1,labels=signif(10^side1,digits=2),side=1,cex.axis=1.5)
  axis(side=2,cex.axis=1.75)
  
  # ##Histogram of discharge
  # h = hist(10^prevdata[,i]/(60*60*24*7)*0.0283168, plot=F, nclass=20)
  # h$density = h$counts/sum(h$counts)*100
  # plot(h, freq=F, xlim = c(0, 100))
}
#beep(sound=3)