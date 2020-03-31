rm(list=ls())
dev.off()
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
##Load models
load("new_svm_dis.RData")
load("new_rvm_dis.RData")
load("new_bay_dis.RData")
load("new_ann_dis.RData")

load("new_svm_tn.RData")
load("new_rvm_tn.RData")
load("new_bay_tn.RData")
load("new_ann_tn.RData")

load("new_svm_tp.RData")
load("new_rvm_tp.RData")
load("new_bay_tp.RData")
load("new_ann_tp.RData")
#######################################################################################################
##Read in entire dataset

HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
x = is.na(HAB$Abundance_cells)


alldata = read.csv("alldata.csv")
alldata$X = NULL

scalefactors = scale(alldata)
alldata = as.data.frame(scale(alldata))

alldata$State = as.factor(states)
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
##Adjust for ann
iso1 = isoreg(sort(ann_prob(ann_dis,alldata))[sort(ann_prob(ann_dis,alldata))<0.5],
            seq(0,0.5,length.out=length(which(sort(ann_prob(ann_dis,alldata))<0.5))))
iso2 = isoreg(sort(ann_prob(ann_dis,alldata))[sort(ann_prob(ann_dis,alldata))>0.5],
            seq(0,0.5,length.out=length(which(sort(ann_prob(ann_dis,alldata))>0.5))))
iso = isoreg(sort(ann_prob(ann_dis,anndata[(0.9*nrow(anndata)):nrow(anndata),])),
             seq(0,1,length.out=length((0.9*nrow(anndata)):nrow(anndata))))
#####################################################################################################
##Create fig5 
tiff("fig5_new.tiff", width = 3150, height = 1446, units = "px",
     pointsize=10, res = 300, compression = "none")
layout(matrix(1:15, byrow=T,nrow=3,ncol=5))
par(mai=c(0.37,0.37,0.2,0.2))
#############################################################################################3
##Columns of discharge

alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(1:9,10:17,31)]##river discharge
prevdata = alldata
scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)

plot_col = c(13, 15, 17, 11, 18)

## Number of rows for values
numrows=1000
rightrows = rep(TRUE, nrow(prevdata))

##Find mean values
prevdata$State = NULL
meanvals = colMeans(prevdata)

##Values for axes ticks
axes_ticks = rbind(log10(c(15,25,50,75,100,150,200,250,300,350,400,450)),
                   log10(c(0.5,1,2.5,5,10,15,20,25,30)),
                   log10(c(0.25,0.5,1,2.5,5,10,15,20,25,30)),
                   log10(c(0.5,1,2.5,5,25,50,100)),
                   log10(c(2.5,5,10,25,50,100,150)))

ylim=c(0.1, 0.8)
conv_val = function(x){
  mid = which.min(abs(iso$x-x))
  return(iso$y[mid])
  # if(x<0.05){
  #   mid = which.min(abs(iso1$x-x))
  #   return(mean(iso1$y[max(1,mid-500):min(length(iso1$x),mid+500)]))
  #   
  # }else if (x>0.95){
  #   mid = which.min(abs(iso2$x-x))
  #   return(mean(iso2$y[max(1,mid-500):min(length(iso2$x),mid+500)]))
  # }
  # return(x)
}
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

  plot(log10(10^xseq/(60*60*24*7)*0.0283168), svm_prob(svm_dis, dat),#predict(svm_dis_lm,k), 
       main=colnames(alldata)[i],
       type="l",
       lty = "dashed",
       col="red",
       xlab="",
       ylab="",
       ylim=ylim,
       lwd=2.5,
       cex.main=1.5,
       cex.lab=1,
       xlim=c(max(c(-1,min(log10(10^xseq/(60*60*24*7)*0.0283168)))),max(log10(10^xseq/(60*60*24*7)*0.0283168))),
       axes=F)
  
  ##NB-blue
  ##Lin. reg needed to adjust probabilities to same scale
  lines(log10(10^xseq/(60*60*24*7)*0.0283168), 
        lty="dotted",lwd=2.5,col="blue",
        bay_prob(bay_dis, dat))
  abline(h=0.5)
  ##RVM-green
  ##Lin. reg needed to adjust probabilities to same scale
  lines(log10(10^xseq/(60*60*24*7)*0.0283168), 
        lty="solid",lwd=4,col="green",
        rvm_prob(rvm_dis, dat))
  print(c(colnames(alldata)[i],rvm_prob(rvm_dis, dat)[1],rvm_prob(rvm_dis, dat)[nrow(dat)]))
  # ##ANN-cyan
  a_p = ann_prob(ann_dis, dat)
  
  # probs = sapply(a_p, FUN = conv_val )

  # lines(log10(10^xseq/(60*60*24*7)*0.0283168),
  #       probs,
  #       lty="dotdash",lwd=2.5)
  
  ##Create axes
  side1 = axes_ticks[j,]
  ylab = "HAB probabililty"
  if(j>=2) ylab= "" 
  
  title(xlab="",
        ylab=ylab,
        cex.lab = 1.5,
        line = 2.25)
  
  title(xlab=expression("Discharge ("~m^{3}/s~")"),
        ylab="",
        cex.lab = 1.5,
        line = 2.65)
  axis(at=side1,labels=signif(10^side1,digits=2),side=1,cex.axis=1.25)
  axis(side=2,cex.axis=1.25,at=c(0.1,0.3,0.5,0.7,0.9), labels=c(0.1,0.3,0.5,0.7,0.9))
}
plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10),
     axes=F, ann=F)

####################################################################################################
##TN

##Get Data
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(1:9,18:23,32)]
alldata = alldata[alldata$Tampa_TN_mg>6 & alldata$Caloosahatchee_TN_mg>5,]
prevdata = alldata

scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)

plot_col = c(7,8,11,12)+4 ##TN
numrows=1000
rightrows = !is.na(match(format(HAB$Date,"%m"),c(9:12,1:2)))
meanvals = colMeans(prevdata)

##Get axes
axes_ticks = rbind(log10(c(0.1,0.25,0.5,1,2.5,5)),
                   log10(c(1,2.5,5,10,11)),
                   log10(c(0.1,0.25,0.5,1,2.5)),
                   log10(c(2,5,10,20)))
for(j in 1:length(plot_col)){
  i = plot_col[j]
  xseq = seq(from=log10(quantile(10^prevdata[,i])[2]),
             to=(quantile(prevdata[,i])[4]+quantile(prevdata[,i])[5])/2,
             length.out=numrows)
  dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
  dat[,i] = xseq
  dat = scale(dat, center = attr(scalefactors,"scaled:center"),
              scale = attr(scalefactors,"scaled:scale"))
  dat = as.data.frame(dat)
  colnames(dat) = colnames(scalefactors)
  
  
  plot(log10(10^xseq/10^9/7), svm_prob(svm_tn, dat),
       main="",
       type="l",
       lty = "dashed",
       col="red",
       xlab="",
       ylab="",
       ylim=c(0.1,0.9),
       lwd=2.5,
       cex.main=1.5,
       cex.lab=1,
       axes=F)
  
  ##NB-blue
  ##Lin. reg needed to adjust probabilities to same scale
  lines(log10(10^xseq/10^9/7), 
        lty="dotted",lwd=2.5,col="blue",
        bay_prob(bay_tn, dat))
  
  ##RVM-green
  ##Lin. reg needed to adjust probabilities to same scale
  lines(log10(10^xseq/10^9/7), 
        lty="solid",lwd=4,col="green",
        rvm_prob(rvm_tn, dat))
  
  # ##ANN-cyan
  # a_p = ann_prob(ann_dis, dat)
  # conv_val = function(x){
  #   mid = which.min(abs(iso$x-x))
  #   s = seq(from=max(1,mid),
  #           to= min(length(iso$y), mid),
  #           by=1)
  #   return(mean(iso$y[s]))
  # }
  # probs = sapply(a_p, FUN = conv_val )
  #   
  # lines(log10(10^xseq/10^9/7),
  #       a_p, 
  #       lty="dotdash",lwd=2.5)
  
  
  
  ##Axes
  side1 = axes_ticks[j,]
  axis(at=axTicks(side=1),labels=signif(10^axTicks(side=1),digits=2),side=1,cex.axis=1.25)
  axis(side=2,cex.axis=1.25,at=c(0.1,0.3,0.5,0.7,0.9),labels=c("0.1","0.3","0.5","0.7","0.9"))
  
  ylab = "HAB probability"
  if(j>=2) ylab= "" 
  title(xlab="TN (tons/day)",
        ylab=ylab,
        cex.lab = 1.5,
        line=2.25)
  
}
plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10),
     axes=F, ann=F)
legend(0, 9, legend=c("SVM", "RVM", "NB"),
       col=c("red", "green", "blue"), 
       lty=c("dashed","solid","dotted"), cex=2, lwd=2.5)

####################################################################################################
##TP
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(1:9,24:29,33)]
prevdata = alldata

scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))

plot_col = c(7,8,11,12)+4 ##TP
meanvals = colMeans(prevdata)
axes_ticks = rbind(log10(c(0.05,0.05,0.5,0.1,0.25,0.5)),
                   log10(c(0.025,0.05,0.1,0.25,0.5,1,1.5)),
                   log10(c(0.25,0.5,0.75,1,2,3,4,5,6)),
                   log10(c(0.15,0.25,0.5,1,1.5,2,2.5)))

for(j in 1:length(plot_col)){
  i = plot_col[j]
  xseq = seq(from=max(c(quantile(prevdata[,i])[1],6.5+log10(10))),
             to=mean(c(quantile(prevdata[,i])[4],quantile(prevdata[,i])[5])),
             length.out=numrows)
  if(j==2){
    xseq = seq(from=max(c(quantile(prevdata[,i])[1],6.5+log10(10))),
               to=quantile(prevdata[,i])[4],
               length.out=numrows)
  }
  dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
  dat[,i] = xseq
  dat = scale(dat, center = attr(scalefactors,"scaled:center"),
              scale = attr(scalefactors,"scaled:scale"))
  dat = as.data.frame(dat)
  colnames(dat) = colnames(scalefactors)
  
  
  plot(log10(10^xseq/10^9/7), svm_prob(svm_tp, dat),
       main="",
       type="l",
       lty = "dashed",
       col="red",
       xlab="",
       ylab="",
       ylim=c(0.1,0.9),
       lwd=2.5,
       cex.main=1.5,
       cex.lab=1,
       axes=F)
  
  ##NB-blue
  ##Lin. reg needed to adjust probabilities to same scale
  lines(log10(10^xseq/10^9/7), 
        lty="dotted",lwd=2.5,col="blue",
        bay_prob(bay_tp, dat))
  
  ##RVM-green
  ##Lin. reg needed to adjust probabilities to same scale
  lines(log10(10^xseq/10^9/7), 
        lty="solid",lwd=4,col="green",
        rvm_prob(rvm_tp, dat))
  
  # ##ANN-cyan
  # a_p = ann_prob(ann_dis, dat)
  # conv_val = function(x){
  #   mid = which.min(abs(iso$x-x))
  #   s = seq(from=max(1,mid),
  #           to= min(length(iso$y), mid),
  #           by=1)
  #   return(mean(iso$y[s]))
  # }
  # probs = sapply(a_p, FUN = conv_val )
  #   
  # lines(log10(10^xseq/10^9/7),
  #       a_p, 
  #       lty="dotdash",lwd=2.5)
  
  ##Axes
  side1 = axes_ticks[j,]
  axis(at=axTicks(side=1),labels=signif(10^axTicks(side=1),digits=2),side=1,cex.axis=1.25)
  axis(side=2,cex.axis=1.25,at=c(0.1,0.3,0.5,0.7,0.9),labels=c("0.1","0.3","0.5","0.7","0.9"))
  
  ylab = "HAB probability"
  if(j>=2) ylab= "" 
  title(xlab="TP (tons/day)",
        ylab=ylab,
        cex.lab = 1.5,
        line=2.25)
  
}

####################################################################################################
dev.off()