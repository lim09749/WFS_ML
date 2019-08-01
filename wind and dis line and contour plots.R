require(e1071)
require("grDevices") # for colors
rm(list=ls())
setwd("D:/Li_Glibert_Backup")
# setwd("C:/Users/xswang/HAB Research/")
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
numrows=1000
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}
getcontour = function(col1, col2,dismeanvals, numrows, alldata, addto1, addto2, HABsvm){
  cmat = matrix(NA, nrow=numrows, ncol=numrows)
  for(i in 1:nrow(cmat)){
    for(j in 1:ncol(cmat)){
      curr = as.data.frame(rbind(dismeanvals,dismeanvals))
      colnames(curr) = colnames(alldata)[1:(ncol(alldata)-1)]
      curr[,col1] = curr[,col1]+addto1[i]
      curr[,col2] = curr[,col2]+addto2[j]
      cmat[i,j] = attr(predict(HABsvm,curr,
                               probability=T),"probabilities")[1,2]
    }
  }
  return(cmat)
}
##################################################################################################
##Plot it
tiff("wind_line_plots.tiff", width = 24, height = 12, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(1,2),mai=c(1.5,2,1.5,2))
# m <- matrix(c(rep(1:4,each=5),
#               rep(5:9,each=4)),nrow = 2,ncol = 20,byrow = TRUE)
# layout(mat = m,heights=c(1.25,1)/2.25)
#layout.show(n=9)
##wind plots
load("wind_svm.RData")
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(1:17)]##wind data
prevdata = alldata
scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)
meanvals = colMeans(alldata[,-ncol(alldata)])

x_axes = c("Normalized Wind","Normalized Wind")

# x_axes = c(expression(paste("Normalized N-S wind (",sigma,")")),
#            expression(paste("Normalized E-W wind (",sigma,")")))
for(j in 1:2){
  i=c(1:2,4:5)[j]
  xseq = seq(from=-2,
             to=2,
             length.out=numrows)
  dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
  dat[,i] = xseq
  dat = as.data.frame(dat)
  colnames(dat) = colnames(scalefactors)
  pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
  plot(xseq,#undoscale(i,xseq,scalefactors),
       pred,lwd=10,
       type="l",
       main="",
       xlab="",
       axes=F,#main="Station venf1",
       #ylab="HAB probability",
       ylab="",
       ylim=c(0.15,0.5),
       cex.main=2,cex.lab=2.5)
  title(ylab = "HAB probability", cex.lab = 3.5,
        line = 9)
  title(xlab = x_axes[j], cex.lab = 3.5,
        line = 7)
  axis(1,cex.axis=3,mgp=c(3, 2, 0))
  axis(2,las=1,cex.axis=3)
  abline(v=0,lwd=4.5,lty=2)
}
dev.off()
##########################################################################################################################
tiff("contour_plots1.tiff", width = 12, height = 12, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
alldata = read.csv("alldata.csv")
alldata$X = NULL
prevdata = read.csv("alldata.csv")
prevdata$X = NULL
prevdata$State = as.factor(states)

scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)
numrows=100

tn_data = alldata[,c(1:9,18:23,32,34)]
tn_scalefactors = scale(prevdata[,c(1:9,18:23,32)])
tn_meanvals = colMeans(tn_data[,1:(ncol(tn_data)-1)])

load("tn_contour_svm.RData")
xcol = 12
ycol = 15
addto1 = seq(-2, 2, length.out=numrows)
addto2 = seq(-2, 2, length.out=numrows)

cmat = getcontour(xcol, ycol, tn_meanvals, numrows, tn_data, addto1, addto2, HABsvm)
filled.contour(x=log10(10^undoscale(xcol, addto1+tn_meanvals[xcol], tn_scalefactors)/10^9/7),
               y=log10(10^undoscale(ycol, addto2+tn_meanvals[ycol], tn_scalefactors)/10^9/7),
               z=cmat,
               nlevels=20,
               color.palette = function(n) palette(grey(n:0/n)),
               xlab="",#Peace TN (tons/day)
               ylab="",#Hillsborough TN (tons/day)
               plot.axes={axis(side = 1, at = seq(-1.0,1.5,by=0.5),
                               labels=signif(10^seq(-1.0,1.5,by=0.5),digits=2),cex.axis=1.75)
                 axis(side = 2, at = seq(-1.5,1,by=0.5),
                      labels=signif(10^seq(-1.5,1,by=0.5),digits=2),cex.axis=1.75)},
               #axes=F,
               cex.lab=2,
               levels=seq(from=0,to=1,by=0.05),
               key.title = "",#title(main="HAB\nProbability"),
               key.axes = axis(4, seq(0, 1, by = 0.1),cex.axis=1.75))
contour(x=log10(10^undoscale(xcol, addto1+tn_meanvals[xcol], tn_scalefactors)/10^9/7),
        y=log10(10^undoscale(ycol, addto2+tn_meanvals[ycol], tn_scalefactors)/10^9/7),
        z=cmat,
        levels=seq(0,1,by=0.1),add=T)
dev.off()
##########################################################################################################################
tiff("contour_plots2.tiff", width = 12, height = 12, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))

tp_data = alldata[,c(1:9,24:29,33,34)]
tp_scalefactors = scale(prevdata[,c(1:9,24:29,33)])
tp_meanvals = colMeans(tp_data[,1:(ncol(tp_data)-1)])

load("tp_contour_svm.RData")
xcol = 15
ycol = 16
addto1 = seq(-2, 2, length.out=numrows)
addto2 = seq(-0.5, 0.8, length.out=numrows)

cmat = getcontour(xcol, ycol, tp_meanvals, numrows, tp_data, addto1, addto2, HABsvm)
filled.contour(x=log10(10^undoscale(xcol, addto1+tp_meanvals[xcol], tp_scalefactors)/10^9/7),
               y=log10(10^undoscale(ycol, addto2+tp_meanvals[ycol], tp_scalefactors)/10^9/7),
               z=cmat,
               nlevels=20,
               color.palette = function(n) palette(grey(n:0/n)),
               xlab= "",#Hillsborough TP (tons/day)",
               ylab="", #"Caloosahatchee TP (tons/day)",
               plot.axes={axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=2),cex.axis=1.5)
                 axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)
               },
               axes=F,
               cex.lab=2,
               levels=seq(from=0,to=1,by=0.05),
               key.title = title(main="HAB\nProbability"),
               key.axes = axis(4, seq(0, 1, by = 0.1)))
dev.off()
###########################################################################################################################
tiff("river_line_plots.tiff", width = 22, height = 12, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(3,5),mai=rep(0.75,4))
##discharge plots
load("dis_svm.RData")
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(1:9,10:17,31)]##river discharge
prevdata = alldata
scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)
plot_col = c(13,15,17,11,18)
numrows=1000
meanvals = colMeans(prevdata)
axes_ticks = rbind(log10(c(15,25,50,75,100,150,200,250,300,350,400,450)),
                   log10(c(0.5,1,2.5,5,10,15,20,25,30)),
                   log10(c(0.25,0.5,1,2.5,5,10,15,20,25,30)),
                   log10(c(0.5,1,2.5,5,25,50,100)),
                   log10(c(2.5,5,10,25,50,100,150)))
for(j in 1:length(plot_col)){
  i = plot_col[j]
  xseq = seq(from=quantile(prevdata[,i])[2]*0.9,
             to=mean(c(quantile(prevdata[,i])[4],quantile(prevdata[,i])[5])),
             length.out=numrows)
  dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
  dat[,i] = xseq
  dat = scale(dat, center = attr(scalefactors,"scaled:center"),
              scale = attr(scalefactors,"scaled:scale"))
  dat = as.data.frame(dat)
  colnames(dat) = colnames(scalefactors)
  pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
  plot(log10(10^xseq/(60*60*24*7)*0.0283168),pred,main=colnames(alldata)[i],lwd=2,
       xlab="",#,expression("Discharge ("~m^{3}/s~")"), 
       ylab="",#,"HAB probability",
       ylim=c(0.1,0.5),
       cex.main=3,cex.lab=2.5,xlim=c(max(c(-1,min(log10(10^xseq/(60*60*24*7)*0.0283168)))),
                                     max(log10(10^xseq/(60*60*24*7)*0.0283168))),axes=F)
  title(ylab = "HAB probability", cex.lab = 2.5,
        line = 5)
  title(xlab = expression("Discharge ("~m^{3}/s~")"), cex.lab = 2.5,
        line = 5)
  side1 = axes_ticks[j,]
  # axis(at=axTicks(side=1),labels=signif(10^axTicks(side=1),digits=2),side=1,cex.axis=2.25)
  axis(at=side1,labels=signif(10^side1,digits=2),side=1,cex.axis=2.25,mgp=c(3, 1.75, 0))
  axis(side=2,las=1,cex.axis=2.25)
}
############################################################################################################################
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(1:9,18:23,32)]##tn
states=states[alldata$Tampa_TN_mg>6 & alldata$Caloosahatchee_TN_mg>5]
alldata = alldata[alldata$Tampa_TN_mg>6 & alldata$Caloosahatchee_TN_mg>5,]
prevdata = alldata
scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)
load("tn_svm.RData")
plot_col = c(11,7,8,12)+4 ##TN
numrows=1000
rightrows = !is.na(match(format(HAB$Date,"%m"),c(9:12,1:2)))
meanvals = colMeans(prevdata)
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
  pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
  plot(log10(10^xseq/10^9/7),pred,main="",#paste0(strsplit(colnames(alldata)[i],split="_")[[1]][1]," TN"),
       xlab="", ylab="",ylim=c(0.1,0.7),
       col="black",cex.main=3,cex.lab=2.5,lwd=2,
       axes=F)
  title(ylab = "HAB probability", cex.lab = 2.5,
        line = 5)
  title(xlab = "TN (tons/day)", cex.lab = 2.5,
        line = 4)
  
  axes_ticks = rbind(log10(c(0.25,0.5,1,2.5)),
                     log10(c(0.1,0.25,0.5,1,2.5,5)),
                     log10(c(1,2.5,5,10,25)),
                     log10(c(2,5,10,20)))
  side1 = axes_ticks[j,]
  # axis(at=axTicks(side=1),labels=signif(10^axTicks(side=1),digits=2),side=1,cex.axis=2.25)
  axis(at=side1,labels=signif(10^side1,digits=2),side=1,cex.axis=2.25)
  axis(side=2,las=1,cex.axis=2.25)
}

############################################################################################################################
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(1:9,24:29,33)]##tp
prevdata = alldata

scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)

# require(beepr)
load("tp_svm.RData")
plot_col = c(11,7,8,12)+4 ##TP
numrows=1000
meanvals = colMeans(prevdata)
for(j in 1:length(plot_col)){
  i = plot_col[j]
  print(quantile(prevdata[,i]))
  xseq = seq(from=max(c(quantile(prevdata[,i])[1],6.5+log10(10))),
             to=mean(c(quantile(prevdata[,i])[4],quantile(prevdata[,i])[5])),
             length.out=numrows)
  if(j==3){
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
  pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
  plot(log10(10^xseq/10^9/7),pred,main="",#paste0(strsplit(colnames(alldata)[i],split="_")[[1]][1]," TP"),
       xlab="", ylab="",ylim=c(0.1,0.7),lwd=2,
       col="black",cex.main=3,cex.lab=2.5,
       axes=F)
  title(ylab = "HAB probability", cex.lab = 2.5,
        line = 5)
  title(xlab = "TP (tons/day)", cex.lab = 2.5,
        line = 4)
  
  axes_ticks = rbind(log10(c(0.025,0.05,0.1,0.25,0.5,1,1.5)),
                     log10(c(0.25,0.5,0.75,1,2,3,4,5,6)),
                     log10(c(0.05,0.05,0.5,0.1,0.25,0.5)),
                     log10(c(0.15,0.25,0.5,1,1.5,2,2.5)))
  side1 = axes_ticks[j,]
  axis(at=axTicks(side=1),labels=signif(10^axTicks(side=1),digits=2),side=1,cex.axis=2.25)
  #axis(at=side1,labels=signif(10^side1,digits=2),side=1,cex.axis=2.25)
  axis(side=2,las=1,cex.axis=2.25)
}

dev.off()