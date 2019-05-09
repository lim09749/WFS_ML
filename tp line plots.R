require(e1071)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Get Data
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
#####################################################################################

##Read in entire dataset
##Read in totaldata
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(1:9,24:29,33)]##tp
# states= states[alldata$Caloosahatchee_TP_mg>=4]
# HAB = HAB[alldata$Caloosahatchee_TP_mg>=4,]
# 
# alldata = alldata[alldata$Caloosahatchee_TP_mg>=4,]
prevdata = alldata

scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)

# require(beepr)
load("tp_svm.RData")
# HABsvm = tune.svm(as.factor(State)~.,
#                   data=alldata,
#                   type="C-classification",
#                   cost = 2^(-5:10),
#                   kernel="radial",
#                   probability=T)
# HABsvm = HABsvm$best.model
# save(HABsvm,file="tp_svm.RData")
# beep(sound=3)

tiff("tp_temp.tiff", width = 6, height = 6, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
plot_col = c(7,8,11,12)+4 ##TP
numrows=1000
meanvals = colMeans(prevdata)
colors = c("cyan","red","green","purple")
par(mfrow=c(2,2))
for(j in 1:length(plot_col)){
  i = plot_col[j]
  print(quantile(prevdata[,i]))
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
  pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
  plot(log10(10^xseq/10^9/7),pred,main=paste0(strsplit(colnames(alldata)[i],split="_")[[1]][1]," TP"),
       xlab="TP (tons/day)", ylab="HAB probability [%]",ylim=c(0.1,0.7),
       col=colors[j],cex.main=2,cex.lab=1.5,
       axes=F)
  axes_ticks = rbind(log10(c(0.025,0.05,0.1,0.25,0.5,1,1.5)),
                     log10(c(0.25,0.5,0.75,1,2,3,4,5,6)),
                     log10(c(0.05,0.05,0.5,0.1,0.25,0.5)),
                     log10(c(0.15,0.25,0.5,1,1.5,2,2.5)))
  side1 = axes_ticks[j,]
  axis(at=axTicks(side=1),labels=signif(10^axTicks(side=1),digits=2),side=1,cex.axis=1.75)
  #axis(at=side1,labels=signif(10^side1,digits=2),side=1,cex.axis=1.75)
  axis(side=2,cex.axis=1.75)
}
dev.off()

