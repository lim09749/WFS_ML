require(e1071)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Get Data
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
dates = as.numeric(format(HAB$Date,format="%m"))
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
#####################################################################################

##Read in entire dataset
##Read in totaldata
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(1:9,10:17,31)]##river discharge
prevdata = alldata
scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)

load("dis_svm.RData")
# require(beepr)
# HABsvm = tune.svm(as.factor(State)~.,
#                   data=alldata,
#                   type="C-classification",
#                   cost = 2^(-5:10),
#                   kernel="radial",
#                   probability=T)
# HABsvm = HABsvm$best.model
# save(HABsvm,file="dis_svm.RData")
# beep(sound=3)
# 

tiff("dis_temp.tiff", width = 5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(3,2))
plot_col = c(11,c(13,15,17,18))
numrows=1000
rightrows = rep(TRUE, nrow(prevdata))
meanvals = colMeans(prevdata)
colors = c("red","blue","green","cyan","purple")
for(j in 1:length(colors)){
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
  # if(j==1){
    plot(log10(10^xseq/(60*60*24*7)*0.0283168),pred,main=colnames(alldata)[i],
         xlab=expression("Discharge ("~m^{3}/s~")"), ylab="HAB probability [%]",ylim=c(0.1,0.5),
         #col=colors[j],
         cex.main=3,cex.lab=1.5,xlim=c(max(c(-1,min(log10(10^xseq/(60*60*24*7)*0.0283168)))),
                                                     max(log10(10^xseq/(60*60*24*7)*0.0283168))),axes=F)
  # }else{
  #   points(log10(10^xseq/(60*60*24*7)*0.0283168),pred,
  #        col=colors[j],cex.main=3,cex.lab=1.5)
  # }
  axes_ticks = rbind(log10(c(0.5,1,2.5,5,25,50,100)),
                     log10(c(15,25,50,75,100,150,200,250,300,350,400,450)),
                     log10(c(0.5,1,2.5,5,10,15,20,25,30)),
                     log10(c(0.25,0.5,1,2.5,5,10,15,20,25,30)),
                     log10(c(2.5,5,10,25,50,100,150)))
  side1 = axes_ticks[j,]
  # axis(at=axTicks(side=1),labels=signif(10^axTicks(side=1),digits=2),side=1,cex.axis=1.75)
  axis(at=side1,labels=signif(10^side1,digits=2),side=1,cex.axis=1.5)
  axis(side=2,cex.axis=1.75)
}
dev.off()
