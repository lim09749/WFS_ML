require(e1071)
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

##Read in entire dataset
##Read in totaldata
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(1:17)]##wind data
prevdata = alldata
scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)

load("wind_svm.RData")
# require(beepr)
# HABsvm = tune.svm(as.factor(State)~.,
#                   data=alldata,
#                   type="C-classification",
#                   cost = 2^(-5:10),
#                   kernel="radial",
#                   probability=T)
# HABsvm = HABsvm$best.model
# save(HABsvm,file="wind_svm.RData")
# beep(sound=3)
#############################################################################################################
par(mfrow=c(2,2))
plot_col = c(1:2,4:5)
numrows=1000#!is.na(match(months,c(9:12,1:2)))
rightrows = rep(TRUE, nrow(prevdata))
meanvals = colMeans(alldata[,-ncol(alldata)])
colors = c("red","blue","green","cyan")
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}

titles = c("NS wind (venf1)","EW wind (venf1)","NS wind (cdrf1)","EW wind (cdrf1)")
tiff(file = "wind_line_plots.tiff", height = 9, width = 9,  units = "in",
     pointsize=10, res = 600, compression = c("lzw"))
par(mfrow=c(2,2),mai=c(0.75,.75,0.75,0.75))
for(j in 1:length(colors)){
  i = plot_col[j]
  xseq = seq(from=-1.5,
             to=0.9,
             length.out=numrows)
  dat = matrix(rep(c(meanvals),numrows),nrow=numrows,byrow=TRUE)
  dat[,i] = xseq
  # dat = scale(dat, center = attr(scalefactors,"scaled:center"),
  #             scale = attr(scalefactors,"scaled:scale"))
  dat = as.data.frame(dat)
  colnames(dat) = colnames(scalefactors)
  pred = attr(predict(HABsvm,dat,probability=T),"probabilities")[,2]
  plot(undoscale(plot_col[j],xseq,scalefactors),pred,
       xlab="Wind speed (m/s)",axes=F,main=titles[j],
       ylab="HAB probability [%]",
       col=colors[j],cex.main=3,cex.lab=2)
  axis(1,cex.axis=2.25)
  axis(2,cex.axis=2.25)
  abline(v=0,col="blue")
}
dev.off()
