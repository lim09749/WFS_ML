#################################################################################################################
##"contour plots ISEF.R"
##plot contour plots
#################################################################################################################
require(e1071)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Get Data
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
months = as.numeric(format(HAB$Date,"%m"))
righttime=!is.na(match(as.numeric(format(HAB$Date,"%m")),c(1:2,9:12)))

states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
#####################################################################################

##Read in entire dataset
alldata = read.csv("alldata.csv")
alldata$X = NULL
prevdata = read.csv("alldata.csv")
prevdata$X = NULL
prevdata$State = as.factor(states)

scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)
numrows=100

#####################################################################################
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
#####################################################################################
##Discharge contour plots
disdata = alldata[,c(1:9,10:17,31,34)]
dismeanvals = colMeans(disdata[!is.na(match(months,c(9:12,1:2))),1:(ncol(disdata)-1)])
# dismeanvals = colMeans(disdata[,1:(ncol(disdata)-1)])
disscalefactors = scale(prevdata[,c(1:9,10:17,31)])
load("dis_contour_svm.RData")
# HABsvm = tune.svm(as.factor(State)~.,
#                   data = disdata,
#                   type = "C-classification",
#                   cost = 2^(-5:10),
#                   kernel="radial",
#                   probability=T)
# HABsvm = HABsvm$best.model
# save(HABsvm,file="dis_contour_svm.RData")
xcol = 1
ycol = 13
addto = seq(-3, 2, length.out=numrows)
cmat = getcontour(xcol, ycol, dismeanvals, numrows, disdata, addto, addto, HABsvm)
tiff(file = "wind and discharge.tiff", width =10, height = 10, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mai=c(1,1,1,1))
contour(x=undoscale(xcol, addto+dismeanvals[xcol], disscalefactors),
        y=log10(10^undoscale(ycol, addto+dismeanvals[ycol], disscalefactors)/(60*60*24*7)*0.0283168),
        z=cmat,
        xlab= "NS wind at venf1 (m/s)",
        ylab= expression("Suwanee Discharge ("~m^{3}/s~")"),
        axes=F,
        lwd=3,
        cex.main=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 1, at = axTicks(1), labels=signif(axTicks(1),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)
dev.off()

xcol = 11
ycol = 13
addto = seq(-2.5, 2.5, length.out=numrows)
dismeanvals = colMeans(disdata[,1:(ncol(disdata)-1)])
cmat = getcontour(xcol, ycol, dismeanvals, numrows, disdata, addto, addto, HABsvm)
tiff(file = "two discharge.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mai=c(1,1,1,1))
contour(x=log10(10^undoscale(xcol, addto+dismeanvals[ycol], disscalefactors)/(60*60*24*7)*0.0283168),
        y=log10(10^undoscale(ycol, addto+dismeanvals[ycol], disscalefactors)/(60*60*24*7)*0.0283168),
        z=cmat,
        xlab= expression("Peace Discharge ("~m^{3}/s~")"),
        ylab= expression("Suwanee Discharge ("~m^{3}/s~")"),
        axes=F,
        lwd=3,
        cex.main=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)
dev.off()

xcol = 18
ycol = 15
addto = seq(-0.25, 1, length.out=numrows)
dismeanvals = colMeans(disdata[,1:(ncol(disdata)-1)])
cmat = getcontour(xcol, ycol, dismeanvals, numrows, disdata, addto, addto, HABsvm)
tiff(file = "two discharge 2.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mai=c(1,1,1,1))
contour(x=log10(10^undoscale(xcol, addto+dismeanvals[ycol], disscalefactors)/(60*60*24*7)*0.0283168),
        y=log10(10^undoscale(ycol, addto+dismeanvals[ycol], disscalefactors)/(60*60*24*7)*0.0283168),
        z=cmat,
        xlab= expression("Hillsborough Discharge ("~m^{3}/s~")"),
        ylab= expression("Caloosahatchee Discharge ("~m^{3}/s~")"),
        axes=F,
        lwd=3,
        cex.main=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)
dev.off()
#####################################################################################################
##Nutrient contour plots
nutrient_data = alldata
nutrient_scalefactors = scale(prevdata[,1:33])#[,c(1:9,18:29,32:33)])
nutrient_meanvals = colMeans(nutrient_data[righttime,1:(ncol(nutrient_data)-1)])
load("nutrient_contour_svm.RData")
# HABsvm = tune.svm(as.factor(State)~.,
#                   data = nutrient_data,
#                   type = "C-classification",
#                   cost = 2^(-5:10),
#                   kernel="radial",
#                   probability=T)
# HABsvm = HABsvm$best.model
# save(HABsvm,file="nutrient_contour_svm.RData")
################################################################################################
xcol = 20
ycol = 26
addto1 = seq(-2, 1.5, length.out=numrows)
addto2 = seq(-2, 0.25, length.out=numrows)
cmat = getcontour(xcol, ycol, nutrient_meanvals, numrows, nutrient_data, addto1, addto2, HABsvm)
tiff(file = "peace.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mai=rep(1,4))
contour(x=log10(10^undoscale(xcol, addto1+nutrient_meanvals[xcol], nutrient_scalefactors)/10^9/7),
        y=log10(10^undoscale(ycol, addto2+nutrient_meanvals[ycol], nutrient_scalefactors)/10^9/7),
        z=cmat,
        xlab= "Peace TN (tons/day)",
        ylab= "Peace TP (tons/day)",
        axes=F,
        lwd=3,
        cex.main=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)
dev.off()
########################################################################################################3
xcol = 23
ycol = 29
addto1 = seq(-2, 2, length.out=numrows)
addto2 = seq(-2, 2, length.out=numrows)
cmat = getcontour(xcol, ycol, nutrient_meanvals, numrows, nutrient_data, addto1, addto2, HABsvm)
tiff(file = "hillsborough.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mai=rep(1,4))
contour(x=log10(10^undoscale(xcol, addto1+nutrient_meanvals[xcol], nutrient_scalefactors)/10^9/7),
        y=log10(10^undoscale(ycol, addto2+nutrient_meanvals[ycol], nutrient_scalefactors)/10^9/7),
        z=cmat[nrow(cmat):1,],
        xlab= "Hillsborough TN (tons/day)",
        ylab= "Hillsborough TP (tons/day)",
        axes=F,
        lwd=3,
        col="blue",
        cex.main=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)
dev.off()
########################################################################################################3

##TN-TN
tn_data = alldata[,c(1:9,18:23,32,34)]
tn_scalefactors = scale(prevdata[,c(1:9,18:23,32)])
tn_meanvals = colMeans(tn_data[,1:(ncol(tn_data)-1)])

load("tn_contour_svm.RData")
# HABsvm = tune.svm(as.factor(State)~.,
#                   data = tn_data,
#                   type = "C-classification",
#                   cost = 2^(-5:10),
#                   kernel="radial",
#                   probability=T)
# HABsvm = HABsvm$best.model
# save(HABsvm,file="tn_contour_svm.RData")
################################################################################################
xcol = 12
ycol = 15
addto1 = seq(-2, 1.75, length.out=numrows)
addto2 = seq(-2, 2, length.out=numrows)

cmat = getcontour(xcol, ycol, tn_meanvals, numrows, tn_data, addto1, addto2, HABsvm)
tiff(file = "tn and tn.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mai=c(1,1,1,1))
contour(x=log10(10^undoscale(xcol, addto1+tn_meanvals[xcol], tn_scalefactors)/10^9/7),
        y=log10(10^undoscale(ycol, addto2+tn_meanvals[ycol], tn_scalefactors)/10^9/7),
        z=cmat,
        col="red",
        xlab="Peace TN (tons/day)",
        ylab="Hillsborough TN (tons/day)",
        axes=F,
        lwd=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)
dev.off()

########################################################################################################3
##TP-TP
tp_data = alldata[,c(1:9,24:29,33,34)]
tp_scalefactors = scale(prevdata[,c(1:9,24:29,33)])
tp_meanvals = colMeans(tp_data[,1:(ncol(tp_data)-1)])

load("tp_contour_svm.RData")
# HABsvm = tune.svm(as.factor(State)~.,
#                   data = tp_data,
#                   type = "C-classification",
#                   cost = 2^(-5:10),
#                   kernel="radial",
#                   probability=T)
# HABsvm = HABsvm$best.model
# save(HABsvm,file="tp_contour_svm.RData")
################################################################################################
xcol = 15
ycol = 16
addto1 = seq(-1, 1.5, length.out=numrows)
addto2 = seq(-1, 0.5, length.out=numrows)

cmat = getcontour(xcol, ycol, tp_meanvals, numrows, tp_data, addto1, addto2, HABsvm)
tiff(file = "tp and tp.tiff", width =7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mai=c(1,1,1,1))
contour(x=log10(10^undoscale(xcol, addto1+tp_meanvals[xcol], tp_scalefactors)/10^9/7),
        y=log10(10^undoscale(ycol, addto2+tp_meanvals[ycol], tp_scalefactors)/10^9/7),
        z=cmat,
        xlab= "Hillsborough TP (tons/day)",
        ylab= "Caloosahatchee TP (tons/day)",
        axes=F,
        lwd=3,
        cex.main=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)
dev.off()
####################################################################################################
######################################################################################################################################
nutrient_data = alldata[prevdata$Tampa_TN_mg>6 & prevdata$Caloosahatchee_TN_mg>5,
                        c(1:9,24:29,33,34)]#
nutrient_scalefactors = scale(prevdata[prevdata$Tampa_TN_mg>6 & prevdata$Caloosahatchee_TN_mg>5,
                                       c(1:9,24:29,33)])
nutrient_meanvals = colMeans(nutrient_data[,1:(ncol(nutrient_data)-1)])
load("two_tp_contour_svm.RData")
# HABsvm = tune.svm(as.factor(State)~.,
#                   data = nutrient_data,
#                   type = "C-classification",
#                   cost = 2^(-5:10),
#                   kernel="radial",
#                   probability=T)
# HABsvm = HABsvm$best.model
# save(HABsvm,file="two_tp_contour_svm.RData")

xcol = which(colnames(nutrient_data)=="Hillsborough_TP_mg")
ycol = which(colnames(nutrient_data)=="Caloosahatchee_TP_mg")
addto1 = seq(-2, 2, length.out=numrows)
addto2 = seq(-0.15, 1.25, length.out=numrows)
cmat = getcontour(xcol, ycol, nutrient_meanvals, numrows, nutrient_data, addto1, addto2, HABsvm)
tiff("contour_tp.tiff", width = 7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mai=c(1,1,1,1))
contour(x=log10(10^undoscale(xcol, addto1+nutrient_meanvals[xcol], nutrient_scalefactors)/10^9/7),
        y=log10(10^undoscale(ycol, addto2+nutrient_meanvals[ycol], nutrient_scalefactors)/10^9/7),
        z=cmat,
        xlab= "Hillsborough TP (tons/day)",
        ylab= "Caloosahatchee TP (tons/day)",
        axes=F,
        lwd=3,
        col="blue",
        cex.main=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)
dev.off()