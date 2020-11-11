#################################################################################################################
##plot NEW contour plots
#################################################################################################################
require(kernlab)
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

##Remove NA values
x = is.na(HAB$Abundance_cells)

alldata = alldata[!x,]
HAB = HAB[!x,]
prevdata = prevdata[!x, ]

#####################################################################################
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}
fixbounds = function(x){
  y=x
  y[x<=0]=0
  y[x>=1]=1
  y
}

rvm_prob = function(HABrvm, data){
  fixbounds(predict(HABrvm, data))[,1]
}
getcontour = function(col1, col2,dismeanvals, numrows, alldata, addto1, addto2, HABrvm){
  cmat = matrix(NA, nrow=numrows, ncol=numrows)
  for(i in 1:nrow(cmat)){
    for(j in 1:ncol(cmat)){
      curr = as.data.frame(rbind(dismeanvals,dismeanvals))
      colnames(curr) = colnames(alldata)[1:(ncol(alldata)-1)]
      curr[,col1] = curr[,col1]+addto1[i]
      curr[,col2] = curr[,col2]+addto2[j]
      cmat[i,j] = rvm_prob(HABrvm, curr)[1]
    }
  }
  return(cmat)
}
#####################################################################################

##TN-TN
tn_data = alldata[,c(1:9,18:23,32,34)]
tn_scalefactors = scale(prevdata[,c(1:9,18:23,32)])
tn_meanvals = colMeans(tn_data[,1:(ncol(tn_data)-1)])

load(file="tn_contour_rvm.RData")

# HABrvm = kernlab::rvm(as.numeric(as.character(State))~., data= tn_data,
#                       kernel = "rbfdot")
# save(HABrvm,file="tn_contour_rvm.RData")
################################################################################################

tiff(file = "rvm new total March 20.tiff", width =15, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mai=c(1,1,1,1))
par(mfrow=c(1,2))
xcol = 12
ycol = 15
addto1 = seq(-1, 2.5, length.out=numrows)
addto2 = seq(-1.5, 0.5, length.out=numrows)

cmat = getcontour(xcol, ycol, tn_meanvals, numrows, tn_data, addto1, addto2, HABrvm)

contour(x=log10(10^undoscale(xcol, addto1+tn_meanvals[xcol], tn_scalefactors)/10^9/7),
        y=log10(10^undoscale(ycol, addto2+tn_meanvals[ycol], tn_scalefactors)/10^9/7),
        z=cmat,
        col="red",
        xlab="Peace TN (tons/day)",
        ylab="Hillsborough TN (tons/day)",
        axes=F,
        frame.plot=T,
        lwd=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)

########################################################################################################3
##TP-TP
tp_data = alldata[,c(1:9,24:29,33,34)]
tp_scalefactors = scale(prevdata[,c(1:9,24:29,33)])
tp_meanvals = colMeans(tp_data[,1:(ncol(tp_data)-1)])

load("tp_contour_rvm.RData")
# HABrvm = kernlab::rvm(as.numeric(as.character(State))~., data= tp_data,
#                       kernel = "rbfdot")
# save(HABrvm,file="tp_contour_rvm.RData")
################################################################################################
xcol = 15
ycol = 16
addto1 = seq(-1, 1.5, length.out=numrows)
addto2 = seq(-1, 0.5, length.out=numrows)

cmat = getcontour(xcol, ycol, tp_meanvals, numrows, tp_data, addto1, addto2, HABrvm)
contour(x=log10(10^undoscale(xcol, addto1+tp_meanvals[xcol], tp_scalefactors)/10^9/7),
        y=log10(10^undoscale(ycol, addto2+tp_meanvals[ycol], tp_scalefactors)/10^9/7),
        z=cmat,
        xlab= "Hillsborough TP (tons/day)",
        ylab= "Caloosahatchee TP (tons/day)",
        axes=F,
        frame.plot=T,
        lwd=3,
        col="blue",
        cex.main=3,cex.lab=2,
        labcex=2,
        levels=seq(from=0,to=1,by=0.05))
axis(side = 1, at = axTicks(1), labels=signif(10^axTicks(1),digits=2),cex.axis=1.5)
axis(side = 2, at = axTicks(2), labels=signif(10^axTicks(2),digits=2),cex.axis=1.5)
dev.off()
