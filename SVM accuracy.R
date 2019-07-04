require(e1071)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Get HAB
HAB = read.csv("C:/Users/xswang/HAB Research/West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Abundance_cells[is.na(HAB$Abundance_cells)] = 0
HAB$Date = as.Date(HAB$Date)

states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1

ndbc = read.csv("NDBC_station_data.csv")
ndbc = ndbc[,!is.na(colSums(ndbc))]
ndbc$X = NULL
usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL
#####################################################################################
##Log scale to adjust discharge (add to avoid negatives)
usgs[,1:(ncol(usgs)-1)] = log10(usgs[,1:(ncol(usgs)-1)]+4)
#####################################################################################
alldata = cbind(ndbc,usgs)
for(i in 1:ncol(alldata)) alldata[,i] = as.numeric(as.character(alldata[,i]))
alldata$State = as.factor(alldata$State)
alldata[,1:(ncol(alldata)-1)] = scale(alldata[,1:(ncol(alldata)-1)])
prevdata = alldata
alldata = alldata[sample(nrow(alldata)),]

#####################################################################################
##Get HAB values (use fact that all of peace river data is unique)
for(i in 1:ncol(alldata)) print(c(colnames(alldata)[i],length(unique(alldata[,i]))))
HABval = rep(NA, nrow(alldata))
for(i in 1:length(HABval)){
  HABval[i] = HAB$Abundance_cells[which.min(abs(alldata$Peace[i]-prevdata$Peace))] 
}
HABstate = rep(0, length(HABval))
HABstate[HABval >= 10^5] = 1
####################################################################################
##Read in predictions
testpred = read.csv("testpred.csv")
testpred$X = NULL
colnames(testpred) = c("Real","SVM","NB","RVM","NET")
HABstate = testpred$Real
##Match up dates with predictions
date = rep(as.Date("1900-01-01"), nrow(testpred))
for(i in 1:nrow(testpred)){
  date[i] = HAB$Date[which.min(abs(alldata$Peace[i]-prevdata$Peace))] 
}
rightsvm = sort(date[testpred[,2]==HABstate])
rightbay = sort(date[testpred[,3]==HABstate])
rightrvm = sort(date[testpred[,4]==HABstate])
rightnet = sort(date[testpred[,5]==HABstate])
####################################################################################
##Histograms
hist(rightsvm,breaks=as.Date(paste0(1998:2019,"-01-01")))
####################################################################################
##K. Brevis cell counts and accuracy
relFreq <- function(counts){
  return(counts/sum(counts))
}

tiff(file = "hist_comparisons.tiff", height = 7.5, width = 7.5,  units = "in",
     pointsize=10, res = 600, compression = c("lzw"))
par(mfrow=c(2,1))
h1 = hist(log10(HAB$Abundance_cells[is.na(match(HAB$Date,rightsvm))]),breaks=1:8,plot=F)
h2 = hist(log10(HAB$Abundance_cells),breaks=1:8,plot=F)
h1$counts <- relFreq(h1$counts)
h2$counts <- relFreq(h2$counts)
plot(h1,ylim=c(0,0.3),xlab="K. Brevis cell counts (log10(c/l))", main="Histogram of SVMs correct predictions")
plot(h2,ylim=c(0,0.3),xlab="K. Brevis cell counts (log10(c/l))", main="Histogram")
dev.off()
