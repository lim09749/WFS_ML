rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Read in data
HAB = read.csv("C:/Users/xswang/HAB Research/West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)

ndbc = read.csv("NDBC_station_data.csv")
ndbc = ndbc[,!is.na(colSums(ndbc))]
ndbc$X = NULL
usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL

##Log scale to adjust discharge (add to avoid negatives)
usgs[,1:(ncol(usgs)-1)] = log10(usgs[,1:(ncol(usgs)-1)]+4)

##Create entire dataset
alldata = cbind(ndbc,usgs)
for(i in 1:ncol(alldata)) alldata[,i] = as.numeric(as.character(alldata[,i]))
alldata$State = as.factor(alldata$State)
alldata[,1:(ncol(alldata)-1)] = scale(alldata[,1:(ncol(alldata)-1)])
#####################################################################################
pca <- prcomp(alldata[,1:(ncol(alldata)-1)], center = TRUE,scale. = TRUE)
summary(pca)
#plot(pca)
library(devtools)
library(ggbiplot)

tiff(file = "PCA Analysis.tiff", width = 7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
ggbiplot(pca, obs.scale = 1, var.scale = 1, circle = T, ellipse=TRUE, groups=alldata$State)
dev.off()