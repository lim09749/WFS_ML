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
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = as.data.frame(scale(alldata))
#####################################################################################
##Peace discharge
tiff("hist_peace_dis.tiff", width = 4, height = 6, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(2,1))
h1=hist(alldata$Peace[states==0],freq=F,main="Peace during weeks without HABs",xlab="Standardized Values")
hist(alldata$Peace[states==1],freq=F,main="Peace during weeks with HABs",breaks=h1$breaks,col="red",xlab="Standardized Values")
dev.off()
#####################################################################################
##Caloosahatchee TN
tiff("hist_caloosa_tn.tiff", width = 4, height = 6, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(2,1))
h1 = hist(alldata$Caloosahatchee_TN_mg[states==0],freq=F,main="Caloosahatchee TN during weeks without HABs",xlab="Standardized Values")
hist(alldata$Caloosahatchee_TN_mg[states==1],freq=F,main="Caloosahatchee TN during weeks with HABs",breaks=h1$breaks,col="red",xlab="Standardized Values")
dev.off()
#####################################################################################
##Caloosahatchee TP
tiff("hist_caloosa_tp.tiff", width = 4, height = 6, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(2,1))
h1 = hist(alldata$Caloosahatchee_TP_mg[states==0],freq=F,main="Caloosahatchee TP during weeks without HABs",xlab="Standardized Values")
hist(alldata$Caloosahatchee_TP_mg[states==1],freq=F,main="Caloosahatchee TP during weeks with HABs",breaks=h1$breaks,col="red",xlab="Standardized Values")
dev.off()