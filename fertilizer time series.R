rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Get Data
n_fert = read.csv("Global Nitrogen Fertilizer.csv")
p_fert = read.csv("Global Phosphorus Fertilizer.csv")
n_fert$Date = as.numeric(substr(as.character(n_fert$Data),start=8,stop=nchar(as.character(n_fert$Data))))
p_fert$Date = as.numeric(substr(as.character(p_fert$Data),start=8,stop=nchar(as.character(p_fert$Data))))
n_fert$Data = NULL
p_fert$Data= NULL
colnames(n_fert)
######################################################################################
##USA, Western Europe, China, Southeastern Asia, Southern Africa
##Nitrogen
par(mfrow=c(2,1))
plot(n_fert$Date, n_fert$USA,type="l",ylim=c(0,4e+07),col="green",xlab="Year",ylab="N Metric Tons")
lines(n_fert$Date, n_fert$Western.Europe,ylim=c(0,4e+07),col="red")
lines(n_fert$Date, n_fert$China,ylim=c(0,4e+07),col="cyan")
lines(n_fert$Date, n_fert$Southeastern.Asia,ylim=c(0,4e+07),col="purple")
lines(n_fert$Date, n_fert$Southern.Africa,ylim=c(0,4e+07),col="orange")
##Phosphorus
plot(p_fert$Date, p_fert$USA,type="l",ylim=c(0,2e+07),col="green",xlab="Year",ylab="P Metric Tons")
lines(p_fert$Date, p_fert$Western.Europe,ylim=c(0,2e+07),col="red")
lines(p_fert$Date, p_fert$China,ylim=c(0,2e+07),col="cyan")
lines(p_fert$Date, p_fert$Southeastern.Asia,ylim=c(0,2e+07),col="purple")
lines(p_fert$Date, p_fert$Southern.Africa,ylim=c(0,2e+07),col="orange")
