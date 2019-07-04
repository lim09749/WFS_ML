rm(list=ls())
################################################################################################################
##Read in data
##Station 02292000
setwd("C:/Users/xswang/HAB Research")

HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)

criver = read.csv("caloosahatchee_discharge.csv",header=T)
criver = criver[2:nrow(criver),]
criver = criver[,c(3:4)]
criver$X171001_00060_00003 = as.numeric(as.character(criver$X171001_00060_00003))
dates = as.Date(as.character(criver$datetime),format="%m/%d/%Y")
################################################################################################################
##Make it weekly
weeks = seq(as.Date("1998-01-01"),to=as.Date("2018-09-27")+7,by=7)
weeklydischarge = rep(NA, length(weeks)-1)
for(i in 1:length(weeklydischarge)){
  weeklydischarge[i] = mean(criver$X171001_00060_00003[dates>=weeks[i] & dates < (weeks[i]+7)], 
                            na.rm=T)
}
weeklydischarge[is.na(weeklydischarge)] = mean(weeklydischarge,na.rm=T)
weeklydischarge = weeklydischarge*60*60*24*7
plot(weeks[1:(length(weeks)-1)],weeklydischarge)
################################################################################################################
##Find nutrient data
files = c("Caloosahatchee_TN.csv","Caloosahatchee_TP.csv")
dat=weeklydischarge
for(i in 1:length(files)){
  r = read.csv(files[i])
  r$X = NULL
  r = r[,c(10,16)]
  r = na.omit(r)
  r = r[!is.na(r$Result_Value),]
  
  dates = as.Date(unlist(strsplit(as.character(r$SampleDate),split=" "))[seq(from=1,to=1+2*(nrow(r)-1),by=2)],format="%m/%d/%Y")
  inds_to_get = sapply(HAB$Date, FUN = function(x) which.min(abs(x-dates)))
  
  nutrient = (weeklydischarge*28.3168)*(0.001*r$Result_Value[inds_to_get])
  print(any(r==0))
  # nutrient[nutrient<quantile(nutrient,na.rm=T)[2]]=0
  dat=cbind(dat,nutrient)
}
plot(dat[,1],dat[,2])
plot(dat[,1],dat[,3])
################################################################################################################
##Write data
colnames(dat) = c("discharge","TN","TP")
dat = as.data.frame(dat)
write.csv(dat,"caloosahatchee_data.csv")
