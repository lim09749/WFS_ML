rm(list=ls())
######################################################################################################################
##Read in data
setwd("C:/Users/xswang/HAB Research")
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL
usgs$State = NULL

usgs[usgs<0]=0
usgs = usgs*60*60*24*7
rivers = c("Tampa",
           "Myakka",
           "Peace",
           "Withlacoochee",
           "Manatee",
           "Hillsborough")
TN_files = paste0(rivers,"_TN.csv")
TP_files = paste0(rivers,"_TP.csv")
######################################################################################################################
nutrientcol = c(paste0(rivers,"_TN_mg"),paste0(rivers,"_TP_mg"))
totaldata = as.data.frame(matrix(NA, nrow = nrow(HAB), ncol=1))
par(mfrow=c(3,4))
for(i in 1:length(TN_files)){
  r = read.csv(TN_files[i])
  r$X = NULL
  r$Result_Value = as.numeric(as.character(r$Result_Value))
  r = r[!is.na(r$Result_Value),]
  
  dates = as.Date(unlist(strsplit(as.character(r$SampleDate),split=" "))[seq(from=1,to=1+2*(nrow(r)-1),by=2)],format="%m/%d/%Y")
  inds_to_get = sapply(HAB$Date, FUN = function(x) which.min(abs(x-dates)))
  nutrient = (usgs[,colnames(usgs)==rivers[i]]*28.3168)*(0.001*r$Result_Value[inds_to_get])
  plot(usgs[,colnames(usgs)==rivers[i]]*28.3168,nutrient)
  # nutrient[nutrient<quantile(nutrient,na.rm=T)[2]]=0
  totaldata = cbind(totaldata,nutrient)
  print(any(nutrient==0))
}
for(i in 1:length(TP_files)){
  r = read.csv(TP_files[i])
  r$X = NULL
  r$Result_Value = as.numeric(as.character(r$Result_Value))
  r = r[!is.na(r$Result_Value),]
  dates = as.Date(unlist(strsplit(as.character(r$SampleDate),split=" "))[seq(from=1,to=1+2*(nrow(r)-1),by=2)],format="%m/%d/%Y")
  inds_to_get = sapply(HAB$Date, FUN = function(x) which.min(abs(x-dates)))
  nutrient = (usgs[,colnames(usgs)==rivers[i]]*28.3168)*(0.001*r$Result_Value[inds_to_get])
  plot(usgs[,colnames(usgs)==rivers[i]]*28.3168,nutrient)
  # nutrient[nutrient<quantile(nutrient,na.rm=T)[2]]=0
  totaldata = cbind(totaldata,nutrient)
  print(any(nutrient==0))
}
totaldata$V1 = NULL
colnames(totaldata)=nutrientcol
totaldata = as.data.frame(totaldata)
write.csv(totaldata,"nutrientdata.csv")
######################################################################################################################
