rm(list=ls())
######################################################################################################################
setwd("C:/Users/xswang/HAB Research")

# Import and prepare necessary points 
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1

ndbc = read.csv("NDBC_station_data.csv")
ndbc = ndbc[,!is.na(colSums(ndbc))]
ndbc$X = NULL
usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL
totaldata = cbind(ndbc, usgs)
totaldata$State = NULL
######################################################################################################################
rivers = c("Tampa",
           "Myakka",
           "Peace",
           "Withlacoochee",
           "Manatee",
           "Hillsborough")
TN_files = paste0(rivers,"_TN.csv")
TP_files = paste0(rivers,"_TP.csv")
######################################################################################################################
nutrientcol = c(colnames(totaldata),paste0(rivers,"_TN_mg"),paste0(rivers,"_TP_mg"),"State")
par(mfrow=c(3,4))
for(i in 1:length(TN_files)){
  r = read.csv(TN_files[i])
  r$X = NULL
  dates = as.Date(unlist(strsplit(as.character(r$SampleDate),split=" "))[seq(from=1,to=1+2*(nrow(r)-1),by=2)],format="%m/%d/%Y")
  inds_to_get = sapply(HAB$Date, FUN = function(x) which.min(abs(x-dates)))
  nutrient = (totaldata[,colnames(totaldata)==rivers[i]]*28.3168)*(0.001*r$Result_Value[inds_to_get])
  plot(totaldata[,colnames(totaldata)==rivers[i]]*28.3168,nutrient)
  totaldata = cbind(totaldata,nutrient)
}
for(i in 1:length(TP_files)){
  r = read.csv(TP_files[i])
  r$X = NULL
  dates = as.Date(unlist(strsplit(as.character(r$SampleDate),split=" "))[seq(from=1,to=1+2*(nrow(r)-1),by=2)],format="%m/%d/%Y")
  inds_to_get = sapply(HAB$Date, FUN = function(x) which.min(abs(x-dates)))
  nutrient = (totaldata[,colnames(totaldata)==rivers[i]]*28.3168)*(0.001*r$Result_Value[inds_to_get])
  plot(totaldata[,colnames(totaldata)==rivers[i]]*28.3168,nutrient)
  totaldata = cbind(totaldata,nutrient)
}
totaldata = cbind(totaldata,states)
colnames(totaldata)=nutrientcol
totaldata = totaldata[,c(1:9,12:13,18:(ncol(totaldata)))]
write.csv(totaldata,"totaldata.csv")
######################################################################################################################
