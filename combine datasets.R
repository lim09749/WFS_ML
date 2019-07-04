rm(list=ls())
#################################################################################################
##ndbc data
ndbc = read.csv("NDBC_station_data.csv")
ndbc$X = NULL

##usgs data
usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL
usgs$State = NULL
usgs[usgs<0]=0
usgs = usgs*60*60*24*7
usgs = log10(usgs+4)

##nutrient data
nutrient = read.csv("nutrientdata.csv")
nutrient$X = NULL
nutrient = log10(nutrient+4)

##altimetry data
ssh = read.csv("SSH.csv")
ssh$X = NULL
ssh$Week = NULL

##Caloosahatchee data
caloosahatchee = read.csv("caloosahatchee_data.csv")
caloosahatchee$X = NULL
caloosahatchee = log10(caloosahatchee+4)
#################################################################################################
##Combine data, save it as alldata
alldata = cbind(ndbc, usgs, nutrient, ssh, caloosahatchee)
colnames(alldata)[31:33] = c("Caloosahatchee","Caloosahatchee_TN_mg","Caloosahatchee_TP_mg")
for(i in 1:ncol(alldata)){
  alldata[is.na(alldata[,i]),i] = colMeans(alldata, na.rm=T)[i]
}
write.csv(alldata,"alldata.csv")
