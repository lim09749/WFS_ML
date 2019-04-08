rm(list=ls())
require(ncdf4)
require(raster)
####################################################################################################
setwd("C:/Users/xswang/HAB Research")

##Get wind station locations
windstations = c("venf1","cdrf1","42039","42003")
windloc = read.csv("NDBC_lat_lon.csv")
windloc$X = NULL
windloc = windloc[match(windstations, windloc$Station),]

ndbc = read.csv("NDBC_station_data.csv")
ndbc = ndbc[,!is.na(colSums(ndbc))]
ndbc$X = NULL
####################################################################################################
##TAS
nc_data = nc_open("C:/Users/xswang/HAB Research/RCM/TAS/tas_RCM3_cgcm3_1996010103.nc")
lat     = (ncvar_get(nc_data, "lat"))
lon     = (360-ncvar_get(nc_data, "lon"))
tas     = (ncvar_get(nc_data,"tas"))

lat1 = windloc$Latitude[1]
lon1 = windloc$Longitude[1]
bestrow = which((lat-lat1)^2+(lon-lon1)^2==min((lat-lat1)^2+(lon-lon1)^2), arr.ind=TRUE)[1]
bestcol = which((lat-lat1)^2+(lon-lon1)^2==min((lat-lat1)^2+(lon-lon1)^2), arr.ind=TRUE)[2]
atmpraw = tas[bestrow,bestcol,1:dim(tas)[3]]

days = as.Date("1996-01-04")+seq(0,by=7,to=1825)
dailyseq = seq(4, by=56, to=14600)
atmp1 = (sapply(dailyseq, FUN = function(x) mean(atmpraw[x:(x+55)]))-270)[105:260]

ratmp = ndbc[1:length(atmp1),3]
plot(atmp1,ratmp)
print(summary(lm(ratmp~atmp1)))
intercept = lm(ratmp~atmp1)$coefficients[1]
slope = lm(ratmp~atmp1)$coefficients[2]
####################################################################################################
#Get free data
files = c("tas_RCM3_cgcm3_2041010103.nc",
          "tas_RCM3_cgcm3_2046010103.nc",
          "tas_RCM3_cgcm3_2051010103.nc",
          "tas_RCM3_cgcm3_2056010103.nc",
          "tas_RCM3_cgcm3_2061010103.nc",
          "tas_RCM3_cgcm3_2066010103.nc")
files = paste0("C:/Users/xswang/HAB Research/RCM/TAS/", files)

values = rep(NA, 260*6)
for(i in 1:length(files)){
  nc_data = nc_open(files[i])
  lat     = (ncvar_get(nc_data, "lat"))
  lon     = (360-ncvar_get(nc_data, "lon"))
  tas     = (ncvar_get(nc_data,"tas"))
  
  lat1 = windloc$Latitude[1]
  lon1 = windloc$Longitude[1]
  bestrow = which((lat-lat1)^2+(lon-lon1)^2==min((lat-lat1)^2+(lon-lon1)^2), arr.ind=TRUE)[1]
  bestcol = which((lat-lat1)^2+(lon-lon1)^2==min((lat-lat1)^2+(lon-lon1)^2), arr.ind=TRUE)[2]
  atmpraw = tas[bestrow,bestcol,1:dim(tas)[3]]
  
  dailyseq = seq(1, by=56, to=14600)
  atmp1 = (sapply(dailyseq, FUN = function(x) mean(atmpraw[x:(x+55)]))-270)[1:(length(dailyseq)-1)]
  values[(i*260-260+1):(i*260)]=atmp1*slope+intercept
}
write.csv(cbind(as.character(seq(as.Date("2041-01-01"), by = 7, length.out=length(values))),values),
          "Future ATMP 2038 to 2071.csv")
