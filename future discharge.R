rm(list=ls())
require(ncdf4)
require(raster)
####################################################################################################
setwd("C:/Users/xswang/HAB Research")

usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL

usgsstations = c("2306647","2296750","2274325","2323500","2319000","2303330","2300500",
                 "2298830")
usgsloc = read.csv("USGS coord.csv")
usgsloc$X = NULL
usgsloc = usgsloc[match(as.numeric(usgsstations),usgsloc$site_no),]
####################################################################################################
##Increases in precipitation at different locations
pastfiles = c("C:/Users/xswang/HAB Research/RCM/pr/pr_RCM3_cgcm3_1996010103.nc")
futfiles = c("pr_RCM3_cgcm3_2041010103.nc",
          "pr_RCM3_cgcm3_2046010103.nc",
          "pr_RCM3_cgcm3_2051010103.nc",
          "pr_RCM3_cgcm3_2056010103.nc",
          "pr_RCM3_cgcm3_2061010103.nc",
          "pr_RCM3_cgcm3_2066010103.nc")
futfiles = paste0("C:/Users/xswang/HAB Research/RCM/pr/", futfiles)

pastvals = rep(NA, 260)
futvals = matrix(NA, nrow=260, 6)
latstations = usgsloc$lat[1]
lonstations = usgsloc$lon[1]

nc_data = nc_open(pastfiles[1])
lat     = (ncvar_get(nc_data, "lat"))
lon     = (360-ncvar_get(nc_data, "lon"))
pr     = (ncvar_get(nc_data,"pr"))

bestrow = which((lat-latstations)^2+(lon-lonstations)^2==min((lat-latstations)^2+(lon-lonstations)^2), arr.ind=TRUE)[1]
bestcol = which((lat-latstations)^2+(lon-lonstations)^2==min((lat-latstations)^2+(lon-lonstations)^2), arr.ind=TRUE)[2]
print(c(bestrow, bestcol))
dailyseq = seq(1, by=56, to=14600)
pr = pr[bestrow,bestcol,1:dim(pr)[3]]
pastvals = (sapply(dailyseq, FUN = function(x) mean(pr[x:(x+55)])))[1:260]
fromvals = c(34:52)
pastmean = mean(pastvals[c(fromvals,52+fromvals,156+fromvals,208+fromvals,104+fromvals)],
                na.rm=T)
# pastmean = mean(sort(pastvals,decreasing=T)[1:130])
# pastmean = mean(pastvals[c(fromvals,52+fromvals,156+fromvals,208+fromvals,104+fromvals)],
#                 na.rm=T)

futmeans = rep(NA, length(futfiles))
for(i in 1:length(futfiles)){
  nc_data = nc_open(futfiles[i])
  lat     = (ncvar_get(nc_data, "lat"))
  lon     = (360-ncvar_get(nc_data, "lon"))
  pr     = (ncvar_get(nc_data,"pr"))
  
  bestrow = which((lat-latstations)^2+(lon-lonstations)^2==min((lat-latstations)^2+(lon-lonstations)^2), arr.ind=TRUE)[1]
  bestcol = which((lat-latstations)^2+(lon-lonstations)^2==min((lat-latstations)^2+(lon-lonstations)^2), arr.ind=TRUE)[2]
  dailyseq = seq(1, by=56, to=14600)
  pr = pr[bestrow,bestcol,1:dim(pr)[3]]
  futvals[,i] = (sapply(dailyseq, FUN = function(x) mean(pr[x:(x+55)])))[1:260]
  futmeans[i] = mean(pr[c(fromvals,52+fromvals,104+fromvals,156+fromvals,208+fromvals)],na.rm=T)
  # futmeans[i] = mean(sort(futvals[,i],decreasing=T)[1:130])
}
print(colMeans(futvals,na.rm=T)/mean(pastvals))
print(futmeans/pastmean)
