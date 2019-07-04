rm(list=ls())
require(raster)
require(rgdal)
setwd("C:/Users/xswang/HAB Research/")
##################################################################################################################

##Get stack and data
r = stack("global-reanalysis-phy-001-030-daily_1554944976012.nc")
days = as.Date(substr(names(r),2,11),"%Y.%m.%d")

ind_above = c(6944, 7077,7210,7343, 7477)
ind_below = c(10148, 10282,10417,10552, 10688)

##Get values for above ones
values_r = values(r)
abovevalues = apply(t(values_r[ind_above,]), 1, max)
belowvalues = apply(t(values_r[ind_below,]), 1, min)
plot(days,abovevalues-belowvalues)

##Get weekly averages
weeks = seq(as.Date("1998-01-01"),to=as.Date("2018-09-27")+7,by=7)
weekly_above_values = sapply(1:(length(weeks)-1), 
                             FUN = function(x) mean(abovevalues[days>=weeks[x]&days<weeks[x+1]]))
weekly_below_values = sapply(1:(length(weeks)-1), 
                             FUN = function(x) mean(belowvalues[days>=weeks[x]&days<weeks[x+1]]))
lines(weeks[1:(length(weeks)-1)],weekly_above_values-weekly_below_values)

##Write dat
dat = cbind(as.character(weeks[1:(length(weeks)-1)]), 
            as.character(weekly_above_values-weekly_below_values))
colnames(dat) = c("Week","SSH")

##Impute data
##averages of same month
for(i in min(which(is.na(weekly_above_values))):length(weekly_above_values)){
 mon=as.numeric(format(weeks[i],"%m"))
 dat[i,2]=mean(abovevalues[as.numeric(format(days,"%m"))==mon]
               -belowvalues[as.numeric(format(days,"%m"))==mon])
}

dat = as.data.frame(dat)
plot(as.Date(dat$Week), as.numeric(as.character(dat$SSH)))
##Write csv
write.csv(dat, "SSH.csv")

