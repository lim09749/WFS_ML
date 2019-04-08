#####################################################################################3
## Convert raw NDBC station data to daily data
#####################################################################################3
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/Wind buoys")
#####################################################################################3
##Get column names
files = list.files()

uniqcolumns = c()
for(f in files){
  station = read.csv(f)
  station$X = NULL
  uniqcolumns = unique(c(uniqcolumns, colnames(station)))
}
#####################################################################################3

