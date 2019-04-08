rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#######################################################
##Get lat lon
NDBC = read.csv("NDBC_lat_lon.csv")
NDBC$X = NULL
NDBC$Latitude = as.numeric(NDBC$Latitude)
NDBC$Longitude = -as.numeric(NDBC$Longitude)
#######################################################
##Ranges for needed region
lonrange = c(-83.1623, -81.09014)
latrange = c(25.8454, 29.1386)
#######################################################

plot(NDBC$Longitude,NDBC$Latitude)
abline(v=lonrange[1])
abline(v=lonrange[2])
abline(h=latrange[1])
abline(h=latrange[2])
NDBC$Station[NDBC$Longitude>=lonrange[1]&NDBC$Longitude<=lonrange[2]&
               NDBC$Latitude>=latrange[1]&NDBC$Latitude<=latrange[2]]

windstations = c("venf1","cdrf1","arpf1",
                 "fhpf1","clbf1","42013",
                 "42023","bgcf1")
#######################################################
