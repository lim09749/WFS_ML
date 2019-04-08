##############################################################################################
##"KBrevis Subset_west_florida_shelf.R"
##Creates boundary for west florida shelf
########################################################################################
##Read in florida in-situ Karenia Brevis cell counts
rm(list=ls())
setwd("C:/Users/xswang/HAB research")

HAB = read.csv(file="Florida_HAB.csv")
HAB$X <- NULL
##Ranges for needed region
lonrange = c(-83.1623, -81.09014)
latrange = c(25.8454, 29.1386)

plot(HAB$Longitude,HAB$Latitude)
abline(v=lonrange[1])
abline(v=lonrange[2])
abline(h=latrange[1])
abline(h=latrange[2])
##############################################################################################

##Keep in bounds
West_HAB = HAB[HAB$Latitude>=latrange[1] & HAB$Latitude<=latrange[2]
               & HAB$Longitude>=lonrange[1] & HAB$Longitude<=lonrange[2],]

##In bounds
points(West_HAB$Longitude,West_HAB$Latitude,col="red")
##############################################################################################
West_HAB_shore = West_HAB[West_HAB$dist_measurement_shore <= 9000,]

##Restrict it too shore
points(West_HAB_shore$Longitude,West_HAB_shore$Latitude,col="green")
##############################################################################################
##Save it
write.csv(West_HAB, file="West_Florida_HAB.csv")
write.csv(West_HAB_shore, file="West_Florida_Shore_HAB.csv")
