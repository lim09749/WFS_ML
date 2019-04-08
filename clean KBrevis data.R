rm(list=ls())
require(rgdal)

setwd("C:/Users/xswang/HAB research")
##########################################################################################
##Load in data
HAB <- read.csv("Karenia & Karlodinium 1998-2018.11.05.csv")
colnames(HAB) <- c("Date", "Time", "Time_Zone",
                   "Depth", "Site", "Latitude",
                   "Longitude", "County", "Taxa",
                   "Abundance_cells", "Abundance_Notes", "Temp",
                   "Salinity", "DO_percent", "DO",
                   "Collection_Agency")

##Only for Karenia Brevis
HAB <- HAB[HAB$Taxa=="Karenia brevis", ]
##########################################################################################
##Format dates (multiple formats)
HAB$Date <- as.character(HAB$Date)
dates <- rep("0000-01-01", length(HAB$Date))
dates <- as.Date(dates)

##Format 1: 2012-08-18
index1 <- which(!is.na(as.Date(HAB$Date, format = "%m/%d/%Y")))
dates[index1] <- as.Date(HAB$Date, format = "%m/%d/%Y")[index1]

##Since year is stored since 2000, we need to add 730485 to go from 0000 to 2000
##Format 2: 19-Aug-12
index2 <- which(!is.na(as.Date(HAB$Date, format = "%d-%b-%Y")+730485))
dates[index2] <- as.Date(HAB$Date, format = "%d-%b-%Y")[index2]+730485

dates <- as.Date(dates)
HAB$Date <- dates

##Remove DO[%], Time, Time_Zone, Abundance_Notes, Site, Taxa, Collection_Agency
HAB$DO_percent        <- NULL
HAB$Time              <- NULL
HAB$Time_Zone         <- NULL
HAB$Abundance_Notes   <- NULL
HAB$Site              <- NULL
HAB$Taxa              <- NULL
HAB$Collection_Agency <- NULL
#HAB$County          <- NULL
#HAB$Site <- as.character(HAB$Site)
HAB$County <- as.character(HAB$County)

#HAB <- HAB[HAB$Abundance_cells != 0,]
HAB <- HAB[!is.na(HAB$Abundance_cells),]

type <- c("not present/background", "very low", "low",
          "medium", "high")
KBrevisLevel <- rep(NA, nrow(HAB))
KBrevisLevel[HAB$Abundance_cells < 1000]  = 1
KBrevisLevel[HAB$Abundance_cells >= 1000 & HAB$Abundance_cells < 10000]  = 2
KBrevisLevel[HAB$Abundance_cells >= 10000 & HAB$Abundance_cells < 100000]  = 3
KBrevisLevel[HAB$Abundance_cells >= 100000 & HAB$Abundance_cells < 1000000]  = 4
KBrevisLevel[HAB$Abundance_cells >= 1000000]  = 5
hist(KBrevisLevel)

HAB <- cbind(HAB, KBrevisLevel)

#####################################################################################################################################
#Convert Distance to meters by converting lat lon to UTM 

shoreline = read.csv("shore.csv")[, 5:6]
HABcoord = cbind(HAB$Longitude%%360, HAB$Latitude)

convertToMeters <- function(lon, lat){
  xy              <- data.frame(X = c(lon), Y = c(lat))#,ID = c(EventId))
  
  coordinates(xy) <- c("X", "Y")
  
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  
  coordStation             <- as.data.frame(spTransform(xy, CRS("+proj=utm +zone=17 ellps=WGS84")))
  #coordStation$ID <- NULL
  
  coordStation <-as.matrix(coordStation, ncol=2) #Coordinates of stations
  return(coordStation)
}

##Closest shoreline
shoreline_meters <- convertToMeters(shoreline[,1],shoreline[,2])
HAB_meters <- convertToMeters(HABcoord[,1],HABcoord[,2])

nearestInd <- apply(HAB_meters, MARGIN =c(1), FUN = function(x) 
  which.min(rowSums(sweep(shoreline_meters[seq(1,nrow(shoreline),by=10),],2,x)^2)))
closest_shoreline <- shoreline_meters[seq(1,nrow(shoreline),by=10)[nearestInd],]
dist_measurement_shore <- sqrt(rowSums((closest_shoreline-HAB_meters)^2))


##Add to HAB
HAB <- cbind(HAB, dist_measurement_shore)

#####################################################################################################################################
##Save data
write.csv(HAB, "Florida_HAB.csv")
##########################################################################################

# dmeans <- rep(NA, 5)
# for(i in 1:5) dmeans[i] = mean(HAB$[KBrevisLevel==i], na.rm=T)
# plot(dmeans, 1:5)
# 
# unique_county <- unique(HAB$County)
# loc_avg <- rep(0, length(unique_county))
# pos_avg <- rep(0, length(unique_county))
# for(i in 1:nrow(HAB)){
#   pos = which(HAB$County[i]==unique_county)
#   loc_avg[pos] <- loc_avg[pos] + HAB$Abundance_cells[i]
#   pos_avg[pos] <- pos_avg[pos] + 1
# }
# 
