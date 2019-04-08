##################################################################################################################################
## PROGRAM: plot HAB stations.R
## SPECIFICATIONS: Plots the location of sampling stations
#####################################################################################################################################
#### Working Directory
rm(list=ls())
setwd("C:/Users/xswang/HAB research")

#####################################################################################################################################
#### Libraries
require(PBSmapping)
require(sp)
require(raster)
######################################################################################################################
#### EDIT THIS

# Import and prepare necessary points 
sites = read.csv(file="Florida_HAB.csv")[,c(4,5)]#Load data
KBrevisLevel = read.csv(file="Florida_HAB.csv")[,11]
# Main Plot Coordinates 
Nlat <- 31.1
Slat <- 23.9
Elon <- 280.2 
Wlon <- 271.9
######################################################################################################################

# Import base and inset maps
zoomMap <- importGSHHS(xlim=c(Wlon, Elon), ylim=c(Slat, Nlat), maxLevel=3)
sites$X <- sites$Longitude %% 360 #Make Longitudegitude positive, and load in value
sites$Y <- sites$Latitude #Load Latitudeitude
sites$EID <- 1:nrow(sites) #Random number as event ID
attr(sites, "projection")="LL" #Add number for projection type
sites = as.EventData(sites) #Create event ID

# Parameterize and plot main map, add points listed above
tiff(file = "Florida HAB stations.tiff", width =3.75, height = 3.75, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
        ylab = "Latitude",cex=1) 
addPoints(sites, col=2, pch=10, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=5,cex=0.25)
legend(x=272, y= 24.75, 
       legend=c("Sampling Locations"), col=c(2), 
       lty=0, lwd=3, pch=c(19),cex=0.75)
scalebar(xy=c(272.5, 25.25), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1)
dev.off()
##Save shoreline
write.csv(zoomMap, "shore.csv")

######################################################################################################################
