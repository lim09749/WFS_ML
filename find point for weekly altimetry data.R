rm(list=ls())
require(PBSmapping)
require(sp)
require(raster)
require(rgdal)
setwd("C:/Users/xswang/HAB Research/")
##################################################################################################################

r = stack("global-reanalysis-phy-001-030-daily_1554944976012.nc")
days = as.Date(substr(names(r),2,11),"%Y.%m.%d")

rightrowsabove = coordinates(r)[,2]>=26.3 & coordinates(r)[,2]<=26.7
rightrowsbelow = coordinates(r)[,2]>=24.3 & coordinates(r)[,2]<=24.7

# Main Plot Coordinates 
Nlat <- 26.75
Slat <- 24.25
Elon <- 280.2 
Wlon <- 271.9
zoomMap <- importGSHHS(xlim=c(Wlon, Elon), ylim=c(Slat, Nlat), maxLevel=3)

######################################################################################################################

# Import base and inset maps
##thanks to https://topex.ucsd.edu/cgi-bin/get_data.cgi
##for the isobath data
isobath = read.table("isobath.txt")
colnames(isobath) = c("lat","lon","isobath")
isobath = isobath[abs(isobath$isobath+300)<10,]
iso  = as.data.frame(cbind(isobath$lat,isobath$lon))
iso$X = isobath$lat
iso$Y = isobath$lon
iso$EID <- 1:nrow(iso) #Random number as event ID
attr(iso, "projection")="LL" #Add number for projection type
iso = as.EventData(iso) #Create event ID


sitesabove = as.data.frame(coordinates(r)[rightrowsabove,])
sitesabove$X = 360+coordinates(r)[rightrowsabove,1]
sitesabove$Y = coordinates(r)[rightrowsabove,2]
sitesabove$EID <- 1:nrow(sitesabove) #Random number as event ID
attr(sitesabove, "projection")="LL" #Add number for projection type
sitesabove = as.EventData(sitesabove) #Create event ID

sitesbelow = as.data.frame(coordinates(r)[rightrowsbelow,])
sitesbelow$X = 360+coordinates(r)[rightrowsbelow,1]
sitesbelow$Y = coordinates(r)[rightrowsbelow,2]
sitesbelow$EID <- 1:nrow(sitesbelow) #Random number as event ID
attr(sitesbelow, "projection")="LL" #Add number for projection type
sitesbelow = as.EventData(sitesbelow) #Create event ID

# Parameterize and plot main map, add points listed above
plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
        ylab = "Latitude",cex=1) 

addPoints(iso, col=2, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=5,cex=0.25)

##How i found the right point
addPoints(sitesabove, col=1, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=5,cex=0.25)
addPoints(sitesbelow, col=1, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=5,cex=0.25)

identify(sitesabove$X, sitesabove$Y)
identify(sitesbelow$X, sitesbelow$Y)

addPoints(sitesabove[c(28,161,294,427,561),], col=3, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=5,cex=0.25)
addPoints(sitesbelow[c(40,174,309,444,580),], col=3, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=5,cex=0.25)
#########################################################################################################
print(which(rightrowsabove)[c(28,161,294,427,561)]) #7210
print(which(rightrowsbelow)[c(40,174,309,444,580)]) #10417
