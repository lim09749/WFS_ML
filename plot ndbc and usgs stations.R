rm(list=ls())
require(PBSmapping)
require(sp)
require(raster)
######################################################################################################################
#### EDIT THIS
setwd("C:/Users/xswang/HAB Research")
# Import and prepare necessary points 
windstations = c("venf1","cdrf1","42039","42003")
usgsstations = c("2306647","2296750","2274325","2323500","2319000","2303330","2300500",
                 "2298830","2292000")
windloc = read.csv("NDBC_lat_lon.csv")
windloc$X = NULL
windloc = windloc[match(windstations, windloc$Station),]

usgsloc = read.csv("USGS coord.csv")
usgsloc$X = NULL
usgsloc = usgsloc[match(as.numeric(usgsstations),usgsloc$site_no),]

# Main Plot Coordinates 
Nlat <- 31.1
Slat <- 23.9
Elon <- 280.2 
Wlon <- 271.9
######################################################################################################################
shift = function(windsites){
  windsites$X = windsites$X-0.35
  windsites$Y = windsites$Y+0.2
  return(windsites)
}

# Import base and inset maps
zoomMap <- importGSHHS(xlim=c(Wlon, Elon), ylim=c(Slat, Nlat), maxLevel=3)
######################################################################################################################
windsites = windloc
windsites$X <- (-windloc$Longitude) %% 360 #Make Longitudegitude positive, and load in value
windsites$Y <- windloc$Latitude #Load Latitudeitude
windsites$EID <- 1:nrow(windsites) #Random number as event ID
attr(windsites, "projection")="LL" #Add number for projection type
windsites$label = windsites$Station
windsites = as.EventData(windsites) #Create event ID

usgssites = usgsloc
usgssites$X <- (usgsloc$lon) %% 360 #Make longitude positive, and load in value
usgssites$Y <- usgsloc$lat #Load latitude
usgssites$EID <- 1:nrow(usgssites) #Random number as event ID
attr(usgssites, "projection")="LL" #Add number for projection type
usgssites$label = usgssites$site_no
usgssites = as.EventData(usgssites) #Create event ID
######################################################################################################################

tiff(file = "NDBC USGS stations.tiff", width =3.75, height = 3.75, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
        ylab = "Latitude",cex=1) 
addPoints(windsites, col=2, pch=10, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=5,cex=0.25)
addLabels(shift(windsites),cex=0.85)
addPoints(usgssites, col=4, pch=10, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=5,cex=0.25)

legend(x=272.5, y= 25, 
       legend=c("NDBC stations", "USGS stations"), col=c(2,4), 
       lty=0, lwd=3, pch=c(10),cex=0.75)
scalebar(xy=c(272.5, 25.5), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1)
dev.off()

