rm(list=ls())
require(PBSmapping)
require(sp)
require(raster)

setwd("C:/Users/xswang/HAB Research/")
########################################################################################################
##HAB months
##Bloom years: 2002 (01-03), 2005 (01/05-01/06), 2012, 2017-2018
bloom_beg = as.Date(c("2002-01-01","2018-07-26","2005-11-01","2012-11-01"))

##no-HAB months
#No bloom years: 1998, 2007(?), 2010 (03/10-10/11), 2013
nobloom_beg = as.Date(c("1998-10-01","2010-10-01","2007-06-01","2018-01-01"))
########################################################################################################
HABtotal = read.csv("Karenia & Karlodinium 1998-2018.11.05.csv")
colnames(HABtotal) <- c("Date", "Time", "Time_Zone",
                        "Depth", "Site", "Latitude",
                        "Longitude", "County", "Taxa",
                        "Abundance_cells", "Abundance_Notes", "Temp",
                        "Salinity", "DO_percent", "DO",
                        "Collection_Agency")

##Only for Karenia Brevis
HABtotal <- HABtotal[HABtotal$Taxa=="Karenia brevis", ]
##########################################################################################
##Format dates (multiple formats)
HABtotal$Date <- as.character(HABtotal$Date)
dates <- rep("0000-01-01", length(HABtotal$Date))
dates <- as.Date(dates)

##Format 1: 2012-08-18
index1 <- which(!is.na(as.Date(HABtotal$Date, format = "%m/%d/%Y")))
dates[index1] <- as.Date(HABtotal$Date, format = "%m/%d/%Y")[index1]

##Since year is stored since 2000, we need to add 730485 to go from 0000 to 2000
##Format 2: 19-Aug-12
index2 <- which(!is.na(as.Date(HABtotal$Date, format = "%d-%b-%Y")+730485))
dates[index2] <- as.Date(HABtotal$Date, format = "%d-%b-%Y")[index2]+730485

dates <- as.Date(dates)
HABtotal$Date <- dates


##Remove DO[%], Time, Time_Zone, Abundance_Notes, Site, Taxa, Collection_Agency
HABtotal$DO_percent        <- NULL
HABtotal$Time              <- NULL
HABtotal$Time_Zone         <- NULL
HABtotal$Abundance_Notes   <- NULL
HABtotal$Site              <- NULL
HABtotal$Taxa              <- NULL
HABtotal$Collection_Agency <- NULL
#HABtotal$County          <- NULL
#HABtotal$Site <- as.character(HABtotal$Site)
HABtotal$County <- as.character(HABtotal$County)

HABtotal$Date = as.Date(HABtotal$Date)

###############################################################################################
bloomtoplot1 = cbind(HABtotal$Latitude,
                     HABtotal$Longitude, HABtotal$Abundance_cells)[HABtotal$Date>=bloom_beg[1] 
                                                                   & HABtotal$Date<bloom_beg[1]+30,]
bloomtoplot2 = cbind(HABtotal$Latitude,
                     HABtotal$Longitude, HABtotal$Abundance_cells)[HABtotal$Date>=bloom_beg[2] 
                                                                   & HABtotal$Date<bloom_beg[2]+30,]
bloomtoplot3 = cbind(HABtotal$Latitude,
                     HABtotal$Longitude, HABtotal$Abundance_cells)[HABtotal$Date>=bloom_beg[3] 
                                                                   & HABtotal$Date<bloom_beg[3]+30,]
bloomtoplot4 = cbind(HABtotal$Latitude,
                     HABtotal$Longitude, HABtotal$Abundance_cells)[HABtotal$Date>=bloom_beg[4] 
                                                                   & HABtotal$Date<bloom_beg[4]+30,]

nobloomtoplot1 = cbind(HABtotal$Latitude,
                     HABtotal$Longitude, HABtotal$Abundance_cells)[HABtotal$Date>=nobloom_beg[1] 
                                                                   & HABtotal$Date<nobloom_beg[1]+30,]
nobloomtoplot2 = cbind(HABtotal$Latitude,
                     HABtotal$Longitude, HABtotal$Abundance_cells)[HABtotal$Date>=nobloom_beg[2] 
                                                                   & HABtotal$Date<nobloom_beg[2]+30,]
nobloomtoplot3 = cbind(HABtotal$Latitude,
                     HABtotal$Longitude, HABtotal$Abundance_cells)[HABtotal$Date>=nobloom_beg[3] 
                                                                   & HABtotal$Date<nobloom_beg[3]+30,]
nobloomtoplot4 = cbind(HABtotal$Latitude,
                       HABtotal$Longitude, HABtotal$Abundance_cells)[HABtotal$Date>=nobloom_beg[4] 
                                                                     & HABtotal$Date<nobloom_beg[4]+30,]

###############################################################################################
##Help get data
Nlat <- 28.5
Slat <- 26
Elon <- 278.5
Wlon <- 276.5
zoomMap <- importGSHHS(xlim=c(Wlon, Elon), ylim=c(Slat, Nlat), maxLevel=3)

plotblooms = function(toplot){
  sites = as.data.frame(matrix(nrow = nrow(toplot), ncol=3))
  colnames(sites) = c("X","Y","EID")
  sites$X = toplot[,2] %% 360
  sites$Y = toplot[,1]
  sites$EID = 1:nrow(sites)
  attr(sites, "projection")="LL"
  sites = as.EventData(sites)
  return(sites)
}
getColors = function(toplot){
  colval = rep(NA, nrow(toplot))
  colval[toplot[,3]<1000] = "black"
  colval[toplot[,3]>=1000 & toplot[,3] < 10000] = "green"
  colval[toplot[,3]>=10000 & toplot[,3] <100000] = "yellow"
  colval[toplot[,3]>=100000 & toplot[,3]<1000000] = "orange"
  colval[toplot[,3]>=1000000] = "red"
  return(colval)
}
###############################################################################################
##During HAB
divide = 1
tiff(file = "During HAB maps.tiff", width =13, height = 15, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(2,2))
plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
        ylab = "Latitude",
        main = "January, 2002",
        cex.main=3,
        cex.axis=1.5,
        cex.lab=2) 
addPoints(plotblooms(bloomtoplot1), col=getColors(bloomtoplot1), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(bloomtoplot1[,3]+1)+1)/divide)
legend(x=276, y= 29, 
       legend=c("Background (0-10^3 cells/L)", 
                "Very Low (10^3-10^4 cells/L)", 
                "Low (10^4-10^5 cells/L)",
                "Medium (10^5-10^6 cells/L)", 
                "High (10^6-10^7 cells/L)"), 
       col=c("black","green","yellow","orange","red"), 
       lty=0, lwd=3, pch=c(19),cex=1.5)
scalebar(xy=c(277, 25.25), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1.5)

plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
        ylab = "Latitude",
        main = "July, 2017",
        cex.main=3,
        cex.axis=1.5,
        cex.lab=2) 
addPoints(plotblooms(bloomtoplot2), col=getColors(bloomtoplot2), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(bloomtoplot2[,3]+1)+1)/divide)
legend(x=276, y= 29, 
       legend=c("Background (0-10^3 cells/L)", 
                "Very Low (10^3-10^4 cells/L)", 
                "Low (10^4-10^5 cells/L)",
                "Medium (10^5-10^6 cells/L)", 
                "High (10^6-10^7 cells/L)"), 
       col=c("black","green","yellow","orange","red"), 
       lty=0, lwd=3, pch=c(19),cex=1.5)
scalebar(xy=c(277, 25.25), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1.5)

plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
        ylab = "Latitude",
        main = "November, 2005",
        cex.main=3,
        cex.axis=1.5,
        cex.lab=2)  
addPoints(plotblooms(bloomtoplot3), col=getColors(bloomtoplot3), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(bloomtoplot3[,3]+1)+1)/divide)
legend(x=276, y= 29, 
       legend=c("Background (0-10^3 cells/L)", 
                "Very Low (10^3-10^4 cells/L)", 
                "Low (10^4-10^5 cells/L)",
                "Medium (10^5-10^6 cells/L)", 
                "High (10^6-10^7 cells/L)"), 
       col=c("black","green","yellow","orange","red"), 
       lty=0, lwd=3, pch=c(19),cex=1.5)
scalebar(xy=c(277, 25.25), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1.5)

plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
        ylab = "Latitude",
        main = "November, 2012",
        cex.main=3,
        cex.axis=1.5,
        cex.lab=2) 
addPoints(plotblooms(bloomtoplot4), col=getColors(bloomtoplot4), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(bloomtoplot4[,3]+1)+1)/divide)
legend(x=276, y= 29, 
       legend=c("Background (0-10^3 cells/L)", 
                "Very Low (10^3-10^4 cells/L)", 
                "Low (10^4-10^5 cells/L)",
                "Medium (10^5-10^6 cells/L)", 
                "High (10^6-10^7 cells/L)"), 
       col=c("black","green","yellow","orange","red"), 
       lty=0, lwd=3, pch=c(19),cex=1.5)
scalebar(xy=c(277, 25.25), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1.5)

dev.off()

##not during HAB
tiff(file = "Not During HAB maps.tiff", width =13, height = 15, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(2,2))
plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
        ylab = "Latitude",
        main = "October, 1998",
        cex.main=3,
        cex.axis=1.5,
        cex.lab=2) 
addPoints(plotblooms(nobloomtoplot1), col=getColors(nobloomtoplot1), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(nobloomtoplot1[,3]+1)+1)/divide)
legend(x=276, y= 29, 
       legend=c("Background (0-10^3 cells/L)", 
                "Very Low (10^3-10^4 cells/L)", 
                "Low (10^4-10^5 cells/L)",
                "Medium (10^5-10^6 cells/L)", 
                "High (10^6-10^7 cells/L)"), 
       col=c("black","green","yellow","orange","red"), 
       lty=0, lwd=3, pch=c(19),cex=1.5)
scalebar(xy=c(277, 25.25), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1.5)

plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
        ylab = "Latitude",
        main = "October, 2007",
        cex.main=3,
        cex.axis=1.5,
        cex.lab=2) 
addPoints(plotblooms(nobloomtoplot3), col=getColors(nobloomtoplot3), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(nobloomtoplot3[,3]+1)+1)/divide)
legend(x=276, y= 29, 
       legend=c("Background (0-10^3 cells/L)", 
                "Very Low (10^3-10^4 cells/L)", 
                "Low (10^4-10^5 cells/L)",
                "Medium (10^5-10^6 cells/L)", 
                "High (10^6-10^7 cells/L)"), 
       col=c("black","green","yellow","orange","red"), 
       lty=0, lwd=3, pch=c(19),cex=1.5)
scalebar(xy=c(277, 25.25), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1.5)

plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
        ylab = "Latitude",
        main = "October, 2010",
        cex.main=3,
        cex.axis=1.5,
        cex.lab=2) 
addPoints(plotblooms(nobloomtoplot2), col=getColors(nobloomtoplot2), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(nobloomtoplot2[,3]+1)+1)/divide)
legend(x=276, y= 29, 
       legend=c("Background (0-10^3 cells/L)", 
                "Very Low (10^3-10^4 cells/L)", 
                "Low (10^4-10^5 cells/L)",
                "Medium (10^5-10^6 cells/L)", 
                "High (10^6-10^7 cells/L)"), 
       col=c("black","green","yellow","orange","red"), 
       lty=0, lwd=3, pch=c(19),cex=1.5)
scalebar(xy=c(277, 25.25), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1.5)

plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
        ylab = "Latitude",
        main = "October, 2013",
        cex.main=3,
        cex.axis=1.5,
        cex.lab=2) 
addPoints(plotblooms(nobloomtoplot4), col=getColors(nobloomtoplot4), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(nobloomtoplot4[,3]+1)+1)/divide)
legend(x=276, y= 29, 
       legend=c("Background (0-10^3 cells/L)", 
                "Very Low (10^3-10^4 cells/L)", 
                "Low (10^4-10^5 cells/L)",
                "Medium (10^5-10^6 cells/L)", 
                "High (10^6-10^7 cells/L)"), 
       col=c("black","green","yellow","orange","red"), 
       lty=0, lwd=3, pch=c(19),cex=1.5)
scalebar(xy=c(277, 25.25), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1.5)

dev.off()
###############################################################################################
tiff(file = "bsf1.tiff", width =17, height = 9, units = "in",
     pointsize=10, res = 600, compression = c("lzw"))
par(mfrow=c(1,3),mai=c(3.5,2.5,3.5,2.5))
plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "", 
        ylab = "",
        main = "October, 1998",
        cex.main=5,
        cex.lab=3,
        cex.axis=2)
mtext(side=2,"Latitude",cex=2.5,line=7)
mtext(side=1,"Longitude",cex=2.5,line=13)
addPoints(plotblooms(nobloomtoplot1), col=getColors(nobloomtoplot1), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(nobloomtoplot1[,3]+1)+1)/divide)
legend(x=276.5, y= 28.5, 
       legend=c("Background (0-10^3 cells/L)", 
                "Very Low (10^3-10^4 cells/L)", 
                "Low (10^4-10^5 cells/L)",
                "Medium (10^5-10^6 cells/L)", 
                "High (10^6-10^7 cells/L)"), 
       col=c("black","green","yellow","orange","red"), 
       lty=0, lwd=3, pch=c(19),cex=2)
scalebar(xy=c(276.625, 26.25), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1.5)

plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "", 
        ylab = "",
        main = "January, 2002",
        cex.main=5,
        cex.lab=3,
        cex.axis=2)
mtext(side=1,"Longitude",cex=2.5,line=13)
addPoints(plotblooms(bloomtoplot1), col=getColors(bloomtoplot1), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(bloomtoplot1[,3]+1)+1)/divide)
plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "", 
        ylab = "",
        main = "November, 2005",
        cex.main=5,
        cex.lab=3,
        cex.axis=2)
mtext(side=1,"Longitude",cex=2.5,line=13)
addPoints(plotblooms(bloomtoplot3), col=getColors(bloomtoplot3), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(bloomtoplot3[,3]+1)+1)/divide)
dev.off()


tiff(file = "bsf2.tiff", width =17, height = 9, units = "in",
     pointsize=10, res = 600, compression = c("lzw"))
par(mfrow=c(1,3),mai=c(1,1,1,1))
plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "", 
        ylab = "",
        main = "October, 2010",
        cex.main=5,
        cex.lab=3,
        cex.axis=2)
mtext(side=2,"Latitude",cex=2.5,line=7)
mtext(side=1,"Longitude",cex=2.5,line=13)
addPoints(plotblooms(nobloomtoplot2), col=getColors(nobloomtoplot2), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(nobloomtoplot2[,3]+1)+1)/divide)
scalebar(xy=c(276.625, 26.25), 
         d=100, divs=2, type='bar', 
         lonlat=T, below="km",cex=1.5)
plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "", 
        ylab = "",
        main = "November, 2012",
        cex.main=5,
        cex.lab=3,
        cex.axis=2)
mtext(side=1,"Longitude",cex=2.5,line=13)
addPoints(plotblooms(bloomtoplot4), col=getColors(bloomtoplot4), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(bloomtoplot4[,3]+1)+1)/divide)

plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "", 
        ylab = "",
        main = "January, 2018",
        cex.main=5,
        cex.lab=3,
        cex.axis=2)
mtext(side=1,"Longitude",cex=2.5,line=13)
addPoints(plotblooms(nobloomtoplot4), col=getColors(nobloomtoplot4), pch=16, xlim=c(Wlon, Elon),
          ylim=c(Slat, Nlat), lwd=1,cex=(log10(nobloomtoplot4[,3]+1)+1)/divide)
dev.off()
