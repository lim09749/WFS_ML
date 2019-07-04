require(quantmod) #just to find peaks
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Read in data
HAB = read.csv("C:/Users/xswang/HAB Research/West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)

ndbc = read.csv("NDBC_station_data.csv")
ndbc = ndbc[,!is.na(colSums(ndbc))]
ndbc$X = NULL
usgs = read.csv("USGS_discharge_total.csv")
usgs$X = NULL

#######################################################
tiff(file = "USGS KB time series.tiff", width = 7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(3,1))
plot(HAB$Date, log10(HAB$Abundance_cells),type="l",ylim=c(2,8), col = "red",
     xlab = "Date",
     ylab="Kb (log10(c/l))",
     main = "Karenia brevis abundance in West Florida Shelf")
abline(a=5,b=0)
hgreater = which(usgs$Hillsborough>(mean(usgs$Hillsborough)+1*sd(usgs$Hillsborough)))
ogreater = which(usgs$Okeechobee> (mean(usgs$Okeechobee)+1*sd(usgs$Okeechobee)))
abline(v=HAB$Date[hgreater],
       col="blue")
abline(v=HAB$Date[ogreater],
       col="green")
points(HAB$Date[log10(HAB$Abundance_cells)>=5],
       log10(HAB$Abundance_cells)[log10(HAB$Abundance_cells)>=5])
       #col="black")
plot(HAB$Date, log10(usgs$Hillsborough), type = "l",
     xlab = "Date",
     ylab="Discharge (log10(cfs))",
     col="blue",
     main = "Hillsborough River")
abline(a=log10(sd(usgs$Hillsborough)+mean(usgs$Hillsborough)), b=0)
points(HAB$Date[usgs$Hillsborough>(sd(usgs$Hillsborough)+mean(usgs$Hillsborough))],
       log10(usgs$Hillsborough[usgs$Hillsborough>(sd(usgs$Hillsborough)+mean(usgs$Hillsborough))]))
plot(HAB$Date, log10(usgs$Okeechobee), type = "l",
     xlab = "Date",
     ylab="Discharge (log10(cfs))",
     col="green",
     main = "Lake Okeechobee")
abline(a=log10(sd(usgs$Okeechobee)+mean(usgs$Okeechobee)), b=0)
points(HAB$Date[usgs$Okeechobee>(sd(usgs$Okeechobee)+mean(usgs$Okeechobee))],
       log10(usgs$Okeechobee[usgs$Okeechobee>(sd(usgs$Okeechobee)+mean(usgs$Okeechobee))]))
dev.off()
#######################################################

tiff(file = "NDBC KB time series.tiff", width = 7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(3,1))
plot(HAB$Date, log10(HAB$Abundance_cells),type="l",ylim=c(2,8), col = "red",
     xlab = "Date",
     ylab="Kb (log10(c/l))",
     main = "Karenia brevis abundance in West Florida Shelf")
abline(a=5,b=0)
hgreater = which(ndbc$X42039.NS>(mean(ndbc$X42039.NS)+1.5*sd(ndbc$X42039.NS)))
ogreater = which(ndbc$venf1.NS> (mean(ndbc$venf1.NS)+1.5*sd(ndbc$venf1.NS)))
abline(v=HAB$Date[hgreater],
       col="blue")
abline(v=HAB$Date[ogreater],
       col="green")
points(HAB$Date[log10(HAB$Abundance_cells)>=5],
       log10(HAB$Abundance_cells)[log10(HAB$Abundance_cells)>=5])
#col="black")
plot(HAB$Date, (ndbc$X42039.NS), type = "l",
     xlab = "Date",
     ylab="Wind Speed (m/s)",
     col="blue",
     main = "NS wind from Station 42039")
abline(a=(1.5*sd(ndbc$X42039.NS)+mean(ndbc$X42039.NS)), b=0)
points(HAB$Date[ndbc$X42039.NS>(1.5*sd(ndbc$X42039.NS)+mean(ndbc$X42039.NS))],
       (ndbc$X42039.NS[ndbc$X42039.NS>(1.5*sd(ndbc$X42039.NS)+mean(ndbc$X42039.NS))]))
plot(HAB$Date, (ndbc$venf1.NS), type = "l",
     xlab = "Date",
     ylab="Wind Speed (m/s)",
     col="green",
     main = "NS wind from venf1")
abline(a=(1.5*sd(ndbc$venf1.NS)+mean(ndbc$venf1.NS)), b=0)
points(HAB$Date[ndbc$venf1.NS>(1.5*sd(ndbc$venf1.NS)+mean(ndbc$venf1.NS))],
       (ndbc$venf1.NS[ndbc$venf1.NS>(1.5*sd(ndbc$venf1.NS)+mean(ndbc$venf1.NS))]))
dev.off()

