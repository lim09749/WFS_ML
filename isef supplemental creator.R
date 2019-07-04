require(e1071)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Get Data
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
#####################################################################################
##Read in entire dataset
alldata = read.csv("alldata.csv")
alldata$X = NULL
col_labels = c("NS wind at venf1 (m/s)","EW wind at venf1 (m/s)","ATMP at venf1 (Celsius)",
               "NS wind at cdrf1 (m/s)","EW wind at cdrf1 (m/s)","NS wind at 42039 (m/s)",
               "EW wind at 42039 (m/s)","NS wind at 42003 (m/s)","EW wind at 42003 (m/s)",
               "Tampa (cfs)","Peace (cfs)","Okeechobee (cfs)","Suwanee (cfs)","Withlacoochee (cfs)",
               "Hillsborough (cfs)",
               "Manatee (cfs)","Myakka (cfs)","Tampa TN (log10(mg/week))", "Myakka TN (log10(mg/week))",
               "Peace TN (log10(mg/week))","Withlacoochee TN (log10(mg/week))","Manatee TN (log10(mg/week))",
               "Hillsborough TN (log10(mg/week))","Tampa TP (log10(mg/week))", "Myakka TP (log10(mg/week))",
               "Peace TP (log10(mg/week))","Withlacoochee TP (log10(mg/week))","Manatee TP (log10(mg/week))",
               "Hillsborough TP (log10(mg/week))","SSH (meters)","Caloosahatchee (cfs)", "Caloosahatchee TN (log10(mg/week))",
               "Caloosahatchee TP (log10(mg/week))")
#####################################################################################
##Histograms
for(i in 1:11){
  tiff(paste0(i,"_hist.tiff"), width = 8.5, height = 11, units = "in",
       pointsize=10, res = 300, compression = c("lzw"))
  par(mfrow=c(3,2),mar=c(5,5,5,5)+0.1)
  for(col in (3*(i-1)+1):(3*i)){
    hist(alldata[states==0,col],main=col_labels[col],xlab="",cex.main=2,cex.axis=2,cex.lab=2)
    hist(alldata[states==1,col],main=col_labels[col],xlab="",col="red",cex.main=2,cex.axis=2,cex.lab=2)
  }
  dev.off()
}
#####################################################################################
##Time series
for(i in 1:11){
  tiff(paste0(i,"_time_series.tiff"), width = 8.5, height = 11, units = "in",
       pointsize=10, res = 300, compression = c("lzw"))
  par(mfrow=c(3,1),mar=c(5,5,5,5)+0.1)
  for(col in (3*(i-1)+1):(3*i)){
    colors = rep("black",nrow(alldata))
    colors[states==1]="red"
    plot(as.Date(HAB$Date),alldata[,col],main=col_labels[col],xlab="Years",
         cex.main=2,cex.axis=2,cex.lab=2,col=colors,pch=24,ylab=col_labels[col])
  }
  dev.off()
}
#####################################################################################
##Scatter
for(i in 1:11){
  tiff(paste0(i,"_scatter.tiff"), width = 8.5, height = 11, units = "in",
       pointsize=10, res = 300, compression = c("lzw"))
  par(mfrow=c(3,1),mar=c(5,5,5,5)+0.1)
  for(col in (3*(i-1)+1):(3*i)){
    colors = rep("black",nrow(alldata))
    colors[states==1]="red"
    plot(alldata[,col],states,main=col_labels[col],xlab=col_labels[col],
         cex.main=2,cex.axis=2,cex.lab=2,col=colors,ylab="HAB or no HAB")
  }
  dev.off()
}
####################################################################3
##Statistics
vals = rep(NA, 33)
for(col in 1:33){
  vals[col]=t.test(x=alldata[states==0,col],y=alldata[states==1,col])$p.value
}
write.csv(cbind(col_labels,vals),file="pvalue_mat.csv")
####################################################################3
##Map/bubble plots
require(PBSmapping)
require(sp)
require(raster)

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

suffix = c("-01-01","-02-01","-03-01","-04-01","-05-01","-06-01",
  "-07-01","-08-01","-09-01","-10-01","-11-01","-12-01")
years = 1999:2018
totallength = c()
for(y in years){
  totallength=c(totallength,paste0(y,suffix))
}

toprint = as.Date(totallength[seq(from=1,to=length(totallength),by=6)])

#################################################################################################
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
######################################################################################################3
for(i in 1:(length(toprint)/4)){
  tiff(paste0(i,"_bubble_plot.tiff"), width = 8.5, height = 11, units = "in",
       pointsize=10, res = 300, compression = c("lzw"))
  par(mfrow=c(2,2))
  for(j in (4*(i-1)+1):(4*i)){
    plotMap(zoomMap, col="floralwhite", bg="lightblue1", xlab = "Longitude", 
            ylab = "Latitude",
            main = "January, 2002",
            cex.main=3,
            cex.axis=1.5,
            cex.lab=2) 
    sampleloc = cbind(HABtotal$Latitude,HABtotal$Longitude, HABtotal$Abundance_cells)[HABtotal$Date>=toprint[j] & HABtotal$Date<toprint[j]+30,]
    
    addPoints(plotblooms(sampleloc), col=getColors(sampleloc), pch=16, xlim=c(Wlon, Elon),
              ylim=c(Slat, Nlat), lwd=1,cex=(log10(sampleloc[,3]+1)+1)/divide)
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
    
  }
  dev.off()
}
######################################################################################################3
##fertilizer plots
n_fert = read.csv("Global Nitrogen Fertilizer.csv")
p_fert = read.csv("Global Phosphorus Fertilizer.csv")
n_fert$Date = as.numeric(substr(as.character(n_fert$Data),start=8,stop=nchar(as.character(n_fert$Data))))
p_fert$Date = as.numeric(substr(as.character(p_fert$Data),start=8,stop=nchar(as.character(p_fert$Data))))
n_fert$Data = NULL
p_fert$Data= NULL
colnames(n_fert)
#####################################################################################################3
##Start Plotting!
tiff("Global Fertilizer Use.tiff", width = 8.5, height = 11, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))

par(mfrow=c(2,1))

##########################################################################################################
##Nitrogen
n_fert$North_America = n_fert$Canada+n_fert$USA
n_fert$Europe = n_fert$Central.Europe+n_fert$Western.Europe+n_fert$Ukraine+n_fert$Russia..
n_fert$Central_South_America = n_fert$Mexico+n_fert$Rest.Central.America+n_fert$Brazil+n_fert$Rest.South.America
n_fert$Africa = n_fert$Northern.Africa+n_fert$Western.Africa+n_fert$Eastern.Africa+n_fert$Southern.Africa
n_fert$Eur_Asia=n_fert$Asia.Stan+n_fert$Middle.East
n_fert$South_Asia = n_fert$Southeastern.Asia+n_fert$India..+n_fert$Indonesia..
n_fert$East_Asia = n_fert$Korea+n_fert$China..
n_fert$Australia = n_fert$Oceania

names = c("North America","Europe","Central and South America","Africa",
          "Eurasia","South Asia","East Asia","Australia")
n_clean = n_fert[,29:36]
##Fill in underneath
fillin = function(Date, Val,color_fill,density){
  for(i in 1:(length(Date)-1)){
    points(seq(Date[i],Date[i+1],length.out=density),
           seq(Val[i],Val[i+1],length.out=density),col=color_fill,type="h")
  }
}


##Start plotting
plot(n_fert$Date,(rowSums(n_clean)),type="h",ylim=c(0,1.4e+08),xlab="Year",ylab="N Applied")
colors=c("red","blue","green","cyan","orange","purple","pink","yellow")
for(i in 1:ncol(n_clean)){
  if(i!=ncol(n_clean)){
    fillin(n_fert$Date,(rowSums(n_clean[,i:ncol(n_clean)])),density=100,color_fill=colors[i])
  }else{
    fillin(n_fert$Date,(n_clean[,i]),density=100,color_fill=colors[i])
  }
}
legend(1969,1.4e+08,
       legend=names, 
       col=colors, 
       lty=0, lwd=3, pch=c(19),cex=1.25)
#####################################################################################################
p_fert$North_America = p_fert$Canada+p_fert$USA
p_fert$Europe = p_fert$Central.Europe+p_fert$Western.Europe+p_fert$Ukraine+p_fert$Russia..
p_fert$Central_South_America = p_fert$Mexico+p_fert$Rest.Central.America+p_fert$Brazil+p_fert$Rest.South.America
p_fert$Africa = p_fert$Northern.Africa+p_fert$Western.Africa+p_fert$Eastern.Africa+p_fert$Southern.Africa
p_fert$Eur_Asia=p_fert$Asia.Stan+p_fert$Middle.East
p_fert$South_Asia = p_fert$Southeastern.Asia+p_fert$India..+p_fert$Indonesia..
p_fert$East_Asia = p_fert$Korea+p_fert$China..
p_fert$Australia = p_fert$Oceania

p_clean = p_fert[,29:36]
##Fill in underneath
fillin = function(Date, Val,color_fill,density){
  for(i in 1:(length(Date)-1)){
    points(seq(Date[i],Date[i+1],length.out=density),
           seq(Val[i],Val[i+1],length.out=density),col=color_fill,type="h")
  }
}


##Start plotting
plot(p_fert$Date,(rowSums(p_clean)),type="h",ylim=c(0,4.0e+7),xlab="Year",ylab="P Applied")
colors=c("red","blue","green","cyan","orange","purple","pink","yellow")
for(i in 1:ncol(p_clean)){
  if(i!=ncol(p_clean)){
    fillin(p_fert$Date,(rowSums(p_clean[,i:ncol(p_clean)])),density=100,color_fill=colors[i])
  }else{
    fillin(p_fert$Date,(p_clean[,i]),density=100,color_fill=colors[i])
  }
}
dev.off()