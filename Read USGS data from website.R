###########################################################################################
## "Read USGS data from website.R"
##Read in all Florida station data with relevant values
###########################################################################################
rm(list = ls())
setwd("C:/Users/xswang/HAB Research/USGS raw data")
###########################################################################################
##Read in data
options(scipen=999)

stations = read.csv("C:/Users/xswang/HAB Research/all florida stations.csv")
stations[,1] = as.character(stations[,1])
stations[,2] = as.character(stations[,2])
stations[,3] = as.character(stations[,3])

##Pad spacing
for(i in 1:nrow(stations)){
  if(nchar(stations[i,2])<8){
    for(l in 1:(8-nchar(stations[i,2]))){
      stations[i,2] = paste0("0",stations[i,2])
    }
  }
}

##Prefix to get station urls
prefix   = "https://waterservices.usgs.gov/nwis/site/?format=rdb&sites="
URLs     = paste0(prefix, stations$site_no)
###########################################################################################
##Read in latitude and longitude data
lat = rep(NA, length(URLs))
lon = rep(NA, length(URLs))
for(i in 1:length(URLs)){
  try(  {
    page       = read.csv(url(URLs[i]))
    ##find the row in the page and split it with the tab
    s          = strsplit(as.character(page[nrow(page),1]), split = "\t")[[1]]
    ##Convert to numeric to find latitude and longitude
    s          = na.omit(as.numeric(s))
    ##Trick to find index of latitude measurement
    latlon     = na.omit(as.numeric(s))[min(which(s>10 & s<40)):(min(which(s>10 & s<40))+1)]
    ##Save the data
    lat[i]     = latlon[1]
    lon[i]     = latlon[2]
    print(latlon)
  }
  )
  # ##Get page
  # page       = read.csv(url(URLs[i]))
  # ##find the row in the page and split it with the tab
  # s          = strsplit(as.character(page[nrow(page),1]), split = "\t")[[1]]
  # ##Convert to numeric to find latitude and longitude
  # s          = na.omit(as.numeric(s))
  # ##Trick to find index of latitude measurement
  # latlon     = na.omit(as.numeric(s))[min(which(s>10 & s<40)):(min(which(s>10 & s<40))+1)]
  # ##Save the data
  # lat[i]     = latlon[1]
  # lon[i]     = latlon[2]
}
coord = as.data.frame(cbind(lat,lon))
colnames(coord) = c("lat","lon")
write.csv(coord,paste0("C:/Users/xswang/HAB Research/USGS coord.csv"))
###########################################################################################
##Find minimum and maximum of in-situ data
HAB       = read.csv("C:/Users/xswang/HAB Research/Florida_HAB.csv")
HAB$Date  = as.Date(as.character(HAB$Date),format="%m/%d/%Y")

startDT   = min(HAB$Date)
endDT     = max(HAB$Date)
endPoints = as.character(c(startDT, seq(from=min(HAB$Date)+365, by=365, length.out=19), endDT))
###########################################################################################
##Find data
for(e in 1:(length(endPoints)-1)){
  for(i in 1:nrow(stations)){
    link  = paste0("https://waterservices.usgs.gov/nwis/",
                   "iv/?format=rdb&sites=",stations$site_no[i],
                   "&startDT=",endPoints[e],"&endDT=",endPoints[e+1])
    try(  {
      page       = read.csv(url(link))
      page       = as.character(page[,1])
      if(length(page)>5 && !paste0("USGS ",stations$site_no[i]," from ", endPoints[e]," to ",endPoints[e+1],".csv")%in%list.files()){
        begin      = max(which(substr(page,0,1)=="#"))+1
        write.csv(page[1:(begin-1)],paste0("USGS ",stations$site_no[i],
                                           " Heading from ", endPoints[e]," to ",endPoints[e+1],".csv"))
        page       = page[c(begin,(begin+2):length(page))]
        df                        = as.data.frame(matrix(NA, nrow = length(page)-1, 
                                                         ncol = length(strsplit(page[1],split="\t")[[1]])))
        colnames(df)              = strsplit(page[1],split="\t")[[1]]
        df[,1:ncol(df)]           = do.call(rbind, strsplit(page[2:length(page)], split = "\t"))
        write.csv(df, paste0("USGS ",stations$site_no[i]," from ", endPoints[e]," to ",endPoints[e+1],".csv"))
      }
    })
  }
  print(e)
}
###########################################################################################
