###########################################################################################
## "USGS convert to daily.R"
## Convert hourly data to daily data
###########################################################################################
rm(list = ls())
setwd("C:/Users/xswang/HAB Research/USGS raw total")
###########################################################################################
files = list.files()
for(i in 1:length(files)){
  dfval       = read.csv(files[i])
  dfval$X     = NULL
  
  for(j in 1:ncol(dfval)){
    if(j <= 4) dfval[,j] = as.character(dfval[,j])
    else dfval[,j] = as.numeric(dfval[,j])
  }
  
  dfval$tz_cd[dfval$tz_cd=="EDT"] = "EST"
  datetime = as.POSIXct(dfval$datetime, tz = dfval$tz_cd[1], format = "%Y-%m-%d %H:%M")
  attributes(datetime)$tzone = "EST"
  
  ##Get days
  days     <- as.character(unlist(strsplit(dfval$datetime, split = " "))[seq(from=1,to=2*length(datetime),by=2)])
  uniqdays <- as.character(unique(days))
  
  tempStations <- data.frame(matrix(NA, nrow = length(uniqdays), ncol = ncol(dfval)))
  colnames(tempStations) <- colnames(dfval)
  
  for(j in 1:length(uniqdays)){
    d <- uniqdays[j]
    if(ncol(dfval)==5) arr <- dfval[which(d==days),5:ncol(dfval)]/length(which(days==d))
    else arr <- colSums(dfval[which(d==days),5:ncol(dfval)], na.rm=T)/length(which(days==d))
    arr[arr==0] = NA
    tempStations[j,] <- c(as.character(c(dfval$agency_cd[j],dfval$site_no[j],
                                         uniqdays[j],"EST")),
                          arr)
  }
  write.csv(tempStations,paste0("C:/Users/xswang/HAB Research/USGS daily data/",strsplit(files[i],split=".csv")[[1]], " daily.csv"))
  if(i%%10==0) print(i)
}

