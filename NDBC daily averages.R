########################################################################################
## "NDBC daily averages.R"
## Averages wind data
########################################################################################
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/NDBC raw data/")
files = list.files()
########################################################################################
##Test example
for(i in 1:length(files)){
  info = file.info(files[i])
  if(info$size>10){
    wind = read.csv(files[i])
    wind$X = NULL
    
    ##Make sure to put old variable names into the new ones
    oldnames = c("WD","DIR","SPD","GSP","GMN","BAR","H0","DOMPD","AVP","SRAD","SRAD2",
                 "LRAD","LRAD1")
    newnames = c("WDIR", "WDIR","WSPD","GST","GTIME","PRES",
                 "WVHT","DPD","APD","SWRAD","SWRAD","LWRAD","LWRAD")
    
    for(j in 1:length(oldnames)){
      old = oldnames[j]
      if(old %in% colnames(wind)){
        oldcolval = which(colnames(wind)==old)
        newcolval = which(colnames(wind)==newnames[j])
        minrows   = min(which(is.na(wind[,oldcolval])))-1
        wind[1:minrows, newcolval] = wind[1:minrows, oldcolval]
        wind[,oldcolval] = NULL
      }
    }
    
    ##update years
    wind$YYYY[wind$YYYY<100] = wind$YYYY[wind$YYYY<100]+1900
    winddays = as.Date(paste0(wind$YYYY,"-",wind$MM,"-",wind$DD))
  
    for(j in 1:ncol(wind)){
      if(!is.numeric(wind[,j])){
        wind[,j] = as.numeric(as.character(wind[,j]))
      }
    }
    wind[wind==99 | wind==999 | wind==9999] = NA

    
    ##Create a df with the daily values
    mat = data.frame(matrix(NA, nrow = length(unique(winddays)), ncol = ncol(wind)-3))
    regular_col = setdiff(colnames(wind)[5:ncol(wind)],c("WSPD","WDIR"))
    colnames(mat) = c("day","WSPD","WDIR",colnames(wind)[match(regular_col,colnames(wind))])
    for(j in 1:length(unique(winddays))){
      day  = unique(winddays)[j]
      rows = which(winddays==day)
      
      ##Find daily averages of windspeed and direction
      # u  = mean(cos(wind[rows,]$WDIR*pi/180)*wind[rows,]$WSPD,na.rm=T)
      # v  = mean(sin(wind[rows,]$WDIR*pi/180)*wind[rows,]$WSPD,na.rm=T)
      # mat[j,] = c(as.character(day), sqrt(u^2+v^2), atan2(v,u), 
      #             colMeans(wind[rows,match(regular_col,colnames(wind))],na.rm=T))
      # 
      angle = atan2(mean(cos(wind[rows,]$WDIR*pi/180), na.rm= T), ##in degrees
                    mean(sin(wind[rows,]$WDIR*pi/180), na.rm= T))*180/pi
      len = mean(wind[rows,]$WSPD, na.rm=T) 
      mat[j, ] = c(as.character(day), len, angle, 
                   colMeans(wind[rows,match(regular_col,colnames(wind))],na.rm=T))
    }
    mat = na.omit(mat)
    write.csv(mat, paste0("C:/Users/xswang/HAB Research/NDBC daily data/",
                          strsplit(files[i],split=" wind ")[[1]][1], " daily.csv"))
  }
  print(i)
}
