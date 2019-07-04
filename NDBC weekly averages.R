########################################################################################
## "NDBC weekly averages.R"
## Averages wind data by week
########################################################################################
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/NDBC daily data/")
files = list.files()
########################################################################################
##Create weekly data
for(i in 1:length(files)){
  info = file.info(files[i])
  if(info$size>10){
    wind = read.csv(files[i])
    wind$X = NULL
    ##Remove all na values
    wind = wind[,colSums(is.na(wind))<nrow(wind)]
    
    ##Remove all NA values
    nacol = c()
    for(j in 1:ncol(wind)){
      if(all(is.na(wind[,j]))) nacol = c(nacol,j)
    }
    if(length(nacol)>0) wind = wind[,-nacol]
    
    ##Get weekly intervals
    datetime = as.Date(as.character(wind$day))
    weeklyintervals = seq(as.Date("1998-01-01"),max(datetime),by=7)
    
    ##Create used db
    db = as.data.frame(matrix(NA, nrow=length(weeklyintervals),ncol=ncol(wind)))
    colnames(db) = colnames(wind)
    regular_col = setdiff(colnames(wind)[2:ncol(wind)],c("WSPD","WDIR"))
    db[,1] = weeklyintervals
    ##Put in new data
    for(w in 1:length(weeklyintervals)){
      rightrows  = datetime>=weeklyintervals[w] & datetime < (weeklyintervals[w]+7)
      if(length(which(rightrows))>0){
        ##Find daily averages of windspeed and direction
        # u  = mean(cos(wind[rightrows,]$WDIR*pi/180)*wind[rightrows,]$WSPD,na.rm=T)
        # v  = mean(sin(wind[rightrows,]$WDIR*pi/180)*wind[rightrows,]$WSPD,na.rm=T)
        # 
        # db[w,] = c(as.character(weeklyintervals[w]), sqrt(u^2+v^2), atan2(v,u), 
        #            colMeans(wind[rightrows,match(regular_col,colnames(wind))],na.rm=T))
        angle = atan2(mean(cos(wind[rightrows,]$WDIR*pi/180), na.rm= T), ##in degrees
                      mean(sin(wind[rightrows,]$WDIR*pi/180), na.rm= T))*180/pi
        len = mean(wind[rightrows,]$WSPD, na.rm=T) 
        db[w, ] = c(as.character(weeklyintervals[w]), len, angle, 
                     colMeans(wind[rightrows,match(regular_col,colnames(wind))],na.rm=T))
        
      }
    }
    
    write.csv(db, paste0("C:/Users/xswang/HAB Research/NDBC weekly data/",
                          strsplit(files[i],split=" daily.csv")[[1]][1], " weekly.csv"))
  }
  print(i)
}
