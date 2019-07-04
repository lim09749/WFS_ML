#######################################################################
##Convert daily to weekly USGS data
#######################################################################
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/USGS daily data/")
#######################################################################
files = list.files()
for(f in files){
  r = read.csv(f)
  r$X <- NULL
  datetime = as.Date(as.character(r$date))
  weeklyintervals = seq(as.Date("1998-01-01"),max(datetime),by=7)
  
  db = as.data.frame(matrix(NA, nrow=length(weeklyintervals),ncol=ncol(r)-3))
  colnames(db) = c("date",colnames(r[,5:ncol(r)]))
  db$date = weeklyintervals
  
  for(w in 1:length(weeklyintervals)){
    rightrows  = datetime>=weeklyintervals[w] & datetime < (weeklyintervals[w]+7)
    if(ncol(db)==2) db[w,2] = NA
    else{
      db[w,2:ncol(db)]    = colMeans(r[rightrows,5:ncol(r)])
    }
  }
  
  write.csv(db, paste0("C:/Users/xswang/HAB Research/USGS weekly data/",
                       strsplit(f,split=" ")[[1]][2],"_weekly.csv"))
}
