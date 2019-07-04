#######################################################################
##"weekly USGS features.R"
##Convert daily to weekly USGS data
#######################################################################
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/USGS daily features/")
#######################################################################
files = list.files()
for(f in files){
  r = read.csv(f)
  r$X <- NULL
  feature = substr(f,0,nchar(f)-4)
  datetime = as.Date(as.character(r$date))
  weeklyintervals = seq(as.Date("1998-01-01"),max(datetime),by=7)
  
  db = as.data.frame(matrix(NA, nrow=length(weeklyintervals),ncol=ncol(r)))
  colnames(db) = colnames(r)
  db$date = weeklyintervals
  
  for(w in 1:length(weeklyintervals)){
    rightrows  = datetime>=weeklyintervals[w] & datetime < (weeklyintervals[w]+7)
    if(ncol(db)==2) db[w,2] = NA
    else{
      db[w,2:ncol(db)]    = colMeans(r[rightrows,2:ncol(db)])
    }
  }
  
  write.csv(db, paste0("C:/Users/xswang/HAB Research/USGS weekly features/",
                       feature,"_weekly_avg.csv"))
}