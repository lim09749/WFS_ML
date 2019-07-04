rm(list=ls())
######################################################################################################################
##Read in data
setwd("C:/Users/xswang/HAB Research")
rivers = c("Tampa",
           "Myakka",
           "Peace",
           "Withlacoochee",
           "Manatee",
           "Hillsborough",
           "Caloosahatchee")
TN_files = paste0(rivers,"_TN.csv")
TP_files = paste0(rivers,"_TP.csv")
######################################################################################################################
par(mfrow=c(4,3))
for(i in 1:length(TN_files)){
  r = read.csv(TN_files[i])
  r$X = NULL
  dates = as.Date(unlist(strsplit(as.character(r$SampleDate),split=" "))[seq(from=1,to=1+2*(nrow(r)-1),by=2)],format="%m/%d/%Y")
  
  plot(dates, as.numeric(as.character(r$Result_Value)),
       main=rivers[i])
}

for(i in 1:length(TP_files)){
  r = read.csv(TP_files[i])
  r$X = NULL
  dates = as.Date(unlist(strsplit(as.character(r$SampleDate),split=" "))[seq(from=1,to=1+2*(nrow(r)-1),by=2)],format="%m/%d/%Y")
  
  plot(dates, as.numeric(as.character(r$Result_Value)),
       main=rivers[i])
}

