rm(list=ls())
setwd("C:/Users/xswang/HAB Research/USGS weekly data")
#######################################################
##Get names
USGSlocs = read.csv("C:/Users/xswang/HAB Research/all florida stations.csv")
USGSlocs$agency_cd = as.character(USGSlocs$agency_cd)
USGSlocs$site_no = as.character(USGSlocs$site_no)
USGSlocs$station_nm = as.character(USGSlocs$station_nm)
#######################################################
rivers = c("SUWANNEE","WACCASASSA","WITHLACOOCHEE","TAMPA",
           "PEACE","CHARLOTTE","CALOOSAHATCHEE","OKEECHOBEE",
           "HILLSBOROUGH","LITTLE MANATEE","MANATEE", "MYAKKA")
##Proxies for Caloosahatchee river
rivers = c("SAN CARLOS","MYERS","CAPE CORAL", "TARPON")
files = list.files()
for(string in rivers){
  print(string)
  stationids = (USGSlocs[which(grepl(string,USGSlocs$station_nm)),])$site_no
  for(id in stationids){
    if(any(grepl(id,files))){
      r = read.csv(files[grepl(id,files)][which.max(nchar(files[grepl(id,files)]))])
      r$X = NULL
      vec = c()
      for(i in 2:ncol(r)){
        vec = c(vec,length(which(!is.na(r[,i]))))
      }
      print(id)
      print((USGSlocs[which(grepl(string,USGSlocs$station_nm)),])$station_nm[stationids==id])
      print(colnames(r[,2:ncol(r)]))
      print(vec)
    }
  }
}


