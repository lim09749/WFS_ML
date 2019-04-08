rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#################################################################################################
##Get HAB
NDBC = read.csv("NDBC_lat_lon.csv")
NDBC$X = NULL
NDBC$Latitude = as.numeric(NDBC$Latitude)
NDBC$Longitude = -as.numeric(NDBC$Longitude)

HAB = read.csv("C:/Users/xswang/HAB Research/West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
# windstations = c("ktnf1","41009",
#                  "spgf1","fwyf1","mlrf1",
#                  "lonf1","smkf1","venf1",
#                  "cdrf1","42036","42039",
#                  "42003")

windstations = c("venf1","cdrf1","42039","42003")
#"42036",

files = paste0("C:/Users/xswang/HAB Research/NDBC weekly data/", windstations, " weekly.csv")

#################################################################################################
dframe = as.data.frame(matrix(NA, nrow=nrow(HAB), ncol = 13))
for(j in 1:length(files)){
  f=files[j]
  r = read.csv(f)
  r$X = NULL
  r$day = as.Date(as.character(r$day))
  print(f)
  r = r[,match(c("day","WSPD","WDIR","ATMP"),colnames(r))]
  print(colnames(r))
  magnitude = as.numeric(as.character(r[,2]))
  angle = as.numeric(as.character(r[,3]))
  if(length(which(is.na(magnitude)))<=100){
    if(length(angle) < length(states)) {
      magnitude = c(magnitude, rep(NA, length(states)-length(magnitude)))
      magnitude[is.na(magnitude)] = mean(magnitude,na.rm=T)
      angle = c(angle, rep(NA, length(states)-length(magnitude)))
      angle[is.na(magnitude)] = mean(angle,na.rm=T)
    }
    dframe[(j-1)*3+1] = magnitude*(cos(angle*pi/180))
    dframe[(j-1)*3+2] = magnitude*(sin(angle*pi/180))
    colnames(dframe)[(j-1)*3+1:2] = paste0(windstations[j],c(" NS", " EW"))
    
    # dframe[(j-1)*3+1] = magnitude*(cos(angle))
    # dframe[(j-1)*3+2] = magnitude*(sin(angle))
    # colnames(dframe)[(j-1)*3+1:2] = paste0(windstations[j],c(" NS", " EW"))
  }
  if("ATMP"%in%colnames(r)){
    ATMP  = as.numeric(as.character(r$ATMP))
    if(length(which(is.na(ATMP))) <=100){
      if(length(ATMP) < length(states)) {
        ATMP = c(ATMP, rep(NA, length(states)-length(ATMP)))
        ATMP[is.na(ATMP)] = mean(ATMP,na.rm=T)
      }
      dframe[(j-1)*3+3] = ATMP
      colnames(dframe)[(j-1)*3+3] = paste0(windstations[j]," ATMP")
    }
  }
}
dframe = dframe[,c(1,2,3,4,5,7,8,10,11)]
for(i in 1:(ncol(dframe))) dframe[is.na(dframe[,i]),i] = mean(dframe[,i],na.rm=T)
write.csv(dframe,"NDBC_station_data.csv")
