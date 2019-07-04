rm(list=ls())
setwd("C:/Users/xswang/HAB Research/USGS weekly data")
#######################################################
##Get HAB
HAB = read.csv("C:/Users/xswang/HAB Research/West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
USGSlocs = read.csv("C:/Users/xswang/HAB Research/all florida stations.csv")
USGSlocs$agency_cd = as.character(USGSlocs$agency_cd)
USGSlocs$site_no = as.character(USGSlocs$site_no)
USGSlocs$station_nm = as.character(USGSlocs$station_nm)

#######################################################
plotRiver = function(stations){
  mat = matrix(NA, nrow=nrow(HAB),ncol=length(stations))
  files = list.files()
  for(i in 1:length(stations)){
    id = stations[i]
    r = read.csv(files[grepl(id,files)][which.max(nchar(files[grepl(id,files)]))])
    r$X = NULL
    discharge = r$Discharge.cubic.feet.per.second
    if(nrow(r) < length(states))   discharge = c(discharge, rep(NA, length(states)-nrow(r)))
    #plot(discharge[!is.na(discharge)], states[!is.na(discharge)])
    mat[,i] = discharge
    mat[is.na(discharge),i]=mean(discharge,na.rm=T)
  }
  plot(HAB$Date, states)
  #plot(HAB$Date, HAB$Abundance_cells/max(HAB$Abundance_cells, na.rm=T))
  peaks = findPeaks(rowSums(mat))
  peaks = peaks[rowSums(mat)[peaks]>=2*mean(rowSums(mat))]
  lines(HAB$Date, rowSums(mat)/max(rowSums(mat)))
  points(HAB$Date[peaks],
         (rowSums(mat)/max(rowSums(mat)))[peaks],col="red")
  abline(h=mean(rowSums(mat))*2/max(rowSums(mat)))
  abline(v=HAB$Date[peaks])
  return(mat)
}

###################################################################################
par(mfrow=c(2,1))
require(quantmod)

USGSinput = matrix(data = NA, nrow = 1083, ncol = 8)

###################################################################################
##Get Tampa station data 
stations = c("2301740",
             "2301738",
             "2301745",
             "2301750",
             "2301793",
             "2306647",
             "2306654")
mat = plotRiver(stations)
USGSinput[,1] = mat[,which.max(colSums(mat,na.rm=T))]
print(stations[which.max(colSums(mat,na.rm=T))])
title("Tampa Bay")
###################################################################################

##Peace discharge
stations = c("2293987",
             "2294650",
             "2294655",
             "2294898",
             "2295637",
             "2296750")
mat = plotRiver(stations)
USGSinput[,2] = mat[,which.max(colSums(mat,na.rm=T))]
print(stations[which.max(colSums(mat,na.rm=T))])
title("Peace River")
###################################################################################

##Charlotte (only gage height)
# stations = c("2299230")
# plotRiver(stations)
#title("Charlotte River")

###################################################################################

##Caloosahatchee
# stations = c("2292010")
# plotRiver(stations)
# title("Caloosahatchee River")

###################################################################################

##Okeechobee
stations = c("2274005",
             "2274325",
             "2274490",
             "2274505")
mat = plotRiver(stations)
USGSinput[,3] = mat[,which.max(colSums(mat,na.rm=T))]
print(stations[which.max(colSums(mat,na.rm=T))])
title("Okeechobee River")

###################################################################################

##Suwannee
stations = c("2315500",
             "2319500",
             "2320000",
             "2319800",
             "2320500",
             "2323500")
mat = plotRiver(stations)
USGSinput[,4] = mat[,which.max(colSums(mat,na.rm=T))]
print(stations[which.max(colSums(mat,na.rm=T))])
title("Suwannee River")

###################################################################################

##Waccasassa
# stations = c("2313700")
# 
# plotRiver(stations)
# title("Waccasassa River")

###################################################################################

##Withlacoochee
stations = c("2310947",
             "2313000",
             "2311500",
             "2319000",
             "2312000",
             "2312200",
             "2312600",
             "2312720")
mat = plotRiver(stations)
USGSinput[,5] = mat[,which.max(colSums(mat,na.rm=T))]
print(stations[which.max(colSums(mat,na.rm=T))])
title("Withlacoochee River")

###################################################################################

##Hillsborough
stations = c("2301990",
             "2302010",
             "2303000",
             "2303330")
mat = plotRiver(stations)
USGSinput[,6] = mat[,which.max(colSums(mat,na.rm=T))]
print(stations[which.max(colSums(mat,na.rm=T))])
title("Hillsborough River")

###################################################################################
# 
# ##Little Manatee
# stations = c("2300100",
#              "2300210",
#              "2300300",
#              "2300500")
# mat = plotRiver(stations)
# USGSinput[,7] = mat[,which.max(colSums(mat,na.rm=T))]
# print(stations[which.max(colSums(mat,na.rm=T))])
# title("Little Manatee River")
# 
# ###################################################################################

##Manatee
stations = c("2299950",
             "2300100",
             "2300300",
             "2300500")
mat = plotRiver(stations)
USGSinput[,7] = mat[,which.max(colSums(mat,na.rm=T))]
print(stations[which.max(colSums(mat,na.rm=T))])
title("Manatee River")

###################################################################################
##Myakka
stations = c("2297155",
             "2298488",
             "2298492",
             "2298495",
             "2298527",
             "2298530",
             "2298554",
             "2298608",
             "2298830",
             "2299410",
             "2299950")

mat = plotRiver(stations)
USGSinput[,8] = mat[,which.max(colSums(mat,na.rm=T))]
print(stations[which.max(colSums(mat,na.rm=T))])
title("Myakka River")
###################################################################################

USGSinput = as.data.frame(USGSinput)
USGSinput = cbind(USGSinput, as.factor(states))
colnames(USGSinput) = c("Tampa",
                        "Peace",
                        "Okeechobee",
                        "Suwanee",
                        "Withlacoochee",
                        "Hillsborough",
#                        "Little Manatee",
                        "Manatee",
                        "Myakka",
                        "State")

write.csv(USGSinput,"C:/Users/xswang/HAB Research/USGS_discharge_total.csv")

###################################################################################

