##############################################################################################
##"Read in wind data.R"
## Read in wind data from National Data Buoy Center (NDBC) in National Oceanic and Atmospheric
## Administration (NOAA), National Ocean Service (NOS), Coastal Ocean Monitoring
## and Prediction System (COMPS), Everglades National Park
##############################################################################################
##Set working data and load in station location
require(testit) ##for a nice function
setwd("C:/Users/xswang/HAB Research")
wind_stations = read.csv("NDBC_names.csv")$Station
wind_stations = tolower(wind_stations)
##############################################################################################
##Find prefix and suffix 
prefix = "https://www.ndbc.noaa.gov/view_text_file.php?filename="
middle = "h"
suffix = ".txt.gz&dir=data/historical/stdmet/"
##############################################################################################
##Get altitude and longitude
mat    = matrix(NA, nrow = length(wind_stations), ncol = 3)
mat[,1]= wind_stations


getLatLon     = function(str){
  f = strsplit(substr(str,start=6, stop=nchar(str)), split = "N")[[1]]
  lat = as.numeric(f[1])
  f = f[2]
  lon = as.numeric(strsplit(f, split = "W")[[1]][1])
  return(c(lat,lon))
}
lonlatFromURL = function(station){
  URL = paste0("https://www.ndbc.noaa.gov/station_page.php?station=",station)
  page = read.csv2(url(URL))
  page = as.character(page[,1])
  for(i in 1:10){
    if(any(!is.na(getLatLon(page[pos+i])))){
      return(getLatLon(page[pos+i]))
    }  
  }
}

for(i in 1:length(wind_stations)){
  mat[i,2:3] = lonlatFromURL(wind_stations[i])
}
mat = as.data.frame(mat)
colnames(mat) = c("Station","Latitude","Longitude")
write.csv(mat, "NDBC_lat_lon.csv")
##############################################################################################
##Find values for wind stations

for(j in 1:length(wind_stations)){
  w = wind_stations[j]
  wind_df = NULL
  first   = TRUE
  for(y in 1998:2018){
    
    ##Get URL
    URL = paste0(prefix, w, middle, y, suffix)
    ##Check if such file exists
    if(!has_error(read.csv2(url(URL), sep = ""), silent = F)){
      ##Read from page
      page = read.csv2(url(URL), sep = "")
      
      ##Check if not wrong
      if(length(page)!=1){
        ##Make sure its the same thing
        colnames(page)[1] = "YYYY"
        ##Save as character
        for(i in 1:ncol(page)) page[,i] = as.character(page[,i])
        
        ##If first row is simply the unit, remove it
        if(all(is.na(as.numeric(page[1,])))){
          page = page[2:nrow(page),]
        }
        
        ##Put it in wind_df
        if(first){
          wind_df = page
          first = FALSE
        }else{
          ##Crete a new dataframe with space for new metrics
          temp = as.data.frame(matrix(NA, nrow = nrow(wind_df)+nrow(page), 
                                      ncol = length(unique(c(colnames(wind_df),colnames(page))))))
          ##Create colnames
          colnames(temp) = unique(c(colnames(wind_df),colnames(page)))
          
          ##Put in old wind station dataframe
          temp[1:nrow(wind_df),match(colnames(wind_df), colnames(temp))] = wind_df
          
          ##Put in new data
          temp[(nrow(wind_df)+1):nrow(temp),match(colnames(page), colnames(temp))] = page
          ##Put back in wind_df
          wind_df = temp
        }
      }
    }
  }
  if(!is.null(wind_df)){
    ##Write into csv 
    write.csv(wind_df, paste0("C:/Users/xswang/HAB Research/Wind buoys/",w, " wind station.csv"))
  }
}