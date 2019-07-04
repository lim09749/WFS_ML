###########################################################################################
## "USGS combine same station files of different years.R"
## Take csv files of stations and combine them and fill in heading
###########################################################################################
rm(list = ls())
setwd("C:/Users/xswang/HAB Research/USGS raw data")
###########################################################################################
##Read in data

stations = read.csv("C:/Users/xswang/HAB Research/all florida stations.csv")
stations[,1] = as.character(stations[,1])
stations[,2] = as.character(stations[,2])
stations[,3] = as.character(stations[,3])
###########################################################################################
##Write as csv

##Iterate through stations
for(z in 1:length(stations$site_no)){
  ##Get files for heading and data
  station_num  = stations$site_no[z]
  headingprefix = paste0(station_num, " Heading")
  dataprefix    = paste0(station_num, " from")
  lfiles        = list.files()
  
  headingfiles  = lfiles[grepl(headingprefix, lfiles)]
  datafiles     = lfiles[grepl(dataprefix, lfiles)]
  
  ##Get headings for all the datafields
  totalheadings = c()
  totalrow      = 0
  
  if(length(headingfiles)!=0){
    for(i in 1:length(headingfiles)){
      ##Read in heading data
      r   = read.csv(headingfiles[i])
      r$X = NULL
      r   = r$x
      r   = as.character(r)
      
      ##Read in station data
      df   = read.csv(datafiles[i])
      df$X = NULL
      
      ##Get index and all headings
      ind = match("#    TS_ID       Parameter Description",r)
      nheadings = c()
      val = 3
      tmp = ""
      while(r[ind] != "#" && val <= ncol(df)){
        if(val+2 <= ncol(df) && grepl(strsplit(colnames(df)[val+2],split="_")[[1]][2], r[ind])){
          nheadings = c(nheadings, tmp,"next")
          val  = val +2
          code = strsplit(colnames(df)[val],split="_")[[1]][2]
          tmp  = trimws(strsplit(r[ind],
                                 split=code)[[1]][2])
        }
        else if(r[ind] != "#    TS_ID       Parameter Description"){
          tmp = paste0(tmp, r[ind])
        }
        ind = ind + 1
      }
      ##Combine with colnames
      nheadings = c(nheadings, tmp,"next")
      nheadings = nheadings[3:length(nheadings)]
      colnames(df)[5:ncol(df)] = nheadings
      ##Keep important data
      df = df[,c(1,2,3,4,seq(from=5,to=ncol(df),by=2))]
      ##Update total row
      totalrow = totalrow + nrow(df)
      ##Add headings
      totalheadings = unique(c(totalheadings, colnames(df)))
    }
    
    ##Total dataset
    totaldf = data.frame(matrix(NA, nrow = totalrow, ncol = length(totalheadings)))
    colnames(totaldf) = totalheadings
    currrow = 1
    for(i in 1:length(headingfiles)){
      ##Header files
      r   = read.csv(headingfiles[i])
      r$X = NULL
      r   = r$x
      r   = as.character(r)
      
      ##Dataframe
      df   = read.csv(datafiles[i])
      df$X = NULL
      
      ##Get index, and total heading
      ind = match("#    TS_ID       Parameter Description",r)
      nheadings = c()
      val = 3
      tmp = ""
      while(r[ind] != "#" && val <= ncol(df)){
        if(val+2 <= ncol(df) && grepl(strsplit(colnames(df)[val+2],split="_")[[1]][2], r[ind])){
          nheadings = c(nheadings, tmp,"next")
          val  = val +2
          code = strsplit(colnames(df)[val],split="_")[[1]][2]
          tmp  = trimws(strsplit(r[ind],
                                 split=code)[[1]][2])
        }
        else if(r[ind] != "#    TS_ID       Parameter Description"){
          tmp = paste0(tmp, r[ind])
        }
        ind = ind + 1
      }
      nheadings = c(nheadings, tmp,"next")
      nheadings = nheadings[3:length(nheadings)]
      colnames(df)[5:ncol(df)] = nheadings
      
      ##Save dataframe important valuse, save first four columns as strings, and place it
      ##in the correct position
      df = df[,c(1,2,3,4,seq(from=5,to=ncol(df),by=2))]
      for(j in 1:4) df[,j] = as.character(df[,j])
      totaldf[currrow:(currrow+nrow(df)-1),match(colnames(df),totalheadings)] = df
      currrow = currrow + nrow(df)
    }
    write.csv(totaldf, paste0("C:/Users/xswang/HAB Research/USGS raw total/",
                              "USGS ", station_num, " total.csv"))
  }
}