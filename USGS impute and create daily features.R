###########################################################################################
## "USGS impute and create daily features.R"
## Convert daily data and organize it by features
###########################################################################################
rm(list = ls())
setwd("C:/Users/xswang/HAB Research/USGS daily data")

features = read.csv("C:\Users\xswang\HAB Research\USGS_features.csv")
files = list.files()
for(feature in features[,1]){
  discharge_files = c()
  minDate = as.Date("2018/01/01")
  maxDate = as.Date("1990/01/01")
  colsize = 0
  cnames  = c()
  for(i in 1:length(files)){
    r = read.csv(files[i])
    if(any(grepl(feature,colnames(r)))) {
      discharge_files = c(discharge_files, i)
      g = grep(feature,colnames(r))
      colsize = colsize + length(g)
      cnames  = c(cnames, paste0(strsplit(files[i], split = " ")[[1]][2]," ",colnames(r)[g]))
      minDate = min(c(minDate, as.Date(r$datetime)))
      maxDate = max(c(maxDate, as.Date(r$datetime)))
    }
  }
  
  ##Create dataframe
  rownam = c(minDate)
  for(i in 1:(maxDate-minDate)){
    rownam = c(rownam, minDate+i)
  }
  dfval  = data.frame(matrix(NA, nrow=length(rownam), ncol=colsize+1))
  colnames(dfval) = c("date",cnames)
  
  dfval[,1] = rownam
  
  colind = 1
  for(i in 1:length(discharge_files)){
    ##Read in data
    r = read.csv(files[discharge_files[i]])
    
    ##Get dates and rowind
    dates = as.Date(r$datetime)
    rowind   = match(dates, rownam)
    
    ##Get needed columns
    g     = grep(feature,colnames(r))
    
    ##Put in dataframe
    dfval[rowind, colind+(1:length(g))] = r[,g]
    ##Adjust column index
    colind = colind+length(g)
  }
  ##################################################################################
  predictFromLM = function(lmodel, pred){
    return(pred*lmodel$coefficients[2]+lmodel$coefficients[1])
  }
  ## Create matrix where its order by the value r-squared
  rsquaremat = matrix(NA, nrow=ncol(dfval)-1, ncol = ncol(dfval)-1)
  for(i in 2:ncol(dfval)){
    for(j in 2:ncol(dfval)){
      if(i != j && any(!is.na(dfval[,i]) && !is.na(dfval[,j]))){
        rsquaremat[i-1,j-1] = summary(lm(dfval[,i]~dfval[,j]))$r.squared
        rsquaremat[j-1,i-1] = summary(lm(dfval[,j]~dfval[,i]))$r.squared
      }
      else{
        rsquaremat[i-1,j-1] = 0
        rsquaremat[j-1,i-1] = 0
      }
    }
  }
  
  switch = function(x){
    if(length(which(x))!=0) return(min(which(x)))
    return(NA)
  }
  for(i in 2:ncol(dfval)){
    ##NA values for imputation
    w = which(is.na(dfval[,i]))
    
    ##Find rsquare values greater than 0.7 and order them from
    ##largest to smallest
    o = order(rsquaremat[i-1,], decreasing = TRUE)
    if(any(rsquaremat[i-1,o]>=0.7)){
      o = o[1:max(which(rsquaremat[i-1,o]>0.7))]
      
      ##Check of the values are not NA
      o_mat = !is.na(dfval[w,o+1])
      
      ##Get linear model
      lmpredict = lapply(o+1, 
                         FUN=function(x) predictFromLM(lm(dfval[,i]~dfval[,x]), dfval[w,x]))
      lmpredict = as.data.frame(lmpredict)
      
      firstTrue = 0
      if(is.null(dim(o_mat))){
        firstTrue = o_mat
        dfval[w,i] = lmpredict
      }else{
        firstTrue = apply(o_mat, MARGIN=1, FUN = switch)
        dfval[w[which(!is.na(firstTrue))],i] = unlist(lapply(which(!is.na(firstTrue)), FUN = function(x) lmpredict[x,firstTrue[x]]))
      }
    }
  }
  colNA2 = rep(0, ncol(dfval)-1)
  for(i in 2:ncol(dfval)){
    colNA2[i-1] = length(which(is.na(dfval[,i])))
  }
  write.csv(dfval, paste0("C:/Users/xswang/HAB Research/USGS daily features/",feature,".csv"))
  #print(colNA)
  #print(c(feature,colNA2))
  print(feature)
  print(length(which(colNA2<10)))
}

unusable = c()