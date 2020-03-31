require(e1071)
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Get Data
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
#####################################################################################

##Read in entire dataset
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = as.data.frame(scale(alldata))

alldata$State = as.factor(states)
alldata$State = as.numeric(as.character(alldata$State))

prevdata = alldata
alldata$SSH=NULL

pred = read.csv("testpred_block_no_SMOTE.csv")
pred$X = NULL
#######################################################################################
pred$Date = as.Date(pred$Date)
years = as.numeric(substr(pred$Date,0,4))
#######################################################################################
##Overall RVM
precision_ = function(x,y){
  length(which(y[x==1]==1))/length(which(x==1))
}
precision_(pred$RVM, pred$Real)

##1998, 2002, 2009, 2010, and 2013 
ty = c(1998, 2002, 2009, 2010, 2013 )
for( i in ty){
  print(c(i, precision_(pred$RVM[years==i], pred$Real[years==i])))
}
#######################################################################################3
recall_ = function(x, y){
  length(which(x[y==1]==1))/length(which(y==1))
}

recall_(pred$RVM, pred$Real)

for( i in ty){
  print(c(i, recall_(pred$RVM[years==i], pred$Real[years==i])))
}


