require(e1071)
require(neuralnet)
require(kernlab)
require(doParallel)
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
x = is.na(HAB$Abundance_cells)
#states[x[9:10]] =1
#####################################################################################
##Read in entire dataset
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = as.data.frame(scale(alldata))

alldata$State = as.factor(states)
#alldata = alldata[sample(nrow(alldata)),]#k-fold cross-validation
prevdata=alldata
alldata$State = as.numeric(as.character(alldata$State))

alldata = alldata[!x,]
HAB = HAB[!x,]
years = format(HAB$Date, "%Y")
#alldata = alldata[is.na(match(years, c(1998,2002,2009,2010,2013))),]##Get rid of years with offshore forcing
#####################################################################################
##Function
predvalues = function(pred, test_data){
  event     = test_data==1
  eventpred = pred[event]
  perf1 = length(which(eventpred==1))/length(eventpred)
  
  no     = test_data== 0
  nopred = pred[no]
  perf2 = length(which(nopred==0))/length(nopred)
  
  perf3 = length(which(pred==test_data))/length(pred)
  return(c(perf1,perf2,perf3))
}
predrvm = function(HABmod, alldata, test){
  p = as.numeric(predict(HABmod,alldata[test,]))
  pred = rep(0, length(p))
  pred[p >= 0.5] = 1
  return(pred)
}
prednnet = function(HABnet, alldata, test){
  p = compute(HABnet, alldata[test,1:(ncol(alldata)-1)])$net.result
  pred = rep(0, length(p))
  pred[p >= 0.5] = 1
  return(pred)
}

#####################################################################################
##Cross validation
folds = 10
marks = round(seq(1, nrow(alldata),length.out=folds+1))
marks[folds+1] = nrow(alldata)

for(i in 1:folds){
  print(c(HAB$Date[marks[i]],
          HAB$Date[marks[i+1]-1]))
}

##Create cluster for parallel processing
cl <- makeCluster(detectCores()-2)#
registerDoParallel(cl)
## Used parallel processing to make it faster
testpred <- foreach(i=1:folds, .combine=rbind) %dopar%{
  require(DMwR)
  require(e1071)
  test = marks[i]:(marks[i+1]-1)
  train = -test
  ##SVM
  trainx = alldata[train,]
  trainx$State = as.factor(trainx$State)
  trainx = rbind(trainx[trainx$State==0,],
                 SMOTE(State~., trainx, perc.over = 100,k=5,perc.under=0))
  # trainx = rbind(trainx[trainx$State==0,],
  #                trainx[as.character(trainx$State)=="1",],
  #                trainx[as.character(trainx$State)=="1",])
  HABsvm = tune.svm(State~.,data=trainx,
                    type="C-classification",
                    cost = 2^(-5:10),
                    kernel="radial",
                    probability=T)
  HABsvm = HABsvm$best.model
  
  ##Naive Bayes
  trainx = alldata[train,]
  trainx$State = as.factor(trainx$State)
  
  HABnbayes = naiveBayes(State~.,data=trainx)
  
  ##RVM
  require(kernlab)
  trainx = alldata[train,]
  trainx$State = as.factor(trainx$State)
  trainx = rbind(trainx[trainx$State==0,],
                 SMOTE(State~., trainx, perc.over = 50,k=5,perc.under=0))
  trainx$State = as.numeric(trainx$State)
  HABrvm = rvm(State~., data= trainx,
               kernel = "rbfdot")
  
  ##ANN
  require(neuralnet)
  HABnet = neuralnet(State~
                       venf1.NS + venf1.EW + venf1.ATMP +
                       cdrf1.NS + cdrf1.EW + X42039.NS +
                       X42039.EW + X42003.NS + X42003.EW +
                       Tampa + Peace + Okeechobee +
                       Suwanee + Withlacoochee + Hillsborough +
                       Manatee + Myakka + Caloosahatchee +
                       Tampa_TN_mg+Myakka_TN_mg+
                       Peace_TN_mg+Withlacoochee_TN_mg+
                       Manatee_TN_mg+Hillsborough_TN_mg+
                       Tampa_TP_mg+Myakka_TP_mg+
                       Peace_TP_mg+Withlacoochee_TP_mg+
                       Manatee_TP_mg+Hillsborough_TP_mg+
                       Caloosahatchee_TN_mg+Caloosahatchee_TP_mg,
                     #SSH,
                     data=alldata[train,],
                     linear.output = F,
                     hidden = c(20, 10))
  returnval = matrix(NA, nrow=length(test),ncol=5)
  returnval[,1:5] = cbind(as.numeric(as.character(alldata$State[test])),
                          as.numeric(predict(HABsvm,alldata[test,]))-1, 
                          as.numeric(predict(HABnbayes,alldata[test,]))-1,
                          predrvm(HABrvm, alldata[test,]),
                          prednnet(HABnet, alldata, test))
  returnval
}
stopCluster(cl)
colnames(testpred) = c("Real","SVM","NB","RVM","NET")
testpred = as.data.frame(testpred)
write.csv(testpred, 'testpred_SMOTE_nutrients.csv')

####################################################################################################
####CV metrics

##set up colnames
perf=as.data.frame(matrix(NA, nrow=folds, ncol=12),
                   row.names=paste0("trial.",1:10))
colnames(perf) = c("SVM HAB Accuracy", "SVM non-HAB Accuracy", "SVM Total Accuracy",
                   "BAY HAB Accuracy", "BAY non-HAB Accuracy", "BAY Total Accuracy",
                   "RVM HAB Accuracy", "RVM non-HAB Accuracy", "RVM Total Accuracy",
                   "ANN HAB Accuracy", "ANN non-HAB Accuracy", "ANN Total Accuracy")
for(i in 1:folds){
  testind = marks[i]:(marks[i+1]-1)
  test_data=alldata[testind,ncol(alldata)]
  
  predsvm = testpred$SVM[testind]
  prednb = testpred$NB[testind]
  predRVM = testpred$RVM[testind]
  predNET = testpred$NET[testind]
  acc =  c(predvalues(predsvm,test_data),predvalues(prednb, test_data),
           predvalues(predRVM,test_data),predvalues(predNET, test_data))
  perf[i,]=acc
}
print(perf)
print(colMeans(perf,na.rm=T))

##Save data
perf = rbind(perf, colMeans(perf,na.rm=T))
print(perf)
rownames(perf)[11]="average"
##Save data
#write.csv(perf, "cv_without_satellite_not_keeping_data.csv")
####################################################################################################
####Time series

testpred = read.csv("testpred.csv")
#testpred$X = NULL
HABstate = testpred$Real
##Match up dates with predictions
date = rep(as.Date("1900-01-01"), nrow(testpred))
for(i in 1:nrow(testpred)){
  date[i] = HAB$Date[which.min(abs(alldata$Peace_TN_mg[i]-prevdata$Peace_TN_mg))] 
}
rightsvm = date[testpred[,2]==HABstate]
rightbay = date[testpred[,3]==HABstate]
rightrvm = date[testpred[,4]==HABstate]
rightnet = date[testpred[,5]==HABstate]
firstpart = rep(FALSE, nrow(HAB))
firstpart[1:522] = TRUE

#####################################################################################