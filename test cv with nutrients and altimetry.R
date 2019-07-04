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
#####################################################################################

##Read in entire dataset
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = as.data.frame(scale(alldata))

alldata$State = as.factor(states)
alldata = alldata[sample(nrow(alldata)),]
alldata$State = as.numeric(as.character(alldata$State))

#####################################################################################
##Function
predmodclass = function(HABmod, alldata, test){
  pred = predict(HABmod,alldata[test,])
  event     = alldata[test,ncol(alldata)]==1
  eventpred = pred[event]
  perf1 = length(which(eventpred==1))/length(eventpred)
  
  no     = alldata[test,ncol(alldata)]== 0
  nopred = pred[no]
  perf2 = length(which(nopred==0))/length(nopred)
  
  perf3 = length(which(pred==alldata[test,ncol(alldata)]))/length(pred)
  return(c(perf1,perf2,perf3))
}
predmodregress = function(HABmod, alldata, test){
  p = as.numeric(predict(HABmod,alldata[test,]))
  pred = rep(0, length(p))
  pred[p >= 0.5] = 1
  event     = alldata[test,ncol(alldata)]==1
  eventpred = pred[event]
  perf1 = length(which(eventpred==1))/length(eventpred)
  
  no     = alldata[test,ncol(alldata)]== 0
  nopred = pred[no]
  perf2 = length(which(nopred==0))/length(nopred)
  
  perf3 = length(which(pred==alldata[test,ncol(alldata)]))/length(pred)
  return(c(perf1,perf2,perf3))
}
prednnet = function(HABnet, alldata, test){
  p = compute(HABnet, alldata[test,1:(ncol(alldata)-1)])$net.result
  pred = rep(0, length(p))
  pred[p >= 0.5] = 1
  event     = alldata[test,ncol(alldata)]==1
  eventpred = pred[event]
  perf1 = length(which(eventpred==1))/length(eventpred)
  
  no     = alldata[test,ncol(alldata)]== 0
  nopred = pred[no]
  perf2 = length(which(nopred==0))/length(nopred)
  
  perf3 = length(which(pred==alldata[test,ncol(alldata)]))/length(pred)
  return(c(perf1,perf2,perf3))
}
#####################################################################################
##Cross validation
marks = round(seq(1, nrow(alldata),length.out=11))
marks[11] = 1084

##Create cluster for parallel processing
cl <- makeCluster(detectCores()-4)#
registerDoParallel(cl)
## Used parallel processing to make it faster
perf <- foreach(i=1:10, .combine=rbind) %dopar%{
  require(e1071)
  test = marks[i]:(marks[i+1]-1)
  train = -test
  
  ##SVM
  trainx = alldata[train,]
  trainx = rbind(trainx,
                 trainx[as.character(trainx$State)=="1",],
                 trainx[as.character(trainx$State)=="1",])
  HABsvm = tune.svm(as.factor(State)~.,data=trainx,
                    type="C-classification",
                    cost = 2^(-5:10),
                    kernel="radial",
                    probability=T)
  HABsvm = HABsvm$best.model

  ##Naive Bayes
  trainx = alldata[train,]
  HABnbayes = naiveBayes(as.factor(State)~.,data=trainx)

  ##RVM
  require(kernlab)
  trainx = alldata[train,]
  trainx = rbind(trainx,
                 trainx[as.character(trainx$State)=="1",])
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
                       Caloosahatchee_TN_mg+Caloosahatchee_TP_mg+
                       SSH,
                     data=alldata[train,],
                     linear.output = F,
                     hidden = c(20, 10))
  #Print results
  c(predmodclass(HABsvm, alldata, test), predmodclass(HABnbayes, alldata, test),
    predmodregress(HABrvm, alldata, test), prednnet(HABnet, alldata, test))
  # c(predmodclass(HABnbayes, alldata, test))
}

##set up colnames
colnames(perf) = c("SVM HAB Accuracy", "SVM non-HAB Accuracy", "SVM Total Accuracy",
                   "BAY HAB Accuracy", "BAY non-HAB Accuracy", "BAY Total Accuracy",
                   "RVM HAB Accuracy", "RVM non-HAB Accuracy", "RVM Total Accuracy",
                   "ANN HAB Accuracy", "ANN non-HAB Accuracy", "ANN Total Accuracy")
print(perf)
print(colMeans(perf))

##Save data
perf = rbind(perf, colMeans(perf))
rownames(perf) = c(rownames(perf)[1:10],"Average")
print(perf)

##Save data
write.csv(perf, "cv_performance.csv")
stopCluster(cl)
#####################################################################################
##Create model trained on entire dataset

##Save model trained on entire dataset standardized and centered at zero
trainx = alldata
# trainx = rbind(trainx,
#                trainx[as.character(trainx$State)=="1",],
#                trainx[as.character(trainx$State)=="1",])
HABsvm = tune.svm(as.factor(State)~.,
                  data=trainx,
                  type="C-classification",
                  cost = 2^(-5:10),
                  kernel="radial",
                  probability=T)
HABsvm = HABsvm$best.model
# as.factor(State)~venf1.NS+venf1.EW+venf1.ATMP+cdrf1.NS+
#   cdrf1.EW+X42039.NS+X42039.EW+X42003.NS+X42003.EW+Tampa+Peace+Hillsborough+
#   Myakka+
#   Peace_TN_mg+Withlacoochee_TN_mg+Manatee_TN_mg+Hillsborough_TN_mg+
#   Tampa_TP_mg+Withlacoochee_TP_mg+Manatee_TP_mg+Hillsborough_TP_mg+
#   SSH+Caloosahatchee+Caloosahatchee_TN_mg+Caloosahatchee_TP_mg
##Naive Bayes
trainx = alldata
HABnbayes = naiveBayes(as.factor(State)~.,data=trainx)

##RVM
require(kernlab)
trainx = alldata
trainx = rbind(trainx,
               trainx[as.character(trainx$State)=="1",])
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
                     Caloosahatchee_TN_mg+Caloosahatchee_TP_mg+
                     SSH,
                   data=alldata,
                   linear.output = F,
                   hidden = c(20, 10))

print(c(predmodclass(HABsvm, alldata, 1:nrow(alldata)), 
        predmodclass(HABnbayes, alldata, 1:nrow(alldata)),
        predmodregress(HABrvm, alldata, 1:nrow(alldata)), 
        prednnet(HABnet, alldata, 1:nrow(alldata))))


##Save data
save(HABsvm, file="radialSVM.RData")
save(HABnbayes, file="naivebayes.RData")
save(HABrvm, file="radialRVM.RData")
save(HABnet, file="nnet.RData")
