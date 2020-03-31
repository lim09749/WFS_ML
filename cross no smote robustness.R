rm(list=ls())

require(e1071)
require(neuralnet)
require(kernlab)
require(doParallel)
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
#####################################################################################
##Read in entire dataset
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = as.data.frame(scale(alldata))

alldata$State = as.factor(states)
alldata$State = as.numeric(as.character(alldata$State))

alldata = alldata[!x,]
HAB = HAB[!x,]
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
fixbounds = function(x){
  y=x
  y[x<=0]=0
  y[x>=1]=1
  y
}
ann_prob = function(HABnet, data){
  fixbounds(predict(HABnet, data))
}
bay_prob = function(HABnbayes, data){
  fixbounds(predict(HABnbayes, data, "raw")[,2])
}
rvm_prob = function(HABrvm, data){
  fixbounds(predict(HABrvm, data))
}
svm_prob = function(HABsvm, data){
  fixbounds(attr(predict(HABsvm, data, probability=T),"probabilities")[,2]) 
}

#####################################################################################
##k-fold cross validation
folds = 10
marks = round(seq(1, nrow(alldata),length.out=folds+1))
marks[folds+1] = nrow(alldata)

for(i in 1:folds){
  print(c(HAB$Date[marks[i]],
          HAB$Date[marks[i+1]-1]))
}

set.seed(42)
randomizer = sample(nrow(alldata))
alldata = alldata[randomizer,]
HAB = HAB[randomizer,]

##Create cluster for parallel processing
cl <- makeCluster(detectCores()-1)#
registerDoParallel(cl)
s=Sys.time()
## Used parallel processing to make it faster
testpred <- foreach(i=1:folds, .combine=rbind) %dopar%{
  require(DMwR)
  require(e1071)
  test = marks[i]:(marks[i+1]-1)
  train = -test
  ##SVM
  trainx = alldata[train,]
  trainx$State = as.factor(trainx$State)
  trainx = rbind(trainx,
                 trainx[as.character(trainx$State)=="1",],
                 trainx[as.character(trainx$State)=="1",])
  HABsvm = tune.svm(State~.,data=trainx,
                    type="C-classification",
                    cost = 2^(-5:10),
                    kernel="radial",
                    probability=T)
  HABsvm = HABsvm$best.model
  
  ##Naive Bayes
  HABnbayes = naiveBayes(as.factor(State)~.,
                         data=trainx)
  
  # #RVM
  require(kernlab)
  trainx$State = as.numeric(trainx$State)-1
  HABrvm = rvm(State~., data= trainx,
               kernel = "rbfdot")
  
  
  # #ANN
  require(neuralnet)
  HABnet = neuralnet(State~.,
                     data=trainx,
                     linear.output = F,
                     likelihood=T,
                     hidden = c(20, 10))
  
  cbind(HABsvm$tot.nSV,
        NA,
        HABrvm@nRV,
        HABnet$result.matrix[which(rownames(HABnet$result.matrix)=="aic")],
        HABnet$result.matrix[which(rownames(HABnet$result.matrix)=="bic")])
}
print(Sys.time()-s)
stopCluster(cl)
colnames(testpred) = c("SVM","NB","RVM","NET_AIC","NET_BIC")
testpred = as.data.frame(testpred)
write.csv(testpred, 'robustness cross no smote.csv')
