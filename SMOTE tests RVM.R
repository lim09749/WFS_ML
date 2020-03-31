cl <- makeCluster(detectCores()-2)#
registerDoParallel(cl)
require(beepr)
for(j in seq(500,700,25)){
  y=Sys.time()
  j=700
  ## Used parallel processing to make it faster
  testpred <- foreach(i=1:folds, .combine=rbind) %dopar%{
    require(DMwR)
    require(e1071)
    test = marks[i]:(marks[i+1]-1)
    train = -test
    
    ##RVM
    require(kernlab)
    trainx = alldata[train,]
    trainx$State = as.factor(trainx$State)
    trainx = rbind(trainx[trainx$State==0,],
                   SMOTE(State~., trainx, perc.over = j, k=5, perc.under=0))
    trainx$State = as.numeric(trainx$State)-1
    HABrvm = rvm(State~., data= trainx,
                 kernel = "rbfdot")
    
    cbind(as.numeric(as.character(alldata$State[test])),
          predrvm(HABrvm, alldata[test,]))
  }
  colnames(testpred) = c("Real","RVM")
  testpred = as.data.frame(testpred)
  perf=as.data.frame(matrix(NA, nrow=folds, ncol=3),
                     row.names=paste0("trial.",1:10))
  colnames(perf) = c("RVM HAB Accuracy", "RVM non-HAB Accuracy", "RVM Total Accuracy")
  
  for(i in 1:folds){
    testind = marks[i]:(marks[i+1]-1)
    test_data=alldata[testind,ncol(alldata)]
    
    predRVM = testpred$RVM[testind]
    # acc =  c(predvalues(predsvm,test_data),predvalues(prednb, test_data),
    #          predvalues(predRVM,test_data),predvalues(predNET, test_data))
    perf[i,]=predvalues(predRVM,test_data)
  }
  print(j)
  print(colMeans(perf,na.rm=T))
  print(Sys.time()-y)
  beep()
}
beep()