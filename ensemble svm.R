
for(i in 1:10){
  testind = sample(nrow(dftotal),nrow(dftotal)/10)
  maxperf = -5
  maxSVM = NULL
  for(hyperC in 1:20){
    HABksvm <- ksvm(as.factor(states[-testind])~ .,
                    data = dftotal[-testind,-1],
                    kernel="rbfdot",
                    type="C-svc",
                    C=hyperC)
    predksvm = predict(HABksvm, dftotal[testind,-1])
    pred = as.character(predksvm)
    real = as.character(states[testind])
    result = length(which(pred==real))/length(real)
    if(result > maxperf){
      maxperf = result
      maxSVM = HABksvm
    }
  }
  predksvm = predict(maxSVM, dftotal[testind,-1])
  print(classification_results(predksvm, as.factor(states[testind])))
  print("#######################################################")
}

ensemblepredict = apply(perf, 2, FUN=function(x) length(which(x=="2")))
print(classification_results(ensemblepredict, as.factor(states[test])))




require(doParallel)
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

test = 1:(nrow(dftotal)/10)
## Used parallel processing to make it faster
perf <- foreach(i=(max(test)+1):nrow(dftotal), .combine=rbind) %dopar%{
  require(kernlab) ##Needs to load libraries in R
  train = (max(test)+1):nrow(dftotal)
  train = train[ train != i]
  ##Train svm
  maxperf = 0
  maxSVM = NULL
  for(hyperC in 1:20){
    HABksvm <- ksvm(as.factor(states[-test])~ .,
                    data = dftotal[-test,-1],
                    kernel="rbfdot",
                    type="C-svc",
                    C=hyperC)
    predksvm = predict(HABksvm, dftotal[test,-1])
    pred = as.character(predksvm)
    real = as.character(states[test])
    result = length(which(pred[real=="2"]==real[real=="2"]))/length(real[real=="2"])
    if(result > maxperf){
      maxperf = result
      maxSVM = HABksvm
    }
  }

  ##Adjoin the information given
  predksvm = predict(maxSVM, dftotal[test,-1])
  c(as.character(predksvm)) ##predicted on test dataset
}

