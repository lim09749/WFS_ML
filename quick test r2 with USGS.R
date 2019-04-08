#######################################################
##"quick test r2 with USGS.R"
##Compare weekly averages
#######################################################
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/NDBC weekly data")
#######################################################
##Read in Karenia Brevis cell counts
HAB = read.csv("C:/Users/xswang/HAB Research/Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))

minTrue = 500
totaldata = matrix(NA, nrow=nrow(HAB), ncol=1)
totaldata[,1] = HAB$Abundance_cells
files = list.files()
for(i in 1:length(files)){
  f = read.csv(files[i])
  f$X = NULL
  print(as.character(f[1,1]))
  for(j in 2:ncol(f)){
    column = f[1:nrow(totaldata), j]
    if(length(which(!is.na(column)))>=minTrue){
      totaldata = cbind(totaldata, column)
    }
  }
}
print(ncol(totaldata))
setwd("C:/Users/xswang/HAB Research/USGS weekly features")
files = list.files()
for(i in 1:length(files)){
  f = read.csv(files[i])
  f$X = NULL
  print(as.character(f[1,1]))
  for(j in 2:ncol(f)){
    column = f[1:nrow(totaldata), j]
    if(length(which(!is.na(column)))>=minTrue){
      totaldata = cbind(totaldata, column)
    }
  }
}
print(ncol(totaldata))

##NDBC
# [1] "day"  "WSPD" "WDIR" "GST"  "WVHT" "DPD"  "APD"  "MWD"  "ATMP" "WTMP" "DEWP" "mm"  
# [13] "PRES" "TIDE"
# uniqcol = c("WSPD", "WDIR", "GST" , "WVHT" ,"DPD",  
#             "APD",  "MWD",  "ATMP", "WTMP" ,"DEWP", "mm"  ,
#             "PRES", "TIDE")
# for(name in uniqcol){
#   par(mfrow=c(3,1))
#   for(i in 1:length(files)){
#     f = read.csv(files[i])
#     f$X = NULL
#     if(name%in%colnames(f) && nrow(f)>1000){
#       #plot(as.Date(as.character(f$day)), f[,colnames(f)==name],type="l")
#       #ind = which(HAB$Abundance_cells>10^6)
#       #points(as.Date(as.character(f$day))[ind],f[ind,colnames(f)==name],col="red")
#       totaldata = cbind(totaldata, c(f[,colnames(f)==name],rep(NA,nrow(totaldata)-nrow(f))))
#     }
#   }
# }

returncol = function(x){
  x[is.na(x)] = mean(x, na.rm=T)
  return(x)
}
totaldata = apply(totaldata, 2, FUN = returncol)
totaldata = scale(totaldata)
totaldata = totaldata[!is.na(totaldata[,1]),]
totaldata = as.data.frame(totaldata)
colnames(totaldata)[1] = "Abundance"

randomizerows = sample(1:nrow(totaldata),replace=F)
states = rep(1, nrow(HAB))
states[HAB$Abundance_cells>=10^6] = 2

dftotal = totaldata[randomizerows,]
states    = states[randomizerows]
colsd = rep(NA, ncol(dftotal))
for(i in 1:ncol(dftotal)){
  colsd[i] = sd(dftotal[,i])
}
dftotal = dftotal[,!is.na(colsd)]
classification_results = function(pred, real){
  pred = as.character(pred)
  real = as.character(real)
  print("Correct")
  print(length(which(pred==real))/length(pred))
  print("False positive")
  print(length(which(pred=="1" & real=="2"))/length(pred))
  print("False negative")
  print(length(which(pred=="2" & real=="1"))/length(pred))
  print("Accuracy of positive")
  print(length(which(pred[real=="2"]==real[real=="2"]))/length(real[real=="2"]))
  print("Accuracy of negative")
  print(length(which(pred[real=="1"]==real[real=="1"]))/length(real[real=="1"]))
  
}

test = 1:(nrow(dftotal)/10)
#################################################################################################
##Naive bayes
require(e1071)
HABnaivebayes = naiveBayes(as.factor(states[-test])~ .,
                           data = dftotal[-test,-1])
prednbayes    = predict(HABnaivebayes, dftotal[test,-1])
print(classification_results(prednbayes, as.factor(states[test])))
# HABnaivebayes = naiveBayes(as.factor(states[-test])~ .,
#                            data = dftotal[-test,-1])
# prednbayes    = predict(HABnaivebayes, dftotal[test,-1])
# print(classification_results(prednbayes, as.factor(states[test])))

#################################################################################################
##tree
library(tree)
tree.carseats = tree(as.factor(states[-test])~ .,
                     data = dftotal[-test,-1])
summary ( tree.carseats )
plot( tree.carseats )
text( tree.carseats , pretty =0)
predtree = predict (tree.carseats , dftotal[test,-1], type ="class")
print(classification_results(predtree,as.factor(states[test])))
#################################################################################################
##knn
require(knn)
HABknn <- train(as.factor(states[-test])~ .,
                data = dftotal[-test,-1], method='knn')
predknn    = predict(HABknn, dftotal[test,-1])
print(classification_results(predknn, as.factor(states[test])))
#################################################################################################
##SVM
for(hyperC in 2^(2:20)){
  HABksvm <- ksvm(as.factor(states[-test])~ .,
                  data = dftotal[-test,-1],
                  kernel="rbfdot",
                  type="C-svc",
                  C=hyperC)
  predksvm = predict(HABksvm, dftotal[test,-1])
  print(hyperC)
  print(classification_results(predksvm, as.factor(states[test])))
  print("#######################################################")
}
# for(hyperC in 2^(1:20)){
#   HABksvm <- ksvm(as.factor(states[-test])~ .,
#                   data = dftotal[-test,-1],
#                   kernel="rbfdot",
#                   type="C-svc",
#                   C=hyperC)
#   predksvm = predict(HABksvm, dftotal[test,-1])
#   print(hyperC)
#   print(classification_results(predksvm, as.factor(states[test])))
#   print("#######################################################")
# }
#################################################################################################
##GLM

HABglm = glm(as.factor(states[-test]-1) ~ .,
            data=dftotal[-test,-1], family='binomial')
predglm = predict(HABglm, dftotal[test,-1])
