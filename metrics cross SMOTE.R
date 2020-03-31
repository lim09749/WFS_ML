rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Read in data
testpred = read.csv('testpred_cross_SMOTE.csv')
testpred$X = NULL

####################################################################################
##Get accuracies
tot_acc = function(x, y){
  length(which(x==y))/length(x)
}
hab_acc = function(x, y){
  length(which(x[y==1]==1))/length(which(y==1))
}
non_hab_acc = function(x,y){
  length(which(x[y==0]==0))/length(which(y==0))
}

recall_ = function(x, y){
  TP = length(which(x==1 & y==1))
  FP = length(which(x==1 & y==0))
  TN = length(which(x==0 & y==0))
  FN = length(which(x==0 & y==1))
  
  return(TP/(TP+FN))
}

precision_ = function(x,y){
  TP = length(which(x==1 & y==1))
  FP = length(which(x==1 & y==0))
  TN = length(which(x==0 & y==0))
  FN = length(which(x==0 & y==1))
  
  return(TP/(TP+FP))
}

f1_ = function(x,y,beta){
  (1+beta^2)*(recall_(x,y)*precision_(x,y))/(recall_(x,y)+beta^2*precision_(x,y))
}



matthews = function(x,y){
  TP = length(which(x==1 & y==1))
  FP = length(which(x==1 & y==0))
  TN = length(which(x==0 & y==0))
  FN = length(which(x==0 & y==1))
  (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN))*sqrt((TN+FP)*(TN+FN)))
}
######################################################################################3
##Get accuracy matrix
acc_df = as.data.frame(matrix(NA, nrow = 4, ncol=7))
rownames(acc_df) = c("SVM","RVM","NB","NET")
colnames(acc_df) = c("Total Accuracy",
                     "HAB Accuracy",
                     "Non-HAB Accuracy",
                     "Recall",
                     "Precision",
                     "F1",
                     "Matthews")
##SVMs
beta=1
acc_df[1,] = c(tot_acc(testpred$SVM, testpred$Real),
               hab_acc(testpred$SVM, testpred$Real),
               non_hab_acc(testpred$SVM, testpred$Real),
               recall_(testpred$SVM, testpred$Real),
               precision_(testpred$SVM, testpred$Real),
               f1_(testpred$SVM, testpred$Real, beta),
               matthews(testpred$SVM, testpred$Real))

##RVM
acc_df[2,] = c(tot_acc(testpred$RVM, testpred$Real),
               hab_acc(testpred$RVM, testpred$Real),
               non_hab_acc(testpred$RVM, testpred$Real),
               recall_(testpred$RVM, testpred$Real),
               precision_(testpred$RVM, testpred$Real),
               f1_(testpred$RVM, testpred$Real, beta),
               matthews(testpred$RVM, testpred$Real))
##NB
acc_df[3,] = c(tot_acc(testpred$NB, testpred$Real),
               hab_acc(testpred$NB, testpred$Real),
               non_hab_acc(testpred$NB, testpred$Real),
               recall_(testpred$NB, testpred$Real),
               precision_(testpred$NB, testpred$Real),
               f1_(testpred$NB, testpred$Real, beta),
               matthews(testpred$NB, testpred$Real))

##NET
acc_df[4,] = c(tot_acc(testpred$NET, testpred$Real),
               hab_acc(testpred$NET, testpred$Real),
               non_hab_acc(testpred$NET, testpred$Real),
               recall_(testpred$NET, testpred$Real),
               precision_(testpred$NET, testpred$Real),
               f1_(testpred$NET, testpred$Real, beta),
               matthews(testpred$NET, testpred$Real))
##################################################################################################3
write.csv(acc_df, "cross smote perf.csv")
