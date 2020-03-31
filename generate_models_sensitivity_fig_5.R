rm(list=ls())
require(e1071)
require(kernlab)
require(neuralnet)
require(beepr)
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Create functions
fixbounds = function(x){
  y=x
  y[x<=0]=0
  y[x>=1]=1
  y
}
ann_prob = function(HABnet, data){
  fixbounds(predict(HABnet, data))[,1]
}
bay_prob = function(HABnbayes, data){
  fixbounds(predict(HABnbayes, data, "raw")[,2])
}
rvm_prob = function(HABrvm, data){
  fixbounds(predict(HABrvm, data))[,1]
}
svm_prob = function(HABsvm, data){
  fixbounds(attr(predict(HABsvm, data, probability=T),"probabilities")[,2]) 
}
undoscale = function(ind, val, scalefactors){
  return(val * attr(scalefactors,"scaled:scale")[ind]
         + attr(scalefactors,"scaled:center")[ind])
}
##############################################################################################################

##Read in data

HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
x = is.na(HAB$Abundance_cells)


alldata = read.csv("alldata.csv")
alldata$X = NULL

scalefactors = scale(alldata)
alldata = as.data.frame(scale(alldata))

alldata$State = as.factor(states)
alldata$State = as.numeric(as.character(alldata$State))

##Create previous values dataset
prevdata = read.csv("alldata.csv")
prevdata$X = NULL
prevdata$State = as.factor(states)

##Remove NA values
alldata = alldata[!x,]
HAB = HAB[!x,]
prevdata = prevdata[!x, ]
#####################################################################################
##Discharge

##SVM
svmdata = alldata[,c(1:9,10:17,31)]
svmdata$State = as.factor(alldata$State)
HABsvm = tune.svm(as.factor(State)~.,
                  data=svmdata,
                  type="C-classification",
                  cost = 2^(-5:10),
                  kernel="radial",
                  probability=T)
svm_dis = HABsvm$best.model
save(svm_dis,file="new_svm_dis.RData")


##Naive Bayes
nbdata = alldata[,c(1:9,10:17,31)]
nbdata$State = as.factor(alldata$State)
bay_dis = naiveBayes(State~.,data=nbdata)

save(bay_dis, file = "new_bay_dis.RData")

##RVM
require(kernlab)
rvmdata = alldata[,c(1:9,10:17,31)]
rvmdata$State = as.numeric(alldata$State)
rvm_dis = kernlab::rvm(State~., data= rvmdata,
             kernel = "rbfdot")
save(rvm_dis, file = "new_rvm_dis.RData")

##ANN
require(neuralnet)
anndata = alldata[,c(1:9,10:17,31)]
anndata$State = as.numeric(alldata$State)
ann_dis = neuralnet(State~.,
                   data=anndata[1:(nrow(anndata)*0.9),],
                   linear.output = F,
                   hidden = c(20, 10))
save(ann_dis, file = "new_ann_dis.RData")
################################################################################################33
##TP
# states= states[alldata$Caloosahatchee_TP_mg>=4]
# HAB = HAB[alldata$Caloosahatchee_TP_mg>=4,]

##SVM
svmdata = alldata[,c(1:9,24:29,33)]
svmdata$State = as.factor(alldata$State)
HABsvm = tune.svm(as.factor(State)~.,
                  data=svmdata,
                  type="C-classification",
                  cost = 2^(-5:10),
                  kernel="radial",
                  probability=T)
svm_tp = HABsvm$best.model
save(svm_tp,file="new_svm_tp.RData")


##Naive Bayes
nbdata = alldata[,c(1:9,24:29,33)]
nbdata$State = as.factor(alldata$State)
bay_tp = naiveBayes(State~.,data=nbdata)
save(bay_tp, file = "new_bay_tp.RData")

##RVM
require(kernlab)
rvmdata = alldata[,c(1:9,24:29,33)]
rvmdata$State = as.numeric(alldata$State)
rvm_tp = kernlab::rvm(State~., data= rvmdata,
             kernel = "rbfdot")
save(rvm_tp, file = "new_rvm_tp.RData")

##ANN
require(neuralnet)
anndata = alldata[,c(1:9,10:17,31)]
anndata$State = as.numeric(alldata$State)
ann_tp = neuralnet(State~.,
                   data=anndata,
                   linear.output = F,
                   hidden = c(20, 10))
save(ann_tp, file = "new_ann_tp.RData")
####################################################################################################
##TN
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = alldata[,c(1:9,18:23,32)]
states=states[alldata$Tampa_TN_mg>6 & alldata$Caloosahatchee_TN_mg>8 & !x]
alldata = alldata[alldata$Tampa_TN_mg>6 & alldata$Caloosahatchee_TN_mg>8 & !x,]
prevdata = alldata
scalefactors=scale(alldata)
alldata=as.data.frame(scale(alldata))
alldata$State=as.factor(states)

##SVM
svmdata = alldata
svmdata$State = as.factor(alldata$State)
HABsvm = tune.svm(as.factor(State)~.,
                  data=svmdata,
                  type="C-classification",
                  cost = 2^(-5:10),
                  kernel="radial",
                  probability=T)
svm_tn = HABsvm$best.model
save(svm_tn,file="new_svm_tn.RData")


##Naive Bayes
nbdata = alldata
nbdata$State = as.factor(alldata$State)
bay_tn = naiveBayes(State~.,data=nbdata)
save(bay_tn, file = "new_bay_tn.RData")

##RVM
require(kernlab)
rvmdata = alldata
rvmdata$State = as.numeric(alldata$State)-1
rvm_tn = kernlab::rvm(State~., data= rvmdata,
             kernel = "rbfdot")
save(rvm_tn, file = "new_rvm_tn.RData")

##ANN
require(neuralnet)
anndata = alldata
anndata$State = as.numeric(alldata$State)
ann_tn = neuralnet(State~.,
                   data=anndata,
                   linear.output = F,
                   hidden = c(20, 10))
save(ann_tn, file = "new_ann_tn.RData")
