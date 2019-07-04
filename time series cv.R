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
alldata = read.csv("totaldata.csv")
alldata$State = NULL
alldata$X = NULL

##Read in SSH data
tempssh = read.csv("SSH.csv")
tempssh$X = NULL
alldata = cbind(alldata, as.numeric(as.character(tempssh$SSH)))
colnames(alldata)[ncol(alldata)] = "SSH"

##Impute and process data
for(i in 1:ncol(alldata)) alldata[is.na(alldata[,i]),i]=mean(alldata[,i],na.rm=T)
alldata[,10:23] = log10(alldata[,10:23]+4)
alldata[,1:(ncol(alldata)-1)] = scale(alldata[,1:(ncol(alldata)-1)])

alldata$State = as.factor(states)
prevdata = alldata
alldata = alldata[sample(nrow(alldata)),]

#####################################################################################
##Get HAB values (use fact that all of peace river data is unique)
for(i in 1:ncol(alldata)) print(c(colnames(alldata)[i],length(unique(alldata[,i]))))
HABval = rep(NA, nrow(alldata))
for(i in 1:length(HABval)){
  HABval[i] = HAB$Abundance_cells[which.min(abs(alldata$Peace_TN_mg[i]-prevdata$Peace_TN_mg))] 
}
HABstate = rep(0, length(HABval))
HABstate[HABval >= 10^5] = 1
#####################################################################################
##Helper functions
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
##Create time series
marks = round(seq(1, nrow(alldata),length.out=11))
marks[11] = 1084

## Used parallel processing to make it faster
cl <- makeCluster(detectCores()-4)#
registerDoParallel(cl)
testpred <- foreach(i=1:10, .combine=rbind) %dopar%{
  require(e1071)
  test = marks[i]:(marks[i+1]-1)
  train = -test
  
  ##SVM
  trainx = alldata[train,]
  trainx = rbind(trainx,
                 trainx[as.character(trainx$State)=="1",])
  HABsvm = tune.svm(as.factor(State)~.,data=trainx,
                    type="C-classification",
                    cost = 2^(-5:10),
                    kernel="radial",
                    probability=T)
  HABsvm = HABsvm$best.model
  
  ##Naive Bayes
  trainx = alldata[train,]
  trainx = rbind(trainx,
                 trainx[as.character(trainx$State)=="1",])
  HABnbayes = naiveBayes(as.factor(State)~.,data=trainx)

  ##RVM
  require(kernlab)
  trainx = alldata[train,]
  trainx = rbind(trainx,
                 trainx[as.character(trainx$State)=="1",])
  HABrvm = rvm(as.numeric(as.character(State))~., data= trainx,
               kernel = "rbfdot")

  ##ANN
  require(neuralnet)
  HABnet = neuralnet(as.numeric(as.character(State))~venf1.NS + venf1.EW + venf1.ATMP + 
                       cdrf1.NS + cdrf1.EW + X42039.NS + 
                       X42039.EW + X42003.NS + X42003.EW + 
                       Okeechobee + Suwanee + 
                       Tampa_TN_mg+Myakka_TN_mg+
                       Peace_TN_mg+Withlacoochee_TN_mg+
                       Manatee_TN_mg+Hillsborough_TN_mg+
                       Tampa_TP_mg+Myakka_TP_mg+
                       Peace_TP_mg+Withlacoochee_TP_mg+
                       Manatee_TP_mg+Hillsborough_TP_mg+SSH,
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
#####################################################################################
##Read in data
colnames(testpred) = c("Real","SVM","NB","RVM","NET")
write.csv(testpred, 'testpred.csv')
testpred = read.csv("testpred.csv")
testpred$X = NULL
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
##Plot green as correct and red as incorrect
tiff(file = "SVM CV time series.tiff", width = 18, height = 3.75, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(1,2))
plot(HAB$Date[firstpart], 
     log10(HAB$Abundance_cells)[firstpart],type="l",ylim=c(2,8), col = "black",
     xlab = "",
     ylab="Kb (log10(c/l))",
     main = expression(paste(italic("K. brevis")," abundance in West Florida Shelf (SVM)")),
     axes=F,
     cex.main=1.5,cex.lab=1.5)
axis.Date(1, at=seq(as.Date("1998-01-01"), as.Date("2008-01-01"), by="year"),
          cex.axis=1.25)
axis(2,cex.axis=1.25)
abline(a=5,b=0, col = "blue")
points(HAB$Date[firstpart & !is.na(match(HAB$Date, rightsvm))],
       log10(HAB$Abundance_cells)[firstpart & !is.na(match(HAB$Date, rightsvm))], 
       col="green")
legend(x=HAB$Date[1]-10, y= 8, 
       legend=c("Observation",
                "SVM prediction"), 
       col=c("black","green"), 
       lty=c(1,0), lwd=3, pch=c(NA,1),cex=1)

plot(HAB$Date[!firstpart], 
     log10(HAB$Abundance_cells)[!firstpart],type="l",ylim=c(2,8), col = "black",
     xlab = "",
     ylab="",
     main = expression(paste(italic("K. brevis")," abundance in West Florida Shelf (SVM)")),
     axes=F,
     cex.main=1.5,cex.lab=1.5)
axis.Date(1, at=seq(as.Date("1998-01-01"), as.Date("2019-01-01"), by="year"),
          cex.axis=1.25)
axis(2,cex.axis=1.25)
abline(a=5,b=0, col = "blue")
points(HAB$Date[!firstpart & !is.na(match(HAB$Date, rightsvm))],
       log10(HAB$Abundance_cells)[!firstpart & !is.na(match(HAB$Date, rightsvm))], 
       col="green")
dev.off()
####################################################################################
tiff(file = "NB CV time series.tiff", width = 7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(2,1))
plot(HAB$Date[firstpart], log10(HAB$Abundance_cells)[firstpart],type="l",ylim=c(2,8), col = "black",
     xlab = "",
     ylab="Kb (log10(c/l))",
     main = paste0("Naive Bayes"),
     axes=F)
axis.Date(1, at=seq(as.Date("1998-01-01"), as.Date("2018-01-01"), by="year"))
axis(2)
abline(a=5,b=0, col = "blue")
points(HAB$Date[firstpart & !is.na(match(HAB$Date, rightbay))],
       log10(HAB$Abundance_cells)[firstpart & !is.na(match(HAB$Date, rightbay))], col="green")
plot(HAB$Date[!firstpart], log10(HAB$Abundance_cells)[!firstpart],type="l",ylim=c(2,8), col = "black",
     xlab = "",
     ylab="Kb (log10(c/l))",
     main = paste0("Naive Bayes"),
     axes=F)
axis.Date(1, at=seq(as.Date("1998-01-01"), as.Date("2018-01-01"), by="year"))
axis(2)
abline(a=5,b=0, col = "blue")
points(HAB$Date[!firstpart & !is.na(match(HAB$Date, rightbay))],
       log10(HAB$Abundance_cells)[!firstpart & !is.na(match(HAB$Date, rightbay))], col="green")
dev.off()
####################################################################################

# plot(HAB$Date, log10(HAB$Abundance_cells),type="l",ylim=c(2,8), col = "black",
#      xlab = "",
#      ylab="Kb (log10(c/l))",
#      main = paste0("RVM"),
#      axes=F)
# axis.Date(1, at=seq(as.Date("1998-01-01"), as.Date("2018-01-01"), by="year"))
# axis(2)
# abline(a=5,b=0, col = "blue")
# points(HAB$Date[!is.na(match(HAB$Date, rightrvm))],
#        log10(HAB$Abundance_cells)[!is.na(match(HAB$Date, rightrvm))], col="green")

####################################################################################
tiff(file = "ANN CV time series.tiff", width = 7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(2,1))
plot(HAB$Date[firstpart], log10(HAB$Abundance_cells)[firstpart],type="l",ylim=c(2,8), col = "black",
     xlab = "",
     ylab="Kb (log10(c/l))",
     main = paste0("ANN"),
     axes=F)
axis.Date(1, at=seq(as.Date("1998-01-01"), as.Date("2018-01-01"), by="year"))
axis(2)
abline(a=5,b=0, col = "blue")
points(HAB$Date[firstpart & !is.na(match(HAB$Date, rightnet))],
       log10(HAB$Abundance_cells)[firstpart & !is.na(match(HAB$Date, rightnet))], col="green")
plot(HAB$Date[!firstpart], log10(HAB$Abundance_cells)[!firstpart],type="l",ylim=c(2,8), col = "black",
     xlab = "",
     ylab="Kb (log10(c/l))",
     main = paste0("ANN"),
     axes=F)
axis.Date(1, at=seq(as.Date("1998-01-01"), as.Date("2018-01-01"), by="year"))
axis(2)
abline(a=5,b=0, col = "blue")
points(HAB$Date[!firstpart & !is.na(match(HAB$Date, rightnet))],
       log10(HAB$Abundance_cells)[!firstpart & !is.na(match(HAB$Date, rightnet))], col="green")
dev.off()
###############################################################################################3
##HAB events
# plot(HAB$Date, log10(HAB$Abundance_cells),type="l",ylim=c(6,8), col = "black",
#      xlab = "",
#      ylab="Kb (log10(c/l))",
#      main = "Karenia brevis abundance with HABs")
# abline(a=6,b=0, col = "blue")
# points(HAB$Date[right],
#        log10(HAB$Abundance_cells)[right], col="green")
# # points(HAB$Date[is.na(match(HAB$Date,dateevent))],
# #        log10(HAB$Abundance_cells)[is.na(match(HAB$Date,dateevent))], col="red")
# ##non-HAB Events
# plot(HAB$Date, log10(HAB$Abundance_cells),type="l",ylim=c(0,6), col = "black",
#      xlab = "",
#      ylab="Kb (log10(c/l))",
#      main = "Karenia brevis abundance without HABs")
# abline(a=6,b=0, col = "blue")
# points(HAB$Date[right],
#        log10(HAB$Abundance_cells)[right], col="green")
# # points(HAB$Date[is.na(match(HAB$Date,dateevent))],
# #        log10(HAB$Abundance_cells)[is.na(match(HAB$Date,dateevent))], col="red")
