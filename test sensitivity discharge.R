numrows = 1000

dis_col   = c(13, 15, 17, 11, 31)
dis_names = colnames(alldata)[dis_col]

tn_col   = c(NA, 23, 19, 20, 32)
tn_names = rep(NA, 5)

tp_col  = c(NA, 20, 25, 26, 33)
tp_col  = rep(NA, 5)

colors = c("black","red","blue","green")

par(mfrow=c(3,4))
meanvals = colMeans(prevdata[prevdata$State=="1",1:(ncol(prevdata)-1)], na.rm=T)
for(i in dis_col){
  data = as.data.frame(matrix(rep(c(meanvals,0),numrows), byrow=T, nrow=numrows, ncol=ncol(alldata)))
  xseq = seq(from=quantile(prevdata[,i])[2],
                  to=mean(c(quantile(prevdata[,i])[4],quantile(prevdata[,i])[5])),
                  length.out=numrows)
  colnames(data) = colnames(alldata)
  data = scale(data,
               center = c(attr(scalefactors, "scaled:center"),0),
               scale  = c(attr(scalefactors, "scaled:scale"),1))
  data = as.data.frame(data)
  data$State = as.factor(data$State)
  predmat = cbind(svm_prob(HABsvm,data),
                  rvm_prob(HABrvm,data),
                  bay_prob(HABnbayes, data),
                  ann_prob(HABnet, data)) 
  
  predmat = as.data.frame(predmat)
  colnames(predmat) = c("SVM", "RVM", "BAY", "ANN")
  plot(xseq,predmat$SVM, 
       xlab=colnames(alldata)[i],
       ylab="HAB Probability [%]",
       main="SVM")
  plot(xseq,predmat$RVM, 
       xlab=colnames(alldata)[i],
       ylab="HAB Probability [%]",
       main="RVM")
  plot(xseq,predmat$BAY,
       xlab=colnames(alldata)[i],
       ylab="HAB Probability [%]",
       main="BAY")
  plot(xseq,predmat$ANN, 
       xlab=colnames(alldata)[i],
       ylab="HAB Probability [%]",
       main="ANN")
  
}
# 
# for(i in 18:23){#c(11, 13, 14, 15)
#   curr = alldata[1,]
#   curr[1:length(meanvals)] = meanvals
#   
#   preds = rep(NA, numrows)
#   sequence = seq(from=minvals[i],
#                  to=maxvals[i],
#                  length.out=numrows)
#   
#   for(j in 1:numrows){
#     #addto = rep(0, length(curr)-1)
#     addto = curr[,1:(length(curr)-1)]
#     addto[i] = sequence[j]
#     preds[j] = attr(predict(HABsvm,curr[,1:(length(curr)-1)]+addto,
#                             probability=T),"probabilities")[1,2]
#   }
#   plot(undoscale(i, curr[1,i]+sequence, scalefactors),
#        preds, type="l",main=colnames(alldata)[i],#main=names[i],
#        xlab = expression(paste('Discharge (m'^3,"s"^-1,")")),
#        ylab = "HAB probability[%]",axes=F,
#        ylim = c(0, 1),
#        cex.main=3,cex.lab=2,
#        lwd=4)
#   axis(1,cex.axis=2,at=axTicks(side=1),labels=signif(10^axTicks(side=1)*0.0283168,digits=1))
#   axis(2,cex.axis=2)
#   
# }
# dev.off()
