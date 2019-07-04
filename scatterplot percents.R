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
alldata$State = as.factor(states)
alldata$State = as.numeric(as.character(alldata$State))

# par(mfrow=c(4,3))
# for(i in 1:ncol(alldata)) {
#   plot(alldata[,i],main=colnames(alldata)[i],type="l")
#   print(c(colnames(alldata)[i],t.test(alldata[states==0,i],alldata[states==1,i])$p.value))
# }
# 
par(mfrow=c(4,3))
for(i in 1:ncol(alldata)){
  plot(alldata[,i],states,main=colnames(alldata)[i])
}

# par(mfrow=c(4,3))
# for(i in 1:ncol(alldata)){
#   print(c(colnames(alldata)[i],sd(alldata[,i])))
# }
# 
# ##TN correlations with TN
# tncol = c(18:23,32)
# for(x in 1:(length(tncol))){
#   for(y in (x+1):length(tncol)){
#     i=tncol[x]
#     j=tncol[y]
#     print(c(colnames(alldata)[i],colnames(alldata)[j],summary(lm(alldata[,j]~alldata[,i]))$r.squared))
#   }
# }
# 
# tpcol = c(24:29,33)
# for(x in 1:(length(tpcol))){
#   for(y in (x+1):length(tpcol)){
#     i=tpcol[x]
#     j=tpcol[y]
#     print(c(colnames(alldata)[i],colnames(alldata)[j],summary(lm(alldata[,j]~alldata[,i]))$r.squared))
#   }
# }

par(mfrow=c(4,4))
for(i in 1:(ncol(alldata)-1)){
  h3 = hist(alldata[,i],plot=F)
  # h1 = hist(alldata[states==0,i],main=colnames(alldata)[i], breaks=h3$breaks)
  h2 = hist(alldata[states==1,i], breaks=h3$breaks,plot=F)
  plot(h3$breaks[1:(length(h3$breaks)-1)],h2$counts/h3$counts, main=paste0(colnames(alldata)[i],":prob")) 
  plot(h3$breaks[1:(length(h3$breaks)-1)],h3$counts, main=paste0(colnames(alldata)[i],":freq")) 
  # val = quantile(alldata[,i])[2]
  # print(colnames(alldata)[i])
  # print(length(alldata[states==1,i])/length(alldata[,i]))
  # print(length(alldata[states==1&alldata[,i]<val,i])/length(alldata[alldata[,i]<val,i]))
}

