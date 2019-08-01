rm(list=ls())
setwd("D:/Li_Glibert_Backup")
# setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Get Data
filenames = paste0(c("psp","cyano","asp","cfp","dsp"),"_international.csv")
dat= rep(1,83)
for(f in filenames){
  r=read.csv(f)
  print(summary(r$eventYear))
  dat=rbind(dat,r[as.numeric(r$eventYear)>=1970,])
}
dat= dat[2:nrow(dat),]
dat= as.data.frame(dat)
View(colnames(dat))
#####################################################################################
##USA, Western Europe, China, Southeastern Asia, Southern Africa
ushab = dat[as.character(dat$countryName)=="UNITED STATES",]
west_europehab = dat[!is.na(match(as.character(dat$countryName),
                                               c("PORTUGAL","FRANCE","NORWAY",
                                                 "UNITED KINGDOM","SPAIN","IRELAND",
                                                 "SWEDEN","GERMANY","NETHERLANDS",
                                                 "BELGUIM","LUXEMBOURG"))),]
canadahab = dat[as.character(dat$countryName)=="CANADA",]
southeast_asia = dat[as.character(dat$countryName)=="PHILIPPINES",]
aushab = dat[as.character(dat$countryName)=="AUSTRALIA",]
#####################################################################################
##Get fertilizer data
n_fert = read.csv("Global Nitrogen Fertilizer.csv")
p_fert = read.csv("Global Phosphorus Fertilizer.csv")
n_fert$Date = as.numeric(substr(as.character(n_fert$Data),start=8,stop=nchar(as.character(n_fert$Data))))
p_fert$Date = as.numeric(substr(as.character(p_fert$Data),start=8,stop=nchar(as.character(p_fert$Data))))
n_fert$Data = NULL
p_fert$Data= NULL
#############################################################################################################
years = seq(1970,to=2019,by=5)
uscounts = rep(0, length(years))
wecounts = rep(0,length(years))
canadacounts = rep(0,length(years))
sacounts = rep(0, length(years))

for(i in 1:length(years)){
  if(i == length(years)){
    uscounts[i] = nrow(ushab[as.numeric(as.character(ushab$eventYear))>=years[i],])
    wecounts[i] = nrow(west_europehab[as.numeric(as.character(west_europehab$eventYear))>=years[i],])
    canadacounts[i] = nrow(canadahab[as.numeric(as.character(canadahab$eventYear))>=years[i],])
    sacounts[i] = nrow(southeast_asia[as.numeric(as.character(southeast_asia$eventYear))>=years[i],])
    
  }else{
    uscounts[i] = nrow(ushab[as.numeric(as.character(ushab$eventYear))>=years[i]
                          &as.numeric(as.character(ushab$eventYear))<years[i+1],])
    wecounts[i] = nrow(west_europehab[as.numeric(as.character(west_europehab$eventYear))>=years[i]
                          &as.numeric(as.character(west_europehab$eventYear))<years[i+1],])
    canadacounts[i] = nrow(canadahab[as.numeric(as.character(canadahab$eventYear))>=years[i]
                                   &as.numeric(as.character(canadahab$eventYear))<years[i+1],])
    sacounts[i] = nrow(southeast_asia[as.numeric(as.character(southeast_asia$eventYear))>=years[i]
                                  &as.numeric(as.character(southeast_asia$eventYear))<years[i+1],])
    
  }
}
######################################################################################################
tiff("hist bloom count.tiff", width = 18, height = 6, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mar = c(5, 4, 4, 4) + 2.24,mfrow=c(1,3),mgp=c(3,0.4,0.4))  # Leave space for z axis
barplot(uscounts/5,names.arg=years,xlab="Years",main="United States",cex.main=4,
        ylab="Annual Number of HABs",cex.lab=2.5,cex.names=2.5,cex.axis=2,
        pch=24)#col="red")
barplot(canadacounts/5,names.arg=years,xlab="Years",main="Canada",cex.main=4,
        ylab="Annual Number of HABs",cex.lab=2.5,cex.names=2.5,cex.axis=2,#
        pch=25)#col="green")
barplot(wecounts/5,names.arg=years,xlab="Years",main="Western Europe",cex.main=4,
        ylab="Annual Number of HABs",cex.lab=2.5,cex.names=2.5,cex.axis=2,
        pch=21)#col="blue")
dev.off()
########################################################################################
globalcounts=rep(NA, length(years)-1)
for(i in 1:length(globalcounts)){
  globalcounts[i]=nrow(dat[as.numeric(as.character(dat$eventYear))>=years[i]&
                               as.numeric(as.character(dat$eventYear))<years[i+1],])
}
uscounts = rep(0, length(years)-1)
wecounts = rep(0,length(years)-1)
canadacounts = rep(0,length(years)-1)
sacounts = rep(0, length(years)-1)
auscounts = rep(0, length(years)-1)
for(i in 1:(length(years)-1)){
  if(i == length(years)){
    uscounts[i] = nrow(ushab[as.numeric(as.character(ushab$eventYear))>=years[i],])
    wecounts[i] = nrow(west_europehab[as.numeric(as.character(west_europehab$eventYear))>=years[i],])
    canadacounts[i] = nrow(canadahab[as.numeric(as.character(canadahab$eventYear))>=years[i],])
    sacounts[i] = nrow(southeast_asia[as.numeric(as.character(southeast_asia$eventYear))>=years[i],])
    auscounts[i] = nrow(aushab[as.numeric(as.character(aushab$eventYear))>=years[i],])
  }else{
    uscounts[i] = nrow(ushab[as.numeric(as.character(ushab$eventYear))>=years[i]
                             &as.numeric(as.character(ushab$eventYear))<years[i+1],])
    wecounts[i] = nrow(west_europehab[as.numeric(as.character(west_europehab$eventYear))>=years[i]
                                      &as.numeric(as.character(west_europehab$eventYear))<years[i+1],])
    canadacounts[i] = nrow(canadahab[as.numeric(as.character(canadahab$eventYear))>=years[i]
                                     &as.numeric(as.character(canadahab$eventYear))<years[i+1],])
    sacounts[i] = nrow(southeast_asia[as.numeric(as.character(southeast_asia$eventYear))>=years[i]
                                      &as.numeric(as.character(southeast_asia$eventYear))<years[i+1],])
    auscounts[i]=nrow(aushab[as.numeric(as.character(aushab$eventYear))>=years[i]
                                     &as.numeric(as.character(aushab$eventYear))<years[i+1],])
  }
}
#########################################################################################################
require(e1071)
mean_global = globalcounts/5
mean_global[length(mean_global)]=5/4*mean_global[length(mean_global)]

totalsvm = svm(y=mean_global,x=cbind(p_fert$Grand.Total+n_fert$Grand.Total),
            type="eps-regression")

years=n_fert$Date

tiff("global tn tp.tiff", width = 30, height = 10, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mar = c(6, 6, 6, 6) + 7.5,mfrow=c(1,3),mgp=c(5,2,0))
# barplot(uscounts,names.arg=years,xlab="Years",ylab="Counts of Algal Blooms",cex.lab=1.5,cex.names=1.5,cex.axis=2)
# par(new = TRUE)
# plot(years,n_fert$USA,type="b",lty="dashed",pch=2,xlab="",ylab="",cex.main=3,cex.lab=2.5,
#      main="United States",
#      col="red",
#      ylim=c(0,1.2e+07),
#      axes=F)http://127.0.0.1:36683/help/library/graphics/help/layout
# axis(side=4,cex.axis=1.75)
# mtext("N Fertilizer Applied (Tons)",side=4,line=3,cex=1.25)

plot(n_fert$Grand.Total,mean_global,pch=19,#bg="darkgrey",#col="red",
     xlab="",cex=8,las=1,
     ylab="",cex.main=3,cex.lab=5,cex.axis=3,axes=F)
title(ylab = "Annual No. of HABs", cex.lab = 5,
      line = 10)
title(xlab = "N Fertilizer (MT)", cex.lab = 5,
      line = 10)
axis(side=1,cex.axis=4,lwd=4,padj=0.5)
axis(side=2,cex.axis=4,lwd=4,las=1)
sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3])))
}
df=data.frame(x=(n_fert$Grand.Total-3e+07)/1e+07,
              y=mean_global/220)
fitmodel <- nls(y~a/(1 + exp(-b * (x-c))), data=df,start=list(a=1,b=0.5,c=3))
params=coef(fitmodel)
y2 <- sigmoid(params,seq(0,10,by=0.1))
lines(seq(0,10,by=0.1)*1e+07+3e+07,y2*220,lwd=4)

plot(p_fert$Grand.Total,mean_global,pch=19,#col="blue",
     xlab="",cex=8,las=1,
     #bg="darkgrey",
     ylab="",main="",cex.main=5,cex.lab=5,axes=F)
title(ylab = "Annual No. of HABs", cex.lab = 5,
      line = 10)
title(xlab = "P Fertilizer (MT)", cex.lab = 5,
      line = 10)
axis(side=1,cex.axis=4,lwd=4,padj=0.5)
axis(side=2,cex.axis=4,las=1,lwd=4)
abline(lm(mean_global~p_fert$Grand.Total),lwd=4)
text(2.5e+07,200,labels=expression(paste("r"^"2","=0.45")),cex=6)

plot(mean_global,
     predict(totalsvm,cbind(n_fert$Grand.Total+p_fert$Grand.Total)),las=1,
     pch=19,xlim=c(0.01,220),ylim=c(0.01,220),cex=8, col="black",#bg="darkgrey",
     ylab="",xlab="",cex.main=3,cex.lab=5,cex.axis=3,axes=F)
title(ylab = "Predicted No. of HABs", cex.lab = 5,
      line = 10)
title(xlab = "Observed No. of HABs", cex.lab = 5,
      line = 10)
axis(side=1,cex.axis=4,lwd=4,padj=0.5)
axis(side=2,cex.axis=4,las=1,lwd=4)
text(50,200,labels=expression(paste("r"^"2","=0.85")),cex=6)
abline(0,1,lwd=4)
dev.off()

