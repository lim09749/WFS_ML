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
alldata = read.csv("alldata.csv")
alldata$X = NULL
alldata = as.data.frame(scale(alldata))

alldata$State = as.factor(states)
alldata$State = as.numeric(as.character(alldata$State))

prevdata = alldata
alldata$SSH=NULL

pred = read.csv("testpred_block_no_SMOTE.csv")
pred$X = NULL
###############################################################################
#####################################################################################
##Get HAB values (use fact that all of peace river data is unique)
##Time series plot func
plot_ts = function(file_name, width, height, date, abundance, model_type, right){
  tiff(file = file_name, width = width, height = height, units = "in",
       pointsize=10, res = 300, compression = c("lzw"))
  plot(date[500:length(date)], 
       log10(abundance)[500:length(date)],type="l",ylim=c(2,8.4), col = "black",
       xlab = "",
       ylab="Kb (log10(c/l))",
       main = "",#model_type,#bquote(paste(italic("K. brevis")," abundance in West Florida Shelf (",.(model_type),")")),
       axes=F,
       cex.main=1.5,cex.lab=1.5)
  axis.Date(1, at=seq(as.Date("2008-01-01"), as.Date("2019-01-01"), by="year"),
            cex.axis=1.25)
  axis(2,cex.axis=1.25)
  abline(a=5,b=0, col = "black", lty='dashed')
  print(date[!is.na(match(date, right))])
  points(date[!is.na(match(date, right))][281:length(date[!is.na(match(date, right))])],
         log10(abundance)[!is.na(match(date, right))][281:length(date[!is.na(match(date, right))])], 
         col="green")
  # points(date[is.na(match(date, right))],
  #        log10(abundance)[is.na(match(date, right))],
  #        col="orange")
  # legend(x=date[1]-100, y= 8.25, 
  #        legend=c("Observation",
  #                 paste0(model_type," correct predictions")
  #                 ), 
  #        col=c("black","green"), ncol=1,
  #        lty=c(1,0), lwd=3, pch=c(NA,1),cex=1)
  # 
  dev.off()
  
}
#####################################################################################
plot_ts("new SVM supplementary Mar 20.tiff", width=18, height=4.25, 
        date=HAB$Date, 
        abundance=HAB$Abundance_cells, 
        model_type="SVM", 
        right=as.Date(pred$Date[pred$Real==pred$SVM]))

plot_ts("new RVM supplementary Mar 20.tiff", width=18, height=4.25, 
        date=HAB$Date, 
        abundance=HAB$Abundance_cells, 
        model_type="RVM", 
        right=as.Date(pred$Date[pred$Real==pred$RVM]))

plot_ts("new NB supplementary Mar 20.tiff", width=18, height=4.25, 
        date=HAB$Date, 
        abundance=HAB$Abundance_cells, 
        model_type="NB", 
        right=as.Date(pred$Date[pred$Real==pred$NB]))

plot_ts("new ANN supplementary Mar 20.tiff", width=18, height=4.25, 
        date=HAB$Date, 
        abundance=HAB$Abundance_cells, 
        model_type="ANN", 
        right=as.Date(pred$Date[pred$Real==pred$NET]))

############################################################################################

plot_ts("new RVM supplementary Mar 20_1.tiff", width=9, height=4.25, 
        date=HAB$Date, 
        abundance=HAB$Abundance_cells, 
        model_type="RVM", 
        right=as.Date(pred$Date[pred$Real==pred$RVM]))


plot_ts("new RVM supplementary Mar 20_2.tiff", width=9, height=4.25, 
        date=HAB$Date, 
        abundance=HAB$Abundance_cells, 
        model_type="RVM", 
        right=as.Date(pred$Date[pred$Real==pred$RVM]))
