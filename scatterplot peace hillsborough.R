rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Read in data
HAB = read.csv("C:/Users/xswang/HAB Research/West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)

alldata = read.csv("totaldata.csv")
alldata$X = NULL
for(i in 1:ncol(alldata)) alldata[is.na(alldata[,i]),i]=mean(alldata[,i],na.rm=T)
colors = rep("black",nrow(HAB))
colors[HAB$Abundance_cells>=10^5]="red"

###########################################################################################
tiff(file = "scatterplot for TP and TN.tiff", width = 7.5, height = 7.5, units = "in",
     pointsize=10, res = 300, compression = c("lzw"))
par(mfrow=c(4,1))
plot(HAB$Date, log10(alldata$Peace_TN_mg),type="p",
     xlab = "Date",
     ylab="Total Nitrogen (log10(mg))",
     col=colors,
     main = "Peace River")
plot(HAB$Date, log10(alldata$Peace_TP_mg),type="p",
     xlab = "Date",
     ylab="Total Phosphorus (log10(mg))",
     col=colors,
     main = "Peace River")

plot(HAB$Date, log10(alldata$Hillsborough_TN_mg),type="p",
     xlab = "Date",
     ylab="Total Nitrogen (log10(mg))",
     col=colors,
     main = "Hillsborough River")
plot(HAB$Date, log10(alldata$Hillsborough_TP_mg),type="p",
     xlab = "Date",
     ylab="Total Phosphorus (log10(mg))",
     col=colors,
     main = "Hillsborough River")
dev.off()
