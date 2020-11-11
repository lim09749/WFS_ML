##new CV metrics.R
rm(list=ls())
setwd("C:/Users/xswang/HAB Research/")
#####################################################################################
##Get Data
cv = read.csv("cv_without_satellite_keeping_data.csv")
cv$X = NULL
##Find proportion of data state = 1
HAB = read.csv("West_Florida_HAB_weekly_avg.csv")
HAB$X = NULL
HAB$Abundance_cells = as.numeric(as.character(HAB$Abundance_cells))
HAB$Date = as.Date(HAB$Date)
states = rep(0, nrow(HAB))
states[HAB$Abundance_cells>=10^5] = 1
prop_1 = length(which(states==1))/length(states)

###############################################################################################3
##Create confusion matrix

recall_ = function(total_acc, non_acc, hab_acc, prop_1){
  coeff = matrix(c(1,0,0,1,
                   1-non_acc,-non_acc,0,0,
                   0,0,-hab_acc,1-hab_acc,
                   0,0,1,1), nrow=4, ncol=4, byrow = T)
  b = c(total_acc, 0, 0, prop_1)
  sol = solve(coeff,b)
  sol[4]/sum(sol[c(3,4)])
}

precision_ = function(total_acc, non_acc, hab_acc, prop_1){
  coeff = matrix(c(1,0,0,1,
                   1-non_acc,-non_acc,0,0,
                   0,0,-hab_acc,1-hab_acc,
                   0,0,1,1), nrow=4, ncol=4, byrow = T)
  b = c(total_acc, 0, 0, prop_1)
  sol = solve(coeff,b)
  sol[4]/sum(sol[c(2,4)])
}

#####################################################################################################
recall = as.data.frame(matrix(NA, nrow = 10, ncol=4))
rownames(recall) = paste0("Trial.", 1:10)
colnames(recall) = c("SVM","RVM","ANN","NB")

precision = as.data.frame(matrix(NA, nrow = 10, ncol=4))
rownames(precision) = paste0("Trial.", 1:10)
colnames(precision) = c("SVM","RVM","ANN","NB")

for(i in c(1:10)){
  print(i)
  if(i!=6){
    recall[i, ]= c(recall_(total_acc = cv$SVM.Total.Accuracy[i], 
                           hab_acc = cv$SVM.HAB.Accuracy[i], 
                           non_acc = cv$SVM.non.HAB.Accuracy[i], prop_1),
                   recall_(total_acc = cv$RVM.Total.Accuracy[i], 
                           hab_acc = cv$RVM.HAB.Accuracy[i], 
                           non_acc = cv$RVM.non.HAB.Accuracy[i], prop_1),
                   recall_(total_acc = cv$ANN.Total.Accuracy[i], 
                           hab_acc = cv$ANN.HAB.Accuracy[i], 
                           non_acc = cv$ANN.non.HAB.Accuracy[i], prop_1),
                   recall_(total_acc = cv$BAY.Total.Accuracy[i], 
                           hab_acc = cv$BAY.HAB.Accuracy[i], 
                           non_acc = cv$BAY.non.HAB.Accuracy[i], prop_1))
    precision[i, ]= c(precision_(total_acc = cv$SVM.Total.Accuracy[i], 
                                 hab_acc = cv$SVM.HAB.Accuracy[i], 
                                 non_acc = cv$SVM.non.HAB.Accuracy[i], prop_1),
                      precision_(total_acc = cv$RVM.Total.Accuracy[i], 
                                 hab_acc = cv$RVM.HAB.Accuracy[i], 
                                 non_acc = cv$RVM.non.HAB.Accuracy[i], prop_1),
                      precision_(total_acc = cv$ANN.Total.Accuracy[i], 
                                 hab_acc = cv$ANN.HAB.Accuracy[i], 
                                 non_acc = cv$ANN.non.HAB.Accuracy[i], prop_1),
                      precision_(total_acc = cv$BAY.Total.Accuracy[i], 
                                 hab_acc = cv$BAY.HAB.Accuracy[i], 
                                 non_acc = cv$BAY.non.HAB.Accuracy[i], prop_1))
  }else{
    recall[i, ]= c(recall_(total_acc = cv$SVM.Total.Accuracy[i], 
                           hab_acc = cv$SVM.HAB.Accuracy[i], 
                           non_acc = cv$SVM.non.HAB.Accuracy[i], prop_1),
                   0,
                   recall_(total_acc = cv$ANN.Total.Accuracy[i], 
                           hab_acc = cv$ANN.HAB.Accuracy[i], 
                           non_acc = cv$ANN.non.HAB.Accuracy[i], prop_1),
                   0)
    precision[i, ]= c(precision_(total_acc = cv$SVM.Total.Accuracy[i], 
                                 hab_acc = cv$SVM.HAB.Accuracy[i], 
                                 non_acc = cv$SVM.non.HAB.Accuracy[i], prop_1),
                      0,
                      precision_(total_acc = cv$ANN.Total.Accuracy[i], 
                                 hab_acc = cv$ANN.HAB.Accuracy[i], 
                                 non_acc = cv$ANN.non.HAB.Accuracy[i], prop_1),
                      0)
  }
}
################################################################################################
## Save data
dfr = rbind(recall,colMeans(recall, na.rm=T))
rownames(dfr)[11] = "Avg"
write.csv(dfr, "recall_block.csv")
dfp = rbind(precision,colMeans(precision, na.rm=T))
rownames(dfp)[11] = "Avg"
write.csv(dfp, "precision_block.csv")

####################################################################################################3
##K-fold cross validation

##Block
cv[11,] = c(0.63, 0.85, 0.78, 
            0.54, 0.71, 0.66,
            0.65, 0.56, 0.59,
            0.42, 0.76, 0.65)
print("SVM")
recall_(
  total_acc = cv$SVM.Total.Accuracy[11],
  hab_acc = cv$SVM.HAB.Accuracy[11],
  non_acc = cv$SVM.non.HAB.Accuracy[11],
  prop_1 = prop_1
)
precision_(
  total_acc = cv$SVM.Total.Accuracy[11],
  hab_acc = cv$SVM.HAB.Accuracy[11],
  non_acc = cv$SVM.non.HAB.Accuracy[11],
  prop_1 = prop_1
)

print("RVM")

recall_(
  total_acc = cv$RVM.Total.Accuracy[11],
  hab_acc = cv$RVM.HAB.Accuracy[11],
  non_acc = cv$RVM.non.HAB.Accuracy[11],
  prop_1 = prop_1
)
precision_(
  total_acc = cv$RVM.Total.Accuracy[11],
  hab_acc = cv$RVM.HAB.Accuracy[11],
  non_acc = cv$RVM.non.HAB.Accuracy[11],
  prop_1 = prop_1
)

print("ANN")

recall_(
  total_acc = cv$ANN.Total.Accuracy[11],
  hab_acc = cv$ANN.HAB.Accuracy[11],
  non_acc = cv$ANN.non.HAB.Accuracy[11],
  prop_1 = prop_1
)
precision_(
  total_acc = cv$ANN.Total.Accuracy[11],
  hab_acc = cv$ANN.HAB.Accuracy[11],
  non_acc = cv$ANN.non.HAB.Accuracy[11],
  prop_1 = prop_1
)

print("NB")
recall_(
  total_acc = cv$BAY.Total.Accuracy[11],
  hab_acc = cv$BAY.HAB.Accuracy[11],
  non_acc = cv$BAY.non.HAB.Accuracy[11],
  prop_1 = prop_1
)

precision_(
  total_acc = cv$BAY.Total.Accuracy[11],
  hab_acc = cv$BAY.HAB.Accuracy[11],
  non_acc = cv$BAY.non.HAB.Accuracy[11],
  prop_1 = prop_1
)
