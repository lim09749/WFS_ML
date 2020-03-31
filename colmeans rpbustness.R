rm(list=ls())

setwd("C:/Users/xswang/HAB Research")
print_means = function(filename){
  testpred = read.csv(filename)
  testpred$X = NULL
  print(filename)
  print(colMeans(testpred))
}

print_means("robustness cross no smote.csv")
print_means("robustness cross smote.csv")
print_means("robustness block no smote.csv")
print_means("robustness block smote.csv")

