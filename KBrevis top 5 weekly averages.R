#############################################################################################
## "top 5 weekly averages.R"
## Creates weekly averages of top 5
#############################################################################################
## Read in data
rm(list=ls())
name = "West_Florida_HAB"
file = read.csv(paste0(name,".csv"))
file$X = NULL
#############################################################################################
datetime = as.Date(file$Date,format="%m/%d/%Y")
weeklyintervals = seq(as.Date("1998-01-01"),max(datetime),by=7)

db = as.data.frame(matrix(NA, nrow=length(weeklyintervals),ncol=2))
colnames(db) = c("Date","Abundance_cells")
db$Date = weeklyintervals

for(w in 1:length(weeklyintervals)){
  rightrows  = datetime>=weeklyintervals[w] & datetime < (weeklyintervals[w]+7)
  abundances = file$Abundance_cells[rightrows]
  top5orAll  = min(c(5,length(abundances)))
  db[w,2]    = mean(sort(abundances, decreasing=TRUE)[1:top5orAll])
}

write.csv(db, paste0(name,"_weekly_avg.csv"))

##Some time series plots
plot(weeklyintervals,log10(db$Abundance_cells),type="l", ylim = c(2,8))
abline(a=6,b=0)

