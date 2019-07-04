############################################################################################
## "USGS list variables.R"
## Features for models
rm(list = ls())
setwd("C:/Users/xswang/HAB Research/USGS weekly features")
############################################################################################


# To create metric names

# columns = c()
# files = list.files()
# for(f in files){
#   r = read.csv(f)
#   r$X <- NULL
#   columns = unique(c(columns, colnames(r)))
# }

##The different features
metric_names = c("Discharge.cubic.feet.per.second", 
                 "Discharge.tidally.filtered.cubic.feet.per.second",
                 "Dissolved.organic.matter.fluorescence..fDOM..water.in.situ.micrograms.per.liter.as.quinine.sulfate.equivalents..QSE.",
                 "Dissolved.oxygen.water.unfiltered.milligrams.per.liter",
                 "Dissolved.solids.water.filtered.estimated.by.regression.equation.milligrams.per.liter",
                 "Elevation.of.reservoir.water.surface.above.datum.feet",
                 "Gage.height.feet" , 
                 "Groundwater.level.above.NAVD",
                 "Groundwater.level.above.NGVD",
                 "Mean.water.velocity.for.discharge.computation.feet.per.second",
                 "Nitrate.plus.nitrite.water.in.situ.milligrams.per.liter.as.nitrogen",
                 "Oxidation.reduction.potential.reference.electrode.not.specified.millivolts",
                 "pH.water.unfiltered.field.standard.units",
                 "Precipitation.total.inches",
                 "Reservoir.storage.thousand.acre.feet",
                 "Salinity.water.unfiltered.parts.per.thousand",
                 "Specific.conductance.water.unfiltered.microsiemens.per.centimeter.at.25.degrees.Celsius",
                 "Stream.water.level.elevation.above.NAVD",
                 "Stream.water.level.elevation.above.NGVD",
                 "Temperature.air.degrees.Celsius",
                 "Temperature.water.degrees.Celsius",
                 "Turbidity.water.unfiltered.monochrome.near.infra.red.LED.light.780.900.nm.detection.angle.90...2.5.degrees.formazin.nephelometric.units",
                 "Wind.speed.miles.per.hour")
##########################################################################
##Keep track of their occurrence
metric_occurs = rep(0, length(metric_names))
files = list.files()
for(f in files){
  r = read.csv(f)
  r$X <- NULL
  add = 0
  for(j in 1:ncol(r)){
    for(m in metric_names){
      if(m==substr(colnames(r)[j],start=0,stop=nchar(m))){
        metric_occurs[metric_names==m] = metric_occurs[metric_names==m]+1
      }
    }
  }
}
##########################################################################
##Write csv files
file = cbind(metric_names, metric_occurs)
colnames(file) = c("feature", "Occurrence_of_each_feature")
write.csv(file, "C:/Users/xswang/HAB Research/USGS_features.csv")
