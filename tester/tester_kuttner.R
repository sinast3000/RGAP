# clear
rm(list=ls())


####### Settings ------------------------------------------------------

# load package
devtools::load_all()

model <- f <- fBayes <- fSTN <- list()

####### data ----------------------------------------------------------
country <- "Netherlands"
# data("gap")
# tsList <- as.list(gap[["Netherlands"]][,c("cpih","gdp")])
tslBase <- fetchAmecoData()
tsList <- as.list(tslBase[["Netherlands"]][,c("cpih","gdp")])

tsList$infl <- diff(tsList$cpih)

####### Kuttner model specification -----------------------------------

# default
model[[1]] <- KuttnerModel(tsl = tsList)
model[[1]]
# cycle lags
model[[2]] <- KuttnerModel(tsl = tsList, cycleLag = 0:3)
model[[2]]
KuttnerModel(tsl = tsList, cycleLag = 0:11) # error: too many
# cycle
model[[3]] <- KuttnerModel(tsl = tsList, cycle = "RAR2")
model[[3]]
model[[4]] <- KuttnerModel(tsl = tsList, cycle = "AR1")
model[[4]]
# trend
model[[5]] <- KuttnerModel(tsl = tsList, trend = "RW2")
model[[5]]
model[[6]] <- KuttnerModel(tsl = tsList, trend = "DT")
model[[6]]
# errorAMRA
model[[6]] <- KuttnerModel(tsl = tsList, inflErrorARMA = c(0, 1))
model[[6]]
KuttnerModel(tsl = tsList, inflErrorARMA = c(2, 1)) # error: AR should be 0

#
# model <- KuttnerModel(tsl = tsList, trend = "RW2",
#                       cycleLag = 1, cycle = "AR2")


####### Kuttner MLE fitting -------------------------------------------

for (k in 1:length(model)) {
  parRestr <- initializeRestr(model = model[[k]], type = "hp")
  # parRestr <- initializeRestr(model = model[[k]])
  f[[k]] <- fit(model = model[[k]], parRestr = parRestr)
  plot(f[[k]])
}

####### Kuttner MLE fitting - signal to noise -------------------------

for (k in 1:length(model)) {
  parRestr <- initializeRestr(model = model[[k]])
  fSTN[[k]] <- fit(model = model[[k]], parRestr = parRestr, signalToNoise = 0.2)
  plot(fSTN[[k]])
}


