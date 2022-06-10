# clear
rm(list=ls())


####### Settings ------------------------------------------------------

# load package
devtools::load_all()

model <- f <- fBayes <- fSTN <- list()

####### data ----------------------------------------------------------
country <- "Italy"
# data("gap")
# tsList <- amecoData2input(gap[[country]], alpha = 0.65)
tslBase <- fetchAmecoData()
tsList <- amecoData2input(tslBase[[country]], alpha = 0.65)

####### TFP model specification ---------------------------------------

# default
model[[1]] <- TFPmodel(tsl = tsList)
model[[1]]
# cycle lags
model[[2]] <- TFPmodel(tsl = tsList, cycleLag = 0:5)
model[[2]]
TFPmodel(tsl = tsList, cycleLag = 0:11) # error: too many
# cycle
model[[3]] <- TFPmodel(tsl = tsList, cycle = "RAR2")
model[[3]]
# trend
model[[4]] <- TFPmodel(tsl = tsList, trend = "RW1")
model[[4]]
model[[5]] <- TFPmodel(tsl = tsList, trend = "RW2")
model[[5]]
# errorAMRA
model[[6]] <- TFPmodel(tsl = tsList, cubsErrorARMA = c(2, 0))
model[[6]]
TFPmodel(tsl = tsList, cubsErrorARMA = c(2, 1)) # error: MA should be 0


# model <- TFPmodel(tsl = tsList, trend = "RW2", cycle = "RAR2",
#                   cycleLag = 0, cubsErrorARMA = c(0,0))

####### TFP MLE fitting -----------------------------------------------

for (k in 1:length(model)) {
  parRestr <- initializeRestr(model = model[[k]], type = "hp")
  # parRestr <- initializeRestr(model = model[[k]])
  f[[k]] <- fit(model = model[[k]], parRestr = parRestr)
  plot(f[[k]])
}

####### TFP bayesian fitting ------------------------------------------

for (k in 1:length(model)) {
  # prior <- initializePrior(model = model[[k]])
  prior <- initializePrior(model = model[[k]], MLEfit = f[[k]], MLE = TRUE)
  fBayes[[k]] <- fit(model = model[[k]], method = "bayesian", prior = prior,
                          R = 1000, thin = 2, MLEfit = f[[k]])
  # error for NKP: bayesian methods not available
  plot(fBayes[[k]])

}

####### TFP MLE fitting - signal to noise -----------------------------

for (k in 1:length(model)) {
  parRestr <- initializeRestr(model = model[[k]])
  fSTN[[k]] <- fit(model = model[[k]], parRestr = parRestr, signalToNoise = 0.05)
  # fSTN[[k]] <- fit(model = model[[k]], parRestr = parRestr, R = 1000, signalToNoise = 0.05, method = "bayesian")
  # # warning that signal to noise not applicable for method "bayesian"
  plot(fSTN[[k]])
}

