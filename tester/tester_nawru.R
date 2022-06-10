# clear
rm(list=ls())


####### Settings ------------------------------------------------------

# load package
devtools::load_all()

model <- f <- fBayes <- fSTN <- list()

####### data ----------------------------------------------------------
country <- "France"
# data("gap")
# tsList <- amecoData2input(gap[[country]], alpha = 0.65)
tslBase <- fetchAmecoData()
tsList <- amecoData2input(tslBase[[country]], alpha = 0.65)

####### NAWRU model specification -------------------------------------
D <- matrix(c(2, 2, 2), 1, 3, byrow = TRUE)
L <- matrix(c(0, 0, 0), 1, 3, byrow = TRUE)
exoType <- initializeExo(varNames = c("ws", "prod","tot"), D = D, L = L)


# default
model[[1]] <- NAWRUmodel(tsl = tsList)
model[[1]]
# exo
model[[2]] <- NAWRUmodel(tsl = tsList, exoType = exoType)
model[[2]]
# cycle lags
model[[3]] <- NAWRUmodel(tsl = tsList, cycleLag = 0:5, exoType = exoType)
model[[3]]
NAWRUmodel(tsl = tsList, cycleLag = 0:11, exoType = exoType) # error: too many
# cycle
NAWRUmodel(tsl = tsList, cycle = "RAR2", exoType = exoType) # error: only AR1 or AR2
# trend
model[[4]] <- NAWRUmodel(tsl = tsList, trend = "DT", exoType = exoType)
model[[4]]
model[[5]] <- NAWRUmodel(tsl = tsList, trend = "RW1", exoType = exoType)
model[[5]]
# errorAMRA
model[[6]] <- NAWRUmodel(tsl = tsList, pcErrorARMA = c(2, 0), exoType = exoType)
model[[6]]
NAWRUmodel(tsl = tsList, pcErrorARMA = c(2, 1), exoType = exoType) # error: MA should be 0
# NKP
model[[7]] <- NAWRUmodel(tsl = tsList, type = "NKP", exoType = exoType) # 2 warnings: cycleLag changed and no exo variables
model[[7]]
model[[8]] <- NAWRUmodel(tsl = tsList, type = "NKP", cycleLag = 3) #  warning: cycleLag changed
model[[8]]

# model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2",
#                     type = "TKP", cycleLag = 0:1, exoType = exoType)


####### NAWRU MLE fitting ---------------------------------------------

for (k in 1:length(model)) {
  parRestr <- initializeRestr(model = model[[k]], type = "hp")
  # parRestr <- initializeRestr(model = model[[k]])
  f[[k]] <- fit(model = model[[k]], parRestr = parRestr)
  plot(f[[k]])
}

####### NAWRU bayesian fitting ----------------------------------------

for (k in 1:length(model)) {
  prior <- initializePrior(model = model[[k]])
  # prior <- initializePrior(model = model[[k]], MLEfit = fit[[k]], MLE = TRUE)
  fBayes[[k]] <- fit(model = model[[k]], method = "bayesian", prior = prior,
                       R = 1000, thin = 2, MLEfit = fit[[k]])
  # error for NKP: bayesian methods not available
  plot(fBayes[[k]])

}

####### NAWRU MLE fitting - signal to noise ---------------------------

for (k in 1:length(model)) {
  parRestr <- initializeRestr(model = model[[k]])
  fSTN[[k]] <- fit(model = model[[k]], parRestr = parRestr, signalToNoise = 0.05)
  # fSTN[[k]] <- fit(model = model[[k]], parRestr = parRestr, R = 1000, signalToNoise = 0.05, method = "bayesian")
  # # warning that signal to noise not applicable for method "bayesian"
  plot(fSTN[[k]])
}

