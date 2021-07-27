# --------------------------------------------------------------------------------------------------- #

# clear
rm(list=ls())

# --------------------------------------------------------------------------------------------------- #
# hidden data

# ----- indicators for cubs
country <- data.frame(abbreviation = c("EU","EA","BE","BG","CZ","DK","DE","EE","IE","EL","ES","FR","HR","IT","CY","LV","LT","LU","HU","MT","NL","AT","PL","PT","RO","SI","SK","FI","SE","UK","ME","MK","AL","RS","TR"),
                      name = c("European Union (current composition)","Euro area","Belgium","Bulgaria","Czechia","Denmark","Germany","Estonia","Ireland","Greece","Spain","France",
                                  "Croatia","Italy","Cyprus","Latvia","Lithuania","Luxembourg","Hungary","Malta","Netherlands","Austria","Poland","Portugal",
                                  "Romania","Slovenia","Slovak Republic","Finland","Sweden","United Kingdom","Montenegro","North Macedonia","Albania","Serbia","Turkey"))
source <- data.frame(file = c("main_indicators_nace2.xlsx",
                              "main_indicators_nace2.xlsx",
                              "industry_total_sa_nace2.xlsx"),
                     sheet = c("MONTHLY",
                               "MONTHLY",
                               "INDUSTRY QUARTERLY"),
                     freq = c(12,
                              12,
                              4),
                     key = c("XXX.SERV",
                             "XXX.BUIL",
                             "INDU.XXX.TOT.13.QPS.Q"),
                     type = c("SERV",
                              "BUIL",
                              "INDU"),
                     datastream =c("TTS99BQ",
                                   "41.COBQ",
                                   "CAPUTLQ"))
indicatorList <- list(country = country, source = source)



# ----- ameco data
# ameco ----> source
# load file
filename <- "ameco_autumn2018.xlsx"
dfameco <- openxlsx::read.xlsx(paste0("data-raw/",filename),
                               sheet = "Sheet1",
                               startRow = 1)
# get country abbreviations and names
tmp_abbreviation <- unique(gsub("\\..*", "", dfameco$CODE))
dfameco$tmp <- gsub("\\..*", "", dfameco$CODE)
tmp_name <- NULL
for (abbr in tmp_abbreviation) {
  tmp_name <- c(tmp_name,dfameco$COUNTRY[dfameco$tmp==abbr][1])
}
tmp_name <- gsub("\"linked\" ","",tmp_name)
dfameco$tmp <- NULL
country <- data.frame(abbreviation = tmp_abbreviation,
                      name = tmp_name)
source <- data.frame(key =  c("OVGD",
                               "UVGD",
                               "OKND",
                               "NETD",
                               "NLHA",
                               "ZUTN",
                               "UWCD",
                               "PVGD",
                               "PCPH",
                               "NPAN",
                               "NLHT",
                               "NETN",
                               "NWTD",
                               "PLCD",
                               "OVGM",
                               "OVG5",
                               "OVG4",
                               "ZCPIH",
                               "ZCPIN"),
                    finalname =  c("gdp",
                                   "ngdp",
                                   "k",
                                   "et",
                                   "ahours",
                                   "ur",
                                   "wtotal",
                                   "gdpdefl",
                                   "pconsp",
                                   "popw",
                                   "l",
                                   "etd",
                                   "eet",
                                   "nulc",
                                   "vaindu",
                                   "vaserv",
                                   "vabuil",
                                   "cpih",
                                   "cpin"),
                    keynumber =  c(".1.1.0.0.OVGD",
                                   ".1.0.0.0.UVGD",
                                   ".1.0.0.0.OKND",
                                   ".1.0.0.0.NETD",
                                   ".1.0.0.0.NLHA",
                                   ".1.0.0.0.ZUTN",
                                   ".1.0.0.0.UWCD",
                                   ".3.1.0.0.PVGD",
                                   ".3.1.0.0.PCPH",
                                   ".1.0.0.0.NPAN",
                                   ".1.0.0.0.NLHT",
                                   ".1.0.0.0.NETN",
                                   ".1.0.0.0.NWTD",
                                   ".3.1.0.0.PLCD",
                                   ".1.1.0.0.OVGM",
                                   ".1.1.0.0.OVG5",
                                   ".1.1.0.0.OVG4",
                                   ".1.0.0.0.ZCPIH",
                                   ".3.0.0.0.ZCPIN"))

source$title <- dfameco$TITLE[dfameco$CODE %in% paste0("FRA",source$keynumber)][sapply(paste0("FRA",source$keynumber),match,as.character(dfameco$CODE[dfameco$CODE %in% paste0("FRA",source$keynumber)]))]
source$unit = dfameco$UNIT[dfameco$CODE %in% paste0("FRA",source$keynumber)][sapply(paste0("FRA",source$keynumber),match,as.character(dfameco$CODE[dfameco$CODE %in% paste0("FRA",source$keynumber)]))]

ameco <- list(country = country, source = source)

# save hidden data
# save(indicatorList,ameco, file = "R/sysdata.rda")

# ----- model options

# initialize
cycle <- trend <- error <- cubs <- pcInd <- infl <- list()

# cycle
cycle$AR1 <- data.frame(equation    = "cycle",
                        variant     = "AR1",
                        sysMatrix   = c("T","Q"),
                        variableRow = c("cycle", "cycle"),
                        varName     = c("cPhi1", "cSigma"),
                        lowerBound  = c(NA, 0),
                        upperBound  = NA,
                        statRestr   = c("AR1", ""),
                        mean        = c(0.5, 5e-4),
                        std         = c(0.1, 5e-4),
                        distribution = c("normal", "invgamma"),
                        stringsAsFactors = FALSE)

cycle$AR2 <- data.frame(equation    = "cycle",
                        variant     = "AR2",
                        sysMatrix   = c("T", "T", "Q"),
                        variableRow = c("cycle", "cycleLag1", "cycle"),
                        varName     = c("cPhi1", "cPhi2", "cSigma"),
                        lowerBound  = c(NA, NA, 0),
                        upperBound  = NA,
                        statRestr   = c("AR2", "AR2", ""),
                        mean        = c(0.5, -0.1, 5e-4),
                        std         = c(0.1, 0.1, 5e-4),
                        distribution = c("normal", "normal", "invgamma"),
                        stringsAsFactors = FALSE)

cycle$RAR2 <- data.frame(equation    = "cycle",
                         variant     = "RAR2",
                         sysMatrix   = c("T", "T", "Q"),
                         variableRow = c("cycle", "cycleLag1", "cycle"),
                         varName     = c("cA", "cTau", "cSigma"),
                         lowerBound  = c(0.01, 2.01, 0),
                         upperBound  = c(0.99, 31.99, NA),
                         statRestr   = c(""),
                         mean        = c(0.42, 8, 5e-4),
                         std         = c(0.17, 3.5, 5e-4),
                         distribution = c("beta", "beta", "invgamma"),
                         stringsAsFactors = FALSE)

cycle <- rbind(cycle[[1]], cycle[[2]], cycle[[3]])

# trend
trend$RW1 <- data.frame(equation    = "trend",
                        variant     = "RW1",
                        sysMatrix   = c("T","Q"),
                        variableRow = c("trendDrift", "trend"),
                        varName     = c("tdConst", "tSigma"),
                        lowerBound  = c(NA, 0),
                        upperBound  = NA,
                        statRestr   = c(""),
                        mean        = c(0.15, 5e-6),
                        std         = c(0.01, 5e-6),
                        distribution = c("normal", "invgamma"),
                        stringsAsFactors = FALSE)
trend$RW2 <- data.frame(equation    = "trend",
                        variant     = "RW2",
                        sysMatrix   = c("Q","Q"),
                        variableRow = c("trend", "trendDrift"),
                        varName     = c("tSigma", "tdSigma"),
                        lowerBound  = 0,
                        upperBound  = NA,
                        statRestr   = c(""),
                        mean        = c(0, 5e-6),
                        std         = c(0, 5e-6),
                        distribution = c("invgamma", "invgamma"),
                        stringsAsFactors = FALSE)
trend$DT <- data.frame(equation    = "trend",
                       variant     = "DT",
                       sysMatrix   = c("Tt", "Tt", "Q", "Q"),
                       variableRow = c("const", "trendDrift", "trend", "trendDrift"),
                       varName     = c("tdOmega", "tdPhi","tSigma", "tdSigma"),
                       # lowerBound  = c(0, 0, 0, 0),
                       # upperBound  = c(0.03, 0.99, NA, NA),
                       lowerBound  = c(NA, 0, 0, 0),
                       upperBound  = c(NA, 0.99, NA, NA),
                       statRestr   = c(""),
                       mean        = c(0.015, 0.8, 0, 5e-6),
                       std         = c(0.01, 0.24, 0, 5e-6),
                       distribution = c("normal", "normal", "invgamma", "invgamma"),
                       stringsAsFactors = FALSE)
trend <- rbind(trend[[1]], trend[[2]], trend[[3]])

# E2 error
error$base <- data.frame(equation     = "E2error",
                          variant     = "base",
                          sysMatrix   = "Q",
                          variableRow = "E2errorL0",
                          varName     = "E2Sigma",
                          lowerBound  = 0,
                          upperBound  = NA,
                          statRestr   = "",
                          mean        = 5e-5,
                          std         = 5e-5,
                         distribution = "invgamma",
                         stringsAsFactors = FALSE)
error$MA <- data.frame(equation    = "E2error",
                       variant     = "errorMA",
                       sysMatrix   = rep("T", 5),
                       variableRow = paste0("E2innoL", 0:4),
                       varName     = paste0("E2ErrGamma", 1:5),
                       lowerBound  = NA,
                       upperBound  = NA,
                       statRestr   = c(""),
                       mean        = rep(0, 5),
                       std         = rep(2, 5),
                       distribution = rep("normal", 5),
                       stringsAsFactors = FALSE)
error$AR1 <- data.frame(equation    = "E2error",
                        variant     = "errorAR1",
                        sysMatrix   = "T",
                        variableRow = "E2errorL0",
                        varName     = "E2ErrPhi1",
                        lowerBound  = NA,
                        upperBound  = NA,
                        statRestr   = "AR1",
                        mean        = 0,
                        std         = 0.4,
                        distribution = "normal",
                        stringsAsFactors = FALSE)
error$AR2 <- data.frame(equation    = "E2error",
                        variant     = "errorAR2",
                        sysMatrix   = c("T", "T"),
                        variableRow = c("E2errorL0", "E2errorL1"),
                        varName     = c("E2ErrPhi1", "E2ErrPhi2"),
                        lowerBound  = NA,
                        upperBound  = NA,
                        statRestr   = "AR2",
                        mean        = c(0, 0),
                        std         = c(0.4, 0.4),
                        distribution = c("normal", "normal"),
                        stringsAsFactors = FALSE)
error <- rbind(error[[1]], error[[2]], error[[3]], error[[4]])


# cubs
cubs$base <- data.frame(equation    = "cubs",
                        variant     = "base",
                        sysMatrix   = c("Z", "Z"),
                        variableRow = c("const", "cycle"),
                        varName     = c("cuConst", "cuC0"),
                        # lowerBound  = c(-0.1, NA),
                        # upperBound  = c(0.1, NA),
                        lowerBound  = c(NA, NA),
                        upperBound  = c(NA, NA),
                        statRestr   = c(""),
                        mean        = c(0, 1.4),
                        std         = c(0.03, 0.7),
                        distribution = c("normal", "normal"),
                        stringsAsFactors = FALSE)
cubs$cycleLag <- data.frame(equation    = "cubs",
                            variant     = "cycleLag",
                            sysMatrix   = rep("Z", 10),
                            variableRow = paste0("cycleLag", 1:10),
                            varName     = paste0("cuC", 1:10),
                            lowerBound  = NA,
                            upperBound  = NA,
                            statRestr   = c(""),
                            mean        = rep(0, 10),
                            std         = rep(2, 10),
                            distribution = rep("normal", 10),
                            stringsAsFactors = FALSE)
cubs$cubsAR1 <- data.frame(equation    = "cubs",
                            variant     = "cubsAR1",
                            sysMatrix   = "exo",
                            variableRow = "cubsAR1",
                            varName     = "cuPhi",
                            lowerBound  = NA,
                            upperBound  = NA,
                            statRestr   = "AR1",
                            mean        = 0,
                            std         = 2,
                            distribution = "normal",
                            stringsAsFactors = FALSE)
cubs$cubsAR2 <- data.frame(equation    = "cubs",
                            variant     = "cubsAR2",
                            sysMatrix   = rep("exo", 2),
                            variableRow = paste0("cubsAR", 1:2),
                            varName     = paste0("cuPhi", 1:2),
                            lowerBound  = NA,
                            upperBound  = NA,
                            statRestr   = "AR2",
                            mean        = rep(0, 2),
                            std         = rep(2, 2),
                            distribution = c("normal", "normal"),
                            stringsAsFactors = FALSE)
cubs <- rbind(cubs[[1]], cubs[[2]], cubs[[3]], cubs[[4]])


# nawru ################ not completely done yet
# e.g. need to check whether exo input can be switched
pcInd$base <- data.frame(equation   = "pcInd",
                          variant     = "base",
                          sysMatrix   = c("Z", "Z"),
                          variableRow = c("const", "cycle"),
                          varName     = c("pcConst", "pcC0"),
                          # lowerBound  = c(-0.1, -10),
                          # upperBound  = c(0.1, 0),
                         lowerBound  = c(NA, NA),
                         upperBound  = c(NA, NA),
                         statRestr   = c(""),
                          mean        = c(0, -1),
                          std         = c(0.1, 0.5),
                         distribution = c("normal", "normal"),
                         stringsAsFactors = FALSE)
pcInd$exo <- data.frame(equation   = "pcInd",
                         variant     = "exo",
                         sysMatrix   = c("exo"),
                         variableRow = c("constExo"),
                         varName     = c("pcXXX"),
                         # lowerBound  = -10,
                         # upperBound  = 10,
                        lowerBound  = NA,
                        upperBound  = NA,
                         statRestr   = c(""),
                         mean        = 0,
                         std         = 2,
                        distribution = "normal",
                        stringsAsFactors = FALSE)
pcInd$cycleLag <- data.frame(equation    = "pcInd",
                            variant     = "cycleLag",
                            sysMatrix   = rep("Z", 10),
                            variableRow = paste0("cycleLag", 1:10),
                            varName     = paste0("pcC", 1:10),
                            lowerBound  = NA,
                            upperBound  = NA,
                            statRestr   = c(""),
                            mean        = rep(0, 10),
                            std         = rep(2, 10),
                            distribution = rep("normal", 10),
                            stringsAsFactors = FALSE)
pcInd$NKP <- data.frame(equation   = "pcInd",
                        variant     = "pcIndAR",
                        sysMatrix   = c("exo"),
                        variableRow = c("pcIndl"),
                        varName     = c("pcPhi"),
                        lowerBound  = NA,
                        upperBound  = NA,
                        statRestr   = c("AR1"),
                        mean        = 0,
                        std         = 2,
                        distribution = "normal",
                        stringsAsFactors = FALSE)
pcInd <- rbind(pcInd[[1]], pcInd[[2]], pcInd[[3]], pcInd[[4]])

# kuttner
infl$base <- data.frame(equation   = "infl",
                         variant     = "base",
                         sysMatrix   = c("Z", "exo", "Z"),
                         variableRow = c("const", "gdpGL1", "cycle"),
                         varName     = c("inflConst","inflGdpGL1", "inflC0"),
                         # lowerBound  = c(-10, NA, -10),
                         # upperBound  = c(10, NA, 10),
                        lowerBound  = c(NA, NA, NA),
                        upperBound  = c(NA, NA, NA),
                         statRestr   = c(""),
                         mean        = c(0, 0, 0),
                         std         = c(5, 5, 2),
                         distribution = c("normal", "normal", "normal"),
                        stringsAsFactors = FALSE)
infl$cycleLag <- data.frame(equation    = "infl",
                             variant     = "cycleLag",
                             sysMatrix   = rep("Z", 10),
                             variableRow = paste0("cycleLag", 1:10),
                             varName     = paste0("inflC", 1:10),
                             lowerBound  = NA,
                             upperBound  = NA,
                             statRestr   = c(""),
                             mean        = rep(0, 10),
                             std         = rep(2, 10),
                             distribution = rep("normal", 10),
                            stringsAsFactors = FALSE)
infl <- rbind(infl[[1]], infl[[2]])

# combine
system <- rbind(cycle, trend, error, cubs, pcInd, infl)
dfSystem <- system

# save hidden data
# save(indicatorList, ameco, dfSystem, file = "R/sysdata.rda")
usethis::use_data(indicatorList, ameco, dfSystem, internal = TRUE, overwrite = TRUE)

