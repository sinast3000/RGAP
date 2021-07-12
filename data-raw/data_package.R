
# --------------------------------------------------------------------------------------------------- #

# clear
rm(list=ls())

library(lubridate)
library(zoo)
library(xts)

# --------------------------------------------------------------------------------------------------- #
# package data

# load hidden data
load(file = "R/sysdata.rda")

# ----- indicators for cubs
source <- indicatorList$source
country <- indicatorList$country

# load all files
for (filename in unique(source$file)) {
  sheet <- source$sheet[source$file==filename]
  for (sheetname in sheet) {
    assign(paste0(filename,sheetname), openxlsx::read.xlsx(paste0("data-raw/",filename),
                                                           sheet = sheetname,
                                                           startRow = 1))
  }
}
# loop and assign ts object
for (key in source$key) {
  filename  <- source$file[source$key==key]
  sheetname <- source$sheet[source$key==key]
  type      <- source$type[source$key==key]
  frequency <- source$freq[source$key==key]
  assign(paste0("key_country",key) ,sapply(country$abbreviation, function(x) gsub("XXX",x,key)))
  # get start time
  time <- get(paste0(filename,sheetname))[,1]
  if (frequency == 4) {
    time <- as.yearqtr(time,format = "%Y-Q%q")
    start_date <- c(year(first(time)), quarter(first(time)))}
  if (frequency == 12) {
    time <- as.yearmon(openxlsx::convertToDateTime(time))
    start_date <- c(year(first(time)), month(first(time)))}
  # ts object
  tmp <- ts(get(paste0(filename,sheetname))[,get(paste0("key_country",key))],
            start = start_date,
            frequency = frequency)
  colnames(tmp) <- country$name
  assign(paste0("indicator",type), tmp)
}
# reorganize data
indicator <- list()
for (c in country$name) {
  tmp <- list()
  for (type in source$type) {
    tmp[[tolower(type)]] <- get(paste0("indicator",type))[,c]
  }
  indicator[[c]] <- tmp
}
# temporal aggregation
library(tempdisagg)
indicatorTA <- list()
for (c in country$name) {
  tmp <- NULL
  mat <- NULL
  for (type in source$type) {
    tmp <- tempdisagg::ta(indicator[[c]][[tolower(type)]], convsersion = "average", to = "annual")
    tmp <- .cubsTa(tsObj = indicator[[c]][[tolower(type)]], conversion = "average", frequency = 1)
    mat <- cbind(mat, tmp)
  }
  colnames(mat) <- tolower(source$type)
  indicatorTA[[c]] <- mat
}


# save indicators
# devtools::use_data(indicator, indicator, overwrite = TRUE)
usethis::use_data(indicator, internal = FALSE, overwrite = TRUE)


# ----- ameco data
source <- ameco$source
country <- ameco$country

# load data set
filename <- "ameco_autumn2018.xlsx"
dfameco <- openxlsx::read.xlsx(paste0("data-raw/",filename),
                               sheet = "Sheet1",
                               startRow = 1)
# index and times
index <- !is.na(as.numeric(colnames(dfameco)))
years <- as.numeric(colnames(dfameco)[index])
start <- first(years)
countrylist <- list()
# loop and assign ts object
for (a in country$abbreviation) {
  keynumber <- paste0(a,source$keynumber)
  name <- as.character(country$name[country$abbreviation==a])
  dftmp <- dfameco
  dftmp <- dftmp[dftmp$CODE %in% keynumber,]
  if (dim(dftmp)[1]>0) {
    tstmp <- ts(apply(dftmp[,index],1, function(x) as.numeric(as.character(x))),start=start)
    index_colnames <- sapply(keynumber,match,as.character(dftmp$CODE))
    colnames(tstmp)[index_colnames[!is.na(index_colnames)]] <- as.character(source$finalname[!is.na(index_colnames)])
    countrylist[[name]] <- tstmp
  }
}

# delete insufficient data
countrylist[lapply(countrylist,function(x) dim(x)[2])<17] <- NULL
gap <- countrylist

# merge indicators for cubs
for (c in names(gap)) {
  namesTmp <- c(colnames(gap[[c]]), colnames(indicatorTA[[c]]))
  gap[[c]] <- cbind(gap[[c]], indicatorTA[[c]])
  colnames(gap[[c]]) <- namesTmp
}

# compute cubs
namesCubs <- c("indu","serv", "buil")
namesVACubs <- paste0("va", namesCubs)
for (x in names(gap)) {
  tryCatch({
    tscubs <- cubs(tsCU = gap[[x]][, namesCubs],
                   tsVA = gap[[x]][, namesVACubs])
    namesTmp <- c(colnames(gap[[x]]), "cubs")
    gap[[x]] <- cbind(gap[[x]], tscubs$cubs)
    colnames(gap[[x]]) <- namesTmp
  },
  error = function(cond) {  }
  )

}


# devtools::use_data(gap, gap, overwrite = TRUE)
usethis::use_data(gap, internal = FALSE, overwrite = TRUE)
