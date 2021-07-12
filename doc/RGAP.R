# clear
rm(list=ls())


####### Settings ------------------------------------------------------

# load package
# library(RGAP)
devtools::load_all()

# file path for plots
path <- "doc"
dir.create(path, recursive = TRUE)
# output format
library("knitr")
options(prompt = "R> ", continue = "+  ", width = 77, useFancyQuotes = FALSE)
opts_chunk$set(fig.align = 'center')
render_sweave()


####### Example: NAWRU model specification ----------------------------
data("gap")
tsList <- amecoData2input(gap[["France"]], alpha = 0.65)
exoType <- initializeExo(varNames = "ws")
exoType[1, , "difference"] <- 2
exoType[1, , "lag"] <- 0
model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2",
                    type = "TKP", cycleLag = 0:1, exoType = exoType)
model


####### Example: TFP model specification ------------------------------
data("gap")
tsList <- amecoData2input(gap[["Italy"]], alpha = 0.65)
model <- TFPmodel(tsl = tsList, trend = "DT", cycle = "RAR2",
                  cycleLag = 0, cubsErrorAR = 1)
model


####### Example: data -------------------------------------------------
# --- internal data ---------------------------------------------------
data("gap")
tsList <- amecoData2input(gap[["Germany"]], alpha = 0.65)
# --- cubs ------------------------------------------------------------
namesCubs <- c("indu","serv", "buil")
namesVACubs <- paste0("va", namesCubs)
tscubs <- cubs(tsCU = gap[["Germany"]][, namesCubs],
               tsVA = gap[["Germany"]][, namesVACubs])
tsList <- c(tsList, tscubs)
# --- AMECO data ------------------------------------------------------
tslBase <- fetchAmecoData(country = "Germany")
tsList <- amecoData2input(tslBase, alpha = 0.65)


####### Example: NAWRU model estimation -------------------------------
country <- "France"
dir.create(file.path(path, country), recursive = TRUE)
# --- data, model -----------------------------------------------------
data("gap")
tsList <- amecoData2input(gap[["France"]], alpha = 0.65)
exoType <- initializeExo(varNames = "ws")
exoType[1, , "difference"] <- 2
exoType[1, , "lag"] <- 0
model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2",
                    type = "TKP", cycleLag = 0, exoType = exoType)
# ---
exoType
# --- MLE -------------------------------------------------------------
parRestr <- initializeRestr(model = model, type = "hp")
parRestr
# ---
fit <- fitNAWRU(model = model, parRestr = parRestr)
plot(fit)
# ---
plot(fit, path = file.path(path, country), prefix = "mle")
plot(fit, path = file.path(path, country), prefix = "mle_wide7", width = 7)
# --- anchor ----------------------------------------------------------
fit <- trendAnchor(fit = fit, anchor = 8, h = 10, returnFit = TRUE)
plot(fit)
# ---
plot(fit, path = file.path(path, country) , prefix = "mle_anchor")
plot(fit, path = file.path(path, country) , prefix = "mle_wide7_anchor", width = 7)


####### Example: TFP model estimation ---------------------------------
country <- "Italy"
dir.create(file.path(path, country), recursive = TRUE)
# --- data, model, MLE ------------------------------------------------
data("gap")
tsList <- amecoData2input(gap[["Italy"]], alpha = 0.65)
model <- TFPmodel(tsl = tsList, trend = "RW2", cycle = "RAR2",
                  cycleLag = 0, cubsErrorARMA = c(0, 0))
parRestr <- initializeRestr(model = model, type = "hp")
fit <- fitTFP(model = model, parRestr = parRestr)
plot(fit)
# ---
plot(fit, path = file.path(path, country), prefix = "mle")
plot(fit, path = file.path(path, country), prefix = "mle_wide7", width = 7)
# --- Bayesian --------------------------------------------------------
prior <- initializePrior(model = model)
prior
# ---
fitBayes <- fitTFP(model = model, method = "bayesian", prior = prior,
                   R = 5000, thin = 2, MLEfit = fit)
plot(fitBayes)
plot(fitBayes, posterior = TRUE)
# ---
plot(fitBayes, posterior = TRUE, path = file.path(path, country), prefix = "bayes")
plot(fitBayes, posterior = FALSE, path = file.path(path, country), prefix = "bayes")


####### Example: Estimating the output gap ----------------------------
country <- "Netherlands"
dir.create(file.path(path, country), recursive = TRUE)
# --- data and models -------------------------------------------------
data("gap")
tsList <- amecoData2input(gap[["Netherlands"]], alpha = 0.65)
model <- parRestr <- prior <- fit <- list()
exoType <- initializeExo(varNames = c("ws", "prod", "tot"))
exoType[1, , "difference"] <- 2
exoType[2, , "difference"] <- 1
exoType[1, , "lag"] <- 0
exoType[2, , "lag"] <- 1
model$nawru <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2",
                          type = "TKP", cycleLag = 0, exoType = exoType)
model$tfp <- TFPmodel(tsl = tsList, trend = "DT", cycle = "RAR2",
                      cycleLag = 0, cubsErrorARMA = c(0,0))
# --- MLE -------------------------------------------------------------
parRestr$nawru <- initializeRestr(model = model$nawru, type = "hp")
fit$nawru <- fitNAWRU(model = model$nawru, parRestr = parRestr$nawru)
# ---
plot(fit$nawru, path = file.path(path, country), prefix = "mle")
# ---
parRestr$tfp <- initializeRestr(model = model$tfp, type = "hp")
fit$tfp <- fitTFP(model = model$tfp, parRestr = parRestr$tfp)
# ---
plot(fit$tfp, path = file.path(path, country), prefix = "mle")
# ---
fit$gap <- gapProd(tsl = tsList, NAWRUfit = fit$nawru,
                   TFPfit = fit$tfp, lambda = 100, alpha = 0.65)
plot(fit$gap)
# ---
plot(fit$gap, path = file.path(path, country), prefix = "mle")
# ---
plot(fit$gap, contribution = TRUE)
# ---
plot(fit$gap, contribution = TRUE, path = file.path(path, country), prefix = "mle_contr")


####### Example: alternative models -----------------------------------
country <- "Netherlands"
dir.create(file.path(path, country), recursive = TRUE)
# --- Kuttner model ---------------------------------------------------
data("gap")
tsList <- as.list(gap[["Netherlands"]][,c("cpih","gdp")])
tsList$infl <- diff(tsList$cpih)
model <- KuttnerModel(tsl = tsList, trend = "RW2",
                      cycleLag = 1, cycle = "AR2")
parRestr <- initializeRestr(model = model, type = "hp")
gapKuttner <- fitKuttner(model, parRestr)
# ---
plot(gapKuttner, path = file.path(path, country), prefix = "kuttner")
# --- HP filter -------------------------------------------------------
gapHPfilter <- gapHP(tsList$gdp, lambda = 10)
# ---
plot(gapHPfilter, path = file.path(path, country), prefix = "hp")
# --- plot gap comparison ---------------------------------------------
{
  tsl_plot <- list("Kuttner" = gapKuttner$tsl$gap, "HP" = gapHPfilter$gap, "EC" = fit$gap$tsl$gap)

  # settings
  set <- list()
  set$title <- "Output gap (in %)"
  set$legend <- names(tsl_plot)
  set$titleFontsize <- 10
  set$labelFontsize <- 8
  set$legendFontsize <- 8
  set$colors <- c("black", "grey12", "grey24", "grey36", "grey48", "grey60")
  set$linetype <- c(2,1,3,4,5,6)
  set$freqYear <- 1 + floor(length(tsl_plot[[1]])/frequency(tsl_plot[[1]])/15)
  # data
  date <- zoo::as.Date(time(tsl_plot[[1]]))
  tsm <- do.call(cbind, tsl_plot)
  tsm <- na.trim(tsm)
  df <- data.frame(date = zoo::as.Date(time(tsm)), tsm)
  # plot
  p0 <- ggplot(df, aes(x = date)) +
    theme_classic() + coord_cartesian(xlim = c(date[1], date[length(date)]), expand = FALSE)  +
    theme(plot.title = element_text(size = set$titleFontsize),
          axis.title = element_text(size = set$labelFontsize),
          panel.grid.major.y = element_line(linetype = "solid")) +
    geom_line(aes(y = Kuttner, col = set$legend[1], linetype = set$legend[1])) +
    geom_line(aes(y = HP, col = set$legend[2], linetype = set$legend[2])) +
    geom_line(aes(y = EC, col = set$legend[3], linetype = set$legend[3])) +
    scale_color_manual(name = NULL, values = rev(set$colors[1:(ncol(df)-1)])) +
    scale_linetype_manual(name = NULL, values = set$linetype) +
    scale_x_date(date_labels = "%Y", date_minor_breaks = "1 year",
                 date_breaks = paste0(set$freqYear, " year")) +
    labs(title = set$title, x = "year", y = "") +
    theme(legend.position = c(1, 0.01),
          legend.justification = c(1, 0),
          legend.margin = margin(0, 0, 0, 0),
          legend.text=element_text(size = set$legendFontsize),
          legend.key.size = unit(0.75, 'lines'),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.spacing.x = unit(0, "pt"),
          legend.spacing.y = unit(0, "pt")) +
    guides(col = guide_legend(ncol = 2, byrow = TRUE))
  p0
  ggsave(filename = file.path(path, country, "gap_comparison.png"), plot = p0, width = 7, height = 3)
}


####### Example: Multiple countries -----------------------------------
library(reshape2)
data("gap")
dir.create(file.path(path, "all_mle"), recursive = TRUE)
## --------------------------------------------------------------------
countries <- c("Belgium", "Denmark", "Germany", "Ireland", "Greece",
               "Spain", "France", "Italy", "Luxembourg", "Netherlands",
               "Austria", "Poland", "Portugal", "Finland", "Sweden",
               "United Kingdom")
tsl <- model <- parRestr <- fit <- list()
for (k in countries) {
  # data
  tsl[[k]] <- amecoData2input(gap[[k]], alpha = 0.65)
  # NAWRU
  exoType <- initializeExo(varNames = c("ws", "prod", "tot"))
  exoType[1, , "difference"] <- 2
  exoType[2, , "difference"] <- 1
  exoType[1, , "lag"] <- 0
  exoType[2, , "lag"] <- 1
  model[[k]]$nawru <- NAWRUmodel(tsl = tsl[[k]], trend = "RW2", cycle = "AR2",
                                 type = "TKP", cycleLag = 0, exoType = exoType)
  parRestr[[k]]$nawru <- initializeRestr(model = model[[k]]$nawru, type = "hp")
  fit[[k]]$nawru <- fitNAWRU(model = model[[k]]$nawru, parRestr = parRestr[[k]]$nawru)
  plot(fit[[k]]$nawru, path = file.path(path,"all_mle"), prefix = gsub(" ","_", k))
  # TFP trend
  model[[k]]$tfp <- TFPmodel(tsl = tsl[[k]], trend = "DT", cycle = "RAR2",
                             cycleLag = 0, cubsErrorAR = 0)
  parRestr[[k]]$tfp <- initializeRestr(model = model[[k]]$tfp, type = "hp")
  fit[[k]]$tfp <- fitTFP(model = model[[k]]$tfp, parRestr = parRestr[[k]]$tfp)
  plot(fit[[k]]$tfp, path = file.path(path,"all_mle"), prefix = gsub(" ","_", k))
  # gap
  fit[[k]]$gap <- gapProd(tsl = tsl[[k]], NAWRUfit = fit[[k]]$nawru,
                          TFPfit = fit[[k]]$tfp, lambda = 100, alpha = 0.65)
  plot(fit[[k]]$gap, path = file.path(path,"all_mle"), prefix = gsub(" ","_", k))
}
## ---- save ---------------------------------------------------------
save(tsl, model, parRestr, fit, file = file.path(path, "all_mle", "results.RData"))
## ---- plot ---------------------------------------------------------
{
  tsl_plot <- lapply(lapply(lapply(fit, "[[", "gap"), "[[", "tsl"), "[[", "gap")

  # settings
  set <- list()
  set$title <- "Output gap (in %)"
  set$legend <- names(tsl_plot)
  set$titleFontsize <- 10
  set$labelFontsize <- 8
  set$legendFontsize <- 8
  set$color <-  rep(c("darkolivegreen4", "darkorange2", "deepskyblue3", "darkred", "darkgrey", "black"), each = 3)
  set$linetype <- rep(c(1,2,3),6)
  set$freqYear <- 1 + floor(length(tsl[[1]])/frequency(tsl[[1]])/15)
  # data
  date <- zoo::as.Date(time(tsl_plot[[1]]))
  tsm <- do.call(cbind, tsl_plot)
  df <- data.frame(date = zoo::as.Date(time(tsm)), tsm)
  df_melt <- melt(df, id.vars = "date")
  df_melt$variable <- as.factor(gsub("\\.", " ", as.character(df_melt$variable )))
  # plot
  p0 <- ggplot(df_melt, aes(x = date, y = value, color = variable, linetype = variable)) + geom_line() +
    theme_classic() + coord_cartesian(xlim = c(date[1], date[length(date)]), expand = FALSE)  +
    theme(plot.title = element_text(size = set$titleFontsize),
          axis.title = element_text(size = set$labelFontsize),
          panel.grid.major.y = element_line(linetype = "solid")) +
    scale_color_manual(name = NULL, values = set$color[1:length(countries)]) +
    scale_linetype_manual(name = NULL, values = set$linetype) +
    scale_x_date(date_labels = "%Y", date_minor_breaks = "1 year",
                 date_breaks = paste0(set$freqYear, " year")) +
    labs(title = set$title, x = "year", y = "") +
    theme(legend.position = c(0, 0.01),
          legend.justification = c(-0.1, 0),
          legend.margin = margin(0, 0, 0, 0),
          legend.text=element_text(size = set$legendFontsize),
          legend.key.size = unit(0.75, 'lines'),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.spacing.x = unit(0, "pt"),
          legend.spacing.y = unit(0, "pt")) +
    guides(col = guide_legend(ncol = 3, byrow = TRUE))
  p0
  ggsave(filename = file.path(path, "all_mle", "gap_comparison_countries.png"), plot = p0, width = 7, height = 3)
}


####### Example: GAP program by EC -----------------------------------------
country <- "France"
dir.create(file.path(path, "EC_gap"), recursive = TRUE)
# --- data, model, fit -----------------------------------------------------
data("gap")
tsList <- amecoData2input(gap[["France"]], alpha = 0.65)
exoType <- initializeExo(varNames = "ws")
exoType[1, , "difference"] <- 2
exoType[1, , "lag"] <- 0
model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2",
                    type = "TKP", cycleLag = 0, exoType = exoType)
parRestr <- initializeRestr(model = model, type = "hp")
fit <- fitNAWRU(model = model, parRestr = parRestr)
plot(fit)
write.csv(do.call(cbind,model$tsl), file.path(path, "EC_gap", "data_gap_ec.csv"))
# --- EC gap result --------------------------------------------------------
resEC <- read.csv(file.path(path, "EC_gap", "data_gap_ec_result.csv"))
nawruEC <- ts(resEC$nawru, start = start(model$tsl$ur))
# --- plot -----------------------------------------------------------------
{
  tsl_plot <- list("RGAP" = fit$tsl$nawru, "ECGAP" = nawruEC)

  # settings
  set <- list()
  set$title <- "NAWRU"
  set$legend <- c("RGAP","EC GAP Version 5.0")
  set$titleFontsize <- 10
  set$labelFontsize <- 8
  set$legendFontsize <- 8
  set$colors <- c("black", "grey12", "grey24", "grey36", "grey48", "grey60")
  set$linetype <- c(2,1,3,4,5,6)
  set$freqYear <- 1 + floor(length(tsl_plot[[1]])/frequency(tsl_plot[[1]])/15)
  # data
  date <- zoo::as.Date(time(tsl_plot[[1]]))
  tsm <- do.call(cbind, tsl_plot)
  tsm <- na.trim(tsm)
  df <- data.frame(date = zoo::as.Date(time(tsm)), tsm)
  # plot
  p0 <- ggplot(df, aes(x = date)) +
    theme_classic() + coord_cartesian(xlim = c(date[1], date[length(date)]), expand = FALSE)  +
    theme(plot.title = element_text(size = set$titleFontsize),
          axis.title = element_text(size = set$labelFontsize),
          panel.grid.major.y = element_line(linetype = "solid")) +
    geom_line(aes(y = ECGAP, col = set$legend[1], linetype = set$legend[1])) +
    geom_line(aes(y = RGAP, col = set$legend[2], linetype = set$legend[2])) +
    scale_color_manual(name = NULL, values = rev(set$colors[1:(ncol(df)-1)])) +
    scale_linetype_manual(name = NULL, values = set$linetype) +
    scale_x_date(date_labels = "%Y", date_minor_breaks = "1 year",
                 date_breaks = paste0(set$freqYear, " year")) +
    labs(title = set$title, x = "year", y = "") +
    theme(legend.position = c(1, 0.01),
          legend.justification = c(1, 0),
          legend.margin = margin(0, 0, 0, 0),
          legend.text=element_text(size = set$legendFontsize),
          legend.key.size = unit(0.75, 'lines'),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.spacing.x = unit(0, "pt"),
          legend.spacing.y = unit(0, "pt")) +
    guides(col = guide_legend(ncol = 2, byrow = TRUE))
  p0
  ggsave(filename = file.path(path, "EC_gap", "rgap_ecgap_comparison.png"), plot = p0, width = 7, height = 3)
}

