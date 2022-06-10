# clear
rm(list=ls())


####### Settings ------------------------------------------------------

# load package
# library(RGAP)
library(stats)
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
exoType <- initializeExo(varNames = "ws", D = 2, L = 0)
model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2",
                    type = "TKP", cycleLag = 0:1, exoType = exoType)
model


####### Example: TFP model specification ------------------------------
data("gap")
tsList <- amecoData2input(gap[["Italy"]], alpha = 0.65)
model <- TFPmodel(tsl = tsList, trend = "DT", cycle = "RAR2",
                  cycleLag = 0, cubsErrorARMA = c(1, 0))
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
exoType <- initializeExo(varNames = "ws", D = 2, L = 0)
model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2",
                    type = "TKP", cycleLag = 0, exoType = exoType)
# ---
exoType
# --- MLE -------------------------------------------------------------
parRestr <- initializeRestr(model = model, type = "hp")
parRestr
# ---
f <- fit(model = model, parRestr = parRestr)
plot(f)
# ---
plot(f, path = file.path(path, country), prefix = "mle")
plot(f, path = file.path(path, country), prefix = "mle_wide7", width = 7)
# --- anchor ----------------------------------------------------------
f <- trendAnchor(fit = f, anchor = 8, h = 10, returnFit = TRUE)
plot(f)
# ---
plot(f, path = file.path(path, country) , prefix = "mle_anchor")
plot(f, path = file.path(path, country) , prefix = "mle_wide7_anchor", width = 7)


####### Example: TFP model estimation ---------------------------------
country <- "Italy"
dir.create(file.path(path, country), recursive = TRUE)
# --- data, model, MLE ------------------------------------------------
data("gap")
tsList <- amecoData2input(gap[["Italy"]], alpha = 0.65)
model <- TFPmodel(tsl = tsList, trend = "DT", cycle = "RAR2",
                  cycleLag = 0, cubsErrorARMA = c(0, 0))
parRestr <- initializeRestr(model = model, type = "hp")
f <- fit(model = model, parRestr = parRestr)
plot(f)
# ---
plot(f, path = file.path(path, country), prefix = "mle")
# --- Prediction ------------------------------------------------------
fPred <- predict(object = f, n.ahead = 10)
plot(fPred, alpha = 0.1)
# ---
plot(fPred, alpha = 0.1, combine = FALSE, path = file.path(path, country), 
     prefix = "prediction_wide5", width = 5)
# --- Bayesian --------------------------------------------------------
prior <- initializePrior(model = model)
prior
# ---
fBayes <- fit(model = model, method = "bayesian", prior = prior,
               R = 5000, thin = 2, MLEfit = f)
plot(fBayes)
plot(fBayes, posterior = TRUE)
# ---
plot(fBayes, posterior = TRUE, path = file.path(path, country), prefix = "bayes")
plot(fBayes, posterior = FALSE, path = file.path(path, country), prefix = "bayes")



####### Example: Estimating the output gap ----------------------------
country <- "Netherlands"
dir.create(file.path(path, country), recursive = TRUE)
# --- data and models -------------------------------------------------
data("gap")
tsList <- amecoData2input(gap[["Netherlands"]], alpha = 0.65)
model <- parRestr <- prior <- fits <- list()
D <- matrix(c(2, 2, 2, 1, 1, 1), 2, 3, byrow = TRUE)
L <- matrix(c(0, 0, 0, 1, 1, 1), 2, 3, byrow = TRUE)
exoType <- initializeExo(varNames = c("ws", "prod","tot"), D = D, L = L)
model$nawru <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2",
                          type = "TKP", cycleLag = 0, exoType = exoType)
model$tfp <- TFPmodel(tsl = tsList, trend = "DT", cycle = "RAR2",
                      cycleLag = 0, cubsErrorARMA = c(0,0))
# --- MLE -------------------------------------------------------------
parRestr$nawru <- initializeRestr(model = model$nawru, type = "hp")
fits$nawru <- fit(model = model$nawru, parRestr = parRestr$nawru)
# ---
plot(fits$nawru, path = file.path(path, country), prefix = "mle")
# ---
parRestr$tfp <- initializeRestr(model = model$tfp, type = "hp")
fits$tfp <- fit(model = model$tfp, parRestr = parRestr$tfp)
# ---
plot(fits$tfp, path = file.path(path, country), prefix = "mle")
# ---
fits$gap <- gapProd(tsl = tsList, NAWRUfit = fits$nawru,
                   TFPfit = fits$tfp, lambda = 100, alpha = 0.65)
plot(fits$gap)
# ---
plot(fits$gap, path = file.path(path, country), prefix = "mle")
# ---
plot(fits$gap, contribution = TRUE)
# ---
plot(fits$gap, contribution = TRUE, path = file.path(path, country), prefix = "mle_contr")


####### Example: alternative models -----------------------------------
country <- "Netherlands"
dir.create(file.path(path, country), recursive = TRUE)
# --- Kuttner model ---------------------------------------------------
data("gap")
tsList <- as.list(gap[["Netherlands"]][,c("cpih","gdp")])
tsList$infl <- diff(tsList$cpih)
model <- KuttnerModel(tsl = tsList, trend = "RW2",
                      cycleLag = 1, cycle = "AR2")
parRestr <- initializeRestr(model = model, type = "hp", q = 0.1)
gapKuttner <- fit(model = model, parRestr = parRestr)
plot(gapKuttner)
# ---
plot(gapKuttner, path = file.path(path, country), prefix = "kuttner")
# --- HP filter -------------------------------------------------------
gapHPfilter <- gapHP(tsList$gdp, lambda = 100)
# ---
plot(gapHPfilter, path = file.path(path, country), prefix = "hp")


# --- plot gap comparison ---------------------------------------------
{
  tsl_plot <- list("Kuttner" = gapKuttner$tsl$gap, "HP" = gapHPfilter$gap, "EC" = fits$gap$tsl$gap)

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
    theme(
      legend.position="bottom",
      legend.box = "horizontal",
      legend.box.margin=margin(-15,0,0,0),
      legend.justification="left",
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
    ) +
    guides(color = guide_legend(order=1, ncol = 6, byrow = TRUE),
           linetype = guide_legend(order=1, ncol = 6, byrow = TRUE),
           fill = guide_legend(order=2, ncol = 6)) 
  
  p0
  ggsave(filename = file.path(path, country, "gap_comparison.png"), plot = p0, width = 7, height = 3)
}

####### Example: Multiple countries -----------------------------------
library(reshape2)
data("gap")
type <- "hp"
path_sub <- paste0("all_mle_", type)
dir.create(file.path(path, path_sub), recursive = TRUE)
## --------------------------------------------------------------------
countries <- c("Belgium", "Denmark", "Germany", "Ireland", "Greece",
               "Spain", "France", "Italy", "Luxembourg", "Netherlands",
               "Austria", "Poland", "Portugal", "Finland", "Sweden",
               "United Kingdom")
tsl <- model <- parRestr <- fits <- list()
for (k in countries) {
  # data
  tsl[[k]] <- amecoData2input(gap[[k]], alpha = 0.65)
  # NAWRU
  D <- matrix(c(2, 2, 2, 1, 1, 1), 2, 3, byrow = TRUE)
  L <- matrix(c(0, 0, 0, 1, 1, 1), 2, 3, byrow = TRUE)
  exoType <- initializeExo(varNames = c("ws", "prod","tot"), D = D, L = L)
  model[[k]]$nawru <- NAWRUmodel(tsl = tsl[[k]], trend = "RW2", cycle = "AR2",
                                 type = "TKP", cycleLag = 0, exoType = exoType)
  parRestr[[k]]$nawru <- initializeRestr(model = model[[k]]$nawru, type = type)
  fits[[k]]$nawru <- fit(model = model[[k]]$nawru, parRestr = parRestr[[k]]$nawru)
  plot(fits[[k]]$nawru, path = file.path(path, path_sub), prefix = gsub(" ","_", k))
  # TFP trend
  model[[k]]$tfp <- TFPmodel(tsl = tsl[[k]], trend = "DT", cycle = "RAR2",
                             cycleLag = 0, cubsErrorAR = 0)
  parRestr[[k]]$tfp <- initializeRestr(model = model[[k]]$tfp, type = type)
  fits[[k]]$tfp <- fit(model = model[[k]]$tfp, parRestr = parRestr[[k]]$tfp)
  plot(fits[[k]]$tfp, path = file.path(path, path_sub), prefix = gsub(" ","_", k))
  # gap
  fits[[k]]$gap <- gapProd(tsl = tsl[[k]], NAWRUfit = fits[[k]]$nawru,
                          TFPfit = fits[[k]]$tfp, lambda = 100, alpha = 0.65)
  plot(fits[[k]]$gap, path = file.path(path, path_sub), prefix = gsub(" ","_", k))
}
## ---- save ---------------------------------------------------------
save(tsl, model, parRestr, fits, file = file.path(path, path_sub, "results.RData"))
## ---- plot ---------------------------------------------------------
{
  tsl_plot <- lapply(lapply(lapply(fits, "[[", "gap"), "[[", "tsl"), "[[", "gap")

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
    theme(
      legend.position="bottom",
      legend.box = "horizontal",
      legend.box.margin=margin(-15,0,0,0),
      legend.justification="left",
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
    ) +
    guides(color = guide_legend(order=1, ncol = 6, byrow = TRUE),
           linetype = guide_legend(order=1, ncol = 6, byrow = TRUE),
           fill = guide_legend(order=2, ncol = 6)) 
  
  p0
  ggsave(filename = file.path(path, path_sub, "gap_comparison_countries.png"), plot = p0, width = 7, height = 3)
}



####### Example: Multiple countries model selection -------------------
library(reshape2)
data("gap")
type <- "hp"
path_sub <- paste0("all_mle_model_selection_", type)
dir.create(file.path(path, path_sub), recursive = TRUE)
## --------------------------------------------------------------------
# exclude Poland, Spain, Ireland, Austria
countries <- c("Belgium", "Denmark", "Finland", "France", 
               "Germany", "Greece", "Italy", "Luxembourg", 
               "Netherlands", "Portugal", "Sweden", "United Kingdom")
tsl <- fits <- list()
for (k in countries) {
  # data
  tsl[[k]] <- amecoData2input(gap[[k]])
  # gap
  fits[[k]] <- autoGapProd(tsl = tsl[[k]], type = type, fast = TRUE, nModels = 5, q = 0.1)

  plot(fits[[k]]$nawru$fit[[1]], path = file.path(path,path_sub), prefix = gsub(" ","_", k))
  plot(fits[[k]]$tfp$fit[[1]], path = file.path(path,path_sub), prefix = gsub(" ","_", k))
  plot(fits[[k]]$gap, path = file.path(path,path_sub), prefix = gsub(" ","_", k))
}
## ---- save ---------------------------------------------------------
save(tsl, fits, file = file.path(path, path_sub, "results.RData"))
## ---- plot ---------------------------------------------------------
{
  tsl_plot <- lapply(lapply(lapply(fits, "[[", "gap"), "[[", "tsl"), "[[", "gap")
  
  # settings
  set <- list()
  set$title <- "Output gap (in %)"
  set$legend <- names(tsl_plot)
  set$titleFontsize <- 10
  set$labelFontsize <- 8
  set$legendFontsize <- 8
  set$color <-  rep(c("darkolivegreen4", "darkorange2", "deepskyblue3", "darkred", "darkgrey", "black"), each = 2)
  set$linetype <- rep(c(1,2),6)
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
    theme(
      legend.position="bottom",
      legend.box = "horizontal",
      legend.box.margin=margin(-15,0,0,0),
      legend.justification="left",
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
    ) +
    guides(color = guide_legend(order=1, ncol = 6, byrow = TRUE),
           linetype = guide_legend(order=1, ncol = 6, byrow = TRUE),
           fill = guide_legend(order=2, ncol = 6)) 
  p0
  ggsave(filename = file.path(path, path_sub, "gap_comparison_countries.png"), plot = p0, width = 7, height = 3)
}


####### Example: GAP program by EC -----------------------------------------
country <- "France"
dir.create(file.path(path, "EC_gap"), recursive = TRUE)
# --- data, model, fit -----------------------------------------------------
data("gap")
tsList <- amecoData2input(gap[["France"]], alpha = 0.65)
exoType <- initializeExo(varNames = "ws", D = 2, L = 0)
model <- NAWRUmodel(tsl = tsList, trend = "RW2", cycle = "AR2",
                    type = "TKP", cycleLag = 0, exoType = exoType)
parRestr <- initializeRestr(model = model, type = "hp")
parRestr$pcInd[,"pcddws"] <- c(-10, 10)
parRestr$pcInd[,"pcConst"] <- c(-0.1, 0.1)
parRestr$pcInd[,"pcC0"] <- c(-10, 0)
fits <- fit(model = model, parRestr = parRestr)
plot(fits)
write.csv(do.call(cbind,model$tsl), file.path(path, "EC_gap", "data_gap_ec.csv"))
# --- EC gap result --------------------------------------------------------
resEC <- read.csv(file.path(path, "EC_gap", "data_gap_ec_result.csv"))
nawruEC <- ts(resEC$nawru, start = start(model$tsl$ur))
# --- plot -----------------------------------------------------------------
{
  tsl_plot <- list("RGAP" = fits$tsl$nawru, "ECGAP" = nawruEC)

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
    theme(
      legend.position="bottom",
      legend.box = "horizontal",
      legend.box.margin=margin(-15,0,0,0),
      legend.justification="left",
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
    ) +
    guides(color = guide_legend(order=1, ncol = 6, byrow = TRUE),
           linetype = guide_legend(order=1, ncol = 6, byrow = TRUE),
           fill = guide_legend(order=2, ncol = 6)) 
  p0
  ggsave(filename = file.path(path, "EC_gap", "rgap_ecgap_comparison.png"), plot = p0, width = 7, height = 3)
}

