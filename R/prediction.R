
# -------------------------------------------------------------------------------------------

#' Predictions
#'
#' @description Computes predictions for an object of class \code{NAWRUfit, TFPfit}, or 
#' \code{KuttnerFit} estimated via MLE or Bayesian methods (objects of class \code{fit}).
#'
#' @param object An object of class \code{NAWRUfit}, \code{TFPfit}, or \code{KuttnerFit} 
#'   (objects of class \code{fit}).
#' @param n.ahead An integer specifying the prediction horizon. 
#' @param exogenous A character string specifying the computation of exogenous variables 
#'   included in the model (if applicable). Valid options are \code{exogenous = "mean"} and 
#'   \code{exogenous = "last".}
#' @param returnFit A logical. If \code{TRUE}, an object of the same class as \code{fit}
#'   where the list entry \code{tsl} is replaced. If \code{FALSE}, only the new time series 
#'   list is returned.
#' @param ... Ignored.
#'
#' @export
#' @importFrom stats predict
#' @aliases predict.fit
#' @return The fitted object with an updated time series list \code{tsl}. If 
#'   \code{returnFit = FALSE}, only the updated time series list is returned.
predict.fit <- function(object, n.ahead = 10, exogenous = "mean", returnFit = TRUE, ...) {

  method <- attr(object, "method")
  if (method != "MLE") {
    predictBayes(fit = object, n.ahead = n.ahead, exogenous = exogenous, returnFit = returnFit)
  } else {
    predictMLE(fit = object, n.ahead = n.ahead, exogenous = exogenous, returnFit = returnFit)
  }
  
}


# -------------------------------------------------------------------------------------------

#' Predictions for MLE
#'
#' @description Computes predictions for an object of class \code{NAWRUfit, TFPfit, KuttnerFit}
#' estimated via MLE.
#'
#' @inheritParams predict.fit
#' 
#' @keywords internal
predictMLE <- function(fit, n.ahead = 10, exogenous = "mean", returnFit = TRUE) {

  # adapt model
  model <- fit$SSMfit$model
  model$y <- window(model$y, end = end(model$y) + c(0, n.ahead), extend = TRUE)
  if (dim(model$Z)[3] > 1) {
    Z <- model$Z
    if (exogenous == "last") {
      Z <- Z[, , dim(Z)[3]]
    } else if (exogenous == "mean") {
      Z <- apply(Z, 1:2, mean)
    }
    model$Z <- array(c(c(model$Z), rep(Z, n.ahead)), dim(model$Z) + c(0, 0, n.ahead))
  }
  n_fc <- attr(model, "n") + as.integer(n.ahead)
  attr(model, "n") <- n_fc
  
  # filter and smoother
  out_fc <- KFS(model, simplify = FALSE, filtering = c("state", "signal"), smoothing = c("state", "signal", "disturbance"))
  
  # get observation equation and add forecast
  y <- fit$model$SSModel$y
  y <- window(y, end = end(y) + c(0, n.ahead), extend = TRUE)
  y[is.na(y)] <- out_fc$m[is.na(y)] # fill NAs with predictions
  out_fc$model$y <- y
  
  # get results
  tsl_fc <- .SSresults(out = out_fc, model = fit$model, prediction = TRUE)
  
  # add standard error for forecast of observation equation
  V_mu <- out_fc$V_mu
  tsl_fc$obsSE <- ts(rbind(
      matrix(0, n_fc - n.ahead, NCOL(y)),
      t(sqrt(apply(out_fc$V_mu[, , (n_fc - n.ahead + 1):n_fc], 3, diag)))
    ), start = start(y), frequency = frequency(y))

  # return
  if (returnFit) {
    fit$tsl <- tsl_fc
    attr(fit, "prediction")$n.ahead <- n.ahead
    attr(fit, "prediction")$exogenous <- exogenous
    return(fit)
  } else {
    attr(tsl_fc, "model") <- class(fit)
    return(tsl_fc)
  }

}


# -------------------------------------------------------------------------------------------

#' Predictions for Bayesian estimation
#'
#' @description Computes predictions for an object of class \code{NAWRUfit, TFPfit}
#' estimated via Bayesian methods.
#'
#' @inheritParams predict.fit
#' 
#' @keywords internal
predictBayes <- function(fit, n.ahead = 10, exogenous = "mean", returnFit = TRUE) {
  
  # Bayesian estimation attributes
  R <- attr(fit, "R")
  thin <-  attr(fit, "thin")
  burnin <-  attr(fit, "burnin")
  HPDIprob <-  attr(fit, "HPDIprob")
  FUN <- attr(fit, "FUN")
  
  # attributes
  model <- fit$model
  E2ind <- ifelse(inherits(model, "NAWRUmodel"), 
                  "phillips curve", "cubs")
  trend <- attr(model, "trend")
  cycle <- attr(model, "cycle")
  cycleLag <- attr(model, E2ind)$cycleLag
  cubsAR <- attr(model, E2ind)$cubsAR
  errorARMA <- attr(model, E2ind)$errorARMA
  loc <- model$loc
  
  # observables
  y <- fit$model$SSModel$y
  y <- window(y, end = end(y) + c(0, n.ahead), extend = TRUE)

  # state and parameters
  index <- (burnin / thin + 1):(R / thin)
  state <- fit$mcmc$states[, , index]
  state_f <- fit$mcmc$states_f[, , index]
  param <- fit$mcmc$parameters[index, ]
  
  # initialize and loop
  state_fc <- array(NA, dim(state) + c(n.ahead, 0, 0))
  state_f_fc <- array(NA, dim(state_f) + c(n.ahead, 0, 0))
  colnames(state_fc) <- colnames(state_f_fc) <- colnames(fit$SSMout$alphahat)
  y_fitted_fc <- array(NA, c(dim(y), dim(state_fc)[3]))
  colnames(y_fitted_fc) <- colnames(fit$model$SSModel$y)
  for (k in 1:dim(state)[3]) {
  
    # update model matrices
    SSModel <- .updateSSSystem(
      pars = param[k ,], SSModel = model$SSModel, loc = loc, cycle = cycle,
      trend = trend, errorARMA = errorARMA, bayes = TRUE
    )
    n <- attr(fit$SSMfit$model, "n")
    TT <- SSModel$T[, , 1]
    Z <- SSModel$Z
    if (exogenous == "last") {
      Z <- Z[, , n]
    } else if (exogenous == "mean") {
      Z <- apply(Z, 1:2, mean)
    }
    
    # smoothed and filtered state
    start <- start(model$SSModel$y)
    freq <- frequency(model$SSModel$y)
    alphahat <- ts(state[, , k], start = start, frequency = freq)
    att <- ts(state_f[, , k], start = start, frequency = freq)
    colnames(alphahat) <- colnames(att) <- colnames(fit$SSMout$alphahat)
    alphahat <- window(alphahat, end = end(alphahat) + c(0, n.ahead), extend = TRUE)
    att <- window(att, end = end(att) + c(0, n.ahead), extend = TRUE)
    
    # predict
    for (tt in 1:n.ahead) {
      # smoothed and filtered state
      alphahat[n + tt, ] <- att[n + tt, ] <- TT %*% alphahat[n + tt - 1, ]
      # alphahat[n, ] equals att[n, ]
    }
    state_fc[, , k] <- alphahat
    state_f_fc[, , k] <- att
    
    # predict observations
    y_fitted <- ts(alphahat %*% t(Z), start = start(alphahat), frequency = frequency(alphahat))
    y_fitted_fc[, , k] <- y_fitted
    
  }
  
  
  # ----- estimated states
  tslRes <- .SSresultsBayesian(model = model$SSModel, HPDIprob = HPDIprob, FUN = FUN,
                               state = state_fc, state_f = state_f_fc, obsFitted = y_fitted_fc)
  tslRes$obsSummary <- tslRes$obsFittedSummary
  tslRes$obsSummary[t(apply(t(!is.na(y)), 2, rep, each = NCOL(tslRes$obsSummary)/2))] <- NA
  
  index_na <- is.na(y)
  y[index_na] <- tslRes$obsFitted[index_na] # fill NAs with predictions
  tslRes$obs <- y
  
  # fill ts list depending on different objects
  if (inherits(fit, "NAWRUfit")) {
    
    # nawru
    tsTmp <- state_fc[, "trend", ]
    tslRes$nawruSummary <- ts(mcmcSummary(x = t(tsTmp), HPDIprob = HPDIprob), start = start, frequency = freq)
    tslRes$nawru <- tslRes$nawruSummary[, firstLetterUp(FUN)]
    # fitted pcInd equation
    tsTmp <- y_fitted_fc[, 2, ]
    tslRes$pcIndFittedSummary <- ts(mcmcSummary(x = t(tsTmp), HPDIprob = HPDIprob), start = start, frequency = freq)
    tslRes$pcIndFitted <- tslRes$pcIndFittedSummary[, firstLetterUp(FUN)]
    tslRes$pcInd[index_na[, 1]] <- exp(y)[, 1][index_na[, 1]]
    tslRes$pcIndSummary[!index_na[, 1]] <- NA
    
  } else if (inherits(fit, "TFPfit")) {
    # tfp
    tsTmp <- exp(y_fitted_fc[, 1, ])
    tslRes$tfpSummary <- ts(mcmcSummary(x = t(tsTmp), HPDIprob = HPDIprob), start = start, frequency = freq)
    tslRes$tfp <- tslRes$tfpSummary[, firstLetterUp(FUN)]  
    tslRes$tfp[index_na[, 1]] <- exp(y)[, 1][index_na[, 1]]
    tslRes$tfpSummary[!index_na[, 1]] <- NA
    # tfp growth rate
    tsTmp <- apply(exp(y_fitted_fc[, 1, ]), 2, function(x) (x[2:length(x)] / x[1:(length(x) - 1)]) - 1)
    tslRes$tfpGrowthSummary <- ts(mcmcSummary(x = t(tsTmp), HPDIprob = HPDIprob), start = start + c(0, 1), frequency = freq)
    tslRes$tfpGrowth <- tslRes$tfpGrowthSummary[, firstLetterUp(FUN)]
    tslRes$tfpGrowth[index_na[2:(NROW(index_na)), 1]] <- growth(exp(y[, 1]))[index_na[2:(NROW(index_na)), 1]]
    tslRes$tfpGrowthSummary[!index_na[2:(NROW(index_na)), 1]] <- NA
    # tfp trend
    tsTmp <- exp(state_fc[, "trend", ])
    tslRes$tfpTrendSummary <- ts(mcmcSummary(x = t(tsTmp), HPDIprob = HPDIprob), start = start, frequency = freq)
    tslRes$tfpTrend <- tslRes$tfpTrendSummary[, firstLetterUp(FUN)]
    # tfp trend growth rate
    tsTmp <- apply(exp(state_fc[, "trend", ]), 2, function(x) (x[2:length(x)] / x[1:(length(x) - 1)]) - 1)
    tslRes$tfpTrendGrowthSummary <- ts(mcmcSummary(x = t(tsTmp), HPDIprob = HPDIprob), start = start + c(0, 1), frequency = freq)
    tslRes$tfpTrendGrowth <- tslRes$tfpTrendGrowthSummary[, firstLetterUp(FUN)]
    # fitted cubs equation
    tsTmp <- y_fitted_fc[, 2, ]
    tslRes$cubsFittedSummary <- ts(mcmcSummary(x = t(tsTmp), HPDIprob = HPDIprob), start = start, frequency = freq)
    tslRes$cubsFitted <- tslRes$cubsFittedSummary[, firstLetterUp(FUN)]
    
  } 
  
  # return
  if (returnFit) {
    fit$tsl <- tslRes
    attr(fit, "prediction")$n.ahead <- n.ahead
    attr(fit, "prediction")$exogenous <- exogenous
    return(fit)
  } else {
    attr(tslRes, "model") <- class(fit)
    return(tslRes)
  }
  
}

