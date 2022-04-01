# -------------------------------------------------------------------------------------------

#' Predictions
#'
#' @description Computes predictions for an object of class \code{fitNAWRU, fitTFP, fitKuttner}
#' estimated via MLE or Bayesian methods.
#'
#' @param n.ahead An integer specifying the prediction horizon. 
#' @param exogenous A character string specifying the computation of exogenous variables 
#'   included in the model. Valid options are \code{exogenous = "mean"} and 
#'   \code{exogenous = "last".}
#' @param returnFit A logical. If \code{TRUE}, an object of the same class as \code{fit}
#'   where the list entry \code{tsl} is replaced. If \code{FALSE}, only the new time series 
#'   list is returned.
#' @inheritParams trendAnchor
#' 
#' @export
predict <- function(fit, n.ahead = 10, exogenous = "mean", returnFit = TRUE) {

  method <- attr(fit, "method")
  if (method != "MLE") {
    predictBayes(fit = fit, n.ahead = n.ahead, exogenous = exogenous, returnFit = returnFit)
  } else {
    predictMLE(fit = fit, n.ahead = n.ahead, exogenous = exogenous, returnFit = returnFit)
  }
  
}

# -------------------------------------------------------------------------------------------

#' Predictions for MLE
#'
#' @description Computes predictions for an object of class \code{fitNAWRU, fitTFP, fitKuttner}
#' estimated via MLE.
#'
#' @inheritParams predict
#' 
#' @details Method only applicable for models estimated via MLE.
#' 
#' @keywords internal
predictMLE <- function(fit, n.ahead = 10, exogenous = "mean", returnFit = TRUE) {

  # adapt model
  model <- fit$SSMfit$model
  model$y <- window(model$y, end = end(model$y) + c(0, n.ahead), extend = TRUE)
  if (dim(model$Z)[3] > 1) {
    Z <- model$Z
    if (exogenous == "last") {
      Z <- Z[, , n]
    } else if (exogenous == "mean") {
      Z <- apply(Z, 1:2, mean)
    }
    # model$Z <- array(c(c(model$Z), rep(model$Z[, , dim(model$Z)[3]], n.ahead)), dim(model$Z) + c(0, 0, n.ahead))
    model$Z <- array(c(c(model$Z), rep(Z, n.ahead)), dim(model$Z) + c(0, 0, n.ahead))
  }
  attr(model, "n") <-   attr(model, "n") + as.integer(n.ahead)
  
  # filter and smoother
  out_fc <- KFS(model, simplify = FALSE, filtering = c("state", "signal"), smoothing = c("state", "signal", "disturbance"))
  
  # get relevant time series
  V <- out_fc$V
  alphahat <- out_fc$alphahat
  att <- out_fc$att
  Ptt <- out_fc$Ptt
  y <- fit$model$SSModel$y
  y <- window(y, end = end(y) + c(0, n.ahead), extend = TRUE)
  
  # state equation parameter matrices
  n <- attr(model, "n")
  TT <- model$T[, , 1]
  R <- model$R[, , 1]
  Q <- model$Q[, , 1]
  Z <- model$Z
  H <- model$H[, , 1]
  
  # # fit$model$SSModel$y <- window(fit$model$SSModel$y, 
  # #                               end = end(fit$model$SSModel$y) + c(0, n.ahead), extend = TRUE)
  # # Z_new <- array(NA, dim(Z) + c(0,0,n.ahead))
  # # Z_new[, , 1:dim(Z)[3]] <- Z
  # # Z_new[, , -c(1:dim(Z)[3])] <- Z[, , n]
  # # fit$model$SSModel$Z <- Z_new
  # # out <- KFS(fit$model$SSModel, simplify = FALSE, filtering = c("state", "signal"), smoothing = c("state", "signal", "disturbance"))
  # # 
  # if (exogenous == "last") {
  #   Z <- Z[, , n]
  # } else if (exogenous == "mean") {
  #   Z <- apply(Z, 1:2, mean)
  # }
  # 
  # # smoothed state and variance
  # alphahat <- fit$SSMout$alphahat
  # V <- fit$SSMout$V
  # alphahat <- window(alphahat, end = end(alphahat) + c(0, n.ahead), extend = TRUE)
  # V <- array(c(V, rep(NA, NROW(V) *  NCOL(V) * n.ahead)), c(NROW(V), NCOL(V),n + n.ahead))
  # 
  # # filtered state and variance
  # att <- fit$SSMout$att
  # Ptt <- fit$SSMout$Ptt
  # att <- window(att, end = end(alphahat) + c(0, n.ahead), extend = TRUE)
  # Ptt <- array(c(Ptt, rep(NA, NROW(V) *  NCOL(Ptt) * n.ahead)), c(NROW(Ptt), NCOL(Ptt),n + n.ahead))

  # # observables
  # y <- fit$model$SSModel$y
  # y <- window(y, end = end(y) + c(0, n.ahead), extend = TRUE)

  # # helper matrices
  # RQR <- R %*% Q %*% t(R)
  # 
  # # predict
  # for (tt in 1:n.ahead) {
  #   # smoothed and filtered state
  #   alphahat[n + tt, ] <- TT %*% alphahat[n + tt - 1, ]
  #   att[n + tt, ] <- TT %*% att[n + tt - 1, ]
  #   
  #   # variances of smoothed and filtered state
  #   # V[, , n + tt] <- Ptt[, , n + tt] <- RQRtt + (TT^tt) %*% V[, , n] %*% t(TT^tt)
  #   V[, , n + tt] <- Ptt[, , n + tt] <- RQR + TT %*% V[, , n + tt - 1] %*% t(TT)
  #   # V[, , n] equals Ptt[, , n]
  #   
  # }

  
  # predict observations
  # y_fitted <- ts(alphahat %*% t(Z), start = start(alphahat), frequency = frequency(alphahat))
  if (dim(model$Z)[3] > 1) {
    y_fitted <- ts(t(sapply(1:dim(Z)[3], function(x) alphahat[x, ] %*% t(Z[, , x]))), 
              start = start(alphahat), frequency = frequency(alphahat))
    V_y <- array(sapply(1:dim(V)[3], function(x) Z[, , x] %*% V[, , x] %*% t(Z[, , x]) + H), c(rep(dim(Z)[1], 2), dim(V)[3]))
  } else {
    y_fitted <- ts(alphahat %*% t(Z[, , 1]), start = start(alphahat), frequency = frequency(alphahat))
    V_y <- array(sapply(1:dim(V)[3], function(x) Z[, , 1] %*% V[, , x] %*% t(Z[, , 1]) + H), c(rep(dim(Z)[1], 2), dim(V)[3]))
  }
  V_y <- array(sapply(1:dim(V)[3], function(x) { # set values for which observations exist to NA
    mat <- V_y[, , x]
    diag(mat)[!is.na(y)[x,]] <- NA
    mat
  }), dim(V_y))
  y[is.na(y)] <- y_fitted[is.na(y)] # fill NAs with predictions
  
  # confidence intervals smoothed state
  alphahat_se <- ts(t(sqrt(apply(V, 3, diag))), start = start(alphahat), frequency = frequency(alphahat))
  colnames(alphahat_se) <- colnames(alphahat)

  # confidence intervals filtered state
  att_se <- ts(t(sqrt(apply(Ptt, 3, diag))), start = start(att), frequency = frequency(att))
  colnames(att_se) <- colnames(att)
  
  # confidence interval obervation equation
  y_se <- ts(t(sqrt(apply(V_y, 3, diag))), start = start(alphahat), frequency = frequency(alphahat))
  colnames(y_se) <- colnames(y)

  # new ts list
  tsl_fc <- within(fit$tsl, {
    obs <- y
    obsSE <- y_se
    stateSmoothed <- alphahat
    stateFiltered <- att
    stateSmoothedSE <- alphahat_se
    stateFilteredSE <- att_se
  })
  
  # fill ts list depending on different objects
  trend <- alphahat[, "trend"]
  trend_se <- alphahat_se[, "trend"]
  cycle <- alphahat[, "cycle"]
  cycle_se <- alphahat_se[, "cycle"]
  
  # out_fc <- fit$SSMout
  # out_fc$alphahat <- alphahat
  # out_fc$V <- V
  out_fc$V_y <- V_y
  # out_fc$dims$n <- fit$SSMout$dims$n + n.ahead
  if (inherits(fit, "NAWRUfit")) {
    
    tsl_fc <- within(tsl_fc, {
      nawru <- trend
      nawruSE <- trend_se
    })
    
  } else if (inherits(fit, "TFPfit")) {
    
    tsl_fc <- within(tsl_fc, {
      tfpTrend <- exp(trend)
      tfpTrendGrowth <- growth(exp(trend))
      tfp <- exp(y[, 1])
      tfpGrowth <- growth(exp(y[, 1]))
    })
    
    tsl_fc[c("tfpTrendSE", "tfpTrendGrowthSE")] <- .deltaMethodState(out = out_fc, nameState = "trend")[-1]
    tsl_fc$tfpSE <- .deltaMethodObs(out = out_fc, y = y, nameObs = "logtfp")$expObsSE
    tsl_fc$tfpGrowthSE <- .deltaMethodObs(out = out_fc, y = y, nameObs = "logtfp")$diffObsSE
    
  } else if (inherits(fit, "KuttnerFit")) {
    
    tsl_fc <- within(tsl_fc, {
      potential <- exp(trend)
      potentialGrowth <- growth(exp(trend))
      gap <- cycle * 100
      gdp <- exp(y[, 1])
      gdpGrowth <- growth(exp(y[, 1]))
    })
    
    tsl_fc[[c("potentialGrowthSE")]] <- .deltaMethodState(out = out_fc, nameState = "trend")$diffStateSE
    tsl_fc[[c("potentialSE")]] <- .deltaMethodState(out = out_fc, nameState = "trend")$expStateSE
    tsl_fc[[c("gapSE")]] <- 100 * .deltaMethodState(out = out_fc, nameState = "cycle")$StateSE
    tsl_fc$gdpSE <- .deltaMethodObs(out = out_fc, y = y, nameObs = "loggdp")$expObsSE
    tsl_fc$gdpGrowthSE <- .deltaMethodObs(out = out_fc, y = y, nameObs = "loggdp")$diffObsSE
  }
  
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
#' @description Computes predictions for an object of class \code{fitNAWRU, fitTFP, fitKuttner}
#' estimated via MLE.
#'
#' @inheritParams predict
#' 
#' @details Method only applicable for models estimated via MLE.
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
    # R <- SSModel$R[, , 1]
    # Q <- SSModel$Q[, , 1]
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

