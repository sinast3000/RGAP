# -------------------------------------------------------------------------------------------

#' Computes confidence interval of Inverse Gamma distributed variable with given mean and 
#' standard deviation.
#'
#' @param mu A \code{k x 1} vector of means.
#' @param sd A \code{k x 1} vector of standard deviations.
#' @param qlower A \code{k x 1} vector of lower quantiles.
#' @param qupper A \code{k x 1} vector of upper quantiles.
#'
#' @return A \code{2 x k} matrix containing the lower and upper bounds of the intervals.
#'
#' @importFrom stats qgamma
#' @keywords internal
.intervalIGamma <- function(mu, sd, qlower = 0.025, qupper = 1 - qlower) {
  # convert to shape and scale parameters alpha and beta
  alpha <- 2 + mu^2 / sd^2
  beta <- mu * (alpha - 1)
  # compute confidence interval
  matrix(c(
    1 / qgamma(1 - qlower, shape = alpha, rate = beta),
    1 / qgamma(1 - qupper, shape = alpha, rate = beta)
  ),
  2, length(mu),
  byrow = TRUE
  )
}

# ---------------------------------------------------------------------------------------------------

#' Capitalizes the first letter of a string.
#'
#' @param x A character.
#' @keywords internal
firstLetterUp <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}
