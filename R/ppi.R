#' Prior Probability Interval for a Normal Distribution
#'
#' Calculates the mean and variance of a normal prior probability interval.
#' Given a lower and upper bound of a prior probability interval
#' this function recovers the implied normal prior parameters.
#' Useful for specifying informative priors
#' on the steady-state parameters (elements of \eqn{\Psi}).
#'
#' @param l Numeric. The lower bound of the prior probability interval.
#' @param u Numeric. The upper bound of the prior probability interval.
#' @param interval Numeric. The prior probability mass within the interval.
#'   Default \code{0.95}, i.e. 95%.
#' @param annualized_growthrate Logical. If \code{TRUE}, treats \code{l} and \code{u} as
#'   bounds on the annualized steady-state growth rate and calculates the implied mean
#'   and variance on the corresponding quarterly/monthly scale. Useful if you are working with
#'   a variable \code{diff(log(x))} (which may be scaled by 100). Default \code{FALSE}.
#' @param freq Integer. The data frequency (e.g. \code{4} for quarterly).
#'   Only used when \code{annualized_growthrate = TRUE}. Default \code{4}.
#'
#' @return A list with two elements: \code{mean} and \code{var}, giving the
#'   mean and variance of the implied normal distribution.
#' @export
#' 
#' @details
#' As an example, say we have a CPI variable \code{CPI <- data$CPI}, where CPI is on the quarterly frequency.
#' Then in the model we work with quarter-on-quarter inflation \code{x <- 100*diff(log(CPI))}.
#' Lets say our prior for annualized steady-state inflation of \code{x} is between 1.7 and 2.3 with 95% probability (with mean at 2).
#' This translates to 0.425 and 0.575 with 95% probability (with mean at 0.5). Clearly it is easier to think of a
#' steady-state prior on the annualized scale, hence the \code{annualized_growthrate} argument.
#' Please see \code{vignette("Homoscedastic-steady-state-BVAR")} on how to use in practise.
#' 
#'
#' @examples
#' ppi(l = 1.7, u = 2.3, interval = 0.95)
#' ppi(l = 1.7, u = 2.3, interval = 0.95, annualized_growthrate = TRUE, freq = 4)
ppi <- function(l, u, interval = 0.95, annualized_growthrate = FALSE, freq = 4) {
  
  if (l >= u) {
    stop("lower bound must be less than upper bound")
  }
  
  if (interval <= 0 || interval >= 1) {
    stop("interval must be between 0 and 1")
  }
  
  if (annualized_growthrate && freq <= 0) {
    stop("freq must be positive")
  }
  
  alpha <- 1 - interval
  z     <- qnorm(1 - alpha / 2)
  
  if (!annualized_growthrate) {
    mu    <- (u + l) / 2
    sigma <- (u - l) / (2 * z)
  } else {
    mu    <- (u + l) / 2 / freq
    sigma <- (u - l) / (2 * z) / freq
  }
  
  return(list(mean = mu, var = sigma^2))
}