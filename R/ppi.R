#' Prior Probability Interval for a Normal Distribution
#'
#' Calculates the mean and variance of a normal prior probability interval.
#' Given a lower and upper bound of a (1-alpha) prior probability interval
#' this function recovers the implied normal prior parameters.
#' Useful for specifying informative priors
#' on the steady-state (Psi) parameters in an intuitive way.
#'
#' @param l Numeric. The lower bound of the prior probability interval.
#' @param u Numeric. The upper bound of the prior probability interval.
#' @param interval Numeric. The prior probability mass within the interval.
#'   Default \code{0.95}, i.e. 95%.
#' @param annualized_growthrate Logical. If \code{TRUE}, converts the interval
#'   from annualized to per-period units by dividing by \code{freq}. Default
#'   \code{FALSE}.
#' @param freq Integer. The data frequency (e.g. \code{4} for quarterly).
#'   Only used when \code{annualized_growthrate = TRUE}. Default \code{4}.
#'
#' @return A list with two elements: \code{mean} and \code{var}, giving the
#'   mean and variance of the implied normal distribution.
#' @export
#'
#' @examples
#' ppi(l = 1, u = 3, interval = 0.95)
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