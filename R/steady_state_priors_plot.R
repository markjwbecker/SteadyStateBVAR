#' Plot steady-state priors
#'
#' Produces time series plots of the data along with the implied prior
#' interval for the steady state \eqn{\mu_t = \Psi d_t}, based on the
#' normal prior \eqn{\mathrm{vec}(\Psi) \sim \mathrm{N}(\theta_\Psi, \Omega_\Psi)}
#' specified in \code{\link{priors}}. Since \eqn{\Omega_\Psi} is diagonal,
#' the elements of \eqn{\Psi} are independent under the prior, so the
#' variance of \eqn{\mu_t = \Psi d_t} is obtained by propagating the
#' marginal prior variances through the linear combination implied by
#' \eqn{d_t}, rather than by combining marginal intervals directly.
#'
#' @param x A steady-state \code{bvar} object that has been passed through
#'   \code{\link{priors}}.
#' @param interval Numeric. The prior interval width. Default \code{0.95},
#'   i.e. a 95\% prior interval for the steady state.
#' @param growth_rate_idx Integer vector. Indices of variables specified as
#'   \code{100*diff(log(x))} for which the historical series and the
#'   steady-state prior are converted to annual terms: the historical series
#'   via the same rolling sum used in \code{\link{forecast}}, and the
#'   steady-state prior mean and interval by scaling
#'   \eqn{\theta_\Psi} and the prior bounds by the frequency of the data
#'   (4 for quarterly, 12 for monthly). Default is \code{NULL}.
#' @param plot_idx Integer vector. Indices of variables to plot. If
#'   \code{NULL} (default), all variables are plotted.
#'
#' @return Invisibly returns a list with three matrices, \code{lower},
#'   \code{mean}, and \code{upper}, each of dimension \code{T x k} giving the
#'   steady-state prior bounds and mean over the historical sample.
#' @export
#'
#' @examples
#' \donttest{
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p = 1, deterministic = "constant")
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1, 2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    Jeffreys = TRUE,
#'                    SV = FALSE,
#'                    SV_type = NULL,
#'                    SV_priors = NULL)
#'
#' steady_state_priors_plot(bvar_obj)
#' }
steady_state_priors_plot <- function(x,
                                     interval = 0.95,
                                     growth_rate_idx = NULL,
                                     plot_idx = NULL) {
  
  if (!inherits(x, "bvar")) stop("x must be a 'bvar' object")
  if (is.null(x$setup))  stop("x must be passed through setup")
  if (is.null(x$priors)) stop("x must be passed through priors")
  if (is.null(x$priors$theta_Psi) || is.null(x$priors$Omega_Psi)) {
    stop("theta_Psi and Omega_Psi must be specified in priors")
  }
  
  Y         <- x$data
  freq      <- frequency(Y)
  theta_Psi <- x$priors$theta_Psi
  Omega_Psi <- x$priors$Omega_Psi
  dt        <- x$setup$dt
  k         <- x$setup$k
  q         <- x$setup$q
  
  if (is.null(plot_idx)) plot_idx <- 1:k
  
  alpha  <- 1 - interval
  z      <- qnorm(1 - alpha / 2)
  sd_Psi <- sqrt(diag(Omega_Psi))
  
  Psi_mean <- matrix(theta_Psi, k, q)
  Psi_sd   <- matrix(sd_Psi, k, q)
  
  ss_mean <- t(Psi_mean %*% t(dt))
  ss_var  <- t((Psi_sd^2) %*% t(dt^2))
  ss_sd   <- sqrt(ss_var)
  
  ss_lower <- ss_mean - z * ss_sd
  ss_upper <- ss_mean + z * ss_sd
  
  time_hist <- as.numeric(time(Y))
  
  legend_labels <- paste0("Prior steady state (", 100 * interval, "%)")
  legend_cols   <- "gray40"
  legend_lty    <- 2
  
  for (i in 1:k) {
    
    smply      <- Y[, i]
    line_mean  <- ss_mean[, i]
    line_lower <- ss_lower[, i]
    line_upper <- ss_upper[, i]
    
    if (!is.null(growth_rate_idx) && i %in% growth_rate_idx) {
      
      annual_hist <- rep(NA, length(smply))
      for (t in freq:length(smply)) {
        annual_hist[t] <- sum(smply[(t - (freq - 1)):t])
      }
      annual_hist <- ts(annual_hist, start = start(Y), frequency = freq)
      
      smply      <- annual_hist
      line_mean  <- line_mean  * freq
      line_lower <- line_lower * freq
      line_upper <- line_upper * freq
    }
    
    if (i %in% plot_idx) {
      
      main_lab <- if (!is.null(growth_rate_idx) && i %in% growth_rate_idx) {
        paste(colnames(Y)[i], "(annual)")
      } else {
        colnames(Y)[i]
      }
      
      plot.ts(smply, main = main_lab, xlab = "Time", ylab = NULL,
              col = "black", lwd = 2,
              ylim = range(c(smply, line_lower, line_upper), na.rm = TRUE))
      
      grid(col = "gray", lty = "dotted")
      
      polygon(c(time_hist, rev(time_hist)), c(line_upper, rev(line_lower)),
              col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
      lines(time_hist, line_mean, col = "gray40", lwd = 2, lty = "dashed")
      
      legend("bottomleft", legend = legend_labels, col = legend_cols,
             lty = legend_lty, lwd = 2, bty = "n", cex = 0.8)
    }
  }
  
  invisible(list(lower = ss_lower, mean = ss_mean, upper = ss_upper))
}