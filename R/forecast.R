#' Forecast from a fitted steady-state BVAR model
#'
#' Computes and plots forecasts from a fitted \code{bvar} object. Draws from the 
#' joint predictive distribution are used to construct point forecasts and
#' prediction intervals. Optionally converts monthly or quarterly growth rate forecasts to annual growth
#' rates for selected variables.
#'
#' @param x A steady-state \code{bvar} object that has been passed through \code{\link{fit}}.
#' @param pi Numeric. The prediction interval width. Default \code{0.95}, i.e. 95% prediction interval.
#' @param fcst_type Character. Whether to use \code{"mean"} or \code{"median"}
#'   as the point forecast. Default \code{"mean"}.
#' @param growth_rate_idx Integer vector. Indices of variables of which to convert forecasts to
#'   annual growth rates \eqn{\ln x_{t} - \ln x_{t-f}}, where \eqn{f} is
#'   the frequency of the data (4 for quarterly, 12 for monthly).
#'   Only suitable for variables specified as \eqn{\ln x_{t} - \ln x_{t-1}}, i.e.
#'   \code{diff(log(x))} or \code{100*diff(log(x))}.
#'   Computed by summing up to \eqn{f} log first differences. Default is \code{NULL}.
#' @param plot_idx Integer vector. Indices of variables to plot. If \code{NULL}
#'   (default), all variables are plotted. Forecasts are always computed and
#'   returned for all variables, regardless of \code{plot_idx}.
#' @param show_all Logical. If \code{FALSE} (default), only the last two years
#'   of history are shown alongside the forecast. If \code{TRUE}, the full
#'   history is shown.
#'
#' @return Invisibly returns a list with three matrices: \code{forecast}, \code{lower}, and
#'   \code{upper}, each of dimension \code{H x k} where \code{H} is the
#'   forecast horizon and \code{k} is the number of variables.
#' @export
#'
#' @examples
#' \donttest{
#' #homoscedastic with Jeffreys prior
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    Jeffreys = TRUE,
#'                    SV = FALSE,
#'                    SV_type = NULL,
#'                    SV_priors = NULL)
#'                    
#' bvar_obj <- fit(bvar_obj,
#'                 H = 8,
#'                 d_pred = matrix(rep(1,8)),
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1)
#'
#' (fcst <- forecast(bvar_obj, pi = 0.90, show_all = TRUE))
#' }
forecast <- function(x, pi = 0.95, fcst_type = c("mean", "median"),
                     growth_rate_idx = NULL, plot_idx = NULL, show_all = FALSE) {
  
  fcst_type <- match.arg(fcst_type)
  
  Y    <- x$data
  freq <- frequency(Y)
  
  if (is.null(plot_idx))
    plot_idx <- 1:ncol(Y)
  
  posterior    <- rstan::extract(x$fit$stan)
  y_pred       <- posterior$y_pred
  alpha        <- 1 - pi
  y_pred_m     <- apply(y_pred, c(2, 3), fcst_type)
  y_pred_lower <- apply(y_pred, c(2, 3), quantile, probs = alpha / 2)
  y_pred_upper <- apply(y_pred, c(2, 3), quantile, probs = 1 - alpha / 2)
  
  T <- nrow(Y)
  H <- nrow(y_pred_m)
  m <- ncol(Y)
  
  forecast_ret <- matrix(NA, H, m)
  lower_ret    <- matrix(NA, H, m)
  upper_ret    <- matrix(NA, H, m)
  colnames(forecast_ret) <- colnames(Y)
  colnames(lower_ret)    <- colnames(Y)
  colnames(upper_ret)    <- colnames(Y)
  
  time_hist <- as.numeric(time(Y))
  time_fore <- seq(tail(time_hist, 1) + 1 / freq, by = 1 / freq, length.out = H)
  
  for (i in 1:m) {
    smply      <- Y[, i]
    fcst_m     <- y_pred_m[, i]
    fcst_lower <- y_pred_lower[, i]
    fcst_upper <- y_pred_upper[, i]
    
    if (!is.null(growth_rate_idx) && i %in% growth_rate_idx) {
      
      annual_hist <- rep(NA, length(smply))
      for (t in freq:length(smply)) {
        annual_hist[t] <- sum(smply[(t - (freq - 1)):t])
      }
      annual_hist <- ts(annual_hist, start = start(Y), frequency = freq)
      
      last_obs   <- tail(smply, (freq - 1))
      all_fcst   <- c(last_obs, fcst_m)
      all_lower  <- c(last_obs, fcst_lower)
      all_upper  <- c(last_obs, fcst_upper)
      
      annual_fcst  <- rep(NA, H)
      annual_lower <- rep(NA, H)
      annual_upper <- rep(NA, H)
      
      for (t_h in 1:H) {
        annual_fcst[t_h]  <- sum(all_fcst[t_h:(t_h + (freq - 1))])
        annual_lower[t_h] <- sum(all_lower[t_h:(t_h + (freq - 1))])
        annual_upper[t_h] <- sum(all_upper[t_h:(t_h + (freq - 1))])
      }
      
      forecast_ret[, i] <- annual_fcst
      lower_ret[, i]    <- annual_lower
      upper_ret[, i]    <- annual_upper
      
      if (i %in% plot_idx) {
        time_full  <- c(tail(time_hist, 1), time_fore)
        m_full     <- c(tail(annual_hist, 1), annual_fcst)
        lower_full <- c(tail(annual_hist, 1), annual_lower)
        upper_full <- c(tail(annual_hist, 1), annual_upper)
        
        if (isFALSE(show_all)) {
          xlim_vals    <- c(head(time_fore, 1) - (freq * 2), tail(time_fore, 1))
          hist_in_plot <- annual_hist[time_hist >= xlim_vals[1] & time_hist <= xlim_vals[2]]
          ylim         <- range(c(hist_in_plot, annual_lower, annual_upper), na.rm = TRUE)
          ymin         <- floor(ylim[1] * 2) / 2
          ymax         <- ceiling(ylim[2] * 2) / 2
          yticks       <- seq(ymin, ymax, by = 0.5)
          
          plot.ts(annual_hist, main = paste(colnames(Y)[i], "(annual)"),
                  xlab = "Time", ylab = NULL,
                  xlim = c(head(time_fore, 1) - (freq * 2), tail(time_fore, 1)),
                  ylim = ylim, col = "black", lwd = 2, yaxt = "n")
          abline(v = seq(head(time_hist, 1), tail(time_fore, 1), by = 1 / freq),
                 col = "gray", lty = 2)
          abline(h = yticks, col = "gray", lty = 2)
          axis(side = 2, at = yticks, labels = yticks, las = 1)
          points(as.numeric(time_hist), annual_hist, pch = 16, col = "black")
          points(time_full[-1], m_full[-1], pch = 16, col = "blue")
        } else {
          plot.ts(annual_hist, main = paste(colnames(Y)[i], "(annual)"),
                  xlab = "Time", ylab = NULL, col = "black", lwd = 2,
                  xlim = c(head(time_hist, 1), tail(time_fore, 1)),
                  ylim = range(upper_full, lower_full, annual_hist, na.rm = TRUE))
        }
        polygon(c(time_full, rev(time_full)), c(upper_full, rev(lower_full)),
                col = rgb(0, 0, 1, 0.2), border = NA)
        lines(time_full, m_full, col = "blue", lwd = 2)
      }
      
    } else {
      
      forecast_ret[, i] <- fcst_m
      lower_ret[, i]    <- fcst_lower
      upper_ret[, i]    <- fcst_upper
      
      if (i %in% plot_idx) {
        time_full  <- c(tail(time_hist, 1), time_fore)
        m_full     <- c(tail(smply, 1), fcst_m)
        lower_full <- c(tail(smply, 1), fcst_lower)
        upper_full <- c(tail(smply, 1), fcst_upper)
        
        if (isFALSE(show_all)) {
          xlim_vals    <- c(head(time_fore, 1) - (freq * 2), tail(time_fore, 1))
          hist_in_plot <- smply[time_hist >= xlim_vals[1] & time_hist <= xlim_vals[2]]
          ylim         <- range(c(hist_in_plot, lower_full, upper_full), na.rm = TRUE)
          ymin         <- floor(ylim[1] * 2) / 2
          ymax         <- ceiling(ylim[2] * 2) / 2
          yticks       <- seq(ymin, ymax, by = 0.5)
          
          plot.ts(smply, main = colnames(Y)[i], xlab = "Time", ylab = NULL,
                  xlim = c(head(time_fore, 1) - (freq * 2), tail(time_fore, 1)),
                  ylim = ylim, col = "black", lwd = 2, yaxt = "n")
          abline(v = seq(head(time_hist, 1), tail(time_fore, 1), by = 1 / freq),
                 col = "gray", lty = 2)
          abline(h = yticks, col = "gray", lty = 2)
          axis(side = 2, at = yticks, labels = yticks, las = 1)
          points(as.numeric(time_hist), smply, pch = 16, col = "black")
          points(time_full[-1], m_full[-1], pch = 16, col = "blue")
        } else {
          plot.ts(smply, main = colnames(Y)[i], xlab = "Time", ylab = NULL,
                  col = "black", lwd = 2,
                  xlim = c(head(time_hist, 1), tail(time_fore, 1)),
                  ylim = range(upper_full, lower_full, smply))
        }
        polygon(c(time_full, rev(time_full)), c(upper_full, rev(lower_full)),
                col = rgb(0, 0, 1, 0.2), border = NA)
        lines(time_full, m_full, col = "blue", lwd = 2)
      }
    }
  }
  
  invisible(list(forecast = forecast_ret, lower = lower_ret, upper = upper_ret))
}