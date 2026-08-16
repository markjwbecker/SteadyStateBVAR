#' Forecast from a fitted steady-state BVAR model
#'
#' Computes and plots forecasts from a fitted \code{bvar} object. Draws from the 
#' joint predictive distribution are used to construct point forecasts and
#' prediction intervals. Optionally converts monthly or quarterly growth rate forecasts to annual growth
#' rates for selected variables. Optionally overlays the posterior steady state
#' \eqn{\mu_t = \Psi d_t}, computed directly from posterior draws of
#' \eqn{\Psi} applied to the full historical \code{d_t} and future
#' \code{d_pred}. For growth-rate variables, the steady state is converted
#' to YoY terms using the same rolling-sum window as the data/forecast,
#' applied per posterior draw so that any correlation across periods
#' (including the perfect correlation within an unchanged regime) is
#' captured automatically rather than assumed.
#'
#' @param x A steady-state \code{bvar} object that has been passed through \code{\link{fit}}.
#' @param pi Numeric. The prediction interval width. Default \code{0.95}, i.e. 95% prediction interval.
#' @param fcst_type Character. Whether to use \code{"mean"} or \code{"median"}
#'   as the point forecast. Default \code{"mean"}.
#' @param growth_rate_idx Integer vector. Indices of variables of which to convert forecasts to
#'   annual growth rates \eqn{100 (\ln x_{t} - \ln x_{t-f})}, where \eqn{f} is
#'   the frequency of the data (4 for quarterly, 12 for monthly).
#'   Only suitable for variables specified as \eqn{100 (\ln x_{t} - \ln x_{t-1})}, i.e. \code{100*diff(log(x))}.
#'   Computed by summing up to \eqn{f} log first differences. Default is \code{NULL}.
#' @param plot_idx Integer vector. Indices of variables to plot. If \code{NULL}
#'   (default), all variables are plotted. Forecasts are always computed and
#'   returned for all variables, regardless of \code{plot_idx}.
#' @param show_all Logical. If \code{FALSE} (default), the last eight years
#'   of history are shown alongside the forecast. If \code{TRUE}, the full
#'   history is shown.
#' @param ss Logical. If \code{TRUE}, overlays the posterior steady state
#'   \eqn{\mu_t} and its posterior interval. Default \code{FALSE}.
#' @param ss_type Character. Whether to use \code{"mean"} or \code{"median"}
#'   as the steady-state point estimate. Default \code{"mean"}.
#' @param ss_ci Numeric. The posterior interval width for the steady state.
#'   Default \code{0.95}.
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
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1)
#'
#' (fcst <- forecast(bvar_obj, pi = 0.90, show_all = TRUE, ss = TRUE))
#' }
forecast <- function(x, pi = 0.95, fcst_type = c("mean", "median"),
                     growth_rate_idx = NULL, plot_idx = NULL, show_all = FALSE,
                     ss = FALSE, ss_type = c("mean", "median"), ss_ci = 0.95) {
  
  fcst_type <- match.arg(fcst_type)
  ss_type   <- match.arg(ss_type)
  
  Y    <- x$data
  freq <- frequency(Y)
  
  alpha    <- 1 - pi
  alpha_ss <- 1 - ss_ci
  
  if (is.null(plot_idx))
    plot_idx <- 1:ncol(Y)
  
  posterior    <- rstan::extract(x$fit$stan)
  y_pred       <- posterior$y_pred
  y_pred_m     <- apply(y_pred, c(2, 3), fcst_type)
  y_pred_lower <- apply(y_pred, c(2, 3), quantile, probs = alpha / 2)
  y_pred_upper <- apply(y_pred, c(2, 3), quantile, probs = 1 - alpha / 2)
  
  if (isTRUE(ss)) {
    Psi_draws <- posterior$Psi          # [S, k, q]
    dt_full   <- x$setup$dt             # [N_original, q]
    d_pred    <- x$predict$d_pred       # [H, q]
    k         <- x$setup$k
    S         <- dim(Psi_draws)[1]
    N_full    <- nrow(dt_full)
    H_ss      <- nrow(d_pred)
    dt_all    <- rbind(dt_full, d_pred)     # [N_full + H_ss, q]
    
    mu_all_draws <- array(NA, dim = c(S, nrow(dt_all), k))
    for (s in 1:S) {
      mu_all_draws[s, , ] <- dt_all %*% t(Psi_draws[s, , ])
    }
    
    mu_m_all     <- apply(mu_all_draws, c(2, 3), ss_type)
    mu_lower_all <- apply(mu_all_draws, c(2, 3), quantile, probs = alpha_ss / 2)
    mu_upper_all <- apply(mu_all_draws, c(2, 3), quantile, probs = 1 - alpha_ss / 2)
  }
  
  T <- nrow(Y)
  H <- nrow(y_pred_m)
  m <- ncol(Y)
  
  forecast_ret <- matrix(NA, H, m)
  lower_ret    <- matrix(NA, H, m)
  upper_ret    <- matrix(NA, H, m)
  colnames(forecast_ret) <- colnames(Y)
  colnames(lower_ret)    <- colnames(Y)
  colnames(upper_ret)    <- colnames(Y)
  
  time_hist    <- as.numeric(time(Y))
  time_fore    <- seq(tail(time_hist, 1) + 1 / freq, by = 1 / freq, length.out = H)
  mu_line_time <- c(time_hist, time_fore)
  
  legend_labels <- paste0("Prediction (", 100 * pi, "%)")
  legend_cols   <- "blue"
  legend_lty    <- 1
  if (isTRUE(ss)) {
    legend_labels <- c(legend_labels, paste0("Posterior steady state (", 100 * ss_ci, "%)"))
    legend_cols   <- c(legend_cols, "gray40")
    legend_lty    <- c(legend_lty, 2)
  }
  
  for (i in 1:m) {
    smply      <- Y[, i]
    fcst_m     <- y_pred_m[, i]
    fcst_lower <- y_pred_lower[, i]
    fcst_upper <- y_pred_upper[, i]
    
    is_growth <- !is.null(growth_rate_idx) && i %in% growth_rate_idx
    
    if (isTRUE(ss)) {
      
      if (is_growth) {
        # Roll up mu draws over the same freq-quarter window used for the
        # data/forecast YoY conversion, per posterior draw - this respects
        # whatever correlation exists across periods (e.g. mu_t is identical
        # across an unchanged regime) automatically, since it never assumes
        # independence or tries to combine summary statistics.
        n_tot        <- dim(mu_all_draws)[2]
        mu_annual    <- matrix(NA, n_tot, S)
        for (s in 1:S) {
          v <- mu_all_draws[s, , i]
          for (tt in freq:n_tot) {
            mu_annual[tt, s] <- sum(v[(tt - freq + 1):tt])
          }
        }
        mu_line       <- apply(mu_annual, 1, ss_type, na.rm = TRUE)
        mu_line_lower <- apply(mu_annual, 1, quantile, probs = alpha_ss / 2, na.rm = TRUE)
        mu_line_upper <- apply(mu_annual, 1, quantile, probs = 1 - alpha_ss / 2, na.rm = TRUE)
      } else {
        mu_line       <- mu_m_all[, i]
        mu_line_lower <- mu_lower_all[, i]
        mu_line_upper <- mu_upper_all[, i]
      }
    }
    
    if (is_growth) {
      
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
          xlim_vals    <- c(head(time_fore, 1) - 8, tail(time_fore, 1))
          hist_in_plot <- annual_hist[time_hist >= xlim_vals[1] & time_hist <= xlim_vals[2]]
          ylim_vals    <- c(hist_in_plot, annual_lower, annual_upper)
          if (isTRUE(ss)) ylim_vals <- c(ylim_vals, mu_line_lower, mu_line_upper)
          ylim <- range(ylim_vals, na.rm = TRUE)
          
          plot.ts(annual_hist, main = paste(colnames(Y)[i], "(annual)"),
                  xlab = "Time", ylab = NULL,
                  xlim = xlim_vals,
                  ylim = ylim, col = "black", lwd = 2)
          grid(col = "gray", lty = "dotted") 
          points(as.numeric(time_hist), annual_hist, pch = 16, col = "black", cex=0.8)
          points(time_full[-1], m_full[-1], pch = 16, col = "blue", cex=0.8)
        } else {
          ylim_vals <- c(upper_full, lower_full, annual_hist)
          if (isTRUE(ss)) ylim_vals <- c(ylim_vals, mu_line_lower, mu_line_upper)
          plot.ts(annual_hist, main = paste(colnames(Y)[i], "(annual)"),
                  xlab = "Time", ylab = NULL, col = "black", lwd = 2,
                  xlim = c(head(time_hist, 1), tail(time_fore, 1)),
                  ylim = range(ylim_vals, na.rm = TRUE))
        }
        polygon(c(time_full, rev(time_full)), c(upper_full, rev(lower_full)),
                col = rgb(0, 0, 1, 0.2), border = NA)
        lines(time_full, m_full, col = "blue", lwd = 2)
        
        if (isTRUE(ss)) {
          polygon(c(mu_line_time, rev(mu_line_time)), c(mu_line_upper, rev(mu_line_lower)),
                  col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
          lines(mu_line_time, mu_line, col = "gray40", lwd = 2, lty = "dashed")
        }
        
        legend("bottomleft", legend = legend_labels, col = legend_cols,
               lty = legend_lty, lwd = 2, bty = "n", cex = 0.8)
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
          xlim_vals    <- c(head(time_fore, 1) - 8, tail(time_fore, 1))
          hist_in_plot <- smply[time_hist >= xlim_vals[1] & time_hist <= xlim_vals[2]]
          ylim_vals    <- c(hist_in_plot, lower_full, upper_full)
          if (isTRUE(ss)) ylim_vals <- c(ylim_vals, mu_line_lower, mu_line_upper)
          ylim <- range(ylim_vals, na.rm = TRUE)
          
          plot.ts(smply, main = colnames(Y)[i], xlab = "Time", ylab = NULL,
                  xlim = xlim_vals,
                  ylim = ylim, col = "black", lwd = 2)
          grid(col = "gray", lty = "dotted") 
          points(as.numeric(time_hist), smply, pch = 16, col = "black", cex=0.8)
          points(time_full[-1], m_full[-1], pch = 16, col = "blue", cex=0.8)
        } else {
          ylim_vals <- c(upper_full, lower_full, smply)
          if (isTRUE(ss)) ylim_vals <- c(ylim_vals, mu_line_lower, mu_line_upper)
          plot.ts(smply, main = colnames(Y)[i], xlab = "Time", ylab = NULL,
                  col = "black", lwd = 2,
                  xlim = c(head(time_hist, 1), tail(time_fore, 1)),
                  ylim = range(ylim_vals, na.rm = TRUE))
        }
        polygon(c(time_full, rev(time_full)), c(upper_full, rev(lower_full)),
                col = rgb(0, 0, 1, 0.2), border = NA)
        lines(time_full, m_full, col = "blue", lwd = 2)
        
        if (isTRUE(ss)) {
          polygon(c(mu_line_time, rev(mu_line_time)), c(mu_line_upper, rev(mu_line_lower)),
                  col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
          lines(mu_line_time, mu_line, col = "gray40", lwd = 2, lty = "dashed")
        }
        
        legend("bottomleft", legend = legend_labels, col = legend_cols,
               lty = legend_lty, lwd = 2, bty = "n", cex = 0.8)
      }
    }
  }
  
  invisible(list(forecast = forecast_ret, lower = lower_ret, upper = upper_ret))
}