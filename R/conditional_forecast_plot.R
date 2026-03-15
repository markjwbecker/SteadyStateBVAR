conditional_forecast_plot <- function(cond_fcst_obj,
                               bvar_obj,
                               ci = 0.95,
                               fcst_type = c("mean", "median"),
                               growth_rate_idx = NULL,
                               plot_idx = NULL,
                               estimation = c("stan", "gibbs")) {
  
  estimation <- match.arg(estimation)
  fcst_type  <- match.arg(fcst_type)
  
  if (estimation == "gibbs") {
    y_pred_uncond <- bvar_obj$fit$gibbs$fcst_draws
    y_pred_m_uncond <- apply(y_pred_uncond, c(1, 2), fcst_type)
  } else {
    posterior <- rstan::extract(bvar_obj$fit$stan)
    y_pred_uncond <- posterior$y_pred
    y_pred_m_uncond <- apply(y_pred_uncond, c(2, 3), fcst_type)
  }
  
  Y         <- bvar_obj$data
  freq      <- frequency(Y)
  m         <- ncol(Y)
  if (is.null(plot_idx)) plot_idx <- 1:m
  
  y_pred_m     <- cond_fcst_obj$cond_point
  y_pred_lower <- cond_fcst_obj$cond_lower
  y_pred_upper <- cond_fcst_obj$cond_upper
  
  T          <- nrow(Y)
  H          <- nrow(y_pred_m)
  time_hist  <- as.numeric(time(Y))
  time_fore  <- seq(tail(time_hist, 1) + 1/freq, by = 1/freq, length.out = H)
  
  forecast_ret <- matrix(NA, H, m)
  lower_ret    <- matrix(NA, H, m)
  upper_ret    <- matrix(NA, H, m)
  colnames(forecast_ret) <- colnames(Y)
  colnames(lower_ret)    <- colnames(Y)
  colnames(upper_ret)    <- colnames(Y)
  
  for (i in plot_idx) {
    smply      <- Y[, i]
    fcst_m     <- y_pred_m[, i]
    fcst_lower <- y_pred_lower[, i]
    fcst_upper <- y_pred_upper[, i]
    uncond_m <- y_pred_m_uncond[, i]
    
    if (!is.null(growth_rate_idx) && i %in% growth_rate_idx) {
      annual_hist <- rep(NA, length(smply))
      for (t in freq:length(smply)) {
        annual_hist[t] <- sum(smply[(t - (freq - 1)):t])
      }
      annual_hist <- ts(annual_hist, start = start(Y), frequency = freq)
      last_obs    <- tail(smply, (freq - 1))
      all_fcst    <- c(last_obs, fcst_m)
      all_fcst_uncond    <- c(last_obs, uncond_m)
      annual_fcst <- rep(NA, H)
      annual_fcst_uncond <- rep(NA, H)
      for (t_h in 1:H) {
        annual_fcst[t_h] <- sum(all_fcst[t_h:(t_h + (freq - 1))])
        annual_fcst_uncond[t_h] <- sum(all_fcst_uncond[t_h:(t_h + (freq - 1))])
      }
      all_lower    <- c(last_obs, fcst_lower)
      all_upper    <- c(last_obs, fcst_upper)
      annual_lower <- rep(NA, H)
      annual_upper <- rep(NA, H)
      for (t_h in 1:H) {
        annual_lower[t_h] <- sum(all_lower[t_h:(t_h + (freq - 1))])
        annual_upper[t_h] <- sum(all_upper[t_h:(t_h + (freq - 1))])
      }
      forecast_ret[, i] <- annual_fcst
      lower_ret[, i]    <- annual_lower
      upper_ret[, i]    <- annual_upper
      time_full  <- c(tail(time_hist, 1), time_fore)
      m_full     <- c(tail(annual_hist, 1), annual_fcst)
      lower_full <- c(tail(annual_hist, 1), annual_lower)
      upper_full <- c(tail(annual_hist, 1), annual_upper)
      uncond_full <- c(tail(annual_hist, 1), annual_fcst_uncond)
      plot.ts(annual_hist, main = paste(colnames(Y)[i], "(annual)"), xlab = "Time", ylab = NULL,
              col = "black", lwd = 2,
              xlim = c(head(time_hist, 1), tail(time_fore, 1)),
              ylim = range(upper_full, lower_full, annual_hist, uncond_full, na.rm = TRUE))
    } else {
      forecast_ret[, i] <- fcst_m
      lower_ret[, i]    <- fcst_lower
      upper_ret[, i]    <- fcst_upper
      time_full  <- c(tail(time_hist, 1), time_fore)
      m_full     <- c(tail(smply, 1), fcst_m)
      lower_full <- c(tail(smply, 1), fcst_lower)
      upper_full <- c(tail(smply, 1), fcst_upper)
      uncond_full <- c(tail(smply, 1), uncond_m)
      plot.ts(smply, main = colnames(Y)[i], xlab = "Time", ylab = NULL,
              col = "black", lwd = 2,
              xlim = c(head(time_hist, 1), tail(time_fore, 1)),
              ylim = range(upper_full, lower_full, smply, uncond_full))
    }
    polygon(c(time_full, rev(time_full)), c(upper_full, rev(lower_full)),
            col = rgb(0, 0, 1, 0.2), border = NA)
    lines(time_full, m_full, col = "blue", lwd = 2)
    lines(time_full, uncond_full, col = "red", lwd = 1, lty=1)
  }
  return(list(forecast = forecast_ret, lower = lower_ret, upper = upper_ret))
}