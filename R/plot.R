plot_forecast <- function(fit, Y, ci=0.95, fcst_type=c("mean", "median"), 
                          growth_rate_idx=NULL, plot_idx=NULL,gibbs=FALSE)
  {
  if (gibbs==FALSE){
  posterior <- rstan::extract(fit)
  Y_pred <- posterior$Y_pred
  } else {
  Y_pred <- fit$fcst_draws
  }
  Y_pred_m <- apply(Y_pred, c(2, 3), fcst_type)
  alpha <- 1 - ci
  Y_pred_lower <- apply(Y_pred, c(2, 3), quantile, probs = alpha/2)
  Y_pred_upper <- apply(Y_pred, c(2, 3), quantile, probs = 1 - alpha/2)
  
  T <- nrow(Y)
  H <- nrow(Y_pred_m)
  m <- ncol(Y)
  
  time_hist <- time(Y)
  time_fore <- seq(tail(time_hist, 1) + 0.25, by = 0.25, length.out = H)
  
  if (ncol(Y) < 4){
    par(mfcol = c(ncol(Y), 1))
  } else if ((ncol(Y) < 9)){
    par(mfcol = c(4,2))
  } else {
    par(mfcol = c(1,1))
  }
  par(mfcol = c(1,1))
  on.exit(par(mfrow = c(1,1)))
  
  if (is.null(plot_idx)) plot_idx <- 1:ncol(Y)
  for (i in plot_idx) {
    smply <- Y[, i]
    fcst_m <- Y_pred_m[, i]
    fcst_lower <- Y_pred_lower[, i]
    fcst_upper <- Y_pred_upper[, i]
    
    if (!is.null(growth_rate_idx) && i %in% growth_rate_idx) {
      
      annual_hist <- rep(NA, length(smply))
      for(t in 4:length(smply)){
        annual_hist[t] <- sum(smply[(t-3):t])
      }
      annual_hist <- annual_hist / 4
      annual_hist <- ts(annual_hist, start = start(Y), frequency = 4)
      
      last_obs <- tail(smply,3)
      all_fcst <- c(last_obs, fcst_m)
      
      annual_fcst <- rep(NA, H)
      for(t_h in 1:H){
        annual_fcst[t_h] <- sum(all_fcst[t_h:(t_h+3)])
      }
      annual_fcst <- annual_fcst / 4
      
      
      all_lower <- c(last_obs, fcst_lower)
      all_upper <- c(last_obs, fcst_upper)
      
      annual_lower <- rep(NA,H)
      annual_upper <- rep(NA,H)
      for(t_h in 1:H){
        annual_lower[t_h] <- sum(all_lower[t_h:(t_h+3)])
        annual_upper[t_h] <- sum(all_upper[t_h:(t_h+3)])
      }
      annual_lower <- annual_lower / 4
      annual_upper <- annual_upper / 4
      
      time_full <- c(tail(time_hist, 1), time_fore)
      
      m_full <- c(tail(annual_hist, 1), annual_fcst)
      lower_full <- c(tail(annual_hist, 1), annual_lower)
      upper_full <- c(tail(annual_hist, 1), annual_upper)
      
      ylim <- range(c(annual_lower, annual_upper),na.rm=TRUE)
      
      plot.ts(annual_hist, main = paste(colnames(Y)[i], "(annual)"), xlab = "Time", ylab = NULL,
              xlim = c(1998,2008),
              ylim = ylim,col = "black", lwd = 2)
      
      polygon(
        c(time_full, rev(time_full)),
        c(upper_full, rev(lower_full)),
        col = rgb(0, 0, 1, 0.2),
        border = NA
      )
      
      lines(time_full, m_full, col = "blue", lwd = 2)
      
    } else {
      
      time_full <- c(tail(time_hist, 1), time_fore)
      m_full <- c(tail(smply, 1), fcst_m)
      lower_full <- c(tail(smply, 1), fcst_lower)
      upper_full <- c(tail(smply, 1), fcst_upper)
      
      ylim <- range(c(lower_full, upper_full))
      
      plot.ts(smply, main = colnames(Y)[i], xlab = "Time", ylab = NULL,
              xlim = c(1998,2008),
              ylim = ylim,col = "black", lwd = 2)
      
      polygon(
        c(time_full, rev(time_full)),
        c(upper_full, rev(lower_full)),
        col = rgb(0, 0, 1, 0.2),
        border = NA
      )
      lines(time_full, m_full, col = "blue", lwd = 2)
    }
  }
}
