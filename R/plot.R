plot_forecast <- function(fit, Y, ci=0.95,fcst_type=c("mean", "median")) {
  posterior <- rstan::extract(fit)
  Y_pred <- posterior$Y_pred
  Y_pred_m <- apply(Y_pred, c(2, 3), fcst_type)
  alpha <- 1 - ci
  Y_pred_lower <- apply(Y_pred, c(2, 3), quantile, probs = alpha/2)
  Y_pred_upper <- apply(Y_pred, c(2, 3), quantile, probs = 1 - alpha/2)
  
  T <- nrow(Y)
  H <- nrow(Y_pred_m)
  m <- ncol(Y)
  
  time_hist <- time(Y)
  time_fore <- seq(tail(time_hist, 1) + 1/4, by = 1/4, length.out = H)
  
  par(mfrow = c(ncol(Y), 1))
  
  for (i in 1:ncol(Y)) {
    smply <- Y[, i]
    fcst_m <- Y_pred_m[, i]
    fcst_lower <- Y_pred_lower[, i]
    fcst_upper <- Y_pred_upper[, i]
    
    time_full <- c(tail(time_hist, 1), time_fore)
    m_full <- c(tail(smply, 1), fcst_m)
    lower_full <- c(tail(smply, 1), fcst_lower)
    upper_full <- c(tail(smply, 1), fcst_upper)
    
    ylim <- range(c(smply, lower_full, upper_full))
    
    plot.ts(smply, main = colnames(Y)[i], xlab = "Time", ylab = NULL,
            xlim = c(time_hist[1], tail(time_fore, 1)),
            ylim = ylim, col = "black", lwd = 2)
    
    polygon(
      c(time_full, rev(time_full)),
      c(upper_full, rev(lower_full)),
      col = rgb(0, 0, 1, 0.2),
      border = NA
    )
    lines(time_full, m_full, col = "blue", lwd = 2)
    
  }
}
