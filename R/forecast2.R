forecast2 <- function(x, fcst_type=c("mean", "median"), plot_idx=NULL, xlim, ylim){
  
  alpha <- 1 - ci
  
  Y_pred <- bvar_obj$fit$gibbs$fcst_draws
  Y_pred_m <- apply(Y_pred, c(1, 2), fcst_type)
  Y_pred_sd <- apply(Y_pred, c(1, 2), sd)
  Y_pred_lower <- Y_pred_m - 1 * Y_pred_sd
  Y_pred_upper <- Y_pred_m + 1 * Y_pred_sd
  
  Y <- bvar_obj$data
  freq <- frequency(Y)
  T <- nrow(Y)
  H <- bvar_obj$H
  m <- ncol(Y)
  time_hist <- time(Y)
  time_fore <- seq(tail(time_hist, 1) + 1/freq, by = 1/freq, length.out = H)
  i = plot_idx
  smply <- Y[, i]
  fcst_m <- Y_pred_m[, i]
  fcst_lower <- Y_pred_lower[, i]
  fcst_upper <- Y_pred_upper[, i]
  
  time_full <- c(tail(time_hist, 1), time_fore)
  m_full <- c(tail(smply, 1), fcst_m)
  lower_full <- c(tail(smply, 1), fcst_lower)
  upper_full <- c(tail(smply, 1), fcst_upper)
  
  ymin <- floor(ylim[1]*2)/2
  ymax <- ceiling(ylim[2]*2)/2
  yticks <- seq(ymin, ymax, by = 1)
  plot.ts(smply, main = colnames(Y)[i],ylab = NULL,
          xlim=xlim,
          ylim=ylim,
          col = "black", lwd = 2, yaxt = "n")
  
  axis(side = 2, at = yticks, labels = yticks, las = 1)
  
  lines(time_full[-1], lower_full[-1], col = "red", lwd = 2,lty=3)
  lines(time_full[-1], upper_full[-1], col = "red", lwd = 2,lty=3)
  lines(time_full, m_full, col = "red", lwd = 2)
  
  posterior_stan <- rstan::extract(x$fit$stan)
  Y_pred2 <- posterior_stan$Y_pred
  Y_pred_m2 <- apply(Y_pred2, c(2, 3), fcst_type)
  Y_pred_sd2 <- apply(Y_pred2, c(2, 3), sd)
  Y_pred_lower2 <- Y_pred_m2 - 1 * Y_pred_sd2
  Y_pred_upper2 <- Y_pred_m2 + 1 * Y_pred_sd2
  
  fcst_m2 <- Y_pred_m2[, i]
  fcst_lower2 <- Y_pred_lower2[, i]
  fcst_upper2 <- Y_pred_upper2[, i]
  
  m_full2 <- c(tail(smply, 1), fcst_m2)
  lower_full2 <- c(tail(smply, 1), fcst_lower2)
  upper_full2 <- c(tail(smply, 1), fcst_upper2)
  lines(time_full[-1], lower_full2[-1], col = "blue", lwd = 2,lty=3)
  lines(time_full[-1], upper_full2[-1], col = "blue", lwd = 2,lty=3)
  lines(time_full, m_full2, col = "blue", lwd = 2)
  
  abline(h = seq(ymin, ymax, by = 0.5), col = "gray", lty = 2)
  legend("bottomleft", legend=c("Gibbs", "Stan"), col=c("red", "blue"), lwd=2, bty="n")
}
