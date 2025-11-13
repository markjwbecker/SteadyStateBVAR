forecast2 <- function(x, ci=0.95, fcst_type=c("mean", "median"), plot_idx=NULL, xlim, ylim){
  
  Y_pred <- bvar_obj$fit$gibbs$fcst_draws
  Y_pred_m <- apply(Y_pred, c(1, 2), fcst_type)
  alpha <- 1 - ci
  Y_pred_lower <- apply(Y_pred, c(1, 2), quantile, probs = alpha/2)
  Y_pred_upper <- apply(Y_pred, c(1, 2), quantile, probs = 1 - alpha/2)
  
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
  
  lines(time_full[-1], lower_full[-1], col = "blue", lwd = 2,lty=3)
  lines(time_full[-1], upper_full[-1], col = "blue", lwd = 2,lty=3)
  lines(time_full, m_full, col = "blue", lwd = 2)
  abline(h = seq(ymin, ymax, by = 0.5), col = "gray", lty = 2)
}