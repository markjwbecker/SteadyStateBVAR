plot_forecast <- function(fit, Y) {
  posterior <- rstan::extract(fit)
  Y_pred <- posterior$Y_pred
  Y_pred_mean <- apply(Y_pred, c(2, 3), mean)
  Y_pred_lower <- apply(Y_pred, c(2, 3), quantile, probs = 0.025)
  Y_pred_upper <- apply(Y_pred, c(2, 3), quantile, probs = 0.975)
  
  T <- nrow(Y)
  H <- nrow(Y_pred_mean)
  m <- ncol(Y)
  
  time_hist <- time(Y)
  time_fore <- seq(tail(time_hist, 1) + 1/4, by = 1/4, length.out = H)
  
  par(mfrow = c(3, 1))
  
  for (i in 1:ncol(Y)) {
    fcsty <- c(rep(NA, T - 1), Y[T,i], Y_pred_mean[, i])
    fcstl <- c(rep(NA, T - 1), Y[T,i], Y_pred_lower[, i])
    fcstu <- c(rep(NA, T - 1), Y[T,i], Y_pred_upper[, i])
    smply <- c(Y[,i], rep(NA,length(Y_pred_mean[, i])))
    min.y <- min(na.omit(c(fcsty, fcstl, fcstu, smply)))
    max.y <- max(na.omit(c(fcsty, fcstl, fcstu, smply)))
    ylim <- c(min.y, max.y)
    fcsty <- ts(fcsty, start = time_hist[1], frequency = 4)
    smply <- ts(smply, start = time_hist[1], frequency = 4)
    fcstl <- ts(fcstl, start = time_hist[1], frequency = 4)
    fcstu <- ts(fcstu, start = time_hist[1], frequency = 4)
    plot.ts(fcsty,main=colnames(Y)[i],xlab="Time", ylim=ylim, col = "blue", lwd = 2,ylab=NULL)
    lines(smply, col = "black", lwd = 2)
    lines(fcstl, col = "blue", lty = 3, lwd=2)
    lines(fcstu, col = "blue", lty = 3, lwd=2)
  }
}
