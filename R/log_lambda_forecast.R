log_lambda_forecast <- function(x, ci=0.95, plot_idx=NULL, xlim=NULL, ylim=NULL) {
  draws <- rstan::extract(x$fit$stan)
  log_lambda_array <- x$fit$stanf$log_lambda_array
  N <- dim(bvar_obj$data)[1]
  p <- bvar_obj$setup$p
  k <- bvar_obj$setup$k
  N_est <- N-p
  H <- bvar_obj$predict$H
  alpha = 1-ci
  
  log_lambda_post_mean <- apply(draws$log_lambda, c(2,3), mean)
  if (is.null(plot_idx)) {
    plot_idx <- k
  }
  
  for (i in plot_idx) {
    
    log_lambda_est <- log_lambda_post_mean[, i]
    log_lambda_pred       <- rep(NA, N_est + H)
    log_lambda_pred_lower <- rep(NA, N_est + H)
    log_lambda_pred_upper <- rep(NA, N_est + H)
    log_lambda_pred[N_est] <- tail(log_lambda_est, 1)
    
    for (h in 1:H) {
      draws_h <- log_lambda_array[, h, i]
      log_lambda_pred[N_est+h]       <- mean(draws_h)
      log_lambda_pred_lower[N_est+h] <- quantile(draws_h, alpha/2)
      log_lambda_pred_upper[N_est+h] <- quantile(draws_h, 1-alpha/2)
    }
    
    
    if (is.null(xlim)) {
      xlim <- c(0, N + H + 10)
    }
    if (is.null(ylim)) {
      ylim <- range(c(log_lambda_est,
                      log_lambda_pred_lower,
                      log_lambda_pred_upper),
                    na.rm = TRUE)
    }
    ts.plot(log_lambda_est,
            col = "red",
            lwd = 2,
            main = paste0("ln(lambda_", i, ")"),
            ylim = ylim,
            xlim = xlim)
    
    lines(log_lambda_pred, col="blue", lwd=2)
    lines(log_lambda_pred_lower, col="blue", lty=2, lwd=2)
    lines(log_lambda_pred_upper, col="blue", lty=2, lwd=2)
    
  }
}