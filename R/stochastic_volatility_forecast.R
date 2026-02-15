stochastic_volatility_forecast <- function(x, ci=0.95, plot_idx=NULL, xlim=NULL, ylim=NULL,vol="log_lambda") {
  draws <- rstan::extract(x$fit$stan)
  log_lambda_array <- x$fit$stanf$log_lambda_array
  Sigma_u_array <- x$fit$stanf$Sigma_u_array
  N <- dim(bvar_obj$data)[1]
  p <- bvar_obj$setup$p
  k <- bvar_obj$setup$k
  N_est <- N-p
  H <- bvar_obj$predict$H
  alpha = 1-ci
  
  log_lambda_post_mean <- apply(draws$log_lambda, c(2,3), mean)
  if (vol == "sd") {
    sigma_posterior_mean <- matrix(NA, N_est, k)
    for(t in 1:N_est){
      Sigma_mean <- apply(draws$Sigma_u[,t,,], c(2,3), mean)
      sigma_posterior_mean[t,] <- sqrt(diag(Sigma_mean))
    }
  }
  
  
  if (is.null(plot_idx)) {
    plot_idx <- k
  }
  
  for (i in plot_idx) {
    
    if (vol == "sd") {
      vol_est <- sigma_posterior_mean[, i]
    } else {
      vol_est <- log_lambda_post_mean[, i]
    }
    vol_pred       <- rep(NA, N_est + H)
    vol_pred_lower <- rep(NA, N_est + H)
    vol_pred_upper <- rep(NA, N_est + H)
    vol_pred[N_est] <- tail(vol_est, 1)
    vol_pred_lower[N_est] <- tail(vol_est, 1)
    vol_pred_upper[N_est] <- tail(vol_est, 1)
    
    for (h in 1:H) {
      if (vol == "sd") {
        sigma2_draws <- Sigma_u_array[, h, i, i]
        draws_h <- sqrt(sigma2_draws)
      } else {
        draws_h <- log_lambda_array[, h, i]
      }
      vol_pred[N_est+h]       <- mean(draws_h)
      vol_pred_lower[N_est+h] <- quantile(draws_h, alpha/2)
      vol_pred_upper[N_est+h] <- quantile(draws_h, 1-alpha/2)
    }
    
    
    if (is.null(xlim)) {
      xlim <- c(0, N + H + 10)
    }
    if (is.null(ylim)) {
      ylim <- range(c(vol_est,
                      vol_pred_lower,
                      vol_pred_upper),
                    na.rm = TRUE)
    }
    ts.plot(vol_est,
            col = "red",
            lwd = 2,
            main = if (vol == "sd") {
              paste0("sd(u_", i, ",t)")
            } else {
              paste0("ln(lambda_", i, ")")
            },
            ylim = ylim,
            xlim = xlim)
    
    x_fore <- (N_est):(N_est + H)
    mean_full  <- vol_pred[x_fore]
    lower_full <- vol_pred_lower[x_fore]
    upper_full <- vol_pred_upper[x_fore]
    polygon(x = c(x_fore, rev(x_fore)),
            y = c(upper_full, rev(lower_full)),
            col = rgb(0, 0, 1, 0.2), border = NA)
    lines(x_fore, mean_full, col = "blue", lwd = 2)
    
  }
}