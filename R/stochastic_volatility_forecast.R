#' Plot stochastic volatility estimates and forecasts
#'
#' Plots the posterior distribution of stochastic volatility over the
#' estimation sample and the forecast horizon for a fitted \code{bvar} object
#' with a stochastic volatility specification. Supports plotting either
#' log-volatility (\code{log_lambda}) or posterior standard deviations
#' (\code{sd}).
#'
#' @param x A \code{bvar} object that has been passed through \code{\link{fit}}
#'   with a stochastic volatility specification.
#' @param ci Numeric. The credible interval width. Default \code{0.95}.
#' @param vol Character. The volatility measure to plot. Either
#'   \code{"log_lambda"} (default) for log-volatility or \code{"sd"} for
#'   posterior standard deviations.
#' @param plot_idx Integer vector. Indices of variables to plot. If \code{NULL}
#'   (default), all variables are plotted.
#' @param xlim Numeric vector of length 2. x-axis limits. If \code{NULL}
#'   (default), limits are set automatically.
#' @param ylim Numeric vector of length 2. y-axis limits. If \code{NULL}
#'   (default), limits are set automatically.
#'
#' @return Invisibly returns \code{NULL}. Called for its plotting side effect.
#' @export
#'
#' @examples
#' \dontrun{
#' model <- bvar(data = my_data)
#' model <- setup(model, p = 2, deterministic = "constant")
#' model <- priors(model)
#' model$SV <- TRUE
#' model$SV_type <- "RW"
#' model <- fit(model)
#' stochastic_volatility_forecast(model, ci = 0.95, vol = "log_lambda")
#' }
stochastic_volatility_forecast <- function(x, ci = 0.95, vol = "log_lambda",
                                           plot_idx = NULL, xlim = NULL, ylim = NULL) {
  
  draws          <- rstan::extract(x$fit$stan)
  log_lambda_pred <- draws$log_lambda_pred
  Sigma_u_pred   <- draws$Sigma_u_pred
  
  N     <- dim(x$data)[1]
  p     <- x$setup$p
  k     <- x$setup$k
  N_est <- N - p
  H     <- x$predict$H
  alpha <- 1 - ci
  
  log_lambda_post_mean  <- apply(draws$log_lambda, c(2, 3), mean)
  log_lambda_post_lower <- apply(draws$log_lambda, c(2, 3), quantile, probs = alpha / 2)
  log_lambda_post_upper <- apply(draws$log_lambda, c(2, 3), quantile, probs = 1 - alpha / 2)
  
  if (vol == "sd") {
    sigma_post_mean  <- matrix(NA, N_est, k)
    sigma_post_lower <- matrix(NA, N_est, k)
    sigma_post_upper <- matrix(NA, N_est, k)
    for (t in 1:N_est) {
      Sigma_draws_t <- draws$Sigma_u[, t, , ]
      sigma_draws_t <- apply(Sigma_draws_t, 1, function(S) sqrt(diag(S)))
      sigma_post_mean[t, ]  <- apply(sigma_draws_t, 1, mean)
      sigma_post_lower[t, ] <- apply(sigma_draws_t, 1, quantile, probs = alpha / 2)
      sigma_post_upper[t, ] <- apply(sigma_draws_t, 1, quantile, probs = 1 - alpha / 2)
    }
  }
  
  if (is.null(plot_idx)) plot_idx <- 1:k
  
  user_xlim <- xlim
  user_ylim <- ylim
  
  for (i in plot_idx) {
    
    if (vol == "sd") {
      vol_mean  <- sigma_post_mean[, i]
      vol_lower <- sigma_post_lower[, i]
      vol_upper <- sigma_post_upper[, i]
    } else {
      vol_mean  <- log_lambda_post_mean[, i]
      vol_lower <- log_lambda_post_lower[, i]
      vol_upper <- log_lambda_post_upper[, i]
    }
    
    vol_pred       <- rep(NA, N_est + H)
    vol_pred_lower <- rep(NA, N_est + H)
    vol_pred_upper <- rep(NA, N_est + H)
    
    vol_pred[N_est]       <- tail(vol_mean,  1)
    vol_pred_lower[N_est] <- tail(vol_lower, 1)
    vol_pred_upper[N_est] <- tail(vol_upper, 1)
    
    for (h in 1:H) {
      draws_h <- if (vol == "sd") sqrt(Sigma_u_pred[, h, i, i]) else log_lambda_pred[, h, i]
      vol_pred[N_est + h]       <- mean(draws_h)
      vol_pred_lower[N_est + h] <- quantile(draws_h, alpha / 2)
      vol_pred_upper[N_est + h] <- quantile(draws_h, 1 - alpha / 2)
    }
    
    xlim_i <- if (is.null(user_xlim)) c(0, N + H + 10) else user_xlim
    ylim_i <- if (is.null(user_ylim)) {
      range(c(vol_mean, vol_lower, vol_upper,
              vol_pred_lower, vol_pred_upper), na.rm = TRUE)
    } else {
      user_ylim
    }
    
    ts.plot(vol_mean,
            col  = "red", lwd = 2,
            main = if (vol == "sd") paste0("sd(u_", i, ")") else paste0("ln(lambda_", i, ")"),
            ylim = ylim_i, xlim = xlim_i, ylab = "")
    
    x_insample <- 1:N_est
    polygon(x   = c(x_insample, rev(x_insample)),
            y   = c(vol_upper, rev(vol_lower)),
            col = rgb(1, 0, 0, 0.15), border = NA)
    
    x_fore <- N_est:(N_est + H)
    polygon(x   = c(x_fore, rev(x_fore)),
            y   = c(vol_pred_upper[x_fore], rev(vol_pred_lower[x_fore])),
            col = rgb(0, 0, 1, 0.2), border = NA)
    lines(x_fore, vol_pred[x_fore], col = "blue", lwd = 2)
  }
  
  invisible(NULL)
}