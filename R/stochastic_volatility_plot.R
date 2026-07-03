#' Plot stochastic volatility estimates and forecasts
#'
#' Produces time series plots of posterior stochastic volatility estimates
#' over the estimation sample together with predictive paths over the forecast
#' horizon for a fitted steady-state \code{bvar} object with stochastic volatility.
#'
#' @param x A steady-state \code{bvar} object that has been passed through
#'   \code{\link{fit}} with stochastic volatility enabled.
#' @param ci Numeric. Interval level. Default is \code{0.95}.
#' @param vol Character. Volatility representation: \code{"log_lambda"}
#'   (default) or \code{"sd"}.
#' @param stat Character. Posterior summary statistic to use as point estimate:
#'   \code{"mean"} (default) or \code{"median"}.
#' @param plot_idx Integer vector. Indices of variables to plot. If \code{NULL}
#'   (default), all variables are plotted.
#' @param xlim Numeric vector of length 2. Optional x-axis limits.
#' @param ylim Numeric vector of length 2. Optional y-axis limits.
#'
#' @return An invisible list with:
#' \itemize{
#'   \item \code{estimated}: List of posterior summaries for the estimation sample:
#'     \itemize{
#'       \item \code{point}: \eqn{T \times k} matrix of posterior point estimates
#'         (mean or median depending on \code{stat})
#'       \item \code{lower}: \eqn{T \times k} matrix of lower credible bounds
#'       \item \code{upper}: \eqn{T \times k} matrix of upper credible bounds
#'     }
#'   \item \code{predicted}: List of posterior summaries for the forecast horizon:
#'     \itemize{
#'       \item \code{point}: \eqn{(T+H) \times k} matrix of predictive point estimates
#'         (mean or median depending on \code{stat})
#'       \item \code{lower}: \eqn{(T+H) \times k} matrix of lower bound of prediction interval
#'       \item \code{upper}: \eqn{(T+H) \times k} matrix of upper bound of prediction interval
#'     }
#' }
#'
#' @details
#' The function supports two volatility representations:
#' \itemize{
#'   \item \code{"log_lambda"}: log-volatility states
#'   \item \code{"sd"}: implied standard deviations
#' }
#'
#' For each selected variable, the function plots the posterior point estimate
#' path and credible bands over the estimation sample, and the predictive point
#' estimate path and prediction intervals over the forecast horizon.
#'
#' @export
#' @examples
#' \donttest{
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p = 1, deterministic = "constant")
#' 
#' k <- bvar_obj$setup$k
#' n_free_params_A <- bvar_obj$setup$n_free_params_A
#' 
#' SV_priors <- list(
#' theta_A             =  rep(0, n_free_params_A),
#' Omega_A             =  diag(1000, n_free_params_A),
#' mu_log_lambda_0     =  rep(0, k),
#' sigma2_log_lambda_0 =  rep(1000, k),
#' alpha_phi           =  rep(5, k),
#' beta_phi            = (rep(5, k) - 1) * rep(0.1, k)
#' )
#'
#' bvar_obj <- priors(bvar_obj,
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    SV = TRUE,
#'                    SV_type = "RW",
#'                    SV_priors = SV_priors)
#'
#' bvar_obj <- fit(bvar_obj,
#'                 H = 8,
#'                 d_pred = matrix(rep(1, 8)),
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1,
#'                 verbose = FALSE)
#'
#' stochastic_volatility_plot(bvar_obj, ci = 0.95)
#' }
stochastic_volatility_plot <- function(x, ci = 0.95, vol = "log_lambda",
                                       stat = "mean", plot_idx = NULL,
                                       xlim = NULL, ylim = NULL) {
  
  stat     <- match.arg(stat, c("mean", "median"))
  point_fn <- if (stat == "mean") mean else median
  
  draws           <- rstan::extract(x$fit$stan)
  log_lambda_pred <- draws$log_lambda_pred
  Sigma_u_pred    <- draws$Sigma_u_pred
  
  N     <- dim(x$data)[1]
  p     <- x$setup$p
  k     <- x$setup$k
  N_est <- N - p
  H     <- x$predict$H
  alpha <- 1 - ci
  
  if (vol == "sd") {
    
    sigma_post_point <- matrix(NA, N_est, k)
    sigma_post_lower <- matrix(NA, N_est, k)
    sigma_post_upper <- matrix(NA, N_est, k)
    
    for (t in 1:N_est) {
      Sigma_draws_t <- draws$Sigma_u[, t, , ]
      sigma_draws_t <- apply(Sigma_draws_t, 1, function(S) sqrt(diag(S)))
      
      sigma_post_point[t, ] <- apply(sigma_draws_t, 1, point_fn)
      sigma_post_lower[t, ] <- apply(sigma_draws_t, 1, quantile, probs = alpha / 2)
      sigma_post_upper[t, ] <- apply(sigma_draws_t, 1, quantile, probs = 1 - alpha / 2)
    }
    
  } else if (vol == "log_lambda") {
    
    log_lambda_post_point <- apply(draws$log_lambda, c(2, 3), point_fn)
    log_lambda_post_lower <- apply(draws$log_lambda, c(2, 3), quantile, probs = alpha / 2)
    log_lambda_post_upper <- apply(draws$log_lambda, c(2, 3), quantile, probs = 1 - alpha / 2)
    
  } else {
    stop("vol must be either 'log_lambda' or 'sd'")
  }
  
  if (is.null(plot_idx)) plot_idx <- 1:k
  
  user_xlim <- xlim
  user_ylim <- ylim
  
  est_point <- matrix(NA, N_est, k)
  est_lower <- matrix(NA, N_est, k)
  est_upper <- matrix(NA, N_est, k)
  
  pred_point <- matrix(NA, N_est + H, k)
  pred_lower <- matrix(NA, N_est + H, k)
  pred_upper <- matrix(NA, N_est + H, k)
  
  for (i in plot_idx) {
    
    if (vol == "sd") {
      vol_point <- sigma_post_point[, i]
      vol_lower <- sigma_post_lower[, i]
      vol_upper <- sigma_post_upper[, i]
    } else {
      vol_point <- log_lambda_post_point[, i]
      vol_lower <- log_lambda_post_lower[, i]
      vol_upper <- log_lambda_post_upper[, i]
    }
    
    est_point[, i] <- vol_point
    est_lower[, i] <- vol_lower
    est_upper[, i] <- vol_upper
    
    vol_pred_point <- rep(NA, N_est + H)
    vol_pred_lower <- rep(NA, N_est + H)
    vol_pred_upper <- rep(NA, N_est + H)
    
    vol_pred_point[N_est] <- tail(vol_point, 1)
    vol_pred_lower[N_est] <- tail(vol_lower, 1)
    vol_pred_upper[N_est] <- tail(vol_upper, 1)
    
    for (h in 1:H) {
      
      draws_h <- if (vol == "sd") {
        sqrt(Sigma_u_pred[, h, i, i])
      } else {
        log_lambda_pred[, h, i]
      }
      
      vol_pred_point[N_est + h] <- point_fn(draws_h)
      vol_pred_lower[N_est + h] <- quantile(draws_h, alpha / 2)
      vol_pred_upper[N_est + h] <- quantile(draws_h, 1 - alpha / 2)
    }
    
    pred_point[, i] <- vol_pred_point
    pred_lower[, i] <- vol_pred_lower
    pred_upper[, i] <- vol_pred_upper
    
    xlim_i <- if (is.null(user_xlim)) c(0, N + H + 10) else user_xlim
    
    ylim_i <- if (is.null(user_ylim)) {
      range(c(vol_point, vol_lower, vol_upper,
              vol_pred_lower, vol_pred_upper), na.rm = TRUE)
    } else user_ylim
    
    ts.plot(vol_point,
            col  = "red", lwd = 2,
            main = if (vol == "sd") paste0("sd(u_", i, ")")
            else paste0("ln(lambda_", i, ")"),
            ylim = ylim_i, xlim = xlim_i, ylab = "")
    
    x_insample <- 1:N_est
    polygon(x = c(x_insample, rev(x_insample)),
            y = c(vol_upper, rev(vol_lower)),
            col = rgb(1, 0, 0, 0.15), border = NA)
    
    x_fore <- N_est:(N_est + H)
    polygon(x = c(x_fore, rev(x_fore)),
            y = c(vol_pred_upper[x_fore], rev(vol_pred_lower[x_fore])),
            col = rgb(0, 0, 1, 0.2), border = NA)
    
    lines(x_fore, vol_pred_point[x_fore], col = "blue", lwd = 2)
    
    legend("topleft",
           legend = c("Estimated volatility", "Predicted volatility"),
           col    = c("red", "blue"),
           lwd    = 2,
           bty    = "n")
  }
  
  invisible(list(
    estimated = list(
      point = est_point,
      lower = est_lower,
      upper = est_upper
    ),
    predicted = list(
      point = pred_point,
      lower = pred_lower,
      upper = pred_upper
    )
  ))
}