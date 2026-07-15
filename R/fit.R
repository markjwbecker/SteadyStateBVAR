#' Estimate the steady-state BVAR model using Stan
#'
#' Estimates the steady-state BVAR model using the No-U-Turn sampler (a variant of Hamiltonian Monte Carlo) via Stan.
#' Also generates draws from the joint predictive distribution.
#' Uses the data, setup, and priors stored in the steady-state \code{bvar} object.
#'
#' @param x A steady-state \code{bvar} object that has been passed through
#'   \code{\link{setup}} and \code{\link{priors}}.
#' @param H Positive Integer. Forecast horizon.
#'   Default is \code{1}.
#' @param d_pred Matrix of size \eqn{H \times q}. Future values of the deterministic variables \eqn{d_t}. Default is \code{NULL}, must be provided by the user.
#' @param ... Additional arguments passed directly to the \code{rstan} function \code{\link[rstan]{sampling}}
#' (e.g. \code{iter}, \code{warmup}, \code{chains}, \code{cores}, \code{control},
#' \code{seed}, \code{init}, \code{thin}, \code{algorithm}, \code{pars},
#' \code{include}, \code{refresh}, \code{verbose}, \code{save_warmup},
#' \code{sample_file}, \code{diagnostic_file}). If \code{pars}/\code{include} is
#' used to exclude model parameters required by \code{fit()} for posterior
#' summaries, an error will be raised.
#'
#' @return A fitted steady-state \code{bvar} object with:
#' \itemize{
#' \item \code{x$fit$stan}: An object of class \code{stanfit}
#' \item \code{x$fit$posterior_means}: List of posterior mean estimates
#' \item \code{x$fit$posterior_medians}: List of posterior median estimates
#' }
#'
#' @details
#' The function selects the appropriate precompiled Stan model based on settings from \link{priors}:
#' \itemize{
#'   \item \code{steady_state_bvar_homoscedastic_jeffreys_prior}
#'   \item \code{steady_state_bvar_homoscedastic_inverse_wishart_prior}
#'   \item \code{steady_state_bvar_RW_stochastic_volatility}
#'   \item \code{steady_state_bvar_AR1_stochastic_volatility}
#' }
#' The function estimates the following parameters (see \link{bvar} for details):
#' \itemize{
#'       \item \code{beta}: \eqn{kp \times k} VAR coefficient matrix
#'       \item \code{Psi}: \eqn{k \times q} steady-state parameter matrix
#'       \item \code{Sigma_u}: innovation covariance matrix (\eqn{k \times k} for homoscedastic,
#'         \eqn{T \times k \times k} for stochastic volatility)
#'       \item If Random Walk stochastic volatility:
#'        \itemize{
#'        \item \code{A}: \eqn{k \times k} lower triangular matrix with ones on the diagonal that describes
#'         the contemporaneous interaction of the endogenous variables
#'        \item \code{phi}: \eqn{k}-dimensional vector of log volatility innovation variances
#'        }
#'       \item If AR1 stochastic volatility:
#'        \itemize{
#'        \item \code{A}: \eqn{k \times k} lower triangular matrix with ones on the diagonal that describes
#'         the contemporaneous interaction of the endogenous variables
#'         \item \code{gamma_0}: \eqn{k}-dimensional vector of log volatility intercepts
#'         \item \code{gamma_1}: \eqn{k}-dimensional vector of log volatility slopes
#'         \item \code{Phi}: \eqn{k \times k} log volatility innovation covariance matrix
#'        }
#'     }
#'
#' @export
#' @examples
#' \donttest{
#' #homoscedastic with Jeffreys prior, d_t = constant
#' 
#' yt <- matrix(rnorm(20), 10, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p = 1, deterministic = "constant")
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    Jeffreys = TRUE)
#'                    
#' H <- 8
#' d_pred <- matrix(rep(1,8))
#' colnames(d_pred) <- c("constant")
#' rownames(d_pred) <- paste("Horizon", 1:H)
#' print(d_pred)
#' 
#' bvar_obj <- fit(bvar_obj,
#'                 H = H,
#'                 d_pred = d_pred,
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1)
#'                 
#'                 
#' #homoscedastic with inverse-Wishart prior, d_t = constant and dummy
#' 
#' yt <- matrix(rnorm(20), 10, 2)
#'
#' bvar_obj <- bvar(data = yt)
#' 
#' dummy_variable <- c(rep(1,5), rep(0,5))
#'
#' bvar_obj <- setup(bvar_obj, p = 1,
#'                   deterministic = "constant_and_dummy",
#'                   dummy = dummy_variable)
#'                   
#' k <- bvar_obj$setup$k
#' q <- bvar_obj$setup$q
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, k*q),
#'                    Omega_Psi = diag(0.1, k*q, k*q),
#'                    Jeffreys = FALSE) #inverse-Wishart
#'
#' H <- 8
#' d_pred <- cbind(rep(1, H), rep(0, H))
#' colnames(d_pred) <- c("constant", "dummy")
#' rownames(d_pred) <- paste("Horizon", 1:H)
#' print(d_pred)
#'
#' bvar_obj <- fit(bvar_obj,
#'                 H = H,
#'                 d_pred = d_pred,
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1)
#'  
#'
#' #RW stochastic volatility
#' 
#' yt <- matrix(rnorm(20), 10, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")
#'
#' k <- bvar_obj$setup$k
#' n_free_params_A <- bvar_obj$setup$n_free_params_A
#'
#' SV_priors_RW <- list(
#' theta_A              =  rep(0, n_free_params_A),
#' Omega_A              =  diag(1000, n_free_params_A),
#' mu_log_lambda_1      =  rep(0, k),
#' sigma2_log_lambda_1  =  rep(1000, k),
#' alpha_phi            =  rep(5, k),
#' beta_phi             = (rep(5, k) - 1) * rep(0.1, k)
#' )
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    SV = TRUE,
#'                    SV_type = "RW",
#'                    SV_priors = SV_priors_RW)
#'
#' bvar_obj <- fit(bvar_obj,
#'                 H = 8,
#'                 d_pred = matrix(rep(1,8)),
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1,
#'                 control = list(max_treedepth = 12, adapt_delta = 0.85)
#'                 )
#'
#'
#' #AR1 stochastic volatility
#' 
#' yt <- matrix(rnorm(20), 10, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")
#'
#' k <- bvar_obj$setup$k
#' n_free_params_A <- bvar_obj$setup$n_free_params_A
#'
#' SV_priors_AR1 <- list(
#' theta_A               =  rep(0, n_free_params_A),
#' Omega_A               =  diag(1000, n_free_params_A),
#' theta_gamma_0         =  rep(0.1, k),
#' Omega_gamma_0         =  diag(1000, k),
#' theta_gamma_1         =  rep(0.9, k),
#' Omega_gamma_1         =  diag(10, k),
#' theta_log_lambda_1    =  rep(0.1, k)/(1-rep(0.9, k)),
#' Omega_log_lambda_1    =  diag(1000, k),
#' V_Phi                 = (10 - k - 1) * diag(k),
#' m_Phi                 =  10
#' )
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    SV = TRUE,
#'                    SV_type = "AR1",
#'                    SV_priors = SV_priors_AR1)
#'
#' bvar_obj <- fit(bvar_obj,
#'                 H = 8,
#'                 d_pred = matrix(rep(1,8)),
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1,
#'                 control = list(max_treedepth = 12, adapt_delta = 0.85)
#'                 )
#' }
fit <- function(x, H = 1, d_pred = NULL, ...) {
  
  if (!inherits(x, "bvar")) stop("must be a 'bvar' object")
  if (is.null(x$setup)) stop("must be passed through setup")
  if (is.null(x$priors)) stop("must be passed through priors")
  if (is.null(d_pred)) stop("d_pred must be supplied")
  
  if (!is.numeric(H) || length(H) != 1 || !is.finite(H) || H < 1 || H != floor(H)) {
    stop("H must be a positive integer")
  }
  if (!is.matrix(d_pred)) stop("d_pred must be a matrix")
  if (nrow(d_pred) != H) stop("nrow(d_pred) must equal H")
  if (!is.null(x$setup$q) && ncol(d_pred) != x$setup$q) {
    stop("ncol(d_pred) must equal q")
  }
  
  x$predict$H <- H
  x$predict$d_pred <- d_pred
  
  setup <- x$setup
  priors <- x$priors
  SV <- isTRUE(priors$SV)
  SV_type <- priors$SV_type
  
  stan_data <- c(setup, list(H = H, d_pred = d_pred), priors)
  
  model_name <- if (isFALSE(priors$Jeffreys)) {
    "steady_state_bvar_homoscedastic_inverse_wishart_prior"
  } else {
    "steady_state_bvar_homoscedastic_jeffreys_prior"
  }
  
  if (isTRUE(SV)) {
    if (is.null(SV_type)) stop("SV_type missing in priors")
    model_name <- if (SV_type == "RW") {
      "steady_state_bvar_RW_stochastic_volatility"
    } else if (SV_type == "AR1") {
      "steady_state_bvar_AR1_stochastic_volatility"
    } else {
      stop("SV_type must be 'RW' or 'AR1'")
    }
    stan_data <- c(priors$SV_priors, stan_data)
    if (setup$k == 2 && !is.null(priors$SV_priors$theta_A)) {
      stan_data$theta_A <- as.array(priors$SV_priors$theta_A[1])
    }
  }
  
  x$fit$stan <- rstan::sampling(stanmodels[[model_name]], data = stan_data, ...)
  
  posterior <- rstan::extract(x$fit$stan)
  
  required_params <- if (!SV) {
    c("beta", "Psi", "Sigma_u")
  } else if (SV_type == "AR1") {
    c("beta", "Psi", "A", "gamma_0", "gamma_1", "Phi", "Sigma_u")
  } else {
    c("beta", "Psi", "A", "phi", "Sigma_u")
  }
  
  missing_params <- setdiff(required_params, names(posterior))
  if (length(missing_params) > 0) {
    stop(sprintf(
      "The following required parameters were not found in the Stan output: %s. This likely happened because 'pars'/'include' was passed via ... and excluded them.",
      paste(missing_params, collapse = ", ")
    ))
  }
  
  x$fit$posterior_means <- list()
  x$fit$posterior_medians <- list()
  
  if (!SV) {
    
    x$fit$posterior_means$beta <- apply(posterior$beta, c(2, 3), mean)
    x$fit$posterior_means$Psi <- apply(posterior$Psi, c(2, 3), mean)
    x$fit$posterior_means$Sigma_u <- apply(posterior$Sigma_u, c(2, 3), mean)
    
    x$fit$posterior_medians$beta <- apply(posterior$beta, c(2, 3), median)
    x$fit$posterior_medians$Psi <- apply(posterior$Psi, c(2, 3), median)
    x$fit$posterior_medians$Sigma_u <- apply(posterior$Sigma_u, c(2, 3), median)
    
  } else if (SV_type == "AR1") {
    
    x$fit$posterior_means$beta <- apply(posterior$beta, c(2, 3), mean)
    x$fit$posterior_means$Psi <- apply(posterior$Psi, c(2, 3), mean)
    x$fit$posterior_means$A <- apply(posterior$A, c(2, 3), mean)
    x$fit$posterior_means$gamma_0 <- apply(posterior$gamma_0, 2, mean)
    x$fit$posterior_means$gamma_1 <- apply(posterior$gamma_1, 2, mean)
    x$fit$posterior_means$Phi <- apply(posterior$Phi, c(2, 3), mean)
    x$fit$posterior_means$Sigma_u <- apply(posterior$Sigma_u, c(2, 3, 4), mean)
    
    x$fit$posterior_medians$beta <- apply(posterior$beta, c(2, 3), median)
    x$fit$posterior_medians$Psi <- apply(posterior$Psi, c(2, 3), median)
    x$fit$posterior_medians$A <- apply(posterior$A, c(2, 3), median)
    x$fit$posterior_medians$gamma_0 <- apply(posterior$gamma_0, 2, median)
    x$fit$posterior_medians$gamma_1 <- apply(posterior$gamma_1, 2, median)
    x$fit$posterior_medians$Phi <- apply(posterior$Phi, c(2, 3), median)
    x$fit$posterior_medians$Sigma_u <- apply(posterior$Sigma_u, c(2, 3, 4), median)
    
  } else if (SV_type == "RW") {
    
    x$fit$posterior_means$beta <- apply(posterior$beta, c(2, 3), mean)
    x$fit$posterior_means$Psi <- apply(posterior$Psi, c(2, 3), mean)
    x$fit$posterior_means$A <- apply(posterior$A, c(2, 3), mean)
    x$fit$posterior_means$phi <- apply(posterior$phi, 2, mean)
    x$fit$posterior_means$Sigma_u <- apply(posterior$Sigma_u, c(2, 3, 4), mean)
    
    x$fit$posterior_medians$beta <- apply(posterior$beta, c(2, 3), median)
    x$fit$posterior_medians$Psi <- apply(posterior$Psi, c(2, 3), median)
    x$fit$posterior_medians$A <- apply(posterior$A, c(2, 3), median)
    x$fit$posterior_medians$phi <- apply(posterior$phi, 2, median)
    x$fit$posterior_medians$Sigma_u <- apply(posterior$Sigma_u, c(2, 3, 4), median)
  }
  
  return(x)
}