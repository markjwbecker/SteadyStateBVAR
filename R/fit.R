#' Estimate the steady-state BVAR model using Stan
#'
#' Estimates the steady-state BVAR model using the  No-U-Turn sampler (a variant of Hamiltonian Monte Carlo) via Stan.
#' Uses the data, setup, and priors stored in the steady-state \code{bvar} object. Supports both
#' homoscedastic and stochastic volatility (RW or AR1) specifications.
#'
#'
#' @param x A steady-state \code{bvar} object that has been passed through
#'   \code{\link{setup}} and \code{\link{priors}}.
#' @param H Positive Integer. Forecast horizon.
#'   Default is \code{1}.
#' @param d_pred Matrix of size H x q. Future values of the deterministic variables d_t. Default is \code{NULL}, must be provided by the user.
#' @param iter Integer. Total number of MCMC iterations per chain. Default is \code{2000}.
#' @param warmup Integer. Number of warmup (burn-in) iterations per chain.
#'   Default is \code{floor(iter/2)}.
#' @param chains Integer. Number of MCMC chains. Default is \code{1}.
#' @param cores Positive integer specifying the number of CPU cores used for
#' sampling. Default is \code{getOption("mc.cores", 1L)}, i.e. the \code{mc.cores}
#' option (if it has been set), otherwise defaults to \code{1} core.
#' @param verbose Logical indicating whether to print intermediate output from Stan on the console,
#' defaults to \code{FALSE}.
#' @param auto_write Logical indicating whether Stan models should be
#' automatically written to the disk cache via \code{rstan}. Default is \code{FALSE}.
#'
#' @return A fitted steady-state \code{bvar} object with:
#' \itemize{
#'   \item \code{fit$stan}: An object of class \code{stanfit}
#'   \item \code{fit$posterior_means}: List of posterior mean estimates
#'   \item \code{fit$posterior_medians}: List of posterior median estimates
#' }
#'
#' @details
#' The function selects the appropriate Stan model based on prior settings:
#' \itemize{
#'   \item Homoscedastic with Jeffreys prior:
#'     \code{steady_state_bvar_homoscedastic_jeffreys_prior.stan}
#'   \item Homoscedastic with uninformative inverse-Wishart prior:
#'     \code{steady_state_bvar_homoscedastic_inverse_wishart_prior.stan}
#'   \item Random Walk stochastic volatility:
#'     \code{steady_state_bvar_RW_stochastic_volatility.stan}
#'   \item AR1 stochastic volatility:
#'     \code{steady_state_bvar_AR1_stochastic_volatility.stan}
#' }
#' The function estimates the following parameters (see \link{bvar} for details):
#' \itemize{
#'       \item \code{beta}: kp×k VAR coefficient matrix
#'       \item \code{Psi}: k×q steady-state parameter matrix
#'       \item \code{Sigma_u}: innovation covariance matrix (k×k for homoscedastic,
#'         T×k×k for stochastic volatility)
#'       \item If Random Walk stochastic volatility: 
#'        \itemize{
#'        \item \code{A}: k×k lower triangular matrix with ones on the diagonal that describes
#'         the contemporaneous interaction of the endogenous variables
#'        \item \code{phi}: k-dimensional vector of log volatility innovation variances 
#'        }
#'       \item If AR1 stochastic volatility:
#'        \itemize{
#'        \item \code{A}: k×k lower triangular matrix with ones on the diagonal that describes
#'         the contemporaneous interaction of the endogenous variables
#'         \item \code{gamma_0}: k-dimensional vector of log volatility intercepts
#'         \item \code{gamma_1}: k-dimensional vector of log volatility slopes
#'         \item \code{Phi}: k×k log volatility innovation covariance matrix
#'        }
#'     }
#'
#'
#' @export
#' @examples
#' \dontrun{
#' #homoscedastic with Jeffreys prior
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    Jeffrey = TRUE,
#'                    SV = FALSE,
#'                    SV_type = NULL,
#'                    SV_priors = NULL)
#'                    
#' bvar_obj <- fit(bvar_obj,
#'                 H = 8,
#'                 d_pred = matrix(rep(1,8)),
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1,
#'                 verbose = FALSE,
#'                 auto_write = FALSE)
#'                    
#' #RW stochastic volatility
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")
#' 
#' k <- bvar_obj$setup$k
#' n_free_params_A <- bvar_obj$setup$n_free_params_A
#' 
#' SV_priors_RW <- list(
#' theta_A             =  rep(0, n_free_params_A),
#' Omega_A             =  diag(1000, n_free_params_A),
#' mu_log_lambda_0     =  rep(0, k),
#' sigma2_log_lambda_0 =  rep(1000, k),
#' alpha_phi           =  rep(5, k),
#' beta_phi            = (rep(5, k) - 1) * rep(0.1, k)
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
#'                 verbose = FALSE,
#'                 auto_write = FALSE)
#'                    
#' #AR1 stochastic volatility
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")
#' 
#' k <- bvar_obj$setup$k
#' n_free_params_A <- bvar_obj$setup$n_free_params_A
#' 
#' SV_priors_AR <- list(
#' theta_A            =  rep(0, n_free_params_A),
#' Omega_A            =  diag(1000, n_free_params_A),
#' theta_gamma_0      =  rep(0.1, k),
#' Omega_gamma_0      =  diag(1000, k),
#' theta_gamma_1      =  rep(0.9, k),
#' Omega_gamma_1      =  diag(10, k),
#' theta_log_lambda_0 =  rep(0.1, k)/(1-rep(0.9, k)),
#' Omega_log_lambda_0 =  diag(1000, k),
#' V_0                = (10 - k - 1) * diag(k),
#' m_0                =  10
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
#'                    SV_priors = SV_priors_AR)
#'
#'
#' bvar_obj <- fit(bvar_obj,
#'                 H = 8,
#'                 d_pred = matrix(rep(1,8)),
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1,
#'                 verbose = FALSE,
#'                 auto_write = FALSE)
#' }
fit <- function(x,
                H = 1,
                d_pred = NULL,
                iter = 2000,
                warmup = floor(iter/2),
                chains = 1,
                cores = getOption("mc.cores", 1L),
                verbose = FALSE,
                auto_write = FALSE) {
  
  if (!inherits(x, "bvar")) stop("must be a 'bvar' object")
  if (is.null(x$setup)) stop("must be passed through setup")
  if (is.null(x$priors)) stop("must be passed through priors")
  
  if (is.null(d_pred)) {
    stop("d_pred must be supplied")
  }
  
  if (!is.logical(verbose) ||
      length(verbose) != 1 ||
      is.na(verbose)) {
    stop("verbose must be TRUE or FALSE")
  }
  
  if (!is.numeric(H) ||
      length(H) != 1 ||
      !is.finite(H) ||
      H < 1 ||
      H != floor(H)) {
    stop("H must be a positive integer")
  }
  
  if (!is.numeric(iter) ||
      length(iter) != 1 ||
      !is.finite(iter) ||
      iter < 1 ||
      iter != floor(iter)) {
    stop("iter must be a positive integer")
  }
  
  if (!is.numeric(warmup) ||
      length(warmup) != 1 ||
      !is.finite(warmup) ||
      warmup < 0 ||
      warmup != floor(warmup)) {
    stop("warmup must be a non-negative integer")
  }
  
  if (!is.numeric(chains) ||
      length(chains) != 1 ||
      !is.finite(chains) ||
      chains < 1 ||
      chains != floor(chains)) {
    stop("chains must be a positive integer")
  }
  
  if (!is.numeric(cores) ||
      length(cores) != 1 ||
      !is.finite(cores) ||
      cores < 1 ||
      cores != floor(cores)) {
    stop("cores must be a positive integer")
  }
  
  if (!is.logical(auto_write) ||
      length(auto_write) != 1 ||
      is.na(auto_write)) {
    stop("auto_write must be TRUE or FALSE")
  }
  
  if (!is.matrix(d_pred)) {
    stop("d_pred must be a matrix")
  }
  
  if (nrow(d_pred) != H) {
    stop("nrow(d_pred) must equal H")
  }
  
  if (!is.null(x$setup$q) &&
      ncol(d_pred) != x$setup$q) {
    stop("ncol(d_pred) must equal q")
  }
  
  x$predict$H <- H
  x$predict$d_pred <- d_pred
  
  setup <- x$setup
  priors <- x$priors
  
  SV <- isTRUE(priors$SV)
  SV_type <- priors$SV_type
  
  stan_data <- c(
    setup,
    list(H = x$predict$H, d_pred = x$predict$d_pred),
    priors
  )
  
  stan_file <- if (isFALSE(priors$Jeffrey)) {
    system.file("steady_state_bvar_homoscedastic_inverse_wishart_prior.stan", package = "SteadyStateBVAR")
  } else {
    system.file("steady_state_bvar_homoscedastic_jeffreys_prior.stan", package = "SteadyStateBVAR")
  }
  
  if (isTRUE(SV)) {
    
    if (is.null(SV_type)) stop("SV_type missing in priors")
    
    if (SV_type == "RW") {
      stan_file <- system.file(
        "steady_state_bvar_RW_stochastic_volatility.stan",
        package = "SteadyStateBVAR"
      )
    } else if (SV_type == "AR1") {
      stan_file <- system.file(
        "steady_state_bvar_AR1_stochastic_volatility.stan",
        package = "SteadyStateBVAR"
      )
    } else {
      stop("SV_type must be 'RW' or 'AR1'")
    }
    
    stan_data <- c(priors$SV_priors, stan_data)
    
    if (setup$k == 2 && !is.null(priors$SV_priors$theta_A)) {
      stan_data$theta_A <- as.array(priors$SV_priors$theta_A[1])
    }
  }
  
  if (!nzchar(stan_file)) {
    stop("Stan model file not found in installed package.")
  }
  
  old_mc_cores <- getOption("mc.cores")
  options(mc.cores = cores)
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  
  old_auto_write <- rstan::rstan_options("auto_write")
  on.exit(rstan::rstan_options(auto_write = old_auto_write), add = TRUE)
  rstan::rstan_options(auto_write = auto_write)
  
  x$fit$stan <- rstan::stan(
    file = stan_file,
    data = stan_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    verbose = verbose
  )
  
  posterior <- rstan::extract(x$fit$stan)
  
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