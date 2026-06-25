#' Estimate a Bayesian VAR model using Stan
#'
#' Runs Stan to estimate the BVAR model using data, setup, and priors stored
#' in a \code{bvar} object. Supports standard homoscedastic steady-state BVAR
#' as well as stochastic volatility specifications (RW or AR(1)).
#'
#' Forecast inputs must be provided in \code{x$predict$H} and
#' \code{x$predict$d_pred} before calling \code{fit()}.
#'
#' @param x A \code{bvar} object that has been passed through \code{\link{setup}} and
#'   \code{\link{priors}}.
#' @param iter Integer. Total number of MCMC iterations per chain. Default \code{5000}.
#' @param warmup Integer. Number of warmup (burn-in) iterations per chain. Default \code{2500}.
#' @param chains Integer. Number of MCMC chains. Default \code{2}.
#' @param cores Integer. Number of CPU cores used for sampling.
#' @param auto_write Logical. Whether to enable \code{rstan} auto-write behavior.
#'
#' @return A \code{bvar} object containing:
#'   \item{fit$stan}{Stan fit object}
#'   \item{posterior_means}{Posterior means of key parameters}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' yt <- matrix(rnorm(40, 0, 1), 20, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1)
#'
#' bvar_obj <- priors(bvar_obj,
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2))
#'
#' bvar_obj$predict$H <- 1
#' bvar_obj$predict$d_pred <- matrix(1)
#'
#' bvar_obj <- fit(bvar_obj,
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1,
#'                 auto_write = FALSE)
#'}
fit <- function(x,
                iter = 5000,
                warmup = 2500,
                chains = 2,
                cores = min(chains, parallel::detectCores()),
                auto_write = TRUE) {
  
  if (!inherits(x, "bvar")) stop("must be a 'bvar' object")
  if (is.null(x$setup)) stop("must be passed through setup")
  if (is.null(x$priors)) stop("must be passed through priors")
  
  if (is.null(x$predict) ||
      is.null(x$predict$H) ||
      is.null(x$predict$d_pred)) {
    stop("predict must contain H and d_pred")
  }
  
  if (!is.numeric(x$predict$H) ||
      length(x$predict$H) != 1 ||
      x$predict$H < 1) {
    stop("H must be a positive integer")
  }
  
  if (!is.matrix(x$predict$d_pred)) {
    stop("d_pred must be a matrix")
  }
  
  if (nrow(x$predict$d_pred) != x$predict$H) {
    stop("nrow(d_pred) must equal H")
  }
  
  if (!is.null(x$setup$q) &&
      ncol(x$predict$d_pred) != x$setup$q) {
    stop("ncol(d_pred) must equal q")
  }
  
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
    system.file("inv_wishart_cov.stan", package = "SteadyStateBVAR")
  } else {
    system.file("diffuse_cov.stan", package = "SteadyStateBVAR")
  }
  
  if (isTRUE(SV)) {
    
    if (is.null(SV_type)) stop("SV_type missing in priors")
    
    if (SV_type == "RW") {
      stan_file <- system.file(
        "stochastic_volatility_RW.stan",
        package = "SteadyStateBVAR"
      )
    } else if (SV_type == "AR") {
      stan_file <- system.file(
        "stochastic_volatility_stationaryAR.stan",
        package = "SteadyStateBVAR"
      )
    } else {
      stop("SV_type must be 'RW' or 'AR'")
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
  
  if (isTRUE(auto_write)) {
    rstan::rstan_options(auto_write = TRUE)
  }
  
  x$fit$stan <- rstan::stan(
    file = stan_file,
    data = stan_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    verbose = FALSE
  )
  
  posterior <- rstan::extract(x$fit$stan)
  
  x$posterior_means <- list()
  
  if (!SV) {
    
    x$posterior_means$beta <- apply(posterior$beta, c(2, 3), mean)
    x$posterior_means$Psi <- apply(posterior$Psi, c(2, 3), mean)
    x$posterior_means$Sigma_u <- apply(posterior$Sigma_u, c(2, 3), mean)
    
  } else if (SV_type == "AR") {
    
    x$posterior_means$beta <- apply(posterior$beta, c(2, 3), mean)
    x$posterior_means$Psi <- apply(posterior$Psi, c(2, 3), mean)
    x$posterior_means$A <- apply(posterior$A, c(2, 3), mean)
    x$posterior_means$gamma_0 <- apply(posterior$gamma_0, 2, mean)
    x$posterior_means$gamma_1 <- apply(posterior$gamma_1, 2, mean)
    x$posterior_means$Phi <- apply(posterior$Phi, c(2, 3), mean)
    
  } else if (SV_type == "RW") {
    
    x$posterior_means$beta <- apply(posterior$beta, c(2, 3), mean)
    x$posterior_means$Psi <- apply(posterior$Psi, c(2, 3), mean)
    x$posterior_means$A <- apply(posterior$A, c(2, 3), mean)
    x$posterior_means$phi <- apply(posterior$phi, 2, mean)
  }
  
  return(x)
}