#' Estimate a BVAR model using Stan
#'
#' Runs Stan to estimate the BVAR model using the data, setup, and priors
#' stored in the \code{bvar} object. Supports the standard homoscedastic steady state BVAR, but also stochastic
#' volatility specifications, either Random Walk or AR(1).
#'
#' @param x A \code{bvar} object that has been passed through \code{\link{setup}}
#'   and \code{\link{priors}}.
#' @param iter Integer. Total number of MCMC iterations per chain. Default \code{5000}.
#' @param warmup Integer. Number of warmup (burn-in) iterations per chain.
#'   Default \code{2500}.
#' @param chains Integer. Number of MCMC chains to run. Default \code{2}.
#'
#' @return The \code{bvar} object with \code{fit$stan} and \code{posterior_means}
#'   appended.
#' @export
#'
#' @examples
#' \dontrun{
#' model <- bvar(data = my_data)
#' model <- setup(model, p = 2, deterministic = "constant")
#' model <- priors(model)
#' model <- fit(model, iter = 5000, warmup = 2500, chains = 2)
#' }
fit <- function(x, iter = 5000, warmup = 2500, chains = 2) {
  
  Jeffrey <- x$priors$Jeffrey
  
  stan_data <- c(
    x$setup,
    list(H = x$predict$H, d_pred = x$predict$d_pred),
    x$priors
  )
  
  stan_file <- if (isFALSE(Jeffrey)) {
    system.file("inv_wishart_cov.stan", package = "SteadyStateBVAR")
  } else {
    system.file("diffuse_cov.stan", package = "SteadyStateBVAR")
  }
  
  if (isTRUE(x$SV)) {
    if (x$SV_type == "RW") {
      stan_file <- system.file("stochastic_volatility_RW.stan", package = "SteadyStateBVAR")
    } else if (x$SV_type == "AR") {
      stan_file <- system.file("stochastic_volatility_stationaryAR.stan", package = "SteadyStateBVAR")
    } else {
      stop("Please specify a valid x$SV_type: either 'RW' or 'AR'.")
    }
    stan_data <- c(x$SV_priors, stan_data)
    k <- x$setup$k
    if (k == 2) {
      stan_data$theta_A <- as.array(x$SV_priors$theta_A[1])
    }
  }
  
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  x$fit$stan <- rstan::stan(
    file    = stan_file,
    data    = stan_data,
    iter    = iter,
    warmup  = warmup,
    chains  = chains,
    verbose = FALSE
  )
  
  posterior <- rstan::extract(x$fit$stan)
  
  if (!isTRUE(x$SV)) {
    x$posterior_means$beta    <- apply(posterior$beta,    c(2, 3), mean)
    x$posterior_means$Psi     <- apply(posterior$Psi,     c(2, 3), mean)
    x$posterior_means$Sigma_u <- apply(posterior$Sigma_u, c(2, 3), mean)
    
  } else if (isTRUE(x$SV) && x$SV_type == "AR") {
    x$posterior_means$beta    <- apply(posterior$beta,    c(2, 3), mean)
    x$posterior_means$Psi     <- apply(posterior$Psi,     c(2, 3), mean)
    x$posterior_means$A       <- apply(posterior$A,       c(2, 3), mean)
    x$posterior_means$gamma_0 <- apply(posterior$gamma_0, 2, mean)
    x$posterior_means$gamma_1 <- apply(posterior$gamma_1, 2, mean)
    x$posterior_means$Phi     <- apply(posterior$Phi,     c(2, 3), mean)
    
  } else if (isTRUE(x$SV) && x$SV_type == "RW") {
    x$posterior_means$beta <- apply(posterior$beta, c(2, 3), mean)
    x$posterior_means$Psi  <- apply(posterior$Psi,  c(2, 3), mean)
    x$posterior_means$A    <- apply(posterior$A,    c(2, 3), mean)
    x$posterior_means$phi  <- apply(posterior$phi,  2, mean)
  }
  
  return(x)
}