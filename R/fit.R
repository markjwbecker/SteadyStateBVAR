fit <- function(x, iter = 5000, warmup = 2500, chains = 2, estimation = c("stan", "gibbs")) {
  
  estimation <- match.arg(estimation)
  Jeffrey <- x$priors$Jeffrey
  
  if (estimation == "stan") {
    
    if (!is.null(x$SV)) {
      stan_data <- c(
        x$setup,
        x$priors
      )
    } else {
      stan_data <- c(
        x$setup,
        list(H = x$predict$H, X_pred = x$predict$X_pred),
        x$priors
      )
    }
    
    stan_file <- if (isFALSE(Jeffrey)) {
      system.file("inv_wishart_cov.stan", package = "SteadyStateBVAR")
    } else {
      system.file("diffuse_cov.stan", package = "SteadyStateBVAR")
    }
    
    if (isTRUE(x$SV)) {
      stan_file <- system.file("stochastic_volatility.stan", package = "SteadyStateBVAR")
      stan_data <- c(x$SV_priors, stan_data)
      k <- x$setup$k
      if (k == 2) {
        stan_data$theta_A <- as.array(x$SV_priors$theta_A[1])
      }
      
    }
    rstan::rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
    
    x$fit$stan <- rstan::stan(
      file = stan_file,
      data = stan_data,
      iter = iter,
      warmup = warmup,
      chains = chains,
      verbose = FALSE
    )
    
  } else {
    
    if (chains != 1) {
      stop("For Gibbs sampling, 'chains' must be equal to 1.")
    }
    
    x$fit$gibbs <- estimate_gibbs(
      x = x,
      iter = iter,
      warmup = warmup,
      H = x$predict$H,
      X_pred = x$predict$X_pred,
      Jeffrey = Jeffrey
    )
  }
  
  return(x)
}