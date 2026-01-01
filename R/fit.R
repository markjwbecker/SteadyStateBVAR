fit <- function(x, iter = 5000, warmup = 2500, chains = 2, estimation = c("stan", "gibbs")) {
  
  estimation <- match.arg(estimation)
  Jeffrey <- x$priors$Jeffrey
  
  if (estimation == "stan") {
    
    stan_data <- c(
      x$setup,
      list(H = x$predict$H, X_pred = x$predict$X_pred),
      x$priors
    )
    
    stan_file <- if (isFALSE(Jeffrey)) {
      system.file("inv_wishart_cov.stan", package = "SteadyStateBVAR")
    } else {
      system.file("diffuse_cov.stan", package = "SteadyStateBVAR")
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
    
  } else if (method == "gibbs") {
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
