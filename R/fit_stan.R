fit_stan <- function(x, iter, warmup) {

  Jeffrey = x$priors$Jeffrey
  stan_data <- x$setup
  
  stan_data <- c(
    stan_data,
    list(
      H = x$H,
      X_pred = x$X_pred
    ),
    x$priors
  )
  
  
  if (isFALSE(Jeffrey)) {
    stan_file <- system.file("STEADYSTATEBVAR2.stan", package = "SteadyStateBVAR")
  } else {
    stan_file <- system.file("STEADYSTATEBVAR3.stan", package = "SteadyStateBVAR")
  }
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  stan_fit <- rstan::stan(
    file = stan_file,
    data = stan_data,
    iter = iter,
    warmup = warmup,
    chains = 2,
    verbose = FALSE
  )
  x$fit$stan <- stan_fit
  return(x)
}

