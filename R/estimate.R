estimate <- function(stan_data, n_chains, iter, warmup, H, X_pred, Jeffrey=FALSE) {
  stan_data$H <- H
  stan_data$X_pred <- X_pred
  if (isFALSE(Jeffrey)){
    stan_file <- system.file("STEADYSTATEBVAR2.stan", package = "SteadyStateBVAR")
  } else {
    stan_file <- system.file("STEADYSTATEBVAR3.stan", package = "SteadyStateBVAR")
  }
  rstan_options(auto_write = TRUE)
  options(mc.cores=parallel::detectCores())
  fit <- stan(stan_file,
              data=stan_data,
              chains=n_chains,
              iter=iter,
              warmup=warmup,
              verbose=FALSE)
  
  return(fit)
}
