estimate <- function(stan_data, n_chains, iter, warmup, H, X_pred) {
  stan_data$H <- H
  stan_data$X_pred <- X_pred
  
  fit <- stan(file="/inst/STEADYSTATEBVAR.stan",
              data=stan_data,
              chains=n_chains,
              iter=iter,
              warmup=warmup,
              verbose=TRUE)
  
  return(fit)
}
