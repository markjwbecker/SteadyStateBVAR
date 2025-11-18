fit_gibbs <- function(x, iter = 5000, warmup = 1000) {
  Jeffrey = x$priors$Jeffrey
  x$fit$gibbs <- estimate_gibbs(
           x = x,
           iter = iter,
           warmup = warmup,
           H = x$forecasts$H,
           X_pred = x$forecasts$X_pred,
           Jeffrey = Jeffrey
           )
  
  return(x)
}
