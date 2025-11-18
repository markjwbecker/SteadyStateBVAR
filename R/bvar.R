bvar <- function(data = NULL, setup = NULL, priors = NULL, fit = list(gibbs = NULL, stan = NULL), forecasts=list()) {
  
  obj <- list(
    data   = data,
    setup  = setup,
    priors = priors,
    fit    = fit,
    forecasts = forecasts
  )
  
  class(obj) <- "bvar"
  return(obj)
}
