bvar <- function(data = NULL, setup = NULL, priors = NULL, fit = list(gibbs = NULL, stan = NULL)) {
  
  obj <- list(
    data   = data,
    setup  = setup,
    priors = priors,
    fit    = fit
  )
  
  class(obj) <- "bvar"
  return(obj)
}
