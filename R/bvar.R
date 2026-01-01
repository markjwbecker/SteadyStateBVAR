bvar <- function(data = NULL, setup = NULL, priors = NULL, fit = list(stan = NULL, gibbs = NULL), predict=list()) {
  
  obj <- list(
    data      = data,
    setup     = setup,
    priors    = priors,
    fit       = fit,
    predict   = predict
  )
  
  class(obj) <- "bvar"
  return(obj)
}
