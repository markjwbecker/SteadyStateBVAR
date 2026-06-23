#' Create a BVAR model object
#'
#' Initialises a Bayesian Vector Autoregression (BVAR) model object. This is
#' the starting point for all models in \code{SteadyStateBVAR}. After creation,
#' pass the object to \code{\link{setup}}, \code{\link{priors}}, and
#' \code{\link{fit}} sequentially to build and estimate the model.
#'
#' @param data A numeric matrix of time series data where each column is a
#'   variable and each row is a time period.
#' @param setup Output from \code{\link{setup}}. Default \code{NULL}.
#' @param priors Output from \code{\link{priors}}. Default \code{NULL}.
#' @param fit Output from \code{\link{fit}}. Default \code{NULL}.
#' @param predict A list of prediction results. Default empty list.
#'
#' @return An object of class \code{bvar}.
#' @export
#'
#' @examples
#' \dontrun{
#' model <- bvar(data = my_data)
#' }
bvar <- function(data = NULL, setup = NULL, priors = NULL, fit = NULL, predict = list()) {
  
  obj <- list(
    data    = data,
    setup   = setup,
    priors  = priors,
    fit     = fit,
    predict = predict
  )
  
  class(obj) <- "bvar"
  return(obj)
}