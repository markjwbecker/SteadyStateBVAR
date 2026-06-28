#' Create a steady-state BVAR model object
#'
#' Initialises a A steady-state \code{bvar} object. This is
#' the starting point for all models in \code{SteadyStateBVAR}. After creation,
#' pass the object sequentially to \code{\link{setup}}, \code{\link{priors}}, 
#' and \code{\link{fit}} to build and estimate the model.
#'
#' @param data A numeric matrix or time series of data where each column is a
#'   variable and each row is a time period.
#'
#' @return An object of class \code{bvar}.
#' @export
#'
#' @examples
#' yt <- matrix(rnorm(50), 25, 2)
#' 
#' bvar_obj <- bvar(data = yt)
bvar <- function(data) {
  
  if (!is.matrix(data) && !is.ts(data)) {
    stop("data must be a matrix or time series object")
  }
  
  obj <- list(
    data    = data,
    setup   = NULL,
    priors  = NULL,
    fit     = NULL,
    predict = list()
  )
  
  class(obj) <- "bvar"
  return(obj)
}