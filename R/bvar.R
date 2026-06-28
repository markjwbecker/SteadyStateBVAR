#' Create a steady-state BVAR model object
#'
#' Initialises a steady-state \code{bvar} object. This is the starting point for all models in \code{SteadyStateBVAR}.
#' After creation, pass the object sequentially to \code{\link{setup}}, \code{\link{priors}},
#' and \code{\link{fit}} to build and estimate the model.
#' 
#' @details
#' The model takes the form
#'
#' \deqn{y_t = \Psi d_t + \Pi_1(y_{t-1}-\Psi d_{t-1})+\dots+\Pi_p(y_{t-p}-\Psi d_{t-p})+u_t,}
#'
#' where \eqn{y_t} is an \eqn{k}-dimensional vector of endogenous variables at time \eqn{t}, and \eqn{d_t} is
#' a \eqn{q}-dimensional vector of deterministic (exogenous) variables at time \eqn{t}.
#' Here \eqn{\Pi_\ell} for \eqn{\ell=1,\dots,p} is a \eqn{(k \times k)} matrix of autoregressive parameters,
#' and \eqn{\Psi} is a \eqn{(k \times q)} matrix of steady-state parameters. Note that
#' \eqn{\mathrm{E}(y_t)=\mu_t=\Psi d_t} is the unconditional mean, or the **steady state** of the process.
#' In the case of the homoscedastic steady-state BVAR, we have \eqn{u_t \sim N_k(0,\Sigma_u)}, and for the models with
#' stochastic volatility we have \eqn{u_t \sim N_k(0,\Sigma_{u,t})}.
#'
#' @param data A numeric matrix or time series of data where each column is a
#'   variable and each row is a time period.
#'
#' @return An object of class \code{bvar}.
#' 
#' @references
#' Villani, M. (2009). Steady-state priors for vector autoregressions.
#' \emph{Journal of Applied Econometrics}. 24(4), pp. 630--649.
#'
#'
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