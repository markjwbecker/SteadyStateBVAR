#' @details
#' The steady-state BVAR model takes the form
#'
#' \deqn{y_t = \Psi d_t + \Pi_1(y_{t-1}-\Psi d_{t-1})+\dots+\Pi_p(y_{t-p}-\Psi d_{t-p})+u_t}
#'
#' where \eqn{y_t} is an \eqn{k}-dimensional vector of endogenous variables at time \eqn{t}, and \eqn{d_t} is
#' a \eqn{q}-dimensional vector of deterministic (exogenous) variables at time \eqn{t}.
#' Here \eqn{\Pi_\ell} for \eqn{\ell=1,\dots,p} is a \eqn{(k \times k)} matrix of autoregressive parameters,
#' and \eqn{\Psi} is a \eqn{(k \times q)} matrix of steady-state parameters. Now
#'
#' \deqn{\mathrm{E}(y_t)=\mu_t=\Psi d_t}
#'
#' is the unconditional mean, or the \strong{steady state} of the process.
#' 
#' For the autoregressive parameters, the Minnesota prior is used. For the steady-state parameters,
#' a normal prior is used (which is hopefully informative, at least for some of the parameters, otherwise there is no point in using this model)
#'
#' The package currently supports three specifications for the (reduced-form) innovations \eqn{u_t}
#' 
#' \enumerate{
#'   \item \eqn{u_t \overset{\text{iid}}{\sim} \mathrm{N_k}(0,\Sigma_u)}
#'   \item \eqn{u_t \sim \mathrm{N_k}(0,\Sigma_{u,t})}, where \eqn{\Sigma_{u,t}} is driven by a latent Random Walk process (stochastic volatility)
#'   \item \eqn{u_t \sim \mathrm{N_k}(0,\Sigma_{u,t})}, where \eqn{\Sigma_{u,t}} is driven by a latent AR(1) process (stochastic volatility)
#' }
#' 
#' For more information, please see the package vignettes
#' 
#' \enumerate{
#'   \item \code{vignette("steady-state-BVAR-intro")}
#'   \item \code{vignette("Homoscedastic-steady-state-BVAR")}
#'   \item \code{vignette("RW-stochastic-volatility-steady-state-BVAR")}
#'   \item \code{vignette("AR1-stochastic-volatility-steady-state-BVAR")}
#' }
#' 
#' @references
#' Villani, M. (2009). Steady-state priors for vector autoregressions.
#' \emph{Journal of Applied Econometrics}, 24(4), pp. 630–650.
#'
#' @keywords internal
"_PACKAGE"
## usethis namespace: start
#' @importFrom grDevices rgb
#' @importFrom graphics abline axis lines par points polygon title legend grid
#' @importFrom stats embed frequency plot.ts qnorm quantile rnorm setNames start time ts ts.plot is.ts median
#' @importFrom utils head tail
#' @import Rcpp
#' @import methods
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib SteadyStateBVAR, .registration = TRUE
## usethis namespace: end
NULL