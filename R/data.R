#' @importFrom grDevices rgb
#' @importFrom graphics abline axis lines par points polygon title
#' @importFrom stats embed frequency plot.ts qnorm quantile rnorm setNames start time ts ts.plot
#' @importFrom utils head tail
#' @keywords internal
NULL

#' Koop and Korobilis (2010) Dataset
#'
#' A macroeconomic dataset used in Koop and Korobilis (2010) for estimating
#' Bayesian VAR models.
#'
#' @format A time series matrix with macroeconomic variables as columns.
#' @source Koop, G. and Korobilis, D. (2010). Bayesian multivariate time
#'   series methods for empirical macroeconomics. \emph{Foundations and
#'   Trends in Econometrics}, 3(4), 267-358.
"KoopKorobilis2010"

#' Villani (2009) Dataset
#'
#' A macroeconomic dataset used in Villani (2009) for estimating Bayesian
#' VAR models with steady-state priors.
#'
#' @format A time series matrix with macroeconomic variables as columns.
#' @source Villani, M. (2009). Steady-state priors for vector autoregressions.
#'   \emph{Journal of Applied Econometrics}, 24(4), 630-650.
"villani2009"