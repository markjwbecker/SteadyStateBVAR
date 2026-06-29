#' Koop and Korobilis (2010) dataset
#'
#' Quarterly US macroeconomic data (1953Q1--2006Q3) from Koop and Korobilis (2010)
#' used for estimating Bayesian VAR models.
#'
#' @format ## `KoopKorobilis2010`
#' A multivariate time series object with 215 rows and 3 columns:
#' \describe{
#'   \item{delta pi}{inflation rate (annual percentage change in a chain-weighted GDP price index)}
#'   \item{u}{unemployment rate (seasonally adjusted civilian unemployment rate, all civilian workers aged 16 years or older)}
#'   \item{r}{interest rate (yield on the three-month Treasury bill rate)}
#' }
#' @source Koop, G. and Korobilis, D. (2010). Bayesian Multivariate Time Series
#'  Methods for Empirical Macroeconomics. \emph{Foundations and Trends in Econometrics}. 3(4), pp. 267-358.
"KoopKorobilis2010"

#' Villani (2009) dataset
#'
#' Swedish macroeconomic data (1980Q1–2005Q4) from Section 4.1 in Villani (2009)
#'
#' @format ## `Villani2009`
#' A multivariate time series object with 104 rows and 7 columns:
#' \describe{
#'   \item{delta y_f}{foreign GDP growth}
#'   \item{pi_f}{foreign CPI inflation}
#'   \item{i_f}{foreign 3-month interest rate}
#'   \item{delta y}{domestic GDP growth}
#'   \item{pi}{domestic CPI inflation}
#'   \item{i}{domestic 3-month interest rate}
#'   \item{q}{level of the real exchange rate defined as q = s + p_f - p, where p_f and p
#'   are the foreign and domestic CPI levels in logs and s is the log of the trade-weighted nominal exchange rate}
#'   ...
#' }
#' @source Villani, M. (2009). Steady-state priors for vector autoregressions.
#'   \emph{Journal of Applied Econometrics}, 24(4), pp. 630-650.
"Villani2009"