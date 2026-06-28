# Forecast from a fitted BVAR model

Computes and plots forecasts from a fitted `bvar` object. Posterior
predictive draws from Stan are used to construct point forecasts and
prediction intervals. Optionally converts monthly or quarterly growth
rate forecasts to annual growth rates for selected variables.

## Usage

``` r
forecast(
  x,
  ci = 0.95,
  fcst_type = c("mean", "median"),
  growth_rate_idx = NULL,
  plot_idx = NULL,
  show_all = FALSE
)
```

## Arguments

- x:

  A steady-state `bvar` object that has been passed through
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md).

- ci:

  Numeric. The prediction interval width. Default `0.95`, i.e. 95%
  prediction interval.

- fcst_type:

  Character. Whether to use `"mean"` or `"median"` as the point
  forecast. Default `"mean"`.

- growth_rate_idx:

  Integer vector. Indices of variables to convert to annual growth
  rates. Default is `NULL`.

- plot_idx:

  Integer vector. Indices of variables to plot. If `NULL` (default), all
  variables are plotted.

- show_all:

  Logical. If `FALSE` (default), only the last two years of history are
  shown alongside the forecast. If `TRUE`, the full history is shown.

## Value

Invisibly returns a list with three matrices: `forecast`, `lower`, and
`upper`, each of dimension `H x k` where `H` is the forecast horizon and
`k` is the number of variables.

## Examples

``` r
if (FALSE) { # \dontrun{
#homoscedastic with Jeffreys prior
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")

bvar_obj <- priors(bvar_obj,
                   lambda_1 = 0.2,
                   lambda_2 = 0.5,
                   lambda_3 = 1,
                   first_own_lag_prior_mean = rep(1,2),
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2),
                   Jeffrey = TRUE,
                   SV = FALSE,
                   SV_type = NULL,
                   SV_priors = NULL)
                   
bvar_obj <- fit(bvar_obj,
                H = 8,
                d_pred = matrix(rep(1,8)),
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1,
                verbose = FALSE,
                auto_write = FALSE)

forecast(bvar_obj, ci = 0.90, show_all = TRUE)
} # }
```
