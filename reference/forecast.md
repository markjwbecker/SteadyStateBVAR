# Forecast from a fitted BVAR model

Computes and plots forecasts from a fitted `bvar` object. Posterior
predictive draws from Stan are used to construct point forecasts and
credible intervals. Optionally converts level forecasts to annual growth
rates for selected variables.

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

  A `bvar` object that has been passed through
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md).

- ci:

  Numeric. The credible interval width. Default `0.95`.

- fcst_type:

  Character. Whether to use `"mean"` or `"median"` as the point
  forecast. Default `"mean"`.

- growth_rate_idx:

  Integer vector. Indices of variables to convert to annual growth
  rates. If `NULL` (default), all variables are plotted in levels.

- plot_idx:

  Integer vector. Indices of variables to plot. If `NULL` (default), all
  variables are plotted.

- show_all:

  Logical. If `FALSE` (default), only the last two years of history are
  shown alongside the forecast. If `TRUE`, the full history is shown.

## Value

A list with three matrices: `forecast`, `lower`, and `upper`, each of
dimension `H x k` where `H` is the forecast horizon and `k` is the
number of variables.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- bvar(data = my_data)
model <- setup(model, p = 2, deterministic = "constant")
model <- priors(model)
model <- fit(model)
fcst <- forecast(model, ci = 0.95, fcst_type = "mean")
} # }
```
