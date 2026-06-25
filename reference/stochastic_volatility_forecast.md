# Plot stochastic volatility estimates and forecasts

Plots the posterior distribution of stochastic volatility over the
estimation sample and the forecast horizon for a fitted `bvar` object
with a stochastic volatility specification. Supports plotting either
log-volatility (`log_lambda`) or posterior standard deviations (`sd`).

## Usage

``` r
stochastic_volatility_forecast(
  x,
  ci = 0.95,
  vol = "log_lambda",
  plot_idx = NULL,
  xlim = NULL,
  ylim = NULL
)
```

## Arguments

- x:

  A `bvar` object that has been passed through
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md)
  with a stochastic volatility specification.

- ci:

  Numeric. The credible interval width. Default `0.95`.

- vol:

  Character. The volatility measure to plot. Either `"log_lambda"`
  (default) for log-volatility or `"sd"` for posterior standard
  deviations.

- plot_idx:

  Integer vector. Indices of variables to plot. If `NULL` (default), all
  variables are plotted.

- xlim:

  Numeric vector of length 2. x-axis limits. If `NULL` (default), limits
  are set automatically.

- ylim:

  Numeric vector of length 2. y-axis limits. If `NULL` (default), limits
  are set automatically.

## Value

A list with two components:

- `estimated`:

  List containing posterior summaries for the estimation sample:

  - `mean`: Matrix of posterior means (T × k)

  - `lower`: Matrix of lower credible bounds (T × k)

  - `upper`: Matrix of upper credible bounds (T × k)

- `predicted`:

  List containing posterior summaries for the forecast horizon:

  - `mean`: Matrix of predictive means ((T+H) × k)

  - `lower`: Matrix of predictive lower credible bounds

  - `upper`: Matrix of predictive upper credible bounds

## Details

The function also returns posterior summaries for both the in-sample
estimates and out-of-sample forecasts.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- bvar(data = my_data)
model <- setup(model, p = 2, deterministic = "constant")
model <- priors(model)
model$SV <- TRUE
model$SV_type <- "RW"
model <- fit(model)

res <- stochastic_volatility_forecast(model, ci = 0.95)
str(res)
} # }
```
