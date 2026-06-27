# Plot stochastic volatility estimates and forecasts

Produces time series plots of posterior stochastic volatility estimates
over the estimation sample together with predictive paths over the
forecast horizon for a fitted `bvar` object.

## Usage

``` r
stochastic_volatility_plot(
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

  A `bvar` object estimated via
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md)
  with stochastic volatility enabled.

- ci:

  Numeric. Credible interval level. Default is `0.95`.

- vol:

  Character. Volatility representation: `"log_lambda"` (default) or
  `"sd"`.

- plot_idx:

  Integer vector. Variables to plot. If `NULL`, all variables are used.

- xlim:

  Numeric vector of length 2. Optional x-axis limits.

- ylim:

  Numeric vector of length 2. Optional y-axis limits.

## Value

An invisible list containing:

- estimated:

  List with posterior summaries for the estimation sample:

  - `mean`: T × k matrix of posterior means

  - `lower`: T × k matrix of lower credible bounds

  - `upper`: T × k matrix of upper credible bounds

- predicted:

  List with posterior summaries for the forecast horizon:

  - `mean`: (T+H) × k matrix of predictive means

  - `lower`: (T+H) × k matrix of lower credible bounds

  - `upper`: (T+H) × k matrix of upper credible bounds

## Details

The function supports two volatility representations:

- `"log_lambda"`: log-volatility states

- `"sd"`: implied posterior standard deviations

For each selected variable, the function plots:

- posterior mean path (in-sample)

- credible bands (in-sample)

- predictive mean path (out-of-sample)

- predictive credible bands (out-of-sample)

In addition, posterior summary matrices for both estimation and forecast
periods are returned invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- bvar(data = my_data)
model <- setup(model, p = 2, deterministic = "constant")
model <- priors(model)
model$SV <- TRUE
model$SV_type <- "RW"
model <- fit(model)

res <- stochastic_volatility_plot(model, ci = 0.95)
str(res)
} # }
```
