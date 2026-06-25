# Impulse Response Functions for a BVAR model

Computes and plots impulse response functions (IRFs) from a fitted
`bvar` object. Supports both orthogonalised (OIRF) and generalised
(GIRF) impulse responses, with optional conversion to annual growth
rates.

## Usage

``` r
IRF(
  x,
  H = 16,
  response = NULL,
  shock = NULL,
  type = c("median", "mean"),
  method = c("OIRF", "GIRF"),
  ci = 0.95,
  t = NULL,
  growth_rate_idx = NULL
)
```

## Arguments

- x:

  A `bvar` object that has been passed through
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md).

- H:

  Integer. The forecast horizon for the IRF. Default `16`.

- response:

  Integer. Index of the response variable to plot. If `NULL` (default),
  all responses are plotted.

- shock:

  Integer. Index of the shock variable to plot. If `NULL` (default), all
  shocks are plotted.

- type:

  Character. Whether to use `"median"` or `"mean"` as the point
  estimate. Default `"median"`.

- method:

  Character. The IRF method: `"OIRF"` for orthogonalised or `"GIRF"` for
  generalised impulse responses. Default `"OIRF"`.

- ci:

  Numeric. The credible interval width. Default `0.95`.

- t:

  Integer. Time index for the covariance matrix when using stochastic
  volatility models. If `NULL` (default), the last time period is used.

- growth_rate_idx:

  Integer vector. Indices of variables to convert to annual growth
  rates. If `NULL` (default), no conversion is applied.

## Value

A list with three arrays: the point estimate IRF, `lower`, and `upper`
credible bounds, each of dimension `k x k x (H+1)`.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- bvar(data = my_data)
model <- setup(model, p = 2, deterministic = "constant")
model <- priors(model)
model <- fit(model)
irf <- IRF(model, H = 16, method = "OIRF")
} # }
```
