# Conditional forecasts from a fitted BVAR model

Computes and plots conditional forecasts from a fitted `bvar` object.
Conditions are imposed on specific variables at specific horizons using
the method of Waggoner and Zha (1999). Both conditional and
unconditional forecasts are plotted for comparison.

## Usage

``` r
conditional_forecast(
  bvar_obj,
  conditions,
  ci = 0.95,
  fcst_type = c("mean", "median"),
  growth_rate_idx = NULL,
  plot_idx = NULL
)
```

## Arguments

- bvar_obj:

  A `bvar` object that has been passed through
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md).

- conditions:

  A data frame with three columns: `var` (integer index of the
  conditioned variable), `horizon` (integer forecast horizon at which
  the condition is imposed), and `value` (numeric value the variable is
  conditioned on).

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

## Value

A list with three matrices: `forecast`, `lower`, and `upper`, each of
dimension `H x k`, as well as `cond_draws`, an array of all posterior
conditional forecast draws.

## References

Waggoner, D. F. and Zha, T. (1999). Conditional forecasts in dynamic
multivariate models. *Review of Economics and Statistics*, 81(4),
639-651.

## Examples

``` r
if (FALSE) { # \dontrun{
# Condition on variable 1 being 2.0 at horizon 4
conditions <- data.frame(var = 1, horizon = 4, value = 2.0)
model <- bvar(data = my_data)
model <- setup(model, p = 2, deterministic = "constant")
model <- priors(model)
model <- fit(model)
conditional_forecast(model, conditions = conditions)
} # }
```
