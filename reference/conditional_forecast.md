# Conditional forecasts from a fitted BVAR model

Computes and plots conditional forecasts from a fitted steady-state
`bvar` object. Conditions are imposed on specific variables at specific
horizons using the method of Dieppe, Legrand, and van Roye (2018). Both
conditional and unconditional forecasts are plotted for comparison.
Please note that for the moment, conditional forecasting is only enabled
for the homoscedastic steady-state BVAR, i.e. when `SV=FALSE` in
[`priors()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md).

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

  A steady-state `bvar` object that has been passed through
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md).

- conditions:

  A data frame with three columns: `var` (integer index of the
  conditioned variable), `horizon` (integer forecast horizon at which
  the condition is imposed), and `value` (numeric, the imposed
  condition). Note that `var` and `horizon` must be integers, as they
  are used for indexing.

- ci:

  Numeric. The prediction interval width. Default `0.95`.

- fcst_type:

  Character. Whether to use `"mean"` or `"median"` as the point
  forecast. Default `"mean"`.

- growth_rate_idx:

  Integer vector. Indices of variables to convert to annual growth
  rates. Default is `NULL`.

- plot_idx:

  Integer vector. Indices of variables to plot. If `NULL` (default), all
  variables are plotted.

## Value

Invisibly returns a list with three matrices: `forecast`, `lower`, and
`upper`, each of dimension `H x k`, as well as `cond_draws`, an array of
all posterior conditional forecast draws.

## References

Dieppe, A., Legrand, R., and van Roye, B. (2018). *The Bayesian
Estimation, Analysis and Regression (BEAR) Toolbox Technical guide*.
European Central Bank.

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
                
conditions <- data.frame(var = rep(2,8),
                         horizon = rep(1:8),
                         value   = rep(1,8))
                         
cond_fcst <- conditional_forecast(bvar_obj,
                                  conditions,
                                  ci=0.68,
                                  fcst_type = "mean")
} # }
```
