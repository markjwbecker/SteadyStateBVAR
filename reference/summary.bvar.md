# Summarise a fitted steady-state BVAR model

Computes and prints posterior summaries from a fitted steady-state
`bvar` object. The printed output depends on whether the model is
homoscedastic or includes stochastic volatility (RW or AR1
specification).

## Usage

``` r
# S3 method for class 'bvar'
summary(object, pars = NULL, stat = "mean", t = NULL, ...)
```

## Arguments

- object:

  A steady-state `bvar` object that has been passed through
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md).

- pars:

  Character vector of parameter names to include. If `NULL` (default),
  all available parameters are displayed.

- stat:

  Character. Posterior summary statistic to display: `"mean"` (default)
  or `"median"`.

- t:

  Integer. Time index for the innovation covariance matrix if stochastic
  volatility was estimated. If `NULL`, the latest available `t` is used.

- ...:

  Additional arguments (currently unused).

## Value

Returns the input object invisibly.

## Details

The function summarises the following estimated parameters:

- `beta`: kp×k VAR coefficient matrix

- `Psi`: k×q steady-state parameter matrix

- `Sigma_u`: innovation covariance matrix (k×k for homoscedastic, T×k×k
  for stochastic volatility)

- If Random Walk stochastic volatility:

  - `A`: k×k lower triangular matrix with ones on the diagonal that
    describes the contemporaneous interaction of the endogenous
    variables

  - `phi`: k-dimensional vector of log volatility innovation variances

- If AR1 stochastic volatility:

  - `A`: k×k lower triangular matrix with ones on the diagonal that
    describes the contemporaneous interaction of the endogenous
    variables

  - `gamma_0`: k-dimensional vector of log volatility intercepts

  - `gamma_1`: k-dimensional vector of log volatility slopes

  - `Phi`: k×k log volatility innovation covariance matrix

## Examples

``` r
# \donttest{
yt <- matrix(rnorm(50), 25, 2)
bvar_obj <- bvar(data = yt)
bvar_obj <- setup(bvar_obj, p = 1, deterministic = "constant")
bvar_obj <- priors(bvar_obj,
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2))

bvar_obj <- fit(bvar_obj,
                H = 1,
                d_pred = matrix(1),
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1,
                verbose = FALSE,
                auto_write = FALSE)
#> Error in stan_model(file, model_name = model_name, model_code = model_code,     stanc_ret = NULL, boost_lib = boost_lib, eigen_lib = eigen_lib,     save_dso = save_dso, verbose = verbose): Boost not found; call install.packages('BH')

summary(bvar_obj)
#> Error in summary.bvar(bvar_obj): object must be passed through fit() first
# }
```
