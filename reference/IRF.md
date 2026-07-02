# Impulse Response Functions for a fitted steady-state BVAR model

Computes and plots impulse response functions (IRFs) from a fitted
steady-state `bvar` object. Supports both orthogonalized (OIRF) and
generalized (GIRF) impulse responses, with optional conversion to annual
growth rates.

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

  A steady-state `bvar` object that has been passed through
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

  Character. The IRF method: `"OIRF"` for orthogonalized or `"GIRF"` for
  generalized impulse responses. Default `"OIRF"`.

- ci:

  Numeric. The credible interval width. Default `0.95`, i.e. 95%.

- t:

  Integer. Time index for the covariance matrix when using stochastic
  volatility models. If `NULL` (default), the last time `t` is used.

- growth_rate_idx:

  Integer vector. Indices of variables to convert to annual growth
  rates. Default is `NULL`.

## Value

Invisibly returns a list with three arrays: the point estimate IRF,
`lower`, and `upper` credible bounds, each of dimension `k x k x (H+1)`.

## Examples

``` r
# \donttest{
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
#> Error in stan_model(file, model_name = model_name, model_code = model_code,     stanc_ret = NULL, boost_lib = boost_lib, eigen_lib = eigen_lib,     save_dso = save_dso, verbose = verbose): Boost not found; call install.packages('BH')
                
IRF(bvar_obj)
#> Error: unable to find an inherited method for function ‘extract’ for signature ‘object = "NULL"’
# }
```
