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

  Integer vector. Indices of variables of which to convert forecasts to
  annual growth rates \\\ln x\_{t} - \ln x\_{t-f}\\, where \\f\\ is the
  frequency of the data (4 for quarterly, 12 for monthly). Suitable for
  variables specified as \\\ln x\_{t} - \ln x\_{t-1}\\, i.e.
  `diff(log(x))` or `100*diff(log(x))`. Computed by summing up to \\f\\
  log first differences. Default is `NULL`.

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
                   Jeffreys = TRUE,
                   SV = FALSE,
                   SV_type = NULL,
                   SV_priors = NULL)
                   
bvar_obj <- fit(bvar_obj,
                H = 8,
                d_pred = matrix(rep(1,8)),
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1)
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_homoscedastic_jeffreys_prior' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000153 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.53 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 7
#> Chain 1:            adapt_window = 38
#> Chain 1:            term_buffer = 5
#> Chain 1: 
#> Chain 1: Iteration:   1 / 200 [  0%]  (Warmup)
#> Chain 1: Iteration:  20 / 200 [ 10%]  (Warmup)
#> Chain 1: Iteration:  40 / 200 [ 20%]  (Warmup)
#> Chain 1: Iteration:  51 / 200 [ 25%]  (Sampling)
#> Chain 1: Iteration:  70 / 200 [ 35%]  (Sampling)
#> Chain 1: Iteration:  90 / 200 [ 45%]  (Sampling)
#> Chain 1: Iteration: 110 / 200 [ 55%]  (Sampling)
#> Chain 1: Iteration: 130 / 200 [ 65%]  (Sampling)
#> Chain 1: Iteration: 150 / 200 [ 75%]  (Sampling)
#> Chain 1: Iteration: 170 / 200 [ 85%]  (Sampling)
#> Chain 1: Iteration: 190 / 200 [ 95%]  (Sampling)
#> Chain 1: Iteration: 200 / 200 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.035 seconds (Warm-up)
#> Chain 1:                0.085 seconds (Sampling)
#> Chain 1:                0.12 seconds (Total)
#> Chain 1: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
                
IRF(bvar_obj)

# }
```
