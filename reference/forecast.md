# Forecast from a fitted steady-state BVAR model

Computes and plots forecasts from a fitted steady-state `bvar` object.
Draws from the joint predictive distribution are used to construct point
forecasts and prediction intervals. Optionally converts monthly or
quarterly growth rate forecasts to annual growth rates for variables
specified as `100*diff(log(x))`. Optionally overlays the posterior
steady states.

## Usage

``` r
forecast(
  x,
  pi = 0.95,
  fcst_type = c("mean", "median"),
  growth_rate_idx = NULL,
  plot_idx = NULL,
  show_all = TRUE,
  ss = TRUE,
  ss_type = c("mean", "median"),
  ss_ci = 0.95
)
```

## Arguments

- x:

  A steady-state `bvar` object that has been passed through
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md).

- pi:

  Numeric. The prediction interval width. Default `0.95`, i.e. 95%
  prediction interval.

- fcst_type:

  Character. Whether to use `"mean"` or `"median"` as the point
  forecast. Default `"mean"`.

- growth_rate_idx:

  Integer vector. Indices of variables of which to convert forecasts to
  annual growth rates \\100 (\ln x\_{t} - \ln x\_{t-f})\\, where \\f\\
  is the frequency of the data (4 for quarterly, 12 for monthly). Only
  suitable for variables specified as \\100 (\ln x\_{t} - \ln
  x\_{t-1})\\, i.e. `100*diff(log(x))`. Computed by summing up to \\f\\
  log first differences. Default is `NULL`.

- plot_idx:

  Integer vector. Indices of variables to plot. If `NULL` (default), all
  variables are plotted. Forecasts are always computed and returned for
  all variables, regardless of `plot_idx`.

- show_all:

  Logical. If `TRUE` (default), the full history is shown. If `FALSE`,
  the last eight years of history are shown.

- ss:

  Logical. If `TRUE`, overlays the posterior steady state \\\mu_t = \Psi
  d_t\\. Default `TRUE`. For variables in `growth_rate_idx`, the
  posterior steady-state is annualized.

- ss_type:

  Character. Whether to use `"mean"` or `"median"` of the posterior as
  the steady-state point estimate. Default `"mean"`.

- ss_ci:

  Numeric. The posterior credible interval width for the steady state.
  Default `0.95`, i.e. 95% credible interval.

## Value

Invisibly returns a list with three matrices: `forecast`, `lower`, and
`upper`, each of dimension `H x k` where `H` is the forecast horizon and
`k` is the number of variables.

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
                   Jeffreys = TRUE)
                   
bvar_obj <- fit(bvar_obj,
                H = 8,
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1)
#> ------------------------------------------------------------
#> Forecast horizon:
#> 8
#> 
#> Future deterministic variables (d_pred):
#>     constant
#> h=1        1
#> h=2        1
#> h=3        1
#> h=4        1
#> h=5        1
#> h=6        1
#> h=7        1
#> h=8        1
#> ------------------------------------------------------------
#> Estimating Stan model:
#> steady_state_bvar_homoscedastic_jeffreys_prior
#> 
#> Also generating draws from the joint predictive distribution
#> 
#> ...
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_homoscedastic_jeffreys_prior' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.57 seconds.
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
#> Chain 1:  Elapsed Time: 0.015 seconds (Warm-up)
#> Chain 1:                0.046 seconds (Sampling)
#> Chain 1:                0.061 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> SAMPLING FINISHED
                
forecast(bvar_obj,
         fcst_type = "mean",
         pi = 0.90,
         show_all = TRUE,
         ss = TRUE,
         ss_type = "mean",
         ss_ci = 0.95)


# }
```
