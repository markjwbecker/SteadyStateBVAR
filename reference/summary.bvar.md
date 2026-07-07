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

- `beta`: \\kp \times k\\ VAR coefficient matrix

- `Psi`: \\k \times q\\ steady-state parameter matrix

- `Sigma_u`: innovation covariance matrix (\\k \times k\\ for
  homoscedastic, \\T \times k \times k\\ for stochastic volatility)

- If Random Walk stochastic volatility:

  - `A`: \\k \times k\\ lower triangular matrix with ones on the
    diagonal that describes the contemporaneous interaction of the
    endogenous variables

  - `phi`: \\k\\-dimensional vector of log volatility innovation
    variances

- If AR1 stochastic volatility:

  - `A`: \\k \times k\\ lower triangular matrix with ones on the
    diagonal that describes the contemporaneous interaction of the
    endogenous variables

  - `gamma_0`: \\k\\-dimensional vector of log volatility intercepts

  - `gamma_1`: \\k\\-dimensional vector of log volatility slopes

  - `Phi`: \\k \times k\\ log volatility innovation covariance matrix

## Examples

``` r
# \donttest{
yt <- matrix(rnorm(20), 10, 2)
bvar_obj <- bvar(data = yt)
bvar_obj <- setup(bvar_obj, p = 1, deterministic = "constant")
bvar_obj <- priors(bvar_obj,
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2))

bvar_obj <- fit(bvar_obj,
                H = 1,
                d_pred = matrix(1),
                iter = 100,
                warmup = 50,
                chains = 1,
                cores = 1,
                verbose = FALSE)
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_homoscedastic_jeffreys_prior' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.47 seconds.
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
#> Chain 1: Iteration:  1 / 100 [  1%]  (Warmup)
#> Chain 1: Iteration: 10 / 100 [ 10%]  (Warmup)
#> Chain 1: Iteration: 20 / 100 [ 20%]  (Warmup)
#> Chain 1: Iteration: 30 / 100 [ 30%]  (Warmup)
#> Chain 1: Iteration: 40 / 100 [ 40%]  (Warmup)
#> Chain 1: Iteration: 50 / 100 [ 50%]  (Warmup)
#> Chain 1: Iteration: 51 / 100 [ 51%]  (Sampling)
#> Chain 1: Iteration: 60 / 100 [ 60%]  (Sampling)
#> Chain 1: Iteration: 70 / 100 [ 70%]  (Sampling)
#> Chain 1: Iteration: 80 / 100 [ 80%]  (Sampling)
#> Chain 1: Iteration: 90 / 100 [ 90%]  (Sampling)
#> Chain 1: Iteration: 100 / 100 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.018 seconds (Warm-up)
#> Chain 1:                0.015 seconds (Sampling)
#> Chain 1:                0.033 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.14, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

summary(bvar_obj)
#> Posterior mean estimates
#> ------------------------
#> 
#> 
#> beta
#> --------------------------------------------------------------------------------         
#>            Var1  Var2
#>   Var1.l1  0.06 -0.05
#>   Var2.l1 -0.02 -0.07
#> --------------------------------------------------------------------------------
#> 
#> 
#> Psi
#> --------------------------------------------------------------------------------      
#>        [,1]
#>   Var1 0.11
#>   Var2 0.11
#> --------------------------------------------------------------------------------
#> 
#> 
#> Sigma_u
#> --------------------------------------------------------------------------------      
#>        Var1 Var2
#>   Var1 1.93 0.86
#>   Var2 0.86 0.90
#> --------------------------------------------------------------------------------
#> 
# }
```
