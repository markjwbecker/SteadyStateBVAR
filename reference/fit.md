# Estimate the steady-state BVAR model using Stan

Estimates the steady-state BVAR model using the No-U-Turn sampler (a
variant of Hamiltonian Monte Carlo) via Stan. Also generates draws from
the joint predictive distribution. Uses the data, setup, and priors
stored in the steady-state `bvar` object.

## Usage

``` r
fit(x, H = 1, d_pred = NULL, ...)
```

## Arguments

- x:

  A steady-state `bvar` object that has been passed through
  [`setup`](https://markjwbecker.github.io/SteadyStateBVAR/reference/setup.md)
  and
  [`priors`](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md).

- H:

  Positive Integer. Forecast horizon. Default is `1`.

- d_pred:

  Matrix of size \\H \times q\\. Future values of the deterministic
  variables \\d_t\\. Default is `NULL` (`d_pred` is automatically
  created). If \\d_t\\ contains only a constant, \\d\_{t+h}=1 \\ \forall
  \\ h\\. If \\d_t\\ contains a constant and a dummy, it is assumed that
  the dummy stays at its last observed value for all future forecast
  periods. However, if this is not the intention, the user may supply
  `d_pred` themselves. Naturally, the constant is equal to one for all
  future periods. If \\d_t\\ contains a constant and a trend, the trend
  is extrapolated from its last observed value (i.e.
  \\trend\_{T+h}=trend\_{T}+h\\). Naturally, the constant is equal to
  one for all future periods.

- ...:

  Additional arguments passed directly to the `rstan` function
  [`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  (e.g. `iter`, `warmup`, `chains`, `cores`, `control`, `seed`, `init`,
  `thin`, `algorithm`, `pars`, `include`, `refresh`, `verbose`,
  `save_warmup`, `sample_file`, `diagnostic_file`). If `pars`/`include`
  is used to exclude model parameters required by `fit()` for posterior
  summaries, an error will be raised.

## Value

A fitted steady-state `bvar` object with:

- `x$fit$stan`: An object of class `stanfit`

- `x$fit$posterior_means`: List of posterior mean estimates

- `x$fit$posterior_medians`: List of posterior median estimates

## Details

The function selects the appropriate precompiled Stan model based on
settings from
[priors](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md):

- `steady_state_bvar_homoscedastic_jeffreys_prior`

- `steady_state_bvar_homoscedastic_inverse_wishart_prior`

- `steady_state_bvar_RW_stochastic_volatility`

- `steady_state_bvar_AR1_stochastic_volatility`

The function estimates the following parameters (see
[bvar](https://markjwbecker.github.io/SteadyStateBVAR/reference/bvar.md)
for details):

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
#homoscedastic with Jeffreys prior, d_t = constant

yt <- matrix(rnorm(20), 10, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p = 1, deterministic = "constant")

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
#> Estimating Stan model:
#> steady_state_bvar_homoscedastic_jeffreys_prior
#> 
#> Also generating draws from the joint predictive distribution
#> 
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
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_homoscedastic_jeffreys_prior' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.42 seconds.
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
#> Chain 1:  Elapsed Time: 0.011 seconds (Warm-up)
#> Chain 1:                0.03 seconds (Sampling)
#> Chain 1:                0.041 seconds (Total)
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
                
                
#homoscedastic with inverse-Wishart prior, d_t = constant and dummy

yt <- matrix(rnorm(20), 10, 2)

bvar_obj <- bvar(data = yt)

dummy_variable <- c(rep(1,5), rep(0,5))

bvar_obj <- setup(bvar_obj, p = 1,
                  deterministic = "constant_and_dummy",
                  dummy = dummy_variable)
                  
k <- bvar_obj$setup$k
q <- bvar_obj$setup$q

bvar_obj <- priors(bvar_obj,
                   lambda_1 = 0.2,
                   lambda_2 = 0.5,
                   lambda_3 = 1,
                   first_own_lag_prior_mean = rep(1,2),
                   theta_Psi = rep(0, k*q),
                   Omega_Psi = diag(0.1, k*q, k*q),
                   Jeffreys = FALSE) #inverse-Wishart


bvar_obj <- fit(bvar_obj,
                H = 8,
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1)
#> NOTE: d_pred not supplied
#> it is assumed that the dummy stays at its last observed value (0) for all 8 forecast periods.
#> ------------------------------------------------------------
#> Estimating Stan model:
#> steady_state_bvar_homoscedastic_inverse_wishart_prior
#> 
#> Also generating draws from the joint predictive distribution
#> 
#> Forecast horizon:
#> 8
#> 
#> Future deterministic variables (d_pred):
#>     constant dummy
#> h=1        1     0
#> h=2        1     0
#> h=3        1     0
#> h=4        1     0
#> h=5        1     0
#> h=6        1     0
#> h=7        1     0
#> h=8        1     0
#> ------------------------------------------------------------
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_homoscedastic_inverse_wishart_prior' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.35 seconds.
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
#> Chain 1:  Elapsed Time: 0.01 seconds (Warm-up)
#> Chain 1:                0.024 seconds (Sampling)
#> Chain 1:                0.034 seconds (Total)
#> Chain 1: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
 

#RW stochastic volatility

yt <- matrix(rnorm(20), 10, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")

k <- bvar_obj$setup$k
n_free_params_A <- bvar_obj$setup$n_free_params_A

SV_priors_RW <- list(
theta_A              =  rep(0, n_free_params_A),
Omega_A              =  diag(1000, n_free_params_A),
theta_log_lambda_1   =  rep(0, k),
Omega_log_lambda_1   =  diag(1000, k),
alpha_phi            =  rep(5, k),
beta_phi             = (rep(5, k) - 1) * rep(0.1, k)
)

bvar_obj <- priors(bvar_obj,
                   lambda_1 = 0.2,
                   lambda_2 = 0.5,
                   lambda_3 = 1,
                   first_own_lag_prior_mean = rep(1,2),
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2),
                   SV = TRUE,
                   SV_type = "RW",
                   SV_priors = SV_priors_RW)

bvar_obj <- fit(bvar_obj,
                H = 8,
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1,
                control = list(max_treedepth = 12, adapt_delta = 0.85)
                )
#> ------------------------------------------------------------
#> Estimating Stan model:
#> steady_state_bvar_RW_stochastic_volatility
#> 
#> Also generating draws from the joint predictive distribution
#> 
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
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_RW_stochastic_volatility' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.51 seconds.
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
#> Chain 1:  Elapsed Time: 0.005 seconds (Warm-up)
#> Chain 1:                0.027 seconds (Sampling)
#> Chain 1:                0.032 seconds (Total)
#> Chain 1: 
#> Warning: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess


#AR1 stochastic volatility

yt <- matrix(rnorm(20), 10, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")

k <- bvar_obj$setup$k
n_free_params_A <- bvar_obj$setup$n_free_params_A

SV_priors_AR1 <- list(
theta_A               =  rep(0, n_free_params_A),
Omega_A               =  diag(1000, n_free_params_A),
theta_gamma_0         =  rep(0.1, k),
Omega_gamma_0         =  diag(1000, k),
theta_gamma_1         =  rep(0.9, k),
Omega_gamma_1         =  diag(10, k),
theta_log_lambda_1    =  rep(0.1, k)/(1-rep(0.9, k)),
Omega_log_lambda_1    =  diag(1000, k),
V_Phi                 = (10 - k - 1) * diag(k),
m_Phi                 =  10
)

bvar_obj <- priors(bvar_obj,
                   lambda_1 = 0.2,
                   lambda_2 = 0.5,
                   lambda_3 = 1,
                   first_own_lag_prior_mean = rep(1,2),
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2),
                   SV = TRUE,
                   SV_type = "AR1",
                   SV_priors = SV_priors_AR1)

bvar_obj <- fit(bvar_obj,
                H = 8,
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1,
                control = list(max_treedepth = 12, adapt_delta = 0.85)
                )
#> ------------------------------------------------------------
#> Estimating Stan model:
#> steady_state_bvar_AR1_stochastic_volatility
#> 
#> Also generating draws from the joint predictive distribution
#> 
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
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_AR1_stochastic_volatility' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 6.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.63 seconds.
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
#> Chain 1:  Elapsed Time: 0.006 seconds (Warm-up)
#> Chain 1:                0.091 seconds (Sampling)
#> Chain 1:                0.097 seconds (Total)
#> Chain 1: 
#> Warning: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
# }
```
