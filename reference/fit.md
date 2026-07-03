# Estimate the steady-state BVAR model using Stan

Estimates the steady-state BVAR model using the No-U-Turn sampler (a
variant of Hamiltonian Monte Carlo) via Stan. Uses the data, setup, and
priors stored in the steady-state `bvar` object. Supports both
homoscedastic and stochastic volatility (RW or AR1) specifications.

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

  Matrix of size H x q. Future values of the deterministic variables
  d_t. Default is `NULL`, must be provided by the user.

- ...:

  Additional arguments passed directly to
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
prior settings:

- Homoscedastic with Jeffreys prior:
  `steady_state_bvar_homoscedastic_jeffreys_prior`

- Homoscedastic with uninformative inverse-Wishart prior:
  `steady_state_bvar_homoscedastic_inverse_wishart_prior`

- Random Walk stochastic volatility:
  `steady_state_bvar_RW_stochastic_volatility`

- AR1 stochastic volatility:
  `steady_state_bvar_AR1_stochastic_volatility`

The function estimates the following parameters (see
[bvar](https://markjwbecker.github.io/SteadyStateBVAR/reference/bvar.md)
for details):

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
                cores = 1)
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_homoscedastic_jeffreys_prior' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 7.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.73 seconds.
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
#> Chain 1:  Elapsed Time: 0.045 seconds (Warm-up)
#> Chain 1:                0.105 seconds (Sampling)
#> Chain 1:                0.15 seconds (Total)
#> Chain 1: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

#RW stochastic volatility
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")

k <- bvar_obj$setup$k
n_free_params_A <- bvar_obj$setup$n_free_params_A

SV_priors_RW <- list(
theta_A             =  rep(0, n_free_params_A),
Omega_A             =  diag(1000, n_free_params_A),
mu_log_lambda_0     =  rep(0, k),
sigma2_log_lambda_0 =  rep(1000, k),
alpha_phi           =  rep(5, k),
beta_phi            = (rep(5, k) - 1) * rep(0.1, k)
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
                d_pred = matrix(rep(1,8)),
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1,
                control = list(max_treedepth = 12, adapt_delta = 0.85)
                )
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_RW_stochastic_volatility' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000212 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.12 seconds.
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
#> Chain 1:  Elapsed Time: 1.45 seconds (Warm-up)
#> Chain 1:                5.683 seconds (Sampling)
#> Chain 1:                7.133 seconds (Total)
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

#AR1 stochastic volatility
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")

k <- bvar_obj$setup$k
n_free_params_A <- bvar_obj$setup$n_free_params_A

SV_priors_AR <- list(
theta_A            =  rep(0, n_free_params_A),
Omega_A            =  diag(1000, n_free_params_A),
theta_gamma_0      =  rep(0.1, k),
Omega_gamma_0      =  diag(1000, k),
theta_gamma_1      =  rep(0.9, k),
Omega_gamma_1      =  diag(10, k),
theta_log_lambda_0 =  rep(0.1, k)/(1-rep(0.9, k)),
Omega_log_lambda_0 =  diag(1000, k),
V_0                = (10 - k - 1) * diag(k),
m_0                =  10
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
                   SV_priors = SV_priors_AR)

bvar_obj <- fit(bvar_obj,
                H = 8,
                d_pred = matrix(rep(1,8)),
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1,
                control = list(max_treedepth = 12, adapt_delta = 0.85)
                )
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_AR1_stochastic_volatility' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.00023 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.3 seconds.
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
#> Chain 1:  Elapsed Time: 2.042 seconds (Warm-up)
#> Chain 1:                0.197 seconds (Sampling)
#> Chain 1:                2.239 seconds (Total)
#> Chain 1: 
#> Warning: There were 150 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
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
