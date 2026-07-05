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
  variables \\d_t\\. Default is `NULL`, must be provided by the user.

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
                   Jeffrey = TRUE)
                   
H <- 8
d_pred <- matrix(rep(1,8))
colnames(d_pred) <- c("constant")
rownames(d_pred) <- paste("Horizon", 1:H)
print(d_pred)
#>           constant
#> Horizon 1        1
#> Horizon 2        1
#> Horizon 3        1
#> Horizon 4        1
#> Horizon 5        1
#> Horizon 6        1
#> Horizon 7        1
#> Horizon 8        1

bvar_obj <- fit(bvar_obj,
                H = H,
                d_pred = d_pred,
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1)
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_homoscedastic_jeffreys_prior' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.5e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.45 seconds.
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
#> Chain 1:  Elapsed Time: 0.023 seconds (Warm-up)
#> Chain 1:                0.061 seconds (Sampling)
#> Chain 1:                0.084 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is 1.05, indicating chains have not mixed.
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
                   Jeffrey = FALSE) #inverse-Wishart

H <- 8
d_pred <- cbind(rep(1, H), rep(0, H))
colnames(d_pred) <- c("constant", "dummy")
rownames(d_pred) <- paste("Horizon", 1:H)
print(d_pred)
#>           constant dummy
#> Horizon 1        1     0
#> Horizon 2        1     0
#> Horizon 3        1     0
#> Horizon 4        1     0
#> Horizon 5        1     0
#> Horizon 6        1     0
#> Horizon 7        1     0
#> Horizon 8        1     0

bvar_obj <- fit(bvar_obj,
                H = H,
                d_pred = d_pred,
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1)
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_homoscedastic_inverse_wishart_prior' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.37 seconds.
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
#> Chain 1:  Elapsed Time: 0.021 seconds (Warm-up)
#> Chain 1:                0.047 seconds (Sampling)
#> Chain 1:                0.068 seconds (Total)
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
#> Chain 1: Gradient evaluation took 7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.7 seconds.
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
#> Chain 1:  Elapsed Time: 0.018 seconds (Warm-up)
#> Chain 1:                0.038 seconds (Sampling)
#> Chain 1:                0.056 seconds (Total)
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
#> Chain 1: Rejecting initial value:
#> Chain 1:   Error evaluating the log probability at the initial value.
#> Chain 1: Exception: multi_normal_lpdf: LDLT_Factor of covariance parameter is not positive definite.  last conditional variance is -1.69407e-21. (in 'steady_state_bvar_AR1_stochastic_volatility', line 102, column 6 to column 54)
#> Chain 1: 
#> Chain 1: Gradient evaluation took 6.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.61 seconds.
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
#> Chain 1:  Elapsed Time: 0.573 seconds (Warm-up)
#> Chain 1:                1.815 seconds (Sampling)
#> Chain 1:                2.388 seconds (Total)
#> Chain 1: 
#> Warning: There were 57 divergent transitions after warmup. See
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
