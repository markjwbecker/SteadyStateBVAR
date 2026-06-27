# Estimate a steady-state BVAR model using Stan

Runs Stan to estimate a Bayesian VAR with steady-state formulation using
data, setup, and priors stored in a `bvar` object. Supports both
homoscedastic specifications and stochastic volatility models (RW or
AR1).

## Usage

``` r
fit(
  x,
  iter = 5000,
  warmup = 2500,
  chains = 2,
  cores = min(chains, parallel::detectCores()),
  auto_write = TRUE
)
```

## Arguments

- x:

  A `bvar` object that has been passed through
  [`setup`](https://markjwbecker.github.io/SteadyStateBVAR/reference/setup.md)
  and
  [`priors`](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md).

- iter:

  Integer. Total number of MCMC iterations per chain. Default is `5000`.

- warmup:

  Integer. Number of warmup (burn-in) iterations per chain. Default is
  `2500`.

- chains:

  Integer. Number of MCMC chains. Default is `2`.

- cores:

  Integer. Number of CPU cores used for sampling. Default is
  `min(chains, available cores)`.

- auto_write:

  Logical. Whether to enable `rstan` auto-write. Default is `TRUE`.

## Value

A `bvar` object with:

- `fit$stan`: Stan fit object containing posterior draws

- `posterior_means`: List of posterior mean estimates:

  - `beta`: k×k VAR coefficient matrix

  - `Psi`: k×q steady-state parameter matrix

  - `Sigma_u`: covariance matrix (k×k for homoscedastic, T×k×k for
    stochastic volatility)

  - RW SV: `A`, `phi`

  - AR1 SV: `A`, `gamma_0`, `gamma_1`, `Phi`

## Details

Forecast inputs must be supplied in `x$predict$H` and `x$predict$d_pred`
prior to calling `fit()`.

The function selects the appropriate Stan model based on prior settings:

- Homoscedastic + Jeffreys prior:
  `steady_state_bvar_homoscedastic_jeffreys_prior.stan`

- Homoscedastic + inverse-Wishart prior:
  `steady_state_bvar_homoscedastic_inverse_wishart_prior.stan`

- Stochastic volatility (RW):
  `steady_state_bvar_RW_stochastic_volatility.stan`

- Stochastic volatility (AR(1)):
  `steady_state_bvar_AR1_stochastic_volatility.stan`

For stochastic volatility models, SV-specific priors in
`priors$SV_priors` are merged into the Stan data list.

## Examples

``` r
if (FALSE) { # \dontrun{
yt <- matrix(rnorm(40), 20, 2)
bvar_obj <- bvar(data = yt)
bvar_obj <- setup(bvar_obj, p = 1)
bvar_obj <- priors(bvar_obj,
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2))
bvar_obj$predict$H <- 1
bvar_obj$predict$d_pred <- matrix(1)

bvar_obj <- fit(bvar_obj, iter = 200, warmup = 50,
                chains = 1, cores = 1, auto_write = FALSE)

summary(bvar_obj)
} # }
```
