# Estimate the steady-state BVAR model using Stan

Estimates the steady-state BVAR model using the No-U-Turn sampler (a
variant of Hamiltonian Monte Carlo) via Stan. Uses the data, setup, and
priors stored in the steady-state `bvar` object. Supports both
homoscedastic and stochastic volatility (RW or AR1) specifications.

## Usage

``` r
fit(
  x,
  H = 1,
  d_pred = NULL,
  iter = 2000,
  warmup = floor(iter/2),
  chains = 1,
  cores = getOption("mc.cores", 1L),
  verbose = FALSE,
  auto_write = FALSE
)
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

- iter:

  Integer. Total number of MCMC iterations per chain. Default is `2000`.

- warmup:

  Integer. Number of warmup (burn-in) iterations per chain. Default is
  `floor(iter/2)`.

- chains:

  Integer. Number of MCMC chains. Default is `1`.

- cores:

  Positive integer specifying the number of CPU cores used for sampling.
  Default is `getOption("mc.cores", 1L)`, i.e. the `mc.cores` option (if
  it has been set), otherwise defaults to `1` core.

- verbose:

  Logical indicating whether to print intermediate output from Stan on
  the console, defaults to `FALSE`.

- auto_write:

  Logical indicating whether Stan models should be automatically written
  to the disk cache via `rstan`. Default is `FALSE`.

## Value

A fitted steady-state `bvar` object with:

- `fit$stan`: An object of class `stanfit`

- `fit$posterior_means`: List of posterior mean estimates

- `fit$posterior_medians`: List of posterior median estimates

## Details

The function selects the appropriate Stan model based on prior settings:

- Homoscedastic with Jeffreys prior:
  `steady_state_bvar_homoscedastic_jeffreys_prior.stan`

- Homoscedastic with uninformative inverse-Wishart prior:
  `steady_state_bvar_homoscedastic_inverse_wishart_prior.stan`

- Random Walk stochastic volatility:
  `steady_state_bvar_RW_stochastic_volatility.stan`

- AR1 stochastic volatility:
  `steady_state_bvar_AR1_stochastic_volatility.stan`

The function estimates the following parameters:

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
if (FALSE) { # \dontrun{
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p = 1, deterministic = "constant")

bvar_obj <- priors(bvar_obj,
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2))

bvar_obj <- fit(bvar_obj,
                H = 8,
                d_pred = matrix(rep(1,8)),
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1,
                verbose = FALSE,
                auto_write = FALSE)
} # }
```
