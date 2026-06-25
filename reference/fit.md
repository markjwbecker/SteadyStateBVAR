# Estimate a Bayesian VAR model using Stan

Runs Stan to estimate the BVAR model using data, setup, and priors
stored in a `bvar` object. Supports standard homoscedastic steady-state
BVAR as well as stochastic volatility specifications (RW or AR(1)).

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

  Integer. Total number of MCMC iterations per chain. Default `5000`.

- warmup:

  Integer. Number of warmup (burn-in) iterations per chain. Default
  `2500`.

- chains:

  Integer. Number of MCMC chains. Default `2`.

- cores:

  Integer. Number of CPU cores used for sampling.

- auto_write:

  Logical. Whether to enable `rstan` auto-write behavior.

## Value

A `bvar` object containing:

- fit\$stan:

  Stan fit object

- posterior_means:

  Posterior means of key parameters

## Details

Forecast inputs must be provided in `x$predict$H` and `x$predict$d_pred`
before calling `fit()`.

## Examples

``` r
if (FALSE) { # \dontrun{
yt <- matrix(rnorm(40, 0, 1), 20, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1)

bvar_obj <- priors(bvar_obj,
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2))

bvar_obj$predict$H <- 1
bvar_obj$predict$d_pred <- matrix(1)

bvar_obj <- fit(bvar_obj,
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1,
                auto_write = FALSE)
} # }
```
