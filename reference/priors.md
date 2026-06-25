# Specify priors for a Bayesian VAR model with optional stochastic volatility

Defines Minnesota-style priors for VAR coefficients and priors for
steady-state parameters. Optionally enables a stochastic volatility
specification (random walk or AR(1)) which augments the prior object
with additional SV-related inputs used by the Stan model. A Jeffrey's
prior can be used for the error covariance matrix or, alternatively, an
inverse-Wishart prior.

## Usage

``` r
priors(
  x,
  lambda_1 = 0.2,
  lambda_2 = 0.5,
  lambda_3 = 1,
  first_own_lag_prior_mean = NULL,
  theta_Psi = NULL,
  Omega_Psi = NULL,
  Jeffrey = TRUE,
  SV = FALSE,
  SV_type = NULL,
  SV_priors = NULL
)
```

## Arguments

- x:

  A `bvar` object that has been passed through
  [`setup`](https://markjwbecker.github.io/SteadyStateBVAR/reference/setup.md).

- lambda_1:

  Numeric. Overall tightness of the Minnesota prior controlling overall
  shrinkage. Default `0.2`.

- lambda_2:

  Numeric. Cross-variable shrinkage parameter controlling shrinkage of
  other variables' lags relative to own lags. Default `0.5`.

- lambda_3:

  Numeric. Lag decay parameter controlling how prior variance decreases
  with lag length. Default `1`.

- first_own_lag_prior_mean:

  Numeric vector of length `k`. Prior means for first own lag of each
  variable. If `NULL`, defaults to a zero vector.

- theta_Psi:

  Numeric vector. Prior mean for steady-state parameters. If `NULL`,
  defaults to a zero vector of length `k*q`.

- Omega_Psi:

  Numeric matrix. Prior covariance for steady-state parameters. If
  `NULL`, defaults to a diagonal matrix with variance 100.

- Jeffrey:

  Logical. If `TRUE`, uses a Jeffrey's prior for the error covariance
  matrix. If `FALSE`, uses an inverse-Wishart prior.

- SV:

  Logical. If `TRUE`, enables stochastic volatility specification.
  Default `FALSE`.

- SV_type:

  Character. Type of stochastic volatility model. Must be either `"RW"`
  or `"AR"`. Required if `SV = TRUE`.

- SV_priors:

  List. User-supplied stochastic volatility prior objects passed
  directly to Stan. Required when `SV = TRUE`.

## Value

The `bvar` object with an appended `priors` list containing:

- theta_beta:

  Prior mean for VAR coefficients

- Omega_beta:

  Prior covariance for VAR coefficients

- theta_Psi:

  Prior mean for steady-state parameters

- Omega_Psi:

  Prior covariance for steady-state parameters

- Jeffrey:

  Indicator for Jeffrey's prior usage

- Sigma_AR:

  Residual variance estimates from univariate AR fits

- m_0:

  Inverse-Wishart degrees of freedom (if `Jeffrey = FALSE`)

- V_0:

  Inverse-Wishart scale matrix (if `Jeffrey = FALSE`)

- SV:

  Logical indicator for stochastic volatility model

- SV_type:

  Stochastic volatility specification type

- SV_priors:

  User-supplied SV prior list (if `SV = TRUE`)

## Validation

The function performs the following checks:

- Ensures `x` is a valid `bvar` object

- Ensures `setup` exists and contains required elements

- Ensures lambda parameters are strictly positive

- Ensures `first_own_lag_prior_mean` has length `k` (if provided)

- Ensures `theta_Psi` and `Omega_Psi` dimensions match (if both
  provided)

- Ensures valid `SV_type` when `SV = TRUE`

- Ensures `SV_priors` is provided when `SV = TRUE`

## Examples

``` r
yt <- matrix(rnorm(40, 0, 1), 20, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1)

bvar_obj <- priors(bvar_obj,
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2))
```
