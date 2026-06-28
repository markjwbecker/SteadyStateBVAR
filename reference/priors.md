# Specify priors for the steady-state BVAR model

The function prepares the priors. The Minnesota prior is used for the
autoregressive parameters, and is determined by the overall tightness,
cross-equation tightness, and the lag decay rate. For the steady-state
parameters, a normal prior is used. For the covariance matrix of the
innovations, the user can choose between Jeffreys prior or an
uninformative inverse-Wishart prior. Optionally enables a stochastic
volatility specification for the covariance matrix of the innovations
(random walk or AR(1)).

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

  A steady-state `bvar` object that has been passed through
  [`setup`](https://markjwbecker.github.io/SteadyStateBVAR/reference/setup.md).

- lambda_1:

  Numeric. Overall tightness of the Minnesota prior. Default `0.2`.

- lambda_2:

  Numeric. Cross-equation tightness of the Minnesota prior. Default
  `0.5`.

- lambda_3:

  Numeric. Lag decay rate of the Minnesota prior. Default `1`.

- first_own_lag_prior_mean:

  Numeric vector of length `k`. Prior means for the first own lag of
  each variable. If `NULL`, defaults to a zero vector.

- theta_Psi:

  Numeric vector. Prior mean for vec(Psi), i.e. the steady-state
  parameters. If `NULL`, defaults to the OLS estimates.

- Omega_Psi:

  Numeric matrix. Prior covariance matrix for vec(Psi), i.e. the
  steady-state parameters. If `NULL`, defaults to a diagonal matrix with
  variances `1000`.

- Jeffrey:

  Logical. If `TRUE`, uses a Jeffreys prior for the innovation
  covariance matrix. If `FALSE`, uses an uninformative inverse-Wishart
  prior.

- SV:

  Logical. If `TRUE`, enables stochastic volatility specification.
  Default `FALSE`.

- SV_type:

  Character. Type of stochastic volatility model. Must be either `"RW"`
  or `"AR1"`. Required if `SV = TRUE`.

- SV_priors:

  List. User-supplied stochastic volatility priors. Required when
  `SV = TRUE`.

## Value

The steady-state `bvar` object with an appended `priors` list
containing:

- theta_beta:

  Prior mean vector for vec(beta) constructed with the Minnesota prior

- Omega_beta:

  Prior covariance matrix for vec(beta) constructed with the Minnesota
  prior

- theta_Psi:

  Prior mean for vec(Psi), i.e. the steady-state parameters

- Omega_Psi:

  Prior covariance matrix for vec(Psi), i.e. the steady-state parameters

- Jeffrey:

  Indicator for Jeffreys prior usage

- Sigma_AR:

  Residual variance estimates from univariate AR fits, which are used by
  the Minnesota prior

- m_0:

  Inverse-Wishart prior degrees of freedom (if `Jeffrey = FALSE`)

- V_0:

  Inverse-Wishart prior scale matrix (if `Jeffrey = FALSE`)

- SV:

  Logical indicator for stochastic volatility specification

- SV_type:

  Stochastic volatility specification type

- SV_priors:

  User-supplied SV prior list (if `SV = TRUE`)

## Examples

``` r
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1)

bvar_obj <- priors(bvar_obj,
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2))
```
