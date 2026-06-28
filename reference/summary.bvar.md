# Summarise a fitted steady-state BVAR model

Computes and prints posterior summaries from a fitted `bvar` object. The
printed output depends on whether the model is homoscedastic or includes
stochastic volatility (RW or AR1 specification).

## Usage

``` r
# S3 method for class 'bvar'
summary(object, pars = NULL, stat = "mean", t = NULL, ...)
```

## Arguments

- object:

  A `bvar` object that has been passed through
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md).

- pars:

  Character vector of parameter names to include. If `NULL` (default),
  all available parameters are displayed.

- stat:

  Character. Posterior summary statistic to display: `"mean"` (default)
  or `"median"`.

- t:

  Integer. Time index for stochastic volatility covariance matrices. If
  `NULL`, the latest available period is used.

- ...:

  Additional arguments (currently unused).

## Value

Returns the input object invisibly.

## Details

Parameters are printed in structured blocks. Matrix-valued parameters
are displayed as matrices, while vector-valued stochastic volatility
parameters (e.g. `phi`, `gamma_0`, `gamma_1`) are printed as named
vectors.

For stochastic volatility models, time-varying covariance matrices
`Sigma_u,t` can be extracted at a specified time index.

The function supports both standard VAR models and stochastic volatility
extensions. For SV models, additional parameters may be displayed
depending on specification:

- RW SV: `phi`

- AR1 SV: `gamma_0`, `gamma_1`, `Phi`

lik

Output is printed in blocks with manual formatting for readability.

## Examples

``` r
if (FALSE) { # \dontrun{
yt <- matrix(rnorm(50), 25, 2)
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
