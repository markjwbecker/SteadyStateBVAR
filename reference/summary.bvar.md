# Summarise a fitted BVAR model

Computes posterior means of the model parameters from a fitted `bvar`
object. The parameters returned depend on the model specification:
standard homoscedastic models return `beta`, `Psi`, and `Sigma`;
stochastic volatility models additionally return volatility parameters.

## Usage

``` r
# S3 method for class 'bvar'
summary(object, pars = NULL, ...)

# S3 method for class 'summary.bvar'
print(x, ...)
```

## Arguments

- object:

  A `bvar` object that has been passed through
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md).

- pars:

  Character vector of parameter names to include. If `NULL` (default),
  all parameters are returned.

- ...:

  Further arguments passed to or from other methods.

- x:

  A `summary.bvar` object returned by `summary.bvar`.

## Value

An object of class `summary.bvar`.

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
                
summary(bvar_obj)
} # }
```
