# Restrict VAR coefficients to zero

Applies zero restrictions to the VAR coefficient matrix (beta) by
setting the corresponding prior variances in `Omega_beta` to a value
near zero (prior means are zero by default in the Minnesota prior).

## Usage

``` r
restrict_beta(x, restriction_matrix)
```

## Arguments

- x:

  A steady-state `bvar` object that has been passed through
  [`priors`](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md).

- restriction_matrix:

  A numeric matrix of dimension `k*p x k` where entries of `0` indicate
  coefficients to be restricted to zero and entries of `1` indicate
  unrestricted coefficients.

## Value

The `bvar` object with the restriction matrix stored in `setup` and
`Omega_beta` updated accordingly.

## Examples

``` r
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=2)

bvar_obj <- priors(bvar_obj,
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2))
                   
p <- bvar_obj$setup$p
k <- bvar_obj$setup$k

restriction_matrix <- matrix(1, k*p, k)

restriction_matrix[1, 1] <- 0
restriction_matrix[4, 2] <- 0

bvar_obj <- restrict_beta(bvar_obj, restriction_matrix)
```
