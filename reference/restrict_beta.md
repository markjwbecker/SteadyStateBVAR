# Restrict VAR coefficients to zero

Applies zero restrictions to the VAR coefficient matrix by setting the
corresponding prior variances in `Omega_beta` to a value near zero. This
enforces restrictions through the prior rather than hard-coding them in
the likelihood, which is compatible with the Stan estimation.

## Usage

``` r
restrict_beta(x, restriction_matrix)
```

## Arguments

- x:

  A `bvar` object that has been passed through
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
yt <- matrix(rnorm(40, 0, 1), 20, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1)

bvar_obj <- priors(bvar_obj,
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2))
                   
p <- bvar_obj$setup$p
k <- bvar_obj$setup$k
restriction_matrix <- matrix(1, k*p, k)
restriction_matrix[1, 1] <- 0

bvar_obj <- restrict_beta(bvar_obj, restriction_matrix)
```
