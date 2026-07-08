# Restrict VAR coefficients to zero

Applies zero restrictions to the VAR coefficient matrix \\\beta\\ by
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

The steady-state `bvar` object with the restriction matrix stored in
`setup` and the prior covariance matrix `Omega_beta` updated
accordingly.

## Details

The steady-state BVAR model takes the form

\$\$y_t = \Psi d_t + \Pi_1(y\_{t-1}-\Psi
d\_{t-1})+\dots+\Pi_p(y\_{t-p}-\Psi d\_{t-p})+u_t\$\$

where \\y_t\\ is an \\k\\-dimensional vector of endogenous variables at
time \\t\\, and \\d_t\\ is a \\q\\-dimensional vector of deterministic
(exogenous) variables at time \\t\\. Here \\\Pi\_\ell\\ for
\\\ell=1,\dots,p\\ is a \\(k \times k)\\ matrix of autoregressive
parameters, and \\\Psi\\ is a \\(k \times q)\\ matrix of steady-state
parameters. One can stack the (transposed) \\\Pi_i\\ matrices in the
\\(kp \times k)\\ matrix \\\beta\\ \$\$\beta=\begin{bmatrix}\Pi'\_1 \\
\vdots \\\Pi'\_p\end{bmatrix}\$\$

This function puts zero restrictions on the elements of \\\beta\\.

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

restriction_matrix <- matrix(1,k*p,k, dimnames=list(NULL,c("y1","y2")))

#restrict beta so y2 does not granger-cause y1

restriction_matrix[1, 2] <- 0
restriction_matrix[3, 2] <- 0

print(restriction_matrix)
#>      y1 y2
#> [1,]  1  0
#> [2,]  1  1
#> [3,]  1  0
#> [4,]  1  1

bvar_obj <- restrict_beta(bvar_obj, restriction_matrix)
```
