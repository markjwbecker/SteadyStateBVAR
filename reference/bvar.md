# Create a steady-state BVAR model object

Initialises a steady-state `bvar` object. This is the starting point for
all models in `SteadyStateBVAR`. After creation, pass the object
sequentially to
[`setup`](https://markjwbecker.github.io/SteadyStateBVAR/reference/setup.md),
[`priors`](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md),
and
[`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md)
to build and estimate the model.

## Usage

``` r
bvar(data)
```

## Arguments

- data:

  A numeric matrix or time series of data where each column is a
  variable and each row is a time period.

## Value

An object of class `bvar`.

## Details

The model takes the form

\$\$y_t = \Psi d_t + \Pi_1(y\_{t-1}-\Psi
d\_{t-1})+\dots+\Pi_p(y\_{t-p}-\Psi d\_{t-p})+u_t,\$\$

where \\y_t\\ is an \\k\\-dimensional vector of endogenous variables at
time \\t\\, and \\d_t\\ is a \\q\\-dimensional vector of deterministic
(exogenous) variables at time \\t\\. Here \\\Pi\_\ell\\ for
\\\ell=1,\dots,p\\ is a \\(k \times k)\\ matrix of autoregressive
parameters, and \\\Psi\\ is a \\(k \times q)\\ matrix of steady-state
parameters. Note that \\\mathrm{E}(y_t)=\mu_t=\Psi d_t\\ is the
unconditional mean, or the **steady state** of the process. In the case
of the homoscedastic steady-state BVAR, we have \\u_t \sim
N_k(0,\Sigma_u)\\, and for the models with stochastic volatility we have
\\u_t \sim N_k(0,\Sigma\_{u,t})\\.

## References

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*. 24(4), pp. 630–649.

## Examples

``` r
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)
```
