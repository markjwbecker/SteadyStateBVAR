# SteadyStateBVAR: Bayesian Vector Autoregressions with Steady-State Priors

Provides estimation of Bayesian vector autoregression (BVAR) models with
steady-state priors via 'Stan', along with functions for unconditional
and conditional forecasting, as well as impulse response analysis. For
details on the steady-state BVAR model see Villani (2009)
[doi:10.1002/jae.1065](https://doi.org/10.1002/jae.1065) .

## Details

The steady-state BVAR model takes the form

\$\$y_t = \Psi d_t + \Pi_1(y\_{t-1}-\Psi
d\_{t-1})+\dots+\Pi_p(y\_{t-p}-\Psi d\_{t-p})+u_t\$\$

where \\y_t\\ is an \\k\\-dimensional vector of endogenous variables at
time \\t\\, and \\d_t\\ is a \\q\\-dimensional vector of deterministic
(exogenous) variables at time \\t\\. Here \\\Pi\_\ell\\ for
\\\ell=1,\dots,p\\ is a \\(k \times k)\\ matrix of autoregressive
parameters, and \\\Psi\\ is a \\(k \times q)\\ matrix of steady-state
parameters. Now

\$\$\mathrm{E}(y_t)=\mu_t=\Psi d_t\$\$

is the unconditional mean, or the **steady state** of the process.

Regarding the (reduced-form) innovations \\u_t\\, this package supports
three specifications

1.  \\u_t \overset{\text{iid}}{\sim} \mathrm{N_k}(0,\Sigma_u)\\

2.  \\u_t \sim \mathrm{N_k}(0,\Sigma\_{u,t})\\, where \\\Sigma\_{u,t}\\
    is driven by a latent Random Walk process (stochastic volatility)

3.  \\u_t \sim \mathrm{N_k}(0,\Sigma\_{u,t})\\, where \\\Sigma\_{u,t}\\
    is driven by a latent AR(1) process (stochastic volatility)

For more information, please see the package vignettes

1.  [`vignette("SteadyStateBVAR-intro")`](https://markjwbecker.github.io/SteadyStateBVAR/articles/SteadyStateBVAR-intro.md)

2.  [`vignette("Homoscedastic-steady-state-BVAR")`](https://markjwbecker.github.io/SteadyStateBVAR/articles/Homoscedastic-steady-state-BVAR.md)

3.  [`vignette("RW-stochastic-volatility-steady-state-BVAR")`](https://markjwbecker.github.io/SteadyStateBVAR/articles/RW-stochastic-volatility-steady-state-BVAR.md)

4.  [`vignette("AR1-stochastic-volatility-steady-state-BVAR")`](https://markjwbecker.github.io/SteadyStateBVAR/articles/AR1-stochastic-volatility-steady-state-BVAR.md)

## References

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630–650.

## See also

Useful links:

- <https://github.com/markjwbecker/SteadyStateBVAR>

- <https://markjwbecker.github.io/SteadyStateBVAR/>

- Report bugs at
  <https://github.com/markjwbecker/SteadyStateBVAR/issues>

## Author

**Maintainer**: Mark Becker <mark.jw.becker@gmail.com> \[copyright
holder\]

Authors:

- Mark Becker <mark.jw.becker@gmail.com> \[copyright holder\]
