# Specify priors for the steady-state BVAR model

This function prepares the priors. The Minnesota prior is used for the
autoregressive parameters, and is determined by the overall tightness,
cross-equation tightness, and the lag decay rate. For the steady-state
parameters, a normal prior is used. For the covariance matrix of the
innovations, the user can choose between Jeffreys prior or an
uninformative inverse-Wishart prior. Optionally enables stochastic
volatility where the covariance matrix varies over time.

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
  Jeffreys = TRUE,
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

  Numeric vector of length `k`. Prior means for the first own lags of
  the variables. If `NULL` (default), a zero vector is used.

- theta_Psi:

  Numeric vector. Prior mean vector for \\\text{vec}(\Psi)\\, i.e. the
  steady-state parameters. If `NULL` (default), the OLS estimates are
  used.

- Omega_Psi:

  Numeric matrix. Prior covariance matrix for \\\text{vec}(\Psi)\\, i.e.
  the steady-state parameters. If `NULL` (default), a diagonal matrix
  with variances `1000` is used.

- Jeffreys:

  Logical. If `TRUE` (default), uses Jeffreys prior for the innovation
  covariance matrix. If `FALSE`, uses an uninformative inverse-Wishart
  prior. Only considered if `SV=FALSE`.

- SV:

  Logical. If `TRUE`, enables stochastic volatility specification.
  Default `FALSE`.

- SV_type:

  Character. Type of stochastic volatility model. Must be either `"RW"`
  or `"AR1"`. Required if `SV = TRUE`.

- SV_priors:

  List. User-supplied stochastic volatility priors. Required when
  `SV = TRUE`. The list must contain the following named elements
  depending on `SV_type`:

  - For `"RW"`: `theta_A`, `Omega_A`, `mu_log_lambda_1`,
    `sigma2_log_lambda_1`, `alpha_phi`, `beta_phi`.

  - For `"AR1"`: `theta_A`, `Omega_A`, `theta_gamma_0`, `Omega_gamma_0`,
    `theta_gamma_1`, `Omega_gamma_1`, `theta_log_lambda_1`,
    `Omega_log_lambda_1`, `V_Phi`, `m_Phi`.

## Value

The steady-state `bvar` object with an appended `priors` list
containing:

- theta_beta:

  Prior mean vector for \\\text{vec}(\beta)\\ constructed with the
  Minnesota prior

- Omega_beta:

  Prior covariance matrix for \\\text{vec}(\beta)\\ constructed with the
  Minnesota prior

- theta_Psi:

  Prior mean vector for \\\text{vec}(\Psi)\\, i.e. the steady-state
  parameters

- Omega_Psi:

  Prior covariance matrix for \\\text{vec}(\Psi)\\, i.e. the
  steady-state parameters

- Jeffreys:

  Indicator for Jeffreys prior usage

- Sigma_AR:

  Residual variance estimates from univariate AR fits, which are used by
  the Minnesota prior

- m:

  Inverse-Wishart prior degrees of freedom (if `Jeffreys = FALSE`)

- V:

  Inverse-Wishart prior scale matrix (if `Jeffreys = FALSE`)

- SV:

  Logical indicator for stochastic volatility specification

- SV_type:

  Stochastic volatility specification type

- SV_priors:

  User-supplied SV prior list (if `SV = TRUE`)

## Details

The goal is to estimate \\\beta, \Psi\\, and \\\Sigma_u\\, so priors are
needed. Following Villani (2009), prior independence between \\\beta,
\Psi\\ and \\\Sigma_u\\ is assumed. For \\\beta\\, i.e. the
autoregressive parameter matrix, the Minnesota prior is used

\$\$\mathrm{vec}(\beta) \sim \mathrm{N}\_{kpk}
\left\[\theta\_\beta,\Omega\_\beta\right\]\$\$

The prior means (the elements of \\\theta\_\beta\\) are set to

\$\$ \begin{aligned} \mathrm{E}\left(\Pi\_{\ell}^{(i,j)}\right)&=
\begin{cases}\kappa & \text{if } \ell = 1 \\ \text{and} \\ i = j \\0 &
\text{otherwise}\end{cases}\\ \kappa&=\begin{cases}\kappa^{level} &
\text{if } \text{variable} \\ i \\ \text{is in level} \\ \kappa^{\Delta}
& \text{if } \text{variable} \\ i \\ \text{is in difference}
\end{cases}\\ \end{aligned} \$\$

Here, the autoregressive coefficient \\\Pi\_{\ell}^{(i,j)}\\ is element
\\\left(i,j\right)\\ of \\\Pi\_{\ell}\\ for \\\ell=1,\dots,p\\. As such,
the Minnesota prior sets all prior means for the elements in \\\beta\\
to \\0\\, except for the elements that relate to the first own lags of
the variables, which are set to \\\kappa\\. If variable \\i\\ is in
level (e.g. nominal interest rate), then \\\kappa=\kappa^{level}\\, and
typical choices for \\\kappa^{level}\\ are \\1\\ or \\0.9\\. Evaluating
the equations at their prior means, equation \\i\\ becomes a random walk
if \\\kappa^{level}=1\\ and a persistent stationary AR(1) process if
\\\kappa^{level}=0.9\\. Since the steady state only exists if the
process is stationary, \\0.9\\ is recommended for the steady-state BVAR.
If variable \\i\\ is in difference (e.g. output growth), then
\\\kappa=\kappa^{\Delta}\\, and the most common choice for
\\\kappa^{\Delta}\\ is \\0\\, i.e. equation \\i\\ becomes (when
evaluating it at its prior means) a random walk expressed in first
differences. If a differenced variable still shows some degree of
persistence (can be examined with an ACF plot), a suitable value for
\\\kappa^{\Delta}\\ can be (for example) \\0.6\\ instead of \\0\\.
Moving on to the prior variances, \\\Omega\_\beta\\ is a diagonal matrix
containing the prior variances for the elements in \\\beta\\. They are
specified as

\$\$\mathrm{Var}\left(\Pi\_{\ell}^{(i,j)}\right)=
\begin{cases}\left(\frac{\lambda_1}{\ell^{\lambda_3}}\right)^2 &
\text{if } i = j \\ \left(\frac{\lambda_1
\lambda_2\sigma_i}{\ell^{\lambda_3}\sigma_j}\right)^2& \text{if } i \neq
j \end{cases}\$\$

Here \\\lambda_1\\, \\\lambda_2\\, and \\\lambda_3\\ are scalar
hyperparameters known as the overall tightness, the cross-equation
tightness and the lag decay rate. Furthermore, \\\sigma_i^2\\ is the
\\(i,i)\\:th element of \\\Sigma_u\\, which is unknown and therefore
replaced with an estimate. In this package, it is replaced by the least
squares residual variance from a univariate autoregression for variable
\\i\\ with \\p\\ lags (including the constant and dummy/trend variable
if applicable). Moving on to \\\Psi\\, the steady-state parameter
matrix, the prior is

\$\$\mathrm{vec}(\Psi) \sim
\mathrm{N}\_{kq}\left\[\theta\_\Psi,\Omega\_\Psi\right\]\$\$

This is the core of the steady-state BVAR model. In \\\theta\_\Psi\\,
one specifies the prior beliefs about the location of the steady state,
and in \\\Omega\_\Psi\\, which is assumed to be a diagonal matrix, one
specifies the degree of certainty in those prior beliefs. The prior for
\\\Sigma_u\\ is either the usual non-informative Jeffreys prior

\$\$p(\Sigma_u) \propto\left\|\Sigma_u \right\|^{-(k+1)/2}\$\$

or a proper uninformative inverse-Wishart prior

\$\$\Sigma_u \sim \mathrm{IW}(V,m)\$\$

where \\V\\ is the scale matrix and \\m\\ is the number of degrees of
freedom. An uninformative prior is specified by setting
\\V=(m-k-1)\hat{\Sigma}\_u\\ where \\\hat{\Sigma}\_u\\ is the least
squares estimate from the VAR(\\p\\) (including the constant and
dummy/trend variable if applicable), and \\m=k+2\\. For the stochastic
volatility specifications, the innovation covariance matrix is now
time-varying \\\Sigma\_{u,t}\\. Therefore, stochastic volatility priors
are needed, see
[bvar](https://markjwbecker.github.io/SteadyStateBVAR/reference/bvar.md)
for more details. For the Random Walk (`RW`) stochastic volatility
specification, the following priors are used

\$\$\begin{aligned}a &\sim \mathrm{N}(\theta_A, \Omega_A) \\ \ln
\lambda\_{i,1} &\sim \mathrm{N}(\mu\_{\ln \lambda\_{i,1}},
\sigma^2\_{\ln \lambda\_{i,1}}) \\ \phi_i &\sim
\mathrm{IG}(\alpha\_{\phi_i},\beta\_{\phi_i})\end{aligned}\$\$

Here \\a\\ is a \\k(k-1)/2\\ vector that collects the free parameters in
\\A\\ in row-major order, and \\\ln \lambda\_{i,1}\\ are the time
\\t=1\\ values (initial conditions) of \\\ln \lambda\_{i,t},
i=1,\dots,k\\. Furthermore, \\\phi_i , i=1,\dots,k\\ are the log
volatility innovation variances. For the AR(1) (`AR1`) stochastic
volatility specification, the following priors are used

\$\$\begin{aligned}a &\sim \mathrm{N}(\theta_A, \Omega_A) \\ \gamma\_{0}
&\sim \mathrm{N}(\theta\_{\gamma_0}, \Omega\_{\gamma_0}) \\ \gamma\_{1}
&\sim \mathrm{N}(\theta\_{\gamma_1}, \Omega\_{\gamma_1}) \\ \ln
\lambda\_{1} &\sim \mathrm{N}(\theta\_{\ln \lambda\_{1}}, \Omega\_{\ln
\lambda\_{1}}) \\ \Phi &\sim
\mathrm{IW}(V\_{\Phi},m\_{\Phi})\end{aligned}\$\$

Here \\a\\ is again the \\k(k-1)/2\\ vector that collects the free
parameters in \\A\\ in row-major order, and \\\ln \lambda_1\\ are the
time \\t=1\\ values (initial conditions) of \\\ln \lambda\_{t}\\.
Furthermore, \\\gamma\_{0}\\ are the log volatility intercepts,
\\\gamma\_{1}\\ are the log volatility slopes, and \\\Phi\\ is the log
volatility innovation covariance matrix.

For details on the homoscedastic steady-state BVAR model, see Villani
(2009). For details on the Random Walk stochastic volatility
steady-state BVAR model, see Clark (2011). See Carriero, Clark, and
Marcellino (2024) for the AR(1) stochastic volatility specification
applied to a conventional BVAR.

## References

Carriero, A., Clark, T. E., and Marcellino, M. (2024). Capturing
macro-economic tail risks with Bayesian vector autoregressions. *Journal
of Money, Credit and Banking*, 56(5), pp. 1099–1127.

Clark, T. E. (2011). Real-time density forecasts from Bayesian vector
autoregressions with stochastic volatility. *Journal of Business &
Economic Statistics*, 29(3), pp. 327-341.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630-650.

## Examples

``` r
#homoscedastic with Jeffreys prior
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1)

bvar_obj <- priors(bvar_obj,
                   lambda_1 = 0.2,
                   lambda_2 = 0.5,
                   lambda_3 = 1,
                   first_own_lag_prior_mean = rep(1,2),
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2),
                   Jeffreys = TRUE,
                   SV = FALSE,
                   SV_type = NULL,
                   SV_priors = NULL)
                   
#RW stochastic volatility
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1)

k <- bvar_obj$setup$k
n_free_params_A <- bvar_obj$setup$n_free_params_A

SV_priors_RW <- list(
theta_A              =  rep(0, n_free_params_A),
Omega_A              =  diag(1000, n_free_params_A),
mu_log_lambda_1      =  rep(0, k),
sigma2_log_lambda_1  =  rep(1000, k),
alpha_phi            =  rep(5, k),
beta_phi             = (rep(5, k) - 1) * rep(0.1, k)
)

bvar_obj <- priors(bvar_obj,
                   lambda_1 = 0.2,
                   lambda_2 = 0.5,
                   lambda_3 = 1,
                   first_own_lag_prior_mean = rep(1,2),
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2),
                   SV = TRUE,
                   SV_type = "RW",
                   SV_priors = SV_priors_RW)
                   
#AR1 stochastic volatility
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1)

k <- bvar_obj$setup$k
n_free_params_A <- bvar_obj$setup$n_free_params_A

SV_priors_AR1 <- list(
theta_A               =  rep(0, n_free_params_A),
Omega_A               =  diag(1000, n_free_params_A),
theta_gamma_0         =  rep(0.1, k),
Omega_gamma_0         =  diag(1000, k),
theta_gamma_1         =  rep(0.9, k),
Omega_gamma_1         =  diag(10, k),
theta_log_lambda_1    =  rep(0.1, k)/(1-rep(0.9, k)),
Omega_log_lambda_1    =  diag(1000, k),
V_Phi                 = (10 - k - 1) * diag(k),
m_Phi                 =  10
)

bvar_obj <- priors(bvar_obj,
                   lambda_1 = 0.2,
                   lambda_2 = 0.5,
                   lambda_3 = 1,
                   first_own_lag_prior_mean = rep(1,2),
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2),
                   SV = TRUE,
                   SV_type = "AR1",
                   SV_priors = SV_priors_AR1)
```
