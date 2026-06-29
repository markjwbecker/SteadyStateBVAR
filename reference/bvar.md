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
d\_{t-1})+\dots+\Pi_p(y\_{t-p}-\Psi d\_{t-p})+u_t\$\$

where \\y_t\\ is an \\k\\-dimensional vector of endogenous variables at
time \\t\\, and \\d_t\\ is a \\q\\-dimensional vector of deterministic
(exogenous) variables at time \\t\\. Here \\\Pi\_\ell\\ for
\\\ell=1,\dots,p\\ is a \\(k \times k)\\ matrix of autoregressive
parameters, and \\\Psi\\ is a \\(k \times q)\\ matrix of steady-state
parameters. Note that \\\mathrm{E}(y_t)=\mu_t=\Psi d_t\\ is the
unconditional mean, or the **steady state** of the process. One can
stack the (transposed) \\\Pi_i\\ matrices in the \\(kp \times k)\\
matrix \\\beta\\ \$\$\beta=\begin{bmatrix}\Pi'\_1 \\ \vdots
\\\Pi'\_p\end{bmatrix}\$\$ Then the model can be rewritten as a
nonlinear regression (Karlsson, 2013) \$\$y_t' =d_t'\Psi' +
\left\[w_t'-q_t'(I_p \otimes \Psi') \right\]\beta +u_t'\$\$ where where
\\w_t'=(y\_{t-1}',\dots,y\_{t-p}')\\ is a \\kp\\-dimensional vector of
lagged endogenous variables and \\q_t'=(d\_{t-1}',\dots,d\_{t-p}')\\ is
a \\qp\\-dimensional vector of lagged deterministic (exogenous)
variables, \\I_p\\ is the \\(p \times p)\\ identity matrix and
\\\otimes\\ denotes the Kronecker product. This is how the likelihood is
written in the Stan code. The goal is to estimate \\\beta, \Psi\\, and
\\\Sigma_u\\.

For the innovations to the model, in the case of the homoscedastic
steady-state BVAR, they are \\u_t \sim \mathrm{N_k}(0,\Sigma_u)\\.
However, for models with stochastic volatility, there is instead a
time-varying innovation covariance matrix \\u_t \sim
\mathrm{N_k}(0,\Sigma\_{u,t})\\. The innovations then take the form

\$\$\begin{aligned} u_t &= A^{-1} \Lambda^{0.5}\_t \epsilon_t \\
\epsilon_t &\sim \mathrm{N}(0, \mathrm{I}\_k)\end{aligned}\$\$

where \\A\\ is a lower triangular matrix with ones on the diagonal that
describes the contemporaneous interaction of the endogenous variables,
and

\$\$\Lambda_t = \mathrm{diag}(\lambda\_{1,t},\dots,\lambda\_{k,t})\$\$

contains the time-varying variances (log volatilities) of conditionally
Gaussian shocks. For the `AR1` stochastic volatility specification, the
log volatilities follow AR(1) processes

\$\$\ln \lambda\_{i,t} = \gamma\_{0,i} + \gamma\_{1,i} \ln
\lambda\_{i,t-1} + \nu\_{i,t}, \\ i=1,\dots,k\$\$

where the log volatility AR(1) processes are restricted to the
stationary region, i.e. \\\|\gamma\_{1,i}\|\<1 \\ \forall i\\. For the
`RW` stochastic volatility specification, the log volatilities follow
Random Walk processes

\$\$\gamma\_{0,i}=0, \\ \gamma\_{1,i}=1 \\ \forall i\$\$

The innovations to the log volatilities follow in the `AR1` case

\$\$\nu\_{t} = (\nu\_{1,t},\dots,\nu\_{k,t})'\sim \mathrm{N}(0,
\Phi)\$\$

where \\\Phi\\ is **not** diagonal and as such the innovations to the
log volatilities are allowed to be correlated across variables. For the
`RW` case, \\\Phi\\ is diagonal with variances \\\phi_i\\ for
\\i=1,\dots,k\\.

Note that under both stochastic volatility specifications, the
time-varying covariance matrix is

\$\$\Sigma\_{u,t} = A^{-1} \Lambda_t (A^{-1})'\$\$

For details on the homoscedastic steady-state BVAR model, see Villani
(2009). For the Random Walk stochastic volatility steady-state BVAR
model, see Clark (2011). For details regarding AR(1) stochastic
volatility, see Carriero, Clark and Marcellino (2024).

## References

Carriero, A., Clark, T. E., and Marcellino, M. (2024). Capturing
Macro‐Economic Tail Risks with Bayesian Vector *Journal of Money, Credit
and Banking*. 56(5), pp. 1099–1127.

Clark, T. E. (2011). Real-Time Density Forecasts from Bayesian Vector
Autoregressions with Stochastic Volatility. *Journal of Business &
Economic Statistics*. 29(3), pp. 327–341.

Karlsson, S. (2013). Forecasting with Bayesian Vector Autoregression.
In: Elliott, G. and Timmerman, A. (eds) *Handbook of Economic
Forecasting*. Elsevier B.V. Vol 2, Part B., pp. 791-897.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*. 24(4), pp. 630–649.

## Examples

``` r
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)
```
