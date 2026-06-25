# SteadyStateBVAR

The goal of SteadyStateBVAR is to …

## Installation

You can install the development version of SteadyStateBVAR from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("markjwbecker/SteadyStateBVAR")
```

## Introduction

With this package, the user can estimate the steady-state BVAR(p) model
of Mattias Villani (Villani, 2009). The goal is to use modern Bayesian
tools (Stan) to: i) estimate the model as specified in the original
paper, and ii) extend the model in many different ways. Back in the day,
extensions of the model seemed to be limited by what Mattias Villani had
time to derive.

See for example, Clark (2011): “*In a methodological sense, this paper
extends the estimator of Villani (2009) to include stochastic
volatility.*” Later on, he writes: “*(Special thanks are due to Mattias
Villani for providing the formulas for posterior means and variances of*
$`\Pi`$*and* $`\Psi`$*, which generalize the constant-variance formulas
of Villani 2009.)*”

For a more recent example of the model being at the mercy of Professor
Villani, we can take a look at the [Technical
Guide](https://github.com/european-central-bank/BEAR-toolbox/blob/master/tbx/doc/Technical%20guide.pdf)
for the BEAR (Bayesian Estimation, Analysis and Regression) toolbox by
the ECB. On page 122, Dieppe, Legrand, and van Roye (2018) write:
“*Villani (2009) only provides derivation in the case of the
normal-diffuse prior distribution, so the incoming analysis will be
restricted to this case.*”

Times are different now, and with the help of Stan, we can basically do
whatever we can imagine. To showcase this, in the last section we build
on the steady-state BVAR model of Clark (2011) \[which itself is an
extension of the original model (Villani, 2009)\], by letting the log
volatilities follow correlated AR(1) processes instead of uncorrelated
driftless random walks. And this is done without asking Professor
Villani for any derivations. I simply wrote down the model equations,
put some priors on the parameters, and Stan did the rest!

# Introduction

The steady-state BVAR($`p`$) model (Villani, 2009) is

``` math
y_t = \Psi d_t + \Pi_1(y_{t-1}-\Psi d_{t-1})+\dots+\Pi_p(y_{t-p}-\Psi d_{t-p})+u_t
```

where $`y_t`$ is a $`k`$-dimensional vector of endogenous variables
(time series) at time $`t`$, $`d_t`$ is a $`q`$-dimensional vector of
deterministic (exogenous) variables at time $`t`$, and the
(reduced-form) innovations are $`u_t \sim N_k(0,\Sigma_u)`$ with
independence between time periods. Here $`\Pi_\ell`$ for
$`\ell=1,\dots,p`$ is a $`(k \times k)`$ matrix, and $`\Psi`$ is a
$`(k \times q)`$ matrix. Now

``` math
\mathrm{E}(y_t)=\mu_t=\Psi d_t
```

is the unconditional mean, or the **steady state**, of the process.
Since long-horizon forecasts from stationary VARs converge to the
unconditional mean (steady state), it is naturally very important from a
forecasting perspective to obtain precise inference on $`\Psi`$.

Note that the current version of this package only allows for $`d_t`$ to
contain either a constant, a constant and a dummy variable, or a
constant and a time trend. We can stack the (transposed) $`\Pi_i`$
matrices in the $`(kp \times k)`$ matrix $`\beta`$

``` math
\beta=
\begin{bmatrix}
\Pi'_1 \\
\vdots  \\
\Pi'_p
\end{bmatrix}
```

We can then rewrite the model as a nonlinear regression (Karlsson, 2013)

``` math
y_t' =d_t'\Psi' + \left[w_t'-q_t'(I_p \otimes \Psi') \right]\beta +u_t'
```

where $`w_t'=(y_{t-1}',\dots,y_{t-p}')`$ is a $`kp`$-dimensional vector
of lagged endogenous variables and $`q_t'=(d_{t-1}',\dots,d_{t-p}')`$ is
a $`qp`$-dimensional vector of lagged deterministic (exogenous)
variables, $`I_p`$ is the $`(p \times p)`$ identity matrix and
$`\otimes`$ denotes the Kronecker product. This is how the likelihood is
written in the Stan code. The goal is to estimate $`\beta, \Psi`$, and
$`\Sigma_u`$, so priors are needed. First, prior independence between
$`\beta, \Psi`$ and $`\Sigma_u`$ is assumed. Starting with $`\beta`$, we
use the Minnesota prior

``` math
\mathrm{vec}(\beta) \sim \mathrm{N}_{kpk} \left[\theta_\beta,\Omega_\beta\right]
```

The prior means (the elements of $`\theta_\beta`$) are set to

``` math
\begin{aligned}
\mathrm{E}\left(\Pi_{\ell}^{(i,j)}\right)&=
\begin{cases}
\kappa & \text{if } \ell = 1 \ \mathrm{and} \ i = j \\
0 & \mathrm{otherwise}
\end{cases}\\
\kappa&=
\begin{cases}
\kappa^{level} & \text{if } \mathrm{variable} \ i \ \mathrm{is in level} \\
\kappa^{\Delta} & \text{if }\mathrm{variable} \ i \ \mathrm{is in difference}
\end{cases}\\
\end{aligned}
```

Here, the autoregressive coefficient $`\Pi_{\ell}^{(i,j)}`$ is element
$`\left(i,j\right)`$ of $`\Pi_{\ell}`$ for $`\ell=1,\dots,p`$. As such,
the Minnesota prior sets all prior means for the elements in $`\beta`$
to $`0`$, except for the elements that relate to the first own lags of
the variables, which are set to $`\kappa`$. If variable $`i`$ is in
level (e.g. nominal interest rate), then $`\kappa=\kappa^{level}`$, and
typical choices for $`\kappa^{level}`$ are $`1`$ or $`0.9`$. Evaluating
the equations at their prior means, equation $`i`$ becomes a random walk
if $`\kappa^{level}=1`$ and a persistent stationary AR(1) process if
$`\kappa^{level}=0.9`$. Since the steady state only exists if the
process is stationary, $`0.9`$ is recommended for the steady-state BVAR.
If variable $`i`$ is in difference (e.g. output growth), then
$`\kappa=\kappa^{\Delta}`$, and the most common choice for
$`\kappa^{\Delta}`$ is $`0`$, i.e. equation $`i`$ becomes (when
evaluating it at its prior means) a random walk expressed in first
differences. If a differenced variable still shows some degree of
persistence (can be examined with an ACF plot), a suitable value for
$`\kappa^{\Delta}`$ can be (for example) $`0.6`$ instead of $`0`$.
Moving on to the prior variances, $`\Omega_\beta`$ is a diagonal matrix
containing the prior variances for the elements in $`\beta`$. They are
specified as

``` math
\mathrm{Var}\left(\Pi_{\ell}^{(i,j)}\right)=
\begin{cases}
\left(\frac{\lambda_1}{\ell^{\lambda_3}}\right)^2 & \text{if } i = j \\
\left(\frac{\lambda_1 \lambda_2\sigma_i}{\ell^{\lambda_3}\sigma_j}\right)^2& \text{if } i \neq j
\end{cases}
```

Here $`\lambda_1`$, $`\lambda_2`$, and $`\lambda_3`$ are scalar
hyperparameters known as the overall tightness, the cross-equation
tightness and the lag decay rate. Furthermore, $`\sigma_i^2`$ is the
$`(i,i)`$:th element of $`\Sigma_u`$, which we do not know and therefore
replace with an estimate. In this package, it is replaced by the least
squares residual variance from a univariate autoregression for variable
$`i`$ with $`p`$ lags (including the constant and dummy/trend variable
if applicable). Moving on to $`\Psi`$, the prior we use is

``` math
\mathrm{vec}(\Psi) \sim \mathrm{N}_{kq}\left[\theta_\Psi,\Omega_\Psi\right]
```

This is really the core of the steady-state BVAR model. In
$`\theta_\Psi`$, we specify our prior beliefs about the location of the
steady state, and in $`\Omega_\Psi`$, which we assume to be a diagonal
matrix, we specify our degree of certainty in those prior beliefs.
Finally, the prior for $`\Sigma_u`$ is the usual noninformative Jeffreys
prior

``` math
p(\Sigma_u) \propto\left|\Sigma_u \right|^{-(k+1)/2}
```

However, if the user wants, an inverse Wishart prior can be used instead

``` math
\Sigma_u \sim \mathrm{IW}(V_0,m_0)
```

where $`V_0`$ is the scale matrix and $`m_0\geq k+2`$ is the number of
degrees of freedom. An uninformative prior can be (and is in this
package) specified by setting $`V_0=(m_0-k-1)\hat{\Sigma}_u`$ where
$`\hat{\Sigma}_u`$ is the least squares estimate from the VAR($`p`$)
(including the constant and dummy/trend variable if applicable), and
$`m_0=k+2`$.

This package also allows for stochastic volatility, where we let the
covariance matrix of the innovations vary over time such that we have a
time-varying covariance matrix $`\Sigma_{u,t}`$ (See other vignettes).
