
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SteadyStateBVAR

<!-- badges: start -->
<!-- badges: end -->

With this package the user can estimate the Steady-State BVAR(p) model
by Mattias Villani.

## Installation

You can install the development version of SteadyStateBVAR with:

``` r
#remotes::install_github("markjwbecker/SteadyStateBVAR", force = TRUE, upgrade = "never")
```

## Introduction

The Steady State BVAR($p$) model is

$$
y_t = \Psi x_t + A_1(y_{t-1}-\Psi x_{t-1})+\dots+A_p(y_{t-p}-\Psi x_{t-p})+u_t
$$

where $y_t$ is a $k$-dimensional vector of time series at time $t$, and
$x_t$ is a $q$-dimensional vector of deterministic variables at time
$t$, and $u_t \sim N_k(0,\Sigma_u)$ with independence between time
periods. Also, $A_i$ for $i=1,\dots,p$ is $(k \times k)$, and $\Psi$ is
$(k \times q)$. Note here that

$$
E(y_t)=\mu_t=\Psi x_t
$$

is the **steady state**. We can stack the $A$’s in the $(kp \times k)$
matrix $\beta$

$$
\beta=
\begin{bmatrix}
A'_1 \\ 
\vdots  \\
A'_p
\end{bmatrix}
$$

We can then rewrite the model as a nonlinear regression

$$
y_t' =x_t'\Psi' + \left[w_t'-q_t'(I_p \otimes \Psi') \right]\beta +u_t'
$$

where $w_t'=(y_{t-1}',\dots,y_{t-p}')$ is a $kp$-dimensional vector of
lagged endogenous variables, and $q_t'=(x_{t-1}',\dots,x_{t-p}')$ is a
$qp$-dimensional vector of lagged exogenous (deterministic) variables,
$I_p$ is the $(p \times p)$ identity matrix and $\otimes$ denotes the
Kronecker product. This is how the likelihood is written in the Stan
code. The goal is to estimate $\beta, \Psi$ and $\Sigma_u$, therefore
priors are needed. Starting with $\beta$, we use the Minnesota prior
where

$$
\textrm{vec}(\beta) \sim N_{kpk}\left[\textrm{vec}(\beta_0),\Sigma_{\textrm{vec}(\beta)}\right]
$$

First for $\beta_0$, the Minnesota prior sets all elements to $0$,
except for the elements that relate to the first own lags of the
variables, which are often set to $0.9$ or $1$ for level variables or
also to $0$ for growth rate variables. Further,
$\Sigma_{\textrm{vec}(\beta)}$ is a diagonal matrix containing the prior
variances for the elements in $\beta$. The prior is constructed such
that for the autoregressive coefficient $A_{\ell}^{(i,j)}$ which is
element $\left(i,j\right)$ of $A_{\ell}$ for $\ell=1,\dots,p,$ the prior
variance is given by

$$
\textrm{Var}\left(A_{l}^{(i,j)}\right)=
\begin{cases}
\left(\frac{\lambda_1}{\ell}\right)^2 & \text{if } i = j \\
\left(\frac{\lambda_1 \lambda_2}{\ell}\right)^2 \frac{\sigma_i^2}{\sigma_j^2}& \text{if } i \neq j
\end{cases}
$$

Here $\lambda_1$ and $\lambda_2$ are scalar hyperparameters where the
former is known as the overall tightness, and the latter as the
cross-equation tightness. Furthermore, $\sigma_i^2$ is the $(i,i)$:th
element of $\Sigma_u$, which we do not know, and therefore replace with
an estimate. In this package, it is replaced by the least squares
residual variance from a univariate autoregression for variable $i$ with
$p$ lags (including any potential deterministic variables). Moving on to
$\Psi$ the prior we use is

$$
\textrm{vec}(\Psi) \sim N_{kq}\left[\textrm{vec}(\Psi_0),\Sigma_{\textrm{vec}(\Psi)}\right]
$$

It is assumed that $\Sigma_{\textrm{vec}(\Psi)}$ is a diagonal matrix.
At last, the prior for $\Sigma_u$ is

$$
\Sigma_u \sim IW(V_0,m_0)
$$

Here $V_0$ is the scale matrix and $m_0\geq k+2$ are the degrees of
freedom. Since Stan does not allow the usual noninformative prior
$\left|\Sigma \right|^{-(k+1)/2}$, we will instead specify an
uninformative prior by setting $V_0=(m_0-k-1)\hat{\Sigma}_u$ where
$\hat{\Sigma}_u$ is the least squares estimate from the VAR($p$)
(including any potential deterministic regressors), and $m_0=k+2$.

## Example

We will now replicate the model in the empirical analysis in Section 4.1
in Villani (2009). First let us load the library and also load the data

``` r
remotes::install_github("markjwbecker/SteadyStateBVAR", force = TRUE, upgrade = "never")
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>       ✔  checking for file 'C:\Users\markj\AppData\Local\Temp\Rtmp4WiPWN\remotes3ad0a801fad\markjwbecker-SteadyStateBVAR-fc646e4/DESCRIPTION'
#>       ─  preparing 'SteadyStateBVAR':
#>    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>      NB: this package now depends on R (>=        NB: this package now depends on R (>= 3.5.0)
#>        WARNING: Added dependency on R >= 3.5.0 because serialized objects in
#>      serialize/load version 3 cannot be read in older versions of R.
#>      File(s) containing such objects:
#>        'SteadyStateBVAR/data/villani2009.rda'
#>        'SteadyStateBVAR/inst/STEADYSTATEBVAR.rds'
#> ─  building 'SteadyStateBVAR_0.1.0.tar.gz'
#>      
#> 
library(SteadyStateBVAR)
data("villani2009")
yt <- xt #villani uses xt to denote the endogenous variables
```

The data set contains quarterly data for Sweden over the time period
1980Q1–2005Q4. The seven variables are: trade-weighted measures of
foreign GDP growth $(\Delta y_f)$, CPI inflation $(\pi_f)$ and the
3-month interest rate $(i_f)$, the corresponding domestic variables
($\Delta y$, $\pi$ and $i$), and the level of the real exchange rate
defined as $q=s+p_f-p$, where $p_f$ and $p$ are the foreign and domestic
CPI levels (in logs) and $s$ is the (log of the) trade weighted nominal
exchange rate. As such we have

$$
y_t=
\begin{pmatrix}
\Delta y_f \\
\pi_f \\
i_f \\
\Delta y \\
\pi \\
i \\
q
\end{pmatrix}
$$

The growth rate variables are specified in terms of quarterly rates of
change, i.e. for a variable $z$, we have
$100 \ln \left(\frac{z_t}{z_{t-1}}\right)$. To simplify the prior
selection for the steady states later, we can multiply the growth rate
series by $4$ to instead get annualized quarterly rates of change. We do
that and then we can plot the data.

``` r
growth_idx <- c(1, 2, 4, 5)
yt[, growth_idx] <- yt[, growth_idx] * 4
plot.ts(yt)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

To model the Swedish financial crisis at the beginning of the 90s and
the subsequent shift in monetary policy to inflation targeting and
flexible exchange rate, $x_t$ (deterministic variables at time $t$)
includes a constant term and a dummy for the pre-crisis period, i.e.

$$
x_{t}' =
\begin{cases}
\begin{pmatrix}1 & 1\end{pmatrix} & \text{if } t \le 1992Q4 \\
\begin{pmatrix}1 & 0\end{pmatrix} & \text{if } t > 1992Q4
\end{cases}
$$

``` r
bp = 52 #breakpoint at 1992Q4
dum_var <- c(rep(1,bp), rep(0,nrow(xt)-bp)) #1 if t<=1992Q4, 0 if t>1992Q4
```

To formulate a prior on $\Psi$, note that the specification of $x_t$
implies the following parametrization of the steady state:

$$
\mu_t =
\begin{cases}
\psi_1 + \psi_2 & \text{if } t \le 1992Q4 \\
\psi_1 & \text{if } t > 1992Q4
\end{cases}
$$

where $\psi_i$ is the $i$:th column of $\Psi$. More on this soon. First
we assume $p=4$ is a good lag length for our model. Then we do some
setup.

``` r
p=4
stan_data <- BVAR_setup(yt, p, deterministic="constant_and_dummy", dummy=dum_var)
```

Now let us specify the priors. First we choose overall tightness
$\lambda_1=0.2$ and cross equation tightness $\lambda_2=0.5$. We then
specify the prior means for the first own lags of the variables. For
variables in growth rates, we set the prior mean to $0$, for variables
in levels, we set the prior mean to 0.9.

``` r
lambda1=0.2 #overall tightness
lambda2=0.5 #cross equation tightness

#first own lag prior mean
fol_pm=c(0,   #delta y_f
         0,   #pi_f
         0.9, #i_f
         0,   #delta y
         0,   #pi
         0.9, #i
         0.9  #q
)
```

Now for the steady state priors, we set them according to the 95% prior
probability intervals in Villani (2009). The ‘ppi()’ function is useful
here, simply input the 95% prior probability interval and then the
function outputs the prior mean and variance.

``` r
Psi_0 <- matrix(c(
                  ppi(2,3)$mean,
                  ppi(1.5,2.5)$mean,
                  ppi(4.5,5.5)$mean, 
                  ppi(2,2.5)$mean,
                  ppi(1.7,2.3)$mean,
                  ppi(4,4.5)$mean,
                  ppi(3.85,4)$mean,
                  ppi(-1,1)$mean,
                  ppi(1.5,2.5)$mean,
                  ppi(1.5,2.5)$mean,
                  ppi(-1,1)$mean,
                  ppi(4.3,5.7)$mean,
                  ppi(3,5.5)$mean,
                  ppi(-0.5,0.5)$mean
                  ), nrow=stan_data$k, ncol=stan_data$q)

vec_Psi_vars <- c(ppi(2,3)$var,
                  ppi(1.5,2.5)$var,
                  ppi(4.5,5.5)$var, 
                  ppi(2,2.5)$var,
                  ppi(1.7,2.3)$var,
                  ppi(4,4.5)$var,
                  ppi(3.85,4)$var,
                  ppi(-1,1)$var,
                  ppi(1.5,2.5)$var,
                  ppi(1.5,2.5)$var,
                  ppi(-1,1)$var,
                  ppi(4.3,5.7)$var,
                  ppi(3,5.5)$var,
                  ppi(-0.5,0.5)$var)
```

Let us check out the prior means in

``` r
Psi_0
#>       [,1] [,2]
#> [1,] 2.500 0.00
#> [2,] 2.000 2.00
#> [3,] 5.000 2.00
#> [4,] 2.250 0.00
#> [5,] 2.000 5.00
#> [6,] 4.250 4.25
#> [7,] 3.925 0.00
```

So if we are “after” 1992Q4, the steady state prior for variable $i$ is
row $i$ column $1$. If we are before or at 1992Q4 the steady state prior
for variable $i$ is row $i$ column $1$ + column $2$. Now lets input our
priors and then append them to ‘stan_data’. We also need to attach the
dummy variable.

``` r
priors <- priors(yt, p, lambda1, lambda2, fol_pm, Psi_0, vec_Psi_vars, dummy=dum_var)
stan_data <- c(stan_data, priors)
```

Like in Villani (2009), to incorporate that Sweden is a small economy
and therefore not likely to affect the foreign economy, we restrict the
upper right submatrix in each $A_i, i =1,\dots,k$, to the zero matrix.

``` r
n1 <- 3 #first 3 variables are foreign in yt
n2 <- 4 #the other 4 are domestic
n <- n1 + n2
p = 4

tmp <- matrix(1, n*p, n)

for(i in 1:p){
  rows <- ((i-1)*n + n1 + 1) : (i*n)
  cols <- 1:n1
  tmp[rows, cols] <- 0
  zero_indices <- which(c(tmp) == 0)
}
diag(stan_data$Sigma_vec_beta[zero_indices, zero_indices]) <- 0.00001

#example for A_1, we restrict the elements with zero entries to have very very small prior variance
#and since prior mean is zero, the posterior will be essentially zero 
t(tmp)[,1:7]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,]    1    1    1    0    0    0    0
#> [2,]    1    1    1    0    0    0    0
#> [3,]    1    1    1    0    0    0    0
#> [4,]    1    1    1    1    1    1    1
#> [5,]    1    1    1    1    1    1    1
#> [6,]    1    1    1    1    1    1    1
#> [7,]    1    1    1    1    1    1    1
```

Now we supply our forecast horizon $H$, and also the deterministic
variables for the future periods and then we fit the model.

``` r
#H <- 8
#X_pred <- cbind(rep(1, H), 0)
#fit <- estimate(stan_data, n_chains=4, iter=2000, warmup=500, H=H, X_pred=X_pred)
```

Lets plot the forecast. Lets select a 95% prediction interval and the
mean of the posterior as the point forecast. For the variables in
annualized quarter on quarter growth rates, we transform the historical
data and predictions to yearly growth rates.

``` r
#plot_forecast(fit, yt, ci=0.95, fcst_type="mean",growth_rate_idx=c(1,2,4,5))
```
