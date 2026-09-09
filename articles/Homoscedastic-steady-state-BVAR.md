# Homoscedastic steady-state BVAR (Villani, 2009)

Here we estimate the original (homoscedastic - i.e. constant innovation
covariance matrix \\\Sigma\_{u}\\) steady-state BVAR model from Section
4.1 of Villani (2009). See
[`?bvar`](https://markjwbecker.github.io/SteadyStateBVAR/reference/bvar.md)
for details on the model.

First, let us attach the package and load the data.

``` r

library(SteadyStateBVAR)
data("Villani2009")
yt <- Villani2009
```

The data set contains quarterly data for Sweden over the time period
1980Q1–2005Q4. The seven variables are: trade-weighted measures of
foreign GDP growth \\(\Delta y_f)\\, CPI inflation \\(\pi_f)\\ and the
3-month interest rate \\(i_f)\\, the corresponding domestic variables
(\\\Delta y\\, \\\pi\\ and \\i\\), and the level of the real exchange
rate defined as \\q=s+p_f-p\\, where \\p_f\\ and \\p\\ are the foreign
and domestic CPI levels (in logs) and \\s\\ is the (log of the)
trade-weighted nominal exchange rate. As such, we have

\\ y_t= \begin{pmatrix} \Delta y_f \\ \pi_f \\ i_f \\ \Delta y \\ \pi \\
i \\ q \end{pmatrix} \\

Also, we will leave out the last two observations, so the user can
compare the forecasts produced here to the last forecasts seen in
Figures 1-3 in Villani (2009) to verify that this implementation works
correctly.

``` r

yt <- ts(yt[1:102, ], start = start(yt), frequency = frequency(yt))
plot.ts(yt)
```

![plot of chunk HOMO-1](figure/HOMO-1-1.png)

plot of chunk HOMO-1

Also, let us create the bvar object which we will use throughout here.

``` r

bvar_obj <- bvar(data = yt)
```

To model the Swedish financial crisis at the beginning of the 90s and
the subsequent shift in monetary policy to inflation targeting and
flexible exchange rate, \\d_t\\ (deterministic variables at time \\t\\)
includes a constant term and a dummy for the pre-crisis period, i.e.

\\ d\_{t}' = \begin{cases} \begin{pmatrix}1 & 1\end{pmatrix} & \text{if
} t \le 1992Q4 \\ \begin{pmatrix}1 & 0\end{pmatrix} & \text{if } t \>
1992Q4 \end{cases} \\

``` r

breakpoint <- which(time(yt) == 1992.75)
dum_var <- c(rep(1,breakpoint), rep(0,nrow(yt)-breakpoint))
```

To formulate a prior on \\\Psi\\, note that the specification of \\d_t\\
implies the following parametrization of the steady state:

\\ \mu_t = \begin{cases} \psi_1 + \psi_2 & \text{if } t \le 1992Q4 \\
\psi_1 & \text{if } t \> 1992Q4 \end{cases} \\

where \\\psi_i\\ denotes the \\i\\:th column of \\\Psi\\. We are now
ready to set up the model. Although it is not mentioned which lag length
is used in Villani (2009), we assume \\p=4\\.

``` r

bvar_obj <- setup(bvar_obj,
                  p=4,
                  deterministic = "constant_and_dummy",
                  dummy = dum_var)
```

Now let us specify the priors. We first consider \\\beta\\ (Minnesota
prior). We choose the same values for the hyperparameters as in Villani
(2009), i.e. an overall tightness of \\\lambda_1=0.2\\, a cross-equation
tightness of \\\lambda_2=0.5\\, and a lag decay rate of \\\lambda_3=1\\.
We then specify the prior means for the first own lags of the variables.
We follow Villani (2009), and as such for variables in growth rates, we
set the prior mean to \\0\\. For variables in levels, we set the prior
mean to \\0.9\\.

``` r

lambda_1 <- 0.2
lambda_2 <- 0.5
lambda_3 <- 1.0

#fol_pm = first own lag prior means
fol_pm=c(0,   #delta y_f
         0,   #pi_f
         0.9, #i_f
         0,   #delta y
         0,   #pi
         0.9, #i
         0.9  #q
         )
```

Now moving on to \\\Psi\\, i.r. the steady-state priors, we set them
according to the 95% prior probability intervals (normal distribution)
in Table I in Villani (2009). We first note that for our data here, the
growth rate variables (\\\Delta y_f, \pi_f, \Delta y, \pi\\) are
specified in terms of quarterly rates of change/quarter-on-quarter
growth, i.e. for a variable \\x\\ which is on a quarterly frequency
(`freq=4`), the quarterly growth rate is \\100(\ln (x_t) - \ln
x\_{t-1})\\. The 95% prior probability intervals in Table I are
specified in terms of annualized quarterly growth rates \\400(\ln
(x_t) - \ln x\_{t-1})\\.

The
[`ppi()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/ppi.md)
function is useful here. Simply input the desired 95% prior probability
interval (normal distribution) on the annualized scale with
`annualized_growthrate=TRUE` with the correct frequency `freq`, and the
function returns the corresponding prior mean and variance on the
original scale (quarter-on-quarter growth). Of course, we could also
just annualize our data beforehand, and set
`annualized_growthrate=FALSE`.

``` r

#psi_1 = Psi col 1
#psi_2 = Psi col 2

theta_Psi <- 
  c(
  ppi( 2.00,  3.00, interval = 0.95, annualized_growthrate=TRUE, freq=4)$mean,   #psi_1: delta y_f
  ppi( 1.50,  2.50, interval = 0.95, annualized_growthrate=TRUE, freq=4)$mean,   #psi_1: pi_f
  ppi( 4.50,  5.50, interval = 0.95                                    )$mean,   #psi_1: i_f
  ppi( 2.00,  2.50, interval = 0.95, annualized_growthrate=TRUE, freq=4)$mean,   #psi_1: delta y
  ppi( 1.70,  2.30, interval = 0.95, annualized_growthrate=TRUE, freq=4)$mean,   #psi_1: pi
  ppi( 4.00,  4.50, interval = 0.95                                    )$mean,   #psi_1: i
  ppi( 3.85,  4.00, interval = 0.95                                    )$mean,   #psi_1: q
  ppi(-1.00,  1.00, interval = 0.95, annualized_growthrate=TRUE, freq=4)$mean,   #psi_2: delta y_f
  ppi( 1.50,  2.50, interval = 0.95, annualized_growthrate=TRUE, freq=4)$mean,   #psi_2: pi_f
  ppi( 1.50,  2.50, interval = 0.95                                    )$mean,   #psi_2: i_f
  ppi(-1.00,  1.00, interval = 0.95, annualized_growthrate=TRUE, freq=4)$mean,   #psi_2: delta y
  ppi( 4.30,  5.70, interval = 0.95, annualized_growthrate=TRUE, freq=4)$mean,   #psi_2: pi
  ppi( 3.00,  5.50, interval = 0.95                                    )$mean,   #psi_2: i
  ppi(-0.50,  0.50, interval = 0.95                                    )$mean    #psi_2: q
  )

Omega_Psi <- 
  diag(
  c(
  ppi( 2.00,  3.00, interval = 0.95, annualized_growthrate=TRUE, freq=4)$var,    #psi_1: delta y_f
  ppi( 1.50,  2.50, interval = 0.95, annualized_growthrate=TRUE, freq=4)$var,    #psi_1: pi_f
  ppi( 4.50,  5.50, interval = 0.95                                    )$var,    #psi_1: i_f
  ppi( 2.00,  2.50, interval = 0.95, annualized_growthrate=TRUE, freq=4)$var,    #psi_1: delta y
  ppi( 1.70,  2.30, interval = 0.95, annualized_growthrate=TRUE, freq=4)$var,    #psi_1: pi
  ppi( 4.00,  4.50, interval = 0.95                                    )$var,    #psi_1: i
  ppi( 3.85,  4.00, interval = 0.95                                    )$var,    #psi_1: q
  ppi(-1.00,  1.00, interval = 0.95, annualized_growthrate=TRUE, freq=4)$var,    #psi_2: delta y_f
  ppi( 1.50,  2.50, interval = 0.95, annualized_growthrate=TRUE, freq=4)$var,    #psi_2: pi_f
  ppi( 1.50,  2.50, interval = 0.95                                    )$var,    #psi_2: i_f
  ppi(-1.00,  1.00, interval = 0.95, annualized_growthrate=TRUE, freq=4)$var,    #psi_2: delta y
  ppi( 4.30,  5.70, interval = 0.95, annualized_growthrate=TRUE, freq=4)$var,    #psi_2: pi
  ppi( 3.00,  5.50, interval = 0.95                                    )$var,    #psi_2: i
  ppi(-0.50,  0.50, interval = 0.95                                    )$var     #psi_2: q
  )
  )
```

Finally for \\\Sigma_u\\ we will use the noninformative Jeffreys prior
\\\left\|\Sigma_u \right\|^{-(k+1)/2}\\, as done in Villani (2009). Now
we simply pass everything to the
[`priors()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md)
function.

``` r

bvar_obj <- priors(bvar_obj,
                   lambda_1,
                   lambda_2,
                   lambda_3,
                   fol_pm,
                   theta_Psi,
                   Omega_Psi,
                   Jeffreys=TRUE)
```

We can plot our steady-state priors now after we have passed them to
[`priors()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md)

``` r

par(mfrow=c(3,3))
steady_state_priors_plot(bvar_obj, interval = 0.95, growth_rate_idx = c(1,2,4,5))
par(mfrow=c(1,1))
```

Continuing, as in Villani (2009), we incorporate the assumption that
Sweden is a small economy and therefore unlikely to affect the foreign
economy by restricting the upper-right submatrix of \\\Pi\_\ell\\ for
\\\ell =1,\dots,p\\ or equivalently restricting the bottom-left
submatrix of \\\Pi\_\ell'\\ to the zero matrix. This technique is called
“block exogeneity” (Dieppe, Legrand, and van Roye, 2016). In essence we
treat the foreign economy as exogenous to the domestic economy, although
it is not exogenous in the strict sense (Karlsson, 2013).

``` r

p <- bvar_obj$setup$p
k <- bvar_obj$setup$k
kf <- 3 #first 3 variables are foreign in yt

restriction_matrix <- matrix(1, k*p, k)

for(i in 1:p){
  rows <- ((i-1)*k + kf + 1) : (i*k)
  cols <- 1:kf
  restriction_matrix[rows, cols] <- 0
}
```

We simply pass our \\(kp \times k)\\ restriction matrix to the
[`restrict_beta()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/restrict_beta.md)
function:

``` r

bvar_obj <- restrict_beta(bvar_obj, restriction_matrix)
#> Restrictions applied using restriction matrix:
#> 
#>              delta y_f pi_f i_f delta y pi i q
#> delta y_f.l1         1    1   1       1  1 1 1
#> pi_f.l1              1    1   1       1  1 1 1
#> i_f.l1               1    1   1       1  1 1 1
#> delta y.l1           0    0   0       1  1 1 1
#> pi.l1                0    0   0       1  1 1 1
#> i.l1                 0    0   0       1  1 1 1
#> q.l1                 0    0   0       1  1 1 1
#> delta y_f.l2         1    1   1       1  1 1 1
#> pi_f.l2              1    1   1       1  1 1 1
#> i_f.l2               1    1   1       1  1 1 1
#> delta y.l2           0    0   0       1  1 1 1
#> pi.l2                0    0   0       1  1 1 1
#> i.l2                 0    0   0       1  1 1 1
#> q.l2                 0    0   0       1  1 1 1
#> delta y_f.l3         1    1   1       1  1 1 1
#> pi_f.l3              1    1   1       1  1 1 1
#> i_f.l3               1    1   1       1  1 1 1
#> delta y.l3           0    0   0       1  1 1 1
#> pi.l3                0    0   0       1  1 1 1
#> i.l3                 0    0   0       1  1 1 1
#> q.l3                 0    0   0       1  1 1 1
#> delta y_f.l4         1    1   1       1  1 1 1
#> pi_f.l4              1    1   1       1  1 1 1
#> i_f.l4               1    1   1       1  1 1 1
#> delta y.l4           0    0   0       1  1 1 1
#> pi.l4                0    0   0       1  1 1 1
#> i.l4                 0    0   0       1  1 1 1
#> q.l4                 0    0   0       1  1 1 1
#> 
#> 1 indicates that the parameter is free
#> 0 indicates that the parameter is restricted to zero
```

The function tells us which autoregressive parameters in \\\beta\\ we
restrict to zero.

Now we are almost ready to fit (estimate) the model. When estimating the
model, we are at the same time generating draws from the joint
predictive distribution. To accomplish the latter, we need the forecast
horizon \\H\\, and also a matrix containing the deterministic variables
(\\d_t\\) for the future periods

\\ d\_{\text{pred}}=\begin{bmatrix}d\_{T+1}' \\ \vdots\\ d\_{T+H}'
\end{bmatrix} \\

Since the deterministic variables are i) a constant and ii) a dummy
indicating whether \\t \leq 1992Q4\\, we simply set

\\ d\_{T+1}'=\ldots=d\_{T+H}'=\begin{pmatrix} 1 & 0 \end{pmatrix} \\

``` r

d_pred <- cbind(rep(1, 12), 0)
```

However,
[`fit()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md)
automatically creates `d_pred`, so we do not need to bother with it.

For the estimation, let us choose 4 markov chains, with each having 4000
iterations, and where 2000 of those 4000 are warmup/burn-in iterations.

``` r

bvar_obj <- fit(bvar_obj,
                H = 12,
                iter = 6000,
                warmup = 2000,
                chains = 4,
                cores = 4)
#> NOTE: d_pred not supplied
#> it is assumed that the dummy stays at its last observed value (0) for all 12 forecast periods.
#> ------------------------------------------------------------
#> Forecast horizon:
#> 12
#> 
#> Future deterministic variables (d_pred):
#>      constant dummy
#> h=1         1     0
#> h=2         1     0
#> h=3         1     0
#> h=4         1     0
#> h=5         1     0
#> h=6         1     0
#> h=7         1     0
#> h=8         1     0
#> h=9         1     0
#> h=10        1     0
#> h=11        1     0
#> h=12        1     0
#> ------------------------------------------------------------
#> Estimating Stan model:
#> steady_state_bvar_homoscedastic_jeffreys_prior
#> 
#> Also generating draws from the joint predictive distribution
#> 
#> ...
#> SAMPLING FINISHED
```

Let us look at the elementwise posterior means of \\\beta\\, \\\Psi\\,
and \\\Sigma_u\\.

``` r

summary(bvar_obj)
#> Posterior mean estimates
#> ------------------------
#> 
#> 
#> beta
#> --------------------------------------------------------------------------------              
#>                delta y_f  pi_f   i_f delta y    pi     i     q
#>   delta y_f.l1      0.18  0.03 -0.01    0.12  0.07 -0.12  0.00
#>   pi_f.l1          -0.02  0.32  0.25    0.12 -0.07  0.01  0.00
#>   i_f.l1            0.00  0.04  0.92   -0.04  0.06  0.05  0.00
#>   delta y.l1        0.00  0.00  0.00    0.23 -0.09 -0.10  0.00
#>   pi.l1             0.00  0.00  0.00    0.00  0.08  0.06  0.00
#>   i.l1              0.00  0.00  0.00    0.00  0.02  0.76  0.00
#>   q.l1              0.00  0.00  0.00    1.22  3.94  0.76  0.93
#>   delta y_f.l2      0.03 -0.01  0.09    0.02 -0.02  0.10  0.00
#>   pi_f.l2           0.01  0.02  0.04    0.00 -0.03 -0.15  0.00
#>   i_f.l2           -0.02 -0.01 -0.01    0.00  0.04  0.07  0.00
#>   delta y.l2        0.00  0.00  0.00    0.11 -0.01  0.15  0.00
#>   pi.l2             0.00  0.00  0.00    0.01 -0.04 -0.05  0.00
#>   i.l2              0.00  0.00  0.00   -0.01  0.01  0.04  0.00
#>   q.l2              0.00  0.00  0.00    0.55 -0.38  0.29 -0.04
#>   delta y_f.l3      0.01 -0.01  0.00    0.02 -0.01  0.00  0.00
#>   pi_f.l3          -0.02  0.06 -0.01    0.00  0.08  0.02  0.00
#>   i_f.l3            0.00  0.00  0.02    0.00  0.00  0.03  0.00
#>   delta y.l3        0.00  0.00  0.00    0.07  0.01 -0.02  0.00
#>   pi.l3             0.00  0.00  0.00    0.00  0.02 -0.02  0.00
#>   i.l3              0.00  0.00  0.00    0.01  0.00  0.01  0.00
#>   q.l3              0.00  0.00  0.00   -0.14 -0.02 -0.58  0.00
#>   delta y_f.l4      0.03 -0.01  0.00   -0.01  0.02  0.02  0.00
#>   pi_f.l4           0.00  0.16 -0.03    0.00  0.01  0.02  0.00
#>   i_f.l4            0.00  0.00 -0.02    0.00  0.00  0.03  0.00
#>   delta y.l4        0.00  0.00  0.00   -0.08  0.01  0.03  0.00
#>   pi.l4             0.00  0.00  0.00    0.00  0.06 -0.01  0.00
#>   i.l4              0.00  0.00  0.00    0.00 -0.01  0.00  0.00
#>   q.l4              0.00  0.00  0.00   -0.15 -0.07 -0.17 -0.01
#> --------------------------------------------------------------------------------
#> 
#> 
#> Psi
#> --------------------------------------------------------------------------------           
#>             [,1]  [,2]
#>   delta y_f 0.58  0.08
#>   pi_f      0.50  0.46
#>   i_f       4.95  2.02
#>   delta y   0.58 -0.04
#>   pi        0.49  1.15
#>   i         4.29  4.45
#>   q         3.92 -0.10
#> --------------------------------------------------------------------------------
#> 
#> 
#> Sigma_u
#> --------------------------------------------------------------------------------           
#>             delta y_f  pi_f  i_f delta y    pi     i     q
#>   delta y_f      0.15 -0.01 0.01    0.07 -0.01  0.00  0.00
#>   pi_f          -0.01  0.09 0.05    0.01  0.13  0.04  0.00
#>   i_f            0.01  0.05 0.51    0.01  0.18  0.11  0.00
#>   delta y        0.07  0.01 0.01    0.19 -0.05 -0.01  0.00
#>   pi            -0.01  0.13 0.18   -0.05  0.59  0.11  0.00
#>   i              0.00  0.04 0.11   -0.01  0.11  1.56 -0.01
#>   q              0.00  0.00 0.00    0.00  0.00 -0.01  0.00
#> --------------------------------------------------------------------------------
```

We can access the elementwise posterior means or medians with
`bvar_obj$fit$posterior_means`/`bvar_obj$fit$posterior_medians` if
needed.

Note that `bvar_obj$fit$stan` is an object of class `stanfit`.

``` r

(stanfit <- bvar_obj$fit$stan)
#> Inference for Stan model: steady_state_bvar_homoscedastic_jeffreys_prior.
#> 4 chains, each with iter=6000; warmup=2000; thin=1; 
#> post-warmup draws per chain=4000, total post-warmup draws=16000.
#> 
#>                  mean se_mean    sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
#> beta[1,1]        0.18    0.00  0.09   0.00   0.12   0.18   0.24   0.36 27973    1
#> beta[1,2]        0.03    0.00  0.05  -0.07   0.00   0.03   0.06   0.12 30693    1
#> beta[1,3]       -0.01    0.00  0.13  -0.27  -0.10  -0.01   0.07   0.24 30646    1
#> beta[1,4]        0.12    0.00  0.08  -0.05   0.06   0.12   0.18   0.29 28176    1
#> beta[1,5]        0.07    0.00  0.14  -0.20  -0.02   0.07   0.17   0.34 27704    1
#> beta[1,6]       -0.12    0.00  0.24  -0.59  -0.29  -0.12   0.04   0.36 31525    1
#> beta[1,7]        0.00    0.00  0.01  -0.02  -0.01   0.00   0.00   0.01 34069    1
#> beta[2,1]       -0.02    0.00  0.09  -0.20  -0.08  -0.02   0.04   0.16 27959    1
#> beta[2,2]        0.32    0.00  0.08   0.16   0.26   0.32   0.37   0.48 22207    1
#> beta[2,3]        0.25    0.00  0.17  -0.08   0.14   0.25   0.37   0.59 28795    1
#> beta[2,4]        0.12    0.00  0.11  -0.09   0.05   0.12   0.19   0.33 26303    1
#> beta[2,5]       -0.07    0.00  0.19  -0.45  -0.20  -0.07   0.06   0.30 21450    1
#> beta[2,6]        0.01    0.00  0.32  -0.62  -0.21   0.01   0.23   0.64 27822    1
#> beta[2,7]        0.00    0.00  0.01  -0.01   0.00   0.00   0.01   0.02 33192    1
#> beta[3,1]        0.00    0.00  0.03  -0.06  -0.02  -0.01   0.02   0.05 23108    1
#> beta[3,2]        0.04    0.00  0.02   0.00   0.03   0.04   0.05   0.08 22721    1
#> beta[3,3]        0.92    0.00  0.07   0.78   0.87   0.92   0.97   1.07 22606    1
#> beta[3,4]       -0.04    0.00  0.04  -0.11  -0.06  -0.04  -0.01   0.03 20949    1
#> beta[3,5]        0.06    0.00  0.06  -0.07   0.02   0.06   0.10   0.18 20721    1
#> beta[3,6]        0.05    0.00  0.11  -0.16  -0.02   0.05   0.12   0.26 24579    1
#> beta[3,7]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.01 32982    1
#> beta[4,1]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16318    1
#> beta[4,2]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16103    1
#> beta[4,3]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15919    1
#> beta[4,4]        0.23    0.00  0.09   0.06   0.17   0.23   0.29   0.40 22004    1
#> beta[4,5]       -0.09    0.00  0.12  -0.32  -0.17  -0.09  -0.01   0.14 26890    1
#> beta[4,6]       -0.10    0.00  0.21  -0.50  -0.24  -0.10   0.04   0.31 28704    1
#> beta[4,7]        0.00    0.00  0.00  -0.01   0.00   0.00   0.00   0.01 31715    1
#> beta[5,1]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15810    1
#> beta[5,2]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16122    1
#> beta[5,3]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15479    1
#> beta[5,4]        0.00    0.00  0.04  -0.08  -0.02   0.00   0.03   0.08 29666    1
#> beta[5,5]        0.08    0.00  0.09  -0.09   0.02   0.08   0.13   0.25 26864    1
#> beta[5,6]        0.06    0.00  0.12  -0.18  -0.03   0.06   0.14   0.30 32197    1
#> beta[5,7]        0.00    0.00  0.00  -0.01   0.00   0.00   0.00   0.00 29923    1
#> beta[6,1]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15553    1
#> beta[6,2]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15681    1
#> beta[6,3]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16190    1
#> beta[6,4]        0.00    0.00  0.02  -0.04  -0.01   0.00   0.01   0.04 26817    1
#> beta[6,5]        0.02    0.00  0.04  -0.05   0.00   0.02   0.04   0.09 25431    1
#> beta[6,6]        0.76    0.00  0.08   0.60   0.70   0.76   0.82   0.93 23770    1
#> beta[6,7]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 18997    1
#> beta[7,1]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15724    1
#> beta[7,2]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16413    1
#> beta[7,3]        0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16854    1
#> beta[7,4]        1.22    0.01  0.87  -0.48   0.63   1.22   1.80   2.91 21874    1
#> beta[7,5]        3.94    0.01  1.45   1.08   2.99   3.95   4.94   6.72 22042    1
#> beta[7,6]        0.76    0.02  2.64  -4.38  -1.03   0.77   2.54   5.89 22210    1
#> beta[7,7]        0.93    0.00  0.08   0.78   0.88   0.93   0.98   1.09 21376    1
#> beta[8,1]        0.03    0.00  0.07  -0.11  -0.02   0.03   0.08   0.17 31317    1
#> beta[8,2]       -0.01    0.00  0.03  -0.07  -0.03  -0.01   0.01   0.05 37393    1
#> beta[8,3]        0.09    0.00  0.08  -0.07   0.04   0.09   0.15   0.25 38019    1
#> beta[8,4]        0.02    0.00  0.05  -0.07  -0.01   0.02   0.06   0.12 32215    1
#> beta[8,5]       -0.02    0.00  0.09  -0.18  -0.08  -0.02   0.04   0.15 35755    1
#> beta[8,6]        0.10    0.00  0.15  -0.20   0.00   0.10   0.20   0.39 37247    1
#> beta[8,7]        0.00    0.00  0.00  -0.01   0.00   0.00   0.00   0.01 34741    1
#> beta[9,1]        0.01    0.00  0.06  -0.12  -0.03   0.01   0.05   0.13 33932    1
#> beta[9,2]        0.02    0.00  0.07  -0.11  -0.02   0.02   0.07   0.15 30084    1
#> beta[9,3]        0.04    0.00  0.12  -0.18  -0.04   0.04   0.12   0.27 31658    1
#> beta[9,4]        0.00    0.00  0.07  -0.14  -0.05   0.00   0.04   0.13 33974    1
#> beta[9,5]       -0.03    0.00  0.12  -0.27  -0.11  -0.03   0.05   0.21 32340    1
#> beta[9,6]       -0.15    0.00  0.21  -0.55  -0.29  -0.15  -0.01   0.26 34418    1
#> beta[9,7]        0.00    0.00  0.00  -0.01   0.00   0.00   0.01   0.01 36677    1
#> beta[10,1]      -0.02    0.00  0.02  -0.06  -0.03  -0.02   0.00   0.03 26464    1
#> beta[10,2]      -0.01    0.00  0.02  -0.04  -0.02  -0.01   0.00   0.02 31888    1
#> beta[10,3]      -0.01    0.00  0.07  -0.15  -0.06  -0.01   0.04   0.14 23559    1
#> beta[10,4]       0.00    0.00  0.03  -0.05  -0.02   0.00   0.02   0.05 30850    1
#> beta[10,5]       0.04    0.00  0.05  -0.05   0.01   0.04   0.08   0.14 33102    1
#> beta[10,6]       0.07    0.00  0.08  -0.08   0.02   0.07   0.12   0.22 28534    1
#> beta[10,7]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 26043    1
#> beta[11,1]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16225    1
#> beta[11,2]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16576    1
#> beta[11,3]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15135    1
#> beta[11,4]       0.11    0.00  0.07  -0.02   0.07   0.11   0.16   0.25 31050    1
#> beta[11,5]      -0.01    0.00  0.08  -0.16  -0.06  -0.01   0.04   0.14 32528    1
#> beta[11,6]       0.15    0.00  0.14  -0.12   0.06   0.15   0.24   0.42 32833    1
#> beta[11,7]       0.00    0.00  0.00  -0.01   0.00   0.00   0.00   0.01 38595    1
#> beta[12,1]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15496    1
#> beta[12,2]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15752    1
#> beta[12,3]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15748    1
#> beta[12,4]       0.01    0.00  0.02  -0.04  -0.01   0.01   0.03   0.06 36929    1
#> beta[12,5]      -0.04    0.00  0.07  -0.17  -0.09  -0.04   0.00   0.09 29593    1
#> beta[12,6]      -0.05    0.00  0.08  -0.20  -0.10  -0.05   0.00   0.10 35660    1
#> beta[12,7]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 21394    1
#> beta[13,1]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15884    1
#> beta[13,2]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15909    1
#> beta[13,3]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16113    1
#> beta[13,4]      -0.01    0.00  0.01  -0.04  -0.02  -0.01   0.00   0.02 33594    1
#> beta[13,5]       0.01    0.00  0.03  -0.04   0.00   0.01   0.03   0.06 30084    1
#> beta[13,6]       0.04    0.00  0.07  -0.10  -0.01   0.04   0.09   0.18 25408    1
#> beta[13,7]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16967    1
#> beta[14,1]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 15742    1
#> beta[14,2]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16646    1
#> beta[14,3]       0.00    0.00  0.00   0.00   0.00   0.00   0.00   0.00 16781    1
#> beta[14,4]       0.55    0.00  0.66  -0.75   0.11   0.56   1.00   1.83 28470    1
#> beta[14,5]      -0.38    0.01  1.16  -2.66  -1.16  -0.38   0.40   1.91 26331    1
#> beta[14,6]       0.29    0.01  1.98  -3.61  -1.02   0.28   1.62   4.20 30078    1
#> beta[14,7]      -0.04    0.00  0.07  -0.18  -0.09  -0.04   0.01   0.11 25495    1
#> beta[15,1]       0.01    0.00  0.05  -0.10  -0.03   0.01   0.05   0.11 36319    1
#> beta[15,2]      -0.01    0.00  0.02  -0.06  -0.03  -0.01   0.00   0.03 33028    1
#>  [ reached 'max' / getOption("max.print") -- omitted 1259 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Wed Sep  9 03:40:49 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

As such, we can do the usual `rstan` inference on our fitted model. Let
us plot the posterior draws of the post-crisis steady-state
\\mu\_{t,i}\\ of inflation at t = 102, i.e. \\mu\_{t=102, i=5}\\ Note
that the index in the code is 98 because we have 102 observations where
the first \\p\\ are lost due to effective sample size \\T=N-p\\.

``` r

rstan::plot(stanfit, pars=c("mu[98,5]"), plotfun="hist")
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![plot of chunk HOMO-2](figure/HOMO-2-1.png) We need to multiply the
steady-state coefficient by 4 to obtain the annualized rate.

``` r

posterior <- rstan::extract(stanfit)
posterior$mu
#> , , 1
#> 
#>           
#> iterations      [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]     [,10]     [,11]
#>       [1,] 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221
#>       [2,] 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973
#>       [3,] 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140
#>       [4,] 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452
#>       [5,] 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882
#>       [6,] 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837
#>       [7,] 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086
#>       [8,] 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645
#>       [9,] 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265
#>      [10,] 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198
#>           
#> iterations     [,12]     [,13]     [,14]     [,15]     [,16]     [,17]     [,18]     [,19]     [,20]     [,21]     [,22]
#>       [1,] 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221
#>       [2,] 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973
#>       [3,] 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140
#>       [4,] 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452
#>       [5,] 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882
#>       [6,] 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837
#>       [7,] 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086
#>       [8,] 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645
#>       [9,] 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265
#>      [10,] 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198
#>           
#> iterations     [,23]     [,24]     [,25]     [,26]     [,27]     [,28]     [,29]     [,30]     [,31]     [,32]     [,33]
#>       [1,] 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221
#>       [2,] 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973
#>       [3,] 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140
#>       [4,] 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452
#>       [5,] 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882
#>       [6,] 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837
#>       [7,] 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086
#>       [8,] 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645
#>       [9,] 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265
#>      [10,] 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198
#>           
#> iterations     [,34]     [,35]     [,36]     [,37]     [,38]     [,39]     [,40]     [,41]     [,42]     [,43]     [,44]
#>       [1,] 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221 0.5595221
#>       [2,] 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973 0.6155973
#>       [3,] 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140 0.6024140
#>       [4,] 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452 0.8602452
#>       [5,] 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882 0.5771882
#>       [6,] 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837 0.5235837
#>       [7,] 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086 0.6052086
#>       [8,] 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645 0.7760645
#>       [9,] 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265 0.6830265
#>      [10,] 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198 0.6084198
#>           
#> iterations     [,45]     [,46]     [,47]     [,48]     [,49]     [,50]     [,51]     [,52]     [,53]     [,54]     [,55]
#>       [1,] 0.5595221 0.5595221 0.5595221 0.5595221 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992
#>       [2,] 0.6155973 0.6155973 0.6155973 0.6155973 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229
#>       [3,] 0.6024140 0.6024140 0.6024140 0.6024140 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572
#>       [4,] 0.8602452 0.8602452 0.8602452 0.8602452 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777
#>       [5,] 0.5771882 0.5771882 0.5771882 0.5771882 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207
#>       [6,] 0.5235837 0.5235837 0.5235837 0.5235837 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700
#>       [7,] 0.6052086 0.6052086 0.6052086 0.6052086 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030
#>       [8,] 0.7760645 0.7760645 0.7760645 0.7760645 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792
#>       [9,] 0.6830265 0.6830265 0.6830265 0.6830265 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774
#>      [10,] 0.6084198 0.6084198 0.6084198 0.6084198 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147
#>           
#> iterations     [,56]     [,57]     [,58]     [,59]     [,60]     [,61]     [,62]     [,63]     [,64]     [,65]     [,66]
#>       [1,] 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992
#>       [2,] 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229
#>       [3,] 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572
#>       [4,] 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777
#>       [5,] 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207
#>       [6,] 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700
#>       [7,] 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030
#>       [8,] 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792
#>       [9,] 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774
#>      [10,] 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147
#>           
#> iterations     [,67]     [,68]     [,69]     [,70]     [,71]     [,72]     [,73]     [,74]     [,75]     [,76]     [,77]
#>       [1,] 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992
#>       [2,] 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229
#>       [3,] 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572
#>       [4,] 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777
#>       [5,] 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207
#>       [6,] 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700
#>       [7,] 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030
#>       [8,] 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792
#>       [9,] 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774
#>      [10,] 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147
#>           
#> iterations     [,78]     [,79]     [,80]     [,81]     [,82]     [,83]     [,84]     [,85]     [,86]     [,87]     [,88]
#>       [1,] 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992
#>       [2,] 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229
#>       [3,] 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572
#>       [4,] 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777
#>       [5,] 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207
#>       [6,] 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700
#>       [7,] 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030
#>       [8,] 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792
#>       [9,] 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774
#>      [10,] 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147
#>           
#> iterations     [,89]     [,90]     [,91]     [,92]     [,93]     [,94]     [,95]     [,96]     [,97]     [,98]
#>       [1,] 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992 0.6076992
#>       [2,] 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229 0.6226229
#>       [3,] 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572 0.5522572
#>       [4,] 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777 0.6477777
#>       [5,] 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207 0.5666207
#>       [6,] 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700 0.6176700
#>       [7,] 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030 0.5352030
#>       [8,] 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792 0.5903792
#>       [9,] 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774 0.6215774
#>      [10,] 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147 0.5746147
#> 
#>  [ reached 'max' / getOption("max.print") -- omitted 6 slices ]
hist(4*posterior$mu[,98,5], col="darkred", breaks=30)
abline(v=mean(4*posterior$mu[,98,5]), col="lightgreen")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

plot of chunk unnamed-chunk-14

We can also look at the model forecasts directly with `rstan`. Remember
that we left out the last two observations/quarters, so let us look at
our forecasts of the domestic interest rate and compare them with the
actual values.

``` r

(Villani2009[103:104,6]) #true values
#> [1] 1.478503 1.563795

rstan::plot(stanfit,
            pars=c("y_pred[1,6]", "y_pred[2,6]"),
            show_density = TRUE,
            ci_level = 0.95,
            fill_color = "blue")
#> ci_level: 0.95 (95% intervals)
#> outer_level: 0.95 (95% intervals)
```

![plot of chunk HOMO-3](figure/HOMO-3-1.png)

plot of chunk HOMO-3

So the model overshot a bit, but the true values are within the 68%
prediction interval. Now let us plot the forecasts along with the
historical data. We will choose a 68% prediction interval (“pi”) and the
mean of the predictive distribution as the point forecast. For variables
in quarter-on-quarter growth rates, we transform the historical data and
predictions to yearly growth rates with ‘growth_rate_idx’ where we
specify the index of the growth rate variables in \\y_t\\. Note that
this is not annualization, but we are now computing \\100(\ln x_t - \ln
x\_{t-4})\\, i.e. the annual growth rate, by summing up to fourth
differences.

``` r

fcst <- forecast(bvar_obj,
                 pi = 0.95,
                 fcst_type = "mean",
                 growth_rate_idx = c(4,5),
                 plot_idx = c(4,5,6),
                 ss = TRUE,
                 ss_type = "mean",
                 ss_ci = 0.95,
                 show_all = FALSE)
```

![plot of chunk HOMO-4](figure/HOMO-4-1.png)![plot of chunk
HOMO-4](figure/HOMO-4-2.png)![plot of chunk HOMO-4](figure/HOMO-4-3.png)

For further inspection, we can print the point forecasts

``` r

print(fcst$forecast)
#>       delta y_f      pi_f      i_f  delta y        pi        i        q
#>  [1,] 0.6210271 0.5212379 2.828049 2.778007 0.9236902 2.021417 3.994970
#>  [2,] 0.6433043 0.4480442 3.022851 3.053265 0.9654955 2.056584 3.991218
#>  [3,] 0.6338076 0.4227043 3.164017 3.364018 1.7593205 2.146578 3.986508
#>  [4,] 0.6339701 0.4937232 3.277124 3.443972 1.7494867 2.266961 3.980849
#>  [5,] 0.6348606 0.4560315 3.408140 3.393316 1.7917306 2.385474 3.975736
#>  [6,] 0.6352717 0.4379055 3.522596 3.316294 1.8607163 2.463094 3.971383
#>  [7,] 0.6329443 0.4270316 3.616697 3.207013 1.8858405 2.561096 3.967736
#>  [8,] 0.6319848 0.4365934 3.708594 3.119617 1.8845718 2.679536 3.963671
#>  [9,] 0.6208898 0.4354921 3.770757 3.035116 1.8869789 2.767437 3.960406
#> [10,] 0.6190550 0.4322034 3.849667 2.964441 1.9036960 2.864581 3.957361
#> [11,] 0.6215051 0.4371410 3.913337 2.909372 1.9176751 2.951722 3.954266
#> [12,] 0.6218366 0.4359885 3.957297 2.861531 1.9109909 3.033778 3.951858
```

We can also perform conditional forecasting by following Algorithm 3.3.1
in Dieppe, Legrand, and van Roye (2016). Note that for the structural
shocks, identification is based on the Cholesky factorisation. Also,
please note the limitations of this method, see the detailed discussion
in Section 5.4 of Dieppe, Legrand, and van Roye (2016).

Now suppose we are interested in the forecasts of the domestic interest
rate \\i\\ conditional on a scenario where domestic inflation \\\pi\\
gets really high (post COVID type scenario). Economic theory says the
short interest rate should rise.

First we set up our conditions/scenarios, i.e., which variables, which
horizons, and which values the variables will take during those
horizons. Our conditions are that \\\pi\\ will follow a specified path,
at forecast horizons \\h=1,\dots,H=12\\.

``` r

conditions <- data.frame(
              var        = rep(5,12),
              horizon    = rep(1:12),
              value      = c(1.0,1.5,2.0,1.8, #Note: QoQ scale for inflation here
                             1.5,1.2,1.0,1.0,
                             rep(0.5,4))
              )
```

We then do the conditional forecasting. We again select a 95% PI
(prediction interval) and the mean of the predictive distribution as the
point forecast.

``` r

cond_fcst <- conditional_forecast(bvar_obj,
                                  conditions,
                                  pi=0.95,
                                  fcst_type = "mean",
                                  plot_idx = c(5,6),
                                  growth_rate_idx = c(5))
```

![plot of chunk HOMO-5](figure/HOMO-5-1.png)![plot of chunk
HOMO-5](figure/HOMO-5-2.png)

The short interest rate rises more dramatically compared to the
unconditional case. Makes sense.

Now for some impulse response analysis. We can choose between the
orthogonalized impulse response function (OIRF) and the generalized
impulse response function (GIRF). Similar to forecasting, we can choose
either the mean or the median (the default is the median), and we can
also transform the IRFs for the quarter-on-quarter growth rate variables
to the annual/yearly scale.

``` r

irf <- IRF(bvar_obj,H=20,response=5,impulse=6,type="median",method="OIRF",ci=0.95,growth_rate_idx=5)
```

![plot of chunk HOMO-6](figure/HOMO-6-1.png)

plot of chunk HOMO-6

``` r

irf <- IRF(bvar_obj,H=20,response=4,impulse=6,type="median",method="GIRF",ci=0.95,growth_rate_idx=4)
```

![plot of chunk HOMO-6](figure/HOMO-6-2.png)

plot of chunk HOMO-6

## References

Dieppe, A., Legrand, R., and van Roye, B. (2016). The BEAR toolbox.
*Working Paper Series*, No. 1934. European Central Bank.

Karlsson, S. (2013). Forecasting with Bayesian vector autoregression.
In: Elliott, G. and Timmermann, A. (eds), *Handbook of Economic
Forecasting*. Elsevier B.V., Vol. 2, Part B, pp. 791–897.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630-650.
