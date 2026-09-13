# Random Walk stochastic volatility steady-state BVAR (Clark, 2011)

Here we estimate the steady-state BVAR model with Random Walk stochastic
volatility from Clark (2011), which is an extension of the original
homoscedastic steady-state BVAR model (Villani, 2009). See
[`?bvar`](https://markjwbecker.github.io/SteadyStateBVAR/reference/bvar.md)
for details.

We will estimate the model on a quarterly US data set from Koop and
Korobilis (2010) on the inflation rate \\\Delta \pi_t\\ (the annual
percentage change in a chain-weighted GDP price index), the unemployment
rate \\u_t\\ (seasonally adjusted civilian unemployment rate, all
civilian workers aged 16 years or older) and the interest rate \\r_t\\
(yield on the three-month Treasury bill rate). The sample is
1953Q1-2006Q3 and we have the data vector

\\ y_t = \begin{pmatrix} \Delta \pi_t \\ u_t \\ r_t \end{pmatrix} \\

First, let’s load the package, then import and plot the data.

``` r

library(SteadyStateBVAR)
data("KoopKorobilis2010")
yt <- KoopKorobilis2010
plot.ts(yt)
```

![plot of chunk RW-1](figure/RW-1-1.png)

plot of chunk RW-1

Let’s create the bvar object which we will use throughout here.

``` r

bvar_obj <- bvar(data = yt)
```

We choose 2 lags and only a constant as the deterministic variable.

``` r

bvar_obj <- setup(bvar_obj,
                  p=2,
                  deterministic = "constant")
```

We set the overall tightness to \\\lambda_1 = 0.20\\, cross-equation
tightness to \\\lambda_2 = 0.50\\ and the lag decay rate to \\\lambda_3
= 1.00\\. For the prior means on the first own lags, we set them to
\\0.6\\ for \\\Delta \pi_t\\ and \\0.9\\ for \\u_t\\ and \\r_t\\. Note
that the prior mean on the first own lag of inflation is set to \\0.6\\
instead of \\0\\ to reflect some degree of persistence in the series
(even though it is a growth rate variable).

``` r

lambda_1 <- 0.20
lambda_2 <- 0.50
lambda_3 <- 1.00

fol_pm=c(0.6, # delta pi
         0.9,  #u
         0.9)  #R
```

Now, for the steady-state priors, we specify our prior beliefs with the
help of 95% prior probability intervals (let us pretend that our
steady-state priors are expert based). Remember that the steady-state is
\\\Psi d_t = \mu_t\\, where \\\Psi\\ is \\k \times q\\ and \\d_t\\ is
\\q \times 1\\. Since we only have a constant now, i.e. \\d_t = 1 \\
\forall \\ t\\, we have \\q=1\\ and as such we can directly interpret
the \\k\\-dimensional vector \\\Psi d_t = \mu_t = \mu = \Psi\\ as the
unconditional mean, i.e. the steady state.

``` r

theta_Psi <- 
  c(
  ppi(1.90, 2.10, interval=0.95)$mean,   #Psi: delta pi
  ppi(3.80, 4.50, interval=0.95)$mean,   #Psi: u
  ppi(2.60, 3.90, interval=0.95)$mean    #Psi: r
  )

Omega_Psi <- 
  diag(
  c(
  ppi(1.90, 2.10, interval=0.95)$var,    #Psi: delta pi
  ppi(3.80, 4.50, interval=0.95)$var,    #Psi: u
  ppi(2.60, 3.90, interval=0.95)$var     #Psi: r
  )
  )
```

Now we need to specify our stochastic volatility priors. See
[`?priors`](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md)
for more information about the prior specification (and
[`?bvar`](https://markjwbecker.github.io/SteadyStateBVAR/reference/bvar.md)
or
[`vignette("SteadyStateBVAR-intro")`](https://markjwbecker.github.io/SteadyStateBVAR/articles/SteadyStateBVAR-intro.md)
for details on the model). I take some inspiration from Clark (2011)
below.

``` r

k <- bvar_obj$setup$k
n_free_params_A <- bvar_obj$setup$n_free_params_A
sigma2 <- diag(bvar_obj$setup$Sigma_AR)

SV_priors_RW <- list(
                     theta_A             =  rep(0, n_free_params_A),
                     Omega_A             =  diag(10, n_free_params_A),
                     theta_log_lambda_1  =  log(sigma2), #initial condition of
                     Omega_log_lambda_1  =  diag(4, k),
                     alpha_phi           =  rep(2.5, k),
                     beta_phi            =  rep(0.0875, k)
                    )
```

Here `sigma2` contains the residual variances from AR(\\p\\) models (the
same ones we used in the Minnesota prior).

Let’s put everything into the
[`priors()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md)
function. Please note that `lambda_1`, `lambda_2`, and `lambda_3`, which
are the hyperparameters in the Minnesota prior, have nothing to do with
the \\\ln \lambda\\’s (log volatilities) from the stochastic volatility
specification.

``` r

bvar_obj <- priors(bvar_obj,
                   lambda_1 = lambda_1,
                   lambda_2 = lambda_2,
                   lambda_3 = lambda_3,
                   first_own_lag_prior_mean =fol_pm,
                   theta_Psi = theta_Psi,
                   Omega_Psi = Omega_Psi,
                   SV = TRUE,
                   SV_type = "RW",
                   SV_priors = SV_priors_RW)
```

We can plot the steady-state priors

``` r

steady_state_priors_plot(bvar_obj, interval = 0.95)
```

![plot of chunk RW-2](figure/RW-2-1.png)

plot of chunk RW-2

![plot of chunk RW-2](figure/RW-2-2.png)

plot of chunk RW-2

![plot of chunk RW-2](figure/RW-2-3.png)

plot of chunk RW-2

Now, let us fit the model. Note that we can use arguments from
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
such as `control` where we can tweak `max_treedepth` and `adapt_delta`.
Let us choose 4 markov chains, with each having 6000 iterations, and
where 2000 of those 6000 are warmup/burn-in iterations.

``` r

bvar_obj <- fit(bvar_obj,
                H = 40,
                iter = 6000,
                warmup = 2000,
                chains = 4,
                cores = 4,
                control = list(max_treedepth = 14, adapt_delta = 0.95))
#> ------------------------------------------------------------
#> Forecast horizon:
#> 40
#> 
#> Future deterministic variables (d_pred):
#>      constant
#> h=1         1
#> h=2         1
#> h=3         1
#> h=4         1
#> h=5         1
#> h=6         1
#> h=7         1
#> h=8         1
#> h=9         1
#> h=10        1
#> h=11        1
#> h=12        1
#> h=13        1
#> h=14        1
#> h=15        1
#> h=16        1
#> h=17        1
#> h=18        1
#> h=19        1
#> h=20        1
#> h=21        1
#> h=22        1
#> h=23        1
#> h=24        1
#> h=25        1
#> h=26        1
#> h=27        1
#> h=28        1
#> h=29        1
#> h=30        1
#> h=31        1
#> h=32        1
#> h=33        1
#> h=34        1
#> h=35        1
#> h=36        1
#> h=37        1
#> h=38        1
#> h=39        1
#> h=40        1
#> ------------------------------------------------------------
#> Estimating Stan model:
#> steady_state_bvar_RW_stochastic_volatility
#> 
#> Also generating draws from the joint predictive distribution
#> 
#> ...
#> Warning: There were 1 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> SAMPLING FINISHED
```

Now let’s see the elemntwise posterior means

``` r

summary(bvar_obj, stat="mean", t = 215) #t = 215 for covariance matrix
#> Posterior mean estimates
#> ------------------------
#> 
#> 
#> beta
#> --------------------------------------------------------------------------------             
#>               delta pi     u     r
#>   delta pi.l1     1.27  0.02  0.15
#>   u.l1           -0.09  1.17 -0.16
#>   r.l1            0.00 -0.01  1.04
#>   delta pi.l2    -0.28  0.02 -0.10
#>   u.l2            0.07 -0.23  0.17
#>   r.l2            0.00  0.02 -0.11
#> --------------------------------------------------------------------------------
#> 
#> 
#> Psi
#> --------------------------------------------------------------------------------          
#>            [,1]
#>   delta pi 2.00
#>   u        4.29
#>   r        3.49
#> --------------------------------------------------------------------------------
#> 
#> 
#> Sigma_u,t (t = 215)
#> --------------------------------------------------------------------------------
#>          delta pi     u     r
#> delta pi     0.09 -0.01  0.02
#> u           -0.01  0.02 -0.01
#> r            0.02 -0.01  0.16
#> --------------------------------------------------------------------------------
#> 
#> 
#> A
#> --------------------------------------------------------------------------------          
#>            delta pi    u r
#>   delta pi     1.00 0.00 0
#>   u            0.12 1.00 0
#>   r           -0.23 0.51 1
#> --------------------------------------------------------------------------------
#> 
#> 
#> phi
#> --------------------------------------------------------------------------------
#> delta pi        u        r 
#>     0.04     0.07     0.10 
#> --------------------------------------------------------------------------------
```

You can always look at the `stanfit` object `bvar_obj$fit$stan` directly
if you want. Note that the `z`’s below are not parameters per se, they
are simply used in a reparameterization trick to sample the log
volatilities more efficiently.

``` r

print(bvar_obj$fit$stan)
#> Inference for Stan model: steady_state_bvar_RW_stochastic_volatility.
#> 4 chains, each with iter=6000; warmup=2000; thin=1; 
#> post-warmup draws per chain=4000, total post-warmup draws=16000.
#> 
#>                        mean se_mean    sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
#> beta[1,1]              1.27    0.00  0.06  1.16  1.23  1.27  1.31  1.38 19262    1
#> beta[1,2]              0.02    0.00  0.04 -0.06 -0.01  0.02  0.04  0.09 19833    1
#> beta[1,3]              0.15    0.00  0.08 -0.01  0.09  0.15  0.20  0.31 20015    1
#> beta[2,1]             -0.09    0.00  0.03 -0.16 -0.12 -0.09 -0.07 -0.03 18995    1
#> beta[2,2]              1.17    0.00  0.06  1.06  1.14  1.17  1.21  1.28 16026    1
#> beta[2,3]             -0.16    0.00  0.08 -0.31 -0.21 -0.16 -0.11 -0.01 17441    1
#> beta[3,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.01  0.04 20316    1
#> beta[3,2]             -0.01    0.00  0.02 -0.04 -0.02 -0.01  0.00  0.02 19952    1
#> beta[3,3]              1.04    0.00  0.06  0.92  1.00  1.04  1.08  1.16 20305    1
#> beta[4,1]             -0.28    0.00  0.06 -0.39 -0.32 -0.28 -0.24 -0.17 19495    1
#> beta[4,2]              0.02    0.00  0.04 -0.06 -0.01  0.02  0.04  0.09 20512    1
#> beta[4,3]             -0.10    0.00  0.08 -0.26 -0.16 -0.10 -0.05  0.06 20278    1
#> beta[5,1]              0.07    0.00  0.03  0.01  0.05  0.07  0.09  0.14 19502    1
#> beta[5,2]             -0.23    0.00  0.05 -0.34 -0.27 -0.23 -0.20 -0.13 16082    1
#> beta[5,3]              0.17    0.00  0.07  0.03  0.12  0.17  0.22  0.31 18363    1
#> beta[6,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.01  0.03 21521    1
#> beta[6,2]              0.02    0.00  0.02 -0.01  0.01  0.02  0.03  0.05 20081    1
#> beta[6,3]             -0.11    0.00  0.06 -0.22 -0.15 -0.11 -0.07  0.01 21161    1
#> Psi[1,1]               2.00    0.00  0.05  1.90  1.96  2.00  2.03  2.10 32519    1
#> Psi[2,1]               4.29    0.00  0.18  3.94  4.17  4.29  4.41  4.63 25576    1
#> Psi[3,1]               3.49    0.00  0.32  2.85  3.28  3.50  3.71  4.13 27551    1
#> z[1,1]                -0.06    0.00  0.24 -0.51 -0.22 -0.07  0.09  0.43 23273    1
#> z[1,2]                 0.64    0.00  0.27  0.13  0.46  0.63  0.81  1.17 22280    1
#> z[1,3]                -0.65    0.00  0.34 -1.29 -0.88 -0.66 -0.43  0.04 20622    1
#> z[2,1]                -0.12    0.01  0.98 -2.06 -0.77 -0.13  0.53  1.77 30547    1
#> z[2,2]                 0.20    0.01  0.98 -1.72 -0.46  0.21  0.87  2.10 27268    1
#> z[2,3]                 0.01    0.01  0.97 -1.88 -0.65  0.02  0.67  1.93 26719    1
#> z[3,1]                -0.02    0.01  0.97 -1.94 -0.67 -0.02  0.63  1.88 28474    1
#> z[3,2]                 0.26    0.01  0.99 -1.67 -0.40  0.26  0.93  2.20 27599    1
#> z[3,3]                -0.04    0.01  0.98 -1.93 -0.70 -0.04  0.62  1.89 31015    1
#> z[4,1]                -0.09    0.01  0.99 -2.04 -0.76 -0.09  0.58  1.83 33604    1
#> z[4,2]                -0.43    0.01  0.97 -2.30 -1.08 -0.43  0.23  1.47 31443    1
#> z[4,3]                -0.26    0.01  0.96 -2.15 -0.91 -0.26  0.39  1.63 31876    1
#> z[5,1]                 0.00    0.01  0.97 -1.91 -0.65  0.01  0.66  1.88 30265    1
#> z[5,2]                -0.45    0.01  0.97 -2.32 -1.10 -0.45  0.21  1.45 33327    1
#> z[5,3]                -0.21    0.01  0.98 -2.10 -0.86 -0.22  0.45  1.70 29480    1
#> z[6,1]                 0.01    0.01  0.96 -1.89 -0.64  0.01  0.65  1.89 30483    1
#> z[6,2]                -0.34    0.01  0.96 -2.22 -0.98 -0.34  0.30  1.56 29699    1
#> z[6,3]                -0.12    0.01  0.96 -1.99 -0.77 -0.11  0.53  1.79 28297    1
#> z[7,1]                 0.10    0.01  0.98 -1.82 -0.57  0.10  0.75  2.03 30465    1
#> z[7,2]                -0.23    0.01  0.98 -2.15 -0.88 -0.23  0.44  1.68 31450    1
#> z[7,3]                -0.03    0.01  0.96 -1.92 -0.68 -0.02  0.62  1.85 31752    1
#> z[8,1]                 0.18    0.01  0.99 -1.77 -0.49  0.19  0.85  2.09 31341    1
#> z[8,2]                -0.24    0.01  0.96 -2.11 -0.89 -0.25  0.41  1.63 30865    1
#> z[8,3]                 0.06    0.01  0.94 -1.80 -0.58  0.06  0.71  1.90 30989    1
#> z[9,1]                 0.17    0.01  0.97 -1.76 -0.49  0.17  0.82  2.07 32711    1
#> z[9,2]                -0.13    0.01  0.97 -2.06 -0.78 -0.13  0.54  1.74 32157    1
#> z[9,3]                 0.20    0.01  0.98 -1.71 -0.46  0.20  0.86  2.14 32542    1
#> z[10,1]               -0.19    0.01  0.97 -2.09 -0.84 -0.20  0.46  1.74 29066    1
#> z[10,2]               -0.11    0.01  0.97 -2.00 -0.77 -0.11  0.56  1.75 31008    1
#> z[10,3]               -0.05    0.01  0.96 -1.91 -0.69 -0.05  0.60  1.83 29214    1
#> z[11,1]               -0.17    0.01  0.96 -2.01 -0.82 -0.16  0.48  1.68 30092    1
#> z[11,2]               -0.13    0.01  0.97 -2.00 -0.79 -0.13  0.52  1.79 29564    1
#> z[11,3]               -0.16    0.01  0.96 -2.04 -0.81 -0.16  0.48  1.72 31904    1
#> z[12,1]               -0.34    0.01  0.96 -2.22 -0.99 -0.33  0.32  1.51 31180    1
#> z[12,2]               -0.04    0.01  0.97 -1.93 -0.69 -0.04  0.61  1.85 30183    1
#> z[12,3]               -0.11    0.01  0.96 -1.99 -0.75 -0.11  0.53  1.77 31099    1
#> z[13,1]               -0.27    0.01  0.96 -2.15 -0.92 -0.27  0.38  1.60 31581    1
#> z[13,2]                0.09    0.01  0.96 -1.79 -0.55  0.09  0.74  1.99 31111    1
#> z[13,3]                0.03    0.01  0.96 -1.84 -0.62  0.03  0.68  1.90 30644    1
#> z[14,1]               -0.33    0.01  0.98 -2.24 -0.99 -0.33  0.34  1.57 30312    1
#> z[14,2]                0.13    0.01  0.97 -1.78 -0.53  0.13  0.79  2.02 34045    1
#> z[14,3]                0.03    0.01  0.98 -1.88 -0.63  0.03  0.70  1.96 30712    1
#> z[15,1]               -0.26    0.01  0.97 -2.18 -0.90 -0.26  0.40  1.67 31216    1
#> z[15,2]                0.07    0.01  0.95 -1.81 -0.55  0.08  0.71  1.93 33463    1
#> z[15,3]                0.12    0.01  0.94 -1.72 -0.52  0.12  0.75  1.97 30643    1
#> z[16,1]               -0.21    0.01  0.98 -2.13 -0.89 -0.21  0.45  1.70 31966    1
#> z[16,2]                0.12    0.01  0.96 -1.75 -0.53  0.11  0.76  2.00 29182    1
#> z[16,3]                0.26    0.01  0.96 -1.63 -0.38  0.27  0.91  2.14 28305    1
#> z[17,1]               -0.21    0.01  0.97 -2.13 -0.87 -0.21  0.44  1.68 33780    1
#> z[17,2]                0.16    0.01  0.98 -1.75 -0.52  0.16  0.83  2.05 31229    1
#> z[17,3]                0.32    0.01  0.97 -1.60 -0.33  0.33  0.97  2.22 33603    1
#> z[18,1]               -0.27    0.01  0.96 -2.15 -0.91 -0.28  0.38  1.61 31085    1
#> z[18,2]                0.26    0.01  0.95 -1.61 -0.39  0.26  0.90  2.11 31355    1
#> z[18,3]                0.42    0.01  0.96 -1.48 -0.24  0.41  1.07  2.31 31658    1
#> z[19,1]               -0.19    0.01  0.98 -2.13 -0.85 -0.20  0.47  1.75 31358    1
#> z[19,2]                0.38    0.01  0.96 -1.51 -0.26  0.38  1.02  2.26 28866    1
#> z[19,3]                0.37    0.01  0.96 -1.51 -0.27  0.37  1.01  2.26 33031    1
#> z[20,1]               -0.19    0.01  0.97 -2.11 -0.83 -0.18  0.46  1.71 33247    1
#> z[20,2]                0.03    0.01  0.94 -1.83 -0.60  0.03  0.66  1.88 25517    1
#> z[20,3]                0.31    0.01  0.94 -1.53 -0.32  0.31  0.95  2.17 31832    1
#> z[21,1]               -0.10    0.01  0.97 -2.02 -0.76 -0.10  0.56  1.80 27198    1
#> z[21,2]               -0.44    0.01  0.95 -2.32 -1.07 -0.43  0.20  1.43 31950    1
#> z[21,3]                0.37    0.01  0.93 -1.47 -0.26  0.37  0.99  2.19 32052    1
#> z[22,1]               -0.10    0.01  0.96 -2.02 -0.73 -0.10  0.55  1.78 34435    1
#> z[22,2]               -0.32    0.01  0.96 -2.17 -0.96 -0.31  0.33  1.55 33773    1
#> z[22,3]                0.51    0.01  0.93 -1.33 -0.12  0.51  1.14  2.33 27665    1
#> z[23,1]               -0.05    0.01  0.97 -1.95 -0.71 -0.05  0.62  1.84 30335    1
#> z[23,2]               -0.36    0.01  0.96 -2.23 -1.00 -0.36  0.28  1.52 31059    1
#> z[23,3]               -0.17    0.01  0.95 -2.04 -0.81 -0.18  0.46  1.71 31927    1
#> z[24,1]               -0.19    0.01  0.98 -2.11 -0.85 -0.18  0.47  1.72 35873    1
#> z[24,2]               -0.28    0.01  0.97 -2.19 -0.93 -0.28  0.37  1.59 31372    1
#> z[24,3]               -0.07    0.01  0.96 -1.94 -0.72 -0.06  0.58  1.79 29360    1
#> z[25,1]               -0.17    0.01  0.97 -2.09 -0.82 -0.17  0.49  1.77 31985    1
#> z[25,2]               -0.30    0.01  0.96 -2.17 -0.94 -0.30  0.34  1.58 31226    1
#> z[25,3]                0.06    0.01  0.93 -1.77 -0.56  0.06  0.68  1.87 27878    1
#> z[26,1]               -0.07    0.01  0.96 -1.95 -0.72 -0.06  0.58  1.80 30509    1
#> z[26,2]               -0.19    0.01  0.96 -2.12 -0.84 -0.19  0.46  1.67 29309    1
#> z[26,3]                0.20    0.01  0.94 -1.65 -0.44  0.19  0.84  2.01 32103    1
#> z[27,1]               -0.08    0.01  0.97 -2.00 -0.74 -0.08  0.57  1.80 31427    1
#>  [ reached 'max' / getOption("max.print") -- omitted 6438 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Sun Sep 13 17:00:01 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Let us plot the posterior steady-state. We can see that \\\mu_t=\Psi\\
in the constant only case, as discussed. Note that in the Stan code the
\\\mu_t\\’s are stacked row-wise in the \\T \times k\\ matrix `mu`.

``` r

stanfit <- bvar_obj$fit$stan

rstan::plot(stanfit, pars=c("mu[1,1]",
                            "mu[1,2]",
                            "mu[1,3]",
                            "Psi[1,1]",
                            "Psi[2,1]",
                            "Psi[3,1]"), plotfun="hist")
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![plot of chunk RW-3](figure/RW-3-1.png)

plot of chunk RW-3

We can forecast

``` r

forecast(bvar_obj, pi = 0.68)
```

![plot of chunk RW-4](figure/RW-4-1.png)

plot of chunk RW-4

![plot of chunk RW-4](figure/RW-4-2.png)

plot of chunk RW-4

![plot of chunk RW-4](figure/RW-4-3.png)

plot of chunk RW-4

Let us plot the log volatility estimates and predictions

``` r

stochastic_volatility_plot(bvar_obj, ci = 0.95, vol = "log_lambda")
```

![plot of chunk RW-5](figure/RW-5-1.png)

plot of chunk RW-5

![plot of chunk RW-5](figure/RW-5-2.png)

plot of chunk RW-5

![plot of chunk RW-5](figure/RW-5-3.png)

plot of chunk RW-5

Let us plot the estimates and predictions of the implied reduced-form
innovation standard deviations

``` r

stochastic_volatility_plot(bvar_obj, vol = "sd")
```

![plot of chunk RW-6](figure/RW-6-1.png)

plot of chunk RW-6

![plot of chunk RW-6](figure/RW-6-2.png)

plot of chunk RW-6

![plot of chunk RW-6](figure/RW-6-3.png)

plot of chunk RW-6

We can also produce orthogonalized IRFs

``` r

IRF(bvar_obj, method = "OIRF", t=215, ci=0.68) #latest t
```

![plot of chunk RW-7](figure/RW-7-1.png)

plot of chunk RW-7

## References

Clark, T. E. (2011). Real-time density forecasts from Bayesian vector
autoregressions with stochastic volatility. *Journal of Business &
Economic Statistics*, 29(3), pp. 327–341.

Koop, G. and Korobilis, D. (2010). Bayesian multivariate time series
methods for empirical macroeconomics. *Foundations and Trends in
Econometrics*, 3(4), pp. 267–358.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630–650.
