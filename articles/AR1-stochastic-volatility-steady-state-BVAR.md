# AR(1) stochastic volatility steady-state BVAR

Here we estimate a steady-state BVAR model with AR(1) stochastic
volatility, see
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

![plot of chunk AR(1)-1](figure/AR(1)-1-1.png)

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
for details on the model). Below I take inspiration from Carriero,
Clark, and Marcellino (2024), which uses the exact same AR(1) stochastic
volatility specification, but for a conventional BVAR.

``` r

k <- bvar_obj$setup$k
n_free_params_A <- bvar_obj$setup$n_free_params_A
sigma2 <- diag(bvar_obj$setup$Sigma_AR)

SV_priors_AR1 <- list(
                      theta_A            =  rep(0, n_free_params_A),
                      Omega_A            =  diag(10, n_free_params_A),
                      theta_gamma_0      =  0.1 * log(sigma2),
                      Omega_gamma_0      =  diag(2, k),
                      theta_gamma_1      =  rep(0.9, k),
                      Omega_gamma_1      =  diag(0.04, k),
                      theta_log_lambda_1 =  log(sigma2),
                      Omega_log_lambda_1 =  diag(2, k),
                      V_Phi              = (10 - k - 1) * 0.03 * diag(k),
                      m_Phi              =  10
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
                   SV_type = "AR1",
                   SV_priors = SV_priors_AR1)
```

We can plot the steady-state priors

``` r

par(mfrow=c(3,1))
steady_state_priors_plot(bvar_obj, interval = 0.95)
```

![plot of chunk AR(1)-2](figure/AR(1)-2-1.png)

``` r

par(mfrow=c(1,1))
```

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
                control = list(max_treedepth = 12, adapt_delta = 0.999))
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
#> steady_state_bvar_AR1_stochastic_volatility
#> 
#> Also generating draws from the joint predictive distribution
#> 
#> ...
#> SAMPLING FINISHED
```

Now lets see the elementwise posterior means

``` r

summary(bvar_obj, stat="mean", t = 215) #t = 215 for covariance matrix
#> Posterior mean estimates
#> ------------------------
#> 
#> 
#> beta
#> --------------------------------------------------------------------------------             
#>               delta pi     u     r
#>   delta pi.l1     1.27  0.02  0.16
#>   u.l1           -0.09  1.17 -0.16
#>   r.l1            0.00 -0.01  1.04
#>   delta pi.l2    -0.28  0.01 -0.11
#>   u.l2            0.07 -0.23  0.17
#>   r.l2           -0.01  0.02 -0.11
#> --------------------------------------------------------------------------------
#> 
#> 
#> Psi
#> --------------------------------------------------------------------------------          
#>            [,1]
#>   delta pi 1.99
#>   u        4.28
#>   r        3.51
#> --------------------------------------------------------------------------------
#> 
#> 
#> Sigma_u,t (t = 215)
#> --------------------------------------------------------------------------------
#>          delta pi     u     r
#> delta pi     0.07 -0.01  0.02
#> u           -0.01  0.03 -0.02
#> r            0.02 -0.02  0.16
#> --------------------------------------------------------------------------------
#> 
#> 
#> A
#> --------------------------------------------------------------------------------          
#>            delta pi    u r
#>   delta pi     1.00 0.00 0
#>   u            0.13 1.00 0
#>   r           -0.24 0.45 1
#> --------------------------------------------------------------------------------
#> 
#> 
#> gamma_0
#> --------------------------------------------------------------------------------
#> delta pi        u        r 
#>    -0.14    -0.15    -0.09 
#> --------------------------------------------------------------------------------
#> 
#> 
#> gamma_1
#> --------------------------------------------------------------------------------
#> delta pi        u        r 
#>     0.95     0.95     0.94 
#> --------------------------------------------------------------------------------
#> 
#> 
#> Phi
#> --------------------------------------------------------------------------------          
#>            delta pi    u    r
#>   delta pi     0.06 0.04 0.06
#>   u            0.04 0.07 0.07
#>   r            0.06 0.07 0.14
#> --------------------------------------------------------------------------------
```

You can always look at the `stanfit` object `bvar_obj$fit$stan` directly
if you want. Note that the `z`’s below are not parameters per se, they
are simply used in a reparameterization trick to sample the log
volatilities more efficiently.

``` r

print(bvar_obj$fit$stan)
#> Inference for Stan model: steady_state_bvar_AR1_stochastic_volatility.
#> 4 chains, each with iter=6000; warmup=2000; thin=1; 
#> post-warmup draws per chain=4000, total post-warmup draws=16000.
#> 
#>                        mean se_mean    sd  2.5%   25%   50%   75%  97.5% n_eff Rhat
#> beta[1,1]              1.27    0.00  0.06  1.16  1.23  1.27  1.31   1.38 11980    1
#> beta[1,2]              0.02    0.00  0.04 -0.05  0.00  0.02  0.05   0.10 12982    1
#> beta[1,3]              0.16    0.00  0.08  0.00  0.11  0.17  0.22   0.33 12305    1
#> beta[2,1]             -0.09    0.00  0.04 -0.16 -0.11 -0.09 -0.06  -0.02 13526    1
#> beta[2,2]              1.17    0.00  0.06  1.06  1.13  1.17  1.20   1.28 12496    1
#> beta[2,3]             -0.16    0.00  0.08 -0.31 -0.21 -0.15 -0.10   0.00 12001    1
#> beta[3,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.02   0.04 14612    1
#> beta[3,2]             -0.01    0.00  0.02 -0.05 -0.02 -0.01  0.00   0.02 14117    1
#> beta[3,3]              1.04    0.00  0.06  0.92  1.00  1.04  1.08   1.16 13548    1
#> beta[4,1]             -0.28    0.00  0.06 -0.39 -0.32 -0.28 -0.24  -0.17 11861    1
#> beta[4,2]              0.01    0.00  0.04 -0.07 -0.02  0.01  0.03   0.08 13328    1
#> beta[4,3]             -0.11    0.00  0.08 -0.27 -0.17 -0.11 -0.05   0.05 12374    1
#> beta[5,1]              0.07    0.00  0.03  0.00  0.05  0.07  0.09   0.13 13719    1
#> beta[5,2]             -0.23    0.00  0.05 -0.33 -0.26 -0.23 -0.19  -0.12 12676    1
#> beta[5,3]              0.17    0.00  0.07  0.03  0.12  0.17  0.22   0.31 12366    1
#> beta[6,1]             -0.01    0.00  0.02 -0.04 -0.02 -0.01  0.00   0.02 14863    1
#> beta[6,2]              0.02    0.00  0.02 -0.01  0.01  0.02  0.04   0.06 14611    1
#> beta[6,3]             -0.11    0.00  0.06 -0.22 -0.15 -0.11 -0.07   0.01 13603    1
#> Psi[1,1]               1.99    0.00  0.05  1.89  1.96  2.00  2.03   2.09 31831    1
#> Psi[2,1]               4.28    0.00  0.17  3.93  4.16  4.28  4.40   4.62 24174    1
#> Psi[3,1]               3.51    0.00  0.33  2.84  3.29  3.51  3.74   4.16 23581    1
#> z[1,1]                 0.01    0.00  0.42 -0.75 -0.28 -0.02  0.28   0.89 18694    1
#> z[1,2]                 1.25    0.00  0.43  0.44  0.96  1.23  1.52   2.14 16625    1
#> z[1,3]                -1.01    0.00  0.53 -2.02 -1.37 -1.03 -0.67   0.07 15833    1
#> z[2,1]                 0.00    0.01  0.95 -1.87 -0.64 -0.01  0.63   1.86 24859    1
#> z[2,2]                 0.23    0.01  0.99 -1.73 -0.44  0.24  0.92   2.16 24788    1
#> z[2,3]                -0.08    0.01  0.98 -2.03 -0.74 -0.07  0.58   1.84 28618    1
#> z[3,1]                 0.11    0.01  0.95 -1.77 -0.54  0.11  0.75   1.99 28257    1
#> z[3,2]                 0.28    0.01  0.99 -1.66 -0.39  0.28  0.95   2.24 26282    1
#> z[3,3]                -0.14    0.01  0.98 -2.05 -0.79 -0.13  0.52   1.77 27484    1
#> z[4,1]                -0.38    0.01  0.97 -2.26 -1.04 -0.38  0.27   1.51 27972    1
#> z[4,2]                -0.22    0.01  0.96 -2.11 -0.88 -0.22  0.44   1.65 28932    1
#> z[4,3]                -0.25    0.01  0.97 -2.16 -0.91 -0.25  0.40   1.66 27232    1
#> z[5,1]                -0.25    0.01  0.97 -2.16 -0.90 -0.26  0.41   1.63 25858    1
#> z[5,2]                -0.22    0.01  0.99 -2.15 -0.90 -0.23  0.44   1.76 29208    1
#> z[5,3]                -0.23    0.01  0.98 -2.17 -0.89 -0.23  0.41   1.70 29541    1
#> z[6,1]                -0.16    0.01  0.96 -2.06 -0.81 -0.15  0.49   1.73 28198    1
#> z[6,2]                -0.10    0.01  0.98 -2.00 -0.77 -0.10  0.55   1.81 25721    1
#> z[6,3]                -0.18    0.01  0.98 -2.13 -0.83 -0.18  0.47   1.76 27216    1
#> z[7,1]                 0.08    0.01  0.95 -1.79 -0.55  0.09  0.72   1.96 29519    1
#> z[7,2]                 0.02    0.01  0.95 -1.83 -0.63  0.02  0.66   1.89 27521    1
#> z[7,3]                -0.14    0.01  0.99 -2.06 -0.82 -0.14  0.56   1.79 28849    1
#> z[8,1]                 0.28    0.01  0.94 -1.56 -0.36  0.27  0.91   2.11 25623    1
#> z[8,2]                 0.03    0.01  0.97 -1.88 -0.62  0.03  0.70   1.94 26831    1
#> z[8,3]                -0.08    0.01  0.97 -1.96 -0.74 -0.07  0.57   1.82 29210    1
#> z[9,1]                 0.44    0.01  0.96 -1.44 -0.20  0.44  1.08   2.32 24838    1
#> z[9,2]                 0.21    0.01  0.96 -1.67 -0.43  0.20  0.86   2.07 27585    1
#> z[9,3]                 0.02    0.01  0.99 -1.95 -0.65  0.01  0.68   1.96 29395    1
#> z[10,1]               -0.16    0.01  0.95 -2.01 -0.79 -0.16  0.48   1.69 30738    1
#> z[10,2]                0.14    0.01  0.98 -1.76 -0.52  0.15  0.80   2.06 29941    1
#> z[10,3]               -0.13    0.01  0.97 -2.05 -0.78 -0.13  0.51   1.76 30649    1
#> z[11,1]               -0.21    0.01  0.95 -2.06 -0.86 -0.21  0.43   1.64 29054    1
#> z[11,2]                0.10    0.01  0.98 -1.82 -0.55  0.10  0.77   2.01 25965    1
#> z[11,3]               -0.21    0.01  0.97 -2.11 -0.86 -0.21  0.45   1.69 27749    1
#> z[12,1]               -0.35    0.01  0.97 -2.24 -1.00 -0.34  0.29   1.57 27819    1
#> z[12,2]                0.21    0.01  0.98 -1.69 -0.44  0.22  0.87   2.13 27857    1
#> z[12,3]               -0.19    0.01  0.98 -2.10 -0.85 -0.19  0.46   1.72 30699    1
#> z[13,1]               -0.12    0.01  0.96 -2.00 -0.77 -0.12  0.53   1.75 29850    1
#> z[13,2]                0.39    0.01  0.98 -1.51 -0.28  0.38  1.05   2.30 26817    1
#> z[13,3]               -0.08    0.01  0.98 -1.97 -0.74 -0.08  0.57   1.83 29240    1
#> z[14,1]               -0.20    0.01  0.96 -2.09 -0.85 -0.19  0.45   1.66 28963    1
#> z[14,2]                0.41    0.01  0.96 -1.47 -0.24  0.42  1.05   2.26 28184    1
#> z[14,3]               -0.13    0.01  0.99 -2.06 -0.79 -0.12  0.53   1.80 29759    1
#> z[15,1]               -0.14    0.01  0.96 -1.99 -0.80 -0.14  0.51   1.70 27709    1
#> z[15,2]                0.36    0.01  0.97 -1.56 -0.30  0.36  1.01   2.26 24546    1
#> z[15,3]               -0.07    0.01  0.98 -1.99 -0.75 -0.07  0.61   1.85 29165    1
#> z[16,1]                0.06    0.01  0.94 -1.78 -0.58  0.06  0.71   1.91 27660    1
#> z[16,2]                0.48    0.01  0.98 -1.44 -0.18  0.49  1.14   2.40 28565    1
#> z[16,3]                0.01    0.01  0.99 -1.91 -0.65  0.01  0.68   1.96 31975    1
#> z[17,1]                0.11    0.01  0.96 -1.79 -0.54  0.11  0.76   2.00 26726    1
#> z[17,2]                0.54    0.01  0.97 -1.34 -0.12  0.53  1.19   2.45 32123    1
#> z[17,3]                0.08    0.01  0.98 -1.84 -0.58  0.08  0.74   2.01 26493    1
#> z[18,1]                0.17    0.01  0.96 -1.74 -0.47  0.17  0.80   2.05 25629    1
#> z[18,2]                0.71    0.01  0.98 -1.21  0.06  0.71  1.39   2.63 23293    1
#> z[18,3]                0.16    0.01  1.00 -1.79 -0.51  0.16  0.83   2.14 25946    1
#> z[19,1]                0.38    0.01  0.96 -1.52 -0.26  0.38  1.03   2.23 22396    1
#> z[19,2]                0.83    0.01  0.97 -1.05  0.19  0.83  1.48   2.72 23028    1
#> z[19,3]                0.15    0.01  1.00 -1.78 -0.53  0.14  0.82   2.09 27214    1
#> z[20,1]                0.12    0.01  0.95 -1.76 -0.51  0.12  0.76   1.98 25875    1
#> z[20,2]                0.51    0.01  0.96 -1.39 -0.15  0.52  1.17   2.38 29838    1
#> z[20,3]                0.13    0.01  0.98 -1.77 -0.52  0.13  0.79   2.05 34035    1
#> z[21,1]               -0.05    0.01  0.95 -1.90 -0.70 -0.06  0.59   1.78 30750    1
#> z[21,2]                0.08    0.01  0.97 -1.81 -0.57  0.08  0.74   1.97 28597    1
#> z[21,3]                0.17    0.01  0.98 -1.74 -0.48  0.18  0.83   2.08 28119    1
#> z[22,1]                0.14    0.01  0.96 -1.73 -0.52  0.13  0.79   2.02 25186    1
#> z[22,2]                0.27    0.01  0.95 -1.61 -0.37  0.27  0.90   2.15 26359    1
#> z[22,3]                0.29    0.01  0.99 -1.66 -0.38  0.29  0.96   2.26 25704    1
#> z[23,1]               -0.24    0.01  0.95 -2.10 -0.88 -0.23  0.40   1.61 26772    1
#> z[23,2]               -0.01    0.01  0.97 -1.89 -0.67 -0.01  0.64   1.91 25603    1
#> z[23,3]               -0.11    0.01  0.98 -2.03 -0.76 -0.10  0.54   1.80 30490    1
#> z[24,1]               -0.22    0.01  0.95 -2.09 -0.85 -0.22  0.42   1.65 28798    1
#> z[24,2]                0.10    0.01  0.96 -1.82 -0.55  0.11  0.74   1.96 31039    1
#> z[24,3]               -0.04    0.01  0.99 -1.99 -0.70 -0.04  0.62   1.94 29854    1
#> z[25,1]               -0.09    0.01  0.94 -1.95 -0.73 -0.08  0.54   1.75 26514    1
#> z[25,2]                0.15    0.01  0.97 -1.73 -0.50  0.15  0.80   2.06 26796    1
#> z[25,3]                0.05    0.01  0.97 -1.87 -0.59  0.06  0.70   1.95 29991    1
#> z[26,1]                0.20    0.01  0.96 -1.66 -0.46  0.20  0.85   2.06 26101    1
#> z[26,2]                0.31    0.01  0.97 -1.60 -0.33  0.30  0.97   2.22 24536    1
#> z[26,3]                0.16    0.01  0.96 -1.76 -0.47  0.16  0.81   2.03 27307    1
#> z[27,1]               -0.17    0.01  0.92 -1.98 -0.79 -0.17  0.45   1.65 27949    1
#>  [ reached 'max' / getOption("max.print") -- omitted 6459 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Mon Sep 14 04:50:28 2026.
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

![plot of chunk AR(1)-3](figure/AR(1)-3-1.png)

We can forecast

``` r

par(mfrow=c(3,1))
forecast(bvar_obj, pi = 0.68)
```

![plot of chunk AR(1)-4](figure/AR(1)-4-1.png)

``` r

par(mfrow=c(1,1))
```

Let us plot the log volatility estimates and predictions

``` r

stochastic_volatility_plot(bvar_obj, ci = 0.95, vol = "log_lambda")
```

![plot of chunk AR(1)-5](figure/AR(1)-5-1.png)

![plot of chunk AR(1)-5](figure/AR(1)-5-2.png)

![plot of chunk AR(1)-5](figure/AR(1)-5-3.png)

Let us plot the estimates and predictions of the implied innovation
standard deviations

``` r

stochastic_volatility_plot(bvar_obj, vol = "sd")
```

![plot of chunk AR(1)-6](figure/AR(1)-6-1.png)

![plot of chunk AR(1)-6](figure/AR(1)-6-2.png)

![plot of chunk AR(1)-6](figure/AR(1)-6-3.png)

We can also produce orthogonalized IRFs

``` r

IRF(bvar_obj, method = "OIRF", t=215, ci=0.68) #using Sigma_u,t=215
```

![plot of chunk AR(1)-7](figure/AR(1)-7-1.png)

## References

Carriero, A., Clark, T. E., and Marcellino, M. (2024). Capturing
macro-economic tail risks with Bayesian vector autoregressions. *Journal
of Money, Credit and Banking*, 56(5), pp. 1099–1127.

Koop, G. and Korobilis, D. (2010). Bayesian multivariate time series
methods for empirical macroeconomics. *Foundations and Trends in
Econometrics*, 3(4), pp. 267–358.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630–650.
