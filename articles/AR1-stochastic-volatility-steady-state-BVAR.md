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

plot of chunk AR(1)-1

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
steady-state priors are expert based). Since we only have a constant
now, we have \\q=1\\ and as such \\\Psi\\ has only one column
\\\psi_1=\Psi\\. Since \\d_t = 1 \\ \forall \\ t\\, we have \\\Psi d_t =
\mu_t\\ which simplifies to \\\Psi = \mu\\ and as such we can directly
interpret \\\Psi d_t = \mu_t = \Psi = \mu\\ as the unconditional mean,
i.e. the steady state.

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

steady_state_priors_plot(bvar_obj, interval = 0.95)
```

![plot of chunk AR(1)-2](figure/AR(1)-2-1.png)![plot of chunk
AR(1)-2](figure/AR(1)-2-2.png)![plot of chunk
AR(1)-2](figure/AR(1)-2-3.png)

Now, let us fit the model. Note that we can use arguments from
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
such as `control` where we can tweak `max_treedepth` and `adapt_delta`.
Let us choose 4 markov chains, with each having 10000 iterations, and
where 2500 of those 10000 are warmup/burn-in iterations.

``` r

bvar_obj <- fit(bvar_obj,
                H = 40,
                iter = 5000,
                warmup = 2500,
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
#>   delta pi.l1     1.27  0.02  0.17
#>   u.l1           -0.09  1.17 -0.15
#>   r.l1            0.00 -0.01  1.04
#>   delta pi.l2    -0.28  0.01 -0.11
#>   u.l2            0.07 -0.23  0.17
#>   r.l2           -0.01  0.03 -0.11
#> --------------------------------------------------------------------------------
#> 
#> 
#> Psi
#> --------------------------------------------------------------------------------          
#>            [,1]
#>   delta pi 2.00
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
#> 4 chains, each with iter=5000; warmup=2500; thin=1; 
#> post-warmup draws per chain=2500, total post-warmup draws=10000.
#> 
#>                        mean se_mean    sd  2.5%   25%   50%   75%  97.5% n_eff Rhat
#> beta[1,1]              1.27    0.00  0.06  1.16  1.23  1.27  1.31   1.38  7370    1
#> beta[1,2]              0.02    0.00  0.04 -0.05  0.00  0.02  0.05   0.10  9114    1
#> beta[1,3]              0.17    0.00  0.08  0.00  0.11  0.17  0.22   0.32  7469    1
#> beta[2,1]             -0.09    0.00  0.04 -0.16 -0.11 -0.09 -0.06  -0.02  8410    1
#> beta[2,2]              1.17    0.00  0.06  1.06  1.13  1.17  1.20   1.28  7253    1
#> beta[2,3]             -0.15    0.00  0.08 -0.31 -0.21 -0.16 -0.10   0.00  8122    1
#> beta[3,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.01   0.04  9926    1
#> beta[3,2]             -0.01    0.00  0.02 -0.05 -0.02 -0.01  0.00   0.02  8889    1
#> beta[3,3]              1.04    0.00  0.06  0.92  1.00  1.04  1.09   1.16  8741    1
#> beta[4,1]             -0.28    0.00  0.06 -0.39 -0.32 -0.28 -0.24  -0.17  7459    1
#> beta[4,2]              0.01    0.00  0.04 -0.06 -0.02  0.01  0.04   0.08  9013    1
#> beta[4,3]             -0.11    0.00  0.08 -0.27 -0.17 -0.11 -0.06   0.05  7716    1
#> beta[5,1]              0.07    0.00  0.03  0.00  0.05  0.07  0.09   0.13  8910    1
#> beta[5,2]             -0.23    0.00  0.05 -0.33 -0.26 -0.23 -0.19  -0.12  7384    1
#> beta[5,3]              0.17    0.00  0.07  0.03  0.12  0.17  0.22   0.31  8275    1
#> beta[6,1]             -0.01    0.00  0.02 -0.04 -0.02 -0.01  0.01   0.03 10100    1
#> beta[6,2]              0.03    0.00  0.02 -0.01  0.01  0.03  0.04   0.06  9505    1
#> beta[6,3]             -0.11    0.00  0.06 -0.23 -0.15 -0.11 -0.07   0.00  9002    1
#> Psi[1,1]               2.00    0.00  0.05  1.89  1.96  2.00  2.03   2.09 18550    1
#> Psi[2,1]               4.28    0.00  0.18  3.93  4.16  4.28  4.40   4.62 16909    1
#> Psi[3,1]               3.51    0.00  0.33  2.85  3.29  3.51  3.73   4.15 16690    1
#> z[1,1]                 0.01    0.00  0.41 -0.75 -0.28 -0.01  0.28   0.85 12304    1
#> z[1,2]                 1.26    0.00  0.43  0.45  0.96  1.25  1.53   2.15 10548    1
#> z[1,3]                -1.01    0.01  0.53 -2.01 -1.36 -1.03 -0.66   0.08 10915    1
#> z[2,1]                -0.01    0.01  0.96 -1.88 -0.67 -0.02  0.65   1.87 14492    1
#> z[2,2]                 0.23    0.01  0.98 -1.68 -0.44  0.23  0.91   2.13 16682    1
#> z[2,3]                -0.07    0.01  1.00 -2.03 -0.74 -0.06  0.59   1.88 17628    1
#> z[3,1]                 0.11    0.01  0.95 -1.74 -0.55  0.12  0.74   1.99 16541    1
#> z[3,2]                 0.28    0.01  0.97 -1.59 -0.38  0.29  0.93   2.14 16478    1
#> z[3,3]                -0.13    0.01  0.96 -2.01 -0.77 -0.13  0.52   1.75 17749    1
#> z[4,1]                -0.38    0.01  0.96 -2.25 -1.04 -0.38  0.27   1.47 17893    1
#> z[4,2]                -0.23    0.01  0.98 -2.18 -0.88 -0.23  0.44   1.70 19556    1
#> z[4,3]                -0.27    0.01  0.98 -2.16 -0.93 -0.26  0.41   1.62 17963    1
#> z[5,1]                -0.24    0.01  0.98 -2.15 -0.91 -0.24  0.42   1.67 18317    1
#> z[5,2]                -0.21    0.01  0.96 -2.13 -0.86 -0.22  0.44   1.66 18755    1
#> z[5,3]                -0.22    0.01  0.98 -2.13 -0.88 -0.23  0.44   1.68 17166    1
#> z[6,1]                -0.16    0.01  0.98 -2.06 -0.83 -0.16  0.52   1.75 19658    1
#> z[6,2]                -0.10    0.01  0.96 -1.96 -0.75 -0.10  0.56   1.76 17882    1
#> z[6,3]                -0.18    0.01  0.99 -2.12 -0.86 -0.19  0.49   1.76 17920    1
#> z[7,1]                 0.08    0.01  0.95 -1.80 -0.56  0.09  0.72   1.97 16350    1
#> z[7,2]                 0.02    0.01  0.99 -1.93 -0.66  0.03  0.69   1.95 18693    1
#> z[7,3]                -0.13    0.01  0.95 -2.01 -0.78 -0.13  0.50   1.76 18211    1
#> z[8,1]                 0.27    0.01  0.96 -1.59 -0.38  0.28  0.92   2.15 17077    1
#> z[8,2]                 0.05    0.01  0.96 -1.87 -0.59  0.04  0.69   1.95 17790    1
#> z[8,3]                -0.09    0.01  0.98 -2.05 -0.74 -0.08  0.56   1.81 18357    1
#> z[9,1]                 0.44    0.01  0.94 -1.42 -0.20  0.43  1.08   2.31 15451    1
#> z[9,2]                 0.21    0.01  0.99 -1.70 -0.46  0.21  0.89   2.15 19271    1
#> z[9,3]                 0.02    0.01  0.97 -1.83 -0.64  0.01  0.69   1.90 18231    1
#> z[10,1]               -0.16    0.01  0.92 -1.94 -0.77 -0.14  0.46   1.67 20028    1
#> z[10,2]                0.14    0.01  0.98 -1.79 -0.52  0.15  0.79   2.06 21187    1
#> z[10,3]               -0.15    0.01  0.95 -2.05 -0.79 -0.16  0.49   1.75 17952    1
#> z[11,1]               -0.21    0.01  0.94 -2.04 -0.83 -0.22  0.42   1.63 20526    1
#> z[11,2]                0.10    0.01  0.99 -1.85 -0.56  0.09  0.76   2.03 18079    1
#> z[11,3]               -0.20    0.01  0.96 -2.10 -0.86 -0.21  0.47   1.66 20079    1
#> z[12,1]               -0.33    0.01  0.96 -2.22 -0.98 -0.33  0.30   1.56 18799    1
#> z[12,2]                0.21    0.01  0.98 -1.71 -0.46  0.21  0.86   2.15 20968    1
#> z[12,3]               -0.18    0.01  0.99 -2.11 -0.84 -0.18  0.50   1.78 17817    1
#> z[13,1]               -0.13    0.01  0.93 -1.90 -0.76 -0.14  0.50   1.69 18555    1
#> z[13,2]                0.38    0.01  0.96 -1.48 -0.26  0.38  1.01   2.27 19315    1
#> z[13,3]               -0.10    0.01  0.96 -2.00 -0.75 -0.10  0.55   1.78 18969    1
#> z[14,1]               -0.19    0.01  0.97 -2.09 -0.84 -0.20  0.46   1.71 18477    1
#> z[14,2]                0.40    0.01  0.96 -1.46 -0.25  0.41  1.04   2.26 18940    1
#> z[14,3]               -0.12    0.01  0.98 -2.03 -0.77 -0.13  0.54   1.76 18906    1
#> z[15,1]               -0.14    0.01  0.96 -1.99 -0.80 -0.14  0.53   1.73 19276    1
#> z[15,2]                0.35    0.01  0.97 -1.56 -0.29  0.35  1.00   2.26 19510    1
#> z[15,3]               -0.08    0.01  0.97 -1.97 -0.75 -0.09  0.58   1.79 19754    1
#> z[16,1]                0.06    0.01  0.96 -1.84 -0.58  0.07  0.69   1.93 16323    1
#> z[16,2]                0.48    0.01  0.99 -1.46 -0.18  0.48  1.15   2.43 17730    1
#> z[16,3]                0.00    0.01  0.97 -1.89 -0.65  0.00  0.65   1.90 17826    1
#> z[17,1]                0.12    0.01  0.96 -1.74 -0.53  0.12  0.78   2.00 18321    1
#> z[17,2]                0.55    0.01  0.99 -1.38 -0.12  0.54  1.22   2.47 20164    1
#> z[17,3]                0.06    0.01  0.97 -1.83 -0.59  0.05  0.71   1.96 17932    1
#> z[18,1]                0.19    0.01  0.97 -1.71 -0.47  0.19  0.85   2.09 15458    1
#> z[18,2]                0.72    0.01  0.98 -1.19  0.05  0.72  1.38   2.65 16942    1
#> z[18,3]                0.15    0.01  0.96 -1.76 -0.48  0.14  0.80   2.05 15946    1
#> z[19,1]                0.37    0.01  0.99 -1.57 -0.30  0.38  1.05   2.26 13709    1
#> z[19,2]                0.84    0.01  0.97 -1.06  0.18  0.84  1.50   2.76 17053    1
#> z[19,3]                0.14    0.01  0.99 -1.80 -0.53  0.15  0.81   2.06 18849    1
#> z[20,1]                0.13    0.01  0.96 -1.75 -0.51  0.14  0.77   2.00 16626    1
#> z[20,2]                0.51    0.01  0.94 -1.36 -0.11  0.51  1.16   2.35 16282    1
#> z[20,3]                0.12    0.01  0.97 -1.75 -0.53  0.12  0.77   2.00 18496    1
#> z[21,1]               -0.06    0.01  0.93 -1.89 -0.68 -0.05  0.58   1.77 17813    1
#> z[21,2]                0.10    0.01  1.00 -1.86 -0.58  0.10  0.79   2.03 18883    1
#> z[21,3]                0.18    0.01  0.98 -1.76 -0.47  0.19  0.84   2.11 17709    1
#> z[22,1]                0.15    0.01  0.94 -1.67 -0.49  0.16  0.79   1.98 16883    1
#> z[22,2]                0.27    0.01  0.98 -1.66 -0.38  0.28  0.93   2.19 19036    1
#> z[22,3]                0.29    0.01  0.99 -1.65 -0.38  0.30  0.94   2.26 20193    1
#> z[23,1]               -0.25    0.01  0.94 -2.11 -0.89 -0.25  0.40   1.60 18205    1
#> z[23,2]               -0.02    0.01  0.95 -1.89 -0.67 -0.01  0.62   1.85 17980    1
#> z[23,3]               -0.12    0.01  0.99 -2.04 -0.81 -0.12  0.56   1.82 20144    1
#> z[24,1]               -0.21    0.01  0.95 -2.09 -0.85 -0.20  0.42   1.66 19149    1
#> z[24,2]                0.10    0.01  0.96 -1.80 -0.54  0.10  0.73   1.99 16405    1
#> z[24,3]               -0.05    0.01  0.98 -1.96 -0.70 -0.04  0.61   1.86 18385    1
#> z[25,1]               -0.08    0.01  0.96 -1.96 -0.74 -0.09  0.56   1.81 17266    1
#> z[25,2]                0.15    0.01  0.95 -1.74 -0.49  0.16  0.78   2.02 18633    1
#> z[25,3]                0.06    0.01  0.97 -1.86 -0.60  0.06  0.71   1.95 19376    1
#> z[26,1]                0.22    0.01  0.94 -1.62 -0.42  0.22  0.85   2.07 15466    1
#> z[26,2]                0.32    0.01  0.95 -1.50 -0.34  0.31  0.96   2.20 16746    1
#> z[26,3]                0.16    0.01  0.99 -1.81 -0.51  0.16  0.82   2.08 21405    1
#> z[27,1]               -0.15    0.01  0.95 -2.02 -0.78 -0.16  0.49   1.74 18791    1
#>  [ reached 'max' / getOption("max.print") -- omitted 6459 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Thu Sep 10 19:59:56 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Let us plot the posterior of the steady-state \\mu_t=mu\\

``` r

stanfit <- bvar_obj$fit$stan

rstan::plot(stanfit, pars=c("mu[1,1]",
                            "mu[1,2]",
                            "mu[1,3]"), plotfun="hist")
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![plot of chunk AR(1)-3](figure/AR(1)-3-1.png)

plot of chunk AR(1)-3

We can forecast

``` r

forecast(bvar_obj, pi = 0.68)
```

![plot of chunk AR(1)-4](figure/AR(1)-4-1.png)![plot of chunk
AR(1)-4](figure/AR(1)-4-2.png)![plot of chunk
AR(1)-4](figure/AR(1)-4-3.png)

Let us plot the log volatility estimates and predictions

``` r

stochastic_volatility_plot(bvar_obj, ci = 0.95, vol = "log_lambda")
```

![plot of chunk AR(1)-5](figure/AR(1)-5-1.png)![plot of chunk
AR(1)-5](figure/AR(1)-5-2.png)![plot of chunk
AR(1)-5](figure/AR(1)-5-3.png)

Let us plot the estimates and predictions of the implied innovation
standard deviations

``` r

stochastic_volatility_plot(bvar_obj, vol = "sd")
```

![plot of chunk AR(1)-6](figure/AR(1)-6-1.png)![plot of chunk
AR(1)-6](figure/AR(1)-6-2.png)![plot of chunk
AR(1)-6](figure/AR(1)-6-3.png)

We can also produce orthogonalized IRFs

``` r

IRF(bvar_obj, method = "OIRF", t=215, ci=0.68) #latest t
```

![plot of chunk AR(1)-7](figure/AR(1)-7-1.png)

plot of chunk AR(1)-7

## References

Carriero, A., Clark, T. E., and Marcellino, M. (2024). Capturing
macro-economic tail risks with Bayesian vector autoregressions. *Journal
of Money, Credit and Banking*, 56(5), pp. 1099–1127.

Koop, G. and Korobilis, D. (2010). Bayesian multivariate time series
methods for empirical macroeconomics. *Foundations and Trends in
Econometrics*, 3(4), pp. 267–358.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630–650.
