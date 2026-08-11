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
steady-state priors are expert based). Remember that we only have a
constant now, so \\q=1\\ and therefore \\\Psi\\ only has one column
\\\psi_1=\Psi\\. Since \\d_t = 1 \\ \forall \\ t\\, we have \\\Psi d_t =
\mu_t\\ which simplifies to \\\Psi = \mu\\ and as such we can directly
interpret \\\Psi\\ as the unconditional mean, i.e. the steady state.

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

Now we can fit the model. Note that we can use arguments from
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
such as `control` where we can tweak `max_treedepth` and `adapt_delta`.

``` r

bvar_obj <- fit(bvar_obj,
                H = 40,
                iter = 2000,
                warmup = 1000,
                chains = 2,
                cores = 2,
                control = list(max_treedepth = 12, adapt_delta = 0.999))
#> ------------------------------------------------------------
#> Estimating Stan model:
#> steady_state_bvar_AR1_stochastic_volatility
#> 
#> Also generating draws from the joint predictive distribution
#> 
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
```

Now lets see the posterior means

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
#>   r.l2           -0.01  0.02 -0.11
#> --------------------------------------------------------------------------------
#> 
#> 
#> Psi
#> --------------------------------------------------------------------------------          
#>            [,1]
#>   delta pi 1.99
#>   u        4.29
#>   r        3.51
#> --------------------------------------------------------------------------------
#> 
#> 
#> Sigma_u,t (t = 215)
#> --------------------------------------------------------------------------------
#>          delta pi     u     r
#> delta pi     0.07 -0.01  0.02
#> u           -0.01  0.03 -0.02
#> r            0.02 -0.02  0.17
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
#>     0.95     0.95     0.93 
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
#> 2 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=2000.
#> 
#>                        mean se_mean    sd  2.5%   25%   50%   75%  97.5% n_eff Rhat
#> beta[1,1]              1.27    0.00  0.05  1.16  1.24  1.27  1.31   1.38  1799 1.00
#> beta[1,2]              0.02    0.00  0.04 -0.05  0.00  0.02  0.05   0.10  2157 1.00
#> beta[1,3]              0.17    0.00  0.08  0.01  0.11  0.17  0.22   0.33  1883 1.00
#> beta[2,1]             -0.09    0.00  0.04 -0.16 -0.11 -0.09 -0.06  -0.02  1868 1.00
#> beta[2,2]              1.17    0.00  0.06  1.05  1.13  1.17  1.21   1.28  1412 1.00
#> beta[2,3]             -0.15    0.00  0.08 -0.30 -0.21 -0.16 -0.10   0.00  1075 1.00
#> beta[3,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.02   0.04  1964 1.00
#> beta[3,2]             -0.01    0.00  0.02 -0.04 -0.02 -0.01  0.00   0.02  1729 1.00
#> beta[3,3]              1.04    0.00  0.06  0.92  1.00  1.04  1.08   1.16  1442 1.00
#> beta[4,1]             -0.28    0.00  0.06 -0.39 -0.32 -0.28 -0.24  -0.17  1801 1.00
#> beta[4,2]              0.01    0.00  0.04 -0.07 -0.02  0.01  0.04   0.09  2123 1.00
#> beta[4,3]             -0.11    0.00  0.08 -0.27 -0.17 -0.11 -0.06   0.04  1932 1.00
#> beta[5,1]              0.07    0.00  0.03  0.00  0.05  0.07  0.09   0.14  1775 1.00
#> beta[5,2]             -0.23    0.00  0.05 -0.33 -0.26 -0.23 -0.19  -0.12  1425 1.00
#> beta[5,3]              0.17    0.00  0.07  0.03  0.12  0.17  0.22   0.31  1075 1.00
#> beta[6,1]             -0.01    0.00  0.02 -0.04 -0.02 -0.01  0.01   0.02  2022 1.00
#> beta[6,2]              0.02    0.00  0.02 -0.01  0.01  0.02  0.04   0.06  1783 1.00
#> beta[6,3]             -0.11    0.00  0.06 -0.22 -0.15 -0.11 -0.07   0.00  1447 1.00
#> Psi[1,1]               1.99    0.00  0.05  1.89  1.96  1.99  2.03   2.10  3266 1.00
#> Psi[2,1]               4.29    0.00  0.18  3.94  4.16  4.29  4.41   4.61  3268 1.00
#> Psi[3,1]               3.51    0.01  0.33  2.86  3.28  3.51  3.74   4.13  3426 1.00
#> z[1,1]                 0.01    0.01  0.42 -0.79 -0.28  0.00  0.28   0.89  2243 1.00
#> z[1,2]                 1.23    0.01  0.43  0.39  0.94  1.22  1.52   2.10  2163 1.00
#> z[1,3]                -1.02    0.01  0.52 -1.99 -1.37 -1.05 -0.68   0.03  2379 1.00
#> z[2,1]                 0.00    0.02  0.93 -1.84 -0.60  0.03  0.62   1.77  3169 1.00
#> z[2,2]                 0.27    0.02  0.96 -1.61 -0.39  0.26  0.93   2.19  2531 1.00
#> z[2,3]                -0.08    0.01  0.96 -1.99 -0.72 -0.07  0.54   1.77  4496 1.00
#> z[3,1]                 0.15    0.02  1.00 -1.85 -0.56  0.13  0.83   2.13  3559 1.00
#> z[3,2]                 0.28    0.02  1.01 -1.68 -0.40  0.30  0.96   2.22  4006 1.00
#> z[3,3]                -0.16    0.02  0.99 -2.02 -0.83 -0.15  0.53   1.75  3418 1.00
#> z[4,1]                -0.39    0.02  0.95 -2.20 -1.06 -0.41  0.28   1.47  3831 1.00
#> z[4,2]                -0.24    0.02  0.98 -2.14 -0.90 -0.25  0.45   1.67  3738 1.00
#> z[4,3]                -0.23    0.02  0.98 -2.22 -0.84 -0.21  0.46   1.64  3260 1.00
#> z[5,1]                -0.28    0.02  0.94 -2.06 -0.92 -0.27  0.34   1.59  2953 1.00
#> z[5,2]                -0.20    0.01  0.96 -2.06 -0.85 -0.24  0.43   1.74  4465 1.00
#> z[5,3]                -0.21    0.01  0.95 -2.08 -0.87 -0.20  0.46   1.62  4674 1.00
#> z[6,1]                -0.13    0.02  0.97 -2.00 -0.83 -0.11  0.54   1.75  4130 1.00
#> z[6,2]                -0.11    0.02  0.98 -2.02 -0.77 -0.10  0.55   1.77  3598 1.00
#> z[6,3]                -0.15    0.02  1.01 -2.11 -0.84 -0.14  0.52   1.84  3162 1.00
#> z[7,1]                 0.08    0.02  1.01 -1.94 -0.60  0.07  0.75   2.12  3687 1.00
#> z[7,2]                 0.03    0.02  0.97 -1.89 -0.61  0.03  0.67   1.87  2928 1.00
#> z[7,3]                -0.14    0.02  0.98 -1.99 -0.82 -0.14  0.54   1.80  3721 1.00
#> z[8,1]                 0.25    0.02  0.92 -1.52 -0.38  0.25  0.86   2.06  3285 1.00
#> z[8,2]                 0.04    0.02  0.98 -1.86 -0.62  0.02  0.69   2.01  3275 1.00
#> z[8,3]                -0.09    0.01  0.92 -1.84 -0.73 -0.08  0.55   1.66  4493 1.00
#> z[9,1]                 0.40    0.02  0.93 -1.42 -0.23  0.40  1.05   2.11  2890 1.00
#> z[9,2]                 0.22    0.02  0.97 -1.55 -0.49  0.19  0.91   2.10  2613 1.00
#> z[9,3]                 0.01    0.02  0.99 -1.88 -0.67  0.00  0.68   1.93  4358 1.00
#> z[10,1]               -0.15    0.01  0.91 -1.89 -0.78 -0.16  0.50   1.67  3770 1.00
#> z[10,2]                0.12    0.02  0.94 -1.73 -0.54  0.14  0.75   1.93  2629 1.00
#> z[10,3]               -0.14    0.02  1.01 -2.18 -0.82 -0.12  0.51   1.91  3614 1.00
#> z[11,1]               -0.19    0.02  0.93 -2.01 -0.83 -0.20  0.45   1.57  3167 1.00
#> z[11,2]                0.10    0.02  0.96 -1.82 -0.55  0.09  0.74   1.98  3519 1.00
#> z[11,3]               -0.20    0.02  0.99 -2.15 -0.87 -0.20  0.50   1.72  3802 1.00
#> z[12,1]               -0.34    0.02  0.98 -2.25 -1.01 -0.32  0.31   1.61  2739 1.00
#> z[12,2]                0.19    0.02  0.98 -1.73 -0.46  0.20  0.82   2.11  3007 1.00
#> z[12,3]               -0.17    0.02  1.00 -2.13 -0.84 -0.15  0.48   1.80  3059 1.00
#> z[13,1]               -0.09    0.01  0.95 -1.88 -0.73 -0.10  0.53   1.75  4830 1.00
#> z[13,2]                0.40    0.01  0.95 -1.47 -0.24  0.40  1.05   2.35  4200 1.00
#> z[13,3]               -0.12    0.02  0.99 -2.00 -0.77 -0.13  0.53   1.82  4361 1.00
#> z[14,1]               -0.18    0.01  0.91 -2.02 -0.81 -0.15  0.42   1.63  3788 1.00
#> z[14,2]                0.40    0.02  0.95 -1.41 -0.25  0.38  1.06   2.21  3064 1.00
#> z[14,3]               -0.12    0.02  1.01 -2.07 -0.84 -0.12  0.58   1.86  3454 1.00
#> z[15,1]               -0.14    0.01  0.95 -2.08 -0.78 -0.11  0.49   1.69  4009 1.00
#> z[15,2]                0.35    0.02  0.93 -1.40 -0.27  0.36  0.95   2.27  3659 1.00
#> z[15,3]               -0.11    0.02  1.02 -2.18 -0.77 -0.10  0.58   1.88  4413 1.00
#> z[16,1]                0.04    0.02  0.97 -1.91 -0.62  0.04  0.72   1.89  3451 1.00
#> z[16,2]                0.48    0.02  0.96 -1.43 -0.15  0.50  1.11   2.35  3832 1.00
#> z[16,3]                0.01    0.02  0.96 -1.87 -0.64  0.00  0.63   1.94  3344 1.00
#> z[17,1]                0.11    0.02  0.96 -1.73 -0.53  0.10  0.78   1.89  2827 1.00
#> z[17,2]                0.53    0.02  0.97 -1.33 -0.14  0.52  1.18   2.48  3055 1.00
#> z[17,3]                0.05    0.02  0.96 -1.74 -0.60  0.06  0.69   1.95  4011 1.00
#> z[18,1]                0.20    0.02  0.97 -1.67 -0.47  0.23  0.87   2.04  3278 1.00
#> z[18,2]                0.75    0.02  0.99 -1.17  0.09  0.74  1.43   2.70  3400 1.00
#> z[18,3]                0.18    0.02  1.01 -1.81 -0.52  0.16  0.85   2.18  3529 1.00
#> z[19,1]                0.41    0.02  0.95 -1.45 -0.23  0.42  1.07   2.26  3165 1.00
#> z[19,2]                0.87    0.02  0.95 -1.00  0.23  0.88  1.52   2.74  2273 1.00
#> z[19,3]                0.11    0.02  0.97 -1.76 -0.53  0.12  0.75   1.99  4015 1.00
#> z[20,1]                0.11    0.02  0.91 -1.70 -0.51  0.09  0.74   1.91  2982 1.00
#> z[20,2]                0.50    0.02  0.97 -1.44 -0.17  0.49  1.18   2.38  3109 1.00
#> z[20,3]                0.13    0.02  0.97 -1.79 -0.51  0.11  0.78   2.10  3015 1.00
#> z[21,1]               -0.08    0.02  0.94 -1.92 -0.69 -0.08  0.54   1.79  3457 1.00
#> z[21,2]                0.09    0.02  0.96 -1.77 -0.53  0.10  0.72   2.00  2770 1.00
#> z[21,3]                0.17    0.02  0.99 -1.76 -0.49  0.16  0.83   2.11  3561 1.00
#> z[22,1]                0.14    0.01  0.93 -1.71 -0.47  0.16  0.76   1.97  4529 1.00
#> z[22,2]                0.28    0.02  0.98 -1.67 -0.38  0.29  0.94   2.26  3351 1.00
#> z[22,3]                0.29    0.02  0.99 -1.66 -0.40  0.29  0.97   2.21  3114 1.00
#> z[23,1]               -0.25    0.02  0.96 -2.07 -0.92 -0.25  0.41   1.64  3062 1.00
#> z[23,2]                0.00    0.02  0.94 -1.80 -0.66  0.00  0.64   1.88  3258 1.00
#> z[23,3]               -0.09    0.02  0.98 -1.95 -0.77 -0.12  0.60   1.84  3412 1.00
#> z[24,1]               -0.25    0.02  0.92 -2.05 -0.90 -0.22  0.41   1.49  3038 1.00
#> z[24,2]                0.10    0.01  0.96 -1.86 -0.55  0.08  0.73   2.00  4305 1.00
#> z[24,3]               -0.05    0.02  0.98 -1.94 -0.71 -0.05  0.60   1.86  3503 1.00
#> z[25,1]               -0.10    0.01  0.96 -1.96 -0.74 -0.11  0.56   1.78  5129 1.00
#> z[25,2]                0.14    0.02  0.99 -1.83 -0.54  0.15  0.79   2.12  2873 1.00
#> z[25,3]                0.08    0.02  0.96 -1.74 -0.57  0.09  0.73   1.96  3612 1.00
#> z[26,1]                0.22    0.02  0.94 -1.60 -0.41  0.23  0.86   1.97  3104 1.00
#> z[26,2]                0.31    0.02  0.98 -1.62 -0.36  0.31  1.00   2.23  3673 1.00
#> z[26,3]                0.15    0.02  1.03 -1.87 -0.54  0.17  0.82   2.19  3224 1.00
#> z[27,1]               -0.15    0.02  0.93 -1.95 -0.79 -0.15  0.47   1.68  2831 1.00
#>  [ reached 'max' / getOption("max.print") -- omitted 5700 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Wed Aug 12 00:31:27 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

We can forecast

``` r

forecast(bvar_obj, pi = 0.68, show_all = TRUE)
```

![plot of chunk AR(1)-2](figure/AR(1)-2-1.png)![plot of chunk
AR(1)-2](figure/AR(1)-2-2.png)![plot of chunk
AR(1)-2](figure/AR(1)-2-3.png)

Let us plot the log volatility estimates and predictions

``` r

stochastic_volatility_plot(bvar_obj, ci = 0.95, vol = "log_lambda")
```

![plot of chunk AR(1)-3](figure/AR(1)-3-1.png)![plot of chunk
AR(1)-3](figure/AR(1)-3-2.png)![plot of chunk
AR(1)-3](figure/AR(1)-3-3.png)

Let us plot the estimates and predictions of the implied innovation
standard deviations

``` r

stochastic_volatility_plot(bvar_obj, vol = "sd")
```

![plot of chunk AR(1)-4](figure/AR(1)-4-1.png)![plot of chunk
AR(1)-4](figure/AR(1)-4-2.png)![plot of chunk
AR(1)-4](figure/AR(1)-4-3.png)

We can also produce orthogonalized IRFs

``` r

IRF(bvar_obj, method = "OIRF", t=215, ci=0.68) #latest t
```

![plot of chunk AR(1)-5](figure/AR(1)-5%20-1.png)

plot of chunk AR(1)-5

## References

Carriero, A., Clark, T. E., and Marcellino, M. (2024). Capturing
macro-economic tail risks with Bayesian vector autoregressions. *Journal
of Money, Credit and Banking*, 56(5), pp. 1099–1127.

Koop, G. and Korobilis, D. (2010). Bayesian multivariate time series
methods for empirical macroeconomics. *Foundations and Trends in
Econometrics*, 3(4), pp. 267–358.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630–650.
