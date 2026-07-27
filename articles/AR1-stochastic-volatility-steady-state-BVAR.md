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

![](figure/AR(1)-1-1.png)

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

Now, for the steady-state coefficients we use some toy values (let us
pretend that they are expert based). Remember that we only have a
constant now, so \\q=1\\ and therefore \\\Psi\\ only has one column
\\\psi_1=\Psi\\. Since \\d_t = 1 \\ \forall \\ t\\, we have \\\Psi d_t =
\mu_t\\ which simplifies to \\\Psi = \mu\\ and as such we can directly
interpret \\\Psi\\ as the unconditional mean, i.e. the steady-state.

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
for more information about the prior specification. Below I take some
inspiration from Carriero, Clark, and Marcellino (2024), which uses the
exact same AR(1) stochastic volatility specification, but for a
conventional BVAR.

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
                      V_Phi              = (7 - k - 1) * 0.5 * diag(k),
                      m_Phi              =  7
                     )
```

Here `sigma2` contains the residual variances from AR(\\p\\) models (the
same ones we used in the Minnesota prior).

Let’s put everything into the
[`priors()`](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md)
function.

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
                d_pred = matrix(rep(1, 40)),
                iter = 4000,
                warmup = 1000,
                chains = 2,
                cores = 2,
                control = list(max_treedepth = 12, adapt_delta = 0.999))
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
#>   delta pi.l1     1.26  0.01  0.15
#>   u.l1           -0.09  1.15 -0.16
#>   r.l1            0.00 -0.01  1.04
#>   delta pi.l2    -0.27  0.02 -0.10
#>   u.l2            0.07 -0.22  0.17
#>   r.l2           -0.01  0.03 -0.10
#> --------------------------------------------------------------------------------
#> 
#> 
#> Psi
#> --------------------------------------------------------------------------------          
#>            [,1]
#>   delta pi 1.99
#>   u        4.27
#>   r        3.53
#> --------------------------------------------------------------------------------
#> 
#> 
#> Sigma_u,t (t = 215)
#> --------------------------------------------------------------------------------
#>          delta pi     u     r
#> delta pi     0.09 -0.01  0.03
#> u           -0.01  0.03 -0.02
#> r            0.03 -0.02  0.21
#> --------------------------------------------------------------------------------
#> 
#> 
#> A
#> --------------------------------------------------------------------------------          
#>            delta pi    u r
#>   delta pi     1.00 0.00 0
#>   u            0.12 1.00 0
#>   r           -0.21 0.47 1
#> --------------------------------------------------------------------------------
#> 
#> 
#> gamma_0
#> --------------------------------------------------------------------------------
#> delta pi        u        r 
#>    -0.27    -0.35    -0.18 
#> --------------------------------------------------------------------------------
#> 
#> 
#> gamma_1
#> --------------------------------------------------------------------------------
#> delta pi        u        r 
#>     0.90     0.88     0.88 
#> --------------------------------------------------------------------------------
#> 
#> 
#> Phi
#> --------------------------------------------------------------------------------          
#>            delta pi    u    r
#>   delta pi     0.18 0.06 0.08
#>   u            0.06 0.26 0.09
#>   r            0.08 0.09 0.30
#> --------------------------------------------------------------------------------
```

You can always look at the `stanfit` object `bvar_obj$fit$stan` directly
if you want. Note that the `z`’s below are not parameters per se, they
are simply used in a reparameterization trick to sample the log
volatilities more efficiently.

``` r

print(bvar_obj$fit$stan)
#> Inference for Stan model: steady_state_bvar_AR1_stochastic_volatility.
#> 2 chains, each with iter=4000; warmup=1000; thin=1; 
#> post-warmup draws per chain=3000, total post-warmup draws=6000.
#> 
#>                        mean se_mean    sd  2.5%   25%   50%   75%  97.5% n_eff Rhat
#> beta[1,1]              1.26    0.00  0.06  1.15  1.22  1.26  1.30   1.37  6481    1
#> beta[1,2]              0.01    0.00  0.04 -0.07 -0.02  0.01  0.03   0.08  5457    1
#> beta[1,3]              0.15    0.00  0.08  0.00  0.10  0.15  0.21   0.31  6271    1
#> beta[2,1]             -0.09    0.00  0.04 -0.16 -0.12 -0.09 -0.07  -0.02  5999    1
#> beta[2,2]              1.15    0.00  0.06  1.03  1.11  1.15  1.19   1.26  4027    1
#> beta[2,3]             -0.16    0.00  0.08 -0.31 -0.21 -0.16 -0.11  -0.01  4383    1
#> beta[3,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.02   0.04  6232    1
#> beta[3,2]             -0.01    0.00  0.02 -0.05 -0.02 -0.01  0.00   0.02  6536    1
#> beta[3,3]              1.04    0.00  0.06  0.92  1.00  1.04  1.08   1.16  4533    1
#> beta[4,1]             -0.27    0.00  0.06 -0.38 -0.30 -0.27 -0.23  -0.15  6333    1
#> beta[4,2]              0.02    0.00  0.04 -0.05 -0.01  0.02  0.04   0.09  5539    1
#> beta[4,3]             -0.10    0.00  0.08 -0.26 -0.16 -0.10 -0.05   0.06  6350    1
#> beta[5,1]              0.07    0.00  0.03  0.01  0.05  0.07  0.10   0.14  6054    1
#> beta[5,2]             -0.22    0.00  0.05 -0.32 -0.25 -0.22 -0.18  -0.11  4063    1
#> beta[5,3]              0.17    0.00  0.07  0.02  0.12  0.17  0.22   0.31  4532    1
#> beta[6,1]             -0.01    0.00  0.02 -0.04 -0.02 -0.01  0.01   0.03  6596    1
#> beta[6,2]              0.03    0.00  0.02  0.00  0.01  0.03  0.04   0.06  6470    1
#> beta[6,3]             -0.10    0.00  0.06 -0.21 -0.14 -0.10 -0.06   0.02  4588    1
#> Psi[1,1]               1.99    0.00  0.05  1.89  1.96  1.99  2.03   2.09 11742    1
#> Psi[2,1]               4.27    0.00  0.17  3.93  4.16  4.28  4.39   4.61  9945    1
#> Psi[3,1]               3.53    0.00  0.33  2.88  3.30  3.53  3.76   4.17  9764    1
#> z[1,1]                 0.06    0.01  0.50 -0.82 -0.29  0.03  0.37   1.12  8429    1
#> z[1,2]                 1.06    0.01  0.61 -0.15  0.66  1.06  1.45   2.28  7321    1
#> z[1,3]                -0.83    0.01  0.65 -2.05 -1.28 -0.85 -0.42   0.46  7148    1
#> z[2,1]                -0.09    0.01  0.98 -2.00 -0.73 -0.10  0.55   1.85 11629    1
#> z[2,2]                 0.60    0.01  0.99 -1.33 -0.06  0.60  1.28   2.52 10050    1
#> z[2,3]                -0.10    0.01  0.94 -1.97 -0.74 -0.08  0.52   1.71 11341    1
#> z[3,1]                 0.16    0.01  0.96 -1.78 -0.47  0.16  0.82   2.02  8023    1
#> z[3,2]                 0.83    0.01  0.94 -1.07  0.19  0.83  1.48   2.64  9018    1
#> z[3,3]                -0.20    0.01  0.98 -2.11 -0.88 -0.21  0.48   1.72 11030    1
#> z[4,1]                -0.33    0.01  0.94 -2.11 -0.96 -0.34  0.31   1.47 10160    1
#> z[4,2]                -0.18    0.01  0.93 -1.98 -0.79 -0.17  0.44   1.65 11034    1
#> z[4,3]                -0.42    0.01  0.94 -2.28 -1.05 -0.42  0.21   1.46 13689    1
#> z[5,1]                -0.13    0.01  0.93 -1.96 -0.75 -0.15  0.47   1.76 11014    1
#> z[5,2]                -0.29    0.01  0.97 -2.17 -0.96 -0.29  0.35   1.63 13627    1
#> z[5,3]                -0.41    0.01  0.96 -2.29 -1.06 -0.41  0.24   1.46 10234    1
#> z[6,1]                -0.08    0.01  0.96 -1.92 -0.73 -0.08  0.56   1.82 11833    1
#> z[6,2]                -0.07    0.01  0.92 -1.85 -0.71 -0.07  0.55   1.75 10373    1
#> z[6,3]                -0.31    0.01  0.94 -2.12 -0.94 -0.31  0.31   1.55  9506    1
#> z[7,1]                 0.22    0.01  0.93 -1.65 -0.41  0.22  0.85   2.01  9358    1
#> z[7,2]                 0.21    0.01  0.90 -1.54 -0.39  0.20  0.81   2.01  9572    1
#> z[7,3]                -0.25    0.01  0.96 -2.11 -0.90 -0.25  0.41   1.65 11163    1
#> z[8,1]                 0.46    0.01  0.94 -1.43 -0.16  0.46  1.09   2.31  9705    1
#> z[8,2]                 0.00    0.01  0.93 -1.87 -0.61 -0.01  0.60   1.86 11565    1
#> z[8,3]                -0.10    0.01  0.95 -1.98 -0.75 -0.10  0.55   1.73 10173    1
#> z[9,1]                 0.61    0.01  0.92 -1.19 -0.02  0.62  1.24   2.38  9039    1
#> z[9,2]                 0.32    0.01  0.92 -1.47 -0.32  0.32  0.95   2.13 10395    1
#> z[9,3]                 0.07    0.01  0.94 -1.79 -0.54  0.08  0.71   1.92 10930    1
#> z[10,1]               -0.06    0.01  0.90 -1.82 -0.65 -0.07  0.53   1.69  8187    1
#> z[10,2]                0.19    0.01  0.92 -1.60 -0.42  0.18  0.79   2.02 10299    1
#> z[10,3]               -0.19    0.01  0.96 -2.09 -0.84 -0.18  0.46   1.65  9363    1
#> z[11,1]               -0.10    0.01  0.93 -1.91 -0.73 -0.11  0.54   1.73 10126    1
#> z[11,2]                0.03    0.01  0.94 -1.79 -0.61  0.02  0.66   1.89 10667    1
#> z[11,3]               -0.39    0.01  0.97 -2.30 -1.05 -0.40  0.26   1.50  8894    1
#> z[12,1]               -0.37    0.01  0.93 -2.17 -0.99 -0.38  0.26   1.44 10343    1
#> z[12,2]                0.20    0.01  0.95 -1.68 -0.43  0.20  0.83   2.05 11388    1
#> z[12,3]               -0.33    0.01  0.94 -2.15 -0.97 -0.34  0.32   1.52  9211    1
#> z[13,1]               -0.11    0.01  0.94 -1.94 -0.73 -0.11  0.51   1.75  9802    1
#> z[13,2]                0.54    0.01  0.91 -1.21 -0.08  0.54  1.16   2.31 11037    1
#> z[13,3]               -0.15    0.01  0.93 -1.97 -0.80 -0.15  0.48   1.69  9315    1
#> z[14,1]               -0.28    0.01  0.93 -2.12 -0.92 -0.28  0.35   1.53 10023    1
#> z[14,2]                0.58    0.01  0.90 -1.15 -0.04  0.59  1.21   2.32 10239    1
#> z[14,3]               -0.25    0.01  0.96 -2.13 -0.92 -0.26  0.41   1.66 10871    1
#> z[15,1]               -0.22    0.01  0.92 -2.03 -0.85 -0.22  0.39   1.57  8420    1
#> z[15,2]                0.33    0.01  0.92 -1.45 -0.30  0.34  0.93   2.13  9760    1
#> z[15,3]               -0.17    0.01  0.95 -2.05 -0.80 -0.15  0.48   1.68  9915    1
#> z[16,1]                0.00    0.01  0.94 -1.84 -0.65  0.00  0.62   1.84 11105    1
#> z[16,2]                0.50    0.01  0.94 -1.28 -0.16  0.50  1.15   2.31 11034    1
#> z[16,3]                0.02    0.01  0.91 -1.76 -0.60  0.02  0.62   1.82 10304    1
#> z[17,1]               -0.02    0.01  0.92 -1.82 -0.62 -0.01  0.60   1.79 10470    1
#> z[17,2]                0.55    0.01  0.95 -1.31 -0.07  0.55  1.18   2.44 10542    1
#> z[17,3]                0.08    0.01  0.96 -1.80 -0.56  0.07  0.72   1.95  9810    1
#> z[18,1]               -0.12    0.01  0.97 -2.04 -0.77 -0.11  0.52   1.74 10770    1
#> z[18,2]                0.88    0.01  0.93 -0.90  0.24  0.87  1.52   2.71 10108    1
#> z[18,3]                0.26    0.01  0.94 -1.56 -0.37  0.27  0.90   2.15  8966    1
#> z[19,1]                0.15    0.01  0.97 -1.74 -0.51  0.16  0.81   2.04  9012    1
#> z[19,2]                1.21    0.01  0.92 -0.60  0.58  1.22  1.83   2.99  8349    1
#> z[19,3]                0.22    0.01  0.95 -1.68 -0.40  0.23  0.85   2.10 10613    1
#> z[20,1]               -0.06    0.01  0.93 -1.87 -0.68 -0.04  0.57   1.74 10416    1
#> z[20,2]                0.69    0.01  0.87 -1.00  0.09  0.69  1.28   2.37 10711    1
#> z[20,3]                0.20    0.01  0.95 -1.64 -0.45  0.20  0.84   2.05  9298    1
#> z[21,1]               -0.04    0.01  0.93 -1.86 -0.66 -0.04  0.57   1.76 11663    1
#> z[21,2]               -0.02    0.01  0.93 -1.86 -0.66 -0.03  0.61   1.81 10112    1
#> z[21,3]                0.29    0.01  0.92 -1.53 -0.34  0.29  0.92   2.11  9245    1
#> z[22,1]                0.12    0.01  0.94 -1.69 -0.51  0.12  0.75   1.94 10845    1
#> z[22,2]                0.31    0.01  0.90 -1.43 -0.32  0.33  0.92   2.07  9477    1
#> z[22,3]                0.54    0.01  0.90 -1.21 -0.06  0.53  1.14   2.31 10347    1
#> z[23,1]               -0.13    0.01  0.94 -2.02 -0.76 -0.13  0.50   1.71  9489    1
#> z[23,2]                0.05    0.01  0.93 -1.74 -0.60  0.05  0.68   1.86  8961    1
#> z[23,3]               -0.24    0.01  0.97 -2.10 -0.90 -0.25  0.41   1.68 13171    1
#> z[24,1]               -0.27    0.01  0.94 -2.11 -0.92 -0.28  0.37   1.56 13025    1
#> z[24,2]                0.25    0.01  0.92 -1.56 -0.38  0.23  0.87   2.04  8839    1
#> z[24,3]               -0.04    0.01  0.94 -1.88 -0.69 -0.06  0.61   1.78  9747    1
#> z[25,1]               -0.18    0.01  0.94 -2.01 -0.81 -0.19  0.46   1.68 11789    1
#> z[25,2]                0.26    0.01  0.92 -1.54 -0.37  0.26  0.88   2.08 11122    1
#> z[25,3]                0.17    0.01  0.94 -1.67 -0.45  0.19  0.80   2.03 11060    1
#> z[26,1]                0.17    0.01  0.94 -1.65 -0.47  0.17  0.79   2.02  9890    1
#> z[26,2]                0.60    0.01  0.90 -1.15 -0.01  0.59  1.20   2.37  9789    1
#> z[26,3]                0.41    0.01  0.93 -1.39 -0.22  0.40  1.03   2.23 10454    1
#> z[27,1]               -0.12    0.01  0.95 -1.98 -0.75 -0.11  0.50   1.74 12057    1
#>  [ reached 'max' / getOption("max.print") -- omitted 3783 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Mon Jul 27 06:18:40 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

We can forecast

``` r

forecast(bvar_obj, pi = 0.68, show_all = TRUE)
```

![](figure/AR(1)-2-1.png)![](figure/AR(1)-2-2.png)![](figure/AR(1)-2-3.png)

Let us plot the log volatility estimates and predictions

``` r

stochastic_volatility_plot(bvar_obj, ci = 0.95, vol = "log_lambda")
```

![](figure/AR(1)-3-1.png)![](figure/AR(1)-3-2.png)![](figure/AR(1)-3-3.png)

Let us plot the estimates and predictions of the implied innovation
standard deviations

``` r

stochastic_volatility_plot(bvar_obj, vol = "sd")
```

![](figure/AR(1)-4-1.png)![](figure/AR(1)-4-2.png)![](figure/AR(1)-4-3.png)

We can also produce orthogonalized IRFs

``` r

IRF(bvar_obj, method = "OIRF", t=215, ci=0.68) #latest t
```

![](figure/AR(1)-5%20-1.png)

## References

Carriero, A., Clark, T. E., and Marcellino, M. (2024). Capturing
macro-economic tail risks with Bayesian vector autoregressions. *Journal
of Money, Credit and Banking*, 56(5), pp. 1099–1127.

Koop, G. and Korobilis, D. (2010). Bayesian multivariate time series
methods for empirical macroeconomics. *Foundations and Trends in
Econometrics*, 3(4), pp. 267–358.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630–650.
