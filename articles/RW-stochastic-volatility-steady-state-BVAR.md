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

![](figure/RW-1-1.png)

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
for more information about the prior specification. I take some
inspiration from Clark (2011) below.

``` r

k <- bvar_obj$setup$k
n_free_params_A <- bvar_obj$setup$n_free_params_A
sigma2 <- diag(bvar_obj$setup$Sigma_AR)

SV_priors_RW <- list(
                     theta_A             =  rep(0, n_free_params_A),
                     Omega_A             =  diag(10, n_free_params_A),
                     mu_log_lambda_1     =  log(sigma2),
                     sigma2_log_lambda_1 =  rep(4, k),
                     alpha_phi           =  rep(5, k),
                     beta_phi            = (rep(5, k) - 1) * rep(0.1, k)
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
                   SV_type = "RW",
                   SV_priors = SV_priors_RW)
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
                control = list(max_treedepth = 14, adapt_delta = 0.95))
```

Now let’s see the posterior means

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
#>   u.l1           -0.09  1.17 -0.16
#>   r.l1            0.00 -0.01  1.04
#>   delta pi.l2    -0.27  0.02 -0.10
#>   u.l2            0.07 -0.23  0.17
#>   r.l2            0.00  0.02 -0.11
#> --------------------------------------------------------------------------------
#> 
#> 
#> Psi
#> --------------------------------------------------------------------------------          
#>            [,1]
#>   delta pi 2.00
#>   u        4.28
#>   r        3.50
#> --------------------------------------------------------------------------------
#> 
#> 
#> Sigma_u,t (t = 215)
#> --------------------------------------------------------------------------------
#>          delta pi     u     r
#> delta pi     0.09 -0.01  0.03
#> u           -0.01  0.02 -0.01
#> r            0.03 -0.01  0.16
#> --------------------------------------------------------------------------------
#> 
#> 
#> A
#> --------------------------------------------------------------------------------          
#>            delta pi   u r
#>   delta pi     1.00 0.0 0
#>   u            0.12 1.0 0
#>   r           -0.22 0.5 1
#> --------------------------------------------------------------------------------
#> 
#> 
#> phi
#> --------------------------------------------------------------------------------
#> delta pi        u        r 
#>     0.06     0.09     0.11 
#> --------------------------------------------------------------------------------
```

You can always look at the `stanfit` object `bvar_obj$fit$stan` directly
if you want. Note that the `z`’s below are not parameters per se, they
are simply used in a reparameterization trick to sample the log
volatilities more efficiently.

``` r

print(bvar_obj$fit$stan)
#> Inference for Stan model: steady_state_bvar_RW_stochastic_volatility.
#> 2 chains, each with iter=4000; warmup=1000; thin=1; 
#> post-warmup draws per chain=3000, total post-warmup draws=6000.
#> 
#>                        mean se_mean    sd  2.5%   25%   50%   75%  97.5% n_eff Rhat
#> beta[1,1]              1.26    0.00  0.06  1.15  1.22  1.26  1.30   1.37  7006    1
#> beta[1,2]              0.01    0.00  0.04 -0.06 -0.01  0.01  0.04   0.09  7865    1
#> beta[1,3]              0.15    0.00  0.08  0.00  0.09  0.15  0.20   0.31  6222    1
#> beta[2,1]             -0.09    0.00  0.03 -0.16 -0.12 -0.09 -0.07  -0.03  7176    1
#> beta[2,2]              1.17    0.00  0.06  1.06  1.13  1.17  1.21   1.28  7552    1
#> beta[2,3]             -0.16    0.00  0.08 -0.32 -0.21 -0.16 -0.10  -0.01  6771    1
#> beta[3,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.01   0.04  8959    1
#> beta[3,2]             -0.01    0.00  0.02 -0.04 -0.02 -0.01  0.00   0.02  9530    1
#> beta[3,3]              1.04    0.00  0.06  0.92  1.00  1.04  1.08   1.16  7421    1
#> beta[4,1]             -0.27    0.00  0.06 -0.38 -0.31 -0.28 -0.23  -0.16  7153    1
#> beta[4,2]              0.02    0.00  0.04 -0.06 -0.01  0.02  0.04   0.09  8015    1
#> beta[4,3]             -0.10    0.00  0.08 -0.26 -0.16 -0.10 -0.05   0.05  6282    1
#> beta[5,1]              0.07    0.00  0.03  0.01  0.05  0.07  0.09   0.13  7279    1
#> beta[5,2]             -0.23    0.00  0.05 -0.33 -0.27 -0.23 -0.19  -0.12  7742    1
#> beta[5,3]              0.17    0.00  0.07  0.03  0.12  0.17  0.22   0.32  6986    1
#> beta[6,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.01   0.03  8316    1
#> beta[6,2]              0.02    0.00  0.02 -0.01  0.01  0.02  0.03   0.05  8810    1
#> beta[6,3]             -0.11    0.00  0.06 -0.22 -0.15 -0.11 -0.07   0.01  7252    1
#> Psi[1,1]               2.00    0.00  0.05  1.89  1.96  2.00  2.03   2.10 16547    1
#> Psi[2,1]               4.28    0.00  0.17  3.94  4.16  4.28  4.40   4.60 13232    1
#> Psi[3,1]               3.50    0.00  0.32  2.87  3.29  3.51  3.72   4.12 12420    1
#> z[1,1]                -0.03    0.00  0.27 -0.52 -0.22 -0.04  0.14   0.51  9674    1
#> z[1,2]                 0.65    0.00  0.29  0.08  0.45  0.64  0.84   1.25  9043    1
#> z[1,3]                -0.66    0.00  0.34 -1.30 -0.90 -0.67 -0.44   0.04  9462    1
#> z[2,1]                -0.13    0.01  0.96 -1.99 -0.79 -0.14  0.50   1.75 15103    1
#> z[2,2]                 0.24    0.01  0.99 -1.69 -0.43  0.23  0.90   2.20 11692    1
#> z[2,3]                 0.00    0.01  0.97 -1.92 -0.65  0.00  0.67   1.95 13873    1
#> z[3,1]                -0.05    0.01  0.97 -1.91 -0.71 -0.05  0.62   1.83 13394    1
#> z[3,2]                 0.31    0.01  0.97 -1.61 -0.33  0.32  0.96   2.23 10142    1
#> z[3,3]                -0.05    0.01  0.96 -1.93 -0.72 -0.05  0.60   1.82 14505    1
#> z[4,1]                -0.09    0.01  0.97 -1.98 -0.74 -0.10  0.55   1.82 14529    1
#> z[4,2]                -0.47    0.01  0.96 -2.33 -1.12 -0.48  0.18   1.44 13354    1
#> z[4,3]                -0.24    0.01  0.93 -2.05 -0.86 -0.23  0.40   1.55 15334    1
#> z[5,1]                 0.01    0.01  1.00 -1.96 -0.65  0.02  0.68   1.99 17285    1
#> z[5,2]                -0.48    0.01  0.97 -2.43 -1.14 -0.46  0.18   1.45 15175    1
#> z[5,3]                -0.20    0.01  0.94 -2.09 -0.82 -0.19  0.42   1.67 14648    1
#> z[6,1]                 0.01    0.01  0.97 -1.86 -0.64  0.01  0.65   1.92 17221    1
#> z[6,2]                -0.35    0.01  0.96 -2.22 -1.00 -0.35  0.31   1.51 13492    1
#> z[6,3]                -0.13    0.01  0.94 -1.99 -0.76 -0.14  0.50   1.70 13985    1
#> z[7,1]                 0.14    0.01  0.98 -1.81 -0.51  0.14  0.79   2.05 14866    1
#> z[7,2]                -0.22    0.01  0.96 -2.08 -0.88 -0.21  0.43   1.67 13625    1
#> z[7,3]                -0.04    0.01  0.95 -1.93 -0.68 -0.03  0.62   1.80 14365    1
#> z[8,1]                 0.23    0.01  0.98 -1.73 -0.44  0.23  0.90   2.13 15389    1
#> z[8,2]                -0.25    0.01  0.94 -2.11 -0.87 -0.23  0.35   1.61 14216    1
#> z[8,3]                 0.07    0.01  0.94 -1.77 -0.57  0.08  0.71   1.86 14645    1
#> z[9,1]                 0.23    0.01  0.97 -1.67 -0.43  0.25  0.86   2.14 16471    1
#> z[9,2]                -0.12    0.01  0.95 -1.96 -0.76 -0.11  0.51   1.73 15187    1
#> z[9,3]                 0.22    0.01  0.94 -1.62 -0.42  0.22  0.83   2.04 13872    1
#> z[10,1]               -0.17    0.01  0.95 -2.02 -0.80 -0.17  0.44   1.71 16319    1
#> z[10,2]               -0.11    0.01  0.95 -1.99 -0.75 -0.11  0.53   1.76 15556    1
#> z[10,3]               -0.06    0.01  0.97 -1.98 -0.71 -0.10  0.60   1.87 15242    1
#> z[11,1]               -0.15    0.01  0.95 -1.99 -0.81 -0.15  0.49   1.73 16400    1
#> z[11,2]               -0.13    0.01  0.96 -2.01 -0.77 -0.13  0.53   1.74 15579    1
#> z[11,3]               -0.16    0.01  0.97 -2.05 -0.83 -0.17  0.49   1.77 14933    1
#> z[12,1]               -0.35    0.01  0.96 -2.28 -0.99 -0.36  0.30   1.50 15723    1
#> z[12,2]               -0.02    0.01  0.96 -1.89 -0.65 -0.02  0.61   1.88 15720    1
#> z[12,3]               -0.10    0.01  0.97 -2.04 -0.74 -0.10  0.55   1.79 14175    1
#> z[13,1]               -0.26    0.01  0.96 -2.11 -0.91 -0.26  0.38   1.62 17321    1
#> z[13,2]                0.12    0.01  0.95 -1.83 -0.51  0.12  0.74   1.98 15008    1
#> z[13,3]                0.04    0.01  0.95 -1.80 -0.59  0.03  0.69   1.91 15647    1
#> z[14,1]               -0.32    0.01  0.96 -2.22 -0.98 -0.33  0.32   1.55 15719    1
#> z[14,2]                0.16    0.01  0.95 -1.67 -0.50  0.16  0.82   2.01 15518    1
#> z[14,3]                0.03    0.01  0.97 -1.86 -0.61  0.04  0.67   1.93 16327    1
#> z[15,1]               -0.26    0.01  0.96 -2.12 -0.93 -0.26  0.39   1.62 15220    1
#> z[15,2]                0.07    0.01  0.96 -1.83 -0.58  0.08  0.73   1.96 15229    1
#> z[15,3]                0.11    0.01  0.95 -1.72 -0.54  0.11  0.78   1.96 15623    1
#> z[16,1]               -0.17    0.01  0.95 -2.03 -0.83 -0.17  0.47   1.67 14146    1
#> z[16,2]                0.15    0.01  0.98 -1.80 -0.52  0.15  0.82   2.06 15368    1
#> z[16,3]                0.26    0.01  0.94 -1.55 -0.38  0.26  0.89   2.08 14018    1
#> z[17,1]               -0.19    0.01  0.95 -2.04 -0.83 -0.18  0.46   1.69 17853    1
#> z[17,2]                0.19    0.01  0.93 -1.61 -0.45  0.20  0.83   1.98 14709    1
#> z[17,3]                0.31    0.01  0.94 -1.54 -0.32  0.31  0.95   2.18 14415    1
#> z[18,1]               -0.30    0.01  0.96 -2.14 -0.96 -0.30  0.34   1.59 14903    1
#> z[18,2]                0.33    0.01  0.96 -1.55 -0.32  0.32  0.96   2.24 16566    1
#> z[18,3]                0.43    0.01  0.96 -1.44 -0.23  0.43  1.08   2.28 13465    1
#> z[19,1]               -0.16    0.01  0.96 -2.08 -0.78 -0.15  0.48   1.73 17488    1
#> z[19,2]                0.44    0.01  0.94 -1.41 -0.19  0.47  1.07   2.24 16186    1
#> z[19,3]                0.39    0.01  0.93 -1.42 -0.22  0.40  1.02   2.20 15204    1
#> z[20,1]               -0.17    0.01  0.97 -2.10 -0.81 -0.18  0.47   1.78 14318    1
#> z[20,2]                0.07    0.01  0.92 -1.75 -0.55  0.06  0.69   1.86 15464    1
#> z[20,3]                0.33    0.01  0.92 -1.46 -0.31  0.33  0.95   2.15 15141    1
#> z[21,1]               -0.06    0.01  0.97 -2.02 -0.71 -0.06  0.60   1.83 13342    1
#> z[21,2]               -0.44    0.01  0.96 -2.32 -1.10 -0.44  0.21   1.41 15241    1
#> z[21,3]                0.37    0.01  0.96 -1.56 -0.27  0.36  1.01   2.27 13459    1
#> z[22,1]               -0.07    0.01  0.99 -2.00 -0.71 -0.07  0.60   1.87 15725    1
#> z[22,2]               -0.31    0.01  0.96 -2.18 -0.95 -0.30  0.35   1.57 13059    1
#> z[22,3]                0.52    0.01  0.94 -1.32 -0.10  0.53  1.17   2.36 16358    1
#> z[23,1]               -0.02    0.01  0.97 -1.95 -0.67 -0.01  0.63   1.86 12680    1
#> z[23,2]               -0.36    0.01  0.98 -2.27 -1.02 -0.36  0.31   1.51 16821    1
#> z[23,3]               -0.18    0.01  0.94 -2.04 -0.81 -0.16  0.46   1.62 16186    1
#> z[24,1]               -0.20    0.01  0.96 -2.10 -0.83 -0.19  0.44   1.69 17884    1
#> z[24,2]               -0.25    0.01  0.94 -2.10 -0.88 -0.26  0.36   1.58 14659    1
#> z[24,3]               -0.07    0.01  0.95 -1.91 -0.72 -0.07  0.59   1.78 14450    1
#> z[25,1]               -0.19    0.01  0.97 -2.06 -0.85 -0.19  0.49   1.67 14051    1
#> z[25,2]               -0.25    0.01  0.92 -2.06 -0.88 -0.25  0.37   1.57 15451    1
#> z[25,3]                0.07    0.01  0.92 -1.73 -0.56  0.07  0.69   1.87 11939    1
#> z[26,1]               -0.04    0.01  0.98 -1.99 -0.68 -0.04  0.60   1.91 15991    1
#> z[26,2]               -0.14    0.01  0.94 -2.02 -0.76 -0.12  0.51   1.68 17109    1
#> z[26,3]                0.20    0.01  0.94 -1.61 -0.44  0.20  0.85   2.02 14011    1
#> z[27,1]               -0.08    0.01  0.97 -1.99 -0.73 -0.08  0.57   1.76 16457    1
#>  [ reached 'max' / getOption("max.print") -- omitted 3759 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Thu Jul 16 04:38:31 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

We can forecast

``` r

forecast(bvar_obj, pi = 0.68, show_all = TRUE)
```

![](figure/RW-2-1.png)![](figure/RW-2-2.png)![](figure/RW-2-3.png)

Let us plot the log volatility estimates and predictions

``` r

stochastic_volatility_plot(bvar_obj, ci = 0.95, vol = "log_lambda")
```

![](figure/RW-3-1.png)![](figure/RW-3-2.png)![](figure/RW-3-3.png)

Let us plot the estimates and predictions of the implied innovation
standard deviations

``` r

stochastic_volatility_plot(bvar_obj, vol = "sd")
```

![](figure/RW-4-1.png)![](figure/RW-4-2.png)![](figure/RW-4-3.png)

We can also produce orthogonalized IRFs

``` r

IRF(bvar_obj, method = "OIRF", t=215, ci=0.68) #latest t
```

![](figure/RW-5-1.png)

## References

Clark, T. E. (2011). Real-time density forecasts from Bayesian vector
autoregressions with stochastic volatility. *Journal of Business &
Economic Statistics*, 29(3), pp. 327-341.

Koop, G. and Korobilis, D. (2010). Bayesian multivariate time series
methods for empirical macroeconomics. *Foundations and Trends in
Econometrics*, 3(4), pp. 267-358.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630-650.
