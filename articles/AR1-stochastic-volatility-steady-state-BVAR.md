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
exact same AR(1) stochastic volatility specification, but for an
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
                      V_Phi              = (5 - k - 1) * 0.1 * diag(k),
                      m_Phi              =  5
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
#>   delta pi.l1     1.27  0.02  0.17
#>   u.l1           -0.09  1.16 -0.15
#>   r.l1            0.00 -0.01  1.04
#>   delta pi.l2    -0.28  0.01 -0.11
#>   u.l2            0.07 -0.22  0.17
#>   r.l2           -0.01  0.02 -0.11
#> --------------------------------------------------------------------------------
#> 
#> 
#> Psi
#> --------------------------------------------------------------------------------          
#>            [,1]
#>   delta pi 1.99
#>   u        4.27
#>   r        3.52
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
#>   r           -0.23 0.44 1
#> --------------------------------------------------------------------------------
#> 
#> 
#> gamma_0
#> --------------------------------------------------------------------------------
#> delta pi        u        r 
#>    -0.16    -0.17    -0.11 
#> --------------------------------------------------------------------------------
#> 
#> 
#> gamma_1
#> --------------------------------------------------------------------------------
#> delta pi        u        r 
#>     0.94     0.94     0.92 
#> --------------------------------------------------------------------------------
#> 
#> 
#> Phi
#> --------------------------------------------------------------------------------          
#>            delta pi    u    r
#>   delta pi     0.08 0.05 0.08
#>   u            0.05 0.10 0.10
#>   r            0.08 0.10 0.20
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
#>                        mean se_mean    sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
#> beta[1,1]              1.27    0.00  0.06  1.16  1.23  1.27  1.31  1.38  4143    1
#> beta[1,2]              0.02    0.00  0.04 -0.06 -0.01  0.02  0.05  0.09  3860    1
#> beta[1,3]              0.17    0.00  0.08  0.01  0.11  0.17  0.22  0.33  4148    1
#> beta[2,1]             -0.09    0.00  0.04 -0.16 -0.11 -0.09 -0.06 -0.02  3526    1
#> beta[2,2]              1.16    0.00  0.06  1.05  1.12  1.16  1.20  1.27  3422    1
#> beta[2,3]             -0.15    0.00  0.08 -0.31 -0.21 -0.15 -0.10  0.00  3855    1
#> beta[3,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.01  0.04  4393    1
#> beta[3,2]             -0.01    0.00  0.02 -0.05 -0.02 -0.01  0.00  0.02  4435    1
#> beta[3,3]              1.04    0.00  0.06  0.93  1.00  1.04  1.09  1.17  3464    1
#> beta[4,1]             -0.28    0.00  0.06 -0.39 -0.32 -0.28 -0.24 -0.16  4046    1
#> beta[4,2]              0.01    0.00  0.04 -0.07 -0.02  0.01  0.04  0.09  3964    1
#> beta[4,3]             -0.11    0.00  0.08 -0.27 -0.17 -0.11 -0.06  0.05  4099    1
#> beta[5,1]              0.07    0.00  0.03  0.00  0.05  0.07  0.09  0.13  3625    1
#> beta[5,2]             -0.22    0.00  0.05 -0.33 -0.26 -0.23 -0.19 -0.12  3463    1
#> beta[5,3]              0.17    0.00  0.07  0.03  0.12  0.17  0.22  0.31  3961    1
#> beta[6,1]             -0.01    0.00  0.02 -0.04 -0.02 -0.01  0.00  0.02  4651    1
#> beta[6,2]              0.02    0.00  0.02 -0.01  0.01  0.02  0.04  0.06  4476    1
#> beta[6,3]             -0.11    0.00  0.06 -0.22 -0.15 -0.11 -0.07  0.00  3464    1
#> Psi[1,1]               1.99    0.00  0.05  1.90  1.96  1.99  2.03  2.10 10586    1
#> Psi[2,1]               4.27    0.00  0.18  3.93  4.16  4.28  4.39  4.61  7486    1
#> Psi[3,1]               3.52    0.00  0.34  2.82  3.30  3.52  3.75  4.18  7764    1
#> z[1,1]                 0.04    0.01  0.43 -0.74 -0.26  0.02  0.31  0.94  4668    1
#> z[1,2]                 1.26    0.01  0.46  0.36  0.95  1.24  1.56  2.17  5532    1
#> z[1,3]                -1.02    0.01  0.56 -2.08 -1.39 -1.03 -0.67  0.14  4554    1
#> z[2,1]                 0.02    0.01  0.96 -1.83 -0.64  0.02  0.66  1.90  7519    1
#> z[2,2]                 0.27    0.01  0.99 -1.64 -0.41  0.27  0.94  2.24  6649    1
#> z[2,3]                -0.08    0.01  0.99 -2.01 -0.74 -0.07  0.61  1.83  9796    1
#> z[3,1]                 0.11    0.01  0.95 -1.75 -0.53  0.11  0.74  1.99  7669    1
#> z[3,2]                 0.31    0.01  0.96 -1.51 -0.35  0.31  0.97  2.20  7788    1
#> z[3,3]                -0.13    0.01  0.96 -2.01 -0.79 -0.14  0.51  1.77  8903    1
#> z[4,1]                -0.44    0.01  0.95 -2.29 -1.08 -0.45  0.17  1.45  7879    1
#> z[4,2]                -0.23    0.01  0.96 -2.12 -0.86 -0.23  0.40  1.67  6911    1
#> z[4,3]                -0.28    0.01  1.00 -2.23 -0.96 -0.28  0.38  1.68  7823    1
#> z[5,1]                -0.29    0.01  0.94 -2.14 -0.94 -0.28  0.35  1.56  7968    1
#> z[5,2]                -0.23    0.01  1.00 -2.22 -0.88 -0.24  0.43  1.73  7802    1
#> z[5,3]                -0.23    0.01  0.99 -2.14 -0.89 -0.25  0.43  1.73  8426    1
#> z[6,1]                -0.17    0.01  0.93 -2.00 -0.79 -0.16  0.44  1.67  9539    1
#> z[6,2]                -0.11    0.01  0.98 -2.02 -0.77 -0.09  0.55  1.78  9178    1
#> z[6,3]                -0.19    0.01  0.98 -2.14 -0.83 -0.19  0.47  1.71 10533    1
#> z[7,1]                 0.09    0.01  0.94 -1.72 -0.56  0.09  0.75  1.94  8671    1
#> z[7,2]                 0.03    0.01  0.98 -1.90 -0.62  0.03  0.67  1.95  9265    1
#> z[7,3]                -0.16    0.01  0.96 -2.10 -0.82 -0.17  0.49  1.69 10401    1
#> z[8,1]                 0.31    0.01  0.95 -1.53 -0.32  0.30  0.97  2.20  8875    1
#> z[8,2]                 0.03    0.01  0.98 -1.89 -0.64  0.05  0.67  1.96  8095    1
#> z[8,3]                -0.10    0.01  0.98 -2.04 -0.77 -0.09  0.56  1.83  9392    1
#> z[9,1]                 0.50    0.01  0.94 -1.30 -0.14  0.51  1.13  2.32  8191    1
#> z[9,2]                 0.21    0.01  0.94 -1.63 -0.42  0.21  0.84  2.06  8865    1
#> z[9,3]                 0.00    0.01  0.97 -1.90 -0.67  0.01  0.65  1.87  9250    1
#> z[10,1]               -0.14    0.01  0.93 -1.97 -0.76 -0.14  0.48  1.75  8134    1
#> z[10,2]                0.14    0.01  0.96 -1.79 -0.52  0.15  0.78  1.99  9674    1
#> z[10,3]               -0.15    0.01  0.98 -2.04 -0.81 -0.15  0.51  1.79 10711    1
#> z[11,1]               -0.22    0.01  0.91 -1.95 -0.84 -0.22  0.40  1.58  8236    1
#> z[11,2]                0.08    0.01  0.98 -1.82 -0.59  0.09  0.75  2.04  8560    1
#> z[11,3]               -0.23    0.01  1.01 -2.19 -0.88 -0.22  0.45  1.73  8555    1
#> z[12,1]               -0.36    0.01  0.91 -2.11 -0.96 -0.37  0.25  1.46  7549    1
#> z[12,2]                0.21    0.01  0.98 -1.72 -0.45  0.21  0.86  2.18 12424    1
#> z[12,3]               -0.21    0.01  0.99 -2.17 -0.88 -0.20  0.46  1.75  9172    1
#> z[13,1]               -0.08    0.01  0.94 -1.93 -0.70 -0.10  0.55  1.76  7763    1
#> z[13,2]                0.41    0.01  0.97 -1.52 -0.24  0.41  1.06  2.32  9488    1
#> z[13,3]               -0.12    0.01  0.98 -2.07 -0.78 -0.11  0.55  1.75  9358    1
#> z[14,1]               -0.20    0.01  0.98 -2.13 -0.88 -0.20  0.48  1.70 11413    1
#> z[14,2]                0.43    0.01  0.97 -1.46 -0.21  0.44  1.09  2.35  8165    1
#> z[14,3]               -0.17    0.01  1.00 -2.14 -0.84 -0.16  0.52  1.74 11178    1
#> z[15,1]               -0.17    0.01  0.96 -2.07 -0.82 -0.17  0.48  1.75  8932    1
#> z[15,2]                0.36    0.01  0.96 -1.52 -0.30  0.35  1.01  2.20  9054    1
#> z[15,3]               -0.12    0.01  1.00 -2.12 -0.80 -0.12  0.54  1.85  8699    1
#> z[16,1]                0.06    0.01  0.94 -1.81 -0.55  0.06  0.69  1.92  8180    1
#> z[16,2]                0.48    0.01  0.99 -1.49 -0.19  0.49  1.15  2.43  8919    1
#> z[16,3]               -0.01    0.01  0.97 -1.86 -0.69 -0.01  0.65  1.86  8815    1
#> z[17,1]                0.11    0.01  0.93 -1.78 -0.50  0.10  0.74  1.91  7399    1
#> z[17,2]                0.55    0.01  0.97 -1.36 -0.10  0.54  1.21  2.47 10740    1
#> z[17,3]                0.02    0.01  0.98 -1.90 -0.65  0.02  0.69  1.96  8155    1
#> z[18,1]                0.21    0.01  0.97 -1.69 -0.45  0.20  0.87  2.12  8813    1
#> z[18,2]                0.73    0.01  0.97 -1.22  0.09  0.74  1.39  2.62  8863    1
#> z[18,3]                0.12    0.01  0.98 -1.84 -0.52  0.12  0.78  2.01  8520    1
#> z[19,1]                0.43    0.01  0.98 -1.49 -0.24  0.43  1.10  2.36  6412    1
#> z[19,2]                0.89    0.01  0.96 -0.99  0.25  0.89  1.54  2.75  7699    1
#> z[19,3]                0.13    0.01  0.97 -1.78 -0.52  0.14  0.77  2.02  8150    1
#> z[20,1]                0.11    0.01  0.93 -1.74 -0.50  0.11  0.73  1.90  8729    1
#> z[20,2]                0.54    0.01  0.95 -1.32 -0.10  0.52  1.19  2.41  9827    1
#> z[20,3]                0.10    0.01  0.98 -1.79 -0.55  0.09  0.76  1.98  9811    1
#> z[21,1]               -0.06    0.01  0.95 -1.92 -0.70 -0.06  0.57  1.79  9958    1
#> z[21,2]                0.09    0.01  0.98 -1.80 -0.56  0.10  0.73  1.99  8200    1
#> z[21,3]                0.16    0.01  0.97 -1.75 -0.49  0.15  0.82  2.04 10677    1
#> z[22,1]                0.19    0.01  0.93 -1.67 -0.43  0.22  0.83  1.97  7851    1
#> z[22,2]                0.27    0.01  0.98 -1.64 -0.40  0.28  0.94  2.20  7553    1
#> z[22,3]                0.28    0.01  0.98 -1.66 -0.37  0.29  0.95  2.19 11927    1
#> z[23,1]               -0.28    0.01  0.92 -2.08 -0.90 -0.28  0.34  1.51 10085    1
#> z[23,2]                0.00    0.01  1.00 -1.96 -0.69 -0.01  0.67  1.96  8562    1
#> z[23,3]               -0.12    0.01  0.98 -2.08 -0.80 -0.12  0.52  1.83  9885    1
#> z[24,1]               -0.22    0.01  0.96 -2.12 -0.87 -0.22  0.44  1.65  8244    1
#> z[24,2]                0.11    0.01  0.95 -1.73 -0.53  0.10  0.75  1.97  7164    1
#> z[24,3]               -0.05    0.01  0.99 -1.99 -0.71 -0.05  0.63  1.89  9442    1
#> z[25,1]               -0.08    0.01  0.94 -1.89 -0.73 -0.07  0.56  1.78  9310    1
#> z[25,2]                0.19    0.01  0.97 -1.73 -0.47  0.19  0.84  2.12  8934    1
#> z[25,3]                0.05    0.01  0.97 -1.82 -0.61  0.05  0.73  1.93 10481    1
#> z[26,1]                0.27    0.01  0.94 -1.57 -0.36  0.27  0.91  2.15  7543    1
#> z[26,2]                0.36    0.01  0.98 -1.58 -0.29  0.36  1.02  2.31  7685    1
#> z[26,3]                0.18    0.01  0.96 -1.71 -0.47  0.19  0.83  2.10  9494    1
#> z[27,1]               -0.15    0.01  0.91 -1.95 -0.76 -0.15  0.46  1.61  9707    1
#>  [ reached 'max' / getOption("max.print") -- omitted 3783 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Thu Jul 16 06:33:12 2026.
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
Econometrics*, 3(4), pp. 267-358.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630-650.
