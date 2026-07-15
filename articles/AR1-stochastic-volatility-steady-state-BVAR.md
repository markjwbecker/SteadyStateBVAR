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
for more information about the prior specification.

``` r

k <- bvar_obj$setup$k
n_free_params_A <- bvar_obj$setup$n_free_params_A

SV_priors_AR1 <- list(
                      theta_A            =  rep(0, n_free_params_A),
                      Omega_A            =  diag(10, n_free_params_A),
                      theta_gamma_0      =  rep(0, k),
                      Omega_gamma_0      =  diag(10, k),
                      theta_gamma_1      =  rep(0.9, k),
                      Omega_gamma_1      =  diag(1, k),
                      theta_log_lambda_1 =  rep(0, k),
                      Omega_log_lambda_1 =  diag(10, k),
                      V_Phi              = (5 - k - 1) * 0.1 * diag(k),
                      m_Phi              =  5
                     )
```

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
#>   r        3.51
#> --------------------------------------------------------------------------------
#> 
#> 
#> Sigma_u,t (t = 215)
#> --------------------------------------------------------------------------------
#>          delta pi     u     r
#> delta pi     0.07 -0.01  0.02
#> u           -0.01  0.03 -0.02
#> r            0.02 -0.02  0.18
#> --------------------------------------------------------------------------------
#> 
#> 
#> A
#> --------------------------------------------------------------------------------          
#>            delta pi    u r
#>   delta pi     1.00 0.00 0
#>   u            0.13 1.00 0
#>   r           -0.24 0.42 1
#> --------------------------------------------------------------------------------
#> 
#> 
#> gamma_0
#> --------------------------------------------------------------------------------
#> delta pi        u        r 
#>    -0.17    -0.18    -0.12 
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
#>   delta pi     0.08 0.05 0.09
#>   u            0.05 0.10 0.10
#>   r            0.09 0.10 0.20
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
#> beta[1,1]              1.27    0.00  0.06  1.16  1.23  1.27  1.31  1.38  6821    1
#> beta[1,2]              0.02    0.00  0.04 -0.05 -0.01  0.02  0.05  0.10  8640    1
#> beta[1,3]              0.17    0.00  0.08  0.01  0.11  0.17  0.23  0.33  7507    1
#> beta[2,1]             -0.09    0.00  0.04 -0.16 -0.11 -0.09 -0.06 -0.02  7593    1
#> beta[2,2]              1.16    0.00  0.06  1.05  1.13  1.16  1.20  1.27  6105    1
#> beta[2,3]             -0.15    0.00  0.07 -0.30 -0.20 -0.15 -0.10 -0.01  6968    1
#> beta[3,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.02  0.04 10269    1
#> beta[3,2]             -0.01    0.00  0.02 -0.05 -0.02 -0.01  0.00  0.02  7367    1
#> beta[3,3]              1.04    0.00  0.06  0.93  1.00  1.04  1.08  1.16  7309    1
#> beta[4,1]             -0.28    0.00  0.06 -0.39 -0.31 -0.28 -0.24 -0.16  6688    1
#> beta[4,2]              0.01    0.00  0.04 -0.06 -0.01  0.01  0.03  0.08  8586    1
#> beta[4,3]             -0.11    0.00  0.08 -0.27 -0.17 -0.11 -0.05  0.05  7193    1
#> beta[5,1]              0.07    0.00  0.03  0.01  0.05  0.07  0.09  0.13  8472    1
#> beta[5,2]             -0.22    0.00  0.05 -0.33 -0.26 -0.22 -0.19 -0.12  6188    1
#> beta[5,3]              0.17    0.00  0.07  0.03  0.12  0.17  0.22  0.31  7472    1
#> beta[6,1]             -0.01    0.00  0.02 -0.04 -0.02 -0.01  0.00  0.02  9445    1
#> beta[6,2]              0.02    0.00  0.02 -0.01  0.01  0.02  0.04  0.06  8249    1
#> beta[6,3]             -0.11    0.00  0.06 -0.22 -0.14 -0.10 -0.07  0.00  7549    1
#> Psi[1,1]               1.99    0.00  0.05  1.89  1.96  2.00  2.03  2.10 15419    1
#> Psi[2,1]               4.27    0.00  0.18  3.93  4.15  4.27  4.39  4.61 11457    1
#> Psi[3,1]               3.51    0.00  0.33  2.85  3.30  3.51  3.73  4.16 12022    1
#> z[1,1]                -0.66    0.00  0.22 -1.05 -0.81 -0.68 -0.52 -0.17  8458    1
#> z[1,2]                 0.05    0.00  0.25 -0.41 -0.12  0.03  0.20  0.58  6509    1
#> z[1,3]                -0.65    0.00  0.29 -1.19 -0.85 -0.66 -0.47 -0.05  8196    1
#> z[2,1]                -0.04    0.01  0.95 -1.91 -0.69 -0.05  0.60  1.82 11924    1
#> z[2,2]                 0.15    0.01  0.99 -1.82 -0.49  0.14  0.80  2.08 11096    1
#> z[2,3]                 0.04    0.01  0.99 -1.89 -0.63  0.05  0.71  1.99 15399    1
#> z[3,1]                 0.03    0.01  0.95 -1.87 -0.61  0.01  0.68  1.89 12389    1
#> z[3,2]                 0.15    0.01  1.00 -1.76 -0.53  0.13  0.82  2.10 12127    1
#> z[3,3]                -0.07    0.01  0.98 -2.01 -0.71 -0.09  0.59  1.84 16432    1
#> z[4,1]                -0.46    0.01  0.95 -2.33 -1.09 -0.48  0.18  1.39 13983    1
#> z[4,2]                -0.27    0.01  0.97 -2.14 -0.91 -0.27  0.39  1.60 12799    1
#> z[4,3]                -0.20    0.01  0.98 -2.12 -0.85 -0.20  0.46  1.70 15001    1
#> z[5,1]                -0.30    0.01  0.95 -2.17 -0.92 -0.29  0.34  1.53 13136    1
#> z[5,2]                -0.27    0.01  0.97 -2.17 -0.93 -0.27  0.39  1.61 14480    1
#> z[5,3]                -0.18    0.01  0.99 -2.11 -0.86 -0.19  0.49  1.78 16042    1
#> z[6,1]                -0.18    0.01  0.93 -2.00 -0.82 -0.17  0.46  1.61 13868    1
#> z[6,2]                -0.13    0.01  0.95 -2.02 -0.75 -0.13  0.49  1.72 14143    1
#> z[6,3]                -0.15    0.01  0.98 -2.06 -0.83 -0.15  0.50  1.77 14733    1
#> z[7,1]                 0.08    0.01  0.96 -1.77 -0.56  0.07  0.73  2.01 13585    1
#> z[7,2]                 0.01    0.01  0.96 -1.89 -0.64  0.02  0.65  1.88 14985    1
#> z[7,3]                -0.13    0.01  0.97 -2.02 -0.80 -0.12  0.53  1.76 14804    1
#> z[8,1]                 0.30    0.01  0.93 -1.54 -0.32  0.30  0.91  2.15 12010    1
#> z[8,2]                 0.03    0.01  0.95 -1.80 -0.62  0.03  0.67  1.87 12452    1
#> z[8,3]                -0.06    0.01  0.99 -1.97 -0.74 -0.07  0.62  1.87 13349    1
#> z[9,1]                 0.52    0.01  0.92 -1.30 -0.10  0.53  1.14  2.26 12144    1
#> z[9,2]                 0.22    0.01  0.95 -1.64 -0.44  0.22  0.88  2.06 13431    1
#> z[9,3]                 0.04    0.01  0.97 -1.83 -0.62  0.05  0.71  1.95 14869    1
#> z[10,1]               -0.14    0.01  0.91 -1.88 -0.77 -0.13  0.48  1.58 12669    1
#> z[10,2]                0.13    0.01  0.92 -1.67 -0.49  0.14  0.76  1.90 11844    1
#> z[10,3]               -0.13    0.01  0.94 -1.99 -0.75 -0.13  0.50  1.71 11799    1
#> z[11,1]               -0.23    0.01  0.94 -2.10 -0.88 -0.23  0.41  1.58 13528    1
#> z[11,2]                0.08    0.01  0.98 -1.86 -0.59  0.08  0.74  1.99 12261    1
#> z[11,3]               -0.21    0.01  0.98 -2.12 -0.88 -0.22  0.45  1.69 16178    1
#> z[12,1]               -0.36    0.01  0.95 -2.22 -0.99 -0.37  0.28  1.49 13863    1
#> z[12,2]                0.20    0.01  0.96 -1.70 -0.42  0.20  0.84  2.04 12106    1
#> z[12,3]               -0.18    0.01  0.99 -2.05 -0.86 -0.19  0.50  1.75 14460    1
#> z[13,1]               -0.08    0.01  0.94 -1.95 -0.73 -0.09  0.54  1.79 12148    1
#> z[13,2]                0.40    0.01  0.99 -1.58 -0.25  0.39  1.07  2.35 13709    1
#> z[13,3]               -0.10    0.01  0.99 -2.00 -0.76 -0.11  0.57  1.84 13725    1
#> z[14,1]               -0.21    0.01  0.94 -2.08 -0.83 -0.22  0.39  1.65 14866    1
#> z[14,2]                0.42    0.01  0.95 -1.44 -0.23  0.43  1.07  2.25 11747    1
#> z[14,3]               -0.16    0.01  0.98 -2.11 -0.81 -0.16  0.49  1.77 13370    1
#> z[15,1]               -0.17    0.01  0.96 -2.05 -0.84 -0.16  0.48  1.71 15110    1
#> z[15,2]                0.34    0.01  0.99 -1.62 -0.31  0.34  0.99  2.28 13797    1
#> z[15,3]               -0.12    0.01  0.98 -2.03 -0.76 -0.12  0.53  1.84 14761    1
#> z[16,1]                0.07    0.01  0.93 -1.77 -0.53  0.08  0.68  1.94 14267    1
#> z[16,2]                0.47    0.01  0.97 -1.46 -0.17  0.48  1.11  2.35 12878    1
#> z[16,3]               -0.02    0.01  0.96 -1.89 -0.67 -0.02  0.63  1.85 15656    1
#> z[17,1]                0.11    0.01  0.95 -1.76 -0.54  0.12  0.74  1.96 12981    1
#> z[17,2]                0.54    0.01  0.98 -1.38 -0.12  0.54  1.20  2.47 14072    1
#> z[17,3]                0.04    0.01  0.97 -1.83 -0.62  0.03  0.70  1.93 12043    1
#> z[18,1]                0.21    0.01  0.96 -1.67 -0.43  0.21  0.85  2.13 13137    1
#> z[18,2]                0.75    0.01  0.96 -1.13  0.11  0.75  1.41  2.61 13769    1
#> z[18,3]                0.13    0.01  0.98 -1.82 -0.53  0.13  0.79  2.04 16029    1
#> z[19,1]                0.45    0.01  0.95 -1.45 -0.20  0.45  1.10  2.26 10038    1
#> z[19,2]                0.91    0.01  0.98 -0.99  0.22  0.91  1.57  2.80 10703    1
#> z[19,3]                0.13    0.01  0.98 -1.80 -0.53  0.13  0.79  2.05 13317    1
#> z[20,1]                0.16    0.01  0.96 -1.76 -0.48  0.17  0.80  2.02 13588    1
#> z[20,2]                0.54    0.01  0.93 -1.30 -0.10  0.56  1.16  2.36 13442    1
#> z[20,3]                0.10    0.01  0.95 -1.75 -0.54  0.10  0.75  1.96 13365    1
#> z[21,1]               -0.06    0.01  0.95 -1.90 -0.72 -0.06  0.58  1.81 13362    1
#> z[21,2]                0.10    0.01  0.97 -1.83 -0.54  0.10  0.76  1.99 13229    1
#> z[21,3]                0.13    0.01  0.94 -1.70 -0.51  0.13  0.75  1.98 13589    1
#> z[22,1]                0.19    0.01  0.93 -1.65 -0.41  0.19  0.83  1.95 11502    1
#> z[22,2]                0.31    0.01  0.96 -1.59 -0.33  0.32  0.95  2.19 12386    1
#> z[22,3]                0.26    0.01  0.98 -1.67 -0.40  0.27  0.92  2.16 11749    1
#> z[23,1]               -0.29    0.01  0.93 -2.16 -0.91 -0.28  0.34  1.57 12764    1
#> z[23,2]                0.00    0.01  0.97 -1.87 -0.68  0.00  0.65  1.94 14757    1
#> z[23,3]               -0.14    0.01  0.98 -2.08 -0.80 -0.13  0.51  1.75 15671    1
#> z[24,1]               -0.25    0.01  0.91 -2.00 -0.89 -0.25  0.37  1.55 15673    1
#> z[24,2]                0.14    0.01  0.97 -1.70 -0.51  0.13  0.79  2.04 14732    1
#> z[24,3]               -0.07    0.01  0.97 -1.97 -0.74 -0.06  0.59  1.84 12854    1
#> z[25,1]               -0.07    0.01  0.95 -1.94 -0.72 -0.07  0.58  1.78 14502    1
#> z[25,2]                0.18    0.01  0.97 -1.71 -0.49  0.18  0.83  2.07 15028    1
#> z[25,3]                0.05    0.01  0.98 -1.86 -0.62  0.05  0.72  1.99 14717    1
#> z[26,1]                0.29    0.01  0.95 -1.57 -0.37  0.29  0.93  2.12 12835    1
#> z[26,2]                0.40    0.01  0.94 -1.44 -0.24  0.42  1.05  2.23 11856    1
#> z[26,3]                0.16    0.01  0.98 -1.72 -0.49  0.17  0.81  2.09 13304    1
#> z[27,1]               -0.13    0.01  0.91 -1.92 -0.75 -0.14  0.48  1.63 13547    1
#>  [ reached 'max' / getOption("max.print") -- omitted 3783 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Wed Jul 15 11:37:52 2026.
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

Koop, G. and Korobilis, D. (2010). Bayesian multivariate time series
methods for empirical macroeconomics. *Foundations and Trends in
Econometrics*, 3(4), pp. 267-358.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630-650.
