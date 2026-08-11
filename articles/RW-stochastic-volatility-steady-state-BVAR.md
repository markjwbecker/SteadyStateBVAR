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
                control = list(max_treedepth = 14, adapt_delta = 0.95))
#> ------------------------------------------------------------
#> Estimating Stan model:
#> steady_state_bvar_RW_stochastic_volatility
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
#>   delta pi.l1     1.27  0.02  0.15
#>   u.l1           -0.09  1.17 -0.16
#>   r.l1            0.00 -0.01  1.04
#>   delta pi.l2    -0.28  0.02 -0.11
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
#> delta pi     0.08 -0.01  0.02
#> u           -0.01  0.02 -0.01
#> r            0.02 -0.01  0.16
#> --------------------------------------------------------------------------------
#> 
#> 
#> A
#> --------------------------------------------------------------------------------          
#>            delta pi   u r
#>   delta pi     1.00 0.0 0
#>   u            0.12 1.0 0
#>   r           -0.23 0.5 1
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
#> 2 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=2000.
#> 
#>                        mean se_mean    sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
#> beta[1,1]              1.27    0.00  0.05  1.16  1.23  1.27  1.30  1.37  1715    1
#> beta[1,2]              0.02    0.00  0.04 -0.06 -0.01  0.02  0.04  0.09  2122    1
#> beta[1,3]              0.15    0.00  0.08 -0.01  0.10  0.15  0.20  0.31  1909    1
#> beta[2,1]             -0.09    0.00  0.04 -0.16 -0.12 -0.09 -0.07 -0.02  2189    1
#> beta[2,2]              1.17    0.00  0.06  1.06  1.13  1.17  1.21  1.28  1353    1
#> beta[2,3]             -0.16    0.00  0.08 -0.32 -0.22 -0.16 -0.11 -0.01  1695    1
#> beta[3,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.01  0.04  2534    1
#> beta[3,2]             -0.01    0.00  0.02 -0.04 -0.02 -0.01  0.00  0.02  2131    1
#> beta[3,3]              1.04    0.00  0.06  0.92  1.00  1.04  1.08  1.16  1776    1
#> beta[4,1]             -0.28    0.00  0.05 -0.39 -0.31 -0.28 -0.24 -0.17  1837    1
#> beta[4,2]              0.02    0.00  0.04 -0.06 -0.01  0.02  0.04  0.09  2193    1
#> beta[4,3]             -0.11    0.00  0.08 -0.26 -0.16 -0.11 -0.05  0.05  1888    1
#> beta[5,1]              0.07    0.00  0.03  0.00  0.05  0.07  0.09  0.14  2092    1
#> beta[5,2]             -0.23    0.00  0.05 -0.34 -0.27 -0.23 -0.20 -0.13  1401    1
#> beta[5,3]              0.17    0.00  0.07  0.02  0.12  0.17  0.22  0.32  1813    1
#> beta[6,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.01  0.03  2368    1
#> beta[6,2]              0.02    0.00  0.02 -0.01  0.01  0.02  0.03  0.05  2179    1
#> beta[6,3]             -0.11    0.00  0.06 -0.22 -0.14 -0.11 -0.07  0.01  1828    1
#> Psi[1,1]               2.00    0.00  0.05  1.90  1.97  2.00  2.03  2.09  5477    1
#> Psi[2,1]               4.29    0.00  0.17  3.95  4.16  4.29  4.41  4.62  3703    1
#> Psi[3,1]               3.49    0.01  0.32  2.87  3.27  3.49  3.72  4.12  3679    1
#> z[1,1]                -0.06    0.00  0.25 -0.52 -0.23 -0.07  0.10  0.45  2569    1
#> z[1,2]                 0.64    0.01  0.27  0.13  0.47  0.64  0.81  1.21  2502    1
#> z[1,3]                -0.66    0.01  0.34 -1.32 -0.88 -0.66 -0.44  0.04  2913    1
#> z[2,1]                -0.11    0.02  0.97 -2.09 -0.78 -0.12  0.56  1.71  3524    1
#> z[2,2]                 0.20    0.02  1.01 -1.72 -0.49  0.20  0.86  2.19  3999    1
#> z[2,3]                 0.02    0.02  1.02 -1.91 -0.65  0.01  0.70  2.07  2926    1
#> z[3,1]                -0.05    0.01  1.00 -2.02 -0.75 -0.01  0.63  1.87  5067    1
#> z[3,2]                 0.25    0.02  0.98 -1.66 -0.41  0.26  0.91  2.20  3608    1
#> z[3,3]                -0.02    0.02  0.98 -1.92 -0.70 -0.04  0.61  1.90  3297    1
#> z[4,1]                -0.08    0.02  0.97 -2.02 -0.71 -0.08  0.56  1.82  4105    1
#> z[4,2]                -0.43    0.02  0.97 -2.33 -1.06 -0.44  0.20  1.52  2510    1
#> z[4,3]                -0.27    0.02  0.97 -2.21 -0.93 -0.27  0.40  1.68  3868    1
#> z[5,1]                 0.00    0.02  0.96 -1.85 -0.64  0.00  0.65  1.83  3611    1
#> z[5,2]                -0.49    0.01  0.98 -2.37 -1.14 -0.50  0.14  1.47  4560    1
#> z[5,3]                -0.20    0.01  0.94 -1.95 -0.85 -0.21  0.42  1.63  5615    1
#> z[6,1]                 0.00    0.01  0.97 -1.87 -0.68 -0.01  0.68  1.90  4522    1
#> z[6,2]                -0.33    0.02  0.93 -2.20 -0.94 -0.30  0.29  1.56  3409    1
#> z[6,3]                -0.12    0.02  0.98 -2.01 -0.78 -0.10  0.50  1.74  3693    1
#> z[7,1]                 0.09    0.01  0.99 -1.88 -0.55  0.09  0.75  2.03  4486    1
#> z[7,2]                -0.20    0.02  1.00 -2.22 -0.83 -0.19  0.47  1.81  3697    1
#> z[7,3]                -0.03    0.01  0.93 -1.84 -0.64 -0.03  0.62  1.77  3976    1
#> z[8,1]                 0.19    0.01  0.99 -1.73 -0.49  0.20  0.85  2.11  4332    1
#> z[8,2]                -0.25    0.02  0.96 -2.18 -0.87 -0.27  0.37  1.66  3486    1
#> z[8,3]                 0.08    0.02  0.93 -1.72 -0.54  0.09  0.70  1.99  3427    1
#> z[9,1]                 0.15    0.02  0.97 -1.82 -0.48  0.14  0.82  2.09  4006    1
#> z[9,2]                -0.16    0.02  1.00 -2.18 -0.83 -0.16  0.50  1.83  3961    1
#> z[9,3]                 0.22    0.01  0.94 -1.67 -0.45  0.20  0.87  2.09  4302    1
#> z[10,1]               -0.16    0.01  0.97 -2.11 -0.82 -0.15  0.48  1.74  4276    1
#> z[10,2]               -0.09    0.02  0.96 -1.94 -0.70 -0.10  0.55  1.77  3422    1
#> z[10,3]               -0.09    0.01  0.94 -1.92 -0.70 -0.11  0.55  1.76  4366    1
#> z[11,1]               -0.20    0.02  1.01 -2.24 -0.90 -0.20  0.48  1.81  3588    1
#> z[11,2]               -0.12    0.02  0.99 -2.02 -0.79 -0.12  0.57  1.82  4291    1
#> z[11,3]               -0.17    0.02  0.96 -2.01 -0.83 -0.17  0.46  1.73  2845    1
#> z[12,1]               -0.33    0.02  0.98 -2.34 -1.00 -0.32  0.34  1.51  3971    1
#> z[12,2]               -0.06    0.01  0.97 -2.01 -0.72 -0.05  0.60  1.85  4406    1
#> z[12,3]               -0.07    0.01  0.99 -1.98 -0.76 -0.07  0.59  1.82  4805    1
#> z[13,1]               -0.29    0.02  0.99 -2.16 -0.98 -0.30  0.41  1.70  4136    1
#> z[13,2]                0.09    0.02  0.94 -1.74 -0.55  0.09  0.70  2.03  3341    1
#> z[13,3]                0.03    0.02  0.96 -1.87 -0.59  0.00  0.67  1.93  3010    1
#> z[14,1]               -0.33    0.02  1.00 -2.24 -1.02 -0.33  0.36  1.63  3690    1
#> z[14,2]                0.13    0.01  0.98 -1.83 -0.52  0.14  0.78  2.03  4565    1
#> z[14,3]                0.06    0.01  0.98 -1.81 -0.63  0.04  0.73  1.99  4385    1
#> z[15,1]               -0.26    0.01  1.00 -2.31 -0.88 -0.25  0.39  1.82  4587    1
#> z[15,2]                0.07    0.02  1.00 -1.81 -0.62  0.05  0.72  2.05  3941    1
#> z[15,3]                0.12    0.02  0.98 -1.83 -0.55  0.12  0.79  2.00  3407    1
#> z[16,1]               -0.21    0.01  0.97 -2.08 -0.83 -0.21  0.41  1.73  4175    1
#> z[16,2]                0.14    0.02  0.95 -1.80 -0.49  0.15  0.75  2.00  3285    1
#> z[16,3]                0.25    0.01  0.98 -1.68 -0.44  0.25  0.91  2.21  4597    1
#> z[17,1]               -0.19    0.02  0.98 -2.09 -0.86 -0.22  0.50  1.66  3577    1
#> z[17,2]                0.16    0.01  0.96 -1.80 -0.47  0.15  0.81  2.09  4284    1
#> z[17,3]                0.30    0.02  0.90 -1.52 -0.31  0.32  0.90  2.13  3604    1
#> z[18,1]               -0.26    0.02  0.94 -2.12 -0.88 -0.26  0.36  1.59  3631    1
#> z[18,2]                0.28    0.02  0.95 -1.57 -0.37  0.28  0.90  2.19  3845    1
#> z[18,3]                0.40    0.02  0.95 -1.36 -0.23  0.39  1.03  2.21  3508    1
#> z[19,1]               -0.17    0.01  0.94 -2.02 -0.83 -0.17  0.50  1.69  4416    1
#> z[19,2]                0.35    0.02  0.97 -1.54 -0.28  0.35  1.01  2.19  3041    1
#> z[19,3]                0.39    0.01  0.99 -1.58 -0.27  0.39  1.07  2.35  4755    1
#> z[20,1]               -0.19    0.02  0.95 -2.05 -0.83 -0.20  0.45  1.67  3849    1
#> z[20,2]                0.00    0.02  0.95 -1.78 -0.66  0.00  0.65  1.87  3892    1
#> z[20,3]                0.33    0.01  0.99 -1.64 -0.30  0.32  0.97  2.25  5095    1
#> z[21,1]               -0.10    0.02  0.94 -1.97 -0.71 -0.08  0.52  1.64  3423    1
#> z[21,2]               -0.44    0.01  0.97 -2.33 -1.08 -0.46  0.20  1.49  4808    1
#> z[21,3]                0.37    0.01  0.93 -1.48 -0.25  0.37  1.02  2.20  4662    1
#> z[22,1]               -0.10    0.01  0.96 -2.02 -0.76 -0.09  0.54  1.85  4300    1
#> z[22,2]               -0.32    0.02  0.96 -2.19 -0.97 -0.34  0.32  1.63  3788    1
#> z[22,3]                0.49    0.02  0.94 -1.31 -0.17  0.48  1.15  2.24  3549    1
#> z[23,1]               -0.09    0.02  0.97 -1.95 -0.76 -0.09  0.57  1.83  3693    1
#> z[23,2]               -0.37    0.02  0.95 -2.27 -0.99 -0.36  0.26  1.53  3592    1
#> z[23,3]               -0.19    0.01  0.93 -1.95 -0.84 -0.20  0.43  1.70  4054    1
#> z[24,1]               -0.19    0.02  0.97 -2.02 -0.86 -0.18  0.49  1.67  3310    1
#> z[24,2]               -0.30    0.02  0.95 -2.09 -0.93 -0.32  0.32  1.61  3328    1
#> z[24,3]               -0.07    0.01  0.93 -1.88 -0.68 -0.05  0.57  1.85  4313    1
#> z[25,1]               -0.15    0.02  0.96 -2.10 -0.80 -0.14  0.47  1.74  3188    1
#> z[25,2]               -0.29    0.01  0.97 -2.25 -0.95 -0.28  0.34  1.63  4363    1
#> z[25,3]                0.08    0.02  0.95 -1.84 -0.57  0.07  0.72  1.98  3440    1
#> z[26,1]               -0.09    0.02  1.01 -2.04 -0.80 -0.08  0.60  1.90  3471    1
#> z[26,2]               -0.18    0.02  0.93 -2.05 -0.80 -0.18  0.46  1.61  3217    1
#> z[26,3]                0.18    0.01  0.89 -1.54 -0.43  0.19  0.81  1.94  3710    1
#> z[27,1]               -0.09    0.01  0.94 -1.97 -0.72 -0.09  0.50  1.73  3916    1
#>  [ reached 'max' / getOption("max.print") -- omitted 5679 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Aug 11 23:14:34 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

We can forecast

``` r

forecast(bvar_obj, pi = 0.68, show_all = TRUE)
```

![plot of chunk RW-2](figure/RW-2-1.png)![plot of chunk
RW-2](figure/RW-2-2.png)![plot of chunk RW-2](figure/RW-2-3.png)

Let us plot the log volatility estimates and predictions

``` r

stochastic_volatility_plot(bvar_obj, ci = 0.95, vol = "log_lambda")
```

![plot of chunk RW-3](figure/RW-3-1.png)![plot of chunk
RW-3](figure/RW-3-2.png)![plot of chunk RW-3](figure/RW-3-3.png)

Let us plot the estimates and predictions of the implied innovation
standard deviations

``` r

stochastic_volatility_plot(bvar_obj, vol = "sd")
```

![plot of chunk RW-4](figure/RW-4-1.png)![plot of chunk
RW-4](figure/RW-4-2.png)![plot of chunk RW-4](figure/RW-4-3.png)

We can also produce orthogonalized IRFs

``` r

IRF(bvar_obj, method = "OIRF", t=215, ci=0.68) #latest t
```

![plot of chunk RW-5](figure/RW-5-1.png)

plot of chunk RW-5

## References

Clark, T. E. (2011). Real-time density forecasts from Bayesian vector
autoregressions with stochastic volatility. *Journal of Business &
Economic Statistics*, 29(3), pp. 327–341.

Koop, G. and Korobilis, D. (2010). Bayesian multivariate time series
methods for empirical macroeconomics. *Foundations and Trends in
Econometrics*, 3(4), pp. 267–358.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630–650.
