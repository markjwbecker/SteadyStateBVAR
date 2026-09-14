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

par(mfrow=c(3,1))
steady_state_priors_plot(bvar_obj, interval = 0.95)
```

![plot of chunk RW-2](figure/RW-2-1.png)

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
#> delta pi     0.08 -0.01  0.02
#> u           -0.01  0.02 -0.01
#> r            0.02 -0.01  0.16
#> --------------------------------------------------------------------------------
#> 
#> 
#> A
#> --------------------------------------------------------------------------------          
#>            delta pi    u r
#>   delta pi     1.00 0.00 0
#>   u            0.11 1.00 0
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
#>                        mean se_mean     sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
#> beta[1,1]              1.27    0.00   0.06  1.16  1.23  1.27  1.31  1.38 16476    1
#> beta[1,2]              0.02    0.00   0.04 -0.06 -0.01  0.02  0.04  0.09 15161    1
#> beta[1,3]              0.15    0.00   0.08 -0.01  0.09  0.15  0.20  0.31 15917    1
#> beta[2,1]             -0.09    0.00   0.03 -0.16 -0.12 -0.09 -0.07 -0.02 15409    1
#> beta[2,2]              1.17    0.00   0.06  1.06  1.13  1.17  1.21  1.28 13476    1
#> beta[2,3]             -0.16    0.00   0.08 -0.32 -0.21 -0.16 -0.11 -0.01 14382    1
#> beta[3,1]              0.00    0.00   0.02 -0.03 -0.01  0.00  0.01  0.03 16445    1
#> beta[3,2]             -0.01    0.00   0.02 -0.04 -0.02 -0.01  0.00  0.02 16618    1
#> beta[3,3]              1.04    0.00   0.06  0.92  1.00  1.04  1.08  1.16 15790    1
#> beta[4,1]             -0.28    0.00   0.06 -0.39 -0.32 -0.28 -0.24 -0.17 16340    1
#> beta[4,2]              0.02    0.00   0.04 -0.06 -0.01  0.02  0.04  0.09 15041    1
#> beta[4,3]             -0.10    0.00   0.08 -0.26 -0.16 -0.10 -0.05  0.06 16209    1
#> beta[5,1]              0.07    0.00   0.03  0.01  0.05  0.07  0.09  0.14 15559    1
#> beta[5,2]             -0.23    0.00   0.05 -0.34 -0.27 -0.23 -0.20 -0.13 13369    1
#> beta[5,3]              0.17    0.00   0.07  0.02  0.12  0.17  0.22  0.32 14638    1
#> beta[6,1]              0.00    0.00   0.02 -0.03 -0.01  0.00  0.01  0.03 17407    1
#> beta[6,2]              0.02    0.00   0.02 -0.01  0.01  0.02  0.03  0.05 16870    1
#> beta[6,3]             -0.11    0.00   0.06 -0.22 -0.15 -0.11 -0.07  0.00 15716    1
#> Psi[1,1]               2.00    0.00   0.05  1.90  1.96  2.00  2.03  2.10 30761    1
#> Psi[2,1]               4.29    0.00   0.18  3.94  4.17  4.29  4.41  4.63 24765    1
#> Psi[3,1]               3.49    0.00   0.32  2.85  3.27  3.50  3.72  4.12 23789    1
#> z[1,1]                -0.06    0.00   0.24 -0.51 -0.23 -0.07  0.10  0.44 22146    1
#> z[1,2]                 0.64    0.00   0.27  0.14  0.46  0.63  0.81  1.18 18139    1
#> z[1,3]                -0.66    0.00   0.33 -1.28 -0.88 -0.67 -0.44  0.03 16521    1
#> z[2,1]                -0.12    0.01   0.97 -2.00 -0.78 -0.13  0.53  1.79 27846    1
#> z[2,2]                 0.20    0.01   0.99 -1.76 -0.48  0.19  0.87  2.17 23937    1
#> z[2,3]                 0.02    0.01   1.00 -1.96 -0.66  0.03  0.69  1.96 24664    1
#> z[3,1]                -0.03    0.01   0.98 -1.94 -0.70 -0.04  0.64  1.87 25559    1
#> z[3,2]                 0.25    0.01   0.97 -1.67 -0.39  0.25  0.91  2.14 27050    1
#> z[3,3]                -0.03    0.01   0.98 -1.96 -0.68 -0.03  0.62  1.86 25219    1
#> z[4,1]                -0.09    0.01   0.98 -2.02 -0.74 -0.09  0.57  1.83 27318    1
#> z[4,2]                -0.43    0.01   0.96 -2.31 -1.07 -0.43  0.23  1.45 27719    1
#> z[4,3]                -0.26    0.01   0.96 -2.14 -0.90 -0.26  0.38  1.62 26686    1
#> z[5,1]                 0.01    0.01   0.97 -1.90 -0.64  0.00  0.67  1.90 30760    1
#> z[5,2]                -0.45    0.01   0.96 -2.32 -1.10 -0.45  0.19  1.42 25173    1
#> z[5,3]                -0.20    0.01   0.97 -2.10 -0.85 -0.19  0.46  1.68 25095    1
#> z[6,1]                 0.00    0.01   0.98 -1.91 -0.65  0.01  0.65  1.92 24233    1
#> z[6,2]                -0.33    0.01   0.97 -2.24 -0.99 -0.33  0.33  1.55 26325    1
#> z[6,3]                -0.12    0.01   0.97 -1.99 -0.78 -0.12  0.53  1.77 24723    1
#> z[7,1]                 0.10    0.01   0.97 -1.82 -0.57  0.11  0.76  2.01 27922    1
#> z[7,2]                -0.21    0.01   0.94 -2.05 -0.85 -0.21  0.43  1.64 27436    1
#> z[7,3]                -0.03    0.01   0.97 -1.95 -0.68 -0.02  0.62  1.89 27116    1
#> z[8,1]                 0.19    0.01   0.98 -1.72 -0.46  0.20  0.86  2.10 26477    1
#> z[8,2]                -0.25    0.01   0.96 -2.15 -0.91 -0.24  0.41  1.62 28646    1
#> z[8,3]                 0.07    0.01   0.96 -1.80 -0.58  0.07  0.72  1.96 27381    1
#> z[9,1]                 0.15    0.01   0.98 -1.76 -0.51  0.15  0.81  2.06 29987    1
#> z[9,2]                -0.13    0.01   0.96 -2.04 -0.76 -0.13  0.50  1.75 27972    1
#> z[9,3]                 0.20    0.01   0.96 -1.64 -0.44  0.21  0.84  2.11 24923    1
#> z[10,1]               -0.19    0.01   0.99 -2.13 -0.85 -0.19  0.48  1.76 28664    1
#> z[10,2]               -0.11    0.01   0.96 -2.01 -0.77 -0.11  0.55  1.75 26042    1
#> z[10,3]               -0.05    0.01   0.95 -1.92 -0.70 -0.06  0.60  1.81 26019    1
#> z[11,1]               -0.17    0.01   0.95 -2.05 -0.82 -0.17  0.46  1.71 26974    1
#> z[11,2]               -0.12    0.01   0.96 -2.01 -0.77 -0.12  0.52  1.75 25376    1
#> z[11,3]               -0.16    0.01   0.97 -2.07 -0.81 -0.16  0.50  1.73 24186    1
#> z[12,1]               -0.33    0.01   0.96 -2.21 -0.98 -0.33  0.33  1.54 26100    1
#> z[12,2]               -0.04    0.01   0.96 -1.94 -0.69 -0.04  0.61  1.85 28753    1
#> z[12,3]               -0.09    0.01   0.96 -1.97 -0.75 -0.10  0.56  1.76 24890    1
#> z[13,1]               -0.27    0.01   0.98 -2.17 -0.93 -0.27  0.41  1.63 26861    1
#> z[13,2]                0.10    0.01   0.98 -1.83 -0.56  0.10  0.76  1.98 26670    1
#> z[13,3]                0.04    0.01   0.97 -1.86 -0.61  0.04  0.69  1.92 25677    1
#> z[14,1]               -0.31    0.01   0.96 -2.18 -0.96 -0.32  0.34  1.56 28207    1
#> z[14,2]                0.12    0.01   0.95 -1.73 -0.52  0.12  0.76  1.97 27832    1
#> z[14,3]                0.02    0.01   0.97 -1.88 -0.62  0.02  0.68  1.91 25876    1
#> z[15,1]               -0.25    0.01   0.96 -2.12 -0.89 -0.25  0.40  1.64 28902    1
#> z[15,2]                0.08    0.01   0.95 -1.78 -0.57  0.08  0.73  1.93 27154    1
#> z[15,3]                0.13    0.01   0.97 -1.81 -0.53  0.13  0.79  2.01 25886    1
#> z[16,1]               -0.21    0.01   0.96 -2.11 -0.85 -0.21  0.45  1.68 26868    1
#> z[16,2]                0.13    0.01   0.97 -1.78 -0.51  0.12  0.77  2.05 26898    1
#> z[16,3]                0.25    0.01   0.96 -1.63 -0.38  0.26  0.88  2.15 28430    1
#> z[17,1]               -0.21    0.01   0.96 -2.10 -0.85 -0.21  0.43  1.65 27073    1
#> z[17,2]                0.16    0.01   0.96 -1.72 -0.48  0.16  0.81  2.06 24423    1
#> z[17,3]                0.30    0.01   0.96 -1.57 -0.36  0.29  0.95  2.16 27675    1
#> z[18,1]               -0.27    0.01   0.97 -2.18 -0.92 -0.27  0.38  1.63 28558    1
#> z[18,2]                0.27    0.01   0.96 -1.63 -0.37  0.27  0.91  2.15 26213    1
#> z[18,3]                0.41    0.01   0.96 -1.48 -0.24  0.42  1.07  2.28 24910    1
#> z[19,1]               -0.19    0.01   0.98 -2.12 -0.85 -0.18  0.47  1.73 27908    1
#> z[19,2]                0.38    0.01   0.95 -1.48 -0.25  0.37  1.01  2.27 25819    1
#> z[19,3]                0.37    0.01   0.96 -1.50 -0.28  0.37  1.02  2.26 23893    1
#> z[20,1]               -0.18    0.01   0.96 -2.06 -0.84 -0.18  0.48  1.70 22746    1
#> z[20,2]                0.02    0.01   0.95 -1.82 -0.63  0.02  0.67  1.89 26998    1
#> z[20,3]                0.32    0.01   0.95 -1.53 -0.31  0.32  0.97  2.17 27698    1
#> z[21,1]               -0.10    0.01   0.97 -2.00 -0.76 -0.10  0.56  1.80 26434    1
#> z[21,2]               -0.45    0.01   0.94 -2.29 -1.09 -0.45  0.18  1.38 29633    1
#> z[21,3]                0.37    0.01   0.94 -1.48 -0.27  0.37  1.00  2.21 24766    1
#> z[22,1]               -0.11    0.01   0.96 -2.03 -0.75 -0.11  0.53  1.77 26218    1
#> z[22,2]               -0.32    0.01   0.95 -2.16 -0.98 -0.33  0.32  1.57 25764    1
#> z[22,3]                0.51    0.01   0.94 -1.34 -0.13  0.51  1.14  2.37 24677    1
#> z[23,1]               -0.05    0.01   0.96 -1.92 -0.69 -0.05  0.59  1.82 24336    1
#> z[23,2]               -0.35    0.01   0.96 -2.22 -1.00 -0.35  0.31  1.52 29899    1
#> z[23,3]               -0.17    0.01   0.95 -2.06 -0.80 -0.17  0.46  1.70 25636    1
#> z[24,1]               -0.19    0.01   0.98 -2.13 -0.83 -0.19  0.46  1.75 27712    1
#> z[24,2]               -0.30    0.01   0.94 -2.15 -0.94 -0.30  0.34  1.55 26602    1
#> z[24,3]               -0.06    0.01   0.94 -1.92 -0.69 -0.05  0.58  1.76 27102    1
#> z[25,1]               -0.16    0.01   0.96 -2.08 -0.82 -0.16  0.49  1.72 26793    1
#> z[25,2]               -0.30    0.01   0.95 -2.18 -0.95 -0.30  0.35  1.56 26469    1
#> z[25,3]                0.06    0.01   0.94 -1.79 -0.57  0.06  0.69  1.89 29444    1
#> z[26,1]               -0.08    0.01   0.96 -1.98 -0.73 -0.07  0.58  1.82 24694    1
#> z[26,2]               -0.21    0.01   0.97 -2.12 -0.86 -0.21  0.46  1.67 26089    1
#> z[26,3]                0.20    0.01   0.94 -1.67 -0.43  0.20  0.83  2.05 28419    1
#> z[27,1]               -0.09    0.01   0.96 -1.97 -0.73 -0.09  0.55  1.80 27863    1
#>  [ reached 'max' / getOption("max.print") -- omitted 6438 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Mon Sep 14 02:11:35 2026.
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

We can forecast

``` r

par(mfrow=c(3,1))
forecast(bvar_obj, pi = 0.68)
```

![plot of chunk RW-4](figure/RW-4-1.png)

``` r

par(mfrow=c(1,1))
```

Let us plot the log volatility estimates and predictions

``` r

stochastic_volatility_plot(bvar_obj, ci = 0.95, vol = "log_lambda")
```

![plot of chunk RW-5](figure/RW-5-1.png)

![plot of chunk RW-5](figure/RW-5-2.png)

![plot of chunk RW-5](figure/RW-5-3.png)

Let us plot the estimates and predictions of the implied reduced-form
innovation standard deviations

``` r

stochastic_volatility_plot(bvar_obj, vol = "sd")
```

![plot of chunk RW-6](figure/RW-6-1.png)

![plot of chunk RW-6](figure/RW-6-2.png)

![plot of chunk RW-6](figure/RW-6-3.png)

We can also produce orthogonalized IRFs

``` r

IRF(bvar_obj, method = "OIRF", t=215, ci=0.68) #using Sigma_u,t=215
```

![plot of chunk RW-7](figure/RW-7-1.png)

## References

Clark, T. E. (2011). Real-time density forecasts from Bayesian vector
autoregressions with stochastic volatility. *Journal of Business &
Economic Statistics*, 29(3), pp. 327–341.

Koop, G. and Korobilis, D. (2010). Bayesian multivariate time series
methods for empirical macroeconomics. *Foundations and Trends in
Econometrics*, 3(4), pp. 267–358.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630–650.
