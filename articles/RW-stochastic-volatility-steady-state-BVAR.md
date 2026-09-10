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

![plot of chunk RW-2](figure/RW-2-1.png)![plot of chunk
RW-2](figure/RW-2-2.png)![plot of chunk RW-2](figure/RW-2-3.png)

Now, let us fit the model. Note that we can use arguments from
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
such as `control` where we can tweak `max_treedepth` and `adapt_delta`.

``` r

bvar_obj <- fit(bvar_obj,
                H = 40,
                iter = 5000,
                warmup = 2500,
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
#>   delta pi.l1     1.27  0.01  0.15
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
#> r            0.02 -0.01  0.15
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
#> 4 chains, each with iter=5000; warmup=2500; thin=1; 
#> post-warmup draws per chain=2500, total post-warmup draws=10000.
#> 
#>                        mean se_mean    sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
#> beta[1,1]              1.27    0.00  0.06  1.16  1.23  1.27  1.31  1.38 11892    1
#> beta[1,2]              0.01    0.00  0.04 -0.06 -0.01  0.01  0.04  0.09 12347    1
#> beta[1,3]              0.15    0.00  0.08 -0.01  0.09  0.15  0.20  0.31 12386    1
#> beta[2,1]             -0.09    0.00  0.03 -0.16 -0.12 -0.09 -0.07 -0.02 13441    1
#> beta[2,2]              1.17    0.00  0.06  1.06  1.13  1.17  1.21  1.28 11519    1
#> beta[2,3]             -0.16    0.00  0.08 -0.32 -0.22 -0.16 -0.11 -0.01 12669    1
#> beta[3,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.01  0.04 14476    1
#> beta[3,2]             -0.01    0.00  0.02 -0.04 -0.02 -0.01  0.00  0.02 14398    1
#> beta[3,3]              1.04    0.00  0.06  0.92  1.00  1.04  1.08  1.16 10977    1
#> beta[4,1]             -0.28    0.00  0.06 -0.39 -0.32 -0.28 -0.24 -0.17 12130    1
#> beta[4,2]              0.02    0.00  0.04 -0.06 -0.01  0.02  0.04  0.09 12696    1
#> beta[4,3]             -0.10    0.00  0.08 -0.26 -0.16 -0.10 -0.05  0.06 12704    1
#> beta[5,1]              0.07    0.00  0.03  0.01  0.05  0.07  0.09  0.13 13670    1
#> beta[5,2]             -0.23    0.00  0.05 -0.34 -0.27 -0.23 -0.20 -0.13 11734    1
#> beta[5,3]              0.17    0.00  0.07  0.02  0.12  0.17  0.22  0.32 12558    1
#> beta[6,1]              0.00    0.00  0.02 -0.03 -0.01  0.00  0.01  0.03 13355    1
#> beta[6,2]              0.02    0.00  0.02 -0.01  0.01  0.02  0.03  0.05 14732    1
#> beta[6,3]             -0.11    0.00  0.06 -0.22 -0.15 -0.11 -0.07  0.01 11335    1
#> Psi[1,1]               2.00    0.00  0.05  1.90  1.96  2.00  2.03  2.10 26226    1
#> Psi[2,1]               4.29    0.00  0.17  3.94  4.17  4.29  4.41  4.63 22511    1
#> Psi[3,1]               3.49    0.00  0.32  2.85  3.28  3.49  3.71  4.11 21213    1
#> z[1,1]                -0.06    0.00  0.25 -0.52 -0.23 -0.07  0.10  0.45 15649    1
#> z[1,2]                 0.64    0.00  0.27  0.12  0.46  0.63  0.81  1.18 16469    1
#> z[1,3]                -0.65    0.00  0.33 -1.28 -0.87 -0.66 -0.44  0.01 14954    1
#> z[2,1]                -0.13    0.01  0.99 -2.07 -0.81 -0.13  0.54  1.81 21560    1
#> z[2,2]                 0.21    0.01  1.01 -1.79 -0.49  0.22  0.89  2.17 17746    1
#> z[2,3]                 0.01    0.01  0.99 -1.92 -0.65  0.02  0.68  1.94 20764    1
#> z[3,1]                -0.04    0.01  0.99 -1.98 -0.72 -0.05  0.63  1.90 21485    1
#> z[3,2]                 0.25    0.01  0.98 -1.69 -0.41  0.23  0.91  2.17 20782    1
#> z[3,3]                -0.02    0.01  0.98 -1.94 -0.69 -0.02  0.64  1.89 24181    1
#> z[4,1]                -0.09    0.01  0.97 -2.00 -0.76 -0.09  0.56  1.80 23805    1
#> z[4,2]                -0.43    0.01  0.96 -2.33 -1.06 -0.42  0.21  1.45 23658    1
#> z[4,3]                -0.25    0.01  0.98 -2.20 -0.91 -0.26  0.40  1.71 21902    1
#> z[5,1]                 0.01    0.01  0.99 -1.90 -0.67  0.02  0.68  1.92 26090    1
#> z[5,2]                -0.45    0.01  0.96 -2.29 -1.10 -0.46  0.21  1.42 20505    1
#> z[5,3]                -0.21    0.01  0.96 -2.08 -0.86 -0.21  0.44  1.71 22716    1
#> z[6,1]                 0.02    0.01  1.01 -1.96 -0.65  0.02  0.69  2.01 26609    1
#> z[6,2]                -0.33    0.01  0.96 -2.17 -0.98 -0.34  0.33  1.57 23067    1
#> z[6,3]                -0.13    0.01  0.97 -2.04 -0.78 -0.14  0.52  1.76 22657    1
#> z[7,1]                 0.10    0.01  0.98 -1.82 -0.56  0.11  0.77  1.99 22424    1
#> z[7,2]                -0.21    0.01  0.96 -2.08 -0.86 -0.21  0.44  1.63 25214    1
#> z[7,3]                -0.03    0.01  0.97 -1.94 -0.69 -0.03  0.62  1.88 25366    1
#> z[8,1]                 0.19    0.01  0.97 -1.71 -0.47  0.19  0.86  2.08 20669    1
#> z[8,2]                -0.26    0.01  0.94 -2.12 -0.89 -0.26  0.36  1.57 20193    1
#> z[8,3]                 0.07    0.01  0.96 -1.82 -0.57  0.08  0.72  1.92 22687    1
#> z[9,1]                 0.15    0.01  0.98 -1.80 -0.50  0.16  0.81  2.09 21445    1
#> z[9,2]                -0.13    0.01  0.96 -2.02 -0.77 -0.14  0.50  1.77 23798    1
#> z[9,3]                 0.20    0.01  0.94 -1.65 -0.43  0.20  0.83  2.06 21499    1
#> z[10,1]               -0.19    0.01  0.97 -2.08 -0.85 -0.18  0.48  1.71 24294    1
#> z[10,2]               -0.12    0.01  0.97 -2.01 -0.78 -0.11  0.53  1.78 23172    1
#> z[10,3]               -0.04    0.01  0.94 -1.85 -0.68 -0.05  0.58  1.81 23630    1
#> z[11,1]               -0.18    0.01  0.94 -2.02 -0.81 -0.18  0.46  1.68 21976    1
#> z[11,2]               -0.11    0.01  0.96 -1.99 -0.75 -0.11  0.53  1.77 24810    1
#> z[11,3]               -0.17    0.01  0.98 -2.07 -0.84 -0.18  0.51  1.72 23164    1
#> z[12,1]               -0.31    0.01  0.96 -2.19 -0.96 -0.32  0.32  1.56 25620    1
#> z[12,2]               -0.03    0.01  0.96 -1.91 -0.68 -0.04  0.60  1.87 26873    1
#> z[12,3]               -0.11    0.01  0.95 -1.96 -0.75 -0.11  0.53  1.77 23155    1
#> z[13,1]               -0.28    0.01  0.98 -2.23 -0.94 -0.28  0.36  1.67 25172    1
#> z[13,2]                0.08    0.01  0.96 -1.80 -0.57  0.07  0.73  1.97 20395    1
#> z[13,3]                0.04    0.01  0.96 -1.85 -0.62  0.06  0.70  1.93 23029    1
#> z[14,1]               -0.33    0.01  0.95 -2.17 -0.98 -0.32  0.31  1.55 24882    1
#> z[14,2]                0.13    0.01  0.96 -1.75 -0.51  0.12  0.77  2.01 23031    1
#> z[14,3]                0.03    0.01  0.98 -1.90 -0.63  0.03  0.69  1.94 26076    1
#> z[15,1]               -0.25    0.01  0.96 -2.15 -0.90 -0.26  0.41  1.61 24749    1
#> z[15,2]                0.06    0.01  0.94 -1.74 -0.57  0.05  0.71  1.87 20847    1
#> z[15,3]                0.12    0.01  0.95 -1.73 -0.52  0.11  0.76  1.97 25490    1
#> z[16,1]               -0.22    0.01  0.97 -2.10 -0.86 -0.21  0.42  1.67 23541    1
#> z[16,2]                0.13    0.01  0.97 -1.78 -0.52  0.13  0.78  2.03 23247    1
#> z[16,3]                0.26    0.01  0.97 -1.66 -0.38  0.27  0.91  2.16 25173    1
#> z[17,1]               -0.21    0.01  0.97 -2.08 -0.87 -0.20  0.45  1.66 23062    1
#> z[17,2]                0.17    0.01  0.96 -1.69 -0.48  0.18  0.82  2.04 24084    1
#> z[17,3]                0.29    0.01  0.95 -1.55 -0.35  0.29  0.94  2.17 22100    1
#> z[18,1]               -0.28    0.01  0.96 -2.20 -0.92 -0.29  0.35  1.59 23699    1
#> z[18,2]                0.26    0.01  0.98 -1.64 -0.40  0.25  0.93  2.21 22581    1
#> z[18,3]                0.41    0.01  0.95 -1.49 -0.22  0.41  1.06  2.27 22083    1
#> z[19,1]               -0.18    0.01  0.97 -2.13 -0.83 -0.19  0.47  1.75 26050    1
#> z[19,2]                0.37    0.01  0.96 -1.52 -0.28  0.37  1.01  2.27 22516    1
#> z[19,3]                0.37    0.01  0.97 -1.51 -0.29  0.38  1.02  2.26 24081    1
#> z[20,1]               -0.18    0.01  0.96 -2.07 -0.84 -0.19  0.47  1.72 25148    1
#> z[20,2]                0.03    0.01  0.96 -1.85 -0.61  0.02  0.68  1.89 25976    1
#> z[20,3]                0.32    0.01  0.97 -1.59 -0.34  0.33  0.97  2.23 27197    1
#> z[21,1]               -0.10    0.01  0.97 -2.02 -0.75 -0.09  0.54  1.82 26111    1
#> z[21,2]               -0.44    0.01  0.94 -2.26 -1.08 -0.45  0.19  1.43 25151    1
#> z[21,3]                0.36    0.01  0.95 -1.51 -0.29  0.37  1.02  2.17 24360    1
#> z[22,1]               -0.11    0.01  0.96 -1.98 -0.75 -0.11  0.55  1.80 23932    1
#> z[22,2]               -0.32    0.01  0.96 -2.22 -0.97 -0.31  0.34  1.52 24630    1
#> z[22,3]                0.50    0.01  0.94 -1.34 -0.13  0.51  1.12  2.33 23799    1
#> z[23,1]               -0.06    0.01  0.96 -1.91 -0.72 -0.06  0.59  1.81 25239    1
#> z[23,2]               -0.36    0.01  0.96 -2.24 -1.00 -0.35  0.29  1.50 21333    1
#> z[23,3]               -0.17    0.01  0.93 -2.00 -0.81 -0.17  0.45  1.63 22604    1
#> z[24,1]               -0.18    0.01  0.96 -2.09 -0.83 -0.19  0.46  1.69 23863    1
#> z[24,2]               -0.29    0.01  0.96 -2.18 -0.93 -0.28  0.35  1.62 26061    1
#> z[24,3]               -0.07    0.01  0.96 -1.94 -0.72 -0.07  0.57  1.83 23637    1
#> z[25,1]               -0.17    0.01  0.98 -2.09 -0.84 -0.17  0.48  1.76 23594    1
#> z[25,2]               -0.31    0.01  0.94 -2.11 -0.95 -0.32  0.33  1.54 24090    1
#> z[25,3]                0.07    0.01  0.92 -1.72 -0.55  0.07  0.70  1.87 20243    1
#> z[26,1]               -0.07    0.01  0.98 -1.95 -0.73 -0.06  0.59  1.85 24714    1
#> z[26,2]               -0.19    0.01  0.97 -2.10 -0.85 -0.18  0.47  1.67 23137    1
#> z[26,3]                0.20    0.01  0.93 -1.62 -0.42  0.21  0.82  2.03 23648    1
#> z[27,1]               -0.08    0.01  0.96 -1.95 -0.74 -0.08  0.59  1.77 20540    1
#>  [ reached 'max' / getOption("max.print") -- omitted 6438 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Thu Sep 10 02:39:21 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Let us plot the posterior of the steady-state \\mu_t=mu=\Psi\\

``` r

stanfit <- bvar_obj$fit$stan

rstan::plot(stanfit, pars=c("mu[1,1]",
                            "mu[1,2]",
                            "mu[1,3]"), plotfun="hist")
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![plot of chunk RW-3](figure/RW-3-1.png)

plot of chunk RW-3

We can forecast

``` r

forecast(bvar_obj, pi = 0.68)
```

![plot of chunk RW-4](figure/RW-4-1.png)![plot of chunk
RW-4](figure/RW-4-2.png)![plot of chunk RW-4](figure/RW-4-3.png)

Let us plot the log volatility estimates and predictions

``` r

stochastic_volatility_plot(bvar_obj, ci = 0.95, vol = "log_lambda")
```

![plot of chunk RW-5](figure/RW-5-1.png)![plot of chunk
RW-5](figure/RW-5-2.png)![plot of chunk RW-5](figure/RW-5-3.png)

Let us plot the estimates and predictions of the implied reduced-form
innovation standard deviations

``` r

stochastic_volatility_plot(bvar_obj, vol = "sd")
```

![plot of chunk RW-6](figure/RW-6-1.png)![plot of chunk
RW-6](figure/RW-6-2.png)![plot of chunk RW-6](figure/RW-6-3.png)

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
