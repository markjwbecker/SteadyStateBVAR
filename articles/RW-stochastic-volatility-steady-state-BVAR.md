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

Now, for the steady-state coefficients we use some toy values (let us
pretend that they are expert based). Remember that we only have a
constant now, so \\q=1\\ and therefore \\\Psi\\ only has one column
\\\psi_1=\Psi\\. Since \\d_t = 1 \\ \forall \\ t\\, we have \\\Psi d_t =
\mu_t\\ which simplifies to \\\Psi = \mu\\ and as such we can directly
interpret \\\Psi\\ as the unconditional mean.

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

SV_priors_RW <- list(
                     theta_A             =  rep(0, n_free_params_A),
                     Omega_A             =  diag(10, n_free_params_A),
                     mu_log_lambda_0     =  rep(0, k),
                     sigma2_log_lambda_0 =  rep(10, k),
                     alpha_phi           =  rep(5, k),
                     beta_phi            = (rep(5, k) - 1) * rep(0.1, k)
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
#> delta pi     0.10 -0.01  0.03
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

Note that you can always look at the `stanfit` object
`bvar_obj$fit$stan` directly if you want

``` r

print(bvar_obj$fit$stan)
#> Inference for Stan model: steady_state_bvar_RW_stochastic_volatility.
#> 2 chains, each with iter=4000; warmup=1000; thin=1; 
#> post-warmup draws per chain=3000, total post-warmup draws=6000.
#> 
#>                        mean se_mean     sd  2.5%   25%   50%   75%  97.5% n_eff Rhat
#> beta[1,1]              1.26    0.00   0.06  1.15  1.22  1.26  1.30   1.37  6895    1
#> beta[1,2]              0.01    0.00   0.04 -0.06 -0.01  0.01  0.04   0.09  7131    1
#> beta[1,3]              0.15    0.00   0.08 -0.01  0.09  0.15  0.20   0.31  6400    1
#> beta[2,1]             -0.09    0.00   0.04 -0.16 -0.12 -0.09 -0.07  -0.02  6833    1
#> beta[2,2]              1.17    0.00   0.06  1.05  1.13  1.17  1.20   1.28  5885    1
#> beta[2,3]             -0.16    0.00   0.08 -0.31 -0.21 -0.16 -0.11  -0.01  5652    1
#> beta[3,1]              0.00    0.00   0.02 -0.03 -0.01  0.00  0.01   0.04  6774    1
#> beta[3,2]             -0.01    0.00   0.02 -0.04 -0.02 -0.01  0.00   0.02  7244    1
#> beta[3,3]              1.04    0.00   0.06  0.93  1.00  1.04  1.08   1.16  6714    1
#> beta[4,1]             -0.27    0.00   0.06 -0.38 -0.31 -0.28 -0.24  -0.16  6773    1
#> beta[4,2]              0.02    0.00   0.04 -0.05 -0.01  0.02  0.04   0.09  7013    1
#> beta[4,3]             -0.10    0.00   0.08 -0.26 -0.16 -0.10 -0.05   0.06  6316    1
#> beta[5,1]              0.07    0.00   0.03  0.01  0.05  0.07  0.09   0.14  6789    1
#> beta[5,2]             -0.23    0.00   0.05 -0.33 -0.26 -0.23 -0.19  -0.12  5890    1
#> beta[5,3]              0.17    0.00   0.07  0.02  0.12  0.17  0.22   0.32  5839    1
#> beta[6,1]              0.00    0.00   0.02 -0.03 -0.01  0.00  0.01   0.03  6697    1
#> beta[6,2]              0.02    0.00   0.02 -0.01  0.01  0.02  0.03   0.05  7342    1
#> beta[6,3]             -0.11    0.00   0.06 -0.22 -0.14 -0.11 -0.07   0.00  6598    1
#> Psi[1,1]               2.00    0.00   0.05  1.90  1.96  2.00  2.03   2.10  9066    1
#> Psi[2,1]               4.28    0.00   0.17  3.94  4.16  4.28  4.40   4.61  8911    1
#> Psi[3,1]               3.50    0.00   0.33  2.85  3.28  3.51  3.73   4.12 10277    1
#> z[1,1]                -0.72    0.00   0.17 -1.04 -0.84 -0.73 -0.61  -0.36  8156    1
#> z[1,2]                -0.22    0.00   0.19 -0.57 -0.35 -0.22 -0.09   0.17  7852    1
#> z[1,3]                -0.56    0.00   0.23 -0.99 -0.71 -0.56 -0.41  -0.10  7357    1
#> z[2,1]                -0.18    0.01   0.98 -2.05 -0.86 -0.17  0.50   1.71 10445    1
#> z[2,2]                 0.12    0.01   0.98 -1.77 -0.53  0.12  0.76   2.08 10198    1
#> z[2,3]                 0.06    0.01   0.97 -1.86 -0.60  0.06  0.73   1.91  7956    1
#> z[3,1]                -0.07    0.01   1.01 -2.11 -0.71 -0.09  0.59   1.88 10943    1
#> z[3,2]                 0.20    0.01   0.97 -1.71 -0.46  0.18  0.86   2.10  9967    1
#> z[3,3]                 0.01    0.01   0.98 -1.89 -0.64  0.00  0.67   1.94 10301    1
#> z[4,1]                -0.14    0.01   0.99 -2.08 -0.83 -0.13  0.55   1.76 10562    1
#> z[4,2]                -0.49    0.01   0.95 -2.34 -1.13 -0.50  0.14   1.37  8577    1
#> z[4,3]                -0.26    0.01   0.95 -2.13 -0.89 -0.26  0.38   1.63 11460    1
#> z[5,1]                -0.03    0.01   0.97 -1.91 -0.69 -0.03  0.64   1.87 12078    1
#> z[5,2]                -0.52    0.01   0.98 -2.45 -1.17 -0.51  0.14   1.40  9936    1
#> z[5,3]                -0.18    0.01   0.96 -2.13 -0.80 -0.16  0.47   1.70 11936    1
#> z[6,1]                 0.00    0.01   0.97 -1.85 -0.68  0.00  0.65   1.89  8108    1
#> z[6,2]                -0.39    0.01   0.93 -2.22 -1.01 -0.40  0.25   1.42 11658    1
#> z[6,3]                -0.12    0.01   0.95 -1.98 -0.76 -0.11  0.52   1.74 10251    1
#> z[7,1]                 0.09    0.01   0.99 -1.82 -0.58  0.07  0.74   2.08  9955    1
#> z[7,2]                -0.24    0.01   0.96 -2.15 -0.87 -0.25  0.39   1.63 10208    1
#> z[7,3]                -0.04    0.01   0.96 -1.87 -0.70 -0.05  0.60   1.84 10731    1
#> z[8,1]                 0.22    0.01   0.95 -1.67 -0.41  0.23  0.86   2.06  8920    1
#> z[8,2]                -0.27    0.01   0.96 -2.18 -0.93 -0.28  0.38   1.62 10401    1
#> z[8,3]                 0.08    0.01   0.98 -1.87 -0.59  0.08  0.74   2.01 11188    1
#> z[9,1]                 0.21    0.01   0.95 -1.69 -0.43  0.21  0.86   2.08  9293    1
#> z[9,2]                -0.13    0.01   0.97 -2.03 -0.79 -0.14  0.52   1.78 10934    1
#> z[9,3]                 0.23    0.01   0.99 -1.70 -0.45  0.22  0.89   2.15  9093    1
#> z[10,1]               -0.18    0.01   0.94 -2.03 -0.80 -0.19  0.45   1.66 10068    1
#> z[10,2]               -0.12    0.01   0.95 -1.99 -0.77 -0.13  0.53   1.71 10766    1
#> z[10,3]               -0.04    0.01   0.93 -1.88 -0.66 -0.03  0.59   1.75 11260    1
#> z[11,1]               -0.16    0.01   0.94 -2.02 -0.81 -0.16  0.48   1.67 10572    1
#> z[11,2]               -0.15    0.01   0.94 -1.99 -0.79 -0.14  0.47   1.69 11174    1
#> z[11,3]               -0.18    0.01   0.96 -2.07 -0.83 -0.18  0.47   1.73 10249    1
#> z[12,1]               -0.34    0.01   0.95 -2.19 -0.98 -0.34  0.31   1.51  8526    1
#> z[12,2]               -0.03    0.01   0.96 -1.86 -0.66 -0.02  0.61   1.82 10875    1
#> z[12,3]               -0.11    0.01   0.95 -1.97 -0.76 -0.12  0.53   1.75  9911    1
#> z[13,1]               -0.29    0.01   0.95 -2.17 -0.93 -0.29  0.35   1.58 11676    1
#> z[13,2]                0.11    0.01   0.95 -1.73 -0.53  0.11  0.75   1.92 12017    1
#> z[13,3]                0.04    0.01   0.94 -1.81 -0.60  0.05  0.68   1.87 10738    1
#> z[14,1]               -0.34    0.01   0.96 -2.23 -1.00 -0.34  0.31   1.52  9817    1
#> z[14,2]                0.18    0.01   0.97 -1.69 -0.49  0.19  0.83   2.08 10687    1
#> z[14,3]                0.03    0.01   1.00 -1.91 -0.65  0.03  0.69   2.01 10555    1
#> z[15,1]               -0.24    0.01   0.99 -2.16 -0.89 -0.24  0.42   1.67 11969    1
#> z[15,2]                0.06    0.01   0.95 -1.78 -0.58  0.06  0.70   1.97 10875    1
#> z[15,3]                0.12    0.01   0.93 -1.68 -0.50  0.11  0.74   1.94 10451    1
#> z[16,1]               -0.19    0.01   0.96 -2.07 -0.85 -0.20  0.48   1.70  9789    1
#> z[16,2]                0.15    0.01   0.95 -1.72 -0.48  0.16  0.81   1.98 10200    1
#> z[16,3]                0.27    0.01   0.95 -1.59 -0.36  0.28  0.89   2.14  9749    1
#> z[17,1]               -0.19    0.01   0.96 -2.10 -0.80 -0.18  0.45   1.70 11333    1
#> z[17,2]                0.19    0.01   0.97 -1.71 -0.47  0.21  0.86   2.10 10313    1
#> z[17,3]                0.31    0.01   0.94 -1.54 -0.32  0.32  0.97   2.09  9666    1
#> z[18,1]               -0.29    0.01   0.98 -2.22 -0.96 -0.30  0.37   1.58 11470    1
#> z[18,2]                0.30    0.01   0.94 -1.55 -0.35  0.31  0.94   2.14 10434    1
#> z[18,3]                0.44    0.01   0.95 -1.46 -0.21  0.44  1.08   2.30  9188    1
#> z[19,1]               -0.16    0.01   0.94 -2.00 -0.80 -0.17  0.47   1.68  9736    1
#> z[19,2]                0.44    0.01   0.95 -1.39 -0.21  0.44  1.09   2.31 10291    1
#> z[19,3]                0.38    0.01   0.94 -1.46 -0.26  0.38  1.02   2.24 11117    1
#> z[20,1]               -0.18    0.01   0.95 -2.03 -0.81 -0.19  0.44   1.68 11557    1
#> z[20,2]                0.06    0.01   0.93 -1.74 -0.58  0.07  0.70   1.85 11933    1
#> z[20,3]                0.31    0.01   0.93 -1.52 -0.31  0.32  0.92   2.13 10493    1
#> z[21,1]               -0.06    0.01   0.98 -1.99 -0.72 -0.06  0.57   1.87 10224    1
#> z[21,2]               -0.45    0.01   0.93 -2.25 -1.08 -0.45  0.19   1.38 11183    1
#> z[21,3]                0.37    0.01   0.94 -1.46 -0.27  0.37  0.99   2.25 10855    1
#> z[22,1]               -0.08    0.01   0.94 -1.90 -0.72 -0.08  0.54   1.76 11004    1
#> z[22,2]               -0.31    0.01   0.92 -2.10 -0.94 -0.33  0.33   1.42  9212    1
#> z[22,3]                0.52    0.01   0.95 -1.36 -0.12  0.52  1.15   2.35 10928    1
#> z[23,1]               -0.01    0.01   0.95 -1.86 -0.65 -0.01  0.63   1.87 10252    1
#> z[23,2]               -0.33    0.01   0.94 -2.19 -0.97 -0.32  0.31   1.52 11729    1
#> z[23,3]               -0.19    0.01   0.95 -2.04 -0.84 -0.19  0.45   1.69 11441    1
#> z[24,1]               -0.20    0.01   0.94 -2.06 -0.86 -0.19  0.44   1.61 10254    1
#> z[24,2]               -0.25    0.01   0.96 -2.12 -0.88 -0.25  0.40   1.60 10151    1
#> z[24,3]               -0.07    0.01   0.93 -1.90 -0.70 -0.08  0.57   1.76  9536    1
#> z[25,1]               -0.17    0.01   0.96 -2.04 -0.83 -0.18  0.47   1.75 10135    1
#> z[25,2]               -0.27    0.01   0.93 -2.06 -0.91 -0.27  0.40   1.52 10356    1
#> z[25,3]                0.07    0.01   0.95 -1.80 -0.57  0.07  0.73   1.97 10269    1
#> z[26,1]               -0.06    0.01   0.94 -1.91 -0.70 -0.05  0.58   1.80  9759    1
#> z[26,2]               -0.13    0.01   0.95 -1.99 -0.78 -0.12  0.51   1.73 10365    1
#> z[26,3]                0.21    0.01   0.93 -1.66 -0.42  0.21  0.83   2.03  9971    1
#> z[27,1]               -0.08    0.01   0.99 -1.99 -0.75 -0.09  0.59   1.89 11061    1
#>  [ reached 'max' / getOption("max.print") -- omitted 3759 rows ]
#> 
#> Samples were drawn using NUTS(diag_e) at Thu Jul  9 04:47:01 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

We can forecast

``` r

forecast(bvar_obj, ci = 0.68, show_all = TRUE)
#> Error in `forecast()`:
#> ! unused argument (ci = 0.68)
```

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
Economic Statistics*, 29(3), pp. 327-341.

Koop, G. and Korobilis, D. (2010). Bayesian multivariate time series
methods for empirical macroeconomics. *Foundations and Trends in
Econometrics*, 3(4), pp. 267-358.

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics*, 24(4), pp. 630-650.
