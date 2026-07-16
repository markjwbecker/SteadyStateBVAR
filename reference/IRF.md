# Impulse Response Functions for a fitted steady-state BVAR model

Computes and plots impulse response functions (IRFs) from a fitted
steady-state `bvar` object. Supports both orthogonalized (OIRF) and
generalized (GIRF) impulse responses, with optional conversion to annual
growth rates.

## Usage

``` r
IRF(
  x,
  H = 16,
  response = NULL,
  shock = NULL,
  type = c("median", "mean"),
  method = c("OIRF", "GIRF"),
  ci = 0.95,
  t = NULL,
  growth_rate_idx = NULL
)
```

## Arguments

- x:

  A steady-state `bvar` object that has been passed through
  [`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md).

- H:

  Integer. The forecast horizon for the IRF. Default `16`.

- response:

  Integer. Index of the response variable to plot. If `NULL` (default),
  all responses are plotted.

- shock:

  Integer. Index of the shock variable to plot. If `NULL` (default), all
  shocks are plotted.

- type:

  Character. Whether to use `"median"` or `"mean"` as the point
  estimate. Default `"median"`.

- method:

  Character. The IRF method: `"OIRF"` for orthogonalized or `"GIRF"` for
  generalized impulse responses. Default `"OIRF"`.

- ci:

  Numeric. The credible interval width. Default `0.95`, i.e. 95%.

- t:

  Integer. Time index for the covariance matrix when using stochastic
  volatility models. If `NULL` (default), the last time `t` is used.

- growth_rate_idx:

  Integer vector. Indices of variables for which the impulse response is
  converted from a quarterly or monthly log first difference to an
  annual growth rate response, i.e. \\\ln x\_{t} - \ln x\_{t-f}\\, where
  \\f\\ is the frequency of the data (4 for quarterly, 12 for monthly).
  Only suitable for variables specified as \\\ln x\_{t} - \ln
  x\_{t-1}\\, i.e. `diff(log(x))` or `100*diff(log(x))`. Computed by
  summing up to \\f\\ periods of the impulse response, treating the
  response in periods prior to the shock as zero. Default is `NULL`.

## Value

Invisibly returns a list with three arrays: the point estimate IRF,
`lower`, and `upper` credible bounds, each of dimension `k x k x (H+1)`.

## Examples

``` r
# \donttest{
#homoscedastic with Jeffreys prior
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")

bvar_obj <- priors(bvar_obj,
                   lambda_1 = 0.2,
                   lambda_2 = 0.5,
                   lambda_3 = 1,
                   first_own_lag_prior_mean = rep(1,2),
                   theta_Psi = rep(0, 2),
                   Omega_Psi = diag(0.1, 2, 2),
                   Jeffreys = TRUE,
                   SV = FALSE,
                   SV_type = NULL,
                   SV_priors = NULL)
                   
bvar_obj <- fit(bvar_obj,
                H = 8,
                d_pred = matrix(rep(1,8)),
                iter = 200,
                warmup = 50,
                chains = 1,
                cores = 1)
#> 
#> SAMPLING FOR MODEL 'steady_state_bvar_homoscedastic_jeffreys_prior' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 9.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.93 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 7
#> Chain 1:            adapt_window = 38
#> Chain 1:            term_buffer = 5
#> Chain 1: 
#> Chain 1: Iteration:   1 / 200 [  0%]  (Warmup)
#> Chain 1: Iteration:  20 / 200 [ 10%]  (Warmup)
#> Chain 1: Iteration:  40 / 200 [ 20%]  (Warmup)
#> Chain 1: Iteration:  51 / 200 [ 25%]  (Sampling)
#> Chain 1: Iteration:  70 / 200 [ 35%]  (Sampling)
#> Chain 1: Iteration:  90 / 200 [ 45%]  (Sampling)
#> Chain 1: Iteration: 110 / 200 [ 55%]  (Sampling)
#> Chain 1: Iteration: 130 / 200 [ 65%]  (Sampling)
#> Chain 1: Iteration: 150 / 200 [ 75%]  (Sampling)
#> Chain 1: Iteration: 170 / 200 [ 85%]  (Sampling)
#> Chain 1: Iteration: 190 / 200 [ 95%]  (Sampling)
#> Chain 1: Iteration: 200 / 200 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.034 seconds (Warm-up)
#> Chain 1:                0.084 seconds (Sampling)
#> Chain 1:                0.118 seconds (Total)
#> Chain 1: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
                
(IRF(bvar_obj))

#> $median_irf
#> , , 1
#> 
#>             [,1]     [,2]
#> [1,]  1.46192096 0.000000
#> [2,] -0.00172196 1.206163
#> 
#> , , 2
#> 
#>             [,1]       [,2]
#> [1,]  0.83918684 -0.0319003
#> [2,] -0.01503758  0.8321868
#> 
#> , , 3
#> 
#>             [,1]        [,2]
#> [1,]  0.47550176 -0.03670799
#> [2,] -0.01008044  0.56180269
#> 
#> , , 4
#> 
#>             [,1]        [,2]
#> [1,]  0.27071675 -0.03443883
#> [2,] -0.01099643  0.38494122
#> 
#> , , 5
#> 
#>              [,1]        [,2]
#> [1,]  0.155277920 -0.02834319
#> [2,] -0.009901752  0.26888653
#> 
#> , , 6
#> 
#>              [,1]       [,2]
#> [1,]  0.091996119 -0.0203872
#> [2,] -0.006070829  0.1776566
#> 
#> , , 7
#> 
#>              [,1]       [,2]
#> [1,]  0.056442674 -0.0124966
#> [2,] -0.004587836  0.1196212
#> 
#> , , 8
#> 
#>              [,1]        [,2]
#> [1,]  0.034480438 -0.00690474
#> [2,] -0.003021792  0.07956163
#> 
#> , , 9
#> 
#>              [,1]         [,2]
#> [1,]  0.020989191 -0.003775322
#> [2,] -0.001950373  0.053207054
#> 
#> , , 10
#> 
#>              [,1]         [,2]
#> [1,]  0.012273960 -0.002194483
#> [2,] -0.001219304  0.036465721
#> 
#> , , 11
#> 
#>               [,1]         [,2]
#> [1,]  0.0065319762 -0.001423404
#> [2,] -0.0007462871  0.023823474
#> 
#> , , 12
#> 
#>               [,1]          [,2]
#> [1,]  0.0036778576 -0.0007980241
#> [2,] -0.0004501086  0.0161092088
#> 
#> , , 13
#> 
#>               [,1]          [,2]
#> [1,]  0.0023648404 -0.0004328076
#> [2,] -0.0002503679  0.0108162092
#> 
#> , , 14
#> 
#>               [,1]          [,2]
#> [1,]  1.578874e-03 -0.0002252567
#> [2,] -7.482414e-05  0.0063916802
#> 
#> , , 15
#> 
#>               [,1]          [,2]
#> [1,]  1.194699e-03 -0.0002222521
#> [2,] -3.162116e-05  0.0043391202
#> 
#> , , 16
#> 
#>               [,1]          [,2]
#> [1,]  7.618632e-04 -0.0002309927
#> [2,] -1.965334e-05  0.0029523692
#> 
#> , , 17
#> 
#>               [,1]          [,2]
#> [1,]  4.640262e-04 -0.0001081337
#> [2,] -9.491287e-06  0.0019239194
#> 
#> 
#> $lower
#> , , 1
#> 
#>            [,1]      [,2]
#> [1,]  1.0707622 0.0000000
#> [2,] -0.6686813 0.8014603
#> 
#> , , 2
#> 
#>            [,1]       [,2]
#> [1,]  0.3347423 -0.2876796
#> [2,] -0.5719690  0.2822722
#> 
#> , , 3
#> 
#>            [,1]       [,2]
#> [1,]  0.1068138 -0.3714979
#> [2,] -0.6184775  0.0796729
#> 
#> , , 4
#> 
#>             [,1]        [,2]
#> [1,]  0.03138235 -0.46490868
#> [2,] -0.81237838  0.02723198
#> 
#> , , 5
#> 
#>              [,1]        [,2]
#> [1,]  0.005330795 -0.50633636
#> [2,] -0.874706346  0.00746849
#> 
#> , , 6
#> 
#>              [,1]         [,2]
#> [1,] -0.002504324 -0.533787814
#> [2,] -0.933664145  0.002360275
#> 
#> , , 7
#> 
#>             [,1]         [,2]
#> [1,] -0.02250002 -0.585682115
#> [2,] -0.99625177 -0.005080609
#> 
#> , , 8
#> 
#>             [,1]         [,2]
#> [1,] -0.03996368 -0.651884740
#> [2,] -1.06859996 -0.009358407
#> 
#> , , 9
#> 
#>            [,1]         [,2]
#> [1,] -0.0334267 -0.719012835
#> [2,] -1.1565645 -0.009966814
#> 
#> , , 10
#> 
#>             [,1]        [,2]
#> [1,] -0.04019535 -0.78999129
#> [2,] -1.26611362 -0.01225547
#> 
#> , , 11
#> 
#>             [,1]        [,2]
#> [1,] -0.03212359 -0.86797326
#> [2,] -1.38925249 -0.01066854
#> 
#> , , 12
#> 
#>             [,1]        [,2]
#> [1,] -0.02802113 -0.95644754
#> [2,] -1.50432546 -0.01028149
#> 
#> , , 13
#> 
#>             [,1]         [,2]
#> [1,] -0.02433876 -1.059352215
#> [2,] -1.63983660 -0.008571817
#> 
#> , , 14
#> 
#>             [,1]         [,2]
#> [1,] -0.01867863 -1.181204309
#> [2,] -1.80601991 -0.006480311
#> 
#> , , 15
#> 
#>             [,1]        [,2]
#> [1,] -0.01387714 -1.32725174
#> [2,] -2.04016300 -0.00621083
#> 
#> , , 16
#> 
#>            [,1]         [,2]
#> [1,] -0.0113278 -1.503653743
#> [2,] -2.3063265 -0.005893544
#> 
#> , , 17
#> 
#>              [,1]         [,2]
#> [1,] -0.009134826 -1.717696370
#> [2,] -2.609320346 -0.005213065
#> 
#> 
#> $upper
#> , , 1
#> 
#>          [,1]    [,2]
#> [1,] 1.997478 0.00000
#> [2,] 0.523475 1.69452
#> 
#> , , 2
#> 
#>           [,1]      [,2]
#> [1,] 1.6670168 0.2068235
#> [2,] 0.5244501 1.4922021
#> 
#> , , 3
#> 
#>           [,1]      [,2]
#> [1,] 1.5125693 0.2648089
#> [2,] 0.5608019 1.5737126
#> 
#> , , 4
#> 
#>           [,1]      [,2]
#> [1,] 1.4012647 0.3093817
#> [2,] 0.5783903 1.5699856
#> 
#> , , 5
#> 
#>           [,1]     [,2]
#> [1,] 1.3227968 0.360725
#> [2,] 0.5720551 1.621630
#> 
#> , , 6
#> 
#>           [,1]      [,2]
#> [1,] 1.2591571 0.3723103
#> [2,] 0.5188606 1.7512005
#> 
#> , , 7
#> 
#>           [,1]      [,2]
#> [1,] 1.1929721 0.3717744
#> [2,] 0.5102253 1.9163306
#> 
#> , , 8
#> 
#>           [,1]      [,2]
#> [1,] 1.1125797 0.3986562
#> [2,] 0.4996484 2.1733524
#> 
#> , , 9
#> 
#>           [,1]      [,2]
#> [1,] 1.0384364 0.4060361
#> [2,] 0.4878854 2.5230883
#> 
#> , , 10
#> 
#>           [,1]      [,2]
#> [1,] 0.9701940 0.4872215
#> [2,] 0.4754788 2.9539133
#> 
#> , , 11
#> 
#>           [,1]     [,2]
#> [1,] 0.9074863 0.521348
#> [2,] 0.4542901 3.481796
#> 
#> , , 12
#> 
#>           [,1]      [,2]
#> [1,] 0.9136274 0.5353846
#> [2,] 0.4249452 4.1268625
#> 
#> , , 13
#> 
#>           [,1]      [,2]
#> [1,] 1.0292011 0.5275165
#> [2,] 0.4172975 4.9143100
#> 
#> , , 14
#> 
#>           [,1]      [,2]
#> [1,] 1.1701312 0.5201342
#> [2,] 0.4060617 5.8755763
#> 
#> , , 15
#> 
#>           [,1]      [,2]
#> [1,] 1.2553516 0.5135174
#> [2,] 0.3875314 7.0498202
#> 
#> , , 16
#> 
#>           [,1]      [,2]
#> [1,] 1.3248850 0.5025517
#> [2,] 0.3861735 8.4857904
#> 
#> , , 17
#> 
#>           [,1]       [,2]
#> [1,] 1.4115900  0.4844836
#> [2,] 0.3919669 10.2441792
#> 
#> 
# }
```
