# Prior Probability Interval for a Normal Distribution

Calculates the mean and variance of a normal prior probability interval.
Given a lower and upper bound of a prior probability interval this
function recovers the implied normal prior parameters. Useful for
specifying informative priors on the steady-state parameters (elements
of \\\Psi\\).

## Usage

``` r
ppi(l, u, interval = 0.95, annualized_growthrate = FALSE, freq = 4)
```

## Arguments

- l:

  Numeric. The lower bound of the prior probability interval.

- u:

  Numeric. The upper bound of the prior probability interval.

- interval:

  Numeric. The prior probability mass within the interval. Default
  `0.95`, i.e. 95%.

- annualized_growthrate:

  Logical. If `TRUE`, treats `l` and `u` as bounds on the annualized
  steady-state growth rate and calculates the implied mean and variance
  on the corresponding quarterly/monthly scale. Useful if you are
  working with a variable `diff(log(x))` (which may be scaled by 100).
  Default `FALSE`.

- freq:

  Integer. The data frequency (e.g. `4` for quarterly). Only used when
  `annualized_growthrate = TRUE`. Default `4`.

## Value

A list with two elements: `mean` and `var`, giving the mean and variance
of the implied normal distribution.

## Details

Consider a CPI variable `CPI <- data$CPI`, observed at quarterly
frequency. In the model, quarter-on-quarter inflation is used
`x <- 100*diff(log(CPI))`. Suppose the prior belief is that annualized
steady-state inflation lies between 1.7 and 2.3 with 95% probability
(mean 2). On the quarterly scale used by `x`, this corresponds to a 95%
interval of 0.425 to 0.575 (mean 0.5). Since it is typically more
natural to elicit a prior on the annualized scale, the
`annualized_growthrate` argument performs this conversion. See
[`vignette("Homoscedastic-steady-state-BVAR")`](https://markjwbecker.github.io/SteadyStateBVAR/articles/Homoscedastic-steady-state-BVAR.md)
for usage in practice.

## Examples

``` r
ppi(l = 1.7, u = 2.3, interval = 0.95)
#> $mean
#> [1] 2
#> 
#> $var
#> [1] 0.0234286
#> 
ppi(l = 1.7, u = 2.3, interval = 0.95, annualized_growthrate = TRUE, freq = 4)
#> $mean
#> [1] 0.5
#> 
#> $var
#> [1] 0.001464287
#> 
```
