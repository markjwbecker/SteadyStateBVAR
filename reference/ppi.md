# Prior Probability Interval for a Normal Distribution

Calculates the mean and variance of a normal prior probability interval.
Given a lower and upper bound of a prior probability interval this
function recovers the implied normal prior parameters. Useful for
specifying informative priors on the steady-state parameters (elements
of Psi).

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

  Logical. If `TRUE`, converts the implied mean and variance from the
  annualized scale to the quarterly/monthly etc. scale by dividing by
  `freq`. Default `FALSE`.

- freq:

  Integer. The data frequency (e.g. `4` for quarterly). Only used when
  `annualized_growthrate = TRUE`. Default `4`.

## Value

A list with two elements: `mean` and `var`, giving the mean and variance
of the implied normal distribution.

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
