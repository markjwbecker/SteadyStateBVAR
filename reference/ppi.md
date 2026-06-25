# Prior Probability Interval for a Normal Distribution

Converts a symmetric prior probability interval into the corresponding
mean and variance of a normal distribution. Given a lower and upper
bound that define a prior probability interval, this function recovers
the implied normal prior parameters. Useful for specifying informative
priors on the steady-state (Psi) parameters in an intuitive way.

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
  `0.95`.

- annualized_growthrate:

  Logical. If `TRUE`, converts the interval from annualized to
  per-period units by dividing by `freq`. Default `FALSE`.

- freq:

  Integer. The data frequency (e.g. `4` for quarterly). Only used when
  `annualized_growthrate = TRUE`. Default `4`.

## Value

A list with two elements: `mean` and `var`, giving the mean and variance
of the implied normal distribution.

## Examples

``` r
# steady-state annualized inflation is between 1% and 3% with 95% probability
ppi(l = 1, u = 3, interval = 0.95)
#> $mean
#> [1] 2
#> 
#> $var
#> [1] 0.2603178
#> 
```
