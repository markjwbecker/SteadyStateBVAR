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

As an example, say we have a CPI variable `CPI <- data$CPI`, where CPI
is on the quarterly frequency. Then in the model we work with
quarter-on-quarter inflation `x <- 100*diff(log(CPI))`. Lets say our
prior for annualized steady-state inflation of `x` is between 1.7 and
2.3 with 95% probability (with mean at 2). This translates to 0.425 and
0.575 with 95% probability (with mean at 0.5). Clearly it is easier to
think of a steady-state prior on the annualized scale, hence the
`annualized_growthrate` argument. Please see
[`vignette("Homoscedastic-steady-state-BVAR")`](https://markjwbecker.github.io/SteadyStateBVAR/articles/Homoscedastic-steady-state-BVAR.md)
on how to use in practise.

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
