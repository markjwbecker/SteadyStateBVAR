
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SteadyStateBVAR

<!-- badges: start -->
<!-- badges: end -->

The goal of SteadyStateBVAR is to …

## Installation

You can install the development version of SteadyStateBVAR from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("markjwbecker/SteadyStateBVAR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
devtools::load_all()
#> ℹ Loading SteadyStateBVAR
#> Loading required package: rstan
#> 
#> Loading required package: StanHeaders
#> 
#> 
#> rstan version 2.32.6 (Stan version 2.32.2)
#> 
#> 
#> For execution on a local, multicore CPU with excess RAM we recommend calling
#> options(mc.cores = parallel::detectCores()).
#> To avoid recompilation of unchanged Stan programs, we recommend calling
#> rstan_options(auto_write = TRUE)
#> For within-chain threading using `reduce_sum()` or `map_rect()` Stan functions,
#> change `threads_per_chain` option:
#> rstan_options(threads_per_chain = 1)
#> 
#> 
#> Do not specify '-march=native' in 'LOCAL_CPPFLAGS' or a Makevars file
#> 
#> Loading required package: bvartools
#> 
#> Loading required package: coda
#> 
#> 
#> Attaching package: 'coda'
#> 
#> 
#> The following object is masked from 'package:rstan':
#> 
#>     traceplot
#> 
#> 
#> Loading required package: Matrix
#library(SteadyStateBVAR)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
