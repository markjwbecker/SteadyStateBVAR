## R CMD check results

0 errors | 0 warnings | 1 note

* Days since last update: 4

The maintainers of StanHeaders and rstan are preparing new submissions to CRAN
and my package version 0.1.0 would fail to build (reverse-dependency failure) due to the
use of old array syntax which is no longer supported. This update provides the new array syntax
which is compatible with both the current CRAN rstan/StanHeaders and the upcoming releases.
See https://github.com/markjwbecker/SteadyStateBVAR/pull/4 and https://github.com/stan-dev/rstantools/pull/154.

## revdepcheck results

There are currently no downstream dependencies for this package.