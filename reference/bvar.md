# Create a steady-state BVAR model object

Initialises a A steady-state `bvar` object. This is the starting point
for all models in `SteadyStateBVAR`. After creation, pass the object
sequentially to
[`setup`](https://markjwbecker.github.io/SteadyStateBVAR/reference/setup.md),
[`priors`](https://markjwbecker.github.io/SteadyStateBVAR/reference/priors.md),
and
[`fit`](https://markjwbecker.github.io/SteadyStateBVAR/reference/fit.md)
to build and estimate the model.

## Usage

``` r
bvar(data)
```

## Arguments

- data:

  A numeric matrix or time series of data where each column is a
  variable and each row is a time period.

## Value

An object of class `bvar`.

## Examples

``` r
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)
```
