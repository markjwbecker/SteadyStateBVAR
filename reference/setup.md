# Set up a BVAR model

A generic function for setting up a BVAR model. Computes OLS estimates
and prepares all matrices needed for prior specification and estimation.

## Usage

``` r
setup(x, ...)

# S3 method for class 'bvar'
setup(
  x,
  p,
  deterministic = c("constant", "constant_and_dummy", "constant_and_trend"),
  dummy = NULL,
  ...
)
```

## Arguments

- x:

  A `bvar` object created by
  [`bvar`](https://markjwbecker.github.io/SteadyStateBVAR/reference/bvar.md).

- ...:

  Further arguments passed to methods.

- p:

  Integer. The lag order of the VAR.

- deterministic:

  Character. The deterministic component to include. One of `"constant"`
  (default), `"constant_and_dummy"`, or `"constant_and_trend"`.

- dummy:

  Optional numeric vector or matrix of dummy variables. Only used when
  `deterministic = "constant_and_dummy"`. Default `NULL`.

## Value

The `bvar` object with a `setup` list appended.

## Examples

``` r
yt <- matrix(rnorm(40, 0, 1), 20, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p = 1)
```
