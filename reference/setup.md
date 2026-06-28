# Set up the steady-state BVAR model

Prepares the matrices needed for prior specification and estimation.
Also computes OLS estimates.

## Usage

``` r
setup(
  x,
  p,
  deterministic = c("constant", "constant_and_dummy", "constant_and_trend"),
  dummy = NULL
)
```

## Arguments

- x:

  A steady-state `bvar` object created by
  [`bvar`](https://markjwbecker.github.io/SteadyStateBVAR/reference/bvar.md).

- p:

  Integer. The lag order of the VAR.

- deterministic:

  Character. The deterministic component to include. One of `"constant"`
  (default), `"constant_and_dummy"`, or `"constant_and_trend"`.

- dummy:

  Numeric vector of a dummy variable. Only used when
  `deterministic = "constant_and_dummy"`. Default `NULL`.

## Value

The `bvar` object with a `setup` list containing the matrices required
for prior specification and estimation, and also the OLS estimates.

## Examples

``` r
yt <- matrix(rnorm(50), 25, 2)

bvar_obj <- bvar(data = yt)

bvar_obj <- setup(bvar_obj, p = 1)
```
