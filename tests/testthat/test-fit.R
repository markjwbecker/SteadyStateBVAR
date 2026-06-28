test_that("fit requires bvar object", {
  expect_error(fit(list()), "must be a 'bvar' object")
})

test_that("fit requires setup", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  expect_error(fit(model), "must be passed through setup")
})

test_that("fit requires priors", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  expect_error(fit(model), "must be passed through priors")
})

test_that("fit requires d_pred", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model),
    "d_pred must be supplied"
  )
})

test_that("fit validates H", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model, H = 0, d_pred = matrix(0, 1, 1)),
    "H must be a positive integer"
  )
})

test_that("fit validates iter", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model, H = 1, d_pred = matrix(0, 1, 1), iter = 0),
    "iter must be a positive integer"
  )
})

test_that("fit validates warmup", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model, H = 1, d_pred = matrix(0, 1, 1), warmup = -1),
    "warmup must be a non-negative integer"
  )
})

test_that("fit validates chains", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model, H = 1, d_pred = matrix(0, 1, 1), chains = 0),
    "chains must be a positive integer"
  )
})

test_that("fit validates cores", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model, H = 1, d_pred = matrix(0, 1, 1), cores = 0),
    "cores must be a positive integer"
  )
})

test_that("fit validates auto_write NA", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model, H = 1, d_pred = matrix(0, 1, 1), auto_write = NA),
    "auto_write must be TRUE or FALSE"
  )
})

test_that("fit validates auto_write type", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model, H = 1, d_pred = matrix(0, 1, 1), auto_write = 1),
    "auto_write must be TRUE or FALSE"
  )
})

test_that("fit validates d_pred type", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model, H = 1, d_pred = 1),
    "d_pred must be a matrix"
  )
})

test_that("fit validates d_pred row dimension", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model, H = 2, d_pred = matrix(0, 1, 3)),
    "nrow\\(d_pred\\) must equal H"
  )
})

test_that("fit validates d_pred column dimension", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model, H = 1, d_pred = matrix(0, 1, 999)),
    "ncol\\(d_pred\\) must equal q"
  )
})