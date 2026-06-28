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

test_that("fit requires H and d_pred", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model),
    "H and d_pred must be supplied"
  )
})

test_that("fit validates H", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model,
        H = 0,
        d_pred = matrix(0, 1, 1),
        cores = 1,
        auto_write = FALSE),
    "H must be a positive integer"
  )
})

test_that("fit validates iter", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model,
        H = 1,
        d_pred = matrix(0, 1, 1),
        iter = 0,
        cores = 1,
        auto_write = FALSE),
    "iter must be a positive integer"
  )
})

test_that("fit validates warmup", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model,
        H = 1,
        d_pred = matrix(0, 1, 1),
        warmup = -1,
        cores = 1,
        auto_write = FALSE),
    "warmup must be a non-negative integer"
  )
})

test_that("fit validates chains", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model,
        H = 1,
        d_pred = matrix(0, 1, 1),
        chains = 0,
        cores = 1,
        auto_write = FALSE),
    "chains must be a positive integer"
  )
})

test_that("fit requires cores", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model,
        H = 1,
        d_pred = matrix(0, 1, 1),
        auto_write = FALSE),
    "Please select how many cores to use"
  )
})

test_that("fit validates cores", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model,
        H = 1,
        d_pred = matrix(0, 1, 1),
        cores = 0,
        auto_write = FALSE),
    "cores must be a positive integer"
  )
})

test_that("fit requires auto_write", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model,
        H = 1,
        d_pred = matrix(0, 1, 1),
        cores = 1),
    "Please select TRUE or FALSE for auto_write"
  )
})

test_that("fit validates auto_write", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(
      model,
      H = 1,
      d_pred = matrix(0, 1, 1),
      cores = 1,
      auto_write = 1
    ),
    "auto_write must be TRUE or FALSE"
  )
})

test_that("fit validates d_pred type", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model,
        H = 1,
        d_pred = 1,
        cores = 1,
        auto_write = FALSE),
    "d_pred must be a matrix"
  )
})

test_that("fit validates d_pred row dimension", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model,
        H = 2,
        d_pred = matrix(0, 1, 3),
        cores = 1,
        auto_write = FALSE),
    "nrow\\(d_pred\\) must equal H"
  )
})

test_that("fit validates d_pred column dimension", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(
    fit(model,
        H = 1,
        d_pred = matrix(0, 1, 999),
        cores = 1,
        auto_write = FALSE),
    "ncol\\(d_pred\\) must equal q"
  )
})