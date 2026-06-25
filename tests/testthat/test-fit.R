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

test_that("fit requires predict", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  expect_error(fit(model), "predict must contain H and d_pred")
})

test_that("fit validates H", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  model$predict <- list(H = 0, d_pred = matrix(0, 1, 1))
  expect_error(fit(model), "H must be a positive integer")
})

test_that("fit validates d_pred type", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  model$predict <- list(H = 1, d_pred = 1)
  expect_error(fit(model), "d_pred must be a matrix")
})

test_that("fit validates d_pred row dimension", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  model$predict <- list(H = 2, d_pred = matrix(0, 1, 3))
  expect_error(fit(model), "nrow\\(d_pred\\) must equal H")
})

test_that("fit validates d_pred column dimension", {
  model <- bvar(data = matrix(rnorm(300), nrow = 100, ncol = 3))
  model <- setup(model, p = 2)
  model <- priors(model)
  
  model$predict <- list(H = 1, d_pred = matrix(0, 1, 999))
  expect_error(fit(model), "ncol\\(d_pred\\) must equal q")
})