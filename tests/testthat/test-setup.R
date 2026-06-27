test_that("setup with constant deterministic", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  result <- setup(model, p = 2, deterministic = "constant")
  
  expect_s3_class(result, "bvar")
  expect_true(!is.null(result$setup))
  expect_equal(result$setup$k, 3)
  expect_equal(result$setup$p, 2)
  expect_equal(result$setup$q, 1)
})

test_that("setup with constant_and_trend", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  result <- setup(model, p = 2, deterministic = "constant_and_trend")
  
  expect_equal(result$setup$q, 2)
  expect_equal(ncol(result$setup$dt), 2)  # constant + trend
})

test_that("setup with constant_and_dummy", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  dummy_var <- c(rep(0, 50), rep(1, 50))
  model <- bvar(data = data)
  result <- setup(model, p = 2, deterministic = "constant_and_dummy", dummy = dummy_var)
  
  expect_equal(result$setup$q, 2)
  expect_equal(result$setup$dummy, dummy_var)
})

test_that("setup computes correct matrix dimensions", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  result <- setup(model, p = 2, deterministic = "constant")
  
  # N = nrow - p = 100 - 2 = 98
  expect_equal(result$setup$N, 98)
  expect_equal(nrow(result$setup$Y), 98)
  expect_equal(nrow(result$setup$X), 98)
  expect_equal(nrow(result$setup$W), 98)
})

test_that("setup computes OLS estimates", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  result <- setup(model, p = 2, deterministic = "constant")
  
  expect_true(!is.null(result$setup$beta_OLS))
  expect_true(!is.null(result$setup$Sigma_u_OLS))
  expect_true(!is.null(result$setup$Psi_OLS))
  expect_equal(nrow(result$setup$Sigma_u_OLS), 3)
})

test_that("setup with different lag orders", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  
  result_p1 <- setup(model, p = 1)
  result_p3 <- setup(model, p = 3)
  
  expect_equal(result_p1$setup$p, 1)
  expect_equal(result_p3$setup$p, 3)
  expect_equal(result_p1$setup$N, 99)  # 100 - 1
  expect_equal(result_p3$setup$N, 97)  # 100 - 3
})

test_that("setup returns updated bvar object", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  result <- setup(model, p = 2)
  
  expect_equal(result$data, model$data)  # data unchanged
  expect_true(!is.null(result$setup))
})

test_that("setup validates p is positive integer", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  
  expect_error(setup(model, p = -1), "p must be a positive integer")
  expect_error(setup(model, p = 0), "p must be a positive integer")
  expect_error(setup(model, p = 2.5), "p must be a positive integer")
})

test_that("setup validates p is less than nrow", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  
  expect_error(setup(model, p = 100), "p must be less than the number of observations")
})

test_that("setup requires data in bvar object", {
  model <- structure(
    list(
      data = NULL,
      setup = NULL,
      priors = NULL,
      fit = NULL,
      predict = list()
    ),
    class = "bvar"
  )
  
  expect_error(
    setup(model, p = 2),
    "x must have data"
  )
})

test_that("setup requires dummy when constant_and_dummy", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  
  expect_error(setup(model, p = 2, deterministic = "constant_and_dummy"),
               "dummy must be provided")
})

test_that("setup validates dummy length", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  dummy_wrong <- c(rep(0, 50))  # Should be 100
  model <- bvar(data = data)
  
  expect_error(setup(model, p = 2, deterministic = "constant_and_dummy", dummy = dummy_wrong),
               "dummy must have length equal to number of observations")
})