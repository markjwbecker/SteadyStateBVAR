test_that("bvar creates object with correct class", {
  model <- bvar()
  expect_s3_class(model, "bvar")
})

test_that("bvar initializes with NULL defaults", {
  model <- bvar()
  expect_null(model$data)
  expect_null(model$setup)
  expect_null(model$priors)
  expect_null(model$fit)
  expect_type(model$predict, "list")
  expect_length(model$predict, 0)
})

test_that("bvar accepts data argument", {
  test_data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  model <- bvar(data = test_data)
  expect_equal(model$data, test_data)
})

test_that("bvar accepts predict list", {
  predict_list <- list(H = 8, d_pred = matrix(1, 8, 1))
  model <- bvar(predict = predict_list)
  expect_equal(model$predict, predict_list)
})

test_that("bvar stores all arguments", {
  test_data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  predict_list <- list(H = 8)
  model <- bvar(data = test_data, predict = predict_list)
  expect_equal(model$data, test_data)
  expect_equal(model$predict, predict_list)
})

test_that("bvar creates proper list structure", {
  model <- bvar()
  expect_type(model, "list")
  expect_named(model, c("data", "setup", "priors", "fit", "predict"))
})

test_that("bvar with realistic data structure", {
  # Create sample time series data: 100 observations, 3 variables
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  colnames(data) <- c("GDP", "Inflation", "Interest_Rate")
  
  model <- bvar(data = data)
  expect_equal(nrow(model$data), 100)
  expect_equal(ncol(model$data), 3)
})
