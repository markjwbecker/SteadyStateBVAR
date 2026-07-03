test_that("bvar requires data argument", {
  expect_error(bvar(), "missing")
})

test_that("bvar validates data input", {
  expect_error(bvar(data = "not_data"), 
               "data must be a matrix or time series object")
})

test_that("bvar creates object with correct class", {
  test_data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  model <- bvar(data = test_data)
  expect_s3_class(model, "bvar")
})

test_that("bvar stores data correctly", {
  test_data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  colnames(test_data) <- c("GDP", "Inflation", "Rate")
  
  model <- bvar(data = test_data)
  expect_equal(model$data, test_data)
  expect_equal(nrow(model$data), 100)
  expect_equal(ncol(model$data), 3)
})

test_that("bvar initializes component slots", {
  test_data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  model <- bvar(data = test_data)
  
  expect_null(model$setup)
  expect_null(model$priors)
  expect_null(model$fit)
  expect_type(model$predict, "list")
})
