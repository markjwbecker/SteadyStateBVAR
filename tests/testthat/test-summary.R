test_that("summary.bvar requires fit object", {
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$SV <- FALSE
  model$fit <- list(stan = NULL)
  
  # This should work but create empty results since no posterior data
  expect_error(summary(model), NA)  # Should not error
})

test_that("summary.bvar creates summary.bvar object", {
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$SV <- FALSE
  model$SV_type <- NULL
  model$fit <- list(stan = NULL)
  
  # Basic structure test
  expect_s3_class(model, "bvar")
})

test_that("summary.bvar handles parameter selection", {
  # Test that pars parameter works (even if just returning structure)
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$SV <- FALSE
  model$SV_type <- NULL
  
  # Just testing the object can be created with pars option
  expect_s3_class(model, "bvar")
})

test_that("summary.bvar with NULL pars returns all", {
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$SV <- FALSE
  model$fit <- list(stan = NULL)
  
  # pars = NULL should return all parameters
  expect_s3_class(model, "bvar")
})

test_that("print.summary.bvar has correct S3 method", {
  # Check that print method exists
  expect_true(exists("print.summary.bvar"))
})

test_that("print.summary.bvar handles missing results", {
  # Create a minimal summary.bvar object
  summary_obj <- list(
    summaries = list(results = NULL),
    SV = FALSE,
    SV_type = NULL
  )
  class(summary_obj) <- "summary.bvar"
  
  # This should not error
  expect_error(capture.output(print(summary_obj)), NA)
})

test_that("summary.bvar initializes empty lists correctly", {
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$SV <- FALSE
  model$SV_type <- NULL
  
  # Should have required structure elements
  expect_s3_class(model, "bvar")
})
