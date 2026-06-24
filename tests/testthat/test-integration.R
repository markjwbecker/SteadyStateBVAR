test_that("basic workflow structure", {
  # Create sample data
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  colnames(data) <- c("Y1", "Y2", "Y3")
  
  # Create bvar object
  model <- bvar(data = data)
  
  expect_s3_class(model, "bvar")
  expect_equal(nrow(model$data), 100)
  expect_equal(ncol(model$data), 3)
})

test_that("ppi can be used for prior specification", {
  # Test realistic scenario: researcher specifies inflation prior
  inflation_prior <- ppi(l = 1.5, u = 3.5, interval = 0.90)
  
  expect_equal(inflation_prior$mean, 2.5)
  expect_gt(inflation_prior$var, 0)
})

test_that("multiple ppi calls with different intervals", {
  priors_list <- list(
    inflation = ppi(l = 1, u = 3),
    growth = ppi(l = 2, u = 4),
    interest = ppi(l = 0.5, u = 2)
  )
  
  expect_length(priors_list, 3)
  expect_true(all(sapply(priors_list, function(x) !is.null(x$mean))))
  expect_true(all(sapply(priors_list, function(x) !is.null(x$var))))
})

test_that("bvar object can chain operations", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  
  # Create and extend object
  model <- bvar(data = data)
  expect_s3_class(model, "bvar")
  
  # Add setup info
  model$setup <- list(k = 3, p = 2)
  expect_equal(model$setup$k, 3)
  
  # Add priors info  
  model$priors <- list(Omega_beta = diag(6))
  expect_equal(dim(model$priors$Omega_beta), c(6, 6))
})

test_that("restriction workflow", {
  # Setup
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$priors <- list(Omega_beta = diag(4))
  
  # Apply restriction
  R <- matrix(c(1, 0, 1, 1, 1, 1, 1, 1), nrow = 4, ncol = 2)
  result <- restrict_beta(model, R)
  
  # Verify restriction applied
  expect_true(!is.null(result$setup$restriction_matrix))
  expect_equal(result$setup$restriction_matrix[1, 1], 0)
})
