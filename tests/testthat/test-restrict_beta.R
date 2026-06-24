test_that("restrict_beta validates dimensions", {
  # Create a minimal bvar object with setup info
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$priors <- list(Omega_beta = diag(8))
  
  # Correct dimensions: k*p x k = 4 x 2
  R_correct <- matrix(1, nrow = 4, ncol = 2)
  result <- restrict_beta(model, R_correct)
  expect_s3_class(result, "bvar")
})

test_that("restrict_beta rejects wrong dimensions", {
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$priors <- list(Omega_beta = diag(8))
  
  # Wrong dimensions
  R_wrong <- matrix(1, nrow = 3, ncol = 2)
  expect_error(restrict_beta(model, R_wrong), 
               "restriction_matrix must have dimension")
})

test_that("restrict_beta stores restriction matrix", {
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$priors <- list(Omega_beta = diag(8))
  
  R <- matrix(c(1, 1, 0, 1, 1, 1, 1, 1), nrow = 4, ncol = 2)
  result <- restrict_beta(model, R)
  expect_equal(result$setup$restriction_matrix, R)
})

test_that("restrict_beta updates Omega_beta diagonal", {
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  omega_orig <- diag(rep(1, 8))
  model$priors <- list(Omega_beta = omega_orig)
  
  R <- matrix(1, nrow = 4, ncol = 2)
  R[1, 1] <- 0  # Restrict first element
  
  result <- restrict_beta(model, R)
  
  # Check that restricted element is very small
  expect_lt(diag(result$priors$Omega_beta)[1], 0.0001)
  # Check other elements unchanged
  expect_equal(diag(result$priors$Omega_beta)[2], 1)
})

test_that("restrict_beta handles missing Omega_beta", {
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$priors <- list()  # No Omega_beta
  
  R <- matrix(1, nrow = 4, ncol = 2)
  expect_warning(restrict_beta(model, R), 
                 "Omega_beta not found")
})

test_that("restrict_beta multiple restrictions", {
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$priors <- list(Omega_beta = diag(8))
  
  R <- matrix(1, nrow = 4, ncol = 2)
  R[c(1, 3), 1] <- 0  # Restrict multiple elements
  
  result <- restrict_beta(model, R)
  
  # Check both restricted elements are small
  expect_lt(diag(result$priors$Omega_beta)[1], 0.0001)
  expect_lt(diag(result$priors$Omega_beta)[3], 0.0001)
})

test_that("restrict_beta returns bvar object", {
  model <- bvar()
  model$setup <- list(k = 2, p = 2)
  model$priors <- list(Omega_beta = diag(8))
  
  R <- matrix(1, nrow = 4, ncol = 2)
  result <- restrict_beta(model, R)
  
  expect_s3_class(result, "bvar")
  expect_type(result, "list")
})
