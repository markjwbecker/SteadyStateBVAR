test_that("restrict_beta stores restriction matrix", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data)
  model <- SteadyStateBVAR::setup(model, p = 2)
  model <- priors(model)
  
  k <- model$setup$k
  p <- model$setup$p
  R <- matrix(1, nrow = k * p, ncol = k)
  R[2, 1] <- 0
  
  result <- restrict_beta(model, R)
  
  expect_equal(result$setup$restriction_matrix, R)
})

test_that("restrict_beta updates Omega_beta for zero restrictions", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data)
  model <- SteadyStateBVAR::setup(model, p = 2)
  model <- priors(model)
  
  k <- model$setup$k
  p <- model$setup$p
  R <- matrix(1, nrow = k * p, ncol = k)
  R[2, 1] <- 0
  
  zero_idx <- which(c(R) == 0)
  result <- restrict_beta(model, R)
  
  expect_equal(diag(result$priors$Omega_beta)[zero_idx], 1e-7)
})

test_that("restrict_beta checks matrix dimensions", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data)
  model <- SteadyStateBVAR::setup(model, p = 2)
  model <- priors(model)
  
  k <- model$setup$k
  p <- model$setup$p
  R_bad <- matrix(1, nrow = k * p - 1, ncol = k)
  
  expect_error(
    restrict_beta(model, R_bad),
    "restriction_matrix must have dimension"
  )
})

test_that("restrict_beta warns if Omega_beta is missing", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data)
  model <- SteadyStateBVAR::setup(model, p = 2)
  model$priors <- list()
  
  k <- model$setup$k
  p <- model$setup$p
  R <- matrix(1, nrow = k * p, ncol = k)
  
  expect_warning(
    restrict_beta(model, R),
    "Omega_beta not found"
  )
})
