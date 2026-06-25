test_that("priors requires bvar object", {
  expect_error(priors(list()), "must be a 'bvar' object")
})

test_that("priors requires setup to be run first", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  
  expect_error(priors(model), "must be passed through setup")
})

test_that("priors validates lambda parameters are positive", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- setup(model, p = 2, deterministic = "constant")
  
  expect_error(priors(model, lambda_1 = -0.1), "must be positive")
  expect_error(priors(model, lambda_2 = 0), "must be positive")
  expect_error(priors(model, lambda_3 = -1), "must be positive")
})

test_that("priors validates first_own_lag_prior_mean length", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- setup(model, p = 2, deterministic = "constant")
  
  expect_error(
    priors(model, first_own_lag_prior_mean = c(1,2)),
    "must have length k"
  )
})

test_that("priors validates theta_Psi and Omega_Psi dimensions match", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- setup(model, p = 2, deterministic = "constant")
  
  expect_error(
    priors(model, theta_Psi = c(1,2,3), Omega_Psi = diag(2)),
    "dimensions must match"
  )
})

test_that("priors default behavior", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- setup(model, p = 2, deterministic = "constant")
  
  result <- priors(model)
  
  expect_s3_class(result, "bvar")
  expect_true(!is.null(result$priors))
})

test_that("theta_beta structure correct", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- setup(model, p = 2, deterministic = "constant")
  
  result <- priors(model)
  
  expect_equal(length(result$priors$theta_beta), 18)
})

test_that("Omega_beta structure correct", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- setup(model, p = 2, deterministic = "constant")
  
  result <- priors(model)
  
  expect_equal(dim(result$priors$Omega_beta), c(18, 18))
})

test_that("Sigma_AR computed", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- setup(model, p = 2, deterministic = "constant")
  
  result <- priors(model)
  
  expect_equal(dim(result$priors$Sigma_AR), c(3,3))
})

test_that("Jeffrey FALSE adds hyperparameters", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- setup(model, p = 2, deterministic = "constant")
  
  r1 <- priors(model, Jeffrey = TRUE)
  r2 <- priors(model, Jeffrey = FALSE)
  
  expect_null(r1$priors$m_0)
  expect_false(is.null(r2$priors$m_0))
})

test_that("SV requires SV_priors", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- setup(model, p = 2, deterministic = "constant")
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "AR"),
    "requires SV_priors"
  )
})

test_that("SV attaches correctly", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- setup(model, p = 2, deterministic = "constant")
  
  sv <- list(test = 1)
  
  result <- priors(model, SV = TRUE, SV_type = "AR", SV_priors = sv)
  
  expect_equal(result$priors$SV_priors, sv)
  expect_true(result$priors$SV)
})