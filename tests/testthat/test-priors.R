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
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  expect_error(priors(model, lambda_1 = -0.1), "must be positive")
  expect_error(priors(model, lambda_2 = 0), "must be positive")
  expect_error(priors(model, lambda_3 = -1), "must be positive")
})

test_that("priors validates first_own_lag_prior_mean length", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  expect_error(
    priors(model, first_own_lag_prior_mean = c(1,2)),
    "must have length k"
  )
})

test_that("priors validates theta_Psi and Omega_Psi dimensions match", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  expect_error(
    priors(model, theta_Psi = c(1,2,3), Omega_Psi = diag(2)),
    "dimensions must match"
  )
})

test_that("priors default behavior", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  result <- priors(model)
  
  expect_s3_class(result, "bvar")
  expect_true(!is.null(result$priors))
})

test_that("theta_beta structure correct", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  result <- priors(model)
  
  expect_equal(length(result$priors$theta_beta), 18)
})

test_that("Omega_beta structure correct", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  result <- priors(model)
  
  expect_equal(dim(result$priors$Omega_beta), c(18, 18))
})

test_that("Sigma_AR computed", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  result <- priors(model)
  
  expect_equal(dim(result$priors$Sigma_AR), c(3,3))
})

test_that("Jeffreys FALSE adds hyperparameters", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  r1 <- priors(model, Jeffreys = TRUE)
  r2 <- priors(model, Jeffreys = FALSE)
  
  expect_null(r1$priors$m)
  expect_false(is.null(r2$priors$m))
})

test_that("Jeffreys FALSE hyperparameters have expected values", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  result <- priors(model, Jeffreys = FALSE)
  
  expect_equal(result$priors$m, k + 2)
  expect_equal(result$priors$V, (result$priors$m - k - 1) * model$setup$Sigma_u_OLS)
})

test_that("SV requires SV_priors", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "AR1"),
    "SV_priors is required"
  )
})

test_that("RW SV_priors validates required names", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A = rep(0, n_free),
    Omega_A = diag(n_free)
    # missing: mu_log_lambda_1, sigma2_log_lambda_1, alpha_phi, beta_phi
  )
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "RW", SV_priors = sv),
    "SV_priors is missing elements"
  )
})

test_that("RW SV_priors validates theta_A length", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A             = rep(0, n_free + 1), # wrong length
    Omega_A             = diag(n_free),
    mu_log_lambda_1     = rep(0, k),
    sigma2_log_lambda_1 = rep(1, k),
    alpha_phi           = rep(1, k),
    beta_phi            = rep(1, k)
  )
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "RW", SV_priors = sv),
    "theta_A must be a vector of length"
  )
})

test_that("RW SV_priors validates Omega_A dimensions", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A             = rep(0, n_free),
    Omega_A             = diag(n_free + 1), # wrong dimensions
    mu_log_lambda_1     = rep(0, k),
    sigma2_log_lambda_1 = rep(1, k),
    alpha_phi           = rep(1, k),
    beta_phi            = rep(1, k)
  )
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "RW", SV_priors = sv),
    "Omega_A must be a"
  )
})

test_that("RW SV_priors validates sigma2_log_lambda_1 positivity", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A             = rep(0, n_free),
    Omega_A             = diag(n_free),
    mu_log_lambda_1     = rep(0, k),
    sigma2_log_lambda_1 = rep(-1, k), # must be positive
    alpha_phi           = rep(1, k),
    beta_phi            = rep(1, k)
  )
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "RW", SV_priors = sv),
    "sigma2_log_lambda_1 must be strictly positive"
  )
})

test_that("RW SV_priors validates alpha_phi positivity", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A             = rep(0, n_free),
    Omega_A             = diag(n_free),
    mu_log_lambda_1     = rep(0, k),
    sigma2_log_lambda_1 = rep(1, k),
    alpha_phi           = rep(-1, k), # must be positive
    beta_phi            = rep(1, k)
  )
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "RW", SV_priors = sv),
    "alpha_phi must be strictly positive"
  )
})

test_that("RW SV_priors validates beta_phi positivity", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A             = rep(0, n_free),
    Omega_A             = diag(n_free),
    mu_log_lambda_1     = rep(0, k),
    sigma2_log_lambda_1 = rep(1, k),
    alpha_phi           = rep(1, k),
    beta_phi            = rep(-1, k) # must be positive
  )
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "RW", SV_priors = sv),
    "beta_phi must be strictly positive"
  )
})

test_that("RW SV attaches correctly with valid priors", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A             = rep(0, n_free),
    Omega_A             = diag(n_free),
    mu_log_lambda_1     = rep(0, k),
    sigma2_log_lambda_1 = rep(1, k),
    alpha_phi           = rep(5, k),
    beta_phi            = rep(0.4, k)
  )
  
  result <- priors(model, SV = TRUE, SV_type = "RW", SV_priors = sv)
  
  expect_true(result$priors$SV)
  expect_equal(result$priors$SV_type, "RW")
  expect_equal(result$priors$SV_priors, sv)
})

test_that("AR1 SV_priors validates required names", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A = rep(0, n_free),
    Omega_A = diag(n_free)
    # missing remaining AR1 elements, including V_Phi and m_Phi
  )
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "AR1", SV_priors = sv),
    "SV_priors is missing elements"
  )
})

test_that("AR1 SV_priors validates Omega_gamma_0 dimensions", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A            = rep(0, n_free),
    Omega_A            = diag(n_free),
    theta_gamma_0      = rep(0, k),
    Omega_gamma_0      = diag(k + 1), # wrong dimensions
    theta_gamma_1      = rep(0, k),
    Omega_gamma_1      = diag(k),
    theta_log_lambda_1 = rep(0, k),
    Omega_log_lambda_1 = diag(k),
    m_Phi              = k + 2,
    V_Phi              = diag(k)
  )
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "AR1", SV_priors = sv),
    "Omega_gamma_0 must be a"
  )
})

test_that("AR1 SV_priors validates m_Phi >= k", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A            = rep(0, n_free),
    Omega_A            = diag(n_free),
    theta_gamma_0      = rep(0, k),
    Omega_gamma_0      = diag(k),
    theta_gamma_1      = rep(0, k),
    Omega_gamma_1      = diag(k),
    theta_log_lambda_1 = rep(0, k),
    Omega_log_lambda_1 = diag(k),
    m_Phi              = k - 1, # too small
    V_Phi              = diag(k)
  )
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "AR1", SV_priors = sv),
    "m_Phi must be a scalar integer >= k"
  )
})

test_that("AR1 SV_priors validates V_Phi dimensions", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A            = rep(0, n_free),
    Omega_A            = diag(n_free),
    theta_gamma_0      = rep(0, k),
    Omega_gamma_0      = diag(k),
    theta_gamma_1      = rep(0, k),
    Omega_gamma_1      = diag(k),
    theta_log_lambda_1 = rep(0, k),
    Omega_log_lambda_1 = diag(k),
    m_Phi              = k + 2,
    V_Phi              = diag(k + 1) # wrong dimensions
  )
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "AR1", SV_priors = sv),
    "V_Phi must be a"
  )
})

test_that("AR1 SV attaches correctly with valid priors", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  k <- model$setup$k
  n_free <- model$setup$n_free_params_A
  
  sv <- list(
    theta_A            = rep(0, n_free),
    Omega_A            = diag(n_free),
    theta_gamma_0      = rep(0, k),
    Omega_gamma_0      = diag(k),
    theta_gamma_1      = rep(0.9, k),
    Omega_gamma_1      = diag(k),
    theta_log_lambda_1 = rep(0, k),
    Omega_log_lambda_1 = diag(k),
    m_Phi              = k + 2,
    V_Phi              = diag(k)
  )
  
  result <- priors(model, SV = TRUE, SV_type = "AR1", SV_priors = sv)
  
  expect_true(result$priors$SV)
  expect_equal(result$priors$SV_type, "AR1")
  expect_equal(result$priors$SV_priors, sv)
})

test_that("SV_type must be RW or AR1", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data = data)
  model <- SteadyStateBVAR::setup(model, p = 2, deterministic = "constant")
  
  expect_error(
    priors(model, SV = TRUE, SV_type = "invalid", SV_priors = list()),
    "SV_type needs to be RW or AR1"
  )
})