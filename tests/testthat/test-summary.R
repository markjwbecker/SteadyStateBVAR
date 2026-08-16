test_that("summary dimensions are correct for homoscedastic Jeffreys model", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data)
  model <- SteadyStateBVAR::setup(model, p = 2)
  model <- priors(model, Jeffreys = TRUE)
  
  k <- model$setup$k
  p <- model$setup$p
  q <- model$setup$q
  
  model <- suppressWarnings(
    fit(model, H = 1, iter = 40, warmup = 20, chains = 1, cores = 1)
  )
  
  posterior <- rstan::extract(model$fit$stan)
  
  expect_equal(dim(model$fit$posterior_means$beta), c(k * p, k))
  expect_equal(dim(model$fit$posterior_means$Psi), c(k, q))
  expect_equal(dim(model$fit$posterior_means$Sigma_u), c(k, k))
  
  expect_null(posterior$A)
  expect_null(posterior$phi)
  expect_null(posterior$gamma_0)
  expect_null(posterior$gamma_1)
})

test_that("summary dimensions are correct for homoscedastic inverse-Wishart model", {
  data <- matrix(rnorm(300), nrow = 100, ncol = 3)
  model <- bvar(data)
  model <- SteadyStateBVAR::setup(model, p = 2)
  model <- priors(model, Jeffreys = FALSE)
  
  k <- model$setup$k
  p <- model$setup$p
  q <- model$setup$q
  
  model <- suppressWarnings(
    fit(model, H = 1, iter = 40, warmup = 20, chains = 1, cores = 1)
  )
  
  expect_equal(dim(model$fit$posterior_means$beta), c(k * p, k))
  expect_equal(dim(model$fit$posterior_means$Psi), c(k, q))
  expect_equal(dim(model$fit$posterior_means$Sigma_u), c(k, k))
})

test_that("summary dimensions are correct for RW stochastic volatility model", {
  data <- matrix(rnorm(200), nrow = 100, ncol = 2)
  model <- bvar(data)
  model <- SteadyStateBVAR::setup(model, p = 1, deterministic = "constant")
  
  k <- model$setup$k
  p <- model$setup$p
  q <- model$setup$q
  n_free_params_A <- model$setup$n_free_params_A
  
  SV_priors_RW <- list(
    theta_A              = rep(0, n_free_params_A),
    Omega_A              = diag(1000, n_free_params_A),
    theta_log_lambda_1   = rep(0, k),
    Omega_log_lambda_1   = diag(1000, k),
    alpha_phi            = rep(5, k),
    beta_phi             = (rep(5, k) - 1) * rep(0.1, k)
  )
  
  model <- priors(model,
                  theta_Psi = rep(0, k),
                  Omega_Psi = diag(0.1, k, k),
                  SV = TRUE,
                  SV_type = "RW",
                  SV_priors = SV_priors_RW)
  
  model <- suppressWarnings(
    fit(model, H = 1, iter = 40, warmup = 20, chains = 1, cores = 1,
        control = list(max_treedepth = 12, adapt_delta = 0.85))
  )
  
  N_est <- nrow(data) - p
  
  expect_equal(dim(model$fit$posterior_means$beta), c(k * p, k))
  expect_equal(dim(model$fit$posterior_means$Psi), c(k, q))
  expect_equal(dim(model$fit$posterior_means$A), c(k, k))
  expect_equal(length(model$fit$posterior_means$phi), k)
  expect_equal(dim(model$fit$posterior_means$Sigma_u), c(N_est, k, k))
  
  expect_null(model$fit$posterior_means$gamma_0)
  expect_null(model$fit$posterior_means$gamma_1)
  expect_null(model$fit$posterior_means$Phi)
})

test_that("summary dimensions are correct for AR1 stochastic volatility model", {
  data <- matrix(rnorm(200), nrow = 100, ncol = 2)
  model <- bvar(data)
  model <- SteadyStateBVAR::setup(model, p = 1, deterministic = "constant")
  
  k <- model$setup$k
  p <- model$setup$p
  q <- model$setup$q
  n_free_params_A <- model$setup$n_free_params_A
  
  SV_priors_AR1 <- list(
    theta_A               = rep(0, n_free_params_A),
    Omega_A               = diag(1000, n_free_params_A),
    theta_gamma_0         = rep(0.1, k),
    Omega_gamma_0         = diag(1000, k),
    theta_gamma_1         = rep(0.9, k),
    Omega_gamma_1         = diag(10, k),
    theta_log_lambda_1    = rep(0.1, k) / (1 - rep(0.9, k)),
    Omega_log_lambda_1    = diag(1000, k),
    V_Phi                 = (10 - k - 1) * diag(k),
    m_Phi                 = 10
  )
  
  model <- priors(model,
                  theta_Psi = rep(0, k),
                  Omega_Psi = diag(0.1, k, k),
                  SV = TRUE,
                  SV_type = "AR1",
                  SV_priors = SV_priors_AR1)
  
  model <- suppressWarnings(
    fit(model, H = 1, iter = 40, warmup = 20, chains = 1, cores = 1,
        control = list(max_treedepth = 12, adapt_delta = 0.85))
  )
  
  N_est <- nrow(data) - p
  
  expect_equal(dim(model$fit$posterior_means$beta), c(k * p, k))
  expect_equal(dim(model$fit$posterior_means$Psi), c(k, q))
  expect_equal(dim(model$fit$posterior_means$A), c(k, k))
  expect_equal(length(model$fit$posterior_means$gamma_0), k)
  expect_equal(length(model$fit$posterior_means$gamma_1), k)
  expect_equal(dim(model$fit$posterior_means$Phi), c(k, k))
  expect_equal(dim(model$fit$posterior_means$Sigma_u), c(N_est, k, k))
  
  expect_null(model$fit$posterior_means$phi)
})