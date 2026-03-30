fit <- function(x, iter = 5000, warmup = 2500, chains = 2, estimation = c("stan", "gibbs")) {
  
  estimation <- match.arg(estimation)
  Jeffrey <- x$priors$Jeffrey
  
  if (estimation == "stan") {
      stan_data <- c(
        x$setup,
        list(H = x$predict$H, d_pred = x$predict$d_pred),
        x$priors
      )

    
    stan_file <- if (isFALSE(Jeffrey)) {
      system.file("inv_wishart_cov.stan", package = "SteadyStateBVAR")
    } else {
      system.file("diffuse_cov.stan", package = "SteadyStateBVAR")
    }
    
    if (isTRUE(x$SV)) {
      if (x$SV_type == "RW") {
        stan_file <- system.file("stochastic_volatility_RW.stan", package = "SteadyStateBVAR")
      } else if (x$SV_type == "AR"){
        stan_file <- system.file("stochastic_volatility_stationaryAR.stan", package = "SteadyStateBVAR")
      } else {
        print("Please specify bvar_obj$SV_type")
      }
      stan_data <- c(x$SV_priors, stan_data)
      k <- x$setup$k
      if (k == 2) {
        stan_data$theta_A <- as.array(x$SV_priors$theta_A[1])
      }
      
    }
    
    if (isTRUE(x$priors$multi_student_t)) {
      stan_file <- system.file("multivariate_t.stan", package = "SteadyStateBVAR")
    }
    rstan::rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
    
    x$fit$stan <- rstan::stan(
      file = stan_file,
      data = stan_data,
      iter = iter,
      warmup = warmup,
      chains = chains,
      verbose = FALSE
    )
    posterior <- rstan::extract(x$fit$stan)
    if (!isTRUE(x$SV)) {
      x$posterior_means$beta <- apply(posterior$beta, c(2, 3), mean)
      x$posterior_means$Psi <- apply(posterior$Psi, c(2, 3), mean)
      x$posterior_means$Sigma_u <- apply(posterior$Sigma_u, c(2, 3), mean)
      
    } else if (isTRUE(x$SV) && x$SV_type == "AR"){
      
      x$posterior_means$beta <- apply(posterior$beta, c(2, 3), mean)
      x$posterior_means$Psi <- apply(posterior$Psi, c(2, 3), mean)
      x$posterior_means$A <- apply(posterior$A, c(2, 3), mean)
      x$posterior_means$gamma_0 <- apply(posterior$gamma_0, 2, mean)
      x$posterior_means$gamma_1 <- apply(posterior$gamma_1, 2, mean)
      x$posterior_means$Phi <- apply(posterior$Phi, c(2, 3), mean)
      
      
    } else if (isTRUE(x$SV) && x$SV_type == "RW"){
      x$posterior_means$beta <- apply(posterior$beta, c(2, 3), mean)
      x$posterior_means$Psi <- apply(posterior$Psi, c(2, 3), mean)
      x$posterior_means$A <- apply(posterior$A, c(2, 3), mean)
      x$posterior_means$phi <- apply(posterior$phi, 2, mean)
    }
  } else {
    
    if (chains != 1) {
      stop("For Gibbs sampling, 'chains' must be equal to 1.")
    }
    
    x$fit$gibbs <- estimate_gibbs(
      x = x,
      iter = iter,
      warmup = warmup,
      H = x$predict$H,
      d_pred = x$predict$d_pred,
      Jeffrey = Jeffrey
    )
    
    x$posterior_means$beta <- x$fit$gibbs$beta_posterior_mean
    x$posterior_means$Psi <- x$fit$gibbs$Psi_posterior_mean
    x$posterior_means$Sigma_u <- x$fit$gibbs$Sigmna_u_posterior_mean
    
  }
  
  return(x)
}