fit <- function(x, iter = 5000, warmup = 2500, chains = 2, estimation = c("stan", "gibbs")) {
  
  estimation <- match.arg(estimation)
  Jeffrey <- x$priors$Jeffrey
  
  if (estimation == "stan") {
    
    if (!is.null(x$SV)) {
      stan_data <- c(
        x$setup,
        x$priors
      )
    } else {
      stan_data <- c(
        x$setup,
        list(H = x$predict$H, X_pred = x$predict$X_pred),
        x$priors
      )
    }
    
    stan_file <- if (isFALSE(Jeffrey)) {
      system.file("inv_wishart_cov.stan", package = "SteadyStateBVAR")
    } else {
      system.file("diffuse_cov.stan", package = "SteadyStateBVAR")
    }
    
    if (isTRUE(x$SV)) {
      stan_file <- system.file("stochastic_volatility.stan", package = "SteadyStateBVAR")
      stan_data <- c(x$SV_priors, stan_data)
      k <- x$setup$k
      if (k == 2) {
        stan_data$theta_A <- as.array(x$SV_priors$theta_A[1])
      }
      
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
    
    ### compute forecasts outside of stan. temporary fix until i get it to work ###
    if (isTRUE(x$SV)) {
      
      draws <- rstan::extract(x$fit$stan)
      n_draws <- dim(draws$beta)[1]
      p <- x$setup$p
      k <- x$setup$k
      q <- x$setup$q
      H <- x$predict$H
      X_pred <- x$predict$X_pred
      forecast_array <- array(NA, dim = c(H, k, n_draws))
      Sigma_u_array <- array(NA, dim = c(n_draws, H, k, k))
      log_lambda_array <- array(NA, dim = c(n_draws, H, k))
      X <- x$setup$X
      Y <- x$setup$Y
      N <- dim(Y)[1]
      
      for (j in seq_len(n_draws)) {
        
        Psi     <- draws$Psi[j,,]
        beta    <- draws$beta[j,,]
        Phi     <- draws$Phi[j,,]
        A       <- draws$A[j,,]
        gamma_0 <- draws$gamma_0[j,]
        gamma_1 <- draws$gamma_1[j,]
        log_lambda_last <- draws$log_lambda[j,N,]
        
        Pi <- vector("list", p)
        for (i in 1:p) {
          rows_idx <- ((i - 1) * k + 1):(i * k)
          Pi[[i]] <- t(beta[rows_idx,])
        }
        
        Y_pred <- matrix(NA, nrow = H, ncol = k)
        log_lambda_pred <- matrix(NA, nrow = H, ncol = k)
        Sigma_u_pred <- vector("list", H)
        u_t <- rep(NA, k)
        yhat_t <- rep(NA, k)
        
        nu_t <- MASS::mvrnorm(1, rep(0, k), Phi)
        log_lambda_pred[1,] = gamma_0 + gamma_1 * log_lambda_last + nu_t;
        for (h in 2:H) {
          nu_t <- MASS::mvrnorm(1, rep(0, k), Phi)
          log_lambda_pred[h,] = gamma_0 + gamma_1 * log_lambda_pred[h-1,] + nu_t;
        }
        log_lambda_array[j,,] <- log_lambda_pred
        Y_pred_mat <- matrix(NA, nrow = H, ncol = k)
        for (h in 1:H) {
          
          Lambda_h = diag(exp(log_lambda_pred[h,]))
          Sigma_u_pred[[h]] = solve(A) %*% Lambda_h %*% t(solve(A))
          u_t <- MASS::mvrnorm(1, rep(0, k), Sigma_u_pred[[h]])
          
          yhat_t <- X_pred[h, ] %*% t(Psi)
          
          if (h > 1) {
            for (i in 1:min(h - 1, p)) {
              term <- (Y_pred_mat[h - i,] - X_pred[h - i,] %*% t(Psi)) %*% t(Pi[[i]])
              yhat_t <- yhat_t + term
            }
          }
          
          if (h <= p) {
            for (i in h:p) {
              term <- (Y[N + h - i,] - X[N + h - i,] %*% t(Psi)) %*% t(Pi[[i]])
              yhat_t <- yhat_t + term
            }
          }
          
          Y_pred_mat[h, ] <- yhat_t + u_t
          Sigma_u_array[j, h, , ] <- Sigma_u_pred[[h]]
        }
        
        forecast_array[,,j] <- Y_pred_mat
      }
      x$fit$stanf$fcst_draws 
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
      X_pred = x$predict$X_pred,
      Jeffrey = Jeffrey
    )
  }
  
  return(x)
}