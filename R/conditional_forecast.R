conditional_forecast <- function(x,
                                 conditions,
                                 ci=0.95,
                                 fcst_type = c("mean", "median"),
                                 estimation = c("stan", "gibbs")){
  fcst_fun <- match.arg(fcst_type)
  estim    <- match.arg(estimation)
  ################################################################################
  p      <- x$setup$p
  k      <- x$setup$k
  Y      <- x$setup$Y
  X      <- x$setup$X
  N      <- x$setup$N
  d_pred <- x$predict$d_pred
  ################################################################################
  # Step 1: define post-warmup draws ("It-Bu"),
  #         forecast horizon h
  #         and the v conditions to be imposed on the series
  
  #post warmup draws
  if (estim == "stan") {
    posterior <- rstan::extract(x$fit$stan)
    n_draws   <- dim(posterior$beta)[1]
  } else {
    n_draws   <- dim(x$fit$gibbs$beta_draws)[3]
  }
  #forecast horizon
  H       <- x$predict$H
  #v conditions
  #these are in 'conditions' argument
  ################################################################################
  cond_forecast_array        <- array(NA, dim = c(H, k, n_draws))
  s                          <- k * H
  
  for (n in 1:n_draws) {
    
    # Step 2: at iteration n, recycle Gibbs draws
    
    if (estim == "stan") {
      
      beta_n    <- posterior$beta[n,,]
      
      if (!is.null(x$SV)) {
      t       <- dim(posterior$Sigma_u)[2]
      Sigma_n <- posterior$Sigma_u[n,t,,]
      } else {
      Sigma_n <- posterior$Sigma_u[n,,]
      }
      
      Psi_n     <- posterior$Psi[n,,]
      D_n       <- t(chol(Sigma_n))
      
    } else {
      
      beta_n    <- x$fit$gibbs$beta_draws[,,n]
      Sigma_n   <- x$fit$gibbs$Sigma_u_draws[,,n]
      Psi_n     <- x$fit$gibbs$Psi_draws[,,n]
      D_n       <- t(chol(Sigma_n))
      
    }
    
    # Step 3: at iteration n, compute the unconditional forecasts excl. shocks
    
    Pi_n <- vector("list", p)
    for (i in 1:p) {
      rows_idx <- ((i-1)*k + 1):(i*k)
      Pi_n[[i]] <- t(beta_n[rows_idx, ])
    }
    
    Y_pred_mat <- matrix(NA, nrow=H, ncol=k)
    
    for (h in 1:H) {
      ytilde_t <- d_pred[h, ] %*% t(Psi_n)
      if (h > 1) {
        for (i in 1:min(h-1, p)) {
          term     <- (Y_pred_mat[h-i, ] - d_pred[h-i, ] %*% t(Psi_n)) %*% t(Pi_n[[i]])
          ytilde_t <- ytilde_t + term
        }
      }
      if (h <= p) {
        for (i in h:p) {
          term     <- (Y[N+h-i, ] - X[N+h-i, ] %*% t(Psi_n)) %*% t(Pi_n[[i]])
          ytilde_t <- ytilde_t + term
        }
      }
      Y_pred_mat[h, ] <- ytilde_t
    }
    # Step 4: compute impulse response matrices Phi and Phi_tilde
    #         NB! my notation: Phi = BEAR notation: Psi
    #         Since we use Psi already for steady-state parameter matrix
    
    Pi_n <- t(beta_n)
    
    Phi_n <- MTS::VARpsi(Pi_n, lag=H)$psi #my notation: Phi = MTS notation:psi
    
    n_Phi <- ncol(Phi_n) / k
    
    Phi_n_tilde <- lapply(1:n_Phi, function(i) {
      cols <- ((i-1)*k + 1):(i*k)
      Phi_n[, cols] %*% D_n
    })
    
    # Step 5: build R (v x s) and r (v x 1)
    
    v <- nrow(conditions)
    R <- matrix(0, v, s)
    r <- rep(0, v)
    
    for (row_idx in 1:v) {
      i     <- conditions$var[row_idx]
      h     <- conditions$horizon[row_idx]
      value <- conditions$value[row_idx]
      
      row <- c()
      for (j in 1:h) {
        row <- c(row, Phi_n_tilde[[h-j+1]][i, ])
      }
      R[row_idx, 1:length(row)] <- row
      r[row_idx] <- value - Y_pred_mat[h, i]
    }
    
    # Step 6: draw constrained shocks
    svd_R  <- svd(R, nu = v, nv = s)
    U      <- svd_R$u
    D      <- diag(svd_R$d)
    V_1     <- svd_R$v[, 1:v]
    V_2     <- svd_R$v[, (v+1):s]
    
    lambda <- rnorm(s - v)
    
    eta    <- V_1%*%solve(D)%*%t(U)%*%r + V_2%*%lambda
    
    # Step 7: calculate the conditional forecasts
    #         with the unconditional forecast values
    #         and the constrained shocks
    
    eta_mat <- matrix(eta, nrow = k, ncol = H)
    Y_cond_mat <- matrix(NA, nrow = H, ncol = k)
    
    for (h in 1:H) {
      sum_shocks <- matrix(0, nrow = k, ncol = 1)
      for (j in 1:h) {
        Psi_tilde_h_minus_j <- Phi_n_tilde[[h - j + 1]]
        eta_t_plus_j <- eta_mat[, j]
        sum_shocks <- sum_shocks + Psi_tilde_h_minus_j %*% eta_t_plus_j
      }
      y_tilde_t_plus_h <- Y_pred_mat[h, ]
      Y_cond_mat[h, ] <- y_tilde_t_plus_h + sum_shocks
    }
    cond_forecast_array[,,n] <- Y_cond_mat
    # Step 8: Repeat the above for all n draws
  }
  
  ################################### END ########################################
  
  alpha <- 1 - ci
  return(list(
    cond_draws = cond_forecast_array,
    cond_point  = apply(cond_forecast_array, c(1,2), fcst_fun),
    cond_lower = apply(cond_forecast_array, c(1,2), quantile, probs = alpha/2),
    cond_upper = apply(cond_forecast_array, c(1,2), quantile, probs = 1 - alpha/2)
  ))
}