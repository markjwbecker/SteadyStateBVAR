estimate_gibbs <- function(x, iter, warmup, H, X_pred, Jeffrey=FALSE){
  
  setup <- x$setup
  priors <- x$priors
  
  Y <- setup$Y
  X <- setup$X
  W <- setup$W
  Q <- setup$Q
  N <- setup$N
  k  <- setup$k
  p  <- setup$p
  q  <- setup$q
  
  beta_hat <- setup$beta_OLS
  SigOLS <- setup$Sigma_OLS
  
  Gamma_d_OLS <- beta_hat[1:(k*p),]
  
  C_hat <- t(beta_hat[(k*p+1):(k*p+q),])
  A <- vector("list", p)
  for (i in 1:p) {
    rows_idx <- ((i - 1) * k + 1):(i * k)
    A[[i]] <- matrix(t(beta_hat[rows_idx, ]), nrow = k, ncol = k)
  }
  A_L <- diag(k)
  for (i in 1:p) {
    A_L <- A_L - A[[i]]
  }
  
  Lambda_OLS <- solve(A_L) %*% C_hat
  
  gamma_d_lbar <-  priors$theta_beta
  Sigma_d_lbar <- priors$Omega_beta
  
  lambda_lbar <- priors$theta_Psi
  Sigma_lambda_lbar <- priors$Omega_Psi
  
  if (isFALSE(Jeffrey)){
    v_ = priors$m_0
    S_ = priors$V_0
  }
  
  Psi <- vector(mode = "list", length = iter)
  gamma_d <- vector(mode = "list", length = iter)
  lambda <- vector(mode = "list", length = iter)
  
  gamma_d[[1]] <- c(Gamma_d_OLS)
  lambda[[1]]  <- c(Lambda_OLS) 
  
  for (j in 2:iter){
    
    Lambda = matrix(lambda[[j-1]],k,q)
    Gamma_d = matrix(gamma_d[[j-1]],k*p,k)
    
    ############ EQ 29 ################
    U = Y - X%*%t(Lambda) - (W-Q%*%(diag(p) %x% t(Lambda))) %*% Gamma_d
    S = t(U)%*%U
    N = nrow(U)
    if (isFALSE(Jeffrey)){
      Psi[[j]] = rinvwishart(N+v_, S+S_)
    } else {
      Psi[[j]] = rinvwishart(N, S)  
    }
    ############ EQ 30 ################
    Y_Lambda = Y-X%*%t(Lambda)
    W_Lambda = (W-Q%*%(diag(p)%x%t(Lambda)))
    
    Sigma_d_bar = solve((solve(Sigma_d_lbar)+solve(Psi[[j]]) %x% (t(W_Lambda) %*% W_Lambda)))
    
    gamma_d_bar = Sigma_d_bar %*% (solve(Sigma_d_lbar)%*%gamma_d_lbar + c(t(W_Lambda)%*%Y_Lambda%*%solve(Psi[[j]])))
    
    gamma_d[[j]] = mvrnorm(1, gamma_d_bar, Sigma_d_bar)
    
    ############ EQ 31 ################
    Gamma_d <- matrix(gamma_d[[j]],k*p,k)
    A_list <- vector("list", p)
    for (i in 1:p) {
      rows_idx <- ((i - 1) * k + 1):(i * k)
      A_list[[i]] <- t(Gamma_d[rows_idx, ])
    }
    blocks <- list(diag(k * q))
    for (i in 1:p) {
      blocks[[i + 1]] <- diag(q) %x% t(A_list[[i]])
    }
    F_prime <- do.call(cbind, blocks)
    F <- t(F_prime)
    
    B = cbind(X,-Q)
    Y_gamma = Y-W%*%Gamma_d
    
    Sigma_lambda_bar = solve(solve(Sigma_lambda_lbar)+t(F) %*% ((t(B)%*%B)%x%solve(Psi[[j]])) %*% F)
    lambda_bar = Sigma_lambda_bar %*% (solve(Sigma_lambda_lbar)%*%lambda_lbar+t(F)%*%c(solve(Psi[[j]])%*%t(Y_gamma)%*%B))
    
    lambda[[j]] <- mvrnorm(1, lambda_bar, Sigma_lambda_bar)
  }
  
  burnin <- warmup
  keep_idx <- seq(burnin + 1, iter)
  
  gamma_keep  <- gamma_d[keep_idx]
  lambda_keep <- lambda[keep_idx]
  Psi_keep    <- Psi[keep_idx]
  
  Gamma_d_array <- array(unlist(gamma_keep), dim = c(k * p, k, length(gamma_keep)))
  lambda_array <- array(unlist(lambda_keep), dim = c(k, q, length(lambda_keep)))
  Psi_array <- simplify2array(Psi_keep)
  
  Gamma_d_post <- apply(Gamma_d_array, c(1, 2), mean)
  lambda_post <- apply(lambda_array, c(1, 2), mean)
  Psi_post <- apply(Psi_array, c(1, 2), mean)
  
  n_draws <- length(lambda_keep)
  
  dummy = setup$dummy
  forecasts_array <- array(NA, dim = c(H, k, n_draws))
  
  for (j in seq_len(n_draws)) {
    Lambda <- matrix(lambda_keep[[j]], nrow = k, ncol = q)
    Gamma_d <- matrix(gamma_keep[[j]], nrow = k * p, ncol = k)
    Psi <- Psi_keep[[j]]                                    
    
    A <- vector("list", p)
    for (i in 1:p) {
      rows_idx <- ((i - 1) * k + 1):(i * k)
      A[[i]] <- t(Gamma_d[rows_idx,])  # k x k
    }
    
    Y_pred_mat <- matrix(NA, nrow = H, ncol = k)
    
    for (h in 1:H) {
      
      u_t <- mvrnorm(1, rep(0, k), Psi)
      yhat_t <- X_pred[h, ] %*% t(Lambda)
      
      if (h > 1) {
        for (i in 1:min(h - 1, p)) {
          term <- (Y_pred_mat[h - i,] - X_pred[h - i,] %*% t(Lambda)) %*% t(A[[i]])
          yhat_t <- yhat_t + term
        }
      }
      
      if (h <= p) {
        for (i in h:p) {
          term <- (Y[N + h - i,] - X[N + h - i,] %*% t(Lambda)) %*% t(A[[i]])
          yhat_t <- yhat_t + term
        }
      }
      
      Y_pred_mat[h, ] <- yhat_t + u_t
    }
    
    forecasts_array[, ,j] <- Y_pred_mat
  }

  res <- list(beta_draws = Gamma_d_array,
              Psi_draws = lambda_array,
              Sigma_u_draws = Psi_array,
              fcst_draws = forecasts_array,
              beta_post_mean=Gamma_d_post,
              Psi_post_mean=lambda_post,
              Sigma_u_post_mean = Psi_post)
  return(res)
}
