estimate_gibbs <- function(x, iter, warmup, H, X_pred, Jeffrey=FALSE){
  #### Algorithm 4 from:
  #### Karlsson, S. (2013). Forecasting with Bayesian Vector Autoregression.
  #### In: Elliott, G. and Timmerman, A. (eds) Handbook of Economic Forecasting.
  #### Elsevier B.V. Vol 2, Part B., pp. 791-897.
  
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
  
  Gamma_d_OLS <- setup$beta_hat[1:(k*p),]
  Lambda_OLS <- setup$Psi_OLS
  
  gamma_d_lbar <-  priors$theta_beta
  Sigma_d_lbar <- priors$Omega_beta
  lambda_lbar <- priors$theta_Psi
  Sigma_lambda_lbar <- priors$Omega_Psi
  
  if (isFALSE(Jeffrey)){
    v_ = priors$m_0
    S_ = priors$V_0
  }
  
  Psi <- vector(mode = "list", length = iter+1)
  gamma_d <- vector(mode = "list", length = iter+1)
  lambda <- vector(mode = "list", length = iter+1)
  
  R <- iter-warmup
  B <- warmup
  forecasts_array <- array(NA, dim = c(H, k, R))
  
  ############### START ###############
  
  #### SELECT STARTING VALUES (OLS ESTIMATES) ####
  
  gamma_d[[1]] <- c(Gamma_d_OLS) #gamma_d^{0}
  lambda[[1]]  <- c(Lambda_OLS) #lambda_d^{0}
  
  for (j in 2:(B+R+1)){
    
    Lambda = matrix(lambda[[j-1]],k,q)
    Gamma_d = matrix(gamma_d[[j-1]],k*p,k)
    
    ############   1   ################
    ############ EQ 29 ################
    U = Y - X%*%t(Lambda) - (W-Q%*%(diag(p) %x% t(Lambda))) %*% Gamma_d
    S = crossprod(U)
    N = nrow(U)
    if (isFALSE(Jeffrey)){
      Psi[[j]] = LaplacesDemon::rinvwishart(N+v_, S+S_)
    } else {
      Psi[[j]] = LaplacesDemon::rinvwishart(N, S)  
    }
    ############   2   ################
    ############ EQ 30 ################
    Y_Lambda = Y-X%*%t(Lambda)
    W_Lambda = (W-Q%*%(diag(p)%x%t(Lambda)))
    
    Sigma_d_bar = solve((solve(Sigma_d_lbar)+solve(Psi[[j]]) %x% (crossprod(W_Lambda))))
    
    gamma_d_bar = Sigma_d_bar %*% (solve(Sigma_d_lbar)%*%gamma_d_lbar + c(t(W_Lambda)%*%Y_Lambda%*%solve(Psi[[j]])))
    
    gamma_d[[j]] = MASS::mvrnorm(1, gamma_d_bar, Sigma_d_bar)
    
    ############   3   ################
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
    
    Sigma_lambda_bar = solve(solve(Sigma_lambda_lbar)+t(F) %*% ((crossprod(B))%x%solve(Psi[[j]])) %*% F)
    lambda_bar = Sigma_lambda_bar %*% (solve(Sigma_lambda_lbar)%*%lambda_lbar+t(F)%*%c(solve(Psi[[j]])%*%t(Y_gamma)%*%B))
    
    lambda[[j]] <- MASS::mvrnorm(1, lambda_bar, Sigma_lambda_bar)
    
    ############   4   ################
    if (j > (B+1)) {
    Lambda <- matrix(lambda[[j]], nrow = k, ncol = q)
    Gamma_d <- matrix(gamma_d[[j]], nrow = k * p, ncol = k)
    Psi <- Psi[[j]]                                    
    
    A <- vector("list", p)
    for (i in 1:p) {
      rows_idx <- ((i - 1) * k + 1):(i * k)
      A[[i]] <- t(Gamma_d[rows_idx,])  # k x k
    }
    
    Y_pred_mat <- matrix(NA, nrow = H, ncol = k)
    
    for (h in 1:H) {
      
      u_t <- MASS::mvrnorm(1, rep(0, k), Psi)
      ytilde_t <- X_pred[h, ] %*% t(Lambda)
      
      if (h > 1) {
        for (i in 1:min(h - 1, p)) {
          term <- (Y_pred_mat[h - i,] - X_pred[h - i,] %*% t(Lambda)) %*% t(A[[i]])
          ytilde_t <- ytilde_t + term
        }
      }
      
      if (h <= p) {
        for (i in h:p) {
          term <- (Y[N + h - i,] - X[N + h - i,] %*% t(Lambda)) %*% t(A[[i]])
          ytilde_t <- ytilde_t + term
        }
      }
      
      Y_pred_mat[h, ] <- ytilde_t + u_t
    }
    
    forecasts_array[, ,j] <- Y_pred_mat
    }
  }
  
  #### The first B draws are discarded as burn-in ###
  keep_idx <- seq(B + 2, B+R+1)
  n_keep <- length(seq(B + 2, B+R+1))
  gamma_d_keep  <- gamma_d[keep_idx]
  lambda_keep <- lambda[keep_idx]
  Psi_keep    <- Psi[keep_idx]
  
  
  ############### FINISH ###############
  
  beta_draws <- array(unlist(gamma_d_keep), dim = c(k * p, k, n_keep))
  Psi_draws <- array(unlist(lambda_keep), dim = c(k, q, n_keep))
  Sigma_u_draws <- simplify2array(Psi_keep)
  
  beta_posterior_mean <- apply(beta_draws, c(1, 2), mean)
  Psi_posterior_mean <- apply(Psi_draws, c(1, 2), mean)
  Sigma_u_posterior_mean <- apply(Sigma_u_draws, c(1, 2), mean)

  res <- list(beta_draws = beta_draws,
              Psi_draws = Psi_draws,
              Sigma_u_draws = Sigma_u_draws,
              fcst_draws = forecasts_array,
              beta_posterior_mean = beta_posterior_mean,
              Psi_posterior_mean = Psi_posterior_mean,
              Sigma_u_posterior_mean = Sigma_u_posterior_mean)
  return(res)
}
