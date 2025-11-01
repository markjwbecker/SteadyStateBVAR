priors <- function(yt, p, lambda1, lambda2, fol_pm, theta_Psi, Omega_Psi, dummy=NULL){
  
  k = ncol(yt)
  dummy <- ts(dummy, start=start(yt), frequency=frequency(yt))
  
  if (is.null(dummy)){
    obj <- bvartools::gen_var(as.ts(yt), p)
  } else {
    obj <- bvartools::gen_var(as.ts(yt), p,
                              exogen=dummy,
                              s=0)
  }
  
  mp <- bvartools::minnesota_prior(
    obj,
    kappa0 = lambda1^2,
    kappa1 = lambda2^2,
    sigma = "AR"
  )
  
  tmp <- diag(solve(mp$v_i)[1:(k*p*k),1:(k*p*k)])
  tmp_mat <- matrix(tmp,k*p,k,byrow=TRUE)
  Omega_beta <- diag(c(tmp_mat))
  mat <- matrix(0, nrow = k*p, ncol = k)
  for (i in 1:k){
    mat[i,i] <- fol_pm[i]
  }
  theta_beta = c(mat)
  
  
  mp2 <- bvartools::minnesota_prior(
    obj,
    kappa0 = lambda1^2,
    kappa1 = lambda2^2,
    sigma = "VAR"
  )
  
  Sigma_u_hat = solve(mp2$sigma_i)
  m_0=k+2
  V_0 = (m_0-k-1)*Sigma_u_hat
  
  priors <- list()
  
  priors$theta_beta <- theta_beta
  priors$Omega_beta <- Omega_beta
  
  priors$theta_Psi <- theta_Psi
  priors$Omega_Psi <- Omega_Psi
  
  priors$V_0 <- V_0
  priors$m_0 <- m_0
  
  priors$Sigma_u_OLS <- Sigma_u_hat
  
  return(priors)
}