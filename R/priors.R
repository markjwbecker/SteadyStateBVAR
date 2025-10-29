priors <- function(yt, p,lambda1, lambda2, fol_pm, Psi_0, vec_Psi_vars, dummy=NULL){
  k = ncol(yt)
  dummy <- ts(dummy, start=start(xt), frequency=frequency(xt))
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
  Sigma_vec_beta <- diag(c(tmp_mat))
  mat <- matrix(0, nrow = k*p, ncol = k)
  for (i in 1:k){
    mat[i,i] <- fol_pm[i]
  }
  vec_beta_0 = c(mat)
  
  
  mp2 <- bvartools::minnesota_prior(
    obj,
    kappa0 = lambda1^2,
    kappa1 = lambda2^2,
    sigma = "VAR"
  )
  
  Sigma_u_hat = solve(mp2$sigma_i)
  m_0=k+2
  V_0 = (m_0-k-1)*Sigma_u_hat
  
  vec_Psi_0 = c(Psi_0)
  Sigma_vec_Psi = diag(vec_Psi_vars)
  
  priors <- list()
  
  priors$vec_beta_0 <- vec_beta_0
  priors$Sigma_vec_beta <- Sigma_vec_beta
  
  priors$vec_Psi_0 <- vec_Psi_0
  priors$Sigma_vec_Psi <- Sigma_vec_Psi
  
  priors$V_0 <- V_0
  priors$m_0 <- m_0
  
  return(priors)
}