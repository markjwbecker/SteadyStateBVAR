priors <- function(Y,p,lambda1,lambda2,fol_pm,Lambda_pr_means,Lambda_pr_vars){
  m = ncol(Y)
  
  obj <- gen_var(as.ts(Y), p)
  mp <- bvartools::minnesota_prior(
  obj,
  kappa0 = lambda1^2,
  kappa1 = lambda2^2,
  sigma = "AR"
  )
  
  tmp <- diag(solve(mp$v_i)[1:(m*p*m),1:(m*p*m)])
  tmp_mat <- matrix(tmp,m*p,m,byrow=TRUE)
  Gamma_d_pr_cov <- diag(c(tmp_mat))
  mat <- matrix(0, nrow = m*p, ncol = m)
  for (i in 1:m){
    mat[i,i] <- fol_pm[i]
  }
  Gamma_d_pr_mean = c(mat)
  
  diag_vars <- matrix(0, m, m)
  for (i in 1:m) {
  arfit <- arima(Y[, i], order = c(p, 0, 0), method = "CSS")
  diag_vars[i, i] <- arfit$sigma2
  }
  
  gamma=m+2
  Psi_pr_scale = (gamma-m-1)*diag_vars
  
  Lambda_pr_mean = c(Lambda_pr_means)
  Lambda_pr_cov = diag(Lambda_pr_vars)
  prior <- list()
  prior$Gamma_d_pr_mean <- Gamma_d_pr_mean
  prior$Gamma_d_pr_cov <- Gamma_d_pr_cov
  prior$Lambda_pr_mean <- Lambda_pr_mean
  prior$Lambda_pr_cov <- Lambda_pr_cov
  prior$Psi_pr_scale <- Psi_pr_scale
  prior$gamma <- 5
  return(prior)
}