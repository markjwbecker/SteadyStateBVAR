priors <- function(x, ...) UseMethod("priors")

priors.bvar <- function(x, lambda_1=0.2, lambda_2=0.5, fol_pm=NULL, theta_Psi=NULL, Omega_Psi=NULL, Jeffrey=FALSE){
  
  setup <- x$setup
  yt <- x$data
  k <- setup$k
  p <- setup$p
  
  dummy <- setup$dummy
  if (!is.null(dummy)) dummy <- ts(dummy, start=start(yt), frequency=frequency(yt))
  
  if (is.null(dummy)){
    obj <- bvartools::gen_var(as.ts(yt), p)
  } else {
    obj <- bvartools::gen_var(as.ts(yt), p, exogen=dummy, s=0)
  }
  
  mp <- bvartools::minnesota_prior(
    obj,
    kappa0 = lambda_1^2,
    kappa1 = lambda_2^2,
    sigma = "AR"
    )
  
  tmp <- diag(solve(mp$v_i)[1:(k*p*k),1:(k*p*k)])
  tmp_mat <- matrix(tmp,k*p,k,byrow=TRUE)
  Omega_beta <- diag(c(tmp_mat))
  
  if (is.null(fol_pm)) fol_pm <- rep(0,k)
  
  mat <- matrix(0, nrow = k*p, ncol = k)
  for (i in 1:k){
    mat[i,i] <- fol_pm[i]
  }
  theta_beta = c(mat)
  if(isFALSE(Jeffrey)){
    m_0=k+2
    V_0 = (m_0-k-1)*setup$Sigma_OLS
    priors$V_0 <- V_0
    priors$m_0 <- m_0
  }
  
  priors <- list()
  
  priors$theta_beta <- theta_beta
  priors$Omega_beta <- Omega_beta
  
  priors$theta_Psi <- theta_Psi
  priors$Omega_Psi <- Omega_Psi
  
  priors$Jeffrey <- Jeffrey
  
  x$priors <- priors
  
  return(x)
}