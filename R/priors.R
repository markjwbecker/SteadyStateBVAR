priors<- function(x, lambda_1=0.2, lambda_2=0.5, lambda_3 = 1, first_own_lag_prior_mean=NULL, theta_Psi=NULL, Omega_Psi=NULL, Jeffrey=FALSE){
  
  priors <- list()
  
  setup <- x$setup
  yt <- x$data
  k <- setup$k
  p <- setup$p
  q <- setup$q
  xt <- setup$xt
  dummy <- setup$dummy
  
  Sigma_AR <- diag(0,k)
  
  for (i in 1:k){
    
    y <- yt[,i]
    N = length(y)-p
    
    Y <- y[-c(1:p)]
    W <- embed(y, dimension = p+1)[, -1]
    X <- xt[-c(1:p), ,drop=F]
    Q <- embed(xt, dimension = p+1)[, -(1:q)]
    
    Z <- cbind(W,X)
    beta_OLS = solve(crossprod(Z,Z),crossprod(Z,Y))
    U = Y-Z%*%beta_OLS
    sigma2 <- crossprod(U,U)/(nrow(Z)-ncol(Z))
    Sigma_AR[i,i] <- sigma2
  }
  
  V <- lapply(1:p, function(x) matrix(0, k, k))
  sigma <- sqrt(diag(Sigma_AR))
  
  for (l in 1:p) {
    for (i in 1:k) {
      for (j in 1:k) {
        if (i == j) {
          V[[l]][i,j] <- (lambda_1/(l^lambda_3))^2
        } else {
          V[[l]][i,j] <- ((lambda_1*lambda_2*sigma[i])/(l^lambda_3*sigma[j]))^2
        }
      }
    }
  }
  V_mat <- do.call(cbind, V)
  Omega_beta <- diag(c(t(V_mat)))
  
  if (is.null(first_own_lag_prior_mean)) first_own_lag_prior_mean <- rep(0,k)
  
  mat <- matrix(0, nrow = k*p, ncol = k)
  for (i in 1:k){
    mat[i,i] <- first_own_lag_prior_mean[i]
  }
  theta_beta = c(mat)
  if(isFALSE(Jeffrey)){
    m_0=k+2
    V_0 = (m_0-k-1)*setup$Sigma_OLS
    priors$V_0 <- V_0
    priors$m_0 <- m_0
  }
  
  priors$theta_beta <- theta_beta
  priors$Omega_beta <- Omega_beta
  
  priors$theta_Psi <- theta_Psi
  priors$Omega_Psi <- Omega_Psi
  
  priors$Jeffrey <- Jeffrey
  priors$Sigma_AR <- Sigma_AR
  x$priors <- priors
  
  return(x)
}