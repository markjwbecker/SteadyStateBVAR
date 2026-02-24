setup <- function(x, ...) UseMethod("setup")

setup.bvar <- function(x, p, deterministic=c("constant", "constant_and_dummy", "constant_and_trend"), dummy=NULL) {
  
  yt <- x$data
  N = nrow(yt)-p
  k = ncol(yt)
  
  deterministic <- match.arg(deterministic)
  
  if (deterministic == "constant") {
    xt <- cbind(rep(1, nrow(yt)))
    q <- 1
  } else if (deterministic == "constant_and_dummy") {
    xt <- cbind(rep(1, nrow(yt)), dummy)
    q <- 2
  } else {
    trend <- 1:nrow(yt)
    xt <- cbind(rep(1, nrow(yt)), trend)
    q <- 2
  }
  
  Y <- yt[-c(1:p), ]
  W <- embed(yt, dimension = p+1)[, -(1:k)]
  X <- xt[-c(1:p), ,drop=F]
  Q <- embed(xt, dimension = p+1)[, -(1:q), drop=F]
  
  Z <- cbind(W,X)
  beta_hat = solve(crossprod(Z),crossprod(Z,Y))
  U = Y-Z%*%beta_hat
  Sigma_u_OLS <- crossprod(U)/(N-k*p-q)
  
  if (q == 1) {
    C_hat <- beta_hat[(k*p+1):(k*p+q),]
  } else {
    C_hat <- t(beta_hat[(k*p+1):(k*p+q),])
  }
  A <- vector("list", p)
  for (i in 1:p) {
    rows_idx <- ((i - 1) * k + 1):(i * k)
    A[[i]] <- matrix(t(beta_hat[rows_idx, ]), nrow = k, ncol = k)
  }
  A_L <- diag(k)
  for (i in 1:p) {
    A_L <- A_L - A[[i]]
  }
  Psi_OLS <- solve(A_L, C_hat)
  beta_OLS = beta_hat[1:(k*p),]
  x$setup <- list(N=N,
                  k=k,
                  p=p,
                  Y=Y,
                  X=X,
                  W=W,
                  Q=Q,
                  q=q,
                  dummy=dummy,
                  beta_OLS=beta_OLS,
                  Sigma_u_OLS=Sigma_u_OLS,
                  Psi_OLS=Psi_OLS,
                  xt=xt)
  return(x)
}
