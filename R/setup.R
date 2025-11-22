setup.bvar <- function(x, p, deterministic=c("constant", "constant_and_dummy"), dummy=NULL) {
  
  deterministic <- match.arg(deterministic)
  yt <- x$data
  N = nrow(yt)-p
  k = ncol(yt)
  
  if (deterministic == "constant") {
    xt <- cbind(rep(1, nrow(yt)))
    q <- 1
  } else {
    xt <- cbind(rep(1, nrow(yt)), dummy)
    q <- 2
  }
  
  Y <- yt[-c(1:p), ]
  W <- embed(yt, dimension = p+1)[, -(1:k)]
  X <- xt[-c(1:p), ,drop=F]
  Q <- embed(xt, dimension = p+1)[, -(1:q)]
  
  Z <- cbind(W,X)
  beta_OLS = solve((t(Z)%*%Z))%*%t(Z)%*%Y
  U = Y-Z%*%beta_OLS
  Sigma_OLS <- t(U)%*%U/(N-k*p-q)
  
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
                Sigma_OLS=Sigma_OLS,
                xt=xt)
  return(x)
}
