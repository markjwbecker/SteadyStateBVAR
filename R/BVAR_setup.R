BVAR_setup <- function(yt, p, deterministic=c("constant", "constant_and_dummy"), dummy=NULL) {
  
  N = dim(yt)[1]-p
  k = dim(yt)[2]
  deterministic <- match.arg(deterministic)
  if (deterministic == "constant") {
    x <- cbind(rep(1, nrow(yt)))
    q <- 1
  } else {
    x <- cbind(rep(1, nrow(yt)), dummy)
    q <- 2
  }
  
  Y <- yt[-c(1:p), ]
  W <- embed(yt, dimension = p+1)[, -(1:k)]
  X <- x[-c(1:p), ]
  Q <- embed(x, dimension = p+1)[, -(1:q)]
  
  stan_data=list(N=N,k=k,p=p,Y=Y,X=X,W=W,Q=Q,q=q)
  return(stan_data)
}
