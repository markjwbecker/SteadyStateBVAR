BVAR_setup <- function(Z, p, det=c("c", "c&d"), dummy=NULL) {

  N = dim(Z)[1]-p
  m = dim(Z)[2]
  det <- match.arg(det)
  if (det == "c") {
    x <- cbind(rep(1, nrow(Z)))
    d <- 1
  } else if (det == "c&t") {
    x <- cbind(rep(1, nrow(Z)), 1:nrow(Z))
    d <- 2
  } else if (det == "c&d"){
    x <- cbind(rep(1, nrow(Z)), dummy)
    d <- 2
  }
  
  Y <- Z[-c(1:p), ]
  W <- embed(Z, dimension = p+1)[, -(1:m)]
  X <- x[-c(1:p), ]
  Q <- embed(x, dimension = p+1)[, -(1:d)]
  
  stan_data=list(N=N,m=m,p=p,Y=Y,X=X,W=W,Q=Q,d=d)
  return(stan_data)
}
