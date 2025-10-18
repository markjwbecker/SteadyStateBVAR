BVAR_setup <- function(Z, p, det=c("c", "c&t", "c&d"), dummy=NULL) {

  N = dim(Z)[1]-p
  m = dim(Z)[2]
  d=2
  det <- match.arg(det)
  if (det == "c") {
    x <- cbind(rep(1, nrow(Z)))
  } else if (det == "c&t") {
    x <- cbind(rep(1, nrow(Z)), 1:nrow(Z))
  } else { # det == "c&d"
    if (is.null(dummy)) {
      print("input dummy variable")
    } else {
      x <- cbind(rep(1, nrow(Z)), dummy)
    }
  }

  Y <- Z[-c(1:p), ]
  W <- embed(Z, dimension = p+1)[, -(1:m)]
  X <- x[-c(1:p), ]
  Q <- embed(x, dimension = p+1)[, -(1:d)]
  obj=list(N=N,m=m,p=p,Y=Y,X=X,W=W,Q=Q)
  return(obj)
}
