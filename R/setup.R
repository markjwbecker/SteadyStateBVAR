#' Set up a BVAR model
#'
#' A generic function for setting up a BVAR model. Computes OLS estimates and
#' prepares all matrices needed for prior specification and estimation.
#'
#' @param x A \code{bvar} object created by \code{\link{bvar}}.
#' @param ... Further arguments passed to methods.
#'
#' @return The \code{bvar} object with a \code{setup} list appended.
#' @export
setup <- function(x, ...) UseMethod("setup")

#' @rdname setup
#'
#' @param p Integer. The lag order of the VAR.
#' @param deterministic Character. The deterministic component to include.
#'   One of \code{"constant"} (default), \code{"constant_and_dummy"}, or
#'   \code{"constant_and_trend"}.
#' @param dummy Optional numeric vector or matrix of dummy variables. Only
#'   used when \code{deterministic = "constant_and_dummy"}. Default \code{NULL}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' model <- bvar(data = my_data)
#' model <- setup(model, p = 2, deterministic = "constant")
#' }
setup.bvar <- function(x, p, deterministic=c("constant", "constant_and_dummy", "constant_and_trend"), dummy=NULL, ...) {
  
  yt <- x$data
  N = nrow(yt)-p
  k = ncol(yt)
  
  deterministic <- match.arg(deterministic)
  
  if (deterministic == "constant") {
    dt <- cbind(rep(1, nrow(yt)))
    q <- 1
  } else if (deterministic == "constant_and_dummy") {
    dt <- cbind(rep(1, nrow(yt)), dummy)
    q <- 2
  } else {
    trend <- 1:nrow(yt)
    dt <- cbind(rep(1, nrow(yt)), trend)
    q <- 2
  }
  
  Y <- yt[-c(1:p), ]
  W <- embed(yt, dimension = p+1)[, -(1:k)]
  X <- dt[-c(1:p), ,drop=F]
  Q <- embed(dt, dimension = p+1)[, -(1:q), drop=F]
  
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
                  dt=dt,
                  D=X)
  return(x)
}
