#' Set up the steady-state BVAR model
#'
#' Prepares the matrices needed for prior specification and estimation. 
#' Also computes OLS estimates.
#'
#' @param x A steady-state \code{bvar} object created by \code{\link{bvar}}.
#' @param p Integer. The lag order of the VAR.
#' @param deterministic Character. The deterministic component to include.
#'   One of \code{"constant"} (default), \code{"constant_and_dummy"}, or
#'   \code{"constant_and_trend"}.
#' @param dummy Numeric vector of a dummy variable. Only
#'   used when \code{deterministic = "constant_and_dummy"}. Default \code{NULL}.
#'
#' @return The \code{bvar} object with a \code{setup} list containing the
#' matrices required for prior specification and estimation, and also the OLS estimates.
#' @export
#'
#' @examples
#' yt <- matrix(rnorm(50), 25, 2)
#' 
#' bvar_obj <- bvar(data = yt)
#' 
#' bvar_obj <- setup(bvar_obj, p = 1)
setup <- function(x, p, deterministic=c("constant", "constant_and_dummy", "constant_and_trend"), dummy=NULL) {
  
  if (!inherits(x, "bvar")) {
    stop("x must be a 'bvar' object")
  }
  
  if (is.null(x$data)) {
    stop("x must have data. Create bvar object with bvar(data = ...)")
  }
  
  if (!is.numeric(p) || p <= 0 || p != floor(p)) {
    stop("p must be a positive integer")
  }
  
  yt <- x$data
  
  if (p >= nrow(yt)) {
    stop("p must be less than the number of observations")
  }
  
  deterministic <- match.arg(deterministic)
  
  if (deterministic == "constant_and_dummy" && is.null(dummy)) {
    stop("dummy must be provided when deterministic = 'constant_and_dummy'")
  }
  
  if (!is.null(dummy) && length(dummy) != nrow(yt)) {
    stop("dummy must have length equal to number of observations")
  }
  
  if (!is.null(dummy) && !is.vector(dummy)) {
    stop("dummy must be a vector")
  }
  
  
  N = nrow(yt) - p
  k = ncol(yt)
  
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
  X <- dt[-c(1:p), , drop=F]
  Q <- embed(dt, dimension = p+1)[, -(1:q), drop=F]
  
  Z <- cbind(W, X)
  beta_hat = solve(crossprod(Z), crossprod(Z, Y))
  U = Y - Z %*% beta_hat
  Sigma_u_OLS <- crossprod(U) / (N - k*p - q)
  
  if (q == 1) {
    C_hat <- beta_hat[(k*p+1):(k*p+q), ]
  } else {
    C_hat <- t(beta_hat[(k*p+1):(k*p+q), ])
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
  beta_OLS = beta_hat[1:(k*p), ]
  
  n_free_params_A <- k*(k-1)/2
  
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
                  D=X,
                  n_free_params_A=n_free_params_A)
  return(x)
}