#' Specify priors for a steady-state BVAR model
#'
#' The Minnesota prior is used for the autoregressive parameters, and is determined by
#' the overall tightness, cross-equation tightness, and the lag decay rate.
#' For the steady-state parameters, a normal prior is used. For the covariance matrix of the innovations,
#' the user can choose between Jeffreys prior or an uninformative inverse-Wishart prior.
#' Optionally enables a stochastic volatility specification for the covariance matrix of the innovations (random walk or AR(1)).
#'
#' @param x A steady-state \code{bvar} object that has been passed through \code{\link{setup}}.
#' @param lambda_1 Numeric. Overall tightness of the Minnesota prior.
#'   Default \code{0.2}.
#' @param lambda_2 Numeric. Cross-equation tightness of the Minnesota prior. Default \code{0.5}.
#' @param lambda_3 Numeric. Lag decay rate of the Minnesota prior. Default \code{1}.
#' @param first_own_lag_prior_mean Numeric vector of length \code{k}. Prior means for the first own lags
#'   of each variable. If \code{NULL}, defaults to a zero vector.
#' @param theta_Psi Numeric vector. Prior mean vector for steady-state parameters. If \code{NULL},
#'   defaults to OLS estimate.
#' @param Omega_Psi Numeric matrix. Prior covariance matrix for steady-state parameters. If \code{NULL},
#'   defaults to a diagonal matrix with variances 1000.
#' @param Jeffrey Logical. If \code{TRUE}, uses a Jeffreys prior for the innovation covariance matrix.
#'   If \code{FALSE}, uses an uninformative inverse-Wishart prior.
#' @param SV Logical. If \code{TRUE}, enables stochastic volatility specification.
#'   Default \code{FALSE}.
#' @param SV_type Character. Type of stochastic volatility model. Must be either \code{"RW"} or \code{"AR1"}.
#'   Required if \code{SV = TRUE}.
#' @param SV_priors List. User-supplied stochastic volatility priors.
#'   Required when \code{SV = TRUE}.
#'
#' @return The \code{bvar} object with an appended \code{priors} list containing:
#'   \item{theta_beta}{Prior mean for vec(beta)}
#'   \item{Omega_beta}{Prior covariance matrix for vec(beta)}
#'   \item{theta_Psi}{Prior mean for vec(Psi), i.e. the steady-state parameters}
#'   \item{Omega_Psi}{Prior covariance matrix for vec(Psi), i.e. the steady-state parameters}
#'   \item{Jeffrey}{Indicator for Jeffreys prior usage}
#'   \item{Sigma_AR}{Residual variance estimates from univariate AR fits, which are used by the Minnesota prior}
#'   \item{m_0}{Inverse-Wishart prior degrees of freedom (if \code{Jeffrey = FALSE})}
#'   \item{V_0}{Inverse-Wishart prior scale matrix (if \code{Jeffrey = FALSE})}
#'   \item{SV}{Logical indicator for stochastic volatility specification}
#'   \item{SV_type}{Stochastic volatility specification type}
#'   \item{SV_priors}{User-supplied SV prior list (if \code{SV = TRUE})}
#'
#' @export
#'
#' @examples
#' yt <- matrix(rnorm(40, 0, 1), 20, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1)
#'
#' bvar_obj <- priors(bvar_obj,
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2))
priors<- function(x,
                  lambda_1=0.2,
                  lambda_2=0.5,
                  lambda_3 = 1,
                  first_own_lag_prior_mean=NULL,
                  theta_Psi=NULL,
                  Omega_Psi=NULL,
                  Jeffrey=TRUE,
                  SV = FALSE,
                  SV_type = NULL,
                  SV_priors = NULL){
  
  if (!inherits(x, "bvar")) stop("must be a 'bvar' object")
  if (is.null(x$setup)) stop("must be passed through setup")
  
  setup <- x$setup
  
  if (!all(c("k","p","q","dt") %in% names(setup))) {
    stop("setup is incomplete")
  }
  
  k <- setup$k
  p <- setup$p
  q <- setup$q
  dt <- setup$dt
  yt <- x$data
  
  if (!is.null(theta_Psi) && !is.null(Omega_Psi)) {
    if (nrow(Omega_Psi) != length(theta_Psi)) {
      stop("dimensions must match")
    }
  }
  
  if (any(c(lambda_1, lambda_2, lambda_3) <= 0)) {
    stop("must be positive")
  }
  
  if (isTRUE(SV)) {
    if (is.null(SV_type)) stop("SV_type is needed")
    if (!SV_type %in% c("RW", "AR1")) stop("SV_type needs to be RW or AR1")
    if (is.null(SV_priors)) stop("SV_priors is required")
  }
  
  priors <- list()
  
  Sigma_AR <- diag(0, k)
  
  for (i in 1:k){
    
    y <- yt[,i]
    
    Y <- y[-c(1:p)]
    W <- embed(y, dimension = p+1)[, -1]
    X <- dt[-c(1:p), ,drop=F]
    
    Z <- cbind(W,X)
    beta_hat = solve(crossprod(Z,Z),crossprod(Z,Y))
    U = Y-Z%*%beta_hat
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
  
  if (!is.null(first_own_lag_prior_mean)) {
    if (length(first_own_lag_prior_mean) != k) {
      stop("must have length k")
    }
  }
  
  if (is.null(first_own_lag_prior_mean)) first_own_lag_prior_mean <- rep(0,k)
  
  mat <- matrix(0, nrow = k*p, ncol = k)
  for (i in 1:k){
    mat[i,i] <- first_own_lag_prior_mean[i]
  }
  theta_beta = c(mat)
  
  if (is.null(theta_Psi)) theta_Psi <- c(x$setup$Psi_OLS)
  if (is.null(Omega_Psi)) Omega_Psi <- diag(1000, k*q, k*q)
  
  if (isFALSE(Jeffrey)){
    m_0=k+2
    V_0 = (m_0-k-1)*setup$Sigma_u_OLS
    priors$V_0 <- V_0
    priors$m_0 <- m_0
  }
  
  priors$theta_beta <- theta_beta
  priors$Omega_beta <- Omega_beta
  priors$theta_Psi <- theta_Psi
  priors$Omega_Psi <- Omega_Psi
  priors$Jeffrey <- Jeffrey
  priors$Sigma_AR <- Sigma_AR
  
  priors$SV <- SV
  priors$SV_type <- SV_type
  
  if (isTRUE(SV)) {
    priors$SV_priors <- SV_priors
  }
  
  x$priors <- priors
  return(x)
}