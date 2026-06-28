#' Specify priors for the steady-state BVAR model
#'
#' The function prepares the priors. The Minnesota prior is used for the autoregressive parameters, and is determined by
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
#' @param first_own_lag_prior_mean Numeric vector of length \code{k}. Prior means for the first own lag
#'   of each variable. If \code{NULL}, defaults to a zero vector.
#' @param theta_Psi Numeric vector. Prior mean for vec(Psi), i.e. the steady-state parameters. If \code{NULL},
#'   defaults to the OLS estimates.
#' @param Omega_Psi Numeric matrix. Prior covariance matrix for vec(Psi), i.e. the steady-state parameters. If \code{NULL},
#'   defaults to a diagonal matrix with variances \code{1000}.
#' @param Jeffrey Logical. If \code{TRUE}, uses a Jeffreys prior for the innovation covariance matrix.
#'   If \code{FALSE}, uses an uninformative inverse-Wishart prior. Only considered if \code{SV=FALSE}.
#' @param SV Logical. If \code{TRUE}, enables stochastic volatility specification.
#'   Default \code{FALSE}.
#' @param SV_type Character. Type of stochastic volatility model. Must be either \code{"RW"} or \code{"AR1"}.
#'   Required if \code{SV = TRUE}.
#' @param SV_priors List. User-supplied stochastic volatility priors. Required when \code{SV = TRUE}.
#'   The list must contain the following named elements depending on \code{SV_type}:
#'   \itemize{
#'     \item For \code{"RW"}: \code{theta_A}, \code{Omega_A}, \code{mu_log_lambda_0},
#'       \code{sigma2_log_lambda_0}, \code{alpha_phi}, \code{beta_phi}.
#'     \item For \code{"AR1"}: \code{theta_A}, \code{Omega_A}, \code{theta_gamma_0},
#'       \code{Omega_gamma_0}, \code{theta_gamma_1}, \code{Omega_gamma_1},
#'       \code{theta_log_lambda_0}, \code{Omega_log_lambda_0}, \code{V_0}, \code{m_0}.
#'   }
#'
#' @return The steady-state \code{bvar} object with an appended \code{priors} list containing:
#'   \item{theta_beta}{Prior mean vector for vec(beta) constructed with the Minnesota prior}
#'   \item{Omega_beta}{Prior covariance matrix for vec(beta) constructed with the Minnesota prior}
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
#' @details
#' 
#' The goal is to estimate \eqn{\beta, \Psi}, and \eqn{\Sigma_u}, so priors are needed.
#' First, prior independence between \eqn{\beta, \Psi} and \eqn{\Sigma_u} is assumed.
#' 
#' Starting with \eqn{\beta}, the Minnesota prior is used
#'
#' \deqn{\mathrm{vec}(\beta) \sim \mathrm{N}_{kpk} \left[\theta_\beta,\Omega_\beta\right].}
#'
#' The prior means (the elements of \eqn{\theta_\beta}) are set to
#'
#' \deqn{
#' \begin{aligned}
#' \mathrm{E}\left(\Pi_{\ell}^{(i,j)}\right)&=
#' \begin{cases}\kappa & \text{if } \ell = 1 \ \text{and} \ i = j \\0 & \text{otherwise}\end{cases}\\
#' \kappa&=\begin{cases}\kappa^{level} & \text{if } \text{variable} \ i \ \text{is in level} \\
#' \kappa^{\Delta} & \text{if } \text{variable} \ i \ \text{is in difference}
#' \end{cases}\\
#' \end{aligned}
#' }
#'
#' Here, the autoregressive coefficient \eqn{\Pi_{\ell}^{(i,j)}} is element \eqn{\left(i,j\right)}
#' of \eqn{\Pi_{\ell}} for \eqn{\ell=1,\dots,p}. As such, the Minnesota prior sets all prior means
#' for the elements in \eqn{\beta} to \eqn{0}, except for the elements that relate to the first
#' own lags of the variables, which are set to \eqn{\kappa}. If variable \eqn{i} is in level
#' (e.g. nominal interest rate), then \eqn{\kappa=\kappa^{level}}, and typical choices for
#' \eqn{\kappa^{level}} are \eqn{1} or \eqn{0.9}. Evaluating the equations at their prior means,
#' equation \eqn{i} becomes a random walk if \eqn{\kappa^{level}=1} and a persistent stationary
#' AR(1) process if \eqn{\kappa^{level}=0.9}. Since the steady state only exists if the process
#' is stationary, \eqn{0.9} is recommended for the steady-state BVAR. If variable \eqn{i} is in difference
#' (e.g. output growth), then \eqn{\kappa=\kappa^{\Delta}}, and the most common choice for
#' \eqn{\kappa^{\Delta}} is \eqn{0}, i.e. equation \eqn{i} becomes (when evaluating it at its prior means)
#' a random walk expressed in first differences. If a differenced variable still shows some degree
#' of persistence (can be examined with an ACF plot), a suitable value for \eqn{\kappa^{\Delta}} can
#' be (for example) \eqn{0.6} instead of \eqn{0}. Moving on to the prior variances, \eqn{\Omega_\beta} is a
#' diagonal matrix containing the prior variances for the elements in \eqn{\beta}. They are specified as
#' 
#' \deqn{\mathrm{Var}\left(\Pi_{\ell}^{(i,j)}\right)=
#' \begin{cases}\left(\frac{\lambda_1}{\ell^{\lambda_3}}\right)^2 & \text{if } i = j \\
#' \left(\frac{\lambda_1 \lambda_2\sigma_i}{\ell^{\lambda_3}\sigma_j}\right)^2& \text{if } i \neq j
#' \end{cases}.}
#' 
#' Here \eqn{\lambda_1}, \eqn{\lambda_2}, and \eqn{\lambda_3} are scalar hyperparameters known as
#' the overall tightness, the cross-equation tightness and the lag decay rate. Furthermore,
#' \eqn{\sigma_i^2} is the \eqn{(i,i)}:th element of \eqn{\Sigma_u}, which is unknown and therefore
#' replaced with an estimate. In this package, it is replaced by the least squares residual variance
#' from a univariate autoregression for variable \eqn{i} with \eqn{p} lags
#' (including the constant and dummy/trend variable if applicable). Moving on to \eqn{\Psi}, the
#' prior is
#' 
#' \deqn{\mathrm{vec}(\Psi) \sim \mathrm{N}_{kq}\left[\theta_\Psi,\Omega_\Psi\right].}
#' 
#' This is the core of the steady-state BVAR model.
#' In \eqn{\theta_\Psi}, one specifies the prior beliefs about the location of the steady state,
#' and in \eqn{\Omega_\Psi}, which is assumed to be a diagonal matrix, one specifies the degree
#' of certainty in those prior beliefs. The prior for \eqn{\Sigma_u} is either the usual noninformative Jeffreys prior
#' 
#' \deqn{p(\Sigma_u) \propto\left|\Sigma_u \right|^{-(k+1)/2},}
#' 
#' or an inverse-Wishart prior
#' 
#' \deqn{\Sigma_u \sim \mathrm{IW}(V_0,m_0),}
#' 
#' where \eqn{V_0} is the scale matrix and \eqn{m_0\geq k+2} is the number of degrees of freedom.
#' An uninformative prior can be (and is in this package) specified by setting
#' \eqn{V_0=(m_0-k-1)\hat{\Sigma}_u} where \eqn{\hat{\Sigma}_u} is the least squares estimate
#' from the VAR(\eqn{p}) (including the constant and dummy/trend variable if applicable), and \eqn{m_0=k+2}.
#' For the stochastic volatility specifications, the innovation covariance matrix is now time-varying \eqn{\Sigma_{u,t}}.
#' Therefore, stochastic volatility priors are needed, see \link{bvar} for more details. For the Random Walk
#' (\code{RW}) stochastic volatility specification, the following priors are used
#' 
#' \deqn{\begin{aligned}a &\sim \mathrm{N}(\theta_A, \Omega_A) \\
#' \ln \lambda_{i,0} &\sim \mathrm{N}(\mu_{\ln \lambda_{i,0}}, \sigma^2_{\ln \lambda_{i,0}}) \\
#' \phi_i &\sim \mathrm{IG}(\alpha_{\phi_i},\beta_{\phi_i})\end{aligned}.}
#' 
#' Here \eqn{a} is a \eqn{k(k-1)/2} vector that collects the free parameters in \eqn{A} in row-major order,
#' and \eqn{\ln \lambda_{i,0}} are the time \eqn{t=0} values (initial conditions) of \eqn{\ln \lambda_{i,t}}.
#' Furthermore, \eqn{\phi_i} are the log volatility innovation variances. Note here that the inverse gamma
#' \eqn{\mathrm{IG}(\alpha, \beta)} distribution is the univariate version of the inverse-Wishart distribution
#' with \eqn{\alpha=m/2, \ \beta = V/2}. For the AR(1) (\code{AR1}) stochastic volatility specification, the following priors are used
#' 
#' \deqn{\begin{aligned}a &\sim \mathrm{N}(\theta_A, \Omega_A) \\
#' \gamma_{0} &\sim \mathrm{N}(\theta_{\gamma_0}, \Omega_{\gamma_0}) \\
#' \gamma_{1} &\sim \mathrm{N}(\theta_{\gamma_1}, \Omega_{\gamma_1}) \\
#' \ln \lambda_{0} &\sim \mathrm{N}(\theta_{\ln \lambda_{0}}, \Omega_{\ln \lambda_{0}}) \\
#' \Phi &\sim \mathrm{IW}(V_0,m_0)\end{aligned}.}
#' 
#' Here \eqn{a} is again the \eqn{k(k-1)/2} vector that collects the free parameters in \eqn{A} in row-major order,
#' and \eqn{\ln \lambda_0} are the time \eqn{t=0} values (initial conditions) of \eqn{\ln \lambda_{t}}.
#' Furthermore, \eqn{\gamma_{0}} are the log volatility intercepts, \eqn{\gamma_{1}} are the log volatility
#' slopes, and \eqn{\Phi} is the log volatility innovation covariance matrix.
#' 
#' @export
#'
#' @examples
#' #homoscedastic with Jeffreys prior
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1)
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    Jeffrey = TRUE,
#'                    SV = FALSE,
#'                    SV_type = NULL,
#'                    SV_priors = NULL)
#'                    
#' #RW stochastic volatility
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1)
#' 
#' k <- bvar_obj$setup$k
#' n_free_params_A <- bvar_obj$setup$n_free_params_A
#' 
#' SV_priors_RW <- list(
#' theta_A             =  rep(0, n_free_params_A),
#' Omega_A             =  diag(1000, n_free_params_A),
#' mu_log_lambda_0     =  rep(0, k),
#' sigma2_log_lambda_0 =  rep(1000, k),
#' alpha_phi           =  rep(5, k),
#' beta_phi            = (rep(5, k) - 1) * rep(0.1, k)
#' )
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    SV = TRUE,
#'                    SV_type = "RW",
#'                    SV_priors = SV_priors_RW)
#'                    
#' #AR1 stochastic volatility
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1)
#' 
#' k <- bvar_obj$setup$k
#' n_free_params_A <- bvar_obj$setup$n_free_params_A
#' 
#' SV_priors_AR <- list(
#' theta_A            =  rep(0, n_free_params_A),
#' Omega_A            =  diag(1000, n_free_params_A),
#' theta_gamma_0      =  rep(0.1, k),
#' Omega_gamma_0      =  diag(1000, k),
#' theta_gamma_1      =  rep(0.9, k),
#' Omega_gamma_1      =  diag(10, k),
#' theta_log_lambda_0 =  rep(0.1, k)/(1-rep(0.9, k)),
#' Omega_log_lambda_0 =  diag(1000, k),
#' V_0                = (10 - k - 1) * diag(k),
#' m_0                =  10
#' )
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    SV = TRUE,
#'                    SV_type = "AR1",
#'                    SV_priors = SV_priors_AR)
#'
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
    
    n_free <- x$setup$n_free_params_A
    
    if (SV_type == "RW") {
      required_names <- c("theta_A", "Omega_A", "mu_log_lambda_0",
                          "sigma2_log_lambda_0", "alpha_phi", "beta_phi")
      missing_names <- setdiff(required_names, names(SV_priors))
      if (length(missing_names) > 0)
        stop(paste("SV_priors is missing elements:", paste(missing_names, collapse = ", ")))
      
      if (length(SV_priors$theta_A) != n_free)
        stop(paste("theta_A must be a vector of length k*(k-1)/2 =", n_free))
      if (!all(dim(SV_priors$Omega_A) == n_free))
        stop(paste("Omega_A must be a", n_free, "x", n_free, "matrix"))
      if (length(SV_priors$mu_log_lambda_0) != k)
        stop(paste("mu_log_lambda_0 must be a vector of length k =", k))
      if (length(SV_priors$sigma2_log_lambda_0) != k)
        stop(paste("sigma2_log_lambda_0 must be a vector of length k =", k))
      if (any(SV_priors$sigma2_log_lambda_0 <= 0))
        stop("sigma2_log_lambda_0 must be strictly positive")
      if (length(SV_priors$alpha_phi) != k)
        stop(paste("alpha_phi must be a vector of length k =", k))
      if (any(SV_priors$alpha_phi <= 0))
        stop("alpha_phi must be strictly positive")
      if (length(SV_priors$beta_phi) != k)
        stop(paste("beta_phi must be a vector of length k =", k))
      if (any(SV_priors$beta_phi <= 0))
        stop("beta_phi must be strictly positive")
      
    } else if (SV_type == "AR1") {
      required_names <- c("theta_A", "Omega_A", "theta_gamma_0", "Omega_gamma_0",
                          "theta_gamma_1", "Omega_gamma_1", "theta_log_lambda_0",
                          "Omega_log_lambda_0", "m_0", "V_0")
      missing_names <- setdiff(required_names, names(SV_priors))
      if (length(missing_names) > 0)
        stop(paste("SV_priors is missing elements:", paste(missing_names, collapse = ", ")))
      
      if (length(SV_priors$theta_A) != n_free)
        stop(paste("theta_A must be a vector of length k*(k-1)/2 =", n_free))
      if (!all(dim(SV_priors$Omega_A) == n_free))
        stop(paste("Omega_A must be a", n_free, "x", n_free, "matrix"))
      if (length(SV_priors$theta_gamma_0) != k)
        stop(paste("theta_gamma_0 must be a vector of length k =", k))
      if (!all(dim(SV_priors$Omega_gamma_0) == k))
        stop(paste("Omega_gamma_0 must be a", k, "x", k, "matrix"))
      if (length(SV_priors$theta_gamma_1) != k)
        stop(paste("theta_gamma_1 must be a vector of length k =", k))
      if (!all(dim(SV_priors$Omega_gamma_1) == k))
        stop(paste("Omega_gamma_1 must be a", k, "x", k, "matrix"))
      if (length(SV_priors$theta_log_lambda_0) != k)
        stop(paste("theta_log_lambda_0 must be a vector of length k =", k))
      if (!all(dim(SV_priors$Omega_log_lambda_0) == k))
        stop(paste("Omega_log_lambda_0 must be a", k, "x", k, "matrix"))
      if (!is.numeric(SV_priors$m_0) || length(SV_priors$m_0) != 1 || SV_priors$m_0 < k)
        stop(paste("m_0 must be a scalar integer >= k =", k))
      if (!all(dim(SV_priors$V_0) == k))
        stop(paste("V_0 must be a", k, "x", k, "matrix"))
    }
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