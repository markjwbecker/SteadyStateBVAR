#' Summarise a fitted BVAR model
#'
#' Computes posterior means of the model parameters from a fitted \code{bvar}
#' object. The parameters returned depend on the model specification:
#' standard homoscedastic models return \code{beta}, \code{Psi}, and
#' \code{Sigma}; stochastic volatility models additionally return volatility
#' parameters.
#'
#' @param object A \code{bvar} object that has been passed through \code{\link{fit}}.
#' @param pars Character vector of parameter names to include. If \code{NULL}
#'   (default), all parameters are returned.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{summary.bvar}.
#' @export
#'
#' @examples
#' \dontrun{
#' yt <- matrix(rnorm(40, 0, 1), 20, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1)
#'
#' bvar_obj <- priors(bvar_obj,
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2))
#'
#' bvar_obj$predict$H <- 1
#' bvar_obj$predict$d_pred <- matrix(1)
#'
#' bvar_obj <- fit(bvar_obj,
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1,
#'                 auto_write = FALSE)
#'                 
#' summary(bvar_obj)
#' }
summary.bvar <- function(object, pars = NULL, ...) {
  
  summary.bvar <- function(object, pars = NULL, ...) {
    
    if (is.null(object$setup)) {
      stop("object must be passed through setup() first")
    }
    
    if (is.null(object$priors)) {
      stop("object must be passed through priors() first")
    }
    
    if (is.null(object$fit)) {
      stop("object must be passed through fit() first")
    }
    
    # Rest of summary function...
  }
  
  to_mat <- function(arr) {
    m <- apply(arr, c(2, 3), mean)
    matrix(m, nrow = nrow(m), ncol = ncol(m))
  }
  
  summaries <- list()
  
  keep_param <- function(lst) {
    if (is.null(pars)) return(lst)
    lst[names(lst) %in% pars]
  }
  
  if (is.null(object$SV)) {
    posterior <- rstan::extract(object$fit$stan)
    summaries$results <- keep_param(list(
      beta  = round(to_mat(posterior$beta), 2),
      Psi   = round(to_mat(posterior$Psi), 2),
      Sigma = round(to_mat(posterior$Sigma_u), 2)
    ))
  } else if (object$SV && object$SV_type == "AR") {
    posterior <- rstan::extract(object$fit$stan)
    summaries$results <- keep_param(list(
      beta    = round(to_mat(posterior$beta), 2),
      Psi     = round(to_mat(posterior$Psi), 2),
      A       = round(to_mat(posterior$A), 2),
      gamma_0 = round(apply(posterior$gamma_0, 2, mean), 2),
      gamma_1 = round(apply(posterior$gamma_1, 2, mean), 2),
      Phi     = round(to_mat(posterior$Phi), 2)
    ))
  } else if (object$SV && object$SV_type == "RW") {
    posterior <- rstan::extract(object$fit$stan)
    summaries$results <- keep_param(list(
      beta = round(to_mat(posterior$beta), 2),
      Psi  = round(to_mat(posterior$Psi), 2),
      A    = round(to_mat(posterior$A), 2),
      phi  = setNames(
        round(apply(posterior$phi, 2, mean), 2),
        paste0("phi_", 1:ncol(posterior$phi))
      )
    ))
  }
  
  out <- list(summaries = summaries, SV = object$SV, SV_type = object$SV_type)
  class(out) <- "summary.bvar"
  return(out)
}

#' @rdname summary.bvar
#' @param x A \code{summary.bvar} object returned by \code{summary.bvar}.
#' @export
print.summary.bvar <- function(x, ...) {
  
  s <- x$summaries$results
  
  for (param_name in setdiff(names(s), c("phi", "gamma_0", "gamma_1"))) {
    cat(param_name, "posterior mean\n")
    print(s[[param_name]])
    cat("\n")
  }
  
  if (!is.null(s$gamma_0)) {
    cat("gamma_0 posterior means\n")
    print(s$gamma_0)
    cat("\n")
  }
  
  if (!is.null(s$gamma_1)) {
    cat("gamma_1 posterior means\n")
    print(s$gamma_1)
    cat("\n")
  }
  
  if (!is.null(s$phi)) {
    cat("phi posterior means\n")
    for (nm in names(s$phi)) {
      cat(" ", nm, ":", s$phi[[nm]], "\n")
    }
    cat("\n")
  }
  
  invisible(x)
}