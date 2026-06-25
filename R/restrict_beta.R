#' Restrict VAR coefficients to zero
#'
#' Applies zero restrictions to the VAR coefficient matrix by setting the
#' corresponding prior variances in \code{Omega_beta} to a value near zero.
#' This enforces restrictions through the prior rather than hard-coding them
#' in the likelihood, which is compatible with the Stan estimation.
#'
#' @param x A \code{bvar} object that has been passed through \code{\link{priors}}.
#' @param restriction_matrix A numeric matrix of dimension \code{k*p x k} where
#'   entries of \code{0} indicate coefficients to be restricted to zero and
#'   entries of \code{1} indicate unrestricted coefficients.
#'
#' @return The \code{bvar} object with the restriction matrix stored in
#'   \code{setup} and \code{Omega_beta} updated accordingly.
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
#'                    
#' p <- bvar_obj$setup$p
#' k <- bvar_obj$setup$k
#' restriction_matrix <- matrix(1, k*p, k)
#' restriction_matrix[1, 1] <- 0
#' 
#' bvar_obj <- restrict_beta(bvar_obj, restriction_matrix)
restrict_beta <- function(x, restriction_matrix) {
  
  k <- x$setup$k
  p <- x$setup$p
  
  if (!all(dim(restriction_matrix) == c(k * p, k))) {
    stop("restriction_matrix must have dimension (k*p x k)")
  }
  
  x$setup$restriction_matrix <- restriction_matrix
  zero_indices <- which(c(restriction_matrix) == 0)
  
  if (!is.null(x$priors$Omega_beta)) {
    diag(x$priors$Omega_beta)[zero_indices] <- 0.0000001
  } else {
    warning("Omega_beta not found in priors: restriction applied but Omega_beta not updated")
  }
  
  return(x)
}