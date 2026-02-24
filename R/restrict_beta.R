restrict_beta <- function(x, restriction_matrix) {
  
  k <- x$setup$k
  p <- x$setup$p
  
  if (!all(dim(restriction_matrix) == c(k * p, k))) {
    stop("restriction_matrix must have dimension (k*p x k)")
  }
  x$setup$restriction_matrix <- restriction_matrix
  
  zero_indices <- which(c(restriction_matrix) == 0)
  
  if (!is.null(x$priors$Omega_beta)) {
    diag(x$priors$Omega_beta[zero_indices, zero_indices]) <- 0.00001
  } else {
    warning("Omega_beta not found in priors: restriction applied but Omega_beta not updated")
  }
  return(x)
}
