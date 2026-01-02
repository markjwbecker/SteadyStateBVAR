summary.bvar <- function(object) {
  
  has_stan  <- !is.null(object$fit$stan)
  has_gibbs <- !is.null(object$fit$gibbs)
  
  if (!has_stan && !has_gibbs)
    stop("No estimation results found in object$fit.")
  
  summaries <- list()
  
  if (has_stan) {
    fit <- object$fit$stan
    posterior <- rstan::extract(fit)
    summaries$stan <- list(
      method = "Stan",
      beta  = round(apply(posterior$beta,  c(2,3), mean), 2),
      Psi   = round(apply(posterior$Psi,   c(2,3), mean), 2),
      Sigma = round(apply(posterior$Sigma_u, c(2,3), mean), 2)
    )
  }
  
  if (has_gibbs) {
    fit <- object$fit$gibbs
    summaries$gibbs <- list(
      method = "Gibbs",
      beta  = round(fit$beta_post_mean, 2),
      Psi   = round(fit$Psi_post_mean, 2),
      Sigma = round(fit$Sigma_u_post_mean, 2)
    )
  }
  
  out <- list(summaries = summaries)
  class(out) <- "summary.bvar"
  return(out)
}


print.summary.bvar <- function(x) {
  
  both_methods <- length(x$summaries) > 1
  
  for (method_name in names(x$summaries)) {
    s <- x$summaries[[method_name]]
    
    if (both_methods) {
      cat("====================================\n")
      cat("Estimation Method:", s$method, "\n")
      cat("====================================\n\n")
    }
    
    cat("beta posterior mean\n"); print(s$beta); cat("\n")
    cat("Psi posterior mean\n"); print(s$Psi); cat("\n")
    cat("Sigma_u posterior mean\n"); print(s$Sigma); cat("\n\n")
  }
  
  invisible(x)
}
