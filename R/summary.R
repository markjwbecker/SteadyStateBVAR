summary.bvar <- function(object) {
  
  has_stan  <- !is.null(object$fit$stan)
  has_gibbs <- !is.null(object$fit$gibbs)
  
  if (!has_stan && !has_gibbs)
    stop("No estimation results found in object$fit.")
  
  summaries <- list()
  
  if (has_stan && is.null(object$SV)) {
    fit <- object$fit$stan
    posterior <- rstan::extract(fit)
    summaries$stan <- list(
      method = "Stan",
      beta  = round(apply(posterior$beta,  c(2,3), mean), 2),
      Psi   = round(apply(posterior$Psi,   c(2,3), mean), 2),
      Sigma = round(apply(posterior$Sigma_u, c(2,3), mean), 2)
    )
  } else if (has_stan && object$SV){
    fit <- object$fit$stan
    posterior <- rstan::extract(fit)
    summaries$stan <- list(
      method = "Stan",
      beta     = round(apply(posterior$beta,  c(2,3), mean), 2),
      Psi      = round(apply(posterior$Psi,   c(2,3), mean), 2),
      A        = round(apply(posterior$A,     c(2,3), mean), 2),
      gamma_0  = round(apply(posterior$gamma_0, 2, mean), 2),
      gamma_1  = round(apply(posterior$gamma_1, 2, mean), 2),
      Phi      = round(apply(posterior$Phi,   c(2,3), mean), 2)
    )
  }
  
  if (has_gibbs) {
    fit <- object$fit$gibbs
    summaries$gibbs <- list(
      method = "Gibbs",
      beta  = round(fit$beta_posterior_mean, 2),
      Psi   = round(fit$Psi_posterior_mean, 2),
      Sigma = round(fit$Sigma_u_posterior_mean, 2)
    )
  }
  
  out <- list(summaries = summaries,
              SV = object$SV)
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
    if (is.null(x$SV)) {
      cat("Sigma_u posterior mean\n"); print(s$Sigma); cat("\n")
    } else {
      cat("A posterior mean\n"); print(s$A); cat("\n")
      cat("gamma_0 posterior mean\n"); print(s$gamma_0); cat("\n")
      cat("gamma_1 posterior mean\n"); print(s$gamma_1); cat("\n")
      cat("Phi posterior mean\n"); print(s$Phi); cat("\n")
    }
  }
  
  invisible(x)
}
