summary.bvar <- function(x, pars = NULL) {
  
  has_stan  <- !is.null(x$fit$stan)
  has_gibbs <- !is.null(x$fit$gibbs)
  
  if (!has_stan && !has_gibbs)
    stop("No estimation results found in bvar_object$fit.")
  
  summaries <- list()
  
  keep_param <- function(lst) {
    if (is.null(pars)) return(lst)
    lst[names(lst) %in% pars]
  }
  
  if (has_stan && is.null(x$SV)) {
    fit <- x$fit$stan
    posterior <- rstan::extract(fit)
    summaries$stan <- keep_param(list(
      method = "Stan",
      beta  = round(apply(posterior$beta,  c(2,3), mean), 2),
      Psi   = round(apply(posterior$Psi,   c(2,3), mean), 2),
      Sigma = round(apply(posterior$Sigma_u, c(2,3), mean), 2)
    ))
  } else if (has_stan && x$SV){
    fit <- x$fit$stan
    posterior <- rstan::extract(fit)
    summaries$stan <- keep_param(list(
      method = "Stan",
      beta     = round(apply(posterior$beta,  c(2,3), mean), 2),
      Psi      = round(apply(posterior$Psi,   c(2,3), mean), 2),
      A        = round(apply(posterior$A,     c(2,3), mean), 2),
      gamma_0  = round(apply(posterior$gamma_0, 2, mean), 2),
      gamma_1  = round(apply(posterior$gamma_1, 2, mean), 2),
      Phi      = round(apply(posterior$Phi,   c(2,3), mean), 2)
    ))
  }
  
  if (has_gibbs) {
    fit <- x$fit$gibbs
    summaries$gibbs <- keep_param(list(
      method = "Gibbs",
      beta  = round(fit$beta_posterior_mean, 2),
      Psi   = round(fit$Psi_posterior_mean, 2),
      Sigma = round(fit$Sigma_u_posterior_mean, 2)
    ))
  }
  
  out <- list(summaries = summaries,SV = x$SV)
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
    
    for (param_name in setdiff(names(s), "method")) {
      cat(param_name, "posterior mean\n"); print(s[[param_name]]); cat("\n")
    }
  }
  
  invisible(x)
}
