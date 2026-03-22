summary.bvar <- function(x, pars = NULL) {
  
  has_stan  <- !is.null(x$fit$stan)
  has_gibbs <- !is.null(x$fit$gibbs)
  
  if (!has_stan && !has_gibbs)
    stop("No estimation results found in bvar_object$fit.")
  
  to_mat <- function(arr) {
    m <- apply(arr, c(2,3), mean)
    matrix(m, nrow = nrow(m), ncol = ncol(m))
  }
  
  summaries <- list()
  
  keep_param <- function(lst) {
    if (is.null(pars)) return(lst)
    lst[names(lst) %in% pars]
  }
  
  if (has_stan && is.null(x$SV)) {
    posterior <- rstan::extract(x$fit$stan)
    summaries$results <- keep_param(list(
      method = "Stan",
      beta  = round(to_mat(posterior$beta), 2),
      Psi   = round(to_mat(posterior$Psi), 2),
      Sigma = round(to_mat(posterior$Sigma_u), 2)
    ))
  } else if (has_stan && x$SV && x$SV_type == "AR") {
    posterior <- rstan::extract(x$fit$stan)
    summaries$results <- keep_param(list(
      method  = "Stan",
      beta    = round(to_mat(posterior$beta), 2),
      Psi     = round(to_mat(posterior$Psi), 2),
      A       = round(to_mat(posterior$A), 2),
      gamma_0 = round(apply(posterior$gamma_0, 2, mean), 2),
      gamma_1 = round(apply(posterior$gamma_1, 2, mean), 2),
      Phi     = round(to_mat(posterior$Phi), 2)
    ))
  } else if (has_stan && x$SV && x$SV_type == "RW") {
    posterior <- rstan::extract(x$fit$stan)
    summaries$results <- keep_param(list(
      method = "Stan",
      beta   = round(to_mat(posterior$beta), 2),
      Psi    = round(to_mat(posterior$Psi), 2),
      A      = round(to_mat(posterior$A), 2),
      phi    = setNames(
        round(apply(posterior$phi, 2, mean), 2),
        paste0("phi_", 1:ncol(posterior$phi))
      )
    ))
  }
  
  if (has_gibbs) {
    fit <- x$fit$gibbs
    summaries$results <- keep_param(list(
      method = "Gibbs",
      beta  = round(fit$beta_posterior_mean, 2),
      Psi   = round(fit$Psi_posterior_mean, 2),
      Sigma = round(fit$Sigma_u_posterior_mean, 2)
    ))
  }
  
  out <- list(summaries = summaries, SV = x$SV, SV_type = x$SV_type)
  class(out) <- "summary.bvar"
  return(out)
}

print.summary.bvar <- function(x) {
  
  s <- x$summaries$results
  
  for (param_name in setdiff(names(s), c("method", "phi", "gamma_0", "gamma_1"))) {
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