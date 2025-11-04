summary_bvar <- function(x, estimation = c("stan", "gibbs")) {
  estimation <- match.arg(estimation)
  
  if (estimation == "gibbs") {
    fit_obj <- x$fit$gibbs
    beta_posterior_mean <- fit_obj$beta_post_mean
    Psi_posterior_mean  <- fit_obj$Psi_post_mean
    Sigma_u_posterior_mean <- fit_obj$Sigma_u_post_mean
  } else {
    fit_obj <- x$fit$stan
    posterior <- rstan::extract(fit_obj)
    beta_posterior_mean <- apply(posterior$beta, c(2,3), mean)
    Psi_posterior_mean  <- apply(posterior$Psi, c(2,3), mean)
    Sigma_u_posterior_mean <- apply(posterior$Sigma_u, c(2,3), mean)
  }
  
  res <- list(
    beta_posterior_mean = round(beta_posterior_mean,2),
    Psi_posterior_mean = round(Psi_posterior_mean,2),
    Sigma_u_posterior_mean = round(Sigma_u_posterior_mean,2)
  )

  return(res)
}