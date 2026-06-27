#' Summarise a fitted steady-state BVAR model
#'
#' Computes and prints posterior summaries from a fitted \code{bvar} object.
#' The printed output depends on whether the model is homoscedastic or includes
#' stochastic volatility (RW or AR1 specification).
#'
#' Parameters are printed in structured blocks. Matrix-valued parameters are
#' displayed as matrices, while vector-valued stochastic volatility parameters
#' (e.g. \code{phi}, \code{gamma_0}, \code{gamma_1}) are printed as named vectors.
#'
#' For stochastic volatility models, time-varying covariance matrices
#' \code{Sigma_u,t} can be extracted at a specified time index.
#'
#' @param object A \code{bvar} object that has been passed through
#'   \code{\link{fit}}.
#' @param pars Character vector of parameter names to include. If
#'   \code{NULL} (default), all available parameters are displayed.
#' @param stat Character. Posterior summary statistic to display:
#'   \code{"mean"} (default) or \code{"median"}.
#' @param t Integer. Time index for stochastic volatility covariance matrices.
#'   If \code{NULL}, the latest available period is used.
#' @param ... Additional arguments (currently unused).
#'
#' @return Returns the input object invisibly.
#'
#' @details
#' The function supports both standard VAR models and stochastic volatility
#' extensions. For SV models, additional parameters may be displayed depending
#' on specification:
#' \itemize{
#'   \item RW  SV: \code{phi}
#'   \item AR1 SV: \code{gamma_0}, \code{gamma_1}, \code{Phi}
#' }
#'
#' Output is printed in blocks with manual formatting for readability.
#'
#' @exportS3Method
summary.bvar <- function(object, pars = NULL, stat = "mean", t = NULL, ...) {
  
  if (is.null(object$setup)) stop("object must be passed through setup() first")
  if (is.null(object$priors)) stop("object must be passed through priors() first")
  if (is.null(object$fit)) stop("object must be passed through fit() first")
  
  stat <- match.arg(stat, c("mean", "median"))
  posterior <- rstan::extract(object$fit$stan)
  
  var_names <- colnames(object$data)
  if (is.null(var_names)) {
    var_names <- paste0("Var", 1:ncol(object$data))
  }
  
  is_sv <- !is.null(object$priors$SV) && object$priors$SV
  
  compute_stat <- function(arr, stat) {
    apply(arr, c(2, 3), if (stat == "mean") mean else median)
  }
  
  results <- list()
  
  # ---------------- beta ----------------
  if (!is.null(posterior$beta)) {
    beta <- compute_stat(posterior$beta, stat)
    
    var_rep <- rep(var_names, object$setup$p)
    lag_names <- paste0("l", rep(1:object$setup$p, each = object$setup$k))
    
    rownames(beta) <- paste0(var_rep, ".", lag_names)
    colnames(beta) <- var_names
    
    results$beta <- round(beta, 2)
  }
  
  # ---------------- Psi ----------------
  if (!is.null(posterior$Psi)) {
    Psi <- compute_stat(posterior$Psi, stat)
    rownames(Psi) <- var_names
    results$Psi <- round(Psi, 2)
  }
  
  # ---------------- Sigma_u ----------------
  if (!is.null(posterior$Sigma_u)) {
    
    if (is_sv) {
      
      N_original <- nrow(object$data)
      N_est <- N_original - object$setup$p
      
      if (is.null(t)) {
        t_est <- N_est
        t_display <- N_original
      } else {
        if (t < object$setup$p + 1 || t > N_original) {
          stop("t must be between ", object$setup$p + 1, " and ", N_original)
        }
        t_est <- t - object$setup$p
        t_display <- t
      }
      
      Sigma <- posterior$Sigma_u[, t_est, , ]
      Sigma <- compute_stat(Sigma, stat)
      
      rownames(Sigma) <- var_names
      colnames(Sigma) <- var_names
      
      results$`Sigma_u,t` <- round(Sigma, 2)
      
    } else {
      
      Sigma <- compute_stat(posterior$Sigma_u, stat)
      rownames(Sigma) <- var_names
      colnames(Sigma) <- var_names
      
      results$Sigma_u <- round(Sigma, 2)
    }
  }
  
  # ---------------- A ----------------
  if (!is.null(posterior$A)) {
    A <- compute_stat(posterior$A, stat)
    rownames(A) <- var_names
    colnames(A) <- var_names
    results$A <- round(A, 2)
  }
  
  # ---------------- RW SV ----------------
  if (!is.null(posterior$phi)) {
    phi <- apply(posterior$phi, 2, if (stat == "mean") mean else median)
    names(phi) <- var_names
    results$phi <- round(phi, 2)
  }
  
  # ---------------- AR1 SV ----------------
  if (!is.null(posterior$gamma_0)) {
    gamma_0 <- apply(posterior$gamma_0, 2, if (stat == "mean") mean else median)
    names(gamma_0) <- var_names
    results$gamma_0 <- round(gamma_0, 2)
  }
  
  if (!is.null(posterior$gamma_1)) {
    gamma_1 <- apply(posterior$gamma_1, 2, if (stat == "mean") mean else median)
    names(gamma_1) <- var_names
    results$gamma_1 <- round(gamma_1, 2)
  }
  
  if (!is.null(posterior$Phi)) {
    Phi <- compute_stat(posterior$Phi, stat)
    rownames(Phi) <- var_names
    colnames(Phi) <- var_names
    results$Phi <- round(Phi, 2)
  }
  
  # ---------------- filter ----------------
  if (!is.null(pars)) {
    results <- results[names(results) %in% pars]
  }
  
  # ---------------- print ----------------
  if (stat == "mean") {
    cat("Posterior", stat, "estimates\n------------------------\n\n")
  } else {
    cat("Posterior", stat, "estimates\n--------------------------\n\n")
  }
  
  for (i in seq_along(results)) {
    
    param_name <- names(results)[i]
    obj <- results[[i]]
    
    # -------------------------
    # Sigma_u,t
    # -------------------------
    if (param_name == "Sigma_u,t") {
      
      cat("Sigma_u,t (t = ", t_display, ")\n\n", sep = "")
      
      print(obj)
      
      # -------------------------
      # beta
      # -------------------------
    } else if (param_name == "beta") {
      
      cat("beta\n")
      print(obj)
      
      # -------------------------
      # Psi
      # -------------------------
    } else if (param_name == "Psi") {
      
      cat("Psi\n")
      print(obj)
      
      # -------------------------
      # A
      # -------------------------
    } else if (param_name == "A") {
      
      cat("A\n")
      print(obj)
      
      # -------------------------
      # phi (RW SV)
      # -------------------------
    } else if (param_name == "phi") {
      
      cat("phi\n\n")
      print(obj)
      
      # -------------------------
      # gamma_0 (AR1 SV)
      # -------------------------
    } else if (param_name == "gamma_0") {
      
      cat("gamma_0\n\n")
      print(obj)
      
      # -------------------------
      # gamma_1 (AR1 SV)
      # -------------------------
    } else if (param_name == "gamma_1") {
      
      cat("gamma_1\n\n")
      print(obj)
      
      # -------------------------
      # Phi (AR1 SV)
      # -------------------------
    } else if (param_name == "Phi") {
      
      cat("Phi\n")
      print(obj)
      
      # -------------------------
      # fallback (future-proofing)
      # -------------------------
    } else {
      
      cat(param_name, "\n\n")
      print(obj)
    }
    
    if (i < length(results)) cat("\n")
  }
  
  invisible(object)
}