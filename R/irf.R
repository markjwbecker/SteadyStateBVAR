IRF <- function(x, H = 16, response, shock, type = c("median", "mean"), method=c("OIRF", "GIRF"), ci=0.95, estimation=c("stan","gibbs"), t=NULL) {
  
  method     <- match.arg(method)
  estimation <- match.arg(estimation)
  type <- match.arg(type)
  
  
  compute_OIRF <- function(A, Sigma, e, N){
    k <- nrow(Sigma)
    OIRF_matrix <- matrix(NA, k*k, N+1)
    count <- 1
    P <- t(chol(Sigma))
    
    for(j in 1:k){
      for(i in 1:k){
        for(n in 0:N){
          OIRF <- A[[n+1]] %*% P %*% e[[j]]   
          OIRF_matrix[count, n+1] <- OIRF[i]
        }
        count <- count + 1
      }
    }
    return(OIRF_matrix)
  }
  
  compute_GIRF <- function(A, Sigma, e, N){
    k <- nrow(Sigma)
    GIRF_matrix <- matrix(NA, k*k, N+1)
    count <- 1
    
    for (j in 1:k){
      for (i in 1:k){
        for (n in 0:N){
          GIRF <- Sigma[j,j]^(-1/2) * A[[n+1]] %*% Sigma %*% e[[j]]
          GIRF_matrix[count, n+1] <- GIRF[i]
        }
        count <- count + 1
      }
    }
    return(GIRF_matrix)
  }
  
  k <- x$setup$k
  p <- x$setup$p
  
  var_names <- colnames(x$data)
  
  if (is.null(var_names)) {
    var_names <- paste0("y_", 1:k)
  }
  
  
  if (estimation == "stan"){
    stan_draws <- rstan::extract(x$fit$stan)
    N_draws <- dim(stan_draws$beta)[1]
  } else {
    gibbs_draws <- x$fit$gibbs
    N_draws <- dim(gibbs_draws$beta_draws)[3]
  }
  
  
  irf_array <- array(NA, dim = c(N_draws, k, k, H+1))
  
  e <- vector("list", k)
  for (nE in 1:k) {
    tmp = rep(0, k)
    tmp[nE] = 1
    e[[nE]] <- tmp
  }
  
  
  for (draw in 1:N_draws) {
    if (estimation == "gibbs") {
      Phi <- t(gibbs_draws$beta_draws[,,draw])
      Sigma <- gibbs_draws$Sigma_u_draws[,,draw]
    } else {
      Phi <- t(stan_draws$beta[draw,,])
      if (!is.null(x$SV)) {
        if (is.null(t)) t = dim(stan_draws$Sigma_u)[2]
        Sigma <- stan_draws$Sigma_u[draw,t,,]
      } else {
        Sigma <- stan_draws$Sigma_u[draw,,]
      }
    }
    
    N <- H
    Psi <- MTS::VARpsi(Phi, lag=N)$psi
    kk <- ncol(Psi)/nrow(Psi)
    
    A <- vector("list", length=kk)
    for (nA in 1:kk) {
      col_idx <- ((nA - 1) * nrow(Psi) + 1):(nA * nrow(Psi))
      A[[nA]] <- Psi[, col_idx]
    }
    
    if (method == "OIRF") {
      tmp_IRF <- compute_OIRF(A,Sigma,e,N)
    } else {
      tmp_IRF <- compute_GIRF(A,Sigma,e,N)
    }
    
    for (j in 1:k){
      for (i in 1:k){
        row_idx <- (j-1)*k + i
        irf_array[draw, i, j, ] <- tmp_IRF[row_idx, ]
      }
    }
  }
  alpha <- 1 - ci
  if (type == "median") {
    m_irf <- apply(irf_array, c(2,3,4), median)
  } else {
    m_irf <- apply(irf_array, c(2,3,4), mean)
  }
  lower_irf  <- apply(irf_array, c(2,3,4), quantile, probs = alpha/2)
  upper_irf  <- apply(irf_array, c(2,3,4), quantile, probs = 1 - alpha/2)
  
  horizon <- 0:H
  max_abs <- max(abs(c(
    lower_irf[response, shock, ],
    upper_irf[response, shock, ],
    m_irf[response, shock, ]
  )))
  
  ylim_range <- c(-max_abs, max_abs)
  
  plot(NA, xlim = c(0, H), ylim = ylim_range, type="n",
       xlab="Horizon", ylab=paste0("Response: ", var_names[response]), font.lab=2)
  
  polygon(c(horizon, rev(horizon)),
          c(upper_irf[response, shock, ], rev(lower_irf[response, shock, ])),
          col=rgb(0,0,1,0.2), border=NA)
  
  lines(horizon, m_irf[response, shock, ], col="blue", lwd=2)
  abline(h=0, col="black", lty=1)
  
  if (type == "median") {
    type_label <- "Median"
  } else {
    type_label <- "Mean"
  }
  
  if (!is.null(x$SV)) {
    title(main = paste0(
      "Posterior ", type_label, " ", method, " (", round(ci*100), "% probability bands)\nt=", t, "\n",
      "Shock: ", var_names[shock]
    ))
    
  } else {
    title(main = paste0(
      "Posterior ", type_label, " ", method, " (", round(ci*100), "% probability bands)\n",
      "\nShock: ", var_names[shock]
    ))
  }

  
  if (type == "median") {
    return(list(median_irf = m_irf, lower = lower_irf, upper = upper_irf))
  } else {
    return(list(mean_irf = m_irf, lower = lower_irf, upper = upper_irf))
  }
}

