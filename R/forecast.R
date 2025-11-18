forecast <- function(x, ci=NULL, fcst_type=c("mean", "median"), 
                          growth_rate_idx=NULL, plot_idx=NULL, estimation = c("stan", "gibbs"),show_all=FALSE,xlim=NULL,ylim=NULL)
  {
  xlim_user <- xlim
  ylim_user <- ylim
  Y <- x$data
  freq <- frequency(Y)
  if (is.null(plot_idx)) plot_idx <- 1:ncol(Y)
  
  if (estimation=="gibbs") {
    Y_pred <- x$fit$gibbs$fcst_draws
    Y_pred_m <- apply(Y_pred, c(1, 2), fcst_type)
    if (!is.null(ci)){
      alpha <- 1 - ci
      Y_pred_lower <- apply(Y_pred, c(1, 2), quantile, probs = alpha/2)
      Y_pred_upper <- apply(Y_pred, c(1, 2), quantile, probs = 1 - alpha/2)
    } else {
      Y_pred_sd <- apply(Y_pred, c(1, 2), sd)
      Y_pred_lower <- Y_pred_m - Y_pred_sd
      Y_pred_upper <- Y_pred_m + Y_pred_sd
    }

  } else {
    posterior <- rstan::extract(x$fit$stan)
    Y_pred <- posterior$Y_pred
    Y_pred_m <- apply(Y_pred, c(2, 3), fcst_type)
    if (!is.null(ci)){
    alpha <- 1 - ci
    Y_pred_lower <- apply(Y_pred, c(2, 3), quantile, probs = alpha/2)
    Y_pred_upper <- apply(Y_pred, c(2, 3), quantile, probs = 1 - alpha/2)
    } else {
      Y_pred_sd <- apply(Y_pred, c(2, 3), sd)
      Y_pred_lower <- Y_pred_m - Y_pred_sd
      Y_pred_upper <- Y_pred_m + Y_pred_sd
    }
  }
  
  
  T <- nrow(Y)
  H <- nrow(Y_pred_m)
  m <- ncol(Y)
  
  time_hist <- time(Y)
  time_fore <- seq(tail(time_hist, 1) + 1/freq, by = 1/freq, length.out = H)
  
  if (is.null(plot_idx)) plot_idx <- 1:ncol(Y)
  for (i in plot_idx) {
    smply <- Y[, i]
    fcst_m <- Y_pred_m[, i]
    fcst_lower <- Y_pred_lower[, i]
    fcst_upper <- Y_pred_upper[, i]
    
    if (!is.null(growth_rate_idx) && i %in% growth_rate_idx) {
      
      annual_hist <- rep(NA, length(smply))
      for(t in freq:length(smply)){
        annual_hist[t] <- sum(smply[(t-(freq-1)):t])
      }
      annual_hist <- annual_hist
      annual_hist <- ts(annual_hist, start = start(Y), frequency = freq)
      
      last_obs <- tail(smply,(freq-1))
      all_fcst <- c(last_obs, fcst_m)
      
      annual_fcst <- rep(NA, H)
      for(t_h in 1:H){
        annual_fcst[t_h] <- sum(all_fcst[t_h:(t_h+(freq-1))])
      }
      annual_fcst <- annual_fcst
      
      
      all_lower <- c(last_obs, fcst_lower)
      all_upper <- c(last_obs, fcst_upper)
      
      annual_lower <- rep(NA,H)
      annual_upper <- rep(NA,H)
      for(t_h in 1:H){
        annual_lower[t_h] <- sum(all_lower[t_h:(t_h+(freq-1))])
        annual_upper[t_h] <- sum(all_upper[t_h:(t_h+(freq-1))])
      }
      annual_lower <- annual_lower
      annual_upper <- annual_upper
      
      time_full <- c(tail(time_hist, 1), time_fore)
      
      m_full <- c(tail(annual_hist, 1), annual_fcst)
      lower_full <- c(tail(annual_hist, 1), annual_lower)
      upper_full <- c(tail(annual_hist, 1), annual_upper)
      
      ylim <- range(c(annual_lower, annual_upper),na.rm=TRUE)
      ymin <- floor(ylim[1]*2)/2
      ymax <- ceiling(ylim[2]*2)/2
      yticks <- seq(ymin-100, ymax+100, by = 0.5)
      
      if (isFALSE(show_all) && is.null(xlim) && is.null(ylim)){
        plot.ts(annual_hist, main = paste(colnames(Y)[i], "(annual)"), xlab = "Time", ylab = NULL,
                xlim = c(head(time_fore,1)-8,tail(time_fore,1)),
                ylim = ylim,col = "black", lwd = 2, yaxt = "n")
        xlim_vals <- c(head(time_hist, 1)-100, tail(time_fore, 1)+100)
        x_quarters <- seq(xlim_vals[1], xlim_vals[2], by = 1/freq)
        #abline(v = x_quarters, col = "gray", lty = 2) 
        abline(h = seq(ymin-100, ymax+100, by = 0.5), col = "gray", lty = 2)
        axis(side = 2, at = yticks, labels = yticks, las = 1)
        points(as.numeric(time_hist), annual_hist, pch=16, col="black")
        points(time_full[-1], m_full[-1], pch=16, col="blue")
      } else if (!is.null(xlim) && !is.null(ylim)) {
        plot.ts(annual_hist, main = paste(colnames(Y)[i], "(annual)"), xlab = "Time", ylab = NULL,
                col = "black", lwd = 2,xlim=xlim_user,
                ylim= ylim_user, yaxt="n")
        xlim_vals <- c(head(time_hist, 1)-100, tail(time_fore, 1)+100)
        x_quarters <- seq(xlim_vals[1], xlim_vals[2], by = 1/freq)
        #abline(v = x_quarters, col = "gray", lty = 2) 
        abline(h = seq(ymin-100, ymax+100, by = 0.5), col = "gray", lty = 2)
        axis(side = 2, at = yticks, labels = yticks, las = 1)
        points(as.numeric(time_hist), annual_hist, pch=16, col="black")
        points(time_full[-1], m_full[-1], pch=16, col="blue")
      } else {
        plot.ts(annual_hist, main = paste(colnames(Y)[i], "(annual)"), xlab = "Time", ylab = NULL,
                col = "black", lwd = 2,xlim=c(head(time_hist, 1), tail(time_fore, 1)),
                ylim=range(upper_full,lower_full,c(annual_hist), na.rm=TRUE))
      }
      
      polygon(
        c(time_full, rev(time_full)),
        c(upper_full, rev(lower_full)),
        col = rgb(0, 0, 1, 0.2),
        border = NA
      )
      lines(time_full, m_full, col = "blue", lwd = 2)
      
    } else {
      
      time_full <- c(tail(time_hist, 1), time_fore)
      m_full <- c(tail(smply, 1), fcst_m)
      lower_full <- c(tail(smply, 1), fcst_lower)
      upper_full <- c(tail(smply, 1), fcst_upper)
      
      ylim <- range(c(lower_full, upper_full))
      ymin <- floor(ylim[1]*2)/2
      ymax <- ceiling(ylim[2]*2)/2
      yticks <- seq(ymin-100, ymax+100, by = 0.5)
      
      if (isFALSE(show_all) && is.null(xlim) && is.null(ylim)){
        plot.ts(smply, main = colnames(Y)[i], xlab = "Time", ylab = NULL,
                xlim = c(head(time_fore,1)-8,tail(time_fore,1)),
                ylim = ylim,col = "black", lwd = 2, yaxt = "n")
        xlim_vals <- c(head(time_hist, 1)-100, tail(time_fore, 1)+100)
        x_quarters <- seq(xlim_vals[1], xlim_vals[2], by = 1/freq)
        abline(v = x_quarters, col = "gray", lty = 2) 
        abline(h = seq(ymin-100, ymax+100, by = 0.5), col = "gray", lty = 2)
        axis(side = 2, at = yticks, labels = yticks, las = 1)
        points(as.numeric(time_hist), smply, pch=16, col="black")
        points(time_full[-1], m_full[-1], pch=16, col="blue")
    } else if (!is.null(xlim) && !is.null(ylim)) {
      plot.ts(smply, main = paste(colnames(Y)[i], "(annual)"), xlab = "Time", ylab = NULL,
              col = "black", lwd = 2,xlim=xlim_user,
              ylim= ylim_user, yaxt="n")
      xlim_vals <- c(head(time_hist, 1)-100, tail(time_fore, 1)+100)
      x_quarters <- seq(xlim_vals[1], xlim_vals[2], by = 1/freq)
      #abline(v = x_quarters, col = "gray", lty = 2) 
      abline(h = seq(ymin-100, ymax+100, by = 0.5), col = "gray", lty = 2)
      axis(side = 2, at = yticks, labels = yticks, las = 1)
      points(as.numeric(time_hist), smply, pch=16, col="black")
      points(time_full[-1], m_full[-1], pch=16, col="blue")
      
     } else {
        plot.ts(smply, main = colnames(Y)[i], xlab = "Time", ylab = NULL,
                col = "black", lwd = 2,xlim=c(head(time_hist, 1), tail(time_fore, 1)),
                ylim=range(upper_full,lower_full,smply))
      }
      
      polygon(
        c(time_full, rev(time_full)),
        c(upper_full, rev(lower_full)),
        col = rgb(0, 0, 1, 0.2),
        border = NA
      )
      lines(time_full, m_full, col = "blue", lwd = 2)
    }
  }
  x$forecasts$forecast = Y_pred_m
  x$forecasts$lower = Y_pred_lower
  x$forecasts$upper = Y_pred_upper
  
  return(x$forecasts)
}
