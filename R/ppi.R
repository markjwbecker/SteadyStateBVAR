ppi <- function(l, u, interval=0.95, annualized_growthrate = FALSE, freq=4) {

  alpha = 1-interval
  z <- qnorm(1 - alpha/2)
  
  if (!annualized_growthrate) {
    sigma <- (u - l) / (2 * z)
    mu <- (u + l) / 2
  } else {
    mu <- (u + l) / 2 / freq
    sigma <- (u - l) / (2 * z) / freq
  }
  
  return(list(mean = mu, var = sigma^2))
}
