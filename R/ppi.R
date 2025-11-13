ppi <- function(l, u, alpha = 0.05, annualized_growthrate = FALSE, freq=4) {

  z <- qnorm(1 - alpha/2)
  
  if (!annualize_growthrate) {
    sigma <- (u - l) / (2 * z)
    mu <- (u + l) / 2
  } else {
    mu <- (u + l) / 2 / freq
    sigma <- (u - l) / (2 * z) / freq
  }
  
  return(list(mean = mu, var = sigma^2))
}
