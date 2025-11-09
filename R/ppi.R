ppi <- function(l, u, annualize_growthrate = FALSE, alpha = 0.05) {

  z <- qnorm(1 - alpha/2)
  
  if (!annualize_growthrate) {
    sigma <- (u - l) / (2 * z)
    mu <- (u + l) / 2
  } else {
    mu <- (u + l) / 2 / 4
    sigma <- (u - l) / (2 * z) / 4
  }
  
  return(list(mean = mu, var = sigma^2))
}
