ppi <- function(l, u, growthrate = FALSE) {

  if (!growthrate) {
    sigma <- (u - l) / (2 * qnorm(0.975))
    mu <- (u + l) / 2
  } else {
    mu <- (u + l) / 2 / 4
    sigma <- (u - l) / (2 * qnorm(0.975)) / 4
  }
  
  return(list(mean = mu, var = sigma^2))
}
