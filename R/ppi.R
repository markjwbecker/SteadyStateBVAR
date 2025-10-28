ppi <- function(l, u){
  sigma <- (u - l) / (2 * qnorm(0.975))
  mu <- (u + l) / 2
  return(list(mean = mu, var = sigma^2))
}