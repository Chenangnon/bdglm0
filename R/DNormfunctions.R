#' @name ddnorm
#' @title The Balanced Discrete Normal Distribution
#'
#' @description Probability mass function, cumulative distribution function, quantile function,
#' and random variate generation, for the balanced discrete normal distribution with parameters
#' \code{mu} (mean) and \code{sigma} (scale).

ddnorm <-  function(x, mu = 0, sigma = 1, log = FALSE) {
  d <- (x-1) * pnorm(x-1, mean = mu, sd = sigma) -
    2 * x * pnorm(x, mean = mu, sd = sigma) +
    (x+1) * pnorm(x+1, mean = mu, sd = sigma) -
  .pmomemtnorm(x, r = 1, mean = mu, sd = sigma) +
    .pmomemtnorm(x-1, r = 1, mean = mu, sd = sigma)
  if(log)
    logb(d)
  else
    d
}

.pmomemtnorm <- function (ymin, ymax = ymin + 1, r = 1, mean = 0, sd = 1) {
  if (r == 0) {
    return(pnorm(ymax, mean = mean, sd = sd) - pnorm(ymin, mean = mean, sd = sd))
  }
  zmin <- (ymin - mean)/sd
  zmax <- (ymax - mean)/sd
  mean * (pnorm(zmax) - pnorm(zmin)) - sd * (dnorm(zmax) - dnorm(zmin))
}

