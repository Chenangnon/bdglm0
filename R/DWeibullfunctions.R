# Probability mass function
# x = x values (quantiles) which the probability masses are desired at
# mu = mean value
# c = shape parameter
# log logical
#' @name ddweibull
#' @aliases pdweibull qdweibull rdweibull
#' @title The Balanced Discrete Weibull Distribution
#' @description Probability mass function, cumulative distribution function, quantile function, random variate generation,
#' and variance function, moment estimation, negative log-likelihood, maximum likelihood estimation for the balanced discrete log normal distribution with parameters \code{mu} (mean) and \code{alpha} (standard deviation at natural logarithmic scale).
#'
#' @usage
#' ddweibull (x, mu, c = 1, log = FALSE)
#' pdweibull (q, mu, c = 1, lower.tail = TRUE, log.p = FALSE)
#' qdweibull (p, mu, c = 1, lower.tail = TRUE, log.p = FALSE)
#' rdweibull (n, mu, c = 1)
#'
#' @param x,q vector of integer/real quantiles.
#' @param p vector of real probability values in the open \code{(0, 1)}.
#' @param n scalar, number of deviates to generate.
#' @param mu vector of positive real mean parameters of the balanced discrete distribution.
#' @param c vector of positive real scale parameters of the balanced discrete distribution.
#' @param log,log.p logical; if TRUE, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#'
#' @details The functions are the direct implementations of formulae in Tovissode et al. (2020).
#'
#' @return \code{ddweibull} gives the probability mass, \code{pdweibull} the distribution function, \code{qdweibull} the quantile function, and \code{rdweibull}
#' generates random deviates. Invalid arguments will result in errors.
#'
#' The length of the result is determined by \code{n} for \code{rdweibull}, and is the maximum of the lengths
#' of the numerical arguments for the other functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result. Only the first
#' elements of the logical arguments are used.
#'
#' @references Tovissode, C.F.; Honfo, S.H.; Doumate, J.T.; Glele Kakai, R. On the Discretization of Continuous Probability Distributions Using a Probabilistic Rounding Mechanism. Mathematics 2021, 9, 555. https://doi.org/10.3390/math9050555
#'
#'  @seealso \link{bdweibull} for generating balanced discrete Weibull family objects for model fitting purposes.

ddweibull <- function(x, mu, c = 1, log = FALSE) {
  ddggamma (x = x, mu = mu, b = 1, c = c, log = log)
}

# Cumulative distribution function
pdweibull <- function(q, mu, c = 1, lower.tail = TRUE, log.p = FALSE) {
  pdggamma (q = q, mu = mu, b = 1, c = c, lower.tail = lower.tail, log.p = log.p)
}

# Quantile function
qdweibull <- function(p, mu, c = 1, lower.tail = TRUE, log.p = FALSE) {
  qdggamma (p = p, mu = mu, b = 1, c = c, lower.tail = lower.tail, log.p = log.p)
}

# Generate random deviates
rdweibull <- function(n, mu, c = 1) {
  rdggamma(n, mu = mu, b = 1, c = c)
}

# Routine for .wzeta function
.wzeta_comps <- function(z, c = 1, a = 1) {
  .ggzeta_comps(z, b = 1, c = c, a = a)
}

# Routine for vardweibull function
.wzeta <- function(mu, c = 1, tol = 1e-10) {
  .ggzeta(mu = mu, b = 1, c = c, tol = tol)
}

# var function
.vardweibull <- function(mu, c = 1, tol = 1e-10) {
  na <- length(a)
  if (length(b) == 1 && na == 1)
    return(as.vector((gamma(1+2/c) / (gamma(1+1/c)^2) - 1) * mu^2 +
                       .wzeta(mu = mu, c = c, tol = tol)))
  AB <- cbind(mu, c)
  mu <- AB[,1]
  c <- AB[,2]
  Wvalues <- apply(AB, MARGIN = 1, FUN = function(ab) {
    .wzeta(mu = ab[1], c = ab[2], tol = tol)
  })
  as.vector((gamma(1+2/c) / (gamma(1+1/c)^2) - 1) * mu^2 + Wvalues)
}

# Find c such that the balanced discrete Weibull distribution with mean mu is equidispersed
.equic_bdweibull <- function (mu, tol = 1e-10) {
  tol <- min(1e-5, max(1e-15, tol))
  lfun <- function(mui) {
    if (mui == 0)
      return(1)
    fun <- function(x) vardweibull(mu = mui, c = exp(x), tol = tol)  - mui
    res <- uniroot(f = fun,
                   lower = logb(.5), upper = logb(5), extendInt = "yes", tol = tol)
    structure(as.vector(exp(res$root)), diffEV = as.vector(res$f.root))
  }
  if (length(mu) == 1)
    lfun(mu)
  else
    sapply(mu, FUN = lfun)
}

# moment estimates
.momente_bdweibull <- function(x, zeta = NULL, tol = 1e-5, maxit = 10) {
  xbar <- mean(x)
  sigma2 <- var(x)
  if (is.null(zeta)) {
    minz <- min(xbar, .25)
    zeta <- if (sigma2 > minz) minz/2 else sigma2/4
  }
  Const <- (sigma2 - zeta)/xbar^2 + 1
  fun = function(lc)
    gamma(1+2*exp(lc))/gamma(1+exp(lc))^2 - Const
  res <- uniroot(f = fun,
                 lower = logb(.75), upper = logb(4), extendInt = "yes",
                 tol = tol, maxiter = maxit)
  c(mu = xbar, c = 1/exp(res$root), sigma2 = sigma2)
}

# Fit a discrete Weibull distribution to an IDD sample x
.mle_bdweibull <- function(x, theta = NULL, alpha = .05, asymp = "norm",
                         method = "BFGS", tol = 1e-10, maxit = 100) {
  if (missing(theta) || is.null(theta))
    start <- logb(.momente_bdweibull(x, maxit = 10)[1:2])
  else {
    start <- logb(theta[1:2])
    names(start) <- c("mu", "c")
  }
  # Negative log likelihood function: theta = c(log(mu), log(alpha))
  nlike.dweibull <- function(theta, x) {
    -sum(ddweibull(x, mu = exp(theta[1]), c = exp(theta[2]), log = TRUE))
  }
  OPTIM <- optim(par = start, fn = nlike.dweibull, x = x,
                 method = method, hessian = TRUE,
                 control = list(maxit = maxit))
  theta <- exp(OPTIM$par)
  log.theta.se <- try(solve(OPTIM$hessian), silent = TRUE)
  if (class(log.theta.se) == "try-error") {
    return(list(theta = theta, CI = NULL, se = NULL, df = length(x)-2,
                  log.se = log.theta.se, cor = NULL, start = exp(start),
                  logLike = -OPTIM$value, AIC = 2 * (OPTIM$value + 2)))
  }
  theta.corr <- cov2cor(log.theta.se)[2]
  names(theta.corr) <- "mu:c"
  log.theta.se <- sqrt(diag(log.theta.se))
  theta.se <- sqrt(theta) * log.theta.se
  df <- length(x) - 2
  if (identical(asymp, "norm"))
    qalpha <- qnorm(alpha/2, lower.tail = FALSE)
  else
    qalpha <- qt(alpha/2, df = df, lower.tail = FALSE)
  half.CI <- qalpha * log.theta.se
  theta.CI <- cbind(lower = exp(OPTIM$par - half.CI), upper = exp(OPTIM$par + half.CI))
  # Add a Test for equidispersion: mu == sigma2
  list(theta = theta, CI = theta.CI, se = theta.se, df = df,
       log.se = log.theta.se, cor = theta.corr, start = exp(start),
       logLike = -OPTIM$value, AIC = 2 * (OPTIM$value + 2))
}
