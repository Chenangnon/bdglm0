#
# Add the survival function S(t)=1-F(t),
# the hazard rate function h(t) = f(t)/S(t)
# the cumulative hazard function sum_{j=0}^t h(j)
#
#' @name ddgamma
#' @aliases pdgamma qdgamma rdgamma
#' @title The Balanced Discrete Gamma Distribution
#' @description Probability mass function, cumulative distribution function, quantile function,
#' and random variate generation, for the balanced discrete gamma distribution with parameters
#' mu (mean) and a (scale).
#'
#' @usage
#' ddgamma (x, mu, a = b/mu, b = mu, log = FALSE)
#'
#' pdgamma (q, mu, a = b/mu, b = mu, lower.tail = TRUE, log.p = FALSE)
#'
#' qdgamma (p, mu, a = b/mu, b = mu, lower.tail = TRUE, log.p = FALSE)
#'
#' rdgamma (n, mu, a = b/mu, b = mu)
#'
#' @param x vector of integer quantiles.
#' @param q vector of real quantiles.
#' @param p vector of real probability values in the open (0, 1).
#' @param mu vector of positive real mean parameters of the balanced discrete distribution.
#' @param a vector of positive real scale parameters of the balanced discrete distribution.
#' @param b vector of positive reals, alternative way to specify the scale of the balanced discrete distribution.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @param n scalar, number of variates to generate.
#'
#' @details The functions are the direct implementations of formulae in Tovissode et al. (2021).
#'
#' @return \code{ddgamma} gives the probability mass, \code{pdgamma} the distribution function,
#' \code{qdgamma} the quantile function, and \code{rdgamma} generates random deviates.
#' Invalid arguments will result in errors.
#'
#' The length of the result is determined by \code{n} for \code{rdgamma}, and is the maximum of the lengths
#' of the numerical arguments for the other functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result. Only the first
#' elements of the logical arguments are used.
#'
#' @seealso \link{bdvariance} for variance computation and the estimation of distribution parameters;
#' and \link{bdgamma} for generating balanced discrete gamma family objects for model fitting purposes.
#'
#' @references Tovissode, C.F.; Honfo, S.H.; Doumate, J.T.; Glele Kakai, R. On the Discretization of Continuous Probability Distributions Using a Probabilistic Rounding Mechanism. Mathematics 2021, 9, 555. https://doi.org/10.3390/math9050555

ddgamma <- function(x, mu, a = b/mu, b = mu, log = FALSE) {
  if (!missing(a))
    b <- a * mu
  zx <- x == 0
  mx_ <- a * (x - 1)
  mx_[zx] <- 0
  d <- tryCatchWE((x - 1) * pgamma(mx_, shape = b) -
                      2 * x * pgamma(a * x, shape = b) +
                      (x + 1) * pgamma(a * (x + 1), shape = b) -
                      mu * (pgamma(a * (x + 1), shape = b + 1) -
                              2 * pgamma(a * x, shape = b + 1) +
                              pgamma(mx_, shape = b + 1)))$value
  d <- pmax(d, exp(-708))
  if (any(zx0 <- ((x < 0) + (x != floor(x))) > 0))
    d[zx0] <- 0
  if (log)
    logb(d)
  else
    d
}

# Cumulative distribution function
pdgamma <- function(q, mu, a = b/mu, b = mu, lower.tail = TRUE, log.p = FALSE) {
  if (!missing(a))
    b <- a * mu
  q <- floor(q)
  pg <- pgamma(a * q, shape = b)
  p <- pg + (q + 1) * (pgamma(a * (q + 1), shape = b) - pg) -
    mu * (pgamma(a * (q + 1), shape = b + 1) - pgamma(a * q, shape = b + 1))
  if (any(zq0 <- q < 0))
    p[zq0] <- 0
  if (!lower.tail)
    p <- 1 - p
  if(log.p)
    logb(p)
  p
}

# Quantile function
qdgamma <- function(p, mu, a = b/mu, b = mu, lower.tail = TRUE, log.p = FALSE) {
  if (!missing(a))
    b <- a * mu
  if(log.p)
    p <- exp(p)
  if (!lower.tail)
    p <- 1 - p
  x <- floor(qgamma(p, shape = b, lower.tail = TRUE, log.p = FALSE) / a)
  u <- pdgamma(x, mu = mu, a = a)
  pu <- u < p
  x[pu] <- x[pu] + 1
  x
}

# Generate random deviates
rdgamma <- function(n, mu, a = b/mu, b = mu) {
  if (!missing(a))
    b <- a * mu
  rg <- rgamma(n, shape = b, rate = a)
  z <- floor(rg)
  r <- rg - z
  nzr = r > 0
  if (any(nzr))
    z[nzr] <- z[nzr] + rbinom(sum(nzr), size = 1, prob = r[nzr])
  z
}

.vardgamma <- function(mu, a = b/mu, b = mu, tol = 1e-10) {
  if (!missing(a))
    b <- a * mu
  na <- length(a)
  if (length(b) == 1 && na == 1)
    return(as.vector(mu / a + .gzeta(b = b, a = a, tol = tol)))
  AB <- cbind(a, b)
  a <- AB[,2]
  Gvalues <- apply(AB, MARGIN = 1, FUN = function(ab) {
    .gzeta(b = ab[1], a = ab[2], tol = tol)
  })
  return(as.vector(mu / a + Gvalues))
}

# Routine for .gzeta function
.gzeta_comps <- function(z, b, a = 1) {
  comp1 <- sum((2 * z + 1) * (pgamma(a * (z + 1), shape = b + 1) -
                                pgamma(a * z, shape = b + 1)))
  comp2 <- sum(pgamma(a * (z + 1), shape = b + 2) -
                 pgamma(a * z, shape = b + 2))
  comp3 <- sum(z * (z + 1) * (pgamma(a * (z + 1), shape = b) -
                                pgamma(a * z, shape = b)))
  c(comp1, comp2, comp3)
}

# Routine for vardgamma function
.gzeta <- function(b, a = 1, tol = 1e-10) {
  z_i <- floor(qgamma(tol/2, shape = b) / a)
  z_f <- floor(qgamma(1 - tol/2, shape = b) / a) + 1
  zcomps <- .gzeta_comps(z_i:z_f, b = b, a = a)
  res <- zcomps[1] * b / a - b * (b + 1) * zcomps[2] / a^2 - zcomps[3]
  structure(res, precision = pgamma(a * z_i, shape = b) + pgamma(a * (z_f + 1), shape = b, lower.tail = FALSE))
}

# approximate moment estimates (approximate for a)
.momente_bdgamma <- function(x, zeta = NULL) {
  xbar <- mean(x)
  sigma2 <- var(x)
  if (is.null(zeta)) {
    minz <- min(xbar/2, .125)
    zeta <- if (sigma2 > minz) minz/2 else sigma2/4
  }
  a <- xbar / (sigma2  - zeta)
  c(mu = xbar, a = a)
}

# An internal version of momente.dgamma
.momente.dgamma <- function(x, zeta = if (sigma2 > minz) minz/2 else sigma2/4) {
  xbar <- mean(x)
  sigma2 <- var(x)
  minz <- min(xbar, .25)
  a <- xbar / (sigma2  - zeta)
  c(mu = xbar, a = a, sigma2 = sigma2)
}

# a for equi-dispersion
.equia_bdgamma <- function (mu, tol = 1e-10) {
  tol <- min(1e-5, max(1e-15, tol))
  lfun <- function(mui) {
    if (mui == 0)
      return(1)
    fun <- function(x) mui / exp(x) + .gzeta(b = exp(x) * mui, a = exp(x), tol = tol) - mui
    res <- uniroot(f = fun,
                   lower = 0, upper = logb(5), extendInt = "yes", tol = tol)
    structure(exp(as.vector(res$root)), diffEV = as.vector(res$f.root))
  }
  if (length(mu) == 1)
    lfun(mu)
  else
    sapply(mu, FUN = lfun)
}

# Fit a discrete gamma distribution to an IDD sample x
.mle_bdgamma <- function(x, theta = NULL, alpha = .05, asymp = "norm",
                       method = "BFGS", tol = 1e-10, maxit = 100) {
  if (missing(theta) || is.null(theta))
    start <- logb(.momente_bdgamma(x)[1:2])
  else {
    start <- logb(theta[1:2])
    names(start) <- c("mu", "a")
  }
  # Negative log likelihood function: theta = c(log(mu), log(a))
  nlike.dgamma <- function(theta, x) {
    -sum(ddgamma(x, mu = exp(theta[1]), a = exp(theta[2]), log = TRUE))
  }
  OPTIM <- optim(par = start, fn = nlike.dgamma, x = x,
                 method = method, hessian = TRUE,
                 control = list(maxit = maxit))
  theta <- exp(OPTIM$par)
  log.theta.se <- try(solve(OPTIM$hessian), silent = TRUE)
  if (any(class(log.theta.se) %in% c("try-error", "error", "condition"))) {
    return(list(theta = theta, CI = NULL, se = NULL, df = length(x)-2,
                log.se = log.theta.se, cor = NULL, start = exp(start),
                logLike = -OPTIM$value, AIC = 2 * (OPTIM$value + 2)))
  }
  theta.corr <- cov2cor(log.theta.se)[2]
  names(theta.corr) <- "mu:a"
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
