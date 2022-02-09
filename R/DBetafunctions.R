#' @name ddbeta
#' @aliases pdbeta qdbeta rdbeta
#' @title The Balanced Discrete Beta Distribution
#' @description Probability mass function, cumulative distribution function, quantile function, and random variate generation for the balanced discrete beta distribution with parameters \code{mu} (mean), \code{a} (scale) and \code{shape}.
#'
#' @usage
#' ddbeta (x, mu, a = 1, shape = mu, log = FALSE)
#' pdbeta (q, mu, a = 1, shape = mu, lower.tail = TRUE, log.p = FALSE)
#' qdbeta (p, mu, a = 1, shape = mu, lower.tail = TRUE, log.p = FALSE)
#' rdbeta (n, mu, a = 1, shape = mu)
#'
#' @param x,q vector of integer/real quantiles.
#' @param p vector of real probability values in the open \code{(0, 1)}.
#' @param mu vector of positive real mean parameters of the balanced discrete distribution.
#' @param a vector of positive real scale parameters of the balanced discrete distribution.
#' @param shape vector of positive reals shape parameters.
#' @param log,log.p logical; if TRUE, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @param n scalar, number of random variates to generate.
#'
#' @details The functions are the direct implementations of formulae in Tovissode et al. (2021).
#'
#' @return \code{ddbeta} gives the probability mass, \code{pdbeta} the distribution function, \code{qdbeta} the quantile function, and \code{rdbeta}
#' generates random deviates. Invalid arguments will result in errors.
#'
#' The length of the result is determined by \code{n} for \code{rdbeta}, and is the maximum of the lengths
#' of the numerical arguments for the other functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result. Only the first
#' elements of the logical arguments are used.
#'
#' @references Tovissode, C.F.; Honfo, S.H.; Doumate, J.T.; Glele Kakai, R. On the Discretization of Continuous Probability Distributions Using a Probabilistic Rounding Mechanism. Mathematics 2021, 9, 555. https://doi.org/10.3390/math9050555
#'
#' @seealso \link{bdbeta} for generating balanced discrete beta family objects for model fitting purposes.
#'
#
ddbeta <- function(x, mu, a = 1, shape = mu, log = FALSE) {
  m <- 1 / (mu + shape)
  b <- a * shape / mu
  zx <- x == 0
  mx_ <- m * (x - 1)
  mx_[zx] <- 0
  mx <- m * x
  mxp <- m * (x + 1)
  d <- (x - 1) * pbeta(mx_, shape1 = a, shape2 = b) -
    2 * x * pbeta(mx, shape1 = a, shape2 = b) +
    (x + 1) * pbeta(mxp, shape1 = a, shape2 = b) -
    mu * (pbeta(mxp, shape1 = a + 1, shape2 = b) -
            2 * pbeta(mx, shape1 = a + 1, shape2 = b) +
            pbeta(mx_, shape1 = a + 1, shape2 = b))
  if (any(zx0 <- ((x < 0) + (x != floor(x)) + (x > ceiling(mu + shape))) > 0))
    d[zx0] <- 0
  if (any(zx0 <- ((mu == 0) + (x != 0)) == 2))
    d[zx0] <- 0
  if (log)
    logb(d)
  else
    d
}

# Cumulative distribution function
pdbeta <- function(q, mu, a = 1, shape = mu, lower.tail = TRUE, log.p = FALSE) {
  q <- floor(q)
  m <- 1 / (mu + shape)
  b <- a * shape / mu
  mq <- m * q
  mq1 <- m * (q + 1)
  pg <- pbeta(mq, shape1 = a, shape2 = b)
  p <- pg + (q + 1) * (pbeta(mq1, shape1 = a, shape2 = b) - pg) -
    mu * (pbeta(mq1, shape1 = a + 1, shape2 = b) -
            pbeta(mq, shape1 = a + 1, shape2 = b))
  if (any(zq0 <- q < 0))
    p[zq0] <- 0
  if (!lower.tail)
    p <- 1 - p
  if(log.p)
    logb(p)
  p
}

# Quantile function
qdbeta <- function(p, mu, a = 1, shape = mu, lower.tail = TRUE, log.p = FALSE) {
  m <- 1 / (mu + shape)
  b <- a * shape / mu
  if (log.p)
    p <- exp(p)
  if (!lower.tail)
    p <- 1 - p
  x <- floor(qbeta(p, shape1 = a, shape2 = b, lower.tail = TRUE, log.p = FALSE) / m)
  validp <- is.na(x)
  if (any(!validp))
    x[!validp] <- -1
  u <- pdbeta(x, mu = mu, shape = shape, a = a)
  pu <- u < p
  x[pu] <- x[pu] + 1
  if (any(!validp))
    x[!validp] <- NaN
  x
}

# Generate random deviates
rdbeta <- function(n, mu, a = 1, shape = mu) {
  rg <- rbeta(n, shape1 = a, shape2 = a * shape / mu) * (mu + shape)
  z <- floor(rg)
  r <- rg - z
  nzr <- r > 0
  if (any(nzr))
    z[nzr] <- z[nzr] + rbinom(sum(nzr), size = 1, prob = r[nzr])
  z
}

# Routine for .bzeta function
.bzeta_comps <- function(z, a, b = a, m = 0.5) {
  mz <- m * z
  mz1 <- m * (z + 1)
  comp1 <- sum((2 * z + 1) * (pbeta(mz1, shape1 = a + 1, shape2 = b) -
                                pbeta(mz, shape1 = a + 1, shape2 = b)))
  comp2 <- sum(pbeta(mz1, shape1 = a + 2, shape2 = b) -
                 pbeta(mz, shape1 = a + 2, shape2 = b))
  comp3 <- sum(z * (z + 1) * (pbeta(mz1, shape1 = a, shape2 = b) -
                                pbeta(mz, shape1 = a, shape2 = b)))
  c(comp1, comp2, comp3)
}

# Routine for .vardbeta function
.bzeta <- function(mu, a = 1, shape = mu, tol = 1e-15) {
  m <- 1 / (mu + shape)
  b <- a * shape / mu
  z_i <- floor(qbeta(tol/2, shape1 = a, shape2 = b) / m)
  z_f <- min(ceiling(1/m), floor(qbeta(1 - tol/2, shape1 = a, shape2 = b) / m) + 1)
  zcomps <- .bzeta_comps(z_i:z_f, a = a, b = b, m = m)
  res <- zcomps[1] * mu - (mu^2) * (a + 1) * zcomps[2] / (a + m * mu) - zcomps[3]
  structure(res, precision = pbeta(m * z_i, shape1 = a, shape2 = b) +
              pbeta(m * (z_f + 1), shape1 = a, shape2 = b, lower.tail = FALSE))
}

# variance function
.vardbeta <- function (mu, a = 1, shape = mu, tol = 1e-15) {
  if (missing(shape))
    shape <- mu
  na <- max(length(a), length(mu))
  if (length(shape) == 1 && na == 1)
    return(as.vector(shape * (mu^2) / ((a+1) * mu + a * shape) +
                         .bzeta(mu = mu, shape = shape, a = a, tol = tol)))
  ABC <- cbind(mu, a, shape)
  mu <- ABC[,1]
  a <- ABC[,2]
  shape <- ABC[,3]
  Bvalues <- apply(ABC, MARGIN = 1, FUN = function(abc) {
    .bzeta(mu = abc[1], a = abc[2], shape = abc[3], tol = tol)
  })
  as.vector(shape * (mu^2) / ((a+1) * mu + a * shape) + Bvalues)
}

# moment estimates (approximate for shape and a)
.momente_bdbeta <- function(x, zeta = NULL, eps = 1e-15) {
  xbar <- mean(x)
  sigma2 <- var(x)
  if (is.null(zeta)) {
    minz <- min(xbar, .25)
    zeta <- if (sigma2 > minz) minz/2 else sigma2/4
  }
  shape <- max(max(x) - xbar, eps + (sigma2  - zeta) / xbar)
  a <- (xbar / (xbar + shape)) * (shape * xbar / (sigma2  - zeta) - 1)
  c(mu = xbar, shape = shape, a = a, sigma2 = sigma2)
}

# Given shape > 1 - .25 * mu, find a such that the balanced
# discrete beta distribution with mean mu is equidispersed
# Not SURE:: The defaults shape = mu + tol ensures that shape > 1 - .25 * mu even for mu <= 4/5
# Take in the paper shape = mu/2, mu, 2 * mu to study equidispersion behavior
.equia_bdbeta <- function (mu, shape = mu + tol, tol = 1e-15) {
  tol <- min(1e-5, max(1e-15, tol))
  lfun <- function(mui, shape) {
    if (any(c(mui, shape) < 0) || shape <= 1 - .25 * mui)
      return(NaN)
    if (mui == 0)
      return(1)
    upper <- if (mui > .25) mui * (shape * mui / (mui - .25) - 1) / (mui + shape)
    else if (shape > .5) mui * (2 * shape - 1) / (mui + shape)
    else return(NaN)
    lower <- min(upper-tol, max(tol, mui * (shape - 1) / (mui + shape)))
    fun <- function(x)
      shape * (mui^2) / ((x+1) * mui + x * shape) + .bzeta(mu = mui, shape = shape, a = x, tol = tol) - mui
    res <- uniroot(f = fun, lower = lower, upper = upper,
                   extendInt = "yes", tol = tol)
    structure(as.vector(res$root), diffEV = as.vector(res$f.root))
  }
  if ((nmu <- length(mu)) == 1)
    lfun(mu, shape = shape[1])
  else if ((nshape <- length(shape)) == 1)
    sapply(mu, FUN = function (x) lfun(x, shape = shape))
  else {
    shape <- rep(shape, length.out = nmu)
    sapply(1:nmu, FUN = function (i) lfun(mu[i], shape = shape[i]))
  }
}

# Given a > , find shape such that the balanced
# discrete beta distribution with mean mu is equidispersed
.equishape_bdbeta <- function (mu, a, tol = 1e-15) {
  tol <- min(1e-5, max(1e-15, tol))
  lfun <- function(mui, a) {
    if (mui == 0)
      return(1)
    if (a > mui)
      return(NaN) # To be studied, often but not always TRUE
    emui <- mui > .25
    maxa <- if (!emui) (mui ^ 2) / (mui - .25)
    else mui / .0000001
    if (any(c(mui, a) < 0) || a >= maxa)
      return(NaN)
    upper <- if (emui) mui / ((a+1) * (mui - .25)) - a / ((a + 1) * mui)
    else mui / ((a+1) * .000001) - a / ((a + 1) * mui)
    lower <- max(tol, min(upper-tol, (mui - a) / ((a+1) * mui)))
    upper <- max(lower+tol, upper, 1)
    fun <- function(x) {
      if (x <= 0)
        return(NaN)
      x * (mui^2) / ((a+1) * mui + a * x) + .bzeta(mu = mui, shape = x, a = a, tol = tol) - mui
    }
    res <- uniroot(f = fun, lower = lower, upper = upper,
                   extendInt = "yes", tol = tol)
    structure(as.vector(res$root), diffEV = as.vector(res$f.root))
  }
  if ((nmu <- length(mu)) == 1)
    lfun(mu, a = a[1])
  else if ((nshape <- length(a)) == 1)
    sapply(mu, FUN = function (x) lfun(x, a = a))
  else {
    a <- rep(a, length.out = nmu)
    sapply(1:nmu, FUN = function (i) lfun(mu[i], a = a[i]))
  }
}

# Fit a discrete beta distribution to an IDD sample x
.mle_bdbeta <- function(x, theta = NULL, alpha = .05, asymp = "norm",
                     method = "BFGS", maxit = 100) {
  if (missing(theta) || is.null(theta))
    start <- logb(.m.estdbeta(x)[1:3])
  else {
    start <- logb(theta[1:3])
    names(start) <- c("mu", "shape", "a")
  }
  # Negative log likelihood function: theta = c(log(mu), log(alpha))
  nlike.dbeta <- function(theta, x) {
    -sum(ddbeta(x, mu = exp(theta[1]), shape = exp(theta[2]), a = exp(theta[3]), log = TRUE))
  }
  OPTIM <- optim(par = start, fn = nlike.dbeta, x = x,
                 method = method, hessian = TRUE,
                 control = list(maxit = maxit))
  theta <- exp(OPTIM$par)
  log.theta.se <- try(solve(OPTIM$hessian), silent = TRUE)
  if (class(log.theta.se) == "try-error") {
    return(  list(theta = theta, CI = NULL, se = NULL, df = length(x)-2,
                  log.se = log.theta.se, cor = NULL, start = exp(start),
                  logLike = -OPTIM$value, AIC = 2 * (OPTIM$value + 2)))
  }
  theta.corr <- cov2cor(log.theta.se)[c(2, 3, 6)]
  names(theta.corr) <- c("mu:shape", "mu:a", "shape:a")
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
       log.se = log.theta.se, cor = theta.corr, start = start,
       logLike = -OPTIM$value, AIC = 2 * (OPTIM$value + 2))
}
