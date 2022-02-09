# See package "rmutil" for Generalized Gamma distribution (continuous) routines

# Probability mass function
# x = x values (quantiles) which the probability masses are desired at
# mu = mean value
# b,c = shape parameters; b = mu (default) yields a = 1 (unit scale) and variance = mu + zeta, i.e. mu < var <= mu + 1/4
# log logical

# Add code to check validity of mu, b and c (all positive) and x (non negative integer)

#' @name ddggamma
#' @aliases pdggamma qdggamma rdggamma
#' @title The Balanced Discrete Generalized Gamma Distribution
#' @description Probability mass function, cumulative distribution function, quantile function, and random variate generation for
#' the balanced discrete generalized gamma distribution with parameters \code{mu} (mean), \code{b} (scale) and \code{c} (shape).
#'
#' @usage
#' ddggamma (x, mu, b = mu, c = 1, log = FALSE)
#' pdggamma (q, mu, b = mu, c = 1, lower.tail = TRUE, log.p = FALSE)
#' qdggamma (p, mu, b = mu, c = 1, lower.tail = TRUE, log.p = FALSE)
#' rdggamma (n, mu, b = mu, c = 1)
#'
#' @param x,q vector of integer/real quantiles.
#' @param p vector of real probability values in the open \code{(0, 1)}.
#' @param n scalar, number of deviates to generate.
#' @param mu vector of positive real mean parameters of the balanced discrete distribution.
#' @param b,c vectors of positive real scale (\code{b}) and shape (\code{c}) parameters of the balanced discrete distribution.
#' @param log.p,log logical; if TRUE, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @details The functions are the direct implementations of formulae in Tovissode et al. (2021).
#'
#' @return \code{ddggamma} gives the probability mass, \code{pdggamma} the distribution function, \code{qdggamma} the quantile function, and \code{rdggamma}
#' generates random deviates. Invalid arguments will result in errors.
#'
#' The length of the result is determined by \code{n} for \code{rdggamma}, and is the maximum of the lengths
#' of the numerical arguments for the other functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result. Only the first
#' elements of the logical arguments are used.
#'
#' @references Tovissode, C.F.; Honfo, S.H.; Doumate, J.T.; Glele Kakai, R. On the Discretization of Continuous Probability Distributions Using a Probabilistic Rounding Mechanism. Mathematics 2021, 9, 555. https://doi.org/10.3390/math9050555
#'
#' @seealso \link{bdggamma} for generating balanced discrete gamma family objects for model fitting purposes.

ddggamma <- function(x, mu, b = mu, c = 1, log = FALSE) {
  a <- gamma(b+1/c) / (mu * gamma(b))
  zx <- x == 0
  ax_ <- (a * (x - 1))^c
  ax_[zx] <- 0
  ax <- (a * x)^c
  axp <- (a * (x + 1))^c
  d <- (x - 1) * pgamma(ax_, shape = b) -
    2 * x * pgamma(ax, shape = b) +
    (x + 1) * pgamma(axp, shape = b) -
    mu * (pgamma(axp, shape = b + 1/c) -
            2 * pgamma(ax, shape = b + 1/c) +
            pgamma(ax_, shape = b + 1/c))
  d <- pmax(d, exp(-708))
  if (any(zx0 <- ((x < 0) + (x != floor(x))) > 0))
    d[zx0] <- 0
  if (log)
    logb(d)
  else
    d
}

# Cumulative distribution function
pdggamma <- function(q, mu, b = mu, c = 1, lower.tail = TRUE, log.p = FALSE) {
  a <- gamma(b+1/c) / (mu * gamma(b))
  q <- floor(q)
  aq <- (a * q)^c
  q1 <- q + 1
  aqp <- (a * q1)^c
  pg <- pgamma(aq, shape = b)
  p <- pg + q1 * (pgamma(aqp, shape = b) - pg) -
    mu * (pgamma(aqp, shape = b + 1/c) - pgamma(aq, shape = b + 1/c))
  if (any(zq0 <- q < 0))
    p[zq0] <- 0
  if (!lower.tail)
    p <- 1 - p
  if(log.p)
    logb(p)
  p
}

# Quantile function
qdggamma <- function(p, mu, b = mu, c = 1, lower.tail = TRUE, log.p = FALSE) {
  a <- gamma(b+1/c) / (mu * gamma(b))
  if(log.p)
    p <- exp(p)
  if (!lower.tail)
    p <- 1 - p
  x <- floor((qgamma(p, shape = b, lower.tail = TRUE, log.p = FALSE)^(1/c)) / a)
  u <- pdggamma(x, mu = mu, b = b, c = c)
  pu <- u < p
  x[pu] <- x[pu] + 1
  x
}

# Generate random deviates
rdggamma <- function(n, mu, b = mu, c = 1) {
  rg <- (rgamma(n, shape = b, rate = 1)^(1/c)) * (mu * gamma(b)) / gamma(b+1/c)
  z <- floor(rg)
  r <- rg - z
  nzr = r > 0
  if (any(nzr))
    z[nzr] <- z[nzr] + rbinom(sum(nzr), size = 1, prob = r[nzr])
  z
}

# Routine for .ggzeta function
.ggzeta_comps <- function(z, b, c = 1, a = 1) {
  nz <- length(z)
  if (nz > 1000) {
    init <- seq(1, nz, by = 1000)
    nt <- length(init)
    final <- init + 1000 - 1
    final[nt] <- nz
    comp1 <- comp2 <- comp3 <- 0
    for (j in 1:nt) {
      zj <- z[init[j]:final[j]]
      comp1 <- comp1 +
        sum((2 * zj + 1) * (pgamma((a * (zj + 1))^c, shape = b + 1/c) -
                             pgamma((a * zj)^c, shape = b + 1/c)))
      comp2 <- comp2 +
        sum(pgamma((a * (zj + 1))^c, shape = b + 2/c) -
              pgamma((a * zj)^c, shape = b + 2/c))
      comp3 <- comp3 +
        sum(zj * (zj + 1) * (pgamma((a * (zj + 1))^c, shape = b) -
                             pgamma((a * zj)^c, shape = b)))
    }
    return(c(comp1, comp2, comp3))
  }
  comp1 <- sum((2 * z + 1) * (pgamma((a * (z + 1))^c, shape = b + 1/c) -
                                pgamma((a * z)^c, shape = b + 1/c)))
  comp2 <- sum(pgamma((a * (z + 1))^c, shape = b + 2/c) -
                 pgamma((a * z)^c, shape = b + 2/c))
  comp3 <- sum(z * (z + 1) * (pgamma((a * (z + 1))^c, shape = b) -
                                pgamma((a * z)^c, shape = b)))
  c(comp1, comp2, comp3)
}

# Routine for vardgamma function
.ggzeta <- function(mu, b = mu, c = 1, tol = 1e-6) {
  a <- exp(lgamma(b + 1/c) - logb(mu) - lgamma(b))
  if (any(c(is.na(a), is.infinite(a))))
    return(NaN)
  z_i <- floor((qgamma(tol/2, shape = b)^(1/c)) / a)
  z_f <- floor((qgamma(1 - tol/2, shape = b)^(1/c)) / a) + 1
  zcomps <- .ggzeta_comps(z_i:z_f, b = b, c = c, a = a)
  res <-  mu * zcomps[1] - (mu^2) * (gamma(b) * gamma(b + 2/c) / (gamma(b + 1/c)^2)) * zcomps[2] - zcomps[3]
  structure(res, precision = pgamma((a * z_i)^c, shape = b) + pgamma((a * (z_f + 1))^c, shape = b, lower.tail = FALSE))
}

# var functions
.vardggamma <- function(mu, b = mu, c = 1, tol = 1e-6) {
  na <- max(length(b), length(mu))
  if (length(c) == 1 && na == 1) {
    tryzval <- try(.ggzeta(mu = mu, b = b, c = c, tol = tol), silent = TRUE)
    if (class(tryzval) == "try-error")
      return(NaN)
    return(as.vector((gamma(b) * gamma(b+2/c) / (gamma(b+1/c)^2) - 1) * mu^2 +
                tryzval))
  }
  ABC <- cbind(mu, b, c)
  mu <- ABC[,1]
  b <- ABC[,2]
  c <- ABC[,3]
  GGvalues <- apply(ABC, MARGIN = 1, FUN = function(abc) {
    tryzval <- try(.ggzeta(mu = abc[1], b = abc[2], c = abc[3], tol = tol), silent = TRUE)
    if (class(tryzval) == "try-error")
      return(NaN)
    tryzval
  })
  as.vector((gamma(b) * gamma(b+2/c) / (gamma(b+1/c)^2) - 1) * mu^2 +
              GGvalues)
}

# Find b such that the balanced discrete generalized gamma distribution with mean mu and c = c is equidispersed
.equib.bdggamma <- function (mu, c = 1, tol = 1e-6) {
  tol <- min(1e-5, max(1e-15, tol))
  lfun <- function(mui, c) {
    if (mui == 0)
      return(tol)
    fun <- function(x) {
      varval <- vardggamma(mu = mui, b = exp(x), c = c, tol = tol)
      if(any(c(is.na(varval), is.infinite(varval))))
        return(100*mui)
      (varval - mui)^2
    }
    res <- try(uniroot(f = fun,
                   lower = logb(mui), upper = logb(mui)+logb(2), extendInt = "yes", tol = tol),
               silent = TRUE)
    if (class(res) == "try-error")
      structure(NA, diffEV = NA, message = res)
    else
      structure(as.vector(exp(res$root)), diffEV = sqrt(as.vector(res$f.root)))
  }
  if ((nmu <- length(mu)) == 1)
    lfun(mu, c[1])
  else {
    if (length(c) == 1)
      sapply(mu, FUN = function(mui) lfun(mui, c))
    else
      sapply(1:nmu, FUN = function(i) lfun(mu[i], c[i]))
  }
}

# Find c such that the balanced discrete generalized gamma distribution with mean mu and b = b is equidispersed
.equic.bdggamma <- function (mu, b = mu, tol = 1e-6) {
  tol <- min(1e-5, max(1e-15, tol))
  lfun <- function(mui, b) {
    if (mui == 0)
      return(1)
    fun <- function(x) vardggamma(mu = mui, b = b, c = exp(x), tol = tol)  - mui
    res <- try(uniroot(f = fun,
                   lower = logb(.5), upper = logb(5), extendInt = "yes", tol = tol),
               silent = TRUE)
    if (class(res) == "try-error")
      structure(NA, diffEV = NA, message = res)
    else
      structure(as.vector(exp(res$root)), diffEV = as.vector(res$f.root))
  }
  if ((nmu <- length(mu)) == 1)
    lfun(mu, b[1])
  else {
    if (length(b) == 1)
      sapply(mu, FUN = function(mui) lfun(mui, b))
    else
      sapply(1:nmu, FUN = function(i) lfun(mu[i], b[i]))
  }
}

# moment estimates very approximate for c, often c = 1 !!!
# Uses P(Y = 0) = proportion of zeros in x ??
.momente_bdggamma <- function(x, zeta = NULL, tol = 1e-5, maxit = 10) {
  xbar <- mean(x)
  sigma2 <- var(x)
  if (is.null(zeta)) {
    minz <- min(xbar, .25)
    zeta <- if (sigma2 > minz) minz/2 else sigma2/4
  }
  Const <- ((sigma2 - zeta)/xbar^2) + 1
  py0 <- mean(x==0)
  fun <- function(lpar) {
    b <- exp(lpar[1])
    invc <- exp(-lpar[2])
    a <- gamma(b + invc) / (xbar * gamma(b))
    ac <- a^(1/invc)
    (gamma(b) * gamma(b+2* invc)/gamma(b+invc)^2 - Const)^2 +
      (py0 - pgamma(ac, shape = b, rate = 1) + xbar * pgamma(ac, shape = b + invc, rate = 1))^2
  }
  b <- xbar^2 / (sigma2  - zeta)
  c <- 1
  res <- optim(par = logb(c(b, c)), fn = fun,
               method = "BFGS", control = list(reltol = tol, maxit = maxit))
  c(mu = xbar, b = exp(res$par[1]), c = exp(res$par[2]), sigma2 = sigma2)
}

# Fit a discrete generalized gamma distribution to an IDD sample x
.mle_bdggamma <- function(x, theta = NULL, alpha = .05, asymp = "norm",
                          method = "BFGS", tol = 1e-10, maxit = 100) {
  if (missing(theta) || is.null(theta))
    start <- logb(.m.estdggamma(x)[1:3])
  else {
    start <- logb(theta[1:3])
  }
  # Negative log likelihood function: theta = c(log(mu), log(b), log(c))
  nlike.dggamma <- function(theta, x) {
    -sum(ddggamma(x, mu = exp(theta[1]), b = exp(theta[2]), c = exp(theta[3]), log = TRUE), na.rm = TRUE)
  }
  names(start) <- c("mu", "b", "c")
  OPTIM <- optim(par = start, fn = nlike.dggamma, x = x,
                 method = method, hessian = TRUE,
                 control = list(maxit = maxit))
  theta <- exp(OPTIM$par)
  log.theta.se <- try(solve(OPTIM$hessian), silent = TRUE)
  if (class(log.theta.se) == "try-error") {
    return(list(theta = theta, CI = NULL, se = NULL, df = length(x)-2,
                  log.se = log.theta.se, cor = NULL, start = exp(start),
                  logLike = -OPTIM$value, AIC = 2 * (OPTIM$value + 2)))
  }
  theta.corr <- cov2cor(log.theta.se)[c(2, 3, 6)]
  names(theta.corr) <- c("mu:b", "mu:c", "b:c")
  log.theta.se <- sqrt(diag(log.theta.se))
  theta.se <- sqrt(theta) * log.theta.se
  df <- length(x) - 3
  if (identical(asymp, "norm"))
    qalpha <- qnorm(alpha/2, lower.tail = FALSE)
  else
    qalpha <- qt(alpha/2, df = df, lower.tail = FALSE)
  half.CI <- qalpha * log.theta.se
  theta.CI <- cbind(lower = exp(OPTIM$par - half.CI), upper = exp(OPTIM$par + half.CI))
  # Add a Test for equidispersion: mu == sigma2
  list(theta = theta, CI = theta.CI, se = theta.se, df = df,
       log.se = log.theta.se, cor = theta.corr, start = exp(start),
       logLike = -OPTIM$value, AIC = 2 * (OPTIM$value + 3))
}
