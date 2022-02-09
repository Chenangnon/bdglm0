#' @name ddlnorm
#' @aliases pdlnorm qdlnorm rdlnorm
#' @title The Balanced Discrete Log Normal Distribution
#' @description Probability mass function, cumulative distribution function, quantile function, and random variate generation,
#' for the balanced discrete log normal distribution with parameters \code{mu} (mean) and \code{alpha} (standard deviation at natural logarithmic scale).
#'
#' @usage
#' ddlnorm (x, mu, alpha = 1, meanlog = logb(mu) - .5 * (alpha^2), log = FALSE)
#' pdlnorm (q, mu, alpha = 1, meanlog = logb(mu) - .5 * (alpha^2), lower.tail = TRUE, log.p = FALSE)
#' qdlnorm (p, mu, alpha = 1, meanlog = logb(mu) - .5 * (alpha^2), lower.tail = TRUE, log.p = FALSE)
#' rdlnorm (n, mu, alpha = 1, meanlog = logb(mu) - .5 * (alpha^2))
#'
#' @param x,q vector of integer/real quantiles.
#' @param p vector of real probability values in the open \code{(0, 1)}.
#' @param mu vector of positive real mean parameters of the balanced discrete log-normal distribution.
#' @param alpha vector of positive real scale parameters of the balanced discrete distribution.
#' @param meanlog vector of positive reals, alternative way to specify the location of the balanced discrete distribution.
#' @param log,log.p logical; if TRUE, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @param n scalar, number of random variates to generate.
#'
#' @details The functions are the direct implementations of expressions derived from formulae in Tovissode et al. (2021).
#'
#' @return \code{ddlnorm} gives the probability mass, \code{pdlnorm} the distribution function, \code{qdlnorm} the quantile function, and \code{rdlnorm}
#' generates random deviates.
#'
#' The length of the result is determined by \code{n} for \code{rdlnorm}, and is the maximum of the lengths
#' of the numerical arguments for the other functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result. Only the first
#' elements of the logical arguments are used.
#'
#' @seealso \link{bdvariance} for variance computation and the estimation of distribution parameters by the method of moments;
#' and \link{bdlnorm} for generating balanced discrete log normal family objects for model fitting purposes.
#'
#' @references Tovissode, C.F.; Honfo, S.H.; Doumate, J.T.; Glele Kakai, R. On the Discretization of Continuous Probability Distributions Using a Probabilistic Rounding Mechanism. Mathematics 2021, 9, 555. https://doi.org/10.3390/math9050555

ddlnorm <- function(x, mu, alpha = 1,
                    meanlog = logb(mu) - .5 * (alpha^2),
                    log = FALSE) {
  if (!missing(meanlog))
    mu <- exp(meanlog + .5 * (alpha^2))
  zx <- x == 0
  x_ <- x - 1
  x_[zx] <- 0
  d <- (x - 1) * plnorm(x_, meanlog = meanlog, sdlog = alpha) -
    2 * x * plnorm(x, meanlog = meanlog, sdlog = alpha) +
    (x + 1) * plnorm(x + 1, meanlog = meanlog, sdlog = alpha) -
    mu * (pgammabar(x + 1, meanlog = meanlog, sdlog = alpha, r = 1) -
            2 * pgammabar(x, meanlog = meanlog, sdlog = alpha, r = 1) +
            pgammabar(x_, meanlog = meanlog, sdlog = alpha, r = 1))
  d <- pmax(d, exp(-708))
  if (any(zx0 <- ((x < 0) + (x != floor(x))) > 0))
    d[zx0] <- 0
  if (log)
    logb(d)
  else
    d
}

#' @name pgammabar

#' @title The pgammabar function
#' @description Compute the \code{pgammabar} function: \code{pnorm((logb(x) - meanlog)/sdlog - r * sdlog)}.
#'
#' @usage pgammabar (x, meanlog, sdlog = 1, r = 0)
#' @param meanlog,sdlog  real location (meanlog) and positive real scale (sdlog) parameters.

#' @details \code{pgammabar} as in Tovissode et al. (2021) except we use meanlog = log(mu) - .5 * (alpha^2) instead of mu.
#'
#' @return a vector, with length the maximum of the lengths of the arguments.
#'
#' @references Tovissode, C.F.; Honfo, S.H.; Doumate, J.T.; Glele Kakai, R. On the Discretization of Continuous Probability Distributions Using a Probabilistic Rounding Mechanism. Mathematics 2021, 9, 555. https://doi.org/10.3390/math9050555
#'
# pgammabar as in the core paper except we use meanlog = log(mu) - .5 * (alpha^2) instead of mu
# sdlog = alpha
pgammabar <- function(x, meanlog, sdlog = 1, r = 0) {
  plnorm(x * exp(- r * (sdlog^2)), meanlog = meanlog, sdlog = sdlog)
#  pnorm((logb(x) - meanlog)/sdlog - r * sdlog)
}

# Cumulative distribution function
pdlnorm <- function(q, mu, alpha = 1, meanlog = logb(mu) - .5 * (alpha^2),
                    lower.tail = TRUE, log.p = FALSE) {
  if (!missing(meanlog))
    mu <- exp(meanlog + .5 * (alpha^2))
  q <- floor(q)
  q1 <- q + 1
  pg <- plnorm(q, meanlog = meanlog, sdlog = alpha)
  p <- pg + q1 * (plnorm(q1, meanlog = meanlog, sdlog = alpha) - pg) -
    mu * (pgammabar(q1, meanlog = meanlog, sdlog = alpha, r = 1) -
            pgammabar(q, meanlog = meanlog, sdlog = alpha, r = 1))
  if (any(zq0 <- q < 0))
    p[zq0] <- 0
  if (!lower.tail)
    p <- 1 - p
  if(log.p)
    logb(p)
  p
}

# Quantile function
qdlnorm <- function(p, mu, alpha = 1, meanlog = logb(mu) - .5 * (alpha^2), lower.tail = TRUE, log.p = FALSE) {
  if(log.p)
    p <- exp(p)
  if (!lower.tail)
    p <- 1 - p
  x <- floor(qlnorm(p, meanlog = meanlog, sdlog = alpha, lower.tail = TRUE, log.p = FALSE))
  u <- pdlnorm(x, meanlog = meanlog, alpha = alpha)
  pu <- u < p
  x[pu] <- x[pu] + 1
  x
}

# Generate random deviates
rdlnorm <- function(n, mu, alpha = 1, meanlog = logb(mu) - .5 * (alpha^2)) {
  rg <- rlnorm (n, meanlog = meanlog, sdlog = alpha)
  z <- floor(rg)
  r <- rg - z
  nzr = r > 0
  if (any(nzr))
    z[nzr] <- z[nzr] + rbinom(sum(nzr), size = 1, prob = r[nzr])
  z
}

# variance function
.vardlnorm <- function(mu, alpha = 1, meanlog = logb(mu) - .5 * (alpha^2), tol = 1e-10) {
  if (!missing(meanlog))
    mu <- exp(meanlog + .5 * (alpha^2))
  na <- length(alpha)
  if (length(mu) == 1 && na == 1)
    return(as.vector((mu^2) * (exp(alpha^2) - 1) +
                       .lnzeta(meanlog = meanlog, alpha = alpha, tol = tol)))
  AB <- cbind(meanlog, alpha)
  meanlog <- AB[,1]
  alpha <- AB[,2]
  Lvalues <- apply(AB, MARGIN = 1, FUN = function(ab) {
    .lnzeta(meanlog = ab[1], alpha = ab[2], tol = tol)
  })
  as.vector((mu^2) * (exp(alpha^2) - 1) + Lvalues)
}

# Routine for .gzeta function a
.lnzeta_comps <- function(z, mu, alpha = 1, meanlog = logb(mu) - .5 * (alpha^2)) {
  comp1 <- sum((2 * z + 1) * (pgammabar(z + 1, meanlog = meanlog, sdlog = alpha, r = 1) -
                                pgammabar(z, meanlog = meanlog, sdlog = alpha, r = 1)))
  comp2 <- sum(pgammabar(z + 1, meanlog = meanlog, sdlog = alpha, r = 2) -
                 pgammabar(z, meanlog = meanlog, sdlog = alpha, r = 2))
  comp3 <- sum(z * (z + 1) * (plnorm(z + 1, meanlog = meanlog, sdlog = alpha) -
                                plnorm(z, meanlog = meanlog, sdlog = alpha)))
  c(comp1, comp2, comp3)
}

# Routine for .vardgamma function
.lnzeta <- function(mu, alpha = 1, meanlog = logb(mu) - .5 * (alpha^2), tol = 1e-10) {
  if (!missing(meanlog))
    mu <- exp(meanlog + .5 * (alpha^2))
  z_i <- floor(qlnorm(tol/2, meanlog = meanlog, sdlog = alpha))
  z_f <- floor(qlnorm(1 - tol/2, meanlog = meanlog, sdlog = alpha)) + 1
  if (is.infinite(z_i)) {
    k <- 0
    tol_i <- tol
    while(is.infinite(z_i) && k < 100) {
      tol_i <- 10 * tol_i
      z_i <- floor(qlnorm(tol_i/2, meanlog = meanlog, sdlog = alpha))
      k <- k + 1
    }
  }
  if (is.infinite(z_i))
    z_i <- max(0, mu - exp(4 * alpha))
  if (is.infinite(z_f)) {
    k <- 0
    tol_i <- tol
    while(is.infinite(z_f) && k < 100) {
      tol_i <- 10 * tol_i
      z_f <- floor(qlnorm(1 - tol_i/2, meanlog = meanlog, sdlog = alpha)) + 1
      k <- k + 1
    }
  }
  if (is.infinite(z_f))
    z_f <- mu + exp(4 * alpha)
  zcomps <- .lnzeta_comps(z_i:z_f, meanlog = meanlog, alpha = alpha)
  res <- zcomps[1] * mu - zcomps[2] * (mu^2) * exp(alpha^2) - zcomps[3]
  structure(res, precision = plnorm(z_i, meanlog = meanlog, sdlog = alpha) +
              plnorm(z_f + 1, meanlog = meanlog, sdlog = alpha, lower.tail = FALSE))
}

.momente_bdlnorm <-
  function(x, zeta = NULL) {
    if (missing(zeta) || is.null(zeta))
      .m.estdln (x)
    else
      .m.estdln (x, zeta = zeta)
  }

# moment estimates (approximate for alpha)
.m.estdln <- function(x, zeta = if (sigma2 > minz) minz/2 else sigma2/4) {
  xbar <- mean(x)
  sigma2 <- var(x)
  minz <- min(xbar, .25)
  alpha <- sqrt(log1p((sigma2  - zeta) / (xbar^2)))
  c(mu = xbar, alpha = alpha, sigma2 = sigma2)
}

# Find alpha such that the balanced discrete log normal distribution with mean mu is equidispersed
.equialpha_bdlnorm <-
  function (mu, tol = 1e-5) {
  tol <- min(1e-5, max(1e-15, tol))
  lfun <- function(mui) {
    if (mui == 0)
      return(1)
    fun <- function(x)
      (mui^2) * (exp(x^2) - 1) + .lnzeta(mu = mui, alpha = x, tol = tol) - mui
    upper <- sqrt(log1p(1/mui))
    zmui <- mui > .25
    if (zmui)
      lower <- sqrt(log1p(1/mui - 1/(4 * mui^2)))
    else
      lower <- upper/4
    res <- uniroot(f = fun, lower = lower, upper = upper,
                   extendInt = if (zmui) "no" else "yes", tol = tol)
    structure(as.vector(res$root), diffEV = as.vector(res$f.root))
  }
  if (length(mu) == 1)
    lfun(mu)
  else
    sapply(mu, FUN = lfun)
}

# Fit a discrete log normal distribution to an IDD sample x
.mle_bdlnorm <- function(x, theta = NULL, alpha = .05, asymp = "norm",
                       method = "BFGS", tol = 1e-10, maxit = 100) {
  if (missing(theta) || is.null(theta))
    start <- logb(.momente_bdlnorm(x)[1:2])
  else {
    start <- logb(theta[1:2])
    names(start) <- c("mu", "alpha")
  }
  # Negative log likelihood function: theta = c(log(mu), log(alpha))
  nlike.dlnorm <- function(theta, x) {
    -sum(ddlnorm(x, mu = exp(theta[1]), alpha = exp(theta[2]), log = TRUE))
  }
  OPTIM <- optim(par = start, fn = nlike.dlnorm, x = x,
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
  names(theta.corr) <- "mu:alpha"
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
