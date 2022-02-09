#' @name dinvgauss
#' @aliases pinvgauss qinvgauss rinvgauss
#' @title The Inverse Gaussian Distribution
#' @description Density, distribution function, quantile function and random generation for the inverse Gaussian distribution with mean equal to \code{mu} and scale parameter \code{a}.
#'
#' @param x,q vector of integer/real quantiles.
#' @param p vector of real probability values in the open \code{(0, 1)}, or the \code{log} of such values when \code{log.p = TRUE}.
#' @param mu vector of positive real mean values.
#' @param a vector of positive real scale parameters.
#' @param log.p,log logical; if TRUE, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @param n integer, number of random deviates to generate.
#'
#' @usage
#' dinvgauss (x, mu, a = mu^2, log = FALSE)
#'
#' pinvgauss (q, mu, a = mu^2, lower.tail = TRUE, log.p = FALSE)
#'
#' qinvgauss (p, mu, a = mu^2, lower.tail = TRUE, log.p = FALSE)
#'
#' rinvgauss (n, mu, a = mu^2)
#'
#' @details
#' The inverse Gaussian distribution has density
#'
#'     \code{f(x) = sqrt(a)/sqrt(2 pi x^3) exp{-(a (x - m)^2/(2 x m^2))}}
#'
#'  where \code{mu} is the mean of the distribution and \code{a} is a scale
#'  parameter such that the variance of the distribution is equal to
#'   \code{mu^3 / a} (Folks & Chhikara, 1978).
#'
#' \code{dinvgauss} gives the probability mass, \code{pinvgauss} the distribution function,
#' \code{qinvgauss} the quantile function, and \code{rinvgauss} generates random variates.
#'
#' The length of the result is determined by \code{n} for \code{rinvgauss}, and is the maximum of the lengths
#' of the numerical arguments for the other functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result. Only the first
#' elements of the logical arguments are used.
#'
#' @source The density and distribution functions (\code{dinvgauss} and \code{pinvgauss}) are based on formulas in Folks & Chhikara (1978).
#'
#' \code{qinvgauss} simply uses the routine \link{uniroot} to reverse \code{pinvgauss}.
#'
#' \code{rinvgauss} is based on Michael et al. (1976).
#'
#' @references Folks, J. L., & Chhikara, R. S. (1978). The inverse Gaussian distribution and its statistical application-a review. Journal of the Royal Statistical Society: Series B (Methodological), 40(3), 263-275.
#'
#' Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating random variates using transformations with multiple roots. The American Statistician, 30(2), 88-90.
#'
#' @seealso \link{ddinvgauss} for the balanced discrete inverse Gaussian distribution.
#'

dinvgauss <-
  function (x, mu, a = mu^2, log = FALSE) {
    x <- cbind(x, mu, a)
    a <- x[,3]
    mu <- x[,2]
    x <- x[,1]
    posx <- x > 0
    invalid0 <- invalid <- !is.finite(mu) | !is.finite(a) | (mu < 0) | (a <= 0) | is.na(mu+a)
    if (any(posx))
      invalid[posx] <- !is.finite(x[posx])
    if (any(!posx))
      invalid[!posx] <- x[!posx] < 0
    invalid <- invalid | invalid0
    logpdf <- x - Inf
    if (any(!invalid)) {
      logpdf[!invalid] <-
        - a[!invalid] * ((x[!invalid] - mu[!invalid])^2) /
        (2 * x[!invalid] * mu[!invalid]^2) +
        .5 * (logb(a[!invalid]) - logb(2 * pi) - 3 * logb(x[!invalid]))
    }
    logpdf[invalid0] <- NaN
    logpdf[(x==0) + (mu == 0) == 2] <- 0
    return(
      if (log) logpdf
      else exp(logpdf))
}

pinvgauss <-
  function (q, mu, a = mu^2, lower.tail = TRUE, log.p = FALSE) {
    q <- cbind(q, mu, a)
    a <- q[,3]
    mu <- q[,2]
    q <- q[,1]
    posx <- q > 0
    invalid0 <- invalid <- !is.finite(mu) | !is.finite(a) | (mu < 0) | (a <= 0) | is.na(mu+a)
    if (any(posx))
      invalid[posx] <- !is.finite(q[posx])
    if (any(!posx))
      invalid[!posx] <- q[!posx] < 0
    invalid <- invalid | invalid0
    cdf <- numeric(length(q))
    if (any(!invalid)) {
      r1 <- q[!invalid] / mu[!invalid]
      r2 <- sqrt(a[!invalid] / q[!invalid])
      cdf[!invalid] <- pnorm(r2 * (r1 - 1)) +
        exp(2 * a[!invalid] / mu[!invalid]) * pnorm(-r2 * (r1 + 1))
    }
    cdf[invalid0] <- NaN
    cdf[(q>=0) + (mu == 0) == 2] <- 1
    if(!lower.tail)
      cdf <- 1 - cdf
    if (log.p)
      cdf <- logb(cdf)
    return(cdf)
  }

qinvgauss <-
  function (p, mu, a = mu^2, lower.tail = TRUE, log.p = FALSE) {
    pmua <- cbind(p, mu, a)
    a <- pmua[,3]
    mu <- pmua[,2]
    p <- pmua[,1]
    if (log.p)
      p <- exp(p)
    if (!lower.tail)
      p <- 1 - p
    invalid0 <- is.na(p)
    if (any(!invalid0))
      invalid0[!invalid0] <- p[!invalid0] <= 0 | p[!invalid0] >= 1
    invalidpars <- !is.finite(mu) | !is.finite(a) | (mu < 0) | (a <= 0) | is.na(mu+a)
    invalid <- invalid0 | invalidpars
    qres <- p + NaN
    if (any(!invalid)) {
      qexpr <- function(q, pmuai) {
        r1 <- q/pmuai[2]
        r2 <- sqrt(pmuai[3] / q)
        pnorm(r2 * (r1 - 1)) +
          exp(2 * pmuai[3]/pmuai[2]) * pnorm(- r2 * (r1 + 1)) -
          pmuai[1]
      }
      b <- min(a[!invalid]/mu[!invalid], 300)
      b <- .5 + exp(2 * b) *
        pnorm(- sqrt(a[!invalid]/mu[!invalid]) * 2)
      b <- p[!invalid] > b
      qres[!invalid] <-
        apply(cbind(pmua[!invalid,, drop = FALSE], b),
              MARGIN = 1,
              FUN = function(pmuai) {
                if (pmuai[4])
                  interval_i <- c(pmuai[2], max(20, pmuai[2]+1))
                else
                  interval_i <- c(.Machine$double.xmin, pmuai[2])
                res <- tryCatchWE(uniroot(f = function(q) qexpr (q, pmuai),
                                   interval = interval_i,
                                   extendInt = "yes",
                                   tol = .Machine$double.eps^0.6))$value
                if(any(class(res) %in% c("try-error", "condition", "error")))
                  return(NaN)
                res$root
              })
    }
    invalid0 <- invalid0 & !is.na(p)
    if (any(invalid0)) {
      if (any(select <- ((p[invalid0] == 0) + (invalidpars[invalid0] == 0)) == 2))
        qres[invalid0][select] <- -Inf
      if (any(select <- (p[invalid0] == 1) + (mu[invalid0] == 0) == 2))
        qres[invalid0][select] <- 0
      if (any(select <- (p[invalid0] == 1 + !invalidpars[invalid0]) == 2))
        qres[invalid0][select] <- Inf
    }
    return(as.vector(qres))
  }

rinvgauss <-
  function(n, mu, a = mu^2) {
    if (length(n) > 1)
      n <- length(n)
    mu <- rep_len(mu, n)
    a <- rep_len(a, n)
    valid0 <- (mu == 0) & (a > 0)
    invalid <- !is.finite(mu) | !is.finite(a) | (mu < 0) | (a <= 0) | is.na(mu+a)
    invalid <- valid0 | invalid
    xres <- mu + NaN
    nvalid <- sum(!invalid)
    if (nvalid > 0) {
      nu0 <- rchisq(nvalid, df = 1) # qchisq(runif(nvalid), df = 1)
      r0 <- mu[!invalid] / a[!invalid]
      x1 <- mu[!invalid] + .5 * nu0 * mu[!invalid] * r0 -
        mu[!invalid] * sqrt(nu0 * r0 + .25 * (r0 * nu0)^2)
      U <- rbinom(nvalid, size = 1, prob = mu[!invalid] / (x1 + mu[!invalid])) == 1
      if (any(U))
        xres[!invalid][U] <- x1[U]
      if (any(!U))
        xres[!invalid][!U] <- (mu[!invalid][!U]^2) / x1[!U]
    }
    if (any(valid0))
      xres[valid0] <- 0
    return(xres)
  }

# Trying an analytical version of 'qinvgauss'
.qinvgauss0 <-
  function (p, mu, a = mu^2, lower.tail = TRUE, log.p = FALSE) {
    pmua <- cbind(p, mu, a)
    a <- pmua[,3]
    mu <- pmua[,2]
    p <- pmua[,1]
    if (log.p)
      p <- exp(p)
    if (!lower.tail)
      p <- 1 - p
    invalid <- p <= 0 | p >= 1
    invalid0 <- !is.finite(mu) | !is.finite(a) | (mu < 0) | (a <= 0)
    invalid <- invalid | invalid0
    qres <- p - Inf
    if (any(!invalid)) {
      nu0 <- qchisq(p[!invalid], df = 1)
      r0 <- mu[!invalid] / a[!invalid]
      x1 <- mu[!invalid] + .5 * nu0 * mu[!invalid] * r0 -
        mu[!invalid] * sqrt(nu0 * r0 + .25 * (r0 * nu0)^2)
      r1 <- x1 / mu[!invalid]
      r2 <- sqrt(a[!invalid] / x1)
      sup1 <- (pnorm(r2 * (r1 - 1)) +
                 exp(2 * a[!invalid] / mu[!invalid]) * pnorm(-r2 * (r1 + 1))) >= p[!invalid]

      print(sup1)
      print(cbind(p[!invalid], x1, sup1))

      if (any(sup1))
        qres[!invalid][sup1] <- x1[sup1]
      if (any(!sup1))
        qres[!invalid][!sup1] <- mu[!invalid][!sup1]^2/x1[!sup1]
    }
    if (any(invalid0))
      qres[invalid0] <- NaN
    if (any(zeromax <- (p == 1) + (mu == 0) == 2))
      qres[zeromax] <- 0
    return(qres)
  }

