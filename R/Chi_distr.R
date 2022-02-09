#' @name dchi
#' @aliases pchi qchi rchi
#' @title The (non-central) Chi distribution
#' @description Density, distribution function, quantile function and random generation for
#' the Chi (KI) distribution with degrees of freedom df and non centrality parameter ncp.
#' @usage  dchi (x, df = 1, ncp = 0, log = FALSE)
#' pchi (q, df = 1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
#' qchi (p, df = 1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
#' rchi (n, df = 1, ncp = 0)
#'
#' @param x,q vector of quantiles.
#' @param p vector of real values in the close (0, 1).
#' @param n scalar; number of deviates to generate.
#' @param df vector of positive real values, degrees of freedom parameter.
#' @param ncp non centrality parameter (non-negative).
#' @param log,log.p logical; if TRUE, probabilities/densities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @seealso See \link{Half-Chi} for the foldded (absolute value of the) chi distribution.

dchi <- function(x, df = 1, ncp = 0, log = FALSE) {
  zerox <- x == 0
  if (any(zerox) && length(df) > 1) {
    n <- length(x)
    df <- rep(df, n)
    d <- numerix(n)
    if (any(df1 <- df[zerox] > 1))
      d[zerox][df1] <-  0
    if (any(df1 <- df[zerox] == 1))
      d[zerox][df1] <-  dnorm(0)
    if (any(df1 <- df[zerox] < 1))
      d[zerox][df1] <- Inf
    if (log)
      d[zerox] <- logb(d[zerox])
    if (any(!zerox)) {
      ncp <- rep(ncp, n)
      d[!zerox] <- dchisq(x = x[!zerox]^2,
                          df = df[!zerox], ncp = ncp[!zerox], log = log)
      if (log)
        d[!zerox] <- d[!zerox] + logb(abs(x[!zerox]))
      else
        d[!zerox] <- d[!zerox] * abs(x[!zerox])
    }
    return(d)
  }
  d <- abs(x) * dchisq(x = x^2, df = df, ncp = ncp)
  if (any(zerox))
    d[zerox] <- if (df == 1) dnorm(0) else if (df < 1) Inf else 0
  if (log)
    logb(d)
  else
    d
}

pchi <- function(q, df = 1, ncp = 0, lower.tail = TRUE, log.p = FALSE) {
  u <- sign(q)
  p <- .5 * (1 + u * pchisq(q^2, df = df, ncp = ncp))
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    logb(p)
  else
    p
}

qchi <- function(p, df = 1, ncp = 0, lower.tail = TRUE, log.p = FALSE) {
  if (log.p)
    p <- exp(p)
  if (!lower.tail)
    p <- 1 - p
  p <- p - .5
  q <- sign(p) * sqrt(qchisq(2 * abs(p), df = df, ncp = ncp))
  q
}

rchi <- function(n, df = 1, ncp = 0) {
  if(df == 1)
    return(rnorm(n, mean = ncp))
  s <- sqrt(rchisq(n, df = df, ncp = ncp))
  u <- runif(n) < 0.5
  if(sum(u) > 0)
    s[u] <- -s[u]
  s
}

#' @name Half-Chi
#' @aliases dhchi phchi qhchi rhchi
#' @title The (non-central) half chi distribution
#' @description Density, distribution function, quantile function and random generation for
#' the half Chi (KI) distribution with degrees of freedom df and non centrality parameter ncp.
#'
#' @usage dhchi (x, df = 1, ncp = 0, log = FALSE)
#'        phchi (q, df = 1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
#'        qhchi (p, df = 1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
#'        rhchi (n, df = 1, ncp = 0)
#'
#' @param x,q vector of non negative quantiles.
#' @param p vector of real values in the close (0, 1).
#' @param n scalar; number of deviates to generate.
#' @param df vector of positive real values, degrees of freedom parameter.
#' @param ncp non centrality parameter (non-negative).
#' @param log,log.p logical; if TRUE, probabilities/densities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @seealso See \link{Chi} for the complete \code{chi} distribution.
#'
dhchi <- function(x, df = 1, ncp = 0, log = FALSE) {
  d <- 2 * x * dchisq(x = x^2, df = df, ncp = ncp)
  d[x < 0] <- 0
  d[x == 0] <- if (df == 1) dnorm(0) else if (df < 1) Inf else 0
  if (log)
    logb(d)
  else
    d
}

phchi <- function(q, df = 1, ncp = 0, lower.tail = TRUE, log.p = FALSE) {
  p <- pchisq(q^2, df = df, ncp = ncp)
  p[q < 0] <- 0
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    logb(p)
  else
    p
}

qhchi <- function(p, df = 1, ncp = 0, lower.tail = TRUE, log.p = FALSE) {
  if (log.p)
    p <- exp(p)
  if (!lower.tail)
    p <- 1 - p
  sqrt(qchisq(p, df = df, ncp = ncp))
}

rhchi <- function(n, df = 1, ncp = 0) {
  sqrt(rchisq(n, df = df, ncp = ncp))
}
