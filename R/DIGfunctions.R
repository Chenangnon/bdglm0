#' @name ddinvgauss
#' @aliases pdinvgauss qdinvgauss rdinvgauss
#' @title The Balanced Discrete Inverse Gaussian Distribution
#' @description Probability mass function, cumulative distribution function, quantile function, and random variate generation,
#' for the balanced discrete inverse Gaussian distribution with parameters \code{mu} (mean) and \code{a} (scale).
#'
#' @usage
#' ddinvgauss (x, mu, a = mu^2, log = FALSE)
#' pdinvgauss (q, mu, a = mu^2, lower.tail = TRUE, log.p = FALSE)
#' qdinvgauss (p, mu, a = mu^2, lower.tail = TRUE, log.p = FALSE)
#' rdinvgauss (n, mu, a = mu^2)
#'
#' @param x,q vector of integer/real quantiles.
#' @param p vector of real probability values in the open \code{(0, 1)}.
#' @param mu vector of positive real mean parameters of the balanced discrete log-normal distribution.
#' @param a vector of positive real scale parameters of the balanced discrete distribution.
#' @param log,log.p logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @param n scalar, number of random variates to generate.
#'
#' @details The functions are the direct implementations of expressions derived from formulae in Tovissode et al. (2021).
#'
#' @return \code{ddinvgauss} gives the probability mass, \code{pdinvgauss} the distribution function, \code{qdinvgauss} the quantile function, and \code{rdinvgauss}
#' generates random deviates.
#'
#' The length of the result is determined by \code{n} for \code{rdlnorm}, and is the maximum of the lengths
#' of the numerical arguments for the other functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result. Only the first
#' elements of the logical arguments are used.
#'
#' @seealso \link{bdvariance} for variance computation and the estimation of distribution parameters,
#' \link{dinvgauss} for the density, cumulative distribution, and quantile functions as well as random variate generation for the continuous inverse Gaussian distribution,
#' and \link{bdinvgauss} for generating balanced discrete log normal family objects for model fitting purposes.
#'
#' @references Tovissode, C.F.; Honfo, S.H.; Doumate, J.T.; Glele Kakai, R. On the Discretization of Continuous Probability Distributions Using a Probabilistic Rounding Mechanism. Mathematics 2021, 9, 555. https://doi.org/10.3390/math9050555
#'
#' @examples
#' # Compare continuous Inverse Gaussian (IG) to discrete IG
#' ## Density and probability mass
#' curve(dinvgauss(x, mu = 6, a = 5), from = 0, to = 20, col = "blue")
#' curve(ddinvgauss(x, mu = 6, a = 5), from = 0, to = 20, add = TRUE, col = "red", type = "h")
#'
#' ## Cumulative probabilities
#' curve(pinvgauss(x, mu = 6, a = 5), from = 0, to = 20, col = "blue")
#' curve(pdinvgauss(x, mu = 6, a = 5), 0, 20, add = TRUE, col = "red", type = "h")
#'

# Probability mass function
ddinvgauss <-
  function(x, mu, a = mu^2, log = FALSE) {
  zx <- x == 0
  x_ <- x - 1
  x_[zx] <- 0
  d <- (x - 1) * pinvgauss(x_, mu = mu, a = a) -
    2 * x * pinvgauss(x, mu = mu, a = a) +
    (x + 1) * pinvgauss(x + 1, mu = mu, a = a) -
    .diffE1_invgauss(y = x, mu = mu, a = a)
  d <- pmax(d, exp(-708))
  if (any(zx0 <- ((x < 0) + (x != floor(x))) > 0))
    d[zx0] <- 0
  if (log)
    logb(d)
  else
    d
}

# Cumulative probabilities
pdinvgauss <-
  function(q, mu, a = mu^2, lower.tail = TRUE, log.p = FALSE) {
    q <- floor(q)
    p <- (q + 1) * pinvgauss(q + 1, mu = mu, a = a) -
      q * pinvgauss(q, mu = mu, a = a) -
      .E1_invgauss (q, mu = mu, a = a)
    p <- pmax(p, exp(-708))
    if (any(zx0 <- ((q < 0) + (q != floor(q))) > 0))
      p[zx0] <- 0
    if (!lower.tail[1])
      p <- 1 - p
    if (log.p[1])
      logb(p)
    else
      p
  }

# Quantile function
qdinvgauss <- function(p, mu, a = mu^2, lower.tail = TRUE, log.p = FALSE) {
  if (log.p[1])
    p <- exp(p)
  if (!lower.tail[1])
    p <- 1 - p
  x <- floor(qinvgauss(p, mu = mu, a = a, lower.tail = TRUE, log.p = FALSE))
  invalidp <- is.na(x)
  if (any(invalidp))
    x[invalidp] <- -1
  u <- pdinvgauss(x, mu = mu, a = a, lower.tail = TRUE, log.p = FALSE)
  pu <- u < p
  x[pu] <- x[pu] + 1
  if (any(invalidp))
    x[invalidp] <- NaN
  x
}

# Generate random deviates
rdinvgauss <- function (n, mu, a = mu^2) {
  rg <- rinvgauss(n, mu = mu, a = a)
  z <- floor(rg)
  r <- rg - z
  nzr <- r > 0
  if (any(nzr))
    z[nzr] <- z[nzr] + rbinom(sum(nzr), size = 1, prob = r[nzr])
  z
}

# Routine for partial expectations
# r might be a vector of positive integers
.E_invgauss <-
  function (r = 1, y, mu, a = mu^2) {
    r <- unique(r)
    mu <- cbind(y, a, mu)
    y <- mu[,1]
    a <- mu[,2]
    mu <- mu[,3]
    a <- a/mu
    c_0 <- y/mu
    s_0 <- a * (c_0 - 1)^2 / (2 * c_0)
    c_1 <- (y + 1)/mu
    s_1 <- a * (c_1 - 1)^2 / (2 * c_1)
    posy <- y >= 0
    case1 <- (c_1 <= 1) & posy
    case2 <- (c_0 >= 1) & posy
    case0 <- (!case1 & !case2) & posy
    if (all(r < 2)) {
      E_0_1 <- numeric(length = length(s_0))
      if (any(case0))
        E_0_1[case0] <- (pgamma(s_1[case0] + 2 * a[case0], shape = .5) -
                           pgamma(s_0[case0] + 2 * a[case0], shape = .5)) * exp(2 * a[case0]) +
          pgamma(s_1[case0], shape = .5) + pgamma(s_0[case0], shape = .5)
      if (any(case1))
        E_0_1[case1] <- exp(2 * a[case1]) *
          (pgamma(s_1[case1] + 2 * a[case1], shape = .5) -
             pgamma(s_0[case1] + 2 * a[case1], shape = .5)) -
          (pgamma(s_1[case1], shape = .5) -
             pgamma(s_0[case1], shape = .5))
      if (any(case2))
        E_0_1[case2] <- exp(2 * a[case2]) *
          (pgamma(s_1[case2] + 2 * a[case2], shape = .5) -
             pgamma(s_0[case2] + 2 * a[case2], shape = .5)) +
          (pgamma(s_1[case2], shape = .5) -
             pgamma(s_0[case2], shape = .5))
      E_0_1 <- .5 * E_0_1
      return(E_0_1 * mu)
    }
    if (all(r == 2)) {
      E_0_2 <- numeric(length = length(s_0))
      if (any(case0)) {
        E_0_2[case0] <- (exp(2 * a[case0]) *
                           (pgamma(s_1[case0] + 2 * a[case0], shape = 1.5) -
                              pgamma(s_0[case0] + 2 * a[case0], shape = 1.5)) +
                           (pgamma(s_1[case0], shape = 1.5) +
                              pgamma(s_0[case0], shape = 1.5))) / a[case0] -
          exp(2 * a[case0]) *
          (pgamma(s_1[case0] + 2 * a[case0], shape = .5) -
             pgamma(s_0[case0] + 2 * a[case0], shape = .5)) +
          pgamma(s_1[case0], shape = .5) + pgamma(s_0[case0], shape = .5)
      }
      if (any(case1)) {
        E_0_2[case1] <- (exp(2 * a[case1]) *
                           (pgamma(s_1[case1] + 2 * a[case1], shape = 1.5) -
                              pgamma(s_0[case1] + 2 * a[case1], shape = 1.5)) -
                           (pgamma(s_1[case1], shape = 1.5) -
                              pgamma(s_0[case1], shape = 1.5))) / a[case1] -
          exp(2 * a[case1]) *
          (pgamma(s_1[case1] + 2 * a[case1], shape = .5) -
             pgamma(s_0[case1] + 2 * a[case1], shape = .5)) -
          (pgamma(s_1[case1], shape = .5) -
             pgamma(s_0[case1], shape = .5))
      }
      if (any(case2)) {
        E_0_2[case2] <- (exp(2 * a[case2]) *
                           (pgamma(s_1[case2] + 2 * a[case2], shape = 1.5) -
                              pgamma(s_0[case2] + 2 * a[case2], shape = 1.5)) +
                           (pgamma(s_1[case2], shape = 1.5) -
                              pgamma(s_0[case2], shape = 1.5))) / a[case2] -
          exp(2 * a[case2]) *
          (pgamma(s_1[case2] + 2 * a[case2], shape = .5) -
             pgamma(s_0[case2] + 2 * a[case2], shape = .5)) +
          (pgamma(s_1[case2], shape = .5) -
             pgamma(s_0[case2], shape = .5))
      }
      E_0_2 <- .5 * E_0_2
      return(E_0_2 * mu^2)
    }
    if (all(r < 3)) {
      E_0_1 <- E_0_2 <- numeric(length = length(s_0))
      if (any(case0)) {
        E_0_1[case0] <- (pgamma(s_1[case0] + 2 * a[case0], shape = .5) -
                           pgamma(s_0[case0] + 2 * a[case0], shape = .5)) * exp(2 * a[case0])
        E_0_2[case0] <- (exp(2 * a[case0]) *
                           (pgamma(s_1[case0] + 2 * a[case0], shape = 1.5) -
                              pgamma(s_0[case0] + 2 * a[case0], shape = 1.5)) +
                           (pgamma(s_1[case0], shape = 1.5) +
                              pgamma(s_0[case0], shape = 1.5))) / a[case0] -
          E_0_1[case0] + pgamma(s_1[case0], shape = .5) +
          pgamma(s_0[case0], shape = .5)
        E_0_1[case0] <- E_0_1[case0] + pgamma(s_1[case0], shape = .5) +
          pgamma(s_0[case0], shape = .5)
      }
      if (any(case1)) {
        E_0_1[case1] <- exp(2 * a[case1]) *
          (pgamma(s_1[case1] + 2 * a[case1], shape = .5) -
             pgamma(s_0[case1] + 2 * a[case1], shape = .5))
        E1add <- (pgamma(s_1[case1], shape = .5) -
                    pgamma(s_0[case1], shape = .5))
        E_0_2[case1] <- (exp(2 * a[case1]) *
                           (pgamma(s_1[case1] + 2 * a[case1], shape = 1.5) -
                              pgamma(s_0[case1] + 2 * a[case1], shape = 1.5)) -
                           (pgamma(s_1[case1], shape = 1.5) -
                              pgamma(s_0[case1], shape = 1.5))) / a[case1] -
          E_0_1[case1] - E1add
        E_0_1[case1] <- E_0_1[case1] - E1add
      }
      if (any(case2)) {
        E_0_1[case2] <- exp(2 * a[case2]) *
          (pgamma(s_1[case2] + 2 * a[case2], shape = .5) -
             pgamma(s_0[case2] + 2 * a[case2], shape = .5))
        E1add <- (pgamma(s_1[case2], shape = .5) -
                    pgamma(s_0[case2], shape = .5))
        E_0_2[case2] <- (exp(2 * a[case2]) *
                           (pgamma(s_1[case2] + 2 * a[case2], shape = 1.5) -
                              pgamma(s_0[case2] + 2 * a[case2], shape = 1.5)) +
                           (pgamma(s_1[case2], shape = 1.5) -
                              pgamma(s_0[case2], shape = 1.5))) / a[case2] -
          E_0_1[case2] + E1add
        E_0_1[case2] <- E_0_1[case2] + E1add
      }
      E_0_1 <- .5 * E_0_1
      E_0_2 <- .5 * E_0_2
      if (r[1] == 2)
        return(cbind(`r=2` = E_0_2 * mu^2, `r=1` = E_0_1 * mu))
      return(cbind(`r=1` = E_0_1 * mu, `r=2` = E_0_2 * mu^2))
    }
    stop("not yet implemented for 'r > 2'")
  }

# Routine '.E_invgauss' specialized for first order partial expectation
.E1_invgauss <-
  function (y, mu, a = mu^2) {
    mu <- cbind(y, a, mu)
    y <- mu[,1]
    a <- mu[,2]
    mu <- mu[,3]
    a <- a/mu
    c_0 <- y/mu
    s_0 <- a * (c_0 - 1)^2 / (2 * c_0)
    c_1 <- (y + 1)/mu
    s_1 <- a * (c_1 - 1)^2 / (2 * c_1)
    posy <- y >= 0
    u <- c_1 <= 1
    v <- c_0 > 1
    E_0_1 <- numeric(length = length(s_0))
    if (any(posy))
      E_0_1[posy] <- exp(2 * a[posy]) *
      (pgamma(s_1[posy] + 2 * a[posy], shape = .5) -
         pgamma(s_0[posy] + 2 * a[posy], shape = .5)) +
      ((-1)^u[posy]) * pgamma(s_1[posy], shape = .5) +
      ((-1)^v[posy]) * pgamma(s_0[posy], shape = .5)
    .5 * E_0_1 * mu
  }

# Routine '.E_invgauss' specialized for second order partial expectation
.E2_invgauss <-
  function (y, mu, a = mu^2) {
    mu <- cbind(y, a, mu)
    y <- mu[,1]
    a <- mu[,2]
    mu <- mu[,3]
    a <- a/mu
    c_0 <- y/mu
    s_0 <- a * (c_0 - 1)^2 / (2 * c_0)
    c_1 <- (y + 1)/mu
    s_1 <- a * (c_1 - 1)^2 / (2 * c_1)
    posy <- y >= 0
    u <- c_1 <= 1
    v <- c_0 > 1
    E_0_2 <- numeric(length = length(s_0))
    if (any(posy)) {
      E_0_2[posy] <- exp(2 * a[posy]) *
      (pgamma(s_1[posy] + 2 * a[posy], shape = 1.5) -
         pgamma(s_0[posy] + 2 * a[posy], shape = 1.5)) +
      ((-1)^u[posy]) * pgamma(s_1[posy], shape = 1.5) +
      ((-1)^v[posy]) * pgamma(s_0[posy], shape = 1.5)
      E_0_2[posy] <- E_0_2[posy] / a[posy] -
      exp(2 * a[posy]) *
        (pgamma(s_1[posy] + 2 * a[posy], shape = .5) -
           pgamma(s_0[posy] + 2 * a[posy], shape = .5)) -
        ((-1)^(1 - u[posy])) * pgamma(s_1[posy], shape = .5) -
        ((-1)^(1 - v[posy])) * pgamma(s_0[posy], shape = .5)
    }
    .5 * E_0_2 * mu^2
  }

# Routine '.E_invgauss' specialized for first and second order partial expectations
.E12_invgauss <-
  function (y, mu, a = mu^2) {
    mu <- cbind(y, a, mu)
    y <- mu[,1]
    a <- mu[,2]
    mu <- mu[,3]
    a <- a/mu
    c_0 <- y/mu
    s_0 <- a * (c_0 - 1)^2 / (2 * c_0)
    c_1 <- (y + 1)/mu
    s_1 <- a * (c_1 - 1)^2 / (2 * c_1)
    posy <- y >= 0
    u <- c_1 <= 1
    v <- c_0 > 1
    E_0_1 <- E_0_2 <- numeric(length = length(s_0))
    if (any(posy)) {
      int1 <- exp(2 * a[posy]) *
        (pgamma(s_1[posy] + 2 * a[posy], shape = .5) -
           pgamma(s_0[posy] + 2 * a[posy], shape = .5))
      int2 <- pgamma(s_1[posy], shape = .5)
      int3 <- pgamma(s_0[posy], shape = .5)
      E_0_1[posy] <- int1 + ((-1)^u[posy]) * int2 + ((-1)^v[posy]) * int3
      E_0_2[posy] <- exp(2 * a[posy]) *
        (pgamma(s_1[posy] + 2 * a[posy], shape = 1.5) -
           pgamma(s_0[posy] + 2 * a[posy], shape = 1.5)) +
        ((-1)^u[posy]) * pgamma(s_1[posy], shape = 1.5) +
        ((-1)^v[posy]) * pgamma(s_0[posy], shape = 1.5)
      E_0_2[posy] <- E_0_2[posy] / a[posy] - int1 -
        ((-1)^(1 - u[posy])) * int2 -
        ((-1)^(1 - v[posy])) * int3
    }
    cbind(`r=1` = .5 * E_0_1 * mu, `r=2` = .5 * E_0_2 * mu^2)
  }

# Routine for difference of successive
# first order partial expectations:
# .E_invgauss (r = 1, y, mu, a = mu^2) - .E_invgauss (r = 1, y - 1, mu, a = mu^2)
.diffE1_invgauss <-
  function (y, mu, a = mu^2) {
    mu <- cbind(y, a, mu)
    y <- mu[,1]
    a <- mu[,2]
    mu <- mu[,3]
    a <- a/mu
    cneg <- (y-1)/mu
    sneg <- a * (cneg - 1)^2 / (2 * cneg)
    c_0 <- y/mu
    s_0 <- a * (c_0 - 1)^2 / (2 * c_0)
    c_1 <- (y + 1)/mu
    s_1 <- a * (c_1 - 1)^2 / (2 * c_1)
    posy <- y > 0
    zeroy <- y == 0
    u0 <- c_0 <= 1
    u1 <- c_1 <= 1
    v_1 <- cneg > 1
    diffE <- numeric(length = length(s_0))
    if (any(posy))
      diffE[posy] <- exp(2 * a[posy]) *
      (pgamma(s_1[posy] + 2 * a[posy], shape = .5) +
         pgamma(sneg[posy] + 2 * a[posy], shape = .5) -
         2 * pgamma(s_0[posy] + 2 * a[posy], shape = .5)) +
      ((-1)^u1[posy]) * pgamma(s_1[posy], shape = .5) -
      2 * ((-1)^u0[posy]) * pgamma(s_0[posy], shape = .5) -
      ((-1)^v_1[posy]) * pgamma(sneg[posy], shape = .5)
    if (any(zeroy))
      diffE[zeroy] <- exp(2 * a[zeroy]) *
      (pgamma(s_1[zeroy] + 2 * a[zeroy], shape = .5) -
         pgamma(s_0[zeroy] + 2 * a[zeroy], shape = .5)) +
      ((-1)^u1[zeroy]) * pgamma(s_1[zeroy], shape = .5) +
      ((-1)^(1 - u0[zeroy])) * pgamma(s_0[zeroy], shape = .5)
    .5 * diffE * mu
  }

##################################
# Routine for .invgzeta function
.invgzeta_comps <- function(z, mu, a = mu^2) {
  E12 <- .E12_invgauss (z, mu = mu, a = a)
  comp1 <- sum((2 * z + 1) * E12[,1])
  comp2 <- sum(E12[,2])
  comp3 <- sum(z * (z + 1) * (pinvgauss(z + 1, mu = mu, a = a) -
                                pinvgauss(z, mu = mu, a = a)))
  c(comp1, comp2, comp3)
}

# Routine for .vardinvgauss function
.invgzeta <- function(mu, a = mu^2, tol = 1e-15) {
  z_i <- floor(qinvgauss(tol/2, mu = mu, a = a))
  if (is.na(z_i))
    z_i <- max(0, floor(mu - 5 * sqrt(mu^3 / a)))
  z_f <- floor(qinvgauss(1 - tol/2, mu = mu, a = a)) + 1
  if (is.na(z_f))
    z_f <- max(2, floor(mu + 5 * sqrt(mu^3 / a)))
  zcomps <- .invgzeta_comps(z_i:z_f, mu = mu, a = a)
  res <- zcomps[1] - zcomps[2] - zcomps[3]
  structure(res, precision = pinvgauss(z_i, mu = mu, a = a) +
              pinvgauss(z_f + 1, mu = mu, a = a, lower.tail = FALSE))
}

# variance function
.vardinvgauss <- function (mu, a = 1, tol = 1e-15) {
  na <- max(length(a), length(mu))
  if (na == 1)
    return(as.vector((mu^3) / a + .invgzeta(mu = mu, a = a, tol = tol)))
  AB <- cbind(mu, a)
  mu <- AB[,1]
  a <- AB[,2]
  IGvalues <- apply(AB, MARGIN = 1, FUN = function(ab) {
    .invgzeta(mu = ab[1], a = ab[2], tol = tol)
  })
  as.vector((mu^3) / a + IGvalues)
}

# moment estimates (approximate for a)
.momente_bdinvgauss <- function(x, zeta = NULL) {
  xbar <- mean(x)
  sigma2 <- var(x)
  if (is.null(zeta)) {
    minz <- min(xbar/2, .125)
    zeta <- if (sigma2 > minz) minz/2 else sigma2/4
  }
  a <- xbar^3 / (sigma2 - zeta)
  c(mu = xbar, a = a, sigma2 = sigma2)
}

# a for equi-dispersion
.equia_bdinvgauss <- function (mu, tol = 1e-10) {
  tol <- min(1e-5, max(1e-15, tol))
  lfun <- function(mui) {
    if (mui == 0)
      return(1)
    fun <- function(x) mui^3 / exp(x) + .invgzeta(mu = mui, a = exp(x), tol = tol) - mui
    res <- tryCatchWE(uniroot(f = fun,
                   lower = 0, upper = logb(5), extendInt = "yes", tol = tol))$value
    if (any(class(res) %in% c("error", "try-error", "condition")))
      return(NaN)
    structure(exp(as.vector(res$root)), diffEV = as.vector(res$f.root))
  }
  if (length(mu) == 1)
    lfun(mu)
  else
    sapply(mu, FUN = lfun)
}

# Fit a discrete gamma distribution to an IDD sample x
.mle_bdinvgauss <- function(x, theta = NULL, alpha = .05, asymp = "norm",
                         method = "BFGS", tol = 1e-10, maxit = 100) {
  if (missing(theta) || is.null(theta))
    start <- logb(.momente_bdinvgauss(x)[1:2])
  else {
    start <- logb(theta[1:2])
    names(start) <- c("mu", "a")
  }
  # Negative log likelihood function: theta = c(log(mu), log(a))
  nlike.dinvgauss <- function(theta, x) {
    -sum(ddinvgauss(x, mu = exp(theta[1]), a = exp(theta[2]), log = TRUE))
  }
  OPTIM <- optim(par = start, fn = nlike.dinvgauss, x = x,
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
