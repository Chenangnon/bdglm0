#' @name pgammadot
#'
#' @title Derivatives of the Incomplete Gamma Integral
#'
#' @description Compute the first two derivatives of the incomplete gamma integral
#'
#' @usage
#' pgammadot (q, shape, tmax = 100, tol = 1.e-6, eps = 1.e-30)
#'
#' @param q,shape numeric, vectors of finite and non-negative values.
#' @param tmax integer, maximum number of iterations in the computations.
#' @param tol,eps control parameters (convergence tolerance and ~zero).
#'
#' @details The Incomplete Gamma Integral (ratio) is considered as given by \code{pgamma(q, shape)}.
#'
#' The computations are based on Moore (1982), i.e. the function uses a series expansion for \code{shape <= q <= 1} or \code{q} < \code{shape}, else otherwise a continued fraction expansion. Machine overflow can occur for large values of \code{q} when \code{q >> shape}.
#'
#' @return A matrix, each line corresponding to the corresponding element of \code{cbind(q, shape)} (element of \code{q} or \code{shape}).
#'
#' The matrix has six columns. The first two contain the gradient of \code{pgamma(q, shape)} (first derivatives with respect to \code{q} and \code{shape}).
#'
#' The columns three to five contain the unique elements of the hessian of \code{pgamma(q, shape)} (second derivatives with respect to \code{q^2}, \code{q * shape}, \code{shape^2}).
#' The column six contain the approximate value of \code{pgamma(q, shape)}.
#'
#' @references Moore, R. J. (1982). Algorithm AS 187: Derivatives of the Incomplete Gamma Integral. Journal of the Royal Statistical Society, Series C (Applied Statistics), 31(3), 330-335.
#'
#' @seealso \link{pgamma}.
#'

pgammadot <-
  function (q, shape, tmax = 100, tol = 1.e-6, eps = 1.e-30) {
    nmq <- names(q)
    shape <- cbind(q, shape)
    q <- shape[,1]
    shape <- shape[,2]
    validq <- is.numeric(q) & is.finite(q) & (q > 0)
    #stop("argument 'q' must be numeric, finite, and positive")
#    if (all(!validq))
#      return(replicate(6, validq))
    validshape <- is.numeric(shape) & is.finite(shape) & (shape > 0)
    valIN <- validq + validshape == 2
    if (!is.numeric(tmax) | !is.finite(tmax) | (tmax <= 0))
      stop("'tmax' must be numeric, a finite positive integer between 20 and 300")
    tmax <- min(max(20, tmax[1]), 300)
    gplog <- lgamma(shape)
    gp1log <- gplog + log(shape)
    psip <- digamma(shape)
    psip1 <- psip + 1/shape
    psidp <- trigamma(shape)
    psidp1 <- psidp - 1/shape^2
    out <- matrix(0, nrow = length(q), ncol = 7)
    if (any(valIN))
    out[valIN,] <- t(apply(cbind(q[valIN], shape[valIN], gplog[valIN],
                                 gp1log[valIN], psip[valIN], psip1[valIN],
                                 psidp[valIN], psidp1[valIN]),
                 MARGIN = 1, FUN = function(x) {
                   resx <- .digami_R(X = x[1], P = x[2],
                                    GPLOG = x[3], GP1LOG = x[4],
                                    PSIP = x[5], PSIP1 = x[6],
                                    PSIDP = x[7], PSIDP1 = x[8],
                                    E = tol, OFLO = 1.E300,
                                    TMAX = tmax, VSMALL = eps)
                   c(resx, attr(resx, "fail"))
                 }))
    res <- out[,c(1, 3, 2, 5, 4, 6), drop = FALSE]
    if (!is.null(nmq))
      nmq <- rep(nmq, length.out = length(q))
    dimnames(res) <- list(nmq,
                          c("q", "shape", "q.2", "q.shape",
                            "shape.2", "pgamma"))
    if (any(out[,7] != 0)) {
      indices <- which(out[,7] != 0)
      warning("convergence problems with elements ",
              indices)
    }
    res
  }
