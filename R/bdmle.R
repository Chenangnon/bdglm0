#' @name bdmle
#' @aliases bdmomente bdvariance bdequiscale bdequishape
#' @title Dispersion and Estimation in Balanced Discrete Distribution Families
#' @description Compute variance, equiscale/bdequishape (scale/shape parameter for equidispersion), approximate moment estimates, and maximum likelihood estimates for parameters of balanced discrete distributions.
#'
#' @usage
#' bdvariance (mu, scale = 1, shape, family = "bdgamma", tol = 1e-10)
#'
#' bdequiscale (mu, shape, family = "bdgamma", tol = 1e-10)
#'
#' bdequishape (mu, scale, family = "bdbeta", tol = 1e-10)
#'
#' bdmomente (x, family = "bdgamma", zeta = NULL)
#'
#' bdmle (mu, family = "bdgamma", tol = 1e-10, ...)
#'
#' @param mu vector of positive real mean parameters of the balanced discrete distribution.
#' @param scale vector of positive real scale parameters of the balanced discrete distribution.
#' @param shape vector of positive real shape parameters of the balanced discrete distribution.
#' @param family character identifying a family of balance discrete distribution.
#' @param tol a tolerance value in the open (0, 1), controls the precision of the variance given by an infinite series which is truncated in the computations.
#' @param x vector of integer quantiles.
#' @param zeta positive scalar in the open (0, 1/4) and representing the \code{zeta} quantity in Tovissode et al. (2021). Defauls to \code{1/8} if \code{sigma2 > 1/8} and \code{sigma2/4} otherwise with sigma2 the sample variance.
#' @param ... additional arguments to be passed to the low level distribution fitting functions (see below).

# @param theta vector of initial values of parameters mu and a
# @param alpha significance level for confidence intervals of parameter estimates.
# @param asymp character, indicates the asymptotic distribution of parameter estimates for building confidence intervals. Defaults to "norm" i.e. normal distribution at log scale. Any other character will use a Student t distribution.
# @param method The method to be used. See \code{Details} of \link{optim}. Can be abbreviated. Defaults to "BFGS".
# @param maxit The maximum number of iterations. Defaults to 100.

# @param tol a tolerance value in the open (0, 1), controls the precision of the variance approximation.
# @param zeta positive scalar in the open (0, 1/4) and representing the zeta quantity in Tovissode et al. (2020). Defauls to 1/8 if sigma2 > 1/8 and sigma2/4 otherwise with sigma2 the sample variance.
# @param theta vector of initial values of parameters mu, b and c
# @param alpha significance level for confidence intervals of parameter estimates.
# @param asymp character, indicates the asymptotic distribution of parameter estimates for building confidence intervals. Defaults to "norm" i.e. normal distribution at log scale. Any other character will use a Student t distribution.
# @param method The method to be used. See 'Details' of \link{optim}. Can be abbreviated. Defaults to "BFGS".
# @param maxit The maximum number of iterations. Defaults to 100.
# @param auxi.par logical, should the shape "c" be considered as an "auxilliary" parameter and calculated by inverting the equidispersion condition? Defaults to TRUE, i.e. mu and b are specified and c is computed as equic.bdggamma(mu, b). The result seems unstable (large asymptotic variance) but stabilizes when the sample size becomes very large (5000-10000). The alternative "b" seems much more unstable, with convergence problems.


#' @details The functions are based on general formulae given in Tovissode et al. (2021).
#'
#' @return \code{bdvariance} computes the variance and \code{bdequiscale} determines for a fixed \code{mu} the scale parameter corresponding to a marginally equidispersed balanced discrete distribution.
#' \code{bdequishape} does the same for the shape, for a fixed \code{mu} and \code{scale} parameter.
#'
#' \code{bdmomente} returns moments estimates for the distribution parameters \code{(mu, a)}, and \code{shape} if applicable.
#' \code{bdmle} computes maximum likelihood estimates.
#'
#' @seealso \link{bdfamilies} for generating balanced discrete family objects for model fitting purposes.
#'
#' @references Tovissode, C.F.; Honfo, S.H.; Doumate, J.T.; Glele Kakai, R. On the Discretization of Continuous Probability Distributions Using a Probabilistic Rounding Mechanism. Mathematics 2021, 9, 555. https://doi.org/10.3390/math9050555

bdvariance <-
  function (mu, scale = 1, shape, family = "bdgamma", tol = 1e-10) {
    switch(family,
           bdgamma = .vardgamma (mu = mu, a = scale, tol = tol),
           bdlnorm = .vardlnorm (mu = mu, alpha = scale, tol = tol),
           bdinvgauss = .vardinvgauss (mu = mu, a = scale, tol = tol),
           bdweibull = .vardweibull (mu = mu, c = scale, tol = tol),
           bdinvgauss = .vardinvgauss (mu = mu, a = scale, tol = tol),
           bdggamma = .vardggamma (mu = mu, b = scale, c = shape, tol = tol),
           bdbeta = .vardbeta (mu = mu, a = scale, shape = shape, tol = tol))
  }

bdequiscale <-
  function (mu, shape, family = "bdgamma", tol = 1e-10) {
    switch(family,
           bdgamma = .equia_bdgamma (mu = mu, tol = tol),
           bdlnorm = .equialpha_bdlnorm (mu = mu, tol = tol),
           bdinvgauss = .equia_bdinvgauss (mu = mu, tol = tol),
           bdweibull = .equic_bdweibull (mu = mu, tol = tol),
           bdinvgauss = .equia_bdinvgauss (mu = mu, tol = tol),
           bdggamma = .equib_bdggamma (mu = mu, c = shape, tol = tol),
           bdbeta = .equia_bdbeta (mu = mu, shape = shape, tol = tol))
  }

bdequishape <-
  function (mu, scale, family = "bdbeta", tol = 1e-10) {
    switch(family,
           bdbeta = .equishape_bdbeta (mu = mu, a = scale, tol = tol),
           bdggamma = .equic_bdggamma (mu = mu, b = scale, tol = tol),
           stop("not recognized three-parameter balanced discrete family"))
  }

bdmomente <-
  function (x, family = "bdgamma", zeta = NULL) {
    switch(family,
           bdgamma = .momente_bdgamma (x = x, zeta = zeta),
           bdlnorm = .momente_bdlnorm (x = x, zeta = zeta),
           bdinvgauss = .momente_bdinvgauss (x = x, zeta = zeta),
           bdweibull = .momente_bdweibull (x = x, zeta = zeta),
           bdggamma = .momente_bdggamma (x = x, zeta = zeta),
           bdinvgauss = .momente_bdinvgauss (x = x, zeta = zeta),
           bdbeta = .momente_bdbeta (x = x, zeta = zeta))
  }

bdmle <-
  function (x, family = "bdgamma", tol = 1e-10, ...) {
    switch(family,
           bdgamma = .mle_bdgamma (x = x, tol = tol),
           bdlnorm = .mle_bdlnorm (x = x, tol = tol),
           bdinvgauss = .mle_bdinvgauss (x = x, tol = tol),
           bdweibull = .mle_bdweibull (x = x, tol = tol),
           bdggamma = .mle_bdggamma (x = x, tol = tol, ...),
           bdbeta = .mle_bdbeta (x = x, tol = tol, ...))
  }
