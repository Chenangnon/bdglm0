#' @name BD-Families
#' @aliases bd.distr.families
#' @title Balanced Discrete Families of Distributions
#' @description Display the list of balanced discrete distribution families available in the package \code{bdglm}.
#' Or check if a continuous distribution has a balanced discrete counterpart implemented in the package.
#'
#' @param family character, (abbreviated) name of a continuous distribution,
#' e.g. 'gamma' for the gamma distribution, 'beta' for the beta distribution.
#'
#'  If \code{family} is not specified, the (abbreviated) names of the available balanced discrete distribution families are \code{NULL} listed and is returned.
#'
#' @usage bd.distr.families(family)
#'
#' @details For a specific family named \code{xxx}, the functions for the probability mass function, cumulative distribution function, quantile function and random variate generation are named in the form ddxxx, pdxxx, qdxxx and rdxxx respectively.
#'
# For the balanced discrete exponential distribution see \link{ddexp}
#
# For the balanced discrete chi-squared distribution see \link{ddchisq}
#
# For the balanced discrete chi distribution see \link{ddchi}
#'
#' For the balanced discrete the Weibull distribution see \link{ddweibull}
#'
#' For the balanced discrete gamma distribution see \link{ddgamma}
#'
#' For the balanced discrete generalized gamma distribution see \link{ddggamma}
#'
#' For the balanced discrete log-normal distribution see \link{ddlnorm}
#'
#' For the balanced discrete normal distribution see \link{ddnorm}
#'
#' For the balanced discrete inverse Gaussian distribution see \link{ddinvgauss}
#'
#' For the balanced discrete beta distribution see \link{ddbeta}
#'
# For the balanced discrete Student's t distribution see \link{ddt}
#'
# For the balanced discrete Fisher's F distribution see \link{ddf}
#
# For the balanced discrete Azalini's skew normal distribution see \link{ddsnorm}
#'
# For the balanced discrete uniform distribution see \link{ddunif}
#
#' @return
#' A logical vector of the same size as the argument \code{family} when the later is specified;
#' otherwise \code{NULL} is returned after displaying a list of the available families with the corresponding (abbreviated) names.
#'
# @seealso See \link{buidbdfamily} for constructing family object/balanced discrete distribution functions from an arbitrary continuous distribution.
# bdfamilies

bd.distr.families <-
  function (family, help = FALSE) {
    if (!missing(family)) {
      if(is.character(family)) {
      OKfamilies <- c('gamma', 'ggamma', 'weibull', 'lnorm',
                      'normal', 'invgauss', 'beta')
      return(family %in% OKfamilies)
      }
      else {
        stop ("argument 'family' must be a character")
      }
    }
    BDFamilies <- c(gamma = "'gamma' for the gamma distribution",
                    #                  exp = "'exp' for the exponential distribution",
                    #                  chisq = "'chisq' for the chi-squared distribution",
                    #                  chi = "'chi' for the chi distribution",
                    ggamma = "'ggamma' for the generalized gamma distribution",
                    weibull = "'weibull' for the Weibull distribution",
                    lnorm = "'lnorm' for the log-normal distribution",
                    normal = "'norm' for the normal distribution",
                    #                  snormal = "'snormal' for the Azalini skew normal distribution",
                    #                  t = "'t' for the Student's t distribution",
                    #                  f = "'f' for the F distribution",
                    #                  unif = "'unif' for the uniform distribution",
                    invgauss = "'invgauss' for the inverse Gaussian distribution",
                    beta = "'beta' for the beta distribution")
                    cat("  Type ?ddxxx for the balanced discrete xxx distribution.\n")
                    cat("  Possible xxx values are:\n")
    cat(paste0("\n      ", BDFamilies))
  }

