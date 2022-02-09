# Review the definition of bdframe.bdbglm, bdframe.bdlnglm, bdframe.bdwglm

#' @name bdframe
#'
#' @title Balanced Discrete Model frame
#'
#' @description Build or Extract the Balanced Discrete Model Frame from a Formula or a Fit. The function and its methods return a \code{list} of the design (response and covariates) and starting values needed to use \code{formula} and any ... arguments.
#'
#' @usage bdframe (formula, family = bdgamma (), ...)
#'
#' @param formula a model \link{formula} or an \code{R} object of a class with a \code{bdframe} method.
#' @param family family of balanced discrete distribution.
#' @param ... further argument passed to or from other methods.
#'
#' @details The default method is \code{bdframe.default} designed for \code{bdglm} class objects.
#'
#' If \code{formula} is a just a \code{formula} object and \code{family} is none of the build-in balanced discrete families,
#' the balanced discrete gamma family if assumed.
#'
#' For completeness, it calls the \link{model.frame} when the argument \code{formula} is not of any bdglm sub-class (e.g. \code{bdgglm}, \code{bdbglm}, ...).
#'
#

bdframe.default <- function(formula, family = bdgamma (), ...) {
  if (inherits(formula, "bdbglm") || identical(family$family, "bdbeta"))
    .bdframe.bdbglm (formula, ...)
  else if (inherits(formula, "bdlnglm") || identical(family$family, "bdlnorm"))
    .bdframe.bdlnglm (formula, ...)
  else if (inherits(formula, "bdwglm") || identical(family$family, "bdweibull"))
    .bdframe.bdwglm (formula, ...)
  else if (any(inherits(formula, "bdgglm"),
               inherits(formula, "formula"),
               identical(family$family, "bdgamma")))
    .bdframe.bdgglm (formula, ...)
  else {
    model.frame(formula(formula), ...)
  }
}

setGeneric(name = "bdframe",
           def = bdframe.default)

.bdframe.bdgglm <- function (formula, formula.scale = ~ 1, data, weights,
                            offset, subset, na.action, start = list(),
                            contrasts = NULL, contrasts.scale = NULL, ...) {
  if (inherits(formula, "bdgglm")) {
    if (!is.null(formula$bdframe))
      return(formula$bdframe)
  }
  bdgglm (formula = formula, formula.scale = formula.scale, data = data,
          weights = weights, offset = offset, subset = subset, na.action = na.action,
          start = start, method = "bdframe", contrasts = contrasts, contrasts.scale = contrasts.scale, ...)
}

.bdframe.bdbglm <- .bdframe.bdlnglm <- .bdframe.bdwglm <- .bdframe.bdgglm
