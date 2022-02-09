#' @name asOneFormula
#' @aliases allVarsRec formula.bdglm
#' @title Formula manipulation
#' @description Operations on formula object
#' @param omit characters to be omitted
#' @param ... formula objects
#' @param object an R object
#'
#'  a \code{bdglm} class object for \code{formula.bdglm}.
#'
#'  a \code{formula} class object or a list of \code{formula} class objects for \code{allVarsRec}.
#'
#' @usage asOneFormula (..., omit = c(".", "pi"))
#'
#' allVarsRec (object)
#'
#' @details \code{asOneFormula} merges the right hand side of two or more \code{formula} objects.
#'
#' \code{allVarsRec} finds all variables in a \code{formula} or a list of \code{formula} objects.
#'
#' \code{formula.bdglm} extracts the \code{formula} component of a \code{bdglm} object.
#'
#' @return \code{asOneFormula} returns one \code{formula} object.
#'
#' \code{allVarsRec} return a character vector with the extracted variable names.
#'
#' \code{formula.bdglm} produces an object of class \code{formula} which contains a symbolic model formula.

asOneFormula <- function (..., omit = c(".", "pi")) {
  names <- unique(allVarsRec(list(...)))
  names <- names[is.na(match(names, omit))]
  if (length(names))
    eval(parse(text = paste("~", paste(names, collapse = "+")))[[1]])
  else ~1
}

allVarsRec <- function (object) {
  if (is.list(object)) {
    unlist(lapply(object, allVarsRec))
  }
  else {
    all.vars(object)
  }
}

formula.bdglm <- function (object,
                           mode = c("mean", "mu", "scale", "a",
                                    "shape", "b", "shape2", "c",
                                    "full"), ...) {
  mode <- match.arg(mode)
  if(inherits(object, "bdgglm"))
    switch(mode[1],
           mu = object$formula,
           mean = object$formula,
           scale =  object$formula.scale,
           a =  object$formula.scale,
           full = object$formula.full)
  else
    switch(mode[1],
           mu = object$formula,
           mean = object$formula,
           scale =  object$formula.scale,
           a =  object$formula.scale,
           shape =  object$formula.shape,
           b =  object$formula.shape,
           shape2 =  object$formula.shape2,
           c =  object$formula.shape2,
           full = object$formula.full)
}
