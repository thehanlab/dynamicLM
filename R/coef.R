#' Get the coefficients of a fitted supermodel in dynamicLM
#'
#' @param object Fitted supermodel
#' @param ... Additional arguments to `coef`
#'
#' @return Vector of coefficients for a Cox landmark supermodel or list of
#'   coefficients for each cause-specific model for a CSC landmark supermodel.
#' @export
coef.dynamicLM <- function(object, ...) {
  return(coef(object$model, ...))
}
