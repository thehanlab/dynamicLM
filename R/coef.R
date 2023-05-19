#' Get the coefficients of a fitted supermodel in dynamicLM
#'
#' @param object
#'
#' @return Vector of coefficients for a Cox landmark supermodel or list of
#'   coefficients for each cause-specific model for a CSC landmark supermodel.
#' @export
coef.dynamicLM <- function(object, ...) {
  return(stats::coef(object$model, ...))
}
