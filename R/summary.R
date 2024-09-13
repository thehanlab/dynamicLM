#' Summarize dynamic_lm objects: not yet implemented
#'
#' @param object dynamicLM object
#' @param ... Additional arguments
#'
#' @export
summary.dynamicLM <- function(object, ...) {
  stop(tidymess("A summary function is not yet implemented for dynamicLM, please
                use print() or plot() or str()."))
}


#' Summarize pen_lm objects: not yet implemented
#'
#' @param object pen_lm object
#' @param ... Additional arguments
#'
#' @export
summary.pen_lm <- function(object, ...) {
  stop(tidymess("A summary function is not yet implemented for dynamicLM, please
                use print() or plot() or str()."))
}


#' Summarize cv.pen_lm objects: not yet implemented
#'
#' @param object cv.pen_lm object
#' @param ... Additional arguments
#'
#' @export
summary.cv.pen_lm <- function(object, ...) {
  stop(tidymess("A summary function is not yet implemented for dynamicLM, please
                use print() or plot() or str()."))
}
