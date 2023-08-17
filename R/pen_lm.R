#' Compute the regularization path of coefficients for a Cox or
#' cause-specific Cox landmark supermodel with lasso or elasticnet penalization.
#'
#' Fit by calling `[glmnet()]`. As in `glmnet`, the model is fit via penalized
#' maximum likelihood to produce a regularization path at a grid of values for
#' the regularization parameter lambda. Input can be as typically done for
#' `glmnet` in the form of `x` and `y` which are a matrix and response object
#'  or with a landmark super dataset specifying the dependent columns in `y`.
#'
#' @param x Either an object of class "LMdataframe" or a matrix. An LMdataframe
#'   can be created by running [stack_data()] and [add_interactions()]. A
#'   matrix should have each row be an observation.
#' @param y If `x` is an LMdataframe, `y` is optional and should be vector of
#'   column names of the data stored in `x` that are to be used as
#'   dependent variables. If not specified, it is assumed that all non-response
#'   variables are the dependent variables.
#'
#'   If `x` is a matrix, `y` should be the response variable: either a
#'   [survival::Surv()] or [prodlim::Hist()] object.
#' @param alpha The elastic net mixing parameter: Lies between 0 and 1. At 1,
#'   the penalty is the LASSO penalty, and at 0, the penalty is the ridge
#'   penalty. The default is 1.
#' @param ... Additional arguments passed to [glmnet()].
#'
#' @return An object of class `pen_lm`. This is a list of `glmnet` objects (one
#'   for each cause-specific Cox model or a list of length one for a regular Cox
#'   model). The object also has attributes `survival.type` (`competing.risk`
#'   or `survival`) and `lmdata` and `xcols` which store the inputs if given.
#'   Functions `print` and `plot` exist for the object. To make predictions,
#'   see [dynamic_lm.pen_lm()] and [predict()].
#' @seealso [print.pen_lm()], [plot.pen_lm()], [dynamic_lm.pen_lm()]
#' @import glmnet
#' @export
#'
#' @examples
#' \dontrun{
#' }
pen_lm <- function(x, y, alpha = 1, ...) {

  checked_input <- match.call()
  checked_input$parent_func <- quote(pen_lm)
  checked_input$CV <- FALSE
  checked_input[[1L]] <- NULL #quote(check_penlm_inputs) #
  checked_input <- do.call("check_penlm_inputs", as.list(checked_input), envir = parent.frame()) #eval(checked_input, parent.frame()) #

  # can call pen_lm again depending on inputs
  if (inherits(checked_input, "pen_lm")) return(checked_input)

  # if not, use checked inputs
  x = checked_input$x
  y = checked_input$y
  lmdata = checked_input$lmdata
  xcols = checked_input$xcols
  alpha = checked_input$alpha

  models <- lapply(y, function(yi) {
    glmnet::glmnet(x = x, y = yi, family = "cox", alpha = alpha, ...)
  })
  if (length(models) > 1){
    attr(models, "survival.type") <- "competing.risk"
  } else {
    attr(models, "survival.type") <- "survival"
  }
  if (!is.null(lmdata)) attr(models, "lmdata") <- lmdata
  if (!is.null(xcols)) attr(models, "xcols") <- xcols
  attr(models, "alpha") <- alpha
  class(models) <- "pen_lm"
  return(models)
}

