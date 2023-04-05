#' Fit a regular or cause-specific Cox landmark supermodel with lasso or
#' elastic-net penalization.
#'
#' Fit by calling `[glmnet()]`. As in `glmnet`, the model is fit via penalized
#' maximum likelihood to produce a regularization path at a grid of values for
#' the regularization parameter lambda. Input can be as typically done for
#' `glmnet` in the form of `x` and `y` or with a landmark super dataset `lmdata`
#' specifying dependent columns in `xcols`.
#'
#' @param x Input matrix with each row an observation.
#' @param y Response variable: either a `Surv` or `Hist` object.
#' @param lmdata An object of class "LMdataframe", this can be created by
#'   running [stack_data()] and [add_interactions()]
#' @param xcols A vector of column names of the data stored in `lmdata` that
#'   are to be used as dependent variables. If not specified, it is assumed that
#'   all non-response variables are the dependent variables.
#' @param alpha The elastic net mixing parameter: Lies betweent 0 and 1. At 1,
#'   the penalty is the LASSO penalty, and at 0, the penalty is the ridge
#'   penalty.
#' @param ... Additional arguments passed to [glmnet()].
#'
#' @return An object of class `penLM`. This is a list of `glmnet` objects (one
#'   for each cause-specific Cox model or a list of length one for a regular Cox
#'   model). The object also has attributes `survival.type` (`competing.risk`
#'   or `survival`) and `lmdata` and `xcols` which store the inputs if given.
#'   Functions `print` and `plot` exist for the object. To make predictions,
#'   see [dynamic_lm.penLM()] and [predict()].
#' @import glmnet
#' @export
#'
#' @examples
#' \dontrun{
#' }
penLM <- function(x, y, lmdata, xcols, alpha = 1, ...) {

  checked_input <- match.call()
  checked_input$parent_func <- quote(penLM)
  checked_input$CV <- FALSE
  checked_input[[1L]] <- quote(check_penLM_inputs)
  checked_input <- eval(checked_input, parent.frame())

  # can call penLM again depending on inputs
  if (class(checked_input) == "penLM") return(checked_input)

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
  class(models) <- "penLM"
  return(models)
}

