#' Compute the regularization path of coefficients for a Cox or
#' cause-specific Cox landmark supermodel with lasso or elasticnet penalization.
#'
#' Fit by calling `[glmnet::glmnet()]`. As in `glmnet`, the model is fit via penalized
#' maximum likelihood to produce a regularization path at a grid of values for
#' the regularization parameter lambda. Input can be as typically done for
#' `glmnet` in the form of `x` and `y` which are a matrix and response object
#'  or with a landmark super dataset specifying the dependent columns in `y`.
#'
#' @param x An "LMdataframe", which can be created by running [stack_data()] and
#'   [add_interactions()].
#' @param y Optional, a vector of column names of the data stored in `lmdata`
#'   that are to be used as dependent variables. If not specified, it is assumed
#'   that all non-response variables are the dependent variables.
#' @param alpha The elastic net mixing parameter: Lies between 0 and 1. At 1,
#'   the penalty is the LASSO penalty, and at 0, the penalty is the ridge
#'   penalty. The default is 1.
#' @param ... Additional arguments passed to `glmnet()`.
#'
#' @return An object of class `pen_lm`. This is a list of `glmnet` objects (one
#'   for each cause-specific Cox model or a list of length one for a regular Cox
#'   model). The object also has attributes `survival.type` (`competing.risk`
#'   or `survival`) and `lmdata` and `xcols` which store the inputs if given.
#'   Functions `print` and `plot` exist for the object. To make predictions,
#'   see [dynamic_lm.pen_lm()] and [predict.dynamicLM()].
#' @seealso [print.pen_lm()], [plot.pen_lm()], [dynamic_lm.pen_lm()]
#' @examples
#' \dontrun{
#' }
#' @import glmnet
#' @export
pen_lm <- function(x, y, alpha = 1, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package \"glmnet\" must be installed to use function pen_lm.",
         call. = FALSE)
  }

  args <- match.call()    # store call
  checked_input <- args   # to be cleaned

  # Save the arguments in 'evaluated' form, so we don't have to find them later
  args$x <- NULL          # remove lmdata (heavy)
  evaluated_args <- as.list(args)
  function_name <- evaluated_args[[1]]
  evaluated_args <- lapply(evaluated_args[-1], eval)
  evaluated_args <- c(list(function_name), evaluated_args)
  args <- as.call(evaluated_args)

  # Clean input
  m <- match(c("x", "y", "alpha"), names(checked_input), 0L)
  checked_input <- as.list(checked_input[m])
  checked_input <- do.call(check_penlm_inputs, checked_input,
                           envir = parent.frame())
  x <- checked_input$x
  y <- checked_input$y
  lmdata <- checked_input$lmdata
  xcols <- checked_input$xcols
  alpha <- checked_input$alpha

  # Get the coefficient path
  # TODO: add parallelization
  models <- lapply(y, function(yi) {
    glmnet::glmnet(x = x, y = yi, family = "cox", alpha = alpha, ...)
  })

  # Additional information for output
  if (length(models) > 1){
    attr(models, "survival.type") <- "competing.risk"
  } else {
    attr(models, "survival.type") <- "survival"
  }
  if (!is.null(lmdata)) attr(models, "lmdata") <- lmdata
  if (!is.null(xcols)) attr(models, "xcols") <- xcols
  attr(models, "alpha") <- alpha
  attr(models, "args") <- args
  class(models) <- "pen_lm"

  return(models)
}

