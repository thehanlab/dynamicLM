#' Cross-validation for a penalized Cox or cause-specific Cox landmark
#' supermodel
#'
#' Fit by calling [glmnet::cv.glmnet()]. As in `cv.glmnet`, k-fold cross validation is
#' performed. This produces a plot and returns optimal values for `lambda`, the
#' penalization parameter. Input can be as typically done for
#' `cv.glmnet` in the form of `x` and `y` which are a matrix and response object
#'  or with a landmark super dataset specifying the dependent columns in `y`.
#'
#' @param x An "LMdataframe", which can be created by running [stack_data()] and
#'   [add_interactions()].
#' @param y Optional, a vector of column names of the data stored in `lmdata`
#'   that are to be used as dependent variables. If not specified, it is assumed
#'   that all non-response variables are the dependent variables.
#' @param id_col Column name or index that identifies individuals in data.
#'   Used to ensure individuals appear in the same cross-validation sets.
#' @param alpha The elastic net mixing parameter: Lies between 0 and 1. At 1,
#'   the penalty is the LASSO penalty, and at 0, the penalty is the ridge
#'   penalty. The default is 1.
#' @param nfolds Number of folds in k-fold cross validation. Default is 10.
#' @param type.measure Loss for cross-validation. Currently the only option is
#'   "deviance" which is the partial-likelihood for the Cox model. If using
#'   cause-specific Cox models, this is evaluated on each model separately.
#' @param seed Set a seed.
#' @param foldid Optional, specify which fold each individual is in.
#' @param ... Additional arguments to `cv.glmnet()`.
#'
#' @return An object of class cv.pen_lm. This is a list of cv.glmnet objects
#'   (one for each cause-specific Cox model or a list of length one for a
#'   regular Cox model). The object also has attributes `survival.type`
#'   (`competing.risk` or `survival`) and `lmdata` and `xcols` which store the
#'   inputs if given.
#'   Functions `print()` and `plot()` exist for the object. To make predictions,
#'   see [dynamic_lm.cv.pen_lm()].
#' @seealso [print.cv.pen_lm()], [plot.cv.pen_lm()], [dynamic_lm.cv.pen_lm()]
#' @examples
#' \dontrun{
#' data(relapse)
#' outcome <- list(time = "Time", status = "event")
#' covars <- list(fixed = c("male", "stage", "bmi"),
#'                varying = c("treatment"))
#' w <- 60; lms <- c(0, 6, 12, 18)
#' lmdata <- stack_data(relapse, outcome, lms, w, covars, format = "long",
#'                      id = "ID", rtime = "T_txgiven")
#' lmdata <- add_interactions(lmdata, func_covars = c("linear", "quadratic"),
#'                            func_lms = c("linear", "quadratic"))
#'
#' # use all covariates
#' cv_model <- cv.pen_lm(lmdata, alpha = 1)
#' print(cv_model, all_causes = TRUE)
#'
#' par(mfrow = c(1, 2))
#' plot(cv_model, all_causes = TRUE)
#'
#' # only use a subset of covariates
#' cv_model1 <- cv.pen_lm(lmdata, y = c("male", "male_LM1", "male_LM2",
#'                                      "stage", "stage_LM1", "stage_LM2"))
#' }
#' @import glmnet riskRegression
#' @export
cv.pen_lm <- function(x, y,
                      id_col,
                      alpha = 1,
                      nfolds = 10,
                      type.measure = "deviance",
                      seed = NULL,
                      foldid = NULL,
                      ...
) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package \"glmnet\" must be installed to use function cv.pen_lm.",
         call. = FALSE)
  }

  args <- match.call()    # store call
  checked_input <- args   # to be cleaned

  # Save the arguments in 'evaluated' form, so we don't have to find them later
  # & remove further arguments we will not lead in later functions
  args$x <- args$id_col <- args$nfolds<- NULL
  args$seed <- args$foldid <- args$type.measure <- NULL
  evaluated_args <- as.list(args)
  function_name <- evaluated_args[[1]]
  evaluated_args <- lapply(evaluated_args[-1], eval)
  evaluated_args <- c(list(function_name), evaluated_args)
  args <- as.call(evaluated_args)

  # Clean input
  m <- match(c("x", "y", "id_col", "alpha"), names(checked_input), 0L)
  checked_input <- as.list(checked_input[m])
  checked_input$CV <- TRUE
  checked_input <- do.call(check_penlm_inputs, checked_input,
                           envir = parent.frame())
  if (inherits(checked_input, "cv.pen_lm")) return(checked_input)
  x <- checked_input$x
  y <- checked_input$y
  lmdata <- checked_input$lmdata
  xcols <- checked_input$xcols
  IDs <- checked_input$IDs
  alpha <- checked_input$alpha
  id_col <- checked_input$id_col
  unique.IDs <- unique(IDs)

  # create foldids (as all instances of an individual must be in the same fold)
  if (missing(foldid)) {
    split.method <- riskRegression::getSplitMethod(paste0("cv", nfolds),
                                                   B = 1,
                                                   N = length(unique.IDs),
                                                   seed = seed)$index
    split.idx <- split.method(1)
    foldid <- split.idx[match(IDs, unique.IDs)]
  }

  # Get the cross-validations
  # TODO: add parallelization!
  models <- lapply(y, function(yi) {
    glmnet::cv.glmnet(x = x, y = yi,
                      family = "cox",
                      alpha = alpha,
                      type.measure = type.measure,
                      foldid = foldid, ...)
  })

  # Additional information for output
  if (length(models) > 1){
    attr(models, "survival.type") <- "competing.risk"
  } else {
    attr(models, "survival.type") <- "survival"
  }
  lmdata$id_col <- id_col
  attr(models, "lmdata") <- lmdata
  attr(models, "xcols") <- xcols
  attr(models, "alpha") <- alpha
  attr(models, "args") <- args
  class(models) <- "cv.pen_lm"

  return(models)
}
