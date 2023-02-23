#' Cross-validation for a  regular or cause-specific Cox landmark supermodel
#' with lasso or elasticnet penalization.
#'
#' Fit by calling `cv.glmnet`. As in `cv.glmnet`, k-fold cross validation is
#' performed. This produces a plot and returns a value for lambda. Input can be
#' as typically done for `cv.glmnet` in the form of `x` and `y` or with a
#' landmark super dataset `LMdata` specifying dependent columns in `xcols`.
#'
#' @param x Input matrix with each row an observation.
#' @param y Response variable: either a `Surv` or `Hist` object.
#' @param LMdata An object of class "LMdataframe", this can be created by
#'   running `cutLMsuper` and `addLMtime`
#' @param xcols A vector of column names of the data stored in `LMdata` that
#'   are to be used as dependent variables. If not specified, it is assumed that
#'   all non-response variables are the dependent variables.
#' @param ID_col Column name or index that identifies individuals in data.
#'   Used to ensure the same individuals appear in the different cross
#'   validation sets.
#' @param alpha The elastic net mixing parameter: Lies betweent 0 and 1. At 1,
#'   the penalty is the LASSO penalty, and at 0, the penalty is the ridge
#'   penalty.
#' @param nfolds Number of folds in k-fold cross validation. Default is 10.
#' @param type.measure Loss for cross-validation. Currently the only option is
#'   "deviance" which is the partial-likelihood for the Cox model. If using
#'   cause-specific Cox models, this is evaluated on each model separately.
#' @param seed Set a seed.
#' @param ... Additional arguments to `cv.glmnet`.
#'
#' @return An object class `cv.penLM`. This is a list of `cv.glmnet` objects
#'   (one for each cause-specific Cox model or a list of length one for a
#'   regular Cox model). The object also has attributes `survival.type`
#'   (`competing.risk` or `survival`) and `LMdata` and `xcols` which store the
#'   inputs if given.
#'   Functions `print` and `plot` exist for the object. To make predictions,
#'   see `fitLM.cv.penLM`.
#' @import glmnet riskRegression
#' @export
#'
#' @examples
#' \dontrun{
#' data(relapse)
#' outcome = list(time="Time", status="event")
#' covars = list(fixed=c("ID","age.at.time.0","male","stage","bmi"),
#'               varying=c("treatment"))
#' w = 60; LMs = c(0,12,24)
#' # Covariate-landmark time interactions
#' func.covars <- list( function(t) t, function(t) t^2)
#' # let hazard depend on landmark time
#' func.LMs <- list( function(t) t, function(t) t^2)
#'
#' # Stack landmark datasets
#' LMdata <- cutLMsuper(relapse, outcome, LMs, w, covars, format="long",
#'                      id="ID", rtime="T_txgiven", right=F)
#' # Update complex LM-varying covariates, note age is in years and LM is in months
#' LMdata$LMdata$age <- LMdata$LMdata$age.at.time.0 + LMdata$LMdata$LM/12
#' # Choose covariates that will have time interaction
#' pred.covars <- c("age","male","stage","bmi","treatment")
#' # Add LM-time interactions
#' LMdata <- addLMtime(LMdata, pred.covars, func.covars, func.LMs)
#'
#' xcols <- c("age","male","stage","bmi","treatment","age_1","age_2","male_1",
#'            "male_2","stage_1","stage_2","bmi_1","bmi_2","treatment_1",
#'            "treatment_2","LM_1","LM_2")
#'
#' pen_supermodel <- cv.penLM(LMdata, xcols)
#' print(pen_supermodel)
#' plot(pen_supermodel)
#' }
cv.penLM <- function(x, y,
                     LMdata, xcols,
                     ID_col,
                     alpha=1,
                     nfolds=10,
                     type.measure="deviance",
                     seed=NULL,
                     foldid=NULL,
                     ...
) {
  checked_input <- match.call()
  checked_input$parent_func = quote(cv.penLM)
  checked_input$CV = TRUE
  checked_input[[1L]] <- quote(check_penLM_inputs)
  checked_input <- eval(checked_input, parent.frame())

  if(class(checked_input) == "cv.penLM") return(checked_input)

  x = checked_input$x
  y = checked_input$y
  LMdata = checked_input$LMdata
  xcols = checked_input$xcols
  IDs = checked_input$IDs
  alpha = checked_input$alpha
  ID_col = checked_input$ID_col
  unique.IDs = unique(IDs)

  # create foldids (as all instances of an individual must be in the same fold)
  if(missing(foldid)){
    split.method <- paste0("cv",nfolds)
    split.idx <- riskRegression::getSplitMethod(split.method, B=1, N=length(unique.IDs), seed=seed)$index
    foldid <- split.idx[match(IDs, unique.IDs)]
  }

  models <- lapply(y, function(yi) {
    glmnet::cv.glmnet(x = x, y = yi,
                      family = "cox",
                      alpha = alpha,
                      type.measure = type.measure,
                      foldid = foldid, ...)
  })
  if (length(models) > 1){
    attr(models, "survival.type") = "competing.risk"
  } else {
    attr(models, "survival.type") = "survival"
  }
  if (!is.null(LMdata)) {
    LMdata$ID_col <- ID_col
    attr(models, "LMdata") = LMdata
  }
  if (!is.null(xcols)) attr(models, "xcols") = xcols
  attr(models, "alpha") = alpha
  class(models) = "cv.penLM"
  return(models)
}
