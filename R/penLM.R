#' Fit a regular or cause-specific Cox landmark supermodel with lasso or
#' elasticnet penalization.
#'
#' Fit by calling `glmnet`. As in `glmnet`, the model is fit via penalized
#' maximum likelihood to produce a regularization path at a grid of values for
#' the regularization parameter lambda. Input can be as typically done for
#' `glmnet` in the form of `x` and `y` or with a landmark super dataset `LMdata`
#' specifying dependent columns in `xcols`.
#'
#' @param x Input matrix with each row an observation.
#' @param y Response variable: either a `Surv` or `Hist` object.
#' @param LMdata An object of class "LMdataframe", this can be created by
#'   running `cutLMsuper` and `addLMtime`
#' @param xcols A vector of column names of the data stored in `LMdata` that
#'   are to be used as dependent variables. If not specified, it is assumed that
#'   all non-response variables are the dependent variables.
#' @param ... Additional arguments to `glmnet`.
#'
#' @return An object class `penLM`. This is a list of `glmnet` objects (one for
#'   each cause-specific Cox model or a list of length one for a regular Cox
#'   model). The object also has attributes `survival.type` (`competing.risk`
#'   or `survival`) and `LMdata` and `xcols` which store the inputs if given.
#'   Functions `print` and `plot` exist for the object. To make predictions,
#'   see `fitLM.penLM`.
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
#' pen_supermodel <- penLM(LMdata, xcols)
#' print(pen_supermodel)
#' plot(pen_supermodel)
#' }
penLM <- function(x, y, LMdata, xcols, ...) {

  checked_input <- match.call()
  checked_input$parent_func = quote(penLM)
  checked_input[[1L]] <- quote(check_penLM_inputs)
  checked_input <- eval(checked_input, parent.frame())

  # can call penLM again depending on inputs
  if(class(checked_input) == "penLM") return(checked_input)

  # if not, use checked inputs
  x = checked_input$x
  y = checked_input$y
  LMdata = checked_input$LMdata
  xcols = checked_input$xcols

  models <- lapply(y, function(yi) {
    glmnet(x = x, y = yi, family = "cox", ...)
  })
  if (length(models) > 1){
    attr(models, "survival.type") = "competing.risk"
  } else {
    attr(models, "survival.type") = "survival"
  }
  if (!is.null(LMdata)) attr(models, "LMdata") = LMdata
  if (!is.null(xcols)) attr(models, "xcols") = xcols
  class(models) = "penLM"
  return(models)
}

