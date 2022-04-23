#' fit a coxph or CSC model to a LM super dataset
#'
#' @param formula The formula to be used, remember to include "+cluster(ID)" for the column that indicates the ID of the individual for robust estimates.
#' @param LMdata  An object of class "LM.data.frame", this can be created by running cutLMsuper and addLMtime
#' @param type "coxph" or "CSC"/"CauseSpecificCox"
#' @param method A character string specifying the method for tie handling. Default is "breslow". More information can be found in coxph.
#' @param ... Arguments given to coxph or CSC.
#'
#' @return An object of class "LMcoxph" or "LMCSC" with components:
#' - superfm: fitted model
#' - type: as input
#' - w, func_covars, func_LMs, LMcovars, allLMcovars, outcome: as in LMdata
#' - LHS: the LHS of the input formula
#' - linear.predictors: the vector of linear predictors, one per subject. Note that this vector has not been centered.
#' - original.landmarks: the LM time point at which prediction was made, one per subject. This has the same order as linear.predictors.
#' @export
#'
fitLM <- function(formula, LMdata, type="coxph", method="breslow", ...){
  # TODO: add FGR?
  # TODO: allow for LMdata to be a normal data.frame and specify all the other parts manually

  LHS = Reduce(paste, deparse(formula[[2]]))

  if(class(LMdata)!="LM.data.frame"){
    stop("data must be of type LM.data.frame")
  }
  data=LMdata$LMdata
  num_preds <- nrow(data)

  if(type=="coxph"){
    superfm <- survival::coxph(formula, data, method=method, ...)
    num_causes <- 1
    models <- list(superfm)
    cl <- "LMcoxph"

  } else if (type=="CauseSpecificCox" | type=="CSC"){
    superfm <- riskRegression::CSC(formula, data, method=method, ...)
    models <- superfm$models
    num_causes <- length(models)
    superfm$call$data <- data
    cl <- "LMCSC"
  }

  func_covars=LMdata$func_covars
  func_LMs=LMdata$func_LMs
  original.landmarks=LMdata$LMdata[[LMdata$LM_col]]

  linear.predictors <- sapply(1:num_preds, function(i){
    sapply(1:num_causes, function(c) {
      coefs <- models[[c]]$coefficients
      sum(
        sapply(names(coefs),
               function(coef_name){ coefs[coef_name] * data[i,coef_name] })
      )
    })
  })

  out=list(superfm=superfm,
           type=type,
           w = LMdata$w,
           end_time=LMdata$end_time,
           func_covars=func_covars,
           func_LMs=func_LMs,
           LMcovars=LMdata$LMcovars,
           allLMcovars=LMdata$allLMcovars,
           outcome=LMdata$outcome,
           LHS=LHS,
           linear.predictors=linear.predictors,
           original.landmarks=original.landmarks
  )
  class(out)=cl

  return(out)
}
