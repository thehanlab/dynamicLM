#' fit a coxph or CSC model to a LM super dataset
#'
#' @param formula The formula to be used, remember to include "+cluster(ID)" for the column that indicates the ID of the individual for robust estimates.
#' @param LMdata  An object of class "LMdataframe", this can be created by running cutLMsuper and addLMtime
#' @param type "coxph" or "CSC"/"CauseSpecificCox"
#' @param method A character string specifying the method for tie handling. Default is "breslow". More information can be found in coxph.
#' @param func_covars A list of functions to use for interactions between LMs and covariates.
#' @param func_LMs A list of functions to use for transformations of the landmark times.
#' @param LM_col Character string specifying the column name that indicates the landmark time point for a row.
#' @param outcome List with items time and status, containing character strings identifying the names of time and status variables, respectively, of the survival outcome
#' @param w Scalar, the value of the prediction window (ie predict w-year/other time period risk from the LM points)
#' @param LMcovars Vector of strings indicating the columns that are to have a LM interaction
#' @param ... Arguments given to coxph or CSC.
#'
#' @return An object of class "LMcoxph" or "LMCSC" with components:
#' - superfm: fitted model
#' - type: as input
#' - w, func_covars, func_LMs, LMcovars, allLMcovars, outcome: as in LMdata
#' - LHS: the LHS of the input formula
#' - linear.predictors: the vector of linear predictors, one per subject. Note that this vector has not been centered.
#' - original.landmarks: the LM time point at which prediction was made, one per subject. This has the same order as linear.predictors.
#' @import survival
#' @export
#'
fitLM <- function(formula, LMdata, type="coxph", method="breslow",
                  func_covars, func_LMs, LM_col, outcome, w, LMcovars, ...){
  # TODO: add FGR?

  LHS = Reduce(paste, deparse(formula[[2]]))

  if(class(LMdata)!="LMdataframe"){
    if(class(LMdata)!="data.frame"){stop("data must be of a data.frame or an object of class LMdataframe")}

    if(missing(func_covars)) stop("For input data that is a data frame, arg func_covars must be specified.")
    if(missing(func_LMs)) stop("For input data that is a data frame, arg func_LMs must be specified.")
    if(missing(LM_col)) stop("For input data that is a data frame, arg LM_col must be specified.")
    if(missing(outcome)) stop("For input data that is a data frame, arg outcome must be specified.")
    if(missing(w)) stop("For input data that is a data frame, arg w must be specified.")
    if(missing(LMcovars)) stop("For input data that is a data frame, arg LMcovars must be specified.")

    allLMcovars <- c(sapply(1:length(func_covars), function(i) paste0(LMcovars,"_",i)),
                     sapply(1:length(func_LMs), function(i) paste0("LM_",i)))
    if (!all(allLMcovars %in% colnames(LMdata))){
      stop(paste0("The data should have all of the following column names: ",paste0(allLMcovars, collapse=", ")))
    }

    data <- LMdata
    original.landmarks <- data[[LM_col]]
    end_time <- max(original.landmarks)

  } else {
    data <- LMdata$LMdata
    func_covars <- LMdata$func_covars
    func_LMs <- LMdata$func_LMs
    original.landmarks <- LMdata$LMdata[[LMdata$LM_col]]
    end_time <- LMdata$end_time
    outcome <- LMdata$outcome
    w <- LMdata$w
    LMcovars <- LMdata$LMcovars
    allLMcovars <- LMdata$allLMcovars
  }
  LMcovars <- intersect(sub("_[^_]+$", "", all.vars(formula[[3]])), LMcovars)

  num_preds <- nrow(data)

  if(type=="coxph"){
    superfm <- coxph(formula, data, method=method, ...)
    num_causes <- 1
    models <- list(superfm)
    cl <- "LMcoxph"

  } else if (type=="CauseSpecificCox" | type=="CSC"){
    if (!requireNamespace("riskRegression", quietly = TRUE)) {
      stop("Package \"riskRegression\" must be installed to use this function.", call. = FALSE)}

    superfm <- riskRegression::CSC(formula, data, method=method, ...)
    models <- superfm$models
    num_causes <- length(models)
    superfm$call$data <- data
    cl <- "LMCSC"
  }

  linear.predictors <-
    t(sapply(1:num_causes, function(c) {
      coefs <- models[[c]]$coefficients
      df <- data[,names(coefs)]
      rowSums(
        data.frame(mapply(`*`,df,coefs))
        )
      })
    )

  out=list(superfm=superfm,
           type=type,
           w=w,
           end_time=end_time,
           func_covars=func_covars,
           func_LMs=func_LMs,
           LMcovars=LMcovars,
           allLMcovars=allLMcovars,
           outcome=outcome,
           LHS=LHS,
           linear.predictors=linear.predictors,
           original.landmarks=original.landmarks
  )
  class(out)=cl

  return(out)
}
