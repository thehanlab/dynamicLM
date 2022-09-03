#' Fit a coxph or CSC model to a LM super dataset
#'
#' @param formula The formula to be used, remember to include "+cluster(ID)" for the column that indicates the ID of the individual for robust error estimates.
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
#'   - model: fitted model
#'   - type: as input
#'   - w, func_covars, func_LMs, LMcovars, allLMcovars, outcome: as in LMdata
#'   - LHS: the LHS of the input formula
#'   - linear.predictors: the vector of linear predictors, one per subject. Note that this vector has not been centered.
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
#' # Choose covariates that will have time interaction
#' pred.covars <- c("age","male","stage","bmi","treatment")
#' # Stack landmark datasets
#' LMdata <- cutLMsuper(relapse, outcome, LMs, w, covs, format="long",
#'                      id="ID", rtime="fup_time", right=F)
#' # Update complex LM-varying covariates, note age is in years and LM is in months
#' LMdata$LMdata$age <- LMdata$LMdata$age.at.time.0 + LMdata$LMdata$LM/12
#' # Add LM-time interactions
#' LMdata <- addLMtime(LMdata, pred.covars, func.covars, func.LMs)
#' formula <- "Hist(Time, event, LM) ~ age + male + stage + bmi + treatment +
#'            age_1 + age_2 + male_1 + male_2 + stage_1 + stage_2 + bmi_1 +
#'            bmi_2 + treatment_1 + treatment_2 + LM_1 + LM_2 + cluster(ID)"
#' supermodel <- fitLM(as.formula(formula), LMdata, "CSC")
#' }
#' @import survival
#' @export
#'
fitLM <- function(formula, LMdata, type="coxph", method="breslow",
                  func_covars, func_LMs, LM_col, outcome, w, LMcovars, ...){

  LHS = Reduce(paste, deparse(formula[[2]]))

  if (!grepl("cluster", as.character(stats::as.formula(formula))[3])){
    message("Did you forget to add a '+ cluster(ID)' term for your ID variable in your formula? No cluster argument was specified in the formula. Standard errors may be estimated incorrectly.")
  }

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
    original.landmarks <- data[[LMdata$LM_col]]
    end_time <- LMdata$end_time
    outcome <- LMdata$outcome
    w <- LMdata$w
    LMcovars <- LMdata$LMcovars
    allLMcovars <- LMdata$allLMcovars
  }
  LMcovars <- intersect(sub("_[^_]+$", "", all.vars(formula[[3]])), LMcovars)

  num_preds <- nrow(data)

  if(type=="coxph"){
    superfm <- survival::coxph(formula, data, method=method, ...)
    models <- list(superfm)
    num_causes <- 1
    superfm$call$data <- data
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

  data = data
  linear.predictors <-
    t(sapply(1:num_causes, function(c) {
      coefs <- models[[c]]$coefficients
      df <- data[,names(coefs)]
      rowSums(
        data.frame(mapply(`*`,df,coefs))
        )
      })
    )

  out=list(model=superfm,
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
