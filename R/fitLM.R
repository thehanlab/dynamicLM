#' Fit a landmark supermodel
#'
#' Fit a Cox or cause-specific Cox landmark supermodel, with or without
#'   regularization.
#'
#' @details To fit coxph or cause-specific Cox models on a landmark super
#' dataset, see `fitLM.LMdataframe`. To fit penalized Cox or penalized
#' cause-specific Cox models see `fitLM.LMpen`.
#'
#' @return An object of class "LMcoxph", "LMCSC", "penLMcoxph" or "penLMCSC"
#'   with components:
#'   - model: fitted model
#'   - type: as input
#'   - w, func_covars, func_LMs, LMcovars, allLMcovars, outcome: as in LMdata
#'   - LHS: the LHS of the input formula
#'   - linear.predictors: the vector of linear predictors, one per subject. Note
#'     that this vector has not been centered.
#' @seealso `fitLM.LMdataframe`, `fitLM.LMpen`
#' @export
#'
# TODO: regularized additions
fitLM <- function(...) {
  UseMethod("fitLM")
}


#' Fit a coxph or CSC model to a LM super dataset
#'
#' @param LMdata  An object of class "LMdataframe", this can be created by
#'   running cutLMsuper and addLMtime
#' @param formula The formula to be used, remember to include "+cluster(ID)" for
#'   the column that indicates the ID of the individual for robust error
#'   estimates.
#' @param type "coxph" or "CSC"/"CauseSpecificCox"
#' @param method A character string specifying the method for tie handling.
#'   Default is "breslow". More information can be found in coxph.
#' @param func_covars A list of functions to use for interactions between LMs
#'   and covariates.
#' @param func_LMs A list of functions to use for transformations of the
#'   landmark times.
#' @param LM_col Character string specifying the column name that indicates the
#'   landmark time point for a row.
#' @param outcome List with items time and status, containing character strings
#'   identifying the names of time and status variables, respectively, of the
#'   survival outcome
#' @param w Scalar, the value of the prediction window (ie predict w-year/other
#'   time period risk from the LM points)
#' @param LMcovars Vector of strings indicating the columns that are to have a
#'   LM interaction
#' @param cluster Variable which clusters the observations (for e.g., identifies
#'   repeated patient IDs), for the purposes of a robust variance.
#' @param x Logical value. If set to true, the LMdata is stored in the returned
#'   object. This is required for internal validation.
#' @param ... Arguments given to coxph or CSC.
#'
#' @return An object of class "LMcoxph" or "LMCSC" with components:
#'   - model: fitted model
#'   - type: as input
#'   - w, func_covars, func_LMs, LMcovars, allLMcovars, outcome: as in LMdata
#'   - LHS: the LHS of the input formula
#'   - linear.predictors: the vector of linear predictors, one per subject. Note
#'     that this vector has not been centered.
#'   - args: arguments used to call model fitting
#'   - ID_col: the cluster argument, usually specifies column with patient ID
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
#' formula <- "Hist(Time, event, LM) ~ age + male + stage + bmi + treatment +
#'            age_1 + age_2 + male_1 + male_2 + stage_1 + stage_2 + bmi_1 +
#'            bmi_2 + treatment_1 + treatment_2 + LM_1 + LM_2 + cluster(ID)"
#' supermodel <- fitLM(LMdata, as.formula(formula), "CSC")
#' }
#' @import survival
#' @export
fitLM.LMdataframe <- function(LMdata,
                              formula,
                              type = "coxph",
                              method = "breslow",
                              func_covars,
                              func_LMs,
                              LM_col,
                              outcome,
                              w,
                              LMcovars,
                              cluster,
                              x = FALSE,
                              ...) {

  # store arguments but not the data (heavy)
  args = match.call()
  args$LMdata = NULL

  # extra LHS of formula
  LHS = getLHS(formula)

  # check if we have a cluster term and extract it
  cluster_check <- as.character(stats::as.formula(formula))[3]
  if (!grepl("cluster", cluster_check)){
    if (missing(cluster)){
      message("Did you forget to specify a cluster argument or add a '+ cluster(ID)' term for your ID variable in your formula? No cluster argument was specified in the formula. Standard errors may be estimated incorrectly.")
    }
  } else {
    cluster <- regmatches(
      cluster_check,
      gregexpr("(?<=cluster\\().*?(?=\\))", cluster_check, perl=T)
    )[[1]]
  }
  ID_col <- cluster

  if(class(LMdata)!="LMdataframe"){
    if(class(LMdata)!="data.frame"){stop("data must be a data frame or an object of class LMdataframe")}

    if(missing(func_covars)) stop("For input data that is a data frame, arg func_covars must be specified.")
    if(missing(func_LMs)) stop("For input data that is a data frame, arg func_LMs must be specified.")
    if(missing(LM_col)) stop("For input data that is a data frame, arg LM_col must be specified.")
    if(missing(outcome)) stop("For input data that is a data frame, arg outcome must be specified.")
    if(missing(w)) stop("For input data that is a data frame, arg w must be specified.")
    if(missing(LMcovars)) stop("For input data that is a data frame, arg LMcovars must be specified.")

    allLMcovars <- c(sapply(1:length(func_covars), function(i) paste0(LMcovars,"_",i)),
                     sapply(1:length(func_LMs), function(i) paste0("LM_",i)))
    if (!all(allLMcovars %in% colnames(LMdata))){
      stop(paste0("The data should have all of the following column names: ",
                  paste0(allLMcovars, collapse=", ")))
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
      stop("Package \"riskRegression\" must be installed to use this function.",
           call. = FALSE)}

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
           ID_col = ID_col,
           linear.predictors=linear.predictors,
           original.landmarks=original.landmarks,
           args=args
  )
  if (x == TRUE) out$LMdata = LMdata
  class(out)=cl

  return(out)
}


#' Fit a penalized coxph or CSC supermodel for a specific coefficient
#'
#' @details The Breslow method is used for handling ties, as we use the `glmnet`
#'   package which does the same.
#'
#' @param object  A fitted object of class "LMpen". This can be created by
#'   calling `penLM` using arguments `LMdata` and `xcols`.
#' @param lambda Value of the penalty parameter `lambda` at which to fit a model. For
#'   cause specific Cox super models, this must be a list or vector of values:
#'   one for each cause.
#'
#' @return An object of class "penLMcoxph" or "penLMCSC" with components:
#'   - model: fitted model
#'   - type: as input
#'   - w, func_covars, func_LMs, LMcovars, allLMcovars, outcome: as in LMdata
#'   - LHS: the LHS of the input formula
#'   - linear.predictors: the vector of linear predictors, one per subject. Note
#'     that this vector has not been centered.
#'   - lambda: the values of lambda for which this model has been fit.
#' @examples
#' \dontrun{
#' }
#' @import survival glmnet
#' @export
fitLM.penLM <- function(object, lambda, ...){
  survival.type <- attr(object, "survival.type")
  LMdata <- attr(object, "LMdata")
  if(is.null(LMdata))
    stop("To fit a penalized model, penLM or cv.penLM must be called with arguments LMdata and xcols. (x,y) is not yet implemented")
  xcols <- attr(object, "xcols")
  data <- LMdata$LMdata
  NC <- length(object)

  # TOOD: check inputs -
  # * is LMdata an LMdataframe etc?
  # * is s the correct length?
  # * allow for (x,y)
  func_covars <- LMdata$func_covars
  func_LMs <- LMdata$func_LMs
  original.landmarks <- data[[LMdata$LM_col]]
  end_time <- LMdata$end_time
  outcome <- LMdata$outcome
  w <- LMdata$w
  LMcovars <- LMdata$LMcovars
  allLMcovars <- LMdata$allLMcovars

  if(survival.type=="survival"){
    if (class(lambda) == "list") lambda <- lambda[[1]]
    glmnet_coefs <- as.vector(coef(object[[1]], s = lambda))

    entry = LMdata$LM_col
    exit = LMdata$outcome$time
    status = LMdata$outcome$status
    LHS_surv <- paste0("Surv(",entry,",",exit,",",status,")")
    formula <- paste0(LHS_surv, "~", paste0(xcols, collapse="+"))
    superfm <- survival::coxph(as.formula(formula), data, method="breslow",
                               iter.max=0, init = glmnet_coefs, ...)

    LHS <- paste0("Hist(",exit,",",status,",",entry,") ~ 1")
    models <- list(superfm)
    num_causes <- 1
    superfm$call$data <- data
    type <- "coxph"
    cl <- "penLMcoxph"
  }

  else if (survival.type=="competing.risk"){
    if (!requireNamespace("riskRegression", quietly = TRUE)) {
      stop("Package \"riskRegression\" must be installed to use this function.",
           call. = FALSE)}

    glmnet_coefs <- lapply(1:NC, function(i){
      as.vector(coef(object[[i]], s = s[[i]]))
    })

    entry = LMdata$LM_col
    exit = LMdata$outcome$time
    status = LMdata$outcome$status
    LHS <- paste0("Hist(",exit,",",status,",",entry,")")
    formula <- paste0(LHS, "~", paste0(xcols, collapse="+"))
    superfm <- CSC.fixed.coefs(as.formula(formula), data, method="breslow",
                               cause.specific.coefs=glmnet_coefs, ...)

    LHS <- paste0(LHS," ~ 1")
    models <- superfm$models
    num_causes <- NC
    superfm$call$data <- data
    type <- "CauseSpecificCox"
    cl <- "penLMCSC"
  }

  # no std errors for penalized Cox
  for(i in 1:num_causes){
    models[[i]]$var <- NULL
  }

  # calculate LPs
  linear.predictors <-
    t(sapply(1:num_causes, function(c) {
      coefs <- models[[c]]$coefficients
      df <- data[,names(coefs)]
      rowSums(data.frame(mapply(`*`,df,coefs)))
    }))

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
           original.landmarks=original.landmarks,
           lambda=lambda
  )
  class(out)=cl
  return(out)
}


#' Fit a penalized cross-validated coxph or CSC super model
#'
#' @details The Breslow method is used for handling ties, as we use the `glmnet`
#'   package which does the same.
#'
#' @param object  A fitted object of class "cv.LMpen". This can be created by
#'   calling `cv.penLM` using arguments `LMdata` and `xcols`.
#' @param lambda Value of the penalty parameter `lambda` at which to fit a model.
#'   Default is "lambda.1se" stored in the object; "lambda.min" can also be used
#'   or a specific value of can be input. For cause-specific Cox super models,
#'   this must be a list or vector of values: one for each cause.
#'
#' @return An object of class "penLMcoxph" or "penLMCSC" with components:
#'   - model: fitted model
#'   - type: as input
#'   - w, func_covars, func_LMs, LMcovars, allLMcovars, outcome: as in LMdata
#'   - LHS: the LHS of the input formula
#'   - linear.predictors: the vector of linear predictors, one per subject. Note
#'     that this vector has not been centered.
#'   - s: the values of lambda for which this model has been fit.
#' @examples
#' \dontrun{
#' }
#' @import survival glmnet
#' @export
fitLM.cv.penLM <- function(object, lambda="lambda.1se", ...){
  if (class(lambda) == "character"){
    if (lambda == "lambda.1se")
      lambda <- lapply(object, function(o) o$lambda.1se)
    else if (lambda == "lambda.min")
      lambda <- lapply(object, function(o) o$lambda.min)
  }

  return(fitLM.penLM(object, lambda=lambda, ...))
}
