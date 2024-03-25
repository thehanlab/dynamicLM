#' Fit a dynamic Cox or cause-specific Cox landmark supermodel with or without
#'   regularization.
#'
#' @description
#' To fit Cox or cause-specific Cox models without regularization see:
#' * [dynamic_lm.LMdataframe()] for use on a stacked landmark dataset
#' * [dynamic_lm.data.frame()] for use on a dataframe
#'
#' To fit penalized Cox or cause-specific Cox models see:
#' * [dynamic_lm.pen_lm()] without cross-validation
#' * [dynamic_lm.cv.pen_lm()] with cross-validation
#'
#' @param ... Arguments to pass to `dynamic_lm()`
#'
#' @return A fitted landmark supermodel object which has components:
#'   - model: fitted model
#'   - type: as input
#'   - w, func_covars, func_lms, lm_covs, all_covs, outcome: as in `lmdata`
#'   - LHS: the survival outcome
#'   - linear.predictors: the vector of linear predictors, one per subject.
#'     Note that this vector has not been centered.
#'
#'  If the model is unpenalized (class "LMcoxph" or "LMCSC") it has additional
#'  components:
#'   - args: arguments used to call model fitting
#'   - id_col: the cluster argument, often specifies the column with patient ID
#'   - lm_col: column name that indicates the landmark time point for a row.
#'
#' If the model is penalized (class "penLMcoxph" or "penLMCSC") it has
#' additional components:
#' - lambda: the values of lambda for which this model has been fit.
#'
#' @seealso [dynamic_lm.LMdataframe()], [dynamic_lm.data.frame()],
#'   [dynamic_lm.pen_lm()], [dynamic_lm.cv.pen_lm()]
#' @export
#'
dynamic_lm <- function(...) {
  UseMethod("dynamic_lm")
}

#' Fit a dynamic Cox or cause-specific Cox landmark supermodel to a stacked
#' landmark dataset
#'
#' @param lmdata  An object of class "LMdataframe", this can be created by
#'   running [dynamicLM::stack_data()] and [dynamicLM::add_interactions()]
#' @param formula The formula to be used, remember to include `+cluster(ID)` for
#'  the column that indicates the ID of the individual for robust error
#'  estimates. See details for further information.
#'  Note that transformations (e.g., `x1*x2`) cannot be used in the formula and
#'  factors/categorical variables must first be made into dummy variables.
#' @param type "coxph" or "CSC"/"CauseSpecificCox"
#' @param method A character string specifying the method for tie handling.
#'   Default is "breslow". More information can be found in [survival::coxph()].
#' @param cluster Variable which clusters the observations (for e.g., identifies
#'   repeated patient IDs), for the purposes of a robust variance. If omitted,
#'   extracted from `formula`.
#' @param x Logical value. If set to true, `lmdata` is stored in the returned
#'   object. This is required for internal validation.
#' @param ... Arguments given to coxph or CSC.
#'
#' @details  For standard survival data (one event and possible censoring), use
#'   `type = "coxph"` and a a formula with left-hand side (LHS) of the form
#'   `Surv(LM, Time, event)`. For competing risks (multiple events and possible
#'   censoring), use `type = "CSC"` and a LHS of the form
#'   `Hist(Time, event, LM)`
#'
#' @return An object of class "LMcoxph" or "LMCSC" with components:
#'   - model: fitted model
#'   - type: as input
#'   - w, func_covars, func_lms, lm_covs, all_covs, outcome: as in `lmdata`
#'   - LHS: the survival outcome
#'   - linear.predictors: the vector of linear predictors, one per subject.
#'     Note that this vector has not been centered.
#'   - args: arguments used to call model fitting
#'   - id_col: the cluster argument, often specifies the column with patient ID
#'   - lm_col: column name that indicates the landmark time point for a row.
#'
#' @examples
#' \dontrun{
#' data(relapse)
#' outcome <- list(time = "Time", status = "event")
#' covars <- list(fixed = c("age.at.time.0", "male", "stage", "bmi"),
#'                varying = c("treatment"))
#' w <- 60; lms <- c(0, 6, 12, 18)
#' # Choose covariates that will have time interaction
#' pred_covars <- c("age", "male", "stage", "bmi", "treatment")
#' # Stack landmark datasets
#' lmdata <- stack_data(relapse, outcome, lms, w, covars, format = "long",
#'                      id = "ID", rtime = "T_txgiven")
#'
#' # Update complex landmark-varying covariates
#' # note age is in years and LM is in months
#' lmdata$data$age <- lmdata$data$age.at.time.0 + lmdata$data$LM/12
#' # Add LM-time interactions
#' lmdata <- add_interactions(lmdata, pred_covars,
#'                            func_covars = c("linear", "quadratic"),
#'                            func_lms = c("linear", "quadratic"))
#'
#' formula <- "Hist(Time, event, LM) ~ age + male + stage + bmi + treatment +
#'            age_1 + age_2 + male_1 + male_2 + stage_1 + stage_2 + bmi_1 +
#'            bmi_2 + treatment_1 + treatment_2 + LM_1 + LM_2 + cluster(ID)"
#' supermodel <- dynamic_lm(lmdata, as.formula(formula), "CSC")
#' print(supermodel)
#'
#' par(mfrow = c(2,3))
#' plot(supermodel)
#' }
#' @import survival
#' @export
#'
dynamic_lm.LMdataframe <-  function(lmdata,
                                    formula,
                                    type = "coxph",
                                    method = "breslow",
                                    cluster,
                                    x = FALSE,
                                    ...) {
  # store arguments but not the data (heavy)
  args <- match.call()
  args$lmdata <- NULL

  # Obtain other inputs
  data <- lmdata$data
  func_covars <- lmdata$func_covars
  func_lms <- lmdata$func_lms
  lm_col <- lmdata$lm_col
  original.landmarks <- data[[lm_col]]
  end_time <- lmdata$end_time
  outcome <- lmdata$outcome
  w <- lmdata$w
  lm_covs <- lmdata$lm_covs
  all_covs <- lmdata$all_covs

  out <- dynamic_lm_helper(formula, type, data, lmdata, method, cluster, x, w,
                           end_time, func_covars, func_lms, lm_covs, all_covs,
                           outcome, lm_col, original.landmarks, args, ...)

  return(out)
}


#' Fit a dynamic Cox or cause-specific Cox landmark supermodel to a dataframe.
#'
#' Note that it is recommended to rather use [stack_data()] and
#' [add_interactions()] to create an object of class LMdataframe rather than
#' directly calling `dynamic_lm()` on a dataframe to ensure the data has the
#' correct form.
#'
#' @param lmdata A dataframe that should be a stacked dataset across landmark
#'   times.
#' @param formula The formula to be used, remember to include `+cluster(ID)` for
#'  the column that indicates the ID of the individual for robust error
#'  estimates. See details for further information.
#'  Note that transformations (e.g., `x1*x2`) cannot be used in the formula and
#'  factors/categorical variables must first be made into dummy variables.
#' @param type "coxph" or "CSC"/"CauseSpecificCox"
#' @param method A character string specifying the method for tie handling.
#'   Default is "breslow". More information can be found in [survival::coxph()].
#' @param func_covars A list of functions to use for interactions between LMs
#'   and covariates.
#' @param func_lms A list of functions to use for transformations of the
#'   landmark times.
#' @param lm_col Character string specifying the column name that indicates the
#'   landmark time point for a row.
#' @param outcome List with items time and status, containing character strings
#'   identifying the names of time and status variables, respectively, of the
#'   survival outcome
#' @param w Scalar, the value of the prediction window (ie predict w-year/other
#'   time period risk from the LM points)
#' @param lm_covs Vector of strings indicating the columns that are to have a
#'   LM interaction
#' @param cluster Variable which clusters the observations (for e.g., identifies
#'   repeated patient IDs), for the purposes of a robust variance. If omitted,
#'   extracted from `formula`.
#' @param x Logical value. If set to true, `lmdata` is stored in the returned
#'   object. This is required for internal validation.
#' @param ... Arguments given to coxph or CSC.
#'
#' @details  For standard survival data (one event and possible censoring), use
#'   `type = "coxph"` and a a formula with left-hand side (LHS) of the form
#'   `Surv(LM, Time, event)`. For competing risks (multiple events and possible
#'   censoring), use `type = "CSC"` and a LHS of the form
#'   `Hist(Time, event, LM)`
#'
#' @return An object of class "LMcoxph" or "LMCSC" with components:
#'   - model: fitted model
#'   - type: as input
#'   - w, func_covars, func_lms, lm_covs, all_covs, outcome: as in input.
#'   - LHS: the survival outcome
#'   - linear.predictors: the vector of linear predictors, one per subject.
#'     Note that this vector has not been centered.
#'   - args: arguments used to call model fitting
#'   - id_col: the cluster argument, often specifies the column with patient ID
#'   - lm_col: column name that indicates the landmark time point for a row.
#'
#' @examples
#' \dontrun{
#' }
#' @import survival
#' @export
#'
dynamic_lm.data.frame <- function(lmdata,
                                  formula,
                                  type = "coxph",
                                  method = "breslow",
                                  func_covars,
                                  func_lms,
                                  lm_col,
                                  outcome,
                                  w,
                                  lm_covs,
                                  cluster,
                                  x = FALSE,
                                  ...) {

  # store arguments but not the data (heavy)
  args <- match.call()
  args$lmdata <- NULL

  # Check arguments are given
  if (missing(func_covars))
    stop(tidymess("For input data that is a dataframe, argument func_covars must
                  be specified."))
  if (missing(func_lms))
    stop(tidymess("For input data that is a dataframe, argument func_lms must be
                  specified."))
  if (missing(lm_col))
    stop(tidymess("For input data that is a dataframe, argument lm_col must be
                  specified."))
  if (missing(outcome))
    stop(tidymess("For input data that is a dataframe, argument outcome must be
                  specified."))
  if (missing(w))
    stop("For input data that is a dataframe, argument w must be specified.")
  if (missing(lm_covs))
    stop(tidymess("For input data that is a dataframe, argument lm_covs must be
                  specified."))

  all_covs <- c(
    sapply(seq_along(func_covars), function(i) paste0(lm_covs, "_", i)),
    sapply(seq_along(func_lms), function(i) paste0("LM_", i))
  )
  if (!all(all_covs %in% colnames(lmdata))) {
    stop(paste0("The data should have all of the following column names: ",
                paste0(all_covs, collapse = ", ")))
  }
  data <- lmdata
  original.landmarks <- data[[lm_col]]
  end_time <- max(original.landmarks)

  out <- dynamic_lm_helper(formula, type, data, lmdata, method, cluster, x, w,
                           end_time, func_covars, func_lms, lm_covs, all_covs,
                           outcome, lm_col, original.landmarks, args, ...)
  return(out)
}


#' Fit a penalized coxph or CSC supermodel for a specific coefficient
#'
#' Use one value of `lambda` to fit a model from which predictions can be made.
#'
#' @details The Breslow method is used for handling ties, as we use the `glmnet`
#'   package which does the same.
#'
#' @param object  A fitted object of class "pen_lm". This can be created by
#'   calling [pen_lm()] using arguments `lmdata` and `xcols`.
#' @param lambda Value of the penalty parameter `lambda` at which to fit a
#'   model. For cause-specific Cox super models, this must be a list or vector
#'   of values: one for each cause.
#' @param ... Additional arguments to pass to [survival::coxph()] or
#'   [riskRegression::CSC()]
#'
#' @return An object of class "penLMcoxph" or "penLMCSC" with components:
#'   - model: fitted model
#'   - type: as input
#'   - w, func_covars, func_lms, lm_covs, all_covs, outcome: as in lmdata.
#'   - LHS: the LHS of the input formula
#'   - linear.predictors: the vector of linear predictors, one per subject. Note
#'     that this vector has not been centered.
#'   - lambda: the values of lambda for which this model has been fit.
#' @examples
#' \dontrun{
#' }
#' @import survival glmnet
#' @export
# TODO: allow for the pen_lm object to have been fit with
# (x,y) = (matrix, response)
dynamic_lm.pen_lm <- function(object, lambda, ...) {
  survival.type <- attr(object, "survival.type")
  lmdata <- attr(object, "lmdata")

  if(is.null(lmdata))
    stop(tidymess("To fit a penalized model, pen_lm() or cv.pen_lm() must be
      called with arguments (x,y) as a LMdataframe and colnames. Matrices and
      responses are not yet implemented"))
  xcols <- attr(object, "xcols")
  data <- lmdata$data
  NC <- length(object)

  if (NC != length(lambda)) stop(tidymess(paste0(
    "One lambda should be specified for each cause. ", length(lambda),
    " lambda(s) were given for ", NC, " cause(s). Please provide ", NC,
    " lambda(s) as a vector or list.")))

  func_covars <- lmdata$func_covars
  func_lms <- lmdata$func_lms
  original.landmarks <- data[[lmdata$lm_col]]
  end_time <- lmdata$end_time
  outcome <- lmdata$outcome
  w <- lmdata$w
  lm_covs <- lmdata$lm_covs
  all_covs <- lmdata$all_covs
  id_col <- lmdata$id_col
  lm_col <- lmdata$lm_col

  if (survival.type == "survival") {
    if (inherits(lambda, "list")) lambda <- lambda[[1]]
    glmnet_coefs <- as.vector(coef(object[[1]], s = lambda))

    entry <- lmdata$lm_col
    exit <- lmdata$outcome$time
    status <- lmdata$outcome$status
    LHS_surv <- paste0("Surv(", entry, ",", exit, ",", status, ")")
    formula <- paste0(LHS_surv, "~", paste0(xcols, collapse = "+"))
    superfm <- survival::coxph(stats::as.formula(formula), data,
                               method = "breslow", iter.max = 0,
                               init = glmnet_coefs, ...)

    LHS <- paste0("Hist(", exit, ",", status, ",", entry, ") ~ 1")
    models <- list(superfm)
    num_causes <- 1
    superfm$call$data <- data
    type <- "coxph"
    cl <- "penLMcoxph"

  } else if (survival.type == "competing.risk") {
    if (!requireNamespace("riskRegression", quietly = TRUE)) {
      stop("Package \"riskRegression\" must be installed to use this function.",
           call. = FALSE)}

    glmnet_coefs <- lapply(1:NC, function(i){
      as.vector(coef(object[[i]], s = lambda[[i]]))
    })

    entry <- lmdata$lm_col
    exit <- lmdata$outcome$time
    status <- lmdata$outcome$status
    LHS <- paste0("Hist(", exit, ",", status, ",", entry, ")")
    formula <- paste0(LHS, "~", paste0(xcols, collapse = "+"))

    superfm <- CSC.fixed.coefs(stats::as.formula(formula), data,
                               method = "breslow",
                               cause.specific.coefs = glmnet_coefs, ...)

    LHS <- paste0(LHS, " ~ 1")
    models <- superfm$models
    num_causes <- NC
    superfm$call$data <- data
    type <- "CauseSpecificCox"
    cl <- "penLMCSC"
  }

  # no std errors for penalized Cox
  for (i in 1:num_causes) {
    models[[i]]$var <- NULL
  }

  # calculate LPs
  linear.predictors <-
    t(sapply(1:num_causes, function(c) {
      coefs <- models[[c]]$coefficients
      df <- data[,names(coefs)]
      rowSums(data.frame(mapply(`*`, df, coefs)))
    }))

  out <- list(model = superfm,
              type = type,
              w = w,
              end_time = end_time,
              func_covars = func_covars,
              func_lms = func_lms,
              lm_covs = lm_covs,
              all_covs = all_covs,
              outcome = outcome,
              LHS = LHS,
              id_col = id_col,
              lm_col = lm_col,
              linear.predictors = linear.predictors,
              original.landmarks = original.landmarks,
              # args = args
              lambda = lambda
  )
  class(out) <- c(cl, "dynamicLM")
  return(out)
}


#' Fit a penalized cross-validated coxph or CSC super model
#'
#' Use one value of `lambda` to fit a model from which predictions can be made.
#'
#' @details The Breslow method is used for handling ties, as we use the `glmnet`
#'   package which does the same.
#'
#' @param object  A fitted object of class "cv.pen_lm". This can be created by
#'   calling `cv.pen_lm` using arguments `lmdata` and `xcols`.
#' @param lambda Value of the penalty parameter `lambda` to fit a model.
#'   Default is "lambda.min"; "lambda.1se" can also be used or a specific value
#'   can be input. For cause-specific Cox super models,
#'   this must be a list or vector of values: one for each cause.
#' @param ... Additional arguments to pass to [survival::coxph()] or
#'   [riskRegression::CSC()]
#'
#' @return An object of class "penLMcoxph" or "penLMCSC" with components:
#'   - model: fitted model
#'   - type: as input
#'   - w, func_covars, func_lms, lm_covs, all_covs, outcome: as in `lmdata`
#'   - LHS: the LHS of the input formula
#'   - linear.predictors: the vector of linear predictors, one per subject. Note
#'     that this vector has not been centered.
#'   - lambda: the values of lambda for which this model has been fit.
#' @examples
#' \dontrun{
#' }
#' @import survival glmnet
#' @export
dynamic_lm.cv.pen_lm <- function(object, lambda = "lambda.min", ...) {
  if (inherits(lambda, "character")) {
    if (lambda == "lambda.1se")
      lambda <- lapply(object, function(o) o$lambda.1se)
    else if (lambda == "lambda.min")
      lambda <- lapply(object, function(o) o$lambda.min)
  }
  return(dynamic_lm.pen_lm(object, lambda = lambda, ...))
}
