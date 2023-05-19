#' Fit a coxph or CSC model to a landmark super dataset, i.e., fit a dynamic
#' landmark supermodel
#'
#' dynamic (dyn) landmark (lm) supermodel --> dynamic_lm
#'
#' @param lmdata  An object of class "LMdataframe", this can be created by
#'   running [dynamicLM::stack_data()] and [dynamicLM::add_interactions()]
#' @param formula The formula to be used, remember to include "+cluster(ID)" for
#'  the column that indicates the ID of the individual for robust error
#'  estimates.
#'  Note that transformations (e.g., `x1*x2`) cannot be used in the formula and
#'  factors/categorical variables must first be made into dummy variables.
#' @param type "coxph" or "CSC"/"CauseSpecificCox"
#' @param method A character string specifying the method for tie handling.
#'   Default is "breslow". More information can be found in coxph.
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
#'   repeated patient IDs), for the purposes of a robust variance.
#' @param x Logical value. If set to true, the `lmdata` is stored in the returned
#'   object. This is required for internal validation.
#' @param ... Arguments given to coxph or CSC.
#'
#' @return An object of class "LMcoxph" or "LMCSC" with components:
#'   - model: fitted model
#'   - type: as input
#'   - w, func_covars, func_lms, lm_covs, all_covs, outcome: as in `lmdata`
#'   - LHS: the LHS of the input formula
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
#' LMs = seq(0, 36, by = 6)
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
#' lmdata <- add_interactions(lmdata, pred_covars, func_covars = "quadratic",
#'                            func_lms = "quadratic")
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
dynamic_lm <-  function(lmdata,
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
  args = match.call()
  args$lmdata = NULL

  # extra LHS of formula
  LHS = getLHS(formula)

  # check if we have a cluster term and extract it
  cluster_check <- as.character(stats::as.formula(formula))[3]
  if (!grepl("cluster", cluster_check)){
    if (missing(cluster)){
      stop("Did you forget to specify a cluster argument or add a '+ cluster(ID)' term for your ID variable in your formula? No cluster argument was specified in the formula. Standard errors may be estimated incorrectly.")
    } else {
      id_col <- cluster
    }
  } else {
    cluster <- regmatches(cluster_check, gregexpr("(?<=cluster\\().*?(?=\\))", cluster_check, perl=T))[[1]]
    id_col <- cluster
  }

  if(!inherits(lmdata,"LMdataframe")){
    if(!inherits(lmdata,"data.frame"))
      stop("data must be of a data.frame or an object of class LMdataframe")

    if (missing(func_covars))
      stop("For input data that is a data frame, arg func_covars must be specified.")
    if (missing(func_lms))
      stop("For input data that is a data frame, arg func_lms must be specified.")
    if (missing(lm_col))
      stop("For input data that is a data frame, arg lm_col must be specified.")
    if (missing(outcome))
      stop("For input data that is a data frame, arg outcome must be specified.")
    if (missing(w))
      stop("For input data that is a data frame, arg w must be specified.")
    if (missing(lm_covs))
      stop("For input data that is a data frame, arg lm_covs must be specified.")

    all_covs <- c(sapply(1:length(func_covars), function(i) paste0(lm_covs,"_",i)),
                     sapply(1:length(func_lms), function(i) paste0("LM_",i)))
    if (!all(all_covs %in% colnames(lmdata))){
      stop(paste0("The data should have all of the following column names: ",paste0(all_covs, collapse=", ")))
    }

    data <- lmdata
    original.landmarks <- data[[lm_col]]
    end_time <- max(original.landmarks)

  } else {
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
  }
  lm_covs <- intersect(sub("_[^_]+$", "", all.vars(formula[[3]])), lm_covs)

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

  out=list(model = superfm,
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
           args = args
  )
  if (x == TRUE) out$data = lmdata
  class(out)=c(cl, "dynamicLM")

  return(out)
}
