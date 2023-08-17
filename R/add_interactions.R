#' Add landmarking time interactions to a super dataset
#'
#' The stacked dataset output is used as input to [dynamic_lm()] to fit a landmark
#' supermodel for dynamic prediction.
#'
#' @param lmdata An object of class "LMdataframe".
#'
#'   This can be created by running [dynamicLM::stack_data()], or creating a
#'   stacked data set and storing it in a list with attributes outcome, w and
#'   end_time (see [dynamicLM::stack_data()] for further description of outcome
#'   and w), end_time is the largest landmarking time.
#'
#' @param lm_covs Vector of strings indicating the columns (covariates) that are
#'   to have an interaction with the landmark times.
#' @param func_covars Either a string (or vector of strings) specifying which
#'   covariate(`x`)-landmark(`t`) interactions to include. One or
#'   multiple of "linear" (`x, x*t`), "quadratic" (`x, x*t^2`),
#'   "log" (`x, log(1 + x)`), or or "exp" (`x, exp(x)`).
#'
#'   Otherwise, a custom list of functions can be specified. For example,
#'   `list( function(t) t, function(t) exp(20*t))` will, for each covariate `x`,
#'   create `x, x*t, exp(20*t)`.
#' @param func_lms Similar to `func_covars`: A list of functions to use for
#'   transformations of the landmark times. Either a string or vector of
#'   strings or a custom list of functions.
#' @param lm_col Character string specifying the column name that indicates the
#'   landmark time point for a row. Obtained from `lmdata` if not input.
#' @param keep Boolean value to indicate whether or not to keep the columns
#'   given by `lm_covs` without the time interactions. Default is TRUE.
#'
#' @return An object of class "LMdataframe" which now also contains LM
#'   time-interactions.
#'   The object has the following components:
#'   - w, outcome: as the input (obtained from lmdata)
#'   - func_covars: as the input
#'   - func_lms: as the input
#'   - lm_covs: as the input
#'   - all_covs: a list of the new columns added. This includes `lm_covs`
#'     if `keep` is TRUE.
#'   - lm_col: as the input
#'
#' @details For each variable "var" in `lm_covs`, new columns var_1,...,var_i
#'   (length(func_covars) == i) are added; one column for each interaction given
#'   in func_covars is added.
#'
#'   Transformations of the LM column are added and labelled as LM_1,...,LM_j
#'   (length(func_lms) == j); one column for each interaction given in func_lms
#'   is added.
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
#' # Update complex landmark-varying covariates
#' # note age is in years and LM is in months
#' lmdata$data$age <- lmdata$data$age.at.time.0 + lmdata$data$LM/12
#' # Add LM-time interactions
#' lmdata <- add_interactions(lmdata, pred_covars,
#'                            func_covars = c("linear", "quadratic"),
#'                            func_lms = c("linear", "quadratic"))
#' head(lmdata$data)
#' }
#'
#' @export
add_interactions <- function(lmdata, lm_covs, func_covars, func_lms, lm_col,
                             keep = T){
  if (missing(lm_col)){
    lm_col <- lmdata$lm_col
  }
  if (lm_col %in% func_covars){
    stop(paste0("arg lm_col (given as/inferred as ",lm_col,
                ") should not be in arg func_covars."))
  }
  data <- lmdata$data

  f1 <- function(t) t
  f2 <- function(t) t^2
  f3 <- function(t) log(1 + t)
  f4 <- function(t) exp(t)

  if (missing(func_covars)) func_covars <- list(f1, f2)
  if (missing(func_lms)) func_lms <- list(f1, f2)
  if (inherits(func_covars, "character")){
    funcs <- list()
    if ("linear" %in% func_covars) funcs <- c(funcs, list(f1))
    if ("quadratic" %in% func_covars) funcs <- c(funcs, list(f2))
    if ("log" %in% func_covars) funcs <- c(funcs, list(f3))
    if ("exp" %in% func_covars) funcs <- c(funcs, list(f4))
    func_covars <- funcs
  }
  if (inherits(func_lms, "character")){
    funcs <- list()
    if ("linear" %in% func_lms) funcs <- c(funcs, list(f1))
    if ("quadratic" %in% func_lms) funcs <- c(funcs, list(f2))
    if ("log" %in% func_lms) funcs <- c(funcs, list(f3))
    if ("exp" %in% func_lms) funcs <- c(funcs, list(f4))
    func_lms <- funcs
  }

  all_covs <- c(lm_covs)
  data_LM <- data[[lm_col]]
  # Add func_covarss: covariate LM interactions
  for(i in 1:length(lm_covs)){
    for (j in 1:length(func_covars)){
      f <- func_covars[[j]]
      name <- paste0(lm_covs[i],"_",j)
      data[[name]]  <- data[[lm_covs[i]]]*f(data_LM)
      all_covs <- c(all_covs, name)
    }
  }
  # Add func_lms: LM interactions
  for (k in 1:length(func_lms)){
    g <- func_lms[[k]]
    name <- paste0("LM_",k)
    data[[name]]  <- g(data_LM)
    all_covs <- c(all_covs, name)
  }

  if(!keep){
    remaining = colnames(data)[! colnames(data)  %in% lm_covs]
    data <- data[remaining]
  }
  lmdata$data <- data
  lmdata$func_covars <- func_covars
  lmdata$func_lms <- func_lms
  lmdata$lm_covs <- lm_covs
  lmdata$all_covs <- unique(all_covs)
  lmdata$lm_col <- lm_col

  return(lmdata)
}


