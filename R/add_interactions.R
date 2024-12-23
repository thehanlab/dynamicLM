#' Add landmarking time interactions to a super dataset
#'
#' The stacked dataset output is used as input to [dynamic_lm()] to fit a landmark
#' supermodel for dynamic prediction.
#'
#' @param lmdata An object of class "LMdataframe"
#'
#'   This can be created by running [dynamicLM::stack_data()], or creating a
#'   stacked data set and storing it in a list with attributes outcome, w and
#'   end_time (see [dynamicLM::stack_data()] for further description of outcome
#'   and w), end_time is the largest landmarking time.
#'
#' @param lm_covs Vector of strings indicating the columns (covariates) that are
#'   to have an interaction with the landmark times.
#' @param func_covars Either a string/vector of strings or list of
#'   functions specifying which covariate-landmark interactions to include.
#'   If `x` are covariates and `t` are landmarks then "linear" (`x, x*t`),
#'   "quadratic" (`x, x*t^2`), "log" (`x, log(1 + x)`), or or "exp"
#'   (`x, exp(x)`) can be specified.
#'
#'   A custom list of functions can be specified. For example,
#'   `list(function(t) t, function(t) exp(20*t))` will, for each covariate,
#'   create `x, x*t, exp(20*t)`.
#' @param func_lms A list of functions to use for transformations of the
#'   landmark times input similarly to `func_covars`, either as a string/
#'   vector of strings or a custom list of functions.
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
#' @details For each variable "var" in `lm_covs`, new columns var_LM1, ...,
#'   var_LMi are added; one column for each interaction given
#'   in func_covars is added (length(func_covars) == i).
#'
#'   Transformations of the LM column are added and labelled as LM1, ..., LMj;
#'   one column for each interaction given in func_lms is added
#'   (length(func_lms) == j).
#'
#' @examples
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
#'
#' @seealso [dynamicLM::stack_data()], [dynamicLM::dynamic_lm()]
#' @export
add_interactions <- function(lmdata, lm_covs,
                             func_covars = c("linear", "quadratic"),
                             func_lms = c("linear", "quadratic"),
                             keep = TRUE) {
  # Validate inputs
  if (!inherits(lmdata, "LMdataframe")) stop("`lmdata` must be of class 'LMdataframe'")
  if (missing(lm_covs)) lm_covs <- lmdata$all_covs
  lm_col <- "LM"
  data <- lmdata$data
  if (is.null(data[[lm_col]])) stop(paste0("`", lm_col, "` column not found in data"))

  # Predefined transformation functions
  predefined_funcs <- function(name) {
    switch(name,
           "linear" = function(t) t,
           "quadratic" = function(t) t^2,
           "log" = function(t) log(1 + t),
           "exp" = function(t) exp(t),
           stop(paste("Unsupported transformation:", name))
    )
  }

  # Process landmark-transformations
  if (inherits(func_covars, "character"))
    func_covars <- lapply(func_covars, predefined_funcs)
  else if (!is.list(funcs))
    stop("`func_covars` must be a character vector or list of functions")

  if (inherits(func_lms, "character"))
    func_lms <- lapply(func_lms, predefined_funcs)
  else if (!is.list(func_lms))
    stop("`func_lms` must be a character vector or list of functions")

  # Initialize output covariates tracker
  all_covs <- lm_covs
  data_lm <- data[[lm_col]]

  # Add covariate-landmark interactions
  for (cov in lm_covs) {
    for (j in seq_along(func_covars)) {
      f <- func_covars[[j]]
      name <- paste0(cov, "_LM", j)
      data[[name]]  <- data[[cov]] * f(data_lm)
      all_covs <- c(all_covs, name)
    }
  }
  # Add landmark-time transformations
  for (j in seq_along(func_lms)){
    g <- func_lms[[j]]
    name <- paste0("LM", j)
    data[[name]]  <- g(data_lm)
    all_covs <- c(all_covs, name)
  }

  # Optionally drop original landmark covariates
  if (!keep) {
    remaining <- colnames(data)[! colnames(data)  %in% lm_covs]
    data <- data[remaining]
  }

  # Update and return `lmdata`
  lmdata$data <- data
  lmdata$func_covars <- func_covars
  lmdata$func_lms <- func_lms
  lmdata$lm_covs <- lm_covs
  lmdata$all_covs <- unique(all_covs)
  lmdata$lm_col <- lm_col

  return(lmdata)
}


