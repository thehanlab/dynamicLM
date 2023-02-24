# -----------------------------------------------------------------------
#' Add landmarking time interactions to a super dataset
#'
#' @param lmdata An object of class "LMdataframe".
#'
#'   This can be created by running [dynamicLM::stack_data()], or creating a
#'   stacked data set and storing it in a list with attributes outcome, w and
#'   end_time (see [dynamicLM::stack_data()] for further description of outcome
#'   and w), end_time is the largest landmarking time.
#' @param lm_covs Vector of strings indicating the columns (covariates) that are
#'   to have a LM interaction
#' @param func_covars A list of functions to use for interactions between LMs
#'   and covariates.
#' @param func_lms A list of functions to use for transformations of the
#'   landmark times.
#' @param LM_col Character string specifying the column name that indicates the
#'   landmark time point for a row.
#' @param keep Boolean value to indicate whether or not to keep the columns
#'   given by `lm_covs` without the time interactions or not. Default=FALSE.
#'
#' @return An object of class "LMdataframe" which now also contains LM
#'   time-interactions.
#'   The object has the following components:
#'   - w, outcome: as the input (obtained from lmdata)
#'   - func_covars: as the input
#'   - func_lms: as the input
#'   - lm_covs: as the input
#'   - allLMcovars: a list of the new columns added
#'   - LM_col: as the input
#'
#' @details For each variable "var" in `lm_covs`, new columns var_1,...,var_i
#'   (length(func_covars) == i) are added; one column for each interaction given
#'   in func_covars is added
#'
#'   Transformations of the LM column are added and labelled as LM_1,...,LM_j
#'   (length(func_lms) == j); one column for each interaction given in func_lms
#'   is added
#'
#' @examples
#' \dontrun{
#' data(relapse)
#' outcome = list(time="Time", status="event")
#' covars = list(fixed=c("ID","age.at.time.0","male","stage","bmi"),
#'               varying=c("treatment"))
#' w = 60; LMs = c(0,12,24)
#' # Covariate-landmark time interactions
#' func_covars <- list( function(t) t, function(t) t^2)
#' # let hazard depend on landmark time
#' func_lms <- list( function(t) t, function(t) t^2)
#' # Choose covariates that will have time interaction
#' pred_covars <- c("age","male","stage","bmi","treatment")
#' # Stack landmark datasets
#' lmdata <- stack_data(relapse, outcome, LMs, w, covars, format="long",
#'                      id="ID", rtime="T_txgiven", right=F)
#' # Update complex LM-varying covariates, note age is in years and LM is in months
#' lmdata$data$age <- lmdata$data$age.at.time.0 + lmdata$data$LM/12
#' # Add LM-time interactions
#' lmdata <- add_interactions(lmdata, pred_covars, func_covars, func_lms)
#' head(lmdata$data)
#' }
#'
#' @export
add_interactions <- function(lmdata, lm_covs, func_covars, func_lms, LM_col="LM",keep=T){
  if (LM_col %in% func_covars){
    stop(paste0("arg LM_col (given as ",LM_col,") should not be in arg func_covars."))
  }
  data <- lmdata$data
  if (missing(func_covars)){
    # f gives covariate-time interactions
    f1 <- function(t) t
    f2 <- function(t) t^2
    func_covars <- list(f1,f2)
  }
  if (missing(func_lms)){
    # g lets the hazard depend on time
    g1 <- function(t) f1(t)
    g2 <- function(t) f2(t)
    func_lms <- list(g1,g2)
  }

  allLMcovars <- c(lm_covs)
  data_LM <- data[[LM_col]]
  # Add func_covarss: covariate LM interactions
  for(i in 1:length(lm_covs)){
    for (j in 1:length(func_covars)){
      f <- func_covars[[j]]
      name <- paste0(lm_covs[i],"_",j)
      data[[name]]  <- data[[lm_covs[i]]]*f(data_LM)
      allLMcovars <- c(allLMcovars, name)
    }
  }
  # Add func_lms: LM interactions
  for (k in 1:length(func_lms)){
    g <- func_lms[[k]]
    name <- paste0("LM_",k)
    data[[name]]  <- g(data_LM)
    allLMcovars <- c(allLMcovars, name)
  }

  if(!keep){
    remaining = colnames(data)[! colnames(data)  %in% lm_covs]
    data <- data[remaining]
  }
  lmdata$data <- data

  lmdata$func_covars <- func_covars
  lmdata$func_lms <- func_lms
  lmdata$lm_covs <- lm_covs
  lmdata$allLMcovars <- unique(allLMcovars)
  lmdata$LM_col <- LM_col

  return(lmdata)
}


