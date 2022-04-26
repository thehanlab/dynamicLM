# -----------------------------------------------------------------------
#' Add landmarking time interations to a dataset
#'
#' @param LMdata An object of class "LMdataframe".
#' This can be created by running cutLMsuper, or creating a stacked data set and storing it in a list with attributes outcome, w and end_time
#' (see cutLMsuper for further description of outcome and w), end_time is the largest landmarking time.
#' @param LMcovars Vector of strings indicating the columns that are to have a LM interaction
#' @param func_covars A list of functions to use for interactions between LMs and covariates.
#' @param func_LMs A list of functions to use for transformations of the landmark times.
#' @param LM_col Character string specifying the column name that indicates the landmark time point for a row.
#' @param keep Boolean value to indicate whether or not to keep the columns given by LMcovars without the time interactions or not. Default=FALSE.
#'
#' @return An object of class "LMdataframe" which now also contains LM time-interactions.
#' For each variable "var" in LMcovars, new columns var_1,...,var_i (length(func_covars) == i) are added. One for each interaction given in func_covars
#' Transformations of the LM column are added and labelled as LM_1,...,LM_j (length(func_LMs) == j). One column for each interaction given in func_LMs
#' The object has the following components:
#' - w, outcome: as the input (obtained from LMdata)
#' - func_covars: as the input
#' - func_LMs: as the input
#' - LMcovars: as the input
#' - allLMcovars: a list of the new columns added
#' - LM_col: as the input
#' @export
addLMtime <- function(LMdata, LMcovars, func_covars, func_LMs, LM_col="LM",keep=F){
  if (LM_col %in% func_covars){
    stop(paste0("arg LM_col (given as ",LM_col,") should not be in arg func_covars."))
  }
  data <- LMdata$LMdata
  if (missing(func_covars)){
    # f gives covariate-time interactions
    f1 <- function(t) 1
    f2 <- function(t) t
    f3 <- function(t) t^2
    func_covars <- list(f1,f2,f3)
  }
  if (missing(func_LMs)){
    # g lets the hazard depend on time
    g1 <- function(t) f2(t)
    g2 <- function(t) f3(t)
    func_LMs <- list(g1,g2)
  }

  allLMcovars <- c()
  data_LM <- data[[LM_col]]
  # Add func_covarss: covariate LM interactions
  for(i in 1:length(LMcovars)){
    for (j in 1:length(func_covars)){
      f <- func_covars[[j]]
      name <- paste0(LMcovars[i],"_",j)
      data[[name]]  <- data[[LMcovars[i]]]*f(data_LM)
      allLMcovars <- c(allLMcovars, name)
    }
  }
  # Add func_LMs: LM interactions
  for (k in 1:length(func_LMs)){
    g <- func_LMs[[k]]
    name <- paste0("LM_",k)
    data[[name]]  <- g(data_LM)
    allLMcovars <- c(allLMcovars, name)
  }

  if(!keep){
    remaining = colnames(data)[! colnames(data)  %in% LMcovars]
    data <- data[remaining]
  }
  LMdata$LMdata <- data

  LMdata$func_covars <- func_covars
  LMdata$func_LMs <- func_LMs
  LMdata$LMcovars <- LMcovars
  LMdata$allLMcovars <- allLMcovars
  LMdata$LM_col <- LM_col

  return(LMdata)
}


