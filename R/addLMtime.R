# -----------------------------------------------------------------------
# addLMtime : Add LM-time interations to a dataset that contains a column LM
# -----------------------------------------------------------------------
# Input:
# - LMdata           : An object of class "LM.data.frame", this can be created by running cutLMsuper, or creating a stacked data set and storing it in a list with attributes outcome and w
# - LMcovars         : List of covariates that are to have a LM interaction
# - func_covars      : A list of functions to use for interactions between LMs and covariates. If fitting a coxph model, the list has length 1, e.g. list(c(f1,f2,f3)). If fitting a CSC model, the list can have length 1 or length=number of causes, if different interactions are desired for different causes.
# - func_LMs         : A list of functions to use for transformations of the LMs. Its form is analogous to func_covars.
# - LM_col           : Character string specifying the column name of the LM value for a data point.
# -----------------------------------------------------------------------
# Output:
# LMdata : An object of class "LM.data.frame" which has the following components:
# - LMdata: The LM super dataset which now also contains LM time-interactions.
#           LM interactions with covariates are newly labelled as var_1,...,var_i if length(func_covars) == i and if var was a variable in LMcovars
#           LM effects are newly labelled as LM_1,...,LM_j if length(func_LMs) == j
# - w, outcome: (already)
# - func_covars: as the input
# - func_LMs: as the input
# - LMcovars: as the input
# - allLMcovars: a list of covariates that include LM-time interactions, i.e. if cov was in LMcovars, allLMcovars will contain cov_1, cov_2, ..., cov_i if there are i func_covars interactions
# - LM_col: as the input
# -----------------------------------------------------------------------
addLMtime <- function(LMdata, LMcovars, func_covars, func_LMs, LM_col="LM"){
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

  data <- data %>% select(-all_of(LMcovars))
  LMdata$LMdata <- data

  LMdata$func_covars <- func_covars
  LMdata$func_LMs <- func_LMs
  LMdata$LMcovars <- LMcovars
  LMdata$allLMcovars <- allLMcovars
  LMdata$LM_col <- LM_col

  return(LMdata)
}


head.LM.data.frame <- function(LMdata){
  print(head(LMdata$LMdata))
}
