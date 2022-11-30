#' Calcutes dynamic risk score at a time for an individual (helper to predLMrisk)
#'
#' @param object A coxph object
#' @param tLM Landmarking time point at which to calculate risk score (time at which the prediction is made)
#' @param data Dataframe (single row) of individual. Must contain the original covariates.
#' @param func_covars A list of functions to use for interactions between LMs and covariates.
#' @param func_LMs A list of functions to use for transformations of the landmark times.
#'
#' @return Numeric risk score
#' @export
#'
riskScore <- function(object, tLM, data, func_covars, func_LMs)
{
  not.null.LM <- !is.null(func_LMs)
  not.null.co <- !is.null(func_covars)

  coefs <- object$coefficients
  pred_covars <- names(coefs)
  idx_LM_covars <- grep("LM",pred_covars, fixed=TRUE)
  if (length(idx_LM_covars)>0){
    LM_covars <- pred_covars[idx_LM_covars]
    bet_covars <- pred_covars[-idx_LM_covars]
  } else {
    not.null.LM <- FALSE
    bet_covars <- pred_covars
  }

  # coef_LM1*g1(t) + coef_LM2*g2(t) + ...
  risk <- 0
  if (not.null.LM){
    risk <- risk + sum(
      sapply(LM_covars, function(coef_name){
        # Get associated function
        idx <-  as.numeric(sub(".*\\D+", "\\1", coef_name))
        return(func_LMs[[idx]](tLM) * coefs[coef_name])
      })
      # X1*coef + X1*t*coef + X1*t^2*coef + ..
    )
  }

  risk <- risk + sum(
    sapply(bet_covars, function(coef_name){
      # Get associated covariate info (remove _i from the name)
      covar <- sub("_\\d+$", "", coef_name)
      # Get associated function & multiply both by coef
      if (coef_name == covar) {
        return(coefs[coef_name] * data[,covar])
      }
      else {
        if (not.null.co){
          idx <- as.numeric(sub(".*\\D+", "\\1", coef_name))
          f <- func_covars[[idx]]
          return(f(tLM) * coefs[coef_name] * data[,covar])
        } else {
          return(0)
        }
      }

    })
  )
  return(risk)
}

