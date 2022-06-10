#' Calculate w-year risk from a landmark time point
#'
#' @param supermodel fitted landmarking supermodel
#' @param newdata dataframe of individuals to make predictions for. Must contain the original covariates (i.e., without landmark interaction)
#' @param tLM time points at which to predict risk of w more years
#' Note tLM must be one value for newdata or must have the same length as the number of rows of newdata
#' (i.e., each datapoint is associated with one LM/prediction time point)
#' @param cause Cause of interest if under competing risks
#' @param extend Argument to allow for predictions at landmark times that are greater than those used in model fitting.
#' Default is FALSE. If set to TRUE, predictions may be unreliable
#' @param silence Silence the warning message when extend is set to TRUE
#' @param complete Only make predictions for data entries with non-NA entries (i.e., non-NA predictions)
#'
#' @return An object of class "LMpred" with components:
#' - preds: a dataframe with columns LM and risk, each entry corresponds to one individual and prediction time point (landmark)
#' - w, type, LHS: as in the fitted super model
#'- data: the newdata given in input
#' @details See the Github for example code
#' @import survival
#' @export
#'
predLMrisk <- function(supermodel, newdata, tLM, cause, extend=F, silence=F, complete=T)
{
  # TODO: allow for different prediction window than w (with a warning)
  # TODO: check supermodel/args are correct = that newdata contains correct covars => check that names of newdata incl. in LMcovars
  # TODO: allow for inputting a fitted model with additional arguments (like fitLM)

  func_covars <- supermodel$func_covars
  func_LM <- supermodel$func_LM
  w <- supermodel$w
  fm <- supermodel$model
  type <- supermodel$type

  if (type == "coxph") {
    models <- list(fm)
    num_causes <- 1
    if (!missing(cause) ){
      if(!is.null(cause)) stop("No cause should be specified for a coxph model.")
    }
    cause<-1

  } else if (type == "CauseSpecificCox" | type =="CSC") {
    models <- fm$models
    num_causes <- length(models)
    if (missing(cause)) { cause <- as.numeric(fm$theCause)
    } else if (is.null(cause)) { cause <- as.numeric(fm$theCause) }
    if (length(cause) > 1) stop(paste0("Can only predict one cause. Provided are: ", paste(cause, collapse = ", "), sep = ""))
    if (!(cause %in% fm$causes)) stop("Error in cause. Not one of the causes in the model.")

  } else {
    stop("Error in fm argument. Input model is of the wrong type.
         \nOnly supported classes are CauseSpecificCox and coxph")
  }

  if (missing(newdata) & !missing(tLM)) {
    stop("newdata must be specified")
  }

  if (!missing(newdata)){
    if (missing(tLM)) {
      tLM = newdata$LM
      if (is.null(tLM)) {
        stop("tLM must be specified")
      }
    }

    ## Check prediction times match with LMs used in training
    if (max(tLM) > supermodel$end_time & !extend){
      stop(paste0("Landmark/prediction time points tLM contains values later than the last LM used in model fitting
                (last LM used in model fitting=",supermodel$end_time," and max tLM value=",max(tLM),").
                If you wish to still make predictions at these times, set arg extend=T but note that results may be unreliable."))
    }
    else if (max(tLM) > supermodel$end_time & extend){
      if (!silence) message(paste0("NOTE:landmark/prediction time points tLM contains values later (max value=",max(tLM),") than the last LM used in model fitting (=",supermodel$end_time,").",
                                   "\nPredictions at times after ",supermodel$end_time," may be unreliable."))
    }
    ## Check prediction times & newdata given are coherent with each other
    num_preds <- nrow(newdata)
    if(!(length(tLM)==num_preds)){
      if (length(tLM) == 1){
        tLM<-rep(tLM,num_preds)
      } else {
        stop("Error in newdata or tLM. Must have length(tLM) == nrow(newdata) or tLM be one landmarking point.")
      }
    }
    ## Check newdata contains the correct covariates
    ## TODO

    ## Get risk scores
    risks <- matrix(sapply(1:num_preds, function(i){
      tLMi <- tLM[i]
      newdatai <- newdata[i,]
      sapply(1:num_causes, function(c) riskScore(models[[c]], tLMi, newdatai, func_covars, func_LM))
    }),nrow=num_causes)
    data <- newdata

  } else {
    ## Get risk scores
    ## Note that linear predictors are centered, so need to un-center them for correct comparison.
    tLM <- supermodel$original.landmarks
    risks <- supermodel$linear.predictors
    num_preds <- ncol(risks)
    # sanity check
    if (num_preds == 0){stop("Newdata must be specified.")}
    if (length(tLM) != num_preds){ stop("Error in newdata or tLM. Must have length(tLM) == nrow(newdata) or tLM be one landmarking point.") }
    data <- fm$call$data
  }

  # Baseline hazards
  sf <- lapply(1:num_causes,function(i) {
    base_data = 0 * models[[i]]$coefficients
    survfit(models[[i]], newdata=base_data)
    })
  sf <- lapply(sf, function(s) data.frame(time=s$time,surv=s$surv,Haz=-log(s$surv)))

  Fw <- rep(0, length(tLM))

  sf1 <-sf[[cause]]
  Fw <- rep(NA, num_preds)
  for(i in 1:num_preds) {
    if (is.na(risks[cause,i])){
      Fw[i] <- NA
    } else {
      tLMi <- tLM[i]
      pred_window <- (tLMi <= sf1$time & sf1$time <= tLMi+w)
      n_times <- sum(pred_window)

      haz <- sf1$Haz[pred_window] * exp(risks[cause,i])
      instHaz <- haz[2:n_times]-haz[1:n_times-1]
      idx <- (instHaz != 0)
      instHaz <- instHaz[idx]

      times <- sf1$time[pred_window][2:n_times]
      times <- times[idx]

      w_adj <- times-tLMi
      surv <- c()
      for (j in 1:length(w_adj)){
        s <- 0
        wj<-w_adj[j]
        for (c in 1:num_causes){
          sfc <- sf[[c]]
          sfc$Haz <- sfc$Haz * exp(risks[c,i])
          f <- stats::stepfun(sfc$time, c(0,sfc$Haz), right=FALSE)
          s <- s + (f(tLMi+wj)-f(tLMi))
        }
        surv <- c(surv,exp(-s))
      }
      Fw[i] <- sum(instHaz*surv)
    }
  }

  if(complete){
    idx = !is.na(Fw)
    preds = data.frame(LM=tLM[idx],risk=Fw[idx])
    data = data[idx,]
  } else {
    preds = data.frame(LM=tLM,risk=Fw)
    data = data
  }

  out = list(preds=preds, w=w, type=type, LHS=supermodel$LHS, data=data, cause=cause)
  class(out) = "LMpred"
  return(out)
}
