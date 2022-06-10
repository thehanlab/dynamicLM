#' Build a stacked super dataset from original dataset (wide or long format)
#'
#' @param data Data frame from which to construct landmark super dataset
#' @param outcome List with items time and status, containing character strings identifying the names of time and status variables, respectively, of the survival outcome
#' @param LMs vector, the value of the landmark time points (points at which prediction is made)
#' @param w Scalar, the value of the prediction window (ie predict w-year/other time period risk from the LM points)
#' @param covs  List with items fixed and varying, containing character strings specifying column names in the data containing time-fixed and time-varying covariates, respectively
#' @param format Character string specifying whether the original data are in wide (default) or in long format
#' @param id Character string specifying the column name in data containing the subject id; only needed if format="long"
#' @param rtime Character string specifying the column name in data containing the (running) time variable associated; only needed if format="long"
#' @param right Boolean (default=TRUE), indicating if the intervals for the time-varying covariates are closed on the right (and open on the left) or vice versa, see cut
#'
#' @return An object of class "LM.data.frame". This the following components:
#'   - LMdata: containing the stacked data set, i.e., the outcome and the values of time-fixed and time-varying covariates taken at the landmark time points. The value of the landmark time point is stored in column LM.
#'   - outcome: same as input
#'   - w: same as input
#'   - end_time: final landmarking point used in training
#' @details This function calls cutLM from the library dynpred, more documentation can be found there. Note that for every landmark tLM given in LMs, there must be at least one patient alive after tLM.
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
#' # Choose covariates that will have time interaction
#' pred.covars <- c("age","male","stage","bmi","treatment")
#' # Stack landmark datasets
#' LMdata <- cutLMsuper(relapse, outcome, LMs, w, covs, format="long", id="ID", rtime="fup_time", right=F)
#' }
#' @export
#'
cutLMsuper <- function(data, outcome, LMs, w, covs, format = c("wide", "long"), id, rtime, right=T){
  if (format == "wide"){
    LMdata <- dynpred::cutLM(data=data,
                    outcome=outcome,
                    LM=LMs[1],
                    horizon=LMs[1]+w,
                    covs=covs,
                    format="wide",
                    right=right)
    if (length(LMs) > 1){
      for (i in 2:length(LMs))
        LMdata <- rbind(LMdata,dynpred::cutLM(data=data,
                                     outcome=outcome,
                                     LM=LMs[i],
                                     horizon=LMs[i]+w,
                                     covs=covs,
                                     format="wide",
                                     right=right))
    }

  } else if (format == "long"){
    LMdata <- dynpred::cutLM(data=data,
                    outcome=outcome,
                    LM=LMs[1],
                    horizon=LMs[1]+w,
                    covs=covs,
                    format, id, rtime, right)
    if (length(LMs) > 1){
      for (i in 2:length(LMs))
        LMdata <- rbind(LMdata,dynpred::cutLM(data=data,
                                     outcome=outcome,
                                     LM=LMs[i],
                                     horizon=LMs[i]+w,
                                     covs=covs,
                                     format, id, rtime, right))
    }
  }
  out=list(LMdata=LMdata, outcome=outcome, w=w, end_time=LMs[length(LMs)])
  class(out)="LMdataframe"
  return(out)
}
