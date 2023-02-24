#' Build a stacked super dataset from original dataset (wide or long format)
#'
#' @param data Data frame from which to construct landmark super dataset
#' @param outcome List with items time and status, containing character strings identifying the names of time and status variables, respectively, of the survival outcome
#' @param lms vector, the value of the landmark time points (points at which prediction is made)
#' @param w Scalar, the value of the prediction window (ie predict w-year/other time period risk from the LM points)
#' @param covs  List with items fixed and varying, containing character strings specifying column names in the data containing time-fixed and time-varying covariates, respectively
#' @param format Character string specifying whether the original data are in wide (default) or in long format
#' @param id Character string specifying the column name in data containing the subject id; only needed if format="long"
#' @param rtime Character string specifying the column name in data containing the (running) time variable associated; only needed if format="long"
#' @param right Boolean (default=TRUE), indicating if the intervals for the time-varying covariates are closed on the right (and open on the left) or vice versa, see cut
#'
#' @return An object of class "LMdataframe". This the following components:
#'   - data: containing the stacked data set, i.e., the outcome and the values of time-fixed and time-varying covariates taken at the landmark time points. The value of the landmark time point is stored in column LM.
#'   - outcome: same as input
#'   - w: same as input
#'   - end_time: final landmarking point used in training
#' @details This function calls cutLM from the library dynpred, more documentation can be found there. Note that for every landmark tLM given in `lms`, there must be at least one patient alive after tLM.
#' @examples
#' \dontrun{
#' data(relapse)
#' outcome = list(time="Time", status="event")
#' covars = list(fixed=c("ID","age.at.time.0","male","stage","bmi"),
#'               varying=c("treatment"))
#' w = 60; lms = c(0,12,24)
#' # Stack landmark datasets
#' LMdata <- stack_data(relapse, outcome, lms, w, covs, format="long",
#'                      id="ID", rtime="T_txgiven", right=F)
#' head(LMdata$data)
#' }
#' @import dynpred
#' @export
#'
stack_data <- function(data, outcome, lms, w, covs, format = c("wide", "long"), id, rtime, right=T){
  if(!all(covs$fixed %in% colnames(data))){
    stop(paste("Fixed column(s): ",
               paste(covs$fixed[!(covs$fixed %in% colnames(data))], collapse=","),
               "are not in the data."))
  }
  if (!is.null(covs$varying)){
    if(!all(covs$varying %in% colnames(data))){
      stop(paste("Varying column(s): ",
                 paste(covs$varying[!(covs$varying %in% colnames(data))], collapse=","),
                 "are not in the data."))
    }
  }
  if (missing(id)){
    if ("ID" %in% colnames(data))
      id = "ID"
  }
  if (!(id %in% colnames(data))){
    stop(paste("ID column ", id,"is not in the data."))
  }

  if (format == "wide"){
    if (!(id %in% covs$fixed)){
      covs$fixed = c(id, covs$fixed)
    }
    LMdata = dynpred::cutLM(data=data,
                    outcome=outcome,
                    LM=lms[1],
                    horizon=lms[1]+w,
                    covs=covs,
                    format="wide",
                    right=right)
    if (length(lms) > 1){
      for (i in 2:length(lms))
        LMdata = rbind(LMdata,dynpred::cutLM(data=data,
                                     outcome=outcome,
                                     LM=lms[i],
                                     horizon=lms[i]+w,
                                     covs=covs,
                                     format="wide",
                                     right=right))
    }
    LMdata = LMdata[,c(which(colnames(LMdata)==id),
                       which(colnames(LMdata)!=id))]

  } else if (format == "long"){
    # guard against errors in cutLM
    if (is.null(covs$varying)) {
      covs$varying = "fake_column1234"
      data["fake_column1234"] = NA
    }

    # call cutLM
    LMdata = dynpred::cutLM(data=data,
                    outcome=outcome,
                    LM=lms[1],
                    horizon=lms[1]+w,
                    covs=covs,
                    format, id, rtime, right)
    if (length(lms) > 1){
      for (i in 2:length(lms))
        LMdata = rbind(LMdata,dynpred::cutLM(data=data,
                                     outcome=outcome,
                                     LM=lms[i],
                                     horizon=lms[i]+w,
                                     covs=covs,
                                     format, id, rtime, right))
    }
  }
  out=list(data=LMdata, outcome=outcome, w=w, end_time=lms[length(lms)])
  class(out)="LMdataframe"
  return(out)
}
