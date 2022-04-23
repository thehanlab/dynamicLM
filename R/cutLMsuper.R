# -----------------------------------------------------------------------
# cutLMsuper: Build stacked super dataset from original dataset
#              (which can be either in wide or long format - see details in cutLM)
# -----------------------------------------------------------------------
# Input:
# - data             : Data frame from which to construct landmark super dataset
# - outcome          : List with items time and status, containing character strings
#                      identifying the names of time and status variables, respectively,
#                      of the survival outcome
# - LMs              : vector, the value of the landmark time points (points at which prediction is made)
# - w                : Scalar, the value of the prediction window (ie predict w-year/other time period risk from the LM points)
# - covs             : List with items fixed and varying, containing character strings
#                      specifying column names in the data containing time-fixed and
#                      time-varying covariates, respectively
# - format	         : Character string specifying whether the original data are in wide (default) or in long format
# - id	             : Character string specifying the column name in data containing the subject id; only needed if format="long"
# - rtime	           : Character string specifying the column name in data containing the (running) time variable associated; only needed if format="long"
# - right	           : Boolean (default=TRUE), indicating if the intervals for the time-varying covariates are closed on the
#                      right (and open on the left) or vice versa, see cut
# -----------------------------------------------------------------------
# Output:
# LMdata             : An object of class "LM.data.frame". This has various attributes
# - $LMdata: containing the stacked data set, i.e., the outcome and the values of time-fixed and time-varying covariates taken at the landmark time points. The value of the landmark time point is stored in column LM.
# - $outcome: same as input
# - $w: same as input
# - $end_time: final landmarking point used in training
# -----------------------------------------------------------------------
cutLMsuper <- function(data, outcome, LMs, w, covs, format = c("wide", "long"), id, rtime, right=T){
  if (format == "wide"){
    LMdata <- cutLM(data=data,
                    outcome=outcome,
                    LM=LMs[1],
                    horizon=LMs[1]+w,
                    covs=covs,
                    format="wide",
                    right=right)
    if (length(LMs) > 1){
      for (i in 2:length(LMs))
        LMdata <- rbind(LMdata,cutLM(data=data,
                                     outcome=outcome,
                                     LM=LMs[i],
                                     horizon=LMs[i]+w,
                                     covs=covs,
                                     format="wide",
                                     right=right))
    }

  } else if (format == "long"){
    LMdata <- cutLM(data=data,
                    outcome=outcome,
                    LM=LMs[1],
                    horizon=LMs[1]+w,
                    covs=covs,
                    format, id, rtime, right)
    if (length(LMs) > 1){
      for (i in 2:length(LMs))
        LMdata <- rbind(LMdata,cutLM(data=data,
                                     outcome=outcome,
                                     LM=LMs[i],
                                     horizon=LMs[i]+w,
                                     covs=covs,
                                     format, id, rtime, right))
    }
  }
  out=list(LMdata=LMdata, outcome=outcome, w=w, end_time=LMs[length(LMs)])
  class(out)="LM.data.frame"
  return(out)
}
