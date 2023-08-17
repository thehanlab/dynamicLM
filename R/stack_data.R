#' Build a stacked dataset from original dataset (wide or long format).
#'
#' This stacked dataset output is used as input to [dynamic_lm()] to fit a
#' landmark supermodel for dynamic prediction. Calling [add_interactions()] on
#' the output before fitting the supermodel allows for landmark time
#' interactions to be included.
#'
#' @param data Data frame from which to construct landmark super dataset
#' @param outcome A list with items time and status, containing strings
#'   identifying the names of time and status variables, respectively, of the
#'   survival outcome
#' @param lms vector, the value of the landmark time points. This should be a
#'   range of points over the interval that prediction will be made. For
#'   example, if 5-year risk predictions are to be made over the first three
#'   years, this could be `c(0, 1.5, 3)`, `c(0, 1, 2, 3)` etc.
#' @param w Scalar, the value of the prediction window (ie predict risk within
#'   time w landmark points)
#' @param covs A list with items fixed and varying, containing character strings
#'   specifying column names in the data containing time-fixed and time-varying
#'   covariates, respectively.
#' @param format Character string specifying whether the original data are in
#'   wide (default) or in long format.
#' @param id Character string specifying the column name in data containing the
#'   subject id.
#' @param rtime Character string specifying the column name in data containing
#'   the (running) time variable associated with the time-varying variables;
#'   only needed if format = "long".
#' @param right Boolean (default = FALSE), indicating if the intervals for the
#'   time-varying covariates are closed on the right (and open on the left) or
#'   vice-versa.
#'
#' @return An object of class "LMdataframe". This the following components:
#'   - data: containing the stacked data set, i.e., the outcome and the values
#'     of time-fixed and time-varying covariates taken at the landmark time
#'     points. The value of the landmark time point is stored in column LM.
#'   - outcome: same as input
#'   - w: same as input
#'   - end_time: final landmarking point used in training
#'   - lm_col: "LM", identifies the landmark time column.
#'
#' @examples
#' \dontrun{
#' data(relapse)
#' outcome <- list(time = "Time", status = "event")
#' covars <- list(fixed = c("age.at.time.0", "male", "stage", "bmi"),
#'                varying = c("treatment"))
#' w <- 60; lms <- c(0, 6, 12, 18)
#' # Stack landmark datasets
#' lmdata <- stack_data(relapse, outcome, lms, w, covars, format = "long",
#'                      id = "ID", rtime = "T_txgiven")
#' head(lmdata$data)
#' }
#' @export
#'
stack_data <- function(data, outcome, lms, w, covs, format = c("wide", "long"),
                       id, rtime, right = FALSE) {
  if (!all(covs$fixed %in% colnames(data))) {
    stop(paste("Fixed column(s): ",
               paste(covs$fixed[!(covs$fixed %in% colnames(data))],
                     collapse = ","),
               "are not in the data."))
  }
  if (!is.null(covs$varying)){
    if (!all(covs$varying %in% colnames(data))) {
      stop(paste("Varying column(s): ",
                 paste(covs$varying[!(covs$varying %in% colnames(data))],
                       collapse = ","),
                 "are not in the data."))
    }
  }
  if (missing(id)) {
    if ("ID" %in% colnames(data))
      id <- "ID"
    else
      stop("argument id must be specified.")
  }
  if (!(id %in% colnames(data))) {
    stop(paste("ID column ", id, "is not in the data."))
  }

  if (format == "wide"){
    if (!(id %in% covs$fixed)) covs$fixed <- c(id, covs$fixed)

    lmdata <- get_lm_data(data = data,
                          outcome = outcome,
                          lm = lms[1],
                          horizon = lms[1] + w,
                          covs = covs,
                          format = "wide",
                          right = right)
    if (length(lms) > 1) {
      for (i in 2:length(lms))
        lmdata <- rbind(lmdata, get_lm_data(data = data,
                                            outcome = outcome,
                                            lm = lms[i],
                                            horizon = lms[i] + w,
                                            covs = covs,
                                            format = "wide",
                                            right = right))
    }
    lmdata <- lmdata[, c(which(colnames(lmdata) == id),
                        which(colnames(lmdata) != id))]

  } else if (format == "long") {
    # call cutLM
    lmdata <- get_lm_data(data = data,
                          outcome = outcome,
                          lm = lms[1],
                          horizon = lms[1] + w,
                          covs = covs,
                          format, id, rtime, right)
    if (length(lms) > 1) {
      for (i in 2:length(lms))
        lmdata <- rbind(lmdata, get_lm_data(data = data,
                                            outcome = outcome,
                                            lm = lms[i],
                                            horizon = lms[i] + w,
                                            covs = covs,
                                            format, id, rtime, right))
    }
  }
  out <- list(
    data = lmdata,
    outcome = outcome,
    w = w,
    end_time = lms[length(lms)],
    lm_col = "LM",
    id_col = id
  )
  class(out) <- "LMdataframe"
  return(out)
}
