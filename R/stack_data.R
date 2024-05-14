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
#'   covariates, respectively. If missing, all columns that are not the outcome,
#'   rtime or id column are set to be time-varying covariates.
#' @param format Character string specifying whether the original data are in
#'   wide (default) or in long format.
#' @param id Character string specifying the column name in data containing the
#'   subject id.
#' @param rtime Character string specifying the column name in data containing
#'   the (running) time variable associated with the time-varying variables;
#'   only needed if format = "long".
#' @param left.open Boolean (default = FALSE), indicating if the intervals for
#'   the time-varying covariates are open on the left (and closed on the right)
#'   or vice-versa.
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
#'
#' @seealso [dynamicLM::add_interactions()], [dynamicLM::dynamic_lm()]
#' @export
#'
stack_data <- function(data, outcome, lms, w, covs, format = c("wide", "long"),
                       id, rtime, left.open = FALSE) {
  ### Check input ###

  # Check all the columns are in the data
  if (!all(covs$fixed %in% colnames(data))) {
    stop(paste("Fixed column(s): ",
               paste(covs$fixed[!(covs$fixed %in% colnames(data))],
                     collapse = ", "),
               "are not in the data."))
  }
  if (!is.null(covs$varying)){
    if (!all(covs$varying %in% colnames(data))) {
      stop(paste("Varying column(s): ",
                 paste(covs$varying[!(covs$varying %in% colnames(data))],
                       collapse = ", "),
                 "are not in the data."))
    }
  }
  if (missing(id)) {
    if ("ID" %in% colnames(data))
      id <- "ID"
    else
      stop("argument id must be specified.")
  }
  if (!(id %in% colnames(data)))
    stop(paste("ID column ", id, "is not in the data."))
  if (!(rtime %in% colnames(data)))
    stop(paste("rtime column ", rtime, "is not in the data."))

  # Make sure to not duplicate columns
  fixed <- setdiff(covs$fixed, c(outcome$time, outcome$status))
  varying <- setdiff(covs$varying, c(outcome$time, outcome$status))

  # Default covs$varying to all covariates if the covs argument is not specified
  if (is.null(fixed) && is.null(varying)) {
    varying <- setdiff(colnames(data), c(id, outcome$time, outcome$status))
    if (!missing(rtime)) {
      varying <- setdiff(varying, rtime)
    }
  }
  covs$fixed <- fixed
  covs$varying <- varying
  all_covs <- c(covs$fixed, covs$varying)

  # Make sure that (ID, rtime) is never duplicated
  if (sum(duplicated(data[, c(id, rtime)])) > 0) {
    stop(tidymess("There are multiple entries for some patients (i.e., id
                  column) at some (running) time points (i.e., rtime column).
                  There can only be one entry for each patient at each time
                  point."))
  }

  # Check there are no variables in the form {name}_{integer}
  matches_pattern <- function(x) grepl("^[^_]+_[0-9]+$", x)
  matching_strings <- all_covs[matches_pattern(all_covs)]
  if (length(matching_strings) > 0) {
    stop(tidymess(paste(
      "Covariates should not be named in the form {name}_{integer} as the
      dynamicLM library reserves this naming convention for time interactions.
      Please rename the following columns:",
      paste(matching_strings, collapse = ", "))))
  }

  # Check that there are no factor variables in the data
  factors <- sapply(data[, all_covs], function(v) inherits(v, "factor"))
  chars <- sapply(data[, all_covs], function(v) inherits(v, "character"))
  if (any(factors)) {
    stop(tidymess(paste(
      "No covariates can be factors. Please convert the following column(s) to
      dummy variables:", paste(all_covs[factors], collapse = ", "), "")))
  }
  if (any(chars)) {
    stop(tidymess(paste(
      "No covariates can be factors or characters. Please convert the following
      column(s) to dummy or numeric variables:",
      paste(all_covs[chars], collapse = ","), "")))
  }


  ### Stack data ###

  if (format == "wide"){
    if (!(id %in% covs$fixed)) covs$fixed <- c(id, covs$fixed)

    lmdata <- lapply(lms, function(lm) {
      get_lm_data(data = data, outcome = outcome, lm = lm, horizon = lm + w,
                  covs = covs, format = "wide", left.open = left.open)
      })
    lmdata <- do.call(rbind, lmdata)

  } else if (format == "long") {
    if (id %in% covs$fixed) covs$fixed <- setdiff(covs$fixed, id)

    data <- data[order(data[[id]], data[[rtime]]), ]
    split.data <- split(data, data[[id]])
    # TODO: can still improve here
    # For example use rle (faster) and then turn into a list
    # run_lengths <- rle(data[[id]])$lengths

    lmdata <- lapply(lms, function(lm) {
      get_lm_data(data = data, outcome = outcome, lm = lm, horizon = lm + w,
                  covs = covs, format = "long", id = id, rtime = rtime,
                  left.open = left.open, split.data = split.data)
    })
    lmdata <- do.call(rbind, lmdata)

  }
  out <- list(
    data = lmdata,
    outcome = outcome,
    w = w,
    end_time = lms[length(lms)],
    lm_col = "LM",
    id_col = id,
    all_covs = all_covs
  )
  class(out) <- "LMdataframe"
  return(out)
}
