#' Build a landmark dataset
#'
#' @param data Data frame from which to construct landmark super dataset
#' @param outcome A list with items time and status, containing character
#'   strings identifying the names of time and status variables, respectively,
#'   of the survival outcome
#' @param lm The value of the landmark time point at which to construct the
#'   landmark dataset.
#' @param horizon Scalar, the value of the prediction window (ie predict risk
#'   within time w landmark points)
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
#' @details This function is based from [dynpred::cutLM()] with minor changes.
#'   The original function was authored by Hein Putter.
#' @references van Houwelingen HC, Putter H (2012). Dynamic Prediction in
#'   Clinical Survival Analysis. Chapman & Hall.
#' @return A landmark dataset.
#' @export
#'
# TODO: add examples
# TODO: add references
get_lm_data <- function(data, outcome, lm, horizon, covs,
                        format = c("wide", "long"), id, rtime, right = TRUE) {
  format <- match.arg(format)
  if (format == "wide") {
    lmdata <- data
    if (!is.null(covs$varying)) {
      for (col in covs$varying)
        lmdata[[col]] <- 1 - as.numeric(lmdata[[col]] > lm)
    }

  } else {
    if (missing(id))
      stop("argument 'id' should be specified for long format data")
    if (missing(rtime))
      stop("argument 'rtime' should be specified for long format data")
    ord <- order(data[[id]], data[[rtime]])
    data <- data[ord, ]
    ids <- unique(data[[id]])
    n <- length(ids)
    lmdata <- data[which(!duplicated(data[[id]])), ]
    for (i in 1:n) {
      wh <- which(data[[id]] == ids[i])
      di <- data[wh, ]
      idx <- cut(lm, c(di[[rtime]], Inf), right = right, labels = FALSE)
      if (!is.na(idx)) {
        lmdata[i, ] <- di[idx, ]
      } else {
        lmdata[i, ] <- di[1, ]
        if (!is.null(covs$varying)) {
          lmdata[i, covs$varying] <- NA
          lmdata[i, rtime] <- NA
        }
      }
    }
  }

  lmdata <- lmdata[lmdata[[outcome$time]] > lm, ]
  if (format == "long")
    lmdata <- lmdata[!is.na(lmdata[[id]]), ]
  lmdata[outcome$status] <- lmdata[[outcome$status]] *
    as.numeric(lmdata[[outcome$time]] <= horizon)
  lmdata[outcome$time] <- pmin(as.vector(lmdata[[outcome$time]]), horizon)
  lmdata$LM <- lm
  if (format == "long")
    cols <- match(c(id, outcome$time, outcome$status, covs$fixed,
                    covs$varying, rtime, "LM"), names(lmdata))
  else cols <- match(c(outcome$time, outcome$status, covs$fixed,
                       covs$varying, "LM"), names(lmdata))
  return(lmdata[, cols])
}
