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
#' @param left.open Boolean (default = FALSE), indicating if the intervals for the
#'   time-varying covariates are open on the left (and closed on the right) or
#'   vice-versa.
#' @param split.data List of data split according to ID. Allows for faster
#'   computation.
#'
#' @details This function is based from [dynpred::cutLM()] with minor changes.
#'   The original function was authored by Hein Putter.
#' @references
#'   - van Houwelingen HC, Putter H (2012). Dynamic Prediction in
#'     Clinical Survival Analysis. Chapman & Hall.
#'   - The dynpred package
#'     (https://cran.r-project.org/web/packages/dynpred/index.html),
#'     in particular, the code for cutLM.
#' @return A landmark dataset.
#'
#' @examples
#' \dontrun{
#' data(relapse)
#' outcome <- list(time = "Time", status = "event")
#' covars <- list(fixed = c("age.at.time.0", "male", "stage", "bmi"),
#'                varying = c("treatment"))
#' lm12 <- get_lm_data(relapse, outcome, lm = 12, horizon = 60, covs = covars,
#'                     format = "long", id = "ID", rtime = "T_txgiven")
#' head(lm12)
#' }
#'
#' @seealso [dynamicLM::stack_data()]
#' @export
get_lm_data <- function(data, outcome, lm, horizon, covs,
                        format = c("wide", "long"), id, rtime,
                        left.open = FALSE, split.data) {
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
    lookup <- FALSE
    if (missing(split.data)) {
      data <- data[order(data[[id]], data[[rtime]]), ]
      ids <- unique(data[[id]])
      n <- length(ids)
      lookup <- TRUE
    } else {
      n <- length(split.data)
    }

    lmdata <- lapply(1:n, function(i) {
      if (lookup) {
        wh <- which(data[[id]] == ids[i])
        di <- data[wh, ]
      } else {
        di <- split.data[[i]]
      }
      t.fups <- di[[rtime]]

      idx <- findInterval(lm, c(t.fups, Inf), left.open = left.open)

      if (idx != 0) {
        return(di[idx, ])
      } else {
        out <- di[1, ]
        if (!is.null(covs$varying)) {
          out[, covs$varying] <- NA
          out[, rtime] <- NA
        }
        return(out)
      }
    })
    lmdata <- do.call(rbind, lmdata)
  }

  lmdata <- lmdata[lmdata[[outcome$time]] > lm, ]
  if (nrow(lmdata) == 0) {
    warning("Landmark dataset for lm = ", lm, " could not be constructed as no individuals are alive after this point.")
    return(NULL)
  }

  lmdata[outcome$status] <- lmdata[[outcome$status]] *
    as.numeric(lmdata[[outcome$time]] <= horizon)
  lmdata[outcome$time] <- pmin(as.vector(lmdata[[outcome$time]]), horizon)
  lmdata$LM <- lm
  if (format == "long")
    cols <- match(c(id, outcome$time, outcome$status, "LM", covs$fixed,
                    covs$varying, rtime), names(lmdata))
  else cols <- match(c(outcome$time, outcome$status, "LM", covs$fixed,
                       covs$varying), names(lmdata))
  return(lmdata[, cols])
}
