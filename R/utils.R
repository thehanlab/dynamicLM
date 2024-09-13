# ----------------------------------------------------------
# tidymess: Helper function to make a string more readable by
# wrapping it
# Reference: https://stackoverflow.com/questions/45693010/how-do-you-format-multiline-r-package-messages
# ----------------------------------------------------------
tidymess <- function(..., prefix = "\n", initial = "") {
  strwrap(..., prefix = prefix, initial = initial)
}


# ----------------------------------------------------------
# find_se_log: Helper function to calculate SE for time varying log HR
# of the form coef[1] + t*coef[2] + t^2*coef[2] + ..
# ----------------------------------------------------------
find_se_log <- function(t, coefs, covar, func_covars) {
  if (!requireNamespace("msm", quietly = TRUE)) {
    stop("Package \"msm\" must be installed to use function find_se_log",
         call. = FALSE)
  }
  form <- "x1"
  if (length(coefs) > 1) {
    for (i in 2:length(coefs)){
      form <- paste(form, sprintf("%s * %f", paste("x", i, sep = ""),
                    func_covars[[i - 1]](t)), sep = " + ")
    }
  }
  form <- paste("~", form)
  se <- msm::deltamethod(stats::as.formula(form), coefs, covar)
  return(se)
}


# ----------------------------------------------------------
# find_se: Helper function to calculate SE for time varying HR
# of the form coef[1] + t*coef[2] + t^2*coef[2] + ..
# ----------------------------------------------------------
find_se <- function(t, coefs, covar, func_covars) {
  if (!requireNamespace("msm", quietly = TRUE)) {
    stop("Package \"msm\" must be installed to use function find_se",
         call. = FALSE)
  }
  form <- "x1"
  if (length(coefs) > 1) {
    for (i in 2:length(coefs)){
      form <- paste(form, sprintf("%s * %f", paste("x",i,sep = ""),
                    func_covars[[i - 1]](t)), sep = " + ")
    }
  }
  form <- paste("~ exp(", form,")")
  se <- msm::deltamethod(stats::as.formula(form), coefs, covar)
  return(se)
}


# ----------------------------------------------------------
# replace_na_with_last: Helper function to update a vector by
# replacing any NA value with the previous non-NA value
# ----------------------------------------------------------
replace_na_with_last <- function(x) {
  c(NA, x)[cummax(seq_along(x)*(!is.na(x)))+1]
}


# ----------------------------------------------------------
# getLHS: Function that takes a formula as an input and obtains
# the LHS, and converts it to a Hist object if it is a Surv object
# ----------------------------------------------------------
getLHS <- function(formula) {
  lhs <- formula[[2]]
  survival_object <- lhs[[1]]

  if (length(lhs) == 3) {
    exit <- lhs[[2]]
    status <- lhs[[3]]
  } else if (survival_object == "Surv") {
    enter <- lhs[[2]]
    exit <- lhs[[3]]
    status <- lhs[[4]]
  } else if (survival_object == "Hist") {
    exit <- lhs[[2]]
    status <- lhs[[3]]
    enter <- lhs[[4]]
  } else {
    stop("Invalid formula: LHS must be an object of type Hist or Surv")
  }

  exit <- deparse(exit)
  status <- deparse(status)

  if (length(lhs) == 3) return(paste0("Hist(", exit, ", ", status, ") ~ 1"))

  enter <- deparse(enter)
  return(paste0("Hist(", exit, ", ", status, ", ", enter, ") ~ 1"))
}


# ----------------------------------------------------------
# update_hist_formula: Function that adjusts a formula's LHS based on a
# specified type ("censor" or "surv") to generate a corresponding Surv formula.
# ----------------------------------------------------------
update_hist_formula <- function(formula, type) {
  lhs <- formula[[2]]
  if (length(lhs) == 3){
    exit <- lhs[[2]]
    status <- lhs[[3]]
  } else {
    exit <- lhs[[2]]
    status <- lhs[[3]]
    enter <- lhs[[4]]
  }

  exit <- deparse(exit)
  status <- deparse(status)

  # if (type == "censor")  outcome <- " == 1"
  if (type == "censor")  outcome <- " == 0"
  else if (type == "surv") outcome <- " != 0"
  else stop("updateHist can only handle censor and surv types.")

  # if (length(LHS) == 3)
  #   out <- paste0("Surv(", exit, ", ", status, outcome, ") ~ 1")
  # else
  #   out <- paste0("Surv(", enter, ", ", exit, ", ", status, outcome, ") ~ 1")
  if (length(lhs) == 3)
    out <- paste0("Hist(", exit, ", ", status, outcome, ") ~ 1")
  else
    out <- paste0("Hist(", exit, ", ", status, outcome,  ", ", enter, ") ~ 1")

  return(out)
}


# ----------------------------------------------------------
# clean_bootstraps: given the outcome from multiple bootstrap replicates,
# clean the table for user readability & to obtain results.
# ----------------------------------------------------------
clean_bootstraps <- function(table, column, alpha, contrasts = FALSE,
                             se.fit = TRUE, summary = FALSE) {
  by_columns <- c("tLM", "model")
  if (summary) by_columns <- c("model")
  if (contrasts) by_columns <- c(by_columns, "reference")

  if (se.fit) {
    out <- table[, data.table::data.table(
      mean(.SD[[column]], na.rm = TRUE),
      se = stats::sd(.SD[[column]], na.rm = TRUE),
      lower = stats::quantile(.SD[[column]], alpha / 2, na.rm = TRUE),
      upper = stats::quantile(.SD[[column]], (1 - alpha / 2), na.rm = TRUE)
    ), by = by_columns, .SDcols = column]
    data.table::setnames(
      out, c(by_columns, column, "se", "lower", "upper"))
    if (contrasts) {
      out[, p := 2 * stats::pnorm(abs(get(column) / se),
                                  lower.tail=FALSE)]
    }
  } else {
    out <- table[, data.table::data.table(
      mean(.SD[[column]], na.rm = TRUE)
    ), by = by_columns, .SDcols = column]
    data.table::setnames(out, c(by_columns, column))
  }
  return(out)
}


# ----------------------------------------------------------
# initialize_df: initializes a data frame with the specified column names and
# one row filled with NA.
# ----------------------------------------------------------
initialize_df <- function(metric, bootstrap) {
  if (metric == "AUC")
    cols <- c("tLM", "model", "times", "AUC", "se", "lower", "upper")
  if (metric == "Brier")
    cols <- c("tLM", "model", "times", "Brier", "se", "lower", "upper")
  if (metric == "delta.AUC")
    cols <- c("tLM", "times", "model", "reference", "delta.AUC", "se",
              "lower", "upper", "p")
  if (metric == "delta.Brier")
    cols <- c("tLM", "times", "model", "reference", "delta.Brier", "se",
              "lower", "upper", "p")
  if (metric == "IF.AUC")
    cols <- c("tLM", "riskRegression_ID", "model", "cause", "times", "IF.AUC")
  if (metric == "IF.Brier")
    cols <- c("tLM", "riskRegression_ID", "model", "cause", "times", "IF.Brier")

  df <- data.frame(matrix(NA, nrow = 1, ncol = length(cols)))
  colnames(df) <- cols
  df$bootstrap <- bootstrap
  df
}


#' Altered code from `riskRegression` of the cause-specific Cox model
#' to fit a CSC model with given coefficients.
#'
#' @param formula Formula to fit the model
#' @param data Data on which to which
#' @param cause Main cause of interest
#' @param cause.specific.coefs Coefficients that each model should be fit with
#' @param ... Additional arguments to coxph.
#'
#' @return CSC model
#' @references
#'   - `riskRegression` package:
#'     <https://cran.r-project.org/web/packages/riskRegression/index.html>
CSC.fixed.coefs <- function(formula, data, cause,
                            cause.specific.coefs, ...) {
  fitter <- "coxph"
  surv.type <- "hazard"

  # {{{ formulae & response
  if (inherits(x = formula, what = "formula")) formula <- list(formula)
  call <- match.call()
  # get outcome information from formula
  Rform <- stats::update(formula[[1]], ".~1")
  response <- eval(Rform[[2]], envir = data)
  if (any(is.na(response)))
    stop("Event history response may not contain missing values")
  time <- response[, "time"]
  status <- response[, "status"]
  event <- prodlim::getEvent(response)
  if ("entry" %in% colnames(response)) {
    entry <- response[, "entry"]
  } else {
    entry <- NULL
  }
  if (any(entry > time)) stop(tidymess("entry > time detected. Entry time into
      the study must be strictly greater than outcome time."))
  ## remove event history variables from data
  if (any((this <- match(all.vars(Rform), names(data), nomatch = 0)) > 0)) {
    if (data.table::is.data.table(data))
      data <- data[, - this, with = FALSE]
    else
      data <- data[, -this]
  }
  # }}}
  # {{{ sorted unique event times
  eventTimes <- unique(sort(as.numeric(time[status != 0])))
  # }}}
  # {{{ causes
  causes <- prodlim::getStates(response)
  NC <- length(causes)
  if (length(cause.specific.coefs) != NC) stop(tidymess(
    "There should be one entry for each cause specific argument but there are ",
     NC, " causes and ", length(cause.specific.coefs),
     " elements in cause.specific.coefs"))

  if (length(formula) != NC[1] && length(formula)>1) stop(tidymess(
    "Wrong number of formulae. Should be one for each cause ", NC, "."))
  if (length(formula) == 1) {
    formula <- lapply(1:NC, function(x) formula[[1]])
  }
  # }}}
  # {{{ find the cause of interest
  if (missing(cause)) {
    theCause <- causes[1]
  } else {
    if ((foundCause <- match(as.character(cause),causes,nomatch=0)) == 0) {
      stop(paste0("Cannot find all requested cause(s) ...\n\n",
                  "Requested cause(s): ", paste0(cause, collapse = ", "),
                  "\n Available causes: ", paste(causes, collapse = ", "),
                  "\n"))
    } else {
      theCause <- causes[foundCause]
    }
  }
  otherCauses <- causes[-match(theCause, causes)]
  # }}}
  # {{{ fit Cox models
  CoxModels <- lapply(1:NC, function(x) {
    if (x == 1)
      causeX <- theCause
    else
      causeX <- otherCauses[x - 1]

    statusX <- as.numeric(event == causeX)

    if (is.null(entry))
      workData <- data.frame(time = time, status = statusX)
    else
      workData <- data.frame(time = time, status = statusX, entry = entry)
    if (any(this <- match(names(data), names(workData), nomatch = 0) > 0)) {
      warning(paste("Variables named",
                    paste(names(data)[this], collapse = ", "),
                    "in data will be ignored."))
      if (data.table::is.data.table(data))
        data <- data[, -this, with = FALSE]
      else
        data <- data[, -this, drop = FALSE]
    }
    workData <- cbind(workData, data)
    if (is.null(entry))
      survresponse <- "survival::Surv(time, status)"
    else
      survresponse <- "survival::Surv(entry, time, status)"
    ## check whether right hand side of formula includes ~.
    allvars <- all.vars(formula[[x]])
    if (any(grepl("^\\.$", allvars))) {
      formulaXX <- stats::as.formula(paste0(survresponse, "~."))
    } else {
      formulaXX <- stats::update(formula[[x]], paste0(survresponse, "~."))
    }

    args <- list(formulaXX, data = workData)
    extra.args <- list(...)
    fit <- do.call("coxph",
                   c(args, list(x = TRUE, y = TRUE, iter.max = 0,
                                init = cause.specific.coefs[[x]]), extra.args))
    fit$call$formula <- formulaXX
    fit$call$data <- workData
    fit
  })
  names(CoxModels) <- paste("Cause", c(theCause, otherCauses))

  # }}}
  out <- list(call = call,
              models = CoxModels,
              response = response,
              eventTimes = eventTimes,
              surv.type = surv.type,
              fitter = fitter,
              theCause = theCause,
              causes = c(theCause, otherCauses))
  class(out) <- "CauseSpecificCox"
  out
}
