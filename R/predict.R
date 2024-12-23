#' Calculate w-year risk from a landmark time point
#'
#' @param object Fitted landmark supermodel
#' @param newdata Either a dataframe of individuals to make predictions for or
#'  an object of class LMdataframe (e.g., created by calling [stack_data()] and
#'  [add_interactions()]). If it is a dataframe, it must contain the
#'  original covariates (i.e., without landmark interaction).
#' @param lms landmark time points that correspond to the entries in `newdata`.
#'  Only required when `newdata` is a data.frame.
#'  `lms` is either a time point, a vector or character string.
#'
#'  * For a single time point, w-year risk is predicted from this time for each
#'    data point.
#'  * For a vector, `lms` must have the same length as the number of rows of
#'    `newdata` (i.e., each data point is associated with one LM/prediction
#'    time point).
#'  * A character string indicates a column in `newdata`.
#'
#' @param cause Cause of interest for competing risks.
#' @param w Prediction window, i.e., predict w-year (/month/..) risk from each
#'   of the `lms`. Defaults to the w used in model fitting.
#'   If `w` > than that used in model fitting, results are unreliable, but can
#'   be produced by setting `extend = T`.
#' @param extend Argument to allow for predictions at landmark times that are
#'   later than those used in model fitting, or prediction windows greater
#'   than the one used in model fitting.
#'   Default is FALSE. If set to TRUE, predictions may be unreliable.
#' @param silence Silence the warning message when extend is set to TRUE.
#' @param complete Only make predictions for data entries with non-NA entries
#'   (i.e., non-NA predictions). Default is TRUE.
#' @param ... Unimplemented for now.
#'
#' @return An object of class "LMpred" with components:
#'   - preds: a dataframe with columns LM and risk, each entry corresponds to
#'     one individual and prediction time point (landmark)
#'   - w, type, LHS: as in the fitted super model
#'   - data: the newdata given in input
#'
#' @references van Houwelingen HC, Putter H (2012). Dynamic Prediction in
#'   Clinical Survival Analysis. Chapman & Hall.
#'
#' @examples
#' data(relapse)
#' outcome <- list(time = "Time", status = "event")
#' covars <- list(fixed = c("male", "stage", "bmi"),
#'                varying = c("treatment"))
#' w <- 60; lms <- c(0, 6, 12, 18)
#' lmdata <- stack_data(relapse, outcome, lms, w, covars, format = "long",
#'                      id = "ID", rtime = "T_txgiven")
#' lmdata <- add_interactions(lmdata, func_covars = c("linear", "quadratic"),
#'                            func_lms = c("linear", "quadratic"))
#'
#' formula <- "Hist(Time, event, LM) ~ male + male_LM1 + male_LM2 +
#'             stage + stage_LM1 + stage_LM2 + bmi + bmi_LM1 + bmi_LM2 +
#'             treatment + treatment_LM1 + treatment_LM2 + LM1 + LM2 + cluster(ID)"
#' supermodel <- dynamic_lm(lmdata, as.formula(formula), "CSC", x = TRUE)
#'
#' p1 <- predict(supermodel)
#' head(p1$preds)
#'
#' @import survival
#' @seealso [stack_data()], [add_interactions()], [dynamic_lm()], [score()],
#'   [calplot()]
#' @export
#'
predict.dynamicLM <- function(object, newdata, lms, cause, w, extend = FALSE,
                              silence = FALSE, complete = TRUE, ...) {
  func_covars <- object$func_covars
  func_lms <- object$func_lms
  model_w <- object$w
  if (missing(w)) {
    w <- model_w
  } else {
    if (w > model_w && !extend) {
      stop(tidymess(paste0(
        "Prediction window w (=", w, ") is larger than the window used in model
        fitting (=", model_w, "). If you wish to still make predictions at these
        times, set arg extend=T but note that results may be unreliable.")))
    } else if (w > model_w & extend) {
      if (!silence)
        message(tidymess(paste0(
          "NOTE: Prediction window w (=", w, ") is larger than the window used
          in model fitting (=", model_w, "). Predictions may be unreliable.")))
    }
  }
  fm <- object$model
  type <- object$type

  if (type == "coxph") {
    models <- list(fm)
    num_causes <- 1
    if (!missing(cause)) {
      if (!is.null(cause))
        stop("No cause should be specified for a coxph model.")
    }
    cause <- 1

  } else if (type == "CauseSpecificCox" || type == "CSC") {
    models <- fm$models
    num_causes <- length(models)

    if (missing(cause)) cause <- as.numeric(fm$theCause)
    else if (is.null(cause)) cause <- as.numeric(fm$theCause)

    if (length(cause) > 1)
      stop(paste0("Can only predict one cause. Provided are: ",
                  paste(cause, collapse = ", "), sep = ""))
    if (!(cause %in% fm$causes))
      stop("Error in cause. Not one of the causes in the model.")

  } else {
    stop("Error in fm argument. Input model is of the wrong type.
         \nOnly supported classes are CauseSpecificCox and coxph")
  }

  if (missing(newdata) & !missing(lms)) {
    stop("Argument newdata must be specified if lms are specified.")
  }

  if (!missing(newdata)) {
    if (inherits(newdata, "LMdataframe")) {
      lms <- newdata$lm_col
      newdata <- newdata$data
    } else if (missing(lms)) {
      lms <- newdata$LM
      if (is.null(lms)) {
        stop("Argument lms must be specified")
      }
    }

    if (inherits(lms, "character")) {
      lms <- newdata[[lms]]
      if (is.null(lms))
        stop("As lms is a string, it must be a column in newdata.")
    }

    ## Check prediction times match with LMs used in training
    if (max(lms) > object$end_time && !extend) {
      stop(tidymess(paste0(
        "Landmark/prediction time points lms contains values later than the last
        LM used in model fitting (last LM used in model fitting=",
        object$end_time, " and max lms value=", max(lms), "). If you wish to
        still make predictions at these times, set arg extend=T but note that
        results may be unreliable.")))
    } else if (max(lms) > object$end_time & extend) {
      if (!silence)
        message(tidymess(paste0(
          "NOTE:landmark/prediction time points lms contains values later (max
          value=", max(lms), ") than the last LM used in model fitting (=",
          object$end_time, "). Predictions at times after ", object$end_time,
          " may be unreliable.")))
    }
    ## Check prediction times & newdata given are coherent with each other
    num_preds <- nrow(newdata)
    if (!(length(lms) == num_preds)) {
      if (length(lms) == 1)
        lms <- rep(lms, num_preds)
      else
        stop(tidymess("Error in newdata or lms. Must have length(lms) ==
                      nrow(newdata) or lms be one landmarking point."))
    }

    ## Get risk scores
    all_combinations <- expand.grid(i = 1:num_preds, c = 1:num_causes)
    risks_values <- apply(all_combinations, 1, function(row) {
      i <- row['i']
      c <- row['c']
      riskScore(models[[c]], lms[i], newdata[i, ], func_covars, func_lms)
    })
    risks <- matrix(risks_values, nrow = num_causes, byrow = TRUE)
    data <- newdata

  } else {
    ## Get risk scores
    ## Note that linear predictors are centered, so need to un-center them for correct comparison.
    lms <- object$original.landmarks
    risks <- object$linear.predictors
    num_preds <- ncol(risks)
    # sanity check
    if (num_preds == 0) stop("Newdata must be specified.")
    if (length(lms) != num_preds)
      stop(tidymess("Error in newdata or lms. Must have length(lms) ==
                    nrow(newdata) or lms be one landmarking point."))
    data <- fm$call$data
  }

  # Baseline hazards
  sf <- lapply(1:num_causes, function(i) {
    base_data <- 0 * models[[i]]$coefficients
    survfit(models[[i]], newdata = base_data)
    })
  sf <- lapply(sf, function(s) {
    data.frame(time = s$time, surv = s$surv, Haz = -log(s$surv))
    })

  Fw <- rep(0, length(lms))

  sf1 <- sf[[cause]]
  Fw <- rep(NA, num_preds)
  # note: using sapply in any of these loops does not increase speed
  for (i in 1:num_preds) {
    if (is.na(risks[cause, i])) {
      Fw[i] <- NA
    } else {
      tLMi <- lms[i]
      pred_window <- (tLMi <= sf1$time & sf1$time <= tLMi + w)
      n_times <- sum(pred_window)

      # cause-specific instant hazard
      haz <- sf1$Haz[pred_window] * exp(risks[cause, i])
      instHaz <- haz[2:n_times] - haz[1:(n_times - 1)]
      idx <- (instHaz != 0)
      instHaz <- instHaz[idx]

      times <- sf1$time[pred_window][2:n_times]
      times <- times[idx]

      # overall survival
      w_adj <- times - tLMi
      surv <- c()
      for (j in seq_along(w_adj)) {
        s <- 0
        wj <- w_adj[j]
        for (c in 1:num_causes){
          sfc <- sf[[c]]
          sfc$Haz <- sfc$Haz * exp(risks[c, i])
          f <- stats::stepfun(sfc$time, c(0, sfc$Haz), right = FALSE)
          s <- s + (f(tLMi + wj) - f(tLMi))
        }
        surv <- c(surv, exp(-s))
      }

      # cumulative incidence from instant hazard and survival
      Fw[i] <- sum(instHaz * surv)
    }
  }

  if (complete) {
    idx <- !is.na(Fw)
    preds <- data.frame(LM = lms[idx], risk = Fw[idx])
    data <- data[idx, ]
  } else {
    preds <- data.frame(LM = lms, risk = Fw)
    data <- data
  }

  out <- list(
    preds = preds,
    w = w,
    type = type,
    LHS = object$LHS,
    data = data,
    cause = cause,
    outcome = object$outcome
  )
  class(out) <- "LMpred"
  return(out)
}
