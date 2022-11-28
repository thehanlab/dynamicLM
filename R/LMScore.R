#' Methods (AUCt, Brier) to score the predictive performance of dynamic risk markers from LM super models
#'
#' @param preds A named list of prediction models, where allowed entries are outputs from predLMrisk
#' @param formula A survival or event history formula. The left hand side is used to compute the expected event status.
#'   If none is given, it is obtained from the prediction object.
#' @param metrics  Character vector specifying which metrics to apply. Choices are "auc" and "brier". Case matters.
#' @param cause Cause of interest if considering competing risks. If left blank, this is inferred from the preds object.
#' @param tLM  Landmark times for which scores must be given. These must be a subset of LM times used during the prediction.
#'    If left blank, all landmark prediction time points are considered.
#' @param unit Time unit for window of prediction, e.g., "year", "month", etc. Only used for printing results.
#' @param split.method Defines the internal validation design as in riskRegression::Score. Options are currently none or bootcv
#' @param B The number of bootstrap steps for cross-validation.
#' @param M The size of the subsamples drawn for bootstrap cross-validation.
#'   If specified it has to be an integer smaller than the size of data.
#' @param ... Additional arguments to pass to Score (riskRegression package).
#'   These arguments have been included for user flexibility but have not been tested and are not necessarily appropriate.
#' @return An object of class "LMScore", which has components:
#'   - auct: dataframe containing time-dependent auc information if "auc" was a metric
#'   - briert: dataframe containing time-dependent brier score if "brier" was a metric
#' @details See the Github for example code
#' @import riskRegression
#' @export
#'
# TODO: update description
LMScore <-
  function(object,
           times,
           metrics = c("auc", "brier"),
           formula,
           data,
           tLM,
           ID_col="ID",
           se.fit = TRUE,
           conf.int = 0.95,
           split.method = "none",
           B = 1,
           M,
           cores = 1,
           seed,
           unit ="year",
           cause,
           ...) {

    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package \"data.table\" must be installed to use function LMScore.",
           call. = FALSE)
    }

    checked_input <- match.call()
    m <- match(c("object", "times", "formula", "data", "tLM", "ID_col",
                 "split.method", "B", "M", "cores", "seed", "cause"), names(checked_input), 0L)
    checked_input <- checked_input[c(1L, m)]
    checked_input[[1L]] <- quote(check_evaluation_inputs)
    checked_input <- eval(checked_input, parent.frame())

    object = checked_input$object
    times = checked_input$times
    preds = checked_input$preds
    pred_LMs = checked_input$pred_LMs
    data = checked_input$data
    NF = checked_input$NF
    w = checked_input$w
    formula = checked_input$formula
    cause = checked_input$cause

    if(is.null(data$b)) data$b <- 1

    B = length(unique(data$b))
    se.fit.b <- se.fit
    if (B > 1) se.fit.b <- FALSE

    auct <- lapply(1:B, function(b) data.table::data.table())
    briert <- lapply(1:B, function(b) data.table::data.table())
    for (b in 1:B){
      for (t in 1:length(times)) {
        tLM = times[t]

        idx = (pred_LMs == tLM) & (data$b == b)
        data_to_test = data[idx,]
        risks_to_test = lapply(1:NF, function(i) {
          preds[idx, i]
        })
        names(risks_to_test) = names(object)

        # TODO: check for CR
        # --> not neccessary for cox
        #
        # if (object[[1]]$type == "coxph") {
        #   risks_to_test = lapply(risks_to_test, function(r) 1-r)
        # }

        if (nrow(data_to_test) != length(risks_to_test[[1]])) {
          stop("nrow(data_to_test)!=length(risks_to_test)")
        }

        # TODO: add a try around this
        # and use NAs if it doesn't work
        score_t <- riskRegression::Score(
          risks_to_test,
          formula = stats::as.formula(formula),
          data = data_to_test,
          metrics = metrics,
          cause = cause,
          times = c(tLM + w - 10e-5),
          se.fit = se.fit.b,
          conf.int = conf.int,
          ...
        )

        if ("auc" %in% metrics) {
          auct[[b]] <- rbind(auct[[b]], cbind(tLM, score_t$AUC$score))
        }

        if ("brier" %in% metrics) {
          briert[[b]] <- rbind(briert[[b]], cbind(tLM, score_t$Brier$score))
        }
      }
    }
    auct <- do.call("rbind", auct)
    briert <- do.call("rbind", briert)
    if (B > 1) {
      if (se.fit==TRUE) {
        alpha = 1-conf.int
        auct <- auct[,
                     data.table::data.table(
                       mean(.SD[["AUC"]],na.rm=TRUE),
                       se=sd(.SD[["AUC"]],na.rm=TRUE),
                       lower=quantile(.SD[["AUC"]],alpha/2,na.rm=TRUE),
                       upper=quantile(.SD[["AUC"]],(1-alpha/2),na.rm=TRUE)),
                     by=c("model","tLM"),.SDcols="AUC"
                     ]
        briert <- briert[,
                         data.table::data.table(
                           mean(.SD[["Brier"]],na.rm=TRUE),
                           se=sd(.SD[["Brier"]],na.rm=TRUE),
                           lower=quantile(.SD[["Brier"]],alpha/2,na.rm=TRUE),
                           upper=quantile(.SD[["Brier"]],(1-alpha/2),na.rm=TRUE)),
                         by=c("model","tLM"),.SDcols="Brier"
                         ]
        data.table::setnames(auct,c("model","tLM","AUC","se","lower","upper"))
        data.table::setnames(briert,c("model","tLM","Brier","se","lower","upper"))
      } else {
        auct <- auct[,data.table::data.table(mean(.SD[["AUC"]],na.rm=TRUE)),by=c("model","tLM"),.SDcols="AUC"]
        briert <- briert[,data.table::data.table(mean(.SD[["Brier"]],na.rm=TRUE)),by=c("model","tLM"),.SDcols="Brier"]
        data.table::setnames(auct,c("model","tLM","AUC"))
        data.table::setnames(briert,c("model","tLM","Brier"))
      }
    }

    outlist <- list(
      auct = auct,
      briert = briert,
      w = w,
      unit = unit
    )
    class(outlist) = "LMScore"
    return(outlist)
  }
