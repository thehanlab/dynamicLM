#' Methods (time-dependent AUC and Brier Score) to score the predictive
#' performance of dynamic risk prediction landmark models.
#'
#' There are three ways to perform assess the predictive performance:
#' apparent/internal, bootstrapped, and external. Accordingly, the named list of
#' prediction models must be as follows:
#' * For both apparent/internal evaluation, objects output from `predLMrisk` or
#'   supermodels fit with `fitLM` may be used as input.
#' * In order to bootstrap, supermodels fit with `fitLM` may be used as input
#'   (note that the argument `x=TRUE` must be specified when fitting the model
#'   in `fitLM`).
#' * For external calibration, supermodels fit with `fitLM` are input along with
#'   new data in the `data` argument. This data can be a LMdataframe or a
#'   dataframe (in which case `tLM` must be specified).
#'
#' For both internal evaluation and bootstrapping, it is assumed that all
#' models in `object` are fit on the same data.
#'
#'
#' @param object A named list of prediction models, where allowed entries are
#'   outputs from `predLMrisk` or supermodels from `fitLM` depending on the type
#'   of calibration.
#' @param times Landmark times for which calibration must be plot. These must be
#'   a subset of LM times used during the prediction
#' @param metrics  Character vector specifying which metrics to apply. Choices
#'   are "auc" and "brier". Case matters.
#' @param formula A survival or event history formula (`Hist(...)`). The left
#'   hand side is used to compute the expected event status.
#'   If none is given, it is obtained from the prediction object.
#' @param data Data for external validation.
#' @param tLM Landmark times corresponding to the patient entries in data. Only
#'   required if data is a dataframe. tLM can be a string (indicating a column
#'   in data), a vector of length nrow(data), or a single value if all patient
#'   entries were obtained at the same landmark time.
#' @param ID_col Column name that identifies individuals in data. If omitted, it
#'   is obtained from the prediction object.
#' @param se.fit If FALSE or 0, no standard errors are calculated.
#' @param conf.int Confidence interval (CI) coverage. Default is 0.95. If
#'   bootstrapping, CIs are calculated from empirical quantiles. If not, for
#'   right censored data, they are calculated by the package `riskRegression` as
#'   in Blanche et al (references).
#' @param split.method Defines the internal validation design as in
#'   `pec::calPlot`. Options are currently "none" or "bootcv".
#'
#'   "none": assess the model in the test data (`data` argument)/data it was trained on.
#'
#'   "bootcv": `B` models are trained on boostrap samples either drawn with replacement of the same size as the original data or without replacement of size `M`. Models are then assessed in observations not in the sample.
#'
#' @param B Number of times bootstrapping is performed.
#' @param M Subsample size for training in cross-validation. Entries not sampled
#'   in the M subsamples are used for validation.
#' @param cores To perform parallel computing, specifies the number of cores.
#' @param seed Optional, integer passed to set.seed. If not given or NA, no seed
#'   is set.
#' @param unit Time unit for window of prediction, e.g., "year", "month", etc.
#'   Only used for printing results.
#' @param cause Cause of interest if considering competing risks. If left blank,
#'   this is inferred from object.
#' @param ... Additional arguments to pass to Score (`riskRegression` package).
#'   These arguments have been included for user flexibility but have not been
#'  tested and should be used with precaution.
#' @param silent Show any error messages when computing `Score` for each
#'  landmark time (and potentially bootstrap iteration)
#' @return An object of class "LMScore", which has components:
#'   - `auct`: dataframe containing time-dependent AUC if "auc" was
#'     included as a metric
#'   - `briert`: dataframe containing time-dependent Brier score if "brier" was
#'     included as a metric
#' @details See the Github for example code
#' @references Paul Blanche, Cecile Proust-Lima, Lucie Loubere, Claudine Berr, Jean- Francois Dartigues, and Helene Jacqmin-Gadda. Quantifying and comparing dynamic predictive accuracy of joint models for longitudinal marker and time-to-event in presence of censoring and competing risks. Biometrics, 71 (1):102–113, 2015.
#'
#' P. Blanche, J-F Dartigues, and H. Jacqmin-Gadda. Estimating and comparing time-dependent areas under receiver operating characteristic curves for censored event times with competing risks. Statistics in Medicine, 32(30):5381–5397, 2013.
#' @import riskRegression
#' @export
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
           silent = T,
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

    # TODO: turn into an lapply?
    # TODO: parallelize
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

        score_t <- try(
          riskRegression::Score(
            risks_to_test,
            formula = stats::as.formula(formula),
            data = data_to_test,
            metrics = metrics,
            cause = cause,
            times = c(tLM + w - 10e-5),
            se.fit = se.fit.b,
            conf.int = conf.int,
            ...
          ), silent = silent
        )

        if (inherits(score_t, "try-error")) {
          auct_b <- data.frame(model=names(object), times=NA, AUC=NA)
          briert_b <- data.frame(model=names(object), times=NA, Brier=NA)
        } else {
          auct_b <- score_t$AUC$score
          briert_b <- score_t$Brier$score
        }
        auct_b$b <- b
        briert_b$b <- b

        if ("auc" %in% metrics) {
          auct[[b]] <- rbind(auct[[b]], cbind(tLM, auct_b))
        }

        if ("brier" %in% metrics) {
          briert[[b]] <- rbind(briert[[b]], cbind(tLM, briert_b))
        }
      }
    }
    auct <- do.call("rbind", auct)
    briert <- do.call("rbind", briert)
    b_na <- auct$b[is.na(auct$AUC)]
    tLM_na <- auct$tLM[is.na(auct$AUC)]
    model_na <- auct$model[is.na(auct$AUC)]

    if(length(b_na) != 0) {
      message(paste0("\nNote that some metrics could not be computed. Set silent=FALSE to potentially see more error messages. Metrics not computed for (model, b, tLM) = ",paste0("(",model_na,", ",b_na,", ",tLM_na,")", collapse=", ")))
    }

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
