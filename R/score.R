#' Methods (time-dependent AUC and Brier Score) to score the predictive
#' performance of dynamic risk prediction landmark models.
#'
#' There are three ways to perform assess the predictive performance:
#' apparent/internal, bootstrapped, and external. Accordingly, the named list of
#' prediction models must be as follows:
#' * For both apparent/internal evaluation, objects output from
#'   [predict.dynamicLM()] or supermodels fit with [dynls()] may be used as input.
#' * In order to bootstrap, supermodels fit with [dynls()] may be used as input
#'   (note that the argument `x=TRUE` must be specified when fitting the model
#'   in `[dynls()]).
#' * For external calibration, supermodels fit with [dynls()] are input along
#'   with new data in the `data` argument. This data can be a LMdataframe or a
#'   dataframe (in which case `tLM` must be specified).
#'
#' For both internal evaluation and bootstrapping, it is assumed that all
#' models in `object` are fit on the same data.
#'
#'
#' @param object A named list of prediction models, where allowed entries are
#'   outputs from [predict.dynamicLM()] or supermodels from [dynls()] depending on the
#'   type of calibration.
#' @param times Landmark times for which calibration must be plot. These must be
#'   a subset of LM times used during the prediction
#' @param metrics  Character vector specifying which metrics to apply. Choices
#'   are "auc" and "brier". Case matters.
#' @param formula A survival or event history formula (`Hist(...)`). The left
#'   hand side is used to compute the expected event status.
#'   If none is given, it is obtained from the prediction object.
#' @param data Data for external validation.
#' @param tLM Landmark times corresponding to the patient entries in data. Only
#'   required if data is specified and is a dataframe.
#'   tLM can be a string (indicating a column in data), a vector of length
#'   nrow(data), or a single value if all patient entries were obtained at the
#'   same landmark time.
#' @param ID_col Column name that identifies individuals in data. If omitted, it
#'   is obtained from the prediction object.
#' @param se.fit If FALSE or 0, no standard errors are calculated.
#' @param conf.int Confidence interval (CI) coverage. Default is 0.95. If
#'   bootstrapping, CIs are calculated from empirical quantiles. If not, for
#'   right censored data, they are calculated by the package [riskRegression] as
#'   in Blanche et al (references).
#' @param split.method Defines the internal validation design. Options are
#'   currently "none" or "bootcv".
#'
#'   "none": assess the model in the test data (`data` argument)/data it was
#'   trained on.
#'
#'   "bootcv": `B` models are trained on boostrap samples either drawn with
#'     replacement of the same size as the original data or without replacement
#'     of size `M`. Models are then assessed in observations not in the sample.
#'
#' @param B Number of times bootstrapping is performed.
#' @param M Subsample size for training in cross-validation. Entries not sampled
#'   in the M subsamples are used for validation.
#' @param cores To perform parallel computing, specifies the number of cores.
#'   (Not yet implemented)
#' @param seed Optional, integer passed to set.seed. If not given or NA, no seed
#'   is set.
#' @param unit Time unit for window of prediction, e.g., "year", "month", etc.
#'   Only used for printing results.
#' @param cause Cause of interest if considering competing risks. If left blank,
#'   this is inferred from object.
#' @param ... Additional arguments to pass to [riskRegression::Score()].
#'   These arguments have been included for user flexibility but have not been
#'   tested and should be used with precaution.
#' @param silent Show any error messages when computing `Score` for each
#'   landmark time (and potentially bootstrap iteration)
#' @param na.rm Ignore bootstraps where there are errors (for example not
#'   enough datasamples) and calculate metrics on remaining values. This is not
#'   recommended. For example, if only one bootstrap sampling has enough data
#'   that live to the prediction window, the standard error will be zero.
#' @return An object of class "LMScore", which has components:
#'   - `auct`: dataframe containing time-dependent AUC if "auc" was
#'     included as a metric
#'   - `briert`: dataframe containing time-dependent Brier score if "brier" was
#'     included as a metric
#' @details See the Github for example code.
#'
#'   If data at late evaluation times is sparse, certain bootstrap samples may
#'   not have patients that live long enough to perform evaluation leading to
#'   the message "Upper limit of followup in bootstrap samples, was too low.
#'   Results at evaluation time(s) beyond these points could not be computed
#'   and are left as NA". In this case, consider only evaluating for earlier
#'   landmarks or performing prediction with a smaller window as data points are
#'   slim. If you wish to see which model/bootstrap/landmark times failed, set
#'   SILENT=FALSE. Set na.rm = TRUE ignores these bootstraps and calculate
#'   metrics from the bootstrap samples that worked (not recommended).
#'
#'   Another message may occur: "Dropping bootstrap b = {X} for model {name} due
#'   to unreliable predictions". As certain approximations are made, numerical
#'   overflow sometimes occurs in predictions for bootstrapped samples. To avoid
#'   potential errors, the whole bootstrap sample is dropped in this case. Note
#'   that input data should be complete otherwise this may occur
#'   unintentionally.
#'
#' @references Paul Blanche, Cecile Proust-Lima, Lucie Loubere, Claudine Berr,
#'   Jean- Francois Dartigues, and Helene Jacqmin-Gadda. Quantifying and
#'   comparing dynamic predictive accuracy of joint models for longitudinal
#'   marker and time-to-event in presence of censoring and competing risks.
#'   Biometrics, 71 (1):102–113, 2015.
#'
#'   P. Blanche, J-F Dartigues, and H. Jacqmin-Gadda. Estimating and comparing
#'   time-dependent areas under receiver operating characteristic curves for
#'   censored event times with competing risks. Statistics in Medicine,
#'   32(30):5381–5397, 2013.
#'
#' @import riskRegression
#' @importFrom data.table .SD
#' @export
Score <-
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
           na.rm = FALSE,
           ...) {

    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package \"data.table\" must be installed to use function Score",
           call. = FALSE)
    }

    get.auc <- FALSE; get.brier <- FALSE
    if ("auc" %in% metrics) get.auc <- TRUE
    if ("brier" %in% metrics) get.brier <- TRUE

    checked_input <- match.call()
    m <- match(c("object", "times", "formula", "data", "tLM", "ID_col",
                 "split.method", "B", "M", "cores", "seed", "cause"), names(checked_input), 0L)
    checked_input <- as.list(checked_input[m])
    checked_input <- do.call(check_evaluation_inputs, checked_input)

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

    num_B = length(unique(data$b))
    se.fit.b <- se.fit
    if (num_B > 1) se.fit.b <- FALSE

    # TODO: parallelize?
    metrics <- lapply(unique(data$b), function(b){
      m_b <- lapply(1:length(times), function(t){
        tLM = times[t]

        idx = (pred_LMs == tLM) & (data$b == b)
        data_to_test = data[idx,]
        risks_to_test = lapply(1:NF, function(i) {
          preds[idx, i]
        })
        names(risks_to_test) = names(object)

        if (nrow(data_to_test) != length(risks_to_test[[1]])) {
          stop("nrow(data_to_test)!=length(risks_to_test)")
        }

        score_t <- suppressMessages(try(
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
          ), silent = T
        ))

        if (inherits(score_t, "try-error")) {
          auct_b <- data.frame(model=names(object), times=NA, AUC=NA)
          briert_b <- data.frame(model=names(object), times=NA, Brier=NA)
        } else {
          auct_b <- score_t$AUC$score
          briert_b <- score_t$Brier$score
        }
        auct_b$b <- b; briert_b$b <- b

        metrics_b_t <- list()
        if (get.auc) metrics_b_t$AUC <- cbind(tLM, auct_b)
        if (get.brier) metrics_b_t$Brier <- cbind(tLM, briert_b)
        metrics_b_t
      })

      m_out <- list()
      if (get.auc) m_out$AUC <- do.call("rbind", lapply(m_b, function(m) m$AUC))
      if (get.brier) m_out$Brier <- do.call("rbind", lapply(m_b, function(m) m$Brier))
      m_out
    })

    auct <- do.call("rbind", lapply(metrics, function(m) m$AUC))
    briert <- do.call("rbind", lapply(metrics, function(m) m$Brier))

    if (B > 1) {
      if (se.fit==TRUE) {
        alpha = 1-conf.int
        auct_out <- auct[,
                     data.table::data.table(
                       mean(.SD[["AUC"]],na.rm=na.rm),
                       se=stats::sd(.SD[["AUC"]],na.rm=T),
                       lower=stats::quantile(.SD[["AUC"]],alpha/2,na.rm=T),
                       upper=stats::quantile(.SD[["AUC"]],(1-alpha/2),na.rm=T)),
                     by=c("model","tLM"),.SDcols="AUC"
                     ]
        briert_out <- briert[,
                         data.table::data.table(
                           mean(.SD[["Brier"]],na.rm=na.rm),
                           se=stats::sd(.SD[["Brier"]],na.rm=T),
                           lower=stats::quantile(.SD[["Brier"]],alpha/2,na.rm=T),
                           upper=stats::quantile(.SD[["Brier"]],(1-alpha/2),na.rm=T)),
                         by=c("model","tLM"),.SDcols="Brier"
                         ]
        data.table::setnames(auct_out,c("model","tLM","AUC","se","lower","upper"))
        data.table::setnames(briert_out,c("model","tLM","Brier","se","lower","upper"))

        if (!silent) {
          b_na <- auct$b[is.na(auct$AUC)]
          tLM_na <- auct$tLM[is.na(auct$AUC)]
          model_na <- auct$model[is.na(auct$AUC)]

          if(length(b_na) != 0) {
            message("Upper limit of followup in bootstrap samples was too low. Results at evaluation time(s) beyond these points could not be computed and are left as NA.")
            message(paste0("Metrics not computed for (model, b, tLM) = ",
                       paste0("(",model_na,", ",b_na,", ",tLM_na,")",
                              collapse=", ")))
          }
        }

        if (na.rm == FALSE) {
          auct_out <- auct_out[c(is.na(auct_out[,3])), `:=` ("se"=NA,"lower"=NA,"upper"=NA)]
          briert_out <- briert_out[c(is.na(briert_out[,3])), `:=` ("se"=NA,"lower"=NA,"upper"=NA)]
        }


      } else {
        auct_out <- auct[,data.table::data.table(mean(.SD[["AUC"]],na.rm=na.rm)),by=c("model","tLM"),.SDcols="AUC"]
        briert_out <- briert[,data.table::data.table(mean(.SD[["Brier"]],na.rm=na.rm)),by=c("model","tLM"),.SDcols="Brier"]
        data.table::setnames(auct_out,c("model","tLM","AUC"))
        data.table::setnames(briert_out,c("model","tLM","Brier"))
      }
    } else {
      auct_out <- auct
      briert_out <- briert
    }

    outlist <- list(
      auct = auct_out,
      briert = briert_out,
      w = w,
      unit = unit
    )
    class(outlist) = "LMScore"
    return(outlist)
  }
