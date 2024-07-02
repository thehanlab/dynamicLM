#' Methods (time-dependent AUC and Brier Score) to score the predictive
#' performance of dynamic risk prediction landmark models.
#'
#' There are three ways to perform assess the predictive performance:
#' apparent/internal, bootstrapped, and external. Accordingly, the named list of
#' prediction models must be as follows:
#' * For both apparent/internal evaluation, objects output from
#'   [predict.dynamicLM()] or supermodels fit with [dynamic_lm()] may be used as
#'   input.
#' * In order to bootstrap, supermodels fit with [dynamic_lm()] may be used as
#'   input (note that the argument `x=TRUE` must be specified when fitting the
#'   model in [dynamic_lm()]).
#' * For external calibration, supermodels fit with [dynamic_lm()] are input
#'   along with new data in the `data` argument. This data can be a LMdataframe
#'   or a dataframe (in which case `lms` must be specified).
#'
#' For both internal evaluation and bootstrapping, it is assumed that all
#' models in `object` are fit on the same data.
#'
#'
#' @param object A named list of prediction models, where allowed entries are
#'   outputs from [predict.dynamicLM()] or supermodels from [dynamic_lm()]
#'   depending on the type of calibration.
#' @param times Landmark times for which calibration must be plot. These must be
#'   a subset of landmark times used during the prediction
#' @param metrics  Character vector specifying which metrics to apply. Choices
#'   are "auc" and "brier".
#' @param formula A survival or event history formula
#'   ([prodlim::Hist()]). The left hand side is used to compute the
#'   expected event status. If none is given, it is obtained from the prediction
#'   object.
#' @param data Data for external validation.
#' @param lms Landmark times corresponding to the patient entries in data. Only
#'   required if data is specified and is a dataframe.
#'   `lms` can be a string (indicating a column in data), a vector of length
#'   nrow(data), or a single value if all patient entries were obtained at the
#'   same landmark time.
#' @param id_col Column name that identifies individuals in data. If omitted, it
#'   is obtained from the prediction object.
#' @param se.fit If FALSE or 0, no standard errors are calculated.
#' @param conf.int Confidence interval (CI) coverage. Default is 0.95. If
#'   bootstrapping, CIs are calculated from empirical quantiles. If not, for
#'   right censored data, they are calculated by the package [riskRegression] as
#'   in Blanche et al (references).
#' @param contrasts If TRUE, perform model comparison tests.
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
#' @param summary TODO
#' @param cores To perform parallel computing, specifies the number of cores.
#'   (Not yet implemented)
#' @param seed Optional, integer passed to set.seed. If not given or NA, no seed
#'   is set.
#' @param cause Cause of interest if considering competing risks. If left blank,
#'   this is inferred from object.
#' @param ... Additional arguments to pass to [riskRegression::Score()].
#'   These arguments have been included for user flexibility but have not been
#'   tested and should be used with precaution.
#' @param silent Show any error messages when computing `score` for each
#'   landmark time (and potentially bootstrap iteration)
#'
#' @return An list with entries `AUC` and `Brier` if "auc" and "brier" were
#'   included as metrics respectively and `AUC_Summary` and/or `Brier_summary`
#'   if `summary` is not null. Each will have entries:
#'   - `score`: data.table containing the metric
#'   - `contrasts`: data.table containing model comparisons
#' @details If data at late landmark times is sparse, some bootstrap samples may
#'   not have patients that live long enough to perform evaluation leading to
#'   the message "Upper limit of followup in bootstrap samples, was too low.
#'   Results at evaluation time(s) beyond these points could not be computed
#'   and are left as NA". In this case, consider only evaluating for earlier
#'   landmarks or performing prediction with a smaller window as data points are
#'   slim. If you wish to see which model/bootstrap/landmark times failed, set
#'   SILENT=FALSE. Currently ignores these bootstraps and calculates
#'   metrics from the bootstrap samples that worked.
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
#' @examples
#' \dontrun{
#' # Internal validation
#' scores <- score(list("Model1" = supermodel),
#'                 times = c(0, 6)) # lms at which to provide calibration plots
#' scores
#'
#' # Bootstrapping
#' # Remember to fit the supermodel with argument 'x = TRUE'
#' scores <- score(list("Model1" = supermodel),
#'                 times = c(0, 6),
#'                 split.method = "bootcv", B = 10) # 10 bootstraps
#' scores
#'
#' par(mfrow=c(1,2))
#' plot(scores)
#'
#' # External validation
#' # Either input an object from predict as the object or a supermodel and
#' # "data" & "lms" argument
#' newdata <- relapse[relapse$T_txgiven == 0, ]
#' newdata$age <- newdata$age.at.time.0
#' newdata$LM <- 0
#' score(list("CSC" = supermodel), cause = 1, data = newdata, lms = "LM")
#' }
#'
#' @import riskRegression
#' @importFrom data.table .SD
#' @export
# TODO: add null.model argument
# TODO: handle na.rm
# TODO: add checks for left-censoring for summary = TRUE
# TODO: add bootstrapping option for summary metric
score <-
  function(object,
           times,
           metrics = c("auc", "brier"),
           formula,
           data,
           lms = "LM",
           id_col,
           se.fit = TRUE,
           conf.int = 0.95,
           contrasts = TRUE,
           split.method = "none",
           B = 1,
           M,
           summary = TRUE,
           cores = 1,
           seed,
           cause,
           silent = TRUE,
           ...) {

    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package 'data.table' must be installed to use function score()",
           call. = FALSE)
    }

    # Get cleaned input
    checked_input <- match.call()
    m <- match(c("object", "times", "formula", "data", "lms", "id_col",
                 "split.method", "B", "M", "cores", "seed", "cause"),
               names(checked_input), 0L)
    checked_input <- as.list(checked_input[m])
    checked_input <- do.call(check_evaluation_inputs, checked_input,
                             envir = parent.frame())
    object <- checked_input$object
    times <- checked_input$times
    preds <- checked_input$preds
    pred_LMs <- checked_input$pred_LMs
    data <- checked_input$data
    NF <- checked_input$NF
    w <- checked_input$w
    formula <- checked_input$formula
    cause <- checked_input$cause
    id_col <- checked_input$id_col

    # Determine which metrics to calculate
    get.auc <- "auc" %in% tolower(metrics)
    get.bs <- "brier" %in% tolower(metrics)
    if (!get.auc && !get.bs)
      stop(tidymess("At least one of the following metrics must be specified:
                    \"auc\", \"brier\"."))
    get.a.iid <- get.auc && summary
    get.b.iid <- get.bs && summary

    if (!all("bootstrap" %in% colnames(data))) data$bootstrap <- 1
    if (length(object) == 1) contrasts <- FALSE

    num_B <- length(unique(data$bootstrap))
    se.fit.b <- if (num_B > 1) FALSE else se.fit

    # TODO: parallelize?
    # Get the metrics for each bootstrap (if not bootstrapping then for the
    # single iteration)
    metrics <- lapply(unique(data$bootstrap), function(b) {
      lapply(seq_along(times), function(t) {
        tLM <- times[t]
        idx <- (pred_LMs == tLM) & (data$bootstrap == b)
        data_to_test <- data[idx, ]

        #       by only including those who live after s, our censoring times
        #       etc align properly but we can also do it as below and set
        #       times = w
        # data_to_test[[time]] <- data_to_test[[time]] - tLM
        # data_to_test[["LM"]] <- data_to_test[["LM"]] - tLM

        risks_to_test <- lapply(1:NF, function(i) preds[idx, i])
        names(risks_to_test) <- names(object)

        if (nrow(data_to_test) != length(risks_to_test[[1]]))
          stop("nrow(data_to_test)!=length(risks_to_test)")

        if (nrow(data_to_test) == 0) {
          warning(tidymess(paste0(
            "Skipping calplot for landmark time ", tLM, " as no data was
            provided for this landmark.")))
          return(NULL)
        }

        score_t <- suppressMessages(try(
          riskRegression::Score(
            risks_to_test,
            formula = stats::as.formula(formula),
            data = data_to_test,
            metrics = metrics,
            cause = cause,
            times = tLM + w - 10e-5,
            # se.fit = se.fit.b,
            conf.int = conf.int,
            keep = c("residuals", "iid"),
            ...
          ), silent = TRUE
        ))
        # print(score_t)

        # Error handing & get what we need from the output
        if (inherits(score_t, "try-error")) {
          # If error: df with one row of NAs of the relevant column names
          auct_b <- initialize_df("AUC")
          briert_b <- initialize_df("Brier")
          auc_contrasts_b <- initialize_df("delta.AUC")
          brier_contrasts_b <- initialize_df("delta.Brier")
          a_iid <- initialize_df("IF.AUC")
          b_iid <- initialize_df("IF.Brier")
        } else {
          auct_b <- cbind(tLM, score_t$AUC$score, bootstrap = b)
          briert_b <- cbind(tLM, score_t$Brier$score, bootstrap = b)
          auc_contrasts_b <- cbind(tLM, score_t$AUC$contrasts, bootstrap = b)
          brier_contrasts_b <- cbind(tLM, score_t$Brier$contrasts,
                                     bootstrap = b)
          if (get.a.iid)
            a_iid <- cbind(tLM, score_t$AUC$iid.decomp, bootstrap = b)
          if (get.b.iid)
            b_iid <- cbind(tLM, score_t$Brier$iid.decomp, bootstrap = b)
        }

        list(
          AUC = if (get.auc) auct_b else NULL,
          Brier = if (get.bs) briert_b else NULL,
          a_contrasts = if (contrasts && get.auc) auc_contrasts_b else NULL,
          b_contrasts = if (contrasts && get.bs) brier_contrasts_b else NULL,
          a_iid = if (get.a.iid) a_iid else NULL,
          b_iid = if (get.b.iid) b_iid else NULL
        )
      })
    })

    # Convert to dataframe
    get_metrics <- function(metrics, metric_name) {
      do.call("rbind", lapply(metrics, function(mb)
        do.call("rbind", lapply(mb, function(mi) mi[[metric_name]]))))
    }
    auct <- get_metrics(metrics, "AUC")
    briert <- get_metrics(metrics, "Brier")
    a_contrasts <- get_metrics(metrics, "a_contrasts")
    b_contrasts <- get_metrics(metrics, "b_contrasts")
    a_iid <- get_metrics(metrics, "a_iid")
    b_iid <- get_metrics(metrics, "b_iid")

    # Clean output and handle bootstrapping results
    if (B > 1) {
      alpha <- 1 - conf.int
      clean_if <- function(condition, df, metric, ...) {
        if (condition) clean_bootstraps(df, metric, alpha, se.fit = se.fit, ...)
        else NULL
      }

      auct_out <- clean_if(get.auc, auct, "AUC")
      a_contrasts_out <- clean_if(contrasts && get.auc, a_contrasts,
                                  "delta.AUC", contrasts = TRUE)
      briert_out <- clean_if(get.bs, briert, "Brier")
      b_contrasts_out <- clean_if(contrasts && get.bs, b_contrasts,
                                  "delta.Brier", contrasts = TRUE)

      if (!silent) {
        b_na <- auct$bootstrap[is.na(auct$AUC)]
        tlm_na <- auct$tLM[is.na(auct$AUC)]
        model_na <- auct$model[is.na(auct$AUC)]
        if (length(b_na) != 0) {
          message(tidymess(paste0(
          "Upper limit of followup in bootstrap samples was too low. Results
          at evaluation time(s) beyond these points could not be computed and
          are left as NA.\nMetrics not computed for (model, b, tLM) = ",
          paste0("(", model_na, ", ", b_na, ", ", tlm_na, ")",
                 collapse = ", "))))
        }
      }
      # if (na.rm == FALSE) {
      #   auct_out <- auct_out[is.na(auct_out[, 3]),
      #                        `:=`("se" = NA, "lower" = NA, "upper" = NA)]
      #   briert_out <- briert_out[is.na(briert_out[, 3]),
      #                            `:=`("se" = NA, "lower" = NA, "upper" = NA)]
      # }

    } else { # B == 1
      auct_out <- auct
      briert_out <- briert
      a_contrasts_out <- a_contrasts
      b_contrasts_out <- b_contrasts
      auct_out$bootstrap <- briert_out$bootstrap <- NULL
      if (contrasts) {
        a_contrasts_out$bootstrap <- b_contrasts_out$bootstrap <- NULL
      }
    }

    outlist <- list()
    if (get.auc)
      outlist$AUC <- list(score = auct_out, contrasts = a_contrasts_out)
    if (get.bs)
      outlist$Brier <- list(score = briert_out, contrasts = b_contrasts_out)
    if (B > 1) {
      outlist$B <- B
      outlist$split.method <- split.method
    }

    # Get summary metrics
    if (summary) {
      if (get.a.iid) {
        outlist$AUC_summary <- summary_metric(
          "AUC", auct, a_contrasts, a_iid, conf.int, object, id_col, B, se.fit)
      }
      if (get.b.iid) {
        outlist$Brier_summary <- summary_metric(
          "Brier", briert, b_contrasts, b_iid, conf.int, object,
          id_col, B, se.fit)
      }
    }

    outlist$w <- w
    class(outlist) <- "LMScore"
    return(outlist)
  }
