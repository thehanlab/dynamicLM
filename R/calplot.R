#' Calibration plots for dynamic risk prediction landmark models.
#'
#' There are three ways to perform calibration: apparent/internal, bootstrapped,
#' and external. Accordingly, the named list of prediction models must be as
#' follows:
#' * For both apparent/internal calbration, objects output from
#'   [predict.dynamicLM()] for supermodels fit with [dynamic_lm()] may be used
#'   as input.
#' * In order to bootstrap, supermodels fit with [dynamic_lm()] may be used as
#'   input (note that the argument `x=TRUE` must be specified when fitting the
#'   model in [dynamic_lm()]).
#' * For external calibration, supermodels fit with [dynamic_lm()] are input
#'   along with new data in the `data` argument. This data can be a LMdataframe
#'   or a dataframe (in which case `lms` must be specified).
#'
#' For both internal calibration and bootstrapping, it is assumed that all
#' models in `object` are fit on the same data.
#'
#' @param object A named list of prediction models, where allowed entries are
#'   outputs from [predict.dynamicLM()] or supermodels from [dynamic_lm()]
#'   depending on the type of calibration.
#' @param times Landmark times for which calibration must be plot. These must be
#'   a subset of landmark times used during the prediction
#' @param formula A survival or event history formula (`Hist(...)`). The left
#"   hand side is used to compute the expected event status.
#'   If none is given, it is obtained from the prediction object.
#' @param data Data for external validation. This can be an object of class
#'   LMdataframe (i.e., created by calling [stack_data()] and
#'   [add_interactions()]), or a data.frame. If it is a data.frame, argument
#'   `lms` must be specified.
#' @param lms Landmark times corresponding to the patient entries in data. Only
#'   required if `data` is specified and is a dataframe.
#'   `lms` can be a string (indicating a column in data), a vector of length
#'   nrow(data), or a single value if all patient entries were obtained at the
#'   same landmark time.
#' @param id_col Column name that identifies individuals in data. If omitted, it
#'   is obtained from the prediction object.
#' @param split.method Defines the internal validation design as in
#'   [pec::calPlot()]. Options are currently "none" or "bootcv".
#'
#'   "none": assess the model in the test data (`data` argument)/data it was
#"   trained on.
#'
#'   "bootcv": `B` models are trained on bootstrap samples either drawn with
#"   replacement of the same size as the original data or without replacement of
#'   size `M`. Models are then assessed in observations not in the sample.
#'
#' @param B Number of times bootstrapping is performed.
#' @param M Subsample size for training in cross-validation. Entries not sampled
#"   in the M subsamples are used for validation.
#' @param cores To perform parallel computing, specifies the number of cores.
#'   (Not yet implemented)
#' @param seed Optional, integer passed to set.seed. If not given or NA, no seed
#"   is set.
#' @param regression_values Default is FALSE. If set to TRUE, the returned list
#'   is appended by another list `regression_values`,
#'   which contains the intercept and slope of a linear regression of each model
#'   for each landmark time (i.e., each calibration plot).
#'   Note that perfect calibration has a slope of 1 and an intercept of 0.
#' @param cause Cause of interest if considering competing risks. If left blank,
#'   this is inferred from object.
#' @param plot If FALSE, do not plot the results, just return a plottable
#'   object. Default is TRUE.
#' @param main Optional title to override default.
#' @param ... Additional arguments to pass to calPlot (`pec` package).
#'   These arguments have been included for user flexibility but have not been
#'   tested and should be used with precaution.
#'
#' @return List of plots of w-year risk, one entry per prediction/landmark time
#'   point. List has a component `$regression_values` (if argument
#'   regression_values is set to TRUE) which is a list of which contains the
#'   intercept and slope of a linear regression of each model
#'   for each landmark time (i.e., each calibration plot).
#'
#' @details When collecting bootstrap samples, the same individuals are
#'   considered across landmarks.
#'   I.e., sample `M` unique individuals, train on the super dataset formed by
#'   these individuals, and validate on the individuals not sampled at the
#'   landmarks they remain alive (or that are given in `times`).
#'
#'  Note that only complete cases of data are considered (whatever type of
#'  calibration is performed).
#'
#'  A comment on the following message:
#'  "Dropping bootstrap b = ... for model ... due
#'  to unreliable predictions". As certain approximations are made, numerical
#'  overflow sometimes occurs in predictions for bootstrapped samples. To avoid
#'  potential errors, the whole bootstrap sample is dropped in this case. Note
#'  that input data should be complete otherwise this may occur
#'  unintentionally. Calibration plots are still produced excluding predictions
#'  made during the bootstrap resampling.
#'
#' @examples
#' \dontrun{
#' # Internal validation
#' par(mfrow = c(2, 2), pty = "s")
#' outlist <- calplot(list("Model1" = supermodel),
#'                    method = "quantile", q = 5,  # method for calibration plot
#'                    regression_values = TRUE,    # output regression values
#'                    ylim = c(0, 0.4), xlim = c(0, 0.4)) # optional
#' outlist$regression_values
#'
#' # Bootstrapping
#' # Remember to fit the supermodel with argument 'x = TRUE'
#' par(mfrow = c(2, 2), pty = "s")
#' outlist <- calplot(list("Model1" = supermodel),
#'                    method = "quantile", q = 5,
#'                    split.method = "bootcv", B = 10, # 10 bootstraps
#'                    ylim = c(0, 0.4), xlim = c(0, 0.4))
#'
#' # External validation
#' # a) newdata is a dataframe
#' newdata <- relapse[relapse$T_txgiven == 0, ]
#' newdata$age <- newdata$age.at.time.0
#' newdata$LM <- 0
#' par(mfrow = c(1, 1))
#' cal <- calplot(list("Model1" = supermodel), data = newdata, lms = "LM",
#'                method = "quantile", q = 5, ylim = c(0, 0.1), xlim = c(0, 0.1))
#'
#' # b) newdata is a landmark dataset
#' par(mfrow = c(2, 2), pty = "s")
#' lmdata_new <- lmdata
#' cal <- calplot(list("Model1" = supermodel), data = lmdata_new,
#'                method = "quantile", q = 10, ylim = c(0, 0.4), xlim = c(0, 0.4))
#' }
#' @import prodlim
#' @export
#' @seealso [dynamicLM::score()], [pec::calPlot()]
#'
calplot <-
  function(object,
           times,
           formula,
           data,
           lms,
           id_col = "ID",
           split.method = "none",
           B = 1,
           M,
           cores = 1,
           seed,
           regression_values = FALSE,
           cause,
           plot = TRUE,
           main,
           ...) {

    ### Check input and set up some initial variables ###

    if (!requireNamespace("pec", quietly = TRUE)) {
      stop("Package \"pec\" must be installed to use function calplot",
           call. = FALSE)
    }
    if (!(inherits(object, "list"))) stop("`object` must be a named list.")


    checked_input <- match.call()
    m <- match(c("object", "times", "formula", "data", "lms", "id_col",
                  "split.method", "B", "M", "cores", "seed", "cause"),
               names(checked_input), 0L)
    checked_input <- as.list(checked_input[m])
    checked_input <- do.call(check_evaluation_inputs, checked_input)

    object <- checked_input$object
    times <- checked_input$times
    preds <- checked_input$preds
    pred_LMs <- checked_input$pred_LMs
    data <- checked_input$data
    NF <- checked_input$NF
    w <- checked_input$w
    formula <- checked_input$formula
    cause <- checked_input$cause
    indicator <- checked_input$indicator

    # Plotting parameters
    add_title <- TRUE
    if (!missing(main)) {
      add_title <- FALSE
    }

    # Calibration
    outlist <- list()
    if (regression_values) reg_values_list <- list()
    for (t in seq_along(times)) {
      tLM <- times[t]

      idx <- pred_LMs == tLM
      data_to_test <- data[idx, ]
      risks_to_test <- lapply(1:NF, function(i) {
        preds[idx, i]
      })
      names(risks_to_test) <- names(object)
      if (object[[1]]$type == "coxph") {
        risks_to_test <- lapply(risks_to_test, function(r) 1 - r)
      }

      if (nrow(data_to_test) != length(risks_to_test[[1]])) {
        stop("nrow(data_to_test)!=length(risks_to_test)")
      }

      if (nrow(data_to_test) == 0) {
        warning(tidymess(paste0(
          "Skipping calplot for landmark time ", tLM, " as no data was provided
          for this landmark.")))
        next
      }

      x <- NULL
      x <- pec::calPlot(
        risks_to_test,
        time = tLM + w,
        formula = stats::as.formula(formula),
        data = data_to_test,
        cause = cause,
        plot = plot,
        ...
      )

      if (regression_values){
        for (i in 1:NF) {
          name <- names(object)[i]
          coefs <- stats::lm(Obs ~ Pred, x$plotFrames[[name]])$coefficients
          coefs <- c(tLM, coefs)
          names(coefs)[1] <- "LM"
          reg_values_list[[name]] <- rbind(reg_values_list[[name]], coefs)
        }
      }

      if (plot) {
        if (add_title) {
          title <- paste("Risk calibration at landmark", tLM)
          graphics::title(main = title)
        } else {
          graphics::title(main = main)
        }

        # Older: Sub gave the number at risk and with an event
        # if (sub) {
        #   num_patients <- length(risks_to_test[[1]])
        #   status <- object[[1]]$outcome$status
        #   num_events <- sum(data_to_test[[status]] == indicator)
        #   subtitle <- paste0("#at risk=", num_patients,
        #                     ",#that undergo event=", num_events)
        #   graphics::title(sub = substitute(paste(italic(subtitle))))
        # }

        # Potential addition: sub shows the regression values
        # @param sub If TRUE and `regression_values` is also set to TRUE, add a
        #    subheading with the regression slope and intercept.
        # @param digits If `sub` and `regression_values` are TRUE, determines
        #    the number of digits to print.
        # if (sub && regression_values) {
        #   slopes <- sapply(1:NF, function(i) {
        #     reg_values_list[[i]][nrow(reg_values_list[[i]]), "Pred"]
        #   })
        #   intercepts <- sapply(1:NF, function(i) {
        #     reg_values_list[[i]][nrow(reg_values_list[[i]]), "(Intercept)"]
        #   })
        #   subtitle <- paste0(names(object),
        #                      ": slope:", round(slopes, digits),
        #                      ", intercept:", round(intercepts, digits), "\n")
        #   graphics::title(sub = substitute(paste(italic(subtitle))))
        # }
      }

      outlist[[t]] <- x

    }
    if (regression_values) {
      for (i in 1:NF) {
        name <- names(object)[i]
        rownames(reg_values_list[[name]]) <- NULL
        colnames(reg_values_list[[name]]) <- c("LM", "Intercept", "Slope")
      }
      outlist[["regression_values"]] <- reg_values_list
    }
    class(outlist) <- "LMcalibrationPlot"
    return(outlist)
  }
