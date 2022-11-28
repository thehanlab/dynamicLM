#' Calibration plots for dynamic risk prediction landmark models.
#'
#' There are three ways to perform calibration: apparent/internal, bootstrapped,
#' and external. Accordingly, the named list of prediction models must be as
#' follows:
#' * For both apparent/internal calbration, objects output from `predLMrisk`
#'   for supermodels fit with `fitLM` may be used as input.
#' * In order to bootstrap, supermodels fit with `fitLM` may be used as input
#'   (note that the argument `x=TRUE` must be specified when fitting the model
#'   in `fitLM`).
#' * For external calibration, supermodels fit with `fitLM` are input along with
#'   new data in the `data` argument. This data can be a LMdataframe or a
#'   dataframe (in which case `tLM` must be specified).
#'
#' For both internal calibration and bootstrapping, it is assumed that all
#' models in `object` are fit on the same data.
#'
#' @param object A named list of prediction models, where allowed entries are outputs from `predLMrisk` or supermodels from `fitLM` depending on the type of calibration.
#' @param times Landmark times for which calibration must be plot. These must be a subset of LM times used during the prediction
#' @param formula A survival or event history formula (`Hist(...)`). The left hand side is used to compute the expected event status.
#'   If none is given, it is obtained from the prediction object.
#' @param data Data for external validation.
#' @param tLM Landmark times corresponding to the patient entries in data. Only required if data is a dataframe.
#'   tLM can be a string (indicating a column in data), a vector of length nrow(data),
#'   or a single value if all patient entries were obtained at the same landmark time.
#' @param ID_col Column name that identifies individuals in data. If omitted, it is obtained from the prediction object.
#' @param split.method Defines the internal validation design as in `pec::calPlot`. Options are currently "none" or "bootcv".
#'
#'   "none": assess the model in the test data (`data` argument)/data it was trained on.
#'
#'   "bootcv": `B` models are trained on boostrap samples either drawn with replacement of the same size as the original data or without replacement of size `M`. Models are then assessed in observations not in the sample.
#'
#' @param B Number of times bootstrapping is performed.
#' @param M Subsample size for training in cross-validation. Entries not sampled in the M subsamples are used for validation.
#' @param cores To perform parallel computing, specifies the number of cores.
#' @param seed Optional, integer passed to set.seed. If not given or NA, no seed is set.
#' @param regression_values Default is FALSE. If set to TRUE, the returned list is appended by a list `regression_values`,
#'   which contains the intercept and slope of a linear regression of each model for each landmark time (i.e., each calibration plot).
#'   Note that perfect calibration has a slope of 1 and an intercept of 0.
#' @param unit Time unit for window of prediction, e.g., "year", "month", etc. Only used for printing results.
#' @param cause Cause of interest if considering competing risks. If left blank, this is inferred from object.
#' @param plot If FALSE, do not plot the results, just return a plottable object. Default is TRUE.
#' @param main Optional title to override default.
#' @param sub If TRUE, add a subheading with the number of individuals at risk, and the number that under the event of interest.
#'   Default is TRUE. Set to FALSE for bootstrapping.
#' @param ... Additional arguments to pass to calPlot (`pec` package).
#'   These arguments have been included for user flexibility but have not been
#'   tested and should be used with precaution.
#'
#' @return List of plots of w-year risk, one entry per prediction/landmark time point
#' @details When collecting bootstrap samples, the same individuals are considered across landmarks.
#'   I.e., sample `M` unique individuals, train on the super dataset formed by these individuals, and validate on the individuals not sampled at the landmarks they remain alive (or that are given in `times`).
#'
#'  Note that only complete cases of data are considered (whatever type of calibration is performed). Furthermore, most errors in plotting occur when a formula is not given. Formulas can look like `Hist(Time,event,LM)~1` / similar...
#'
#'  See the [github](https://github.com/thehanlab/dynamicLM) for detailed example code.
#' @examples
#' \dontrun{
#' par(mfrow=c(2,2),pty="s")
#' outlist = LMcalPlot(list("Model1"=p1),
#'                     unit="month",            # for the title
#'                     times=c(6,12,18,24),     # landmarks at which to provide calibration plots
#'                     formula="Hist(event,Time,LM)~1",
#'                     method="quantile", q=10, # method for calibration plot
#'                     ylim=c(0,0.4), xlim=c(0,0.4))
#' }
#' @import prodlim
#' @export
#'
LMcalPlot <-
  function(object,
           times,
           formula,
           data,
           tLM,
           ID_col="ID",
           split.method = "none",
           B = 1,
           M,
           cores = 1,
           seed,
           regression_values = FALSE,
           unit = "year",
           cause,
           plot = T,
           main,
           sub = T,
           ...) {

    ### Check input and set up some initial variables ###

    if (!requireNamespace("pec", quietly = T)) {
      stop("Package \"pec\" must be installed to use function LMcalPlot.", call. = F)
    }
    if (!(class(object)=="list")) stop("object must be a named list.")


    checked_input <- match.call()
    m <- match(c("object", "times", "formula", "data", "tLM", "ID_col",
                  "split.method", "B", "M", "cores", "seed", "cause"), names(checked_input), 0L)
    checked_input <- checked_input[c(1L, m)]
    checked_input[[1L]] <- quote(check_evaluation_inputs)
    # print(checked_input)
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
    indicator = checked_input$indicator

    # Plotting parameters
    add_title = T
    if (!missing(main)) {
      add_title = F
    }

    # Calibration
    outlist = list()
    if (regression_values) reg_values_list = list()
    for (t in 1:length(times)) {
      tLM = times[t]

      idx = pred_LMs == tLM
      data_to_test = data[idx,]
      risks_to_test = lapply(1:NF, function(i) {
        preds[idx, i]
      })
      names(risks_to_test) = names(object)
      if (object[[1]]$type == "coxph") {
        risks_to_test = lapply(risks_to_test, function(r) 1-r)
      }

      if (nrow(data_to_test) != length(risks_to_test[[1]])) {
        stop("nrow(data_to_test)!=length(risks_to_test)")
      }

      x = NULL
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
          name = names(object)[i]
          coefs = lm(Pred ~ Obs, x$plotFrames[[name]])$coefficients
          coefs = c(tLM, coefs)
          names(coefs)[1] = "LM"
          reg_values_list[[name]] = rbind(reg_values_list[[name]], coefs)
        }
      }

      if (plot) {
        if (add_title) {
          title = paste0("Calibration of ",w,"-",unit," risk \n measured at LM time ",tLM)
          graphics::title(main = title)
        }
        else {
          graphics::title(main = main)
        }

        if (sub) {
          num_patients = length(risks_to_test[[1]])
          status = object[[1]]$outcome$status
          num_events = sum(data_to_test[[status]] == indicator)
          subtitle = paste0("#at risk=", num_patients,",#that undergo event=", num_events)
          graphics::title(sub = substitute(paste(italic(subtitle))))
        }
      }


      outlist[[t]] <- x

    }
    if (regression_values) {
      for (i in 1:NF){
        name = names(object)[i]
        rownames(reg_values_list[[name]]) <- NULL
        colnames(reg_values_list[[name]]) <- c("LM", "Intercept", "Slope")
      }
      outlist[["regression_values"]] <- reg_values_list
    }
    return(outlist)
  }
