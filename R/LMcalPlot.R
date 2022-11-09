#' Calibration plots for dynamic risk prediction landmark models
#'
#' @param object A named list of prediction models, where allowed entries are outputs from predLMrisk or supermodels from fitLM
#' @param times Landmark times for which calibration must be plot. These must be a subset of LM times used during the prediction
#' @param formula A survival or event history formula (Hist(..)). The left hand side is used to compute the expected event status.
#'   If none is given, it is obtained from the prediction object.
#' @param data Data on which to predict if the object is a fitted model. If bootstrapping, must be an object
#'   of class LMdataframe. If not bootstrapping, it can be a dataframe.
#' @param tLM Landmark times corresponding to the patient entries in data. Only required if data is a dataframe.
#'   tLM can be a string (indicating a column in data), a vector of length nrow(data),
#'   or a single value if all patient entries were obtained at the same landmark time.
#' @param ID_col Column name that identifies individuals in data. Only required if bootstrapping.
#' @param split.method Defines the internal validation design as in pec::calPlot. Options are currently none or bootcv
#' @param B Number of times bootstrapping is performed.
#' @param M Subsample size for training in cross-validation. Entries not sampled in the M subsamples are used for validation.
#' @param unit The unit of w, i.e. w-unit prediction ("year","month", etc...). Only used to label the plot.
#' @param cause Cause of interest if considering competing risks
#' @param plot If FALSE, do not plot the results, just return a plottable object. Default is TRUE.
#' @param main Optional title to override default.
#' @param sub If TRUE, add a subheading with the number of individuals at risk, and the number that under the event of interest.
#'   Default is TRUE. Set to FALSE for bootstrapping.
#' @param ... Additional arguments to pass to calPlot
#'
#' @return List of plots of w-year risk, one entry per prediction/landmark time point
#' @details Most errors in plotting occur when a formula is not given. Formulas can look like `Surv(LM,Time,event)~1` / `Surv(LM,Time,event==1)~1` / `Hist(Time,event,LM)~1` / similar...
#'
#'  See the Github for example code on using LMcalPlot in general.
#' @examples
#' \dontrun{
#' par(mfrow=c(2,2),pty="s")
#' outlist = LMcalPlot(list("Model1"=p1),
#'                     unit="month",            # for the title
#'                     times=c(6,12,18,24),       # landmarks at which to provide calibration plots
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
           tLM, ##
           ID_col="ID", ##
           split.method = "none", ##
           B = 1, ##
           M, ##
           unit = "year",
           cause,
           plot = T,
           main,
           sub = T,
           ...) {
    if (!requireNamespace("pec", quietly = T)) {
      stop("Package \"pec\" must be installed to use function LMcalPlot.", call. = F)
    }
    if (!(class(object)=="list")) stop("object must be a named list.")

    # first checks
    split.method <- tolower(split.method)
    if (split.method == "none") B <- 1

    # check naming of objects
    NF <- length(object)
    if (is.null(names(object))) {
      names(object) <- paste("Model", seq(NF))
    } else {
      names(object)[(names(object)=="")] <- paste("Model", seq(NF))[(names(object)=="")]
    }
    names(object) <- make.unique(names(object))

    # check bootstrapping and inputs are coherent
    for (i in 1:NF){
      name = names(object)[[i]]
      # check compatibility of object type with bootstrapping and predictions
      if (class(object[[i]]) == "LMpred"){
        if (split.method != "none") {
          stop(paste("Cannot bootstrap with deterministic risk predictions:", name))
        }
        if (!missing(data)) {
          stop(paste0("Cannot use deterministic predictions of type LMpred (",name,") on new data, data argument should not be specified."))
        }


      } else if (!(class(object[[i]]) %in% c("LMCSC", "LMcoxph"))) {
        stop(paste("all prediction models in object must be of class LMCSC, LMcoxph (i.e., output from fitLM) or LMpred (i.e., output from predLMrisk) but", name, "is of class",class(object[[i]])))


      } else {
        # check if we have data if bootstrapping
        if (missing(data)){
          model.class = class(object[[i]])
          stop(paste("Models of class ", model.class, "do not store original data, data must be specified when bootstrapping."))
        }
        if ((split.method == "bootcv") & (class(data) != "LMdataframe")){
          stop("In order to bootstrap, data must be of class LMdataframe in order to refit models")
        } else if ((split.method == "bootcv") & (class(data) == "LMdataframe")){
          sub <- FALSE

          LMdata <- data
          LMdata$LMdata <- LMdata$LMdata[complete.cases(LMdata$LMdata), ]
          data <- LMdata$LMdata

          LM_col <- LMdata$LM_col

          if(!(ID_col %in% colnames(data))){
            stop(paste("An ID column that is in the data must be provided;", ID_col, "is not a column."))
          }
        }
      }
    }

    # more checks
    w <- object[[1]]$w
    if (NF > 1) {
      for (i in 2:NF) {
        if (object[[i]]$w != w)
          stop("prediction window w is not the same for all prediction models.")
      }
    }
    if (missing(cause))
      cause <- object[[1]]$cause
    if (is.null(cause))
      cause <- 1
    if (missing(formula))
      formula <- object[[1]]$LHS

    # split method
    perform.boot <- F
    if (split.method == "bootcv"){
      if (!missing(tLM)){
        warning("LMdataframe objects have LM columns, tLM is ignored.")
      }
      unique.inds = unique(data[[ID_col]])
      split.method <- riskRegression::getSplitMethod(split.method,B=B,N=length(unique.inds),M=M)
      B <- split.method$B
      split.idx <- split.method$index
      perform.boot <- !is.null(split.idx)
      if (perform.boot) split.idx <- replace(split.idx, seq_along(split.idx), unique.inds[split.idx])# get IDs of the individuals
    }
    else if (!(split.method %in% c("none", "noplan"))) {
      paste("split.method of type", split.method, "is not supported.")
    }


    ### Create data we need: preds, preds_LM, data

    if (class(object[[1]]) == "LMpred") {
      if (NF > 1) {
        for (i in 2:NF) {
          if (class(object[[i]]) != "LMpred")
            stop("All prediction models must either be supermodels or of type LMpred")
        }
      }

      # check data
      preds <- data.frame(lapply(object, function(o) o$preds$risk))
      colnames(preds) <- names(object)
      pred_LMs <- object[[1]]$preds$LM
      data <- object[[1]]$data
      num_preds <- nrow(object[[1]]$preds)

      type <- lapply(object, function(o) ifelse(o$type == "coxph", "coxph", "CSC"))

      if (NF > 1) {
        for (i in 2:NF) {
          if (nrow(object[[i]]$preds) != num_preds)
            stop("number of predictions is not the same for all prediction models.")
          if (!isTRUE(all.equal(object[[i]]$preds$LM, pred_LMs)))
            stop("LM points for individuals across prediction models are not the name.")
        }
      }
    }


    else if (!perform.boot) {
      # TODO: consider including w, extend, silence, complete as args to predLMrisk
      # TODO: consider allowing not LMdata
      # TODO: test what happens if data is a dataframe
      preds = lapply(object, function(o) predLMrisk(o, data, tLM, cause))
      args = match.call()
      args$data = NULL
      args$tLM = NULL
      args$object = preds
      return(eval(args))
    }


    else if (perform.boot){
      # TODO: handle cores
      cores = 1

      pred.list <- parallel::mclapply(1:B,function(b){
        outcome <- object[[1]]$outcome

        id_train_b <- split.method$index[,b]
        id_train_b <- data[[ID_col]] %in% id_train_b

        data_val_b <- data[!id_train_b, ]
        tLMs_b <- data_val_b[[LM_col]]
        outcomes_val_b <- data_val_b[c(outcome$time, outcome$status, LM_col)]

        data_train_b <- LMdata
        data_train_b$LMdata <- data[id_train_b, ]

        preds.b <- do.call("cbind",lapply(1:NF,function(f){
          original_model <- object[[f]]
          args = original_model$args
          args$LMdata = data_train_b
          model.b <- eval(args)

          pred.b <- try(
            # TODO: consider including w, extend, silence, complete as args to predLMrisk
            # TODO: use complete=T
            # TODO: consider allowing not LMdata
            predLMrisk(model.b,newdata=data_val_b,tLM=tLMs_b,cause=model.b$cause, complete=F),
            silent = F
          )
          if (inherits(pred.b, "try-error")){
            rep(NA,NROW(data_val_b))
          } else {
            pred.b$preds$risk
          }
        }))
        colnames(preds.b) <- names(object)
        cbind(outcomes_val_b, preds.b)

      }, mc.cores=cores)

      pred.df <- do.call("rbind",pred.list)
      pred.df <- pred.df[complete.cases(pred.df), ] # TODO: remove/fix
      preds <- pred.df[names(object)]
      pred_LMs <- pred.df[[LM_col]]
      data <- pred.df[c(outcome$time, outcome$status, LM_col)]
      num_preds <- nrow(data)
      type <- lapply(object, function(o) ifelse(class(o) == "LMcoxph", "coxph", "CSC"))

      rm(pred.df)
    }

    # check LM times
    if (missing(tLM))
      times <- unique(pred_LMs)
    else {
      if (all(tLM %in% pred_LMs))
        times <- tLM
      else {
        tLM = paste(tLM, collapse = ",")
        pred_LMs = paste(pred_LMs, collapse = ",")
        stop(paste("arg tLM (= ",tLM,") must be a subset of landmark prediction times (= ",pred_LMs,")"))
      }
    }

    # Plotting parameters
    add_title = T
    if (!missing(main)) {
      add_title = F
    }

    # Calibration
    outlist = list()
    for (t in 1:length(times)) {
      tLM = times[t]

      idx = pred_LMs == tLM
      data_to_test = data[idx,]
      risks_to_test = lapply(1:NF, function(i) preds[idx, i])
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
          num_events = sum(data_to_test[[status]] == cause)
          subtitle = paste0("#at risk=", num_patients,",#that undergo event=", num_events)
          graphics::title(sub = substitute(paste(italic(subtitle))))
        }
      }

      outlist[[t]] <- x
    }
    return(outlist)
  }
