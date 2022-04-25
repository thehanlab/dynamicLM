#' Calibration plots for dynamic risk prediction landmark models
#'
#' @param preds A named list of prediction models, where allowed entries are outputs from predLMrisk
#' @param cause Cause of interest if considering competing risks
#' @param tLM Landmark times for which calibration must be plot. These must be a subset of LM times used during the prediction
#' @param formula A survival or event history formula. The left hand side is used to compute the expected event status.
#' If none is given, it is obtained from the prediction object.
#' @param plot If FALSE, do not plot the results, just return a plottable object.
#' @param main Optional title to override default.
#' @param sub If TRUE, add a subheading with the number of individuals at risk, and the number that under the event of interest
#' @param ... Additional arguments to pass to calPlot
#'
#' @return List of plots of w-year risk, one entry per prediction/landmark time point
#' @import prodlim
#' @export
#'
LMcalPlot <- function(preds,cause,tLM,formula,plot=T,main,sub=T,...){
  # TODO: Add option to show AUCt, Brier on plot

  if (!requireNamespace("pec", quietly = TRUE)) {
    stop("Package \"pec\" must be installed to use function LMcalPlot.", call. = FALSE)}

  pred_LMs <- unique(preds[[1]]$preds$LM)
  if(missing(tLM)) times <- pred_LMs
  else {
    if(all(tLM %in% pred_LMs)) times <- tLM
    else {
      tLM = paste(tLM, collapse = ",")
      pred_LMs = paste(pred_LMs, collapse = ",")
      stop(paste("arg tLM (= ",tLM,") must be a subset of landmark prediction times (= ",pred_LMs,")"))
    }
  }

  w = preds[[1]]$w
  data = eval(preds[[1]]$data)
  num_preds = nrow(preds[[1]]$preds)

  if (missing(cause)) cause <- preds[[1]]$cause
  if (missing(formula)) formula <- paste0(preds[[1]]$LHS,"~1.")

  if (length(preds)>1){
    for (i in 2:length(preds)){
      if (preds[[i]]$w != w) stop("prediction window w is not the same for all prediction models.")
      if (nrow(preds[[i]]$preds) != num_preds) stop("number of predictions is not the same for all prediction models.")
      if (!isTRUE(all.equal(preds[[i]]$preds$LM,preds[[1]]$preds$LM))) stop("LM points used for each data point is not the same for all prediction models.")
    }
  }

  add_title=T
  if(!missing(main)){add_title=F}

  outlist = list()
  for(t in 1:length(times)){
    tLM = times[t]

    idx = preds[[1]]$preds$LM == tLM
    data_to_test = data[idx, ]
    risks_to_test = lapply(preds, function(p) p$preds$risk[idx])

    if (nrow(data_to_test)!=length(risks_to_test[[1]])){
      stop("nrow(data_to_test)!=length(risks_to_test)")
    }

    x = NULL
    x <- pec::calPlot(
      risks_to_test,
      time=tLM+w,
      formula=stats::as.formula(formula),
      data=data_to_test,
      cause=cause,
      plot=plot,
      ...
    )

    if (plot){
      if(add_title){ graphics::title(main = paste0("Calibration of ",w, "-yr risk \n measured at LM time ",tLM)) }
      else { graphics::title(main=main) }

      if(sub){
        num_patients = length( risks_to_test[[1]] )
        num_events = sum( data_to_test$event==cause)
        subtitle = paste0("#at risk=",num_patients,", #that undergo event=",num_events)
        graphics::title(sub = substitute(paste(italic(subtitle) )))
      }
    }


    outlist[[t]] <- x
  }
  return(outlist)
}
