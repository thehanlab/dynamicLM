#' Methods (AUCt, Brier) to score the predictive performance of dynamic risk markers from LM super models
#'
#' @param preds A named list of prediction models, where allowed entries are outputs from predLMrisk
#' @param formula A survival or event history formula. The left hand side is used to compute the expected event status.
#'   If none is given, it is obtained from the prediction object.
#' @param metrics  Character vector specifying which metrics to apply. Choices are "auc" and "brier". Case matters.
#' @param cause Cause of interest if considering competing risks
#' @param tLM  Landmark times for which scores must be given. These must be a subset of LM times used during the prediction
#' @param unit Time unit for window of prediction, e.g., "year", "month", etc. Used for printing results.
#' @param split.method Method for cross-validation. Right now, as in riskRegression, the only option is bootcv.
#' @param B The number of bootstap steps for cross-validation.
#' @param ... Additional arguments to pass to Score (riskRegression package)
#'
#' @return An object of class "LMScore", which has components:
#'   - auct: dataframe containing time-dependent auc information if "auc" was a metric
#'   - briert: dataframe containing time-dependent brier score if "brier" was a metric
#' @details See the Github for example code
#' @export
#'
LMScore <- function(preds,formula,metrics=c("auc","brier"),cause,tLM, unit, split.method, B,...){

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package \"data.table\" must be installed to use function LMScore.", call. = FALSE)}
  if (!requireNamespace("riskRegression", quietly = TRUE)) {
    stop("Package \"riskRegression\" must be installed to use function LMScore.", call. = FALSE)}
  if(missing(unit)) unit = "year"

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

  auct <- data.table::data.table()
  briert <- data.table::data.table()
  for(t in 1:length(times)){
    tLM = times[t]

    idx = preds[[1]]$preds$LM == tLM
    data_to_test = preds[[1]]$data[idx, ]
    risks_to_test = lapply(preds, function(p) p$preds$risk[idx])

    if (nrow(data_to_test)!=length(risks_to_test[[1]])){
      stop("nrow(data_to_test)!=length(risks_to_test)")
    }

    score_t <- riskRegression::Score(risks_to_test,
                     formula=stats::as.formula(formula),
                     data=data_to_test,
                     metrics=metrics,
                     cause=cause,
                     times=c(tLM+w),
                     split.method=split.method,
                     B=B,
                     ...)

    if ("auc" %in% metrics) auct <- rbind(auct, cbind(tLM,score_t$AUC$score))
    if ("brier" %in% metrics) briert <- rbind(briert, cbind(tLM,score_t$Brier$score))

  }
  outlist = list(auct=auct,briert=briert,w=w,unit=unit)
  class(outlist) = "LMScore"
  return(outlist)
}
