#' Plots the dynamic log-hazard ratio of a cox or CSC supermodel
#'
#' @param object An object of class "LMcoxph" or "LMCSC", i.e. a fitted
#'   supermodel
#' @param covars Vector or list of strings indicating the variables to plot
#'   (note these must be given without time interaction label, for e.g., as in
#'   the argument `lm_covs` in [add_interactions()]).
#' @param conf_int Include confidence intervals or not, default is TRUE
#' @param cause Cause of interest if considering competing risks
#' @param end_time Final time point to plot HR, defaults to the last landmark
#'   point used in model fitting.
#' @param logHR Boolean, if true plots the log of the hazard ratio, if false
#'   plots the hazard ratio. Default is TRUE.
#' @param extend Argument to allow for HR to be plot at landmark times that are
#'   later than the LMs used in model fitting.
#'   Default is FALSE. If set to TRUE, the HR may be unreliable.
#' @param silence silence the warning message when end_time > LMs used in
#'   fitting the model
#' @param xlab x label for the plots
#' @param ylab y label for the plots
#' @param ylim y limit for the plots
#' @param main Vector of strings indicating the title of each plot. Must be in
#'   the same order as covars.
#' @param ... Additional arguments passed to plot
#'
#' @return Plots for each variable in covars showing the dynamic hazard ratio
#' @details See our [GitHub](https://github.com/thehanlab/dynamicLM) for example
#'   code
#' @export
#'
plot.dynamicLM <- function(object, covars, conf_int=T, cause, end_time, logHR=T,
                           extend=F, silence=F,
                           xlab="LM time", ylab, ylim, main,  ...){
  fm = object$model

  if(conf_int){
    if (!requireNamespace("msm", quietly = TRUE)) {
      stop("Package \"msm\" must be installed to use this function.", call. = FALSE)}
  }

  if (missing(covars)){ covars <- object$lm_covs }
  if(missing(main)){
    if(is.null(names(covars))) { main <- covars }
    else { main <- names(covars) }
  } else {
    if(length(main)!=length(covars)) stop("# of titles given must equal the # of covariates to plot")
  }
  if (missing(ylab)){
    if(logHR){ylab <- "log HR"}
    else {ylab <- "HR"}
  }
  if (missing(end_time)){
    end_time <- object$end_time

  } else if (end_time > object$end_time & !extend){
    if (!silence) message(paste0("NOTE: arg end_time (=",end_time,") is later than the last LM used in model fitting (=",object$end_time,")",
                                 "\nand has been set back to the last LM used in model fitting. (=",object$end_time,")",
                                 "\nIf you wish to still plot until ",end_time, ", set arg extend=T but note that results after time ",object$end_time," may be unreliable."))
    end_time <- object$end_time

  } else if (end_time > object$end_time & extend){
    if (!silence) warning(paste0("NOTE: arg end_time (=",end_time,") is later than the last LM used in model fitting (=",object$end_time,")",
                                 "\nResults after time ",object$end_time," may be unreliable."))
  }

  if (object$type == "coxph") {
    if (!missing(cause)) {stop("no cause should be input for coxph supermodels.")}
    bet <- fm$coefficients
    func_covars <- object$func_covars
    if(conf_int){ sig <- stats::vcov(fm) }

  } else if (object$type == "CauseSpecificCox" | object$type == "CSC") {
    if (missing(cause)) { cause <- as.numeric(fm$theCause) }
    if (length(cause) > 1) stop(paste0("Can only predict one cause. Provided are: ", paste(cause, collapse = ", "), sep = ""))

    bet <- fm$models[[cause]]$coefficients
    func_covars <- object$func_covars

    if(conf_int){ sig <- stats::vcov(fm$models[[cause]]) }
  }
  if (missing(ylim)) { set_ylim <- T }

  t <- seq(0, end_time, by=0.1)

  for (i in 1:length(covars)){
    idx <- startsWith(names(bet), covars[i])
    bet_var <- bet[idx]
    HR <- sapply(t, function(x){ # eval HR over times x in t
      sum(sapply(1:length(bet_var), function(j){
        var = bet_var[j]
        name = names(bet_var)[j]
        if (name == covars[i]) { return(var) }
        else {
          idx <- as.numeric(sub(".*\\D+", "\\1", name))
          return(func_covars[[idx]](x) * var)
        }
      })) # bet0 + bet1*x + bet2*x^2 + ...
    })
    if(!logHR){
      HR <- exp(HR)
    }

    if (set_ylim) ylim <- c(min(HR),max(HR))
    if(conf_int){
      if(logHR){ se <- sapply(t, find_se_log, bet_var, sig[idx,idx], func_covars) }
      else{ se <- sapply(t, find_se, bet_var, sig[idx,idx], func_covars) }
      lower <- HR - 1.96*se
      upper <- HR + 1.96*se
      if(set_ylim){
        ylim[1] <- min(min(lower), ylim[1])
        ylim[2] <- max(max(upper), ylim[1])
      }
    }

    plot(t, HR, xlab=xlab, ylab=ylab, main=main[i], type="l", ylim=ylim, ...)
    if(logHR){ graphics::lines(t, rep(0,end_time/0.1+1), col="grey") }
    else { graphics::lines(t, rep(1,end_time/0.1+1), col="grey") }
    if(conf_int){
      graphics::lines(t, lower, lty=2)
      graphics::lines(t, upper, lty=2)
    }
  }
}