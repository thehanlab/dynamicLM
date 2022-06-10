#' Plots the dynamic hazard ratio of a cox or CSC supermodel
#'
#' @param supermodel An object of class "LMcoxph" or "LMCSC", i.e. a fitted supermodel
#' @param covars Vector or list of strings indicating the variables to plot the HR of
#'   (note these must be given without time interaction label, for e.g., as in LMcovars).
#' @param CI Include confidence intervals or not, default is TRUE
#' @param cause Cause of interest if considering competing risks
#' @param end_time Final time point to plot HR, defaults to the last landmark point used in model fitting.
#' @param extend Argument to allow for HR to be plot at landmark times that are later than the LMs used in model fitting.
#'   Default is FALSE. If set to TRUE, the HR may be unreliable.
#' @param silence silence the warning message when end_time > LMs used in fitting the model
#' @param xlab x label for the plots
#' @param ylab y label for the plots
#' @param ylim y limit for the plots
#' @param main Vector of strings indicating the title of each plot. Must be in the same order as covars.
#' @param logHR Boolean, if true plots the log of the hazard ratio, if false plots the hazard ratio
#' @param ... Additional arguments passed to plot
#'
#' @return Plots for each variable in covars showing the dynamic hazard ratio
#' @details See the Github for example code
#' @export
#'
plot_dynamic_HR <- function(supermodel, covars, CI=T, cause, end_time, extend=F, silence=F,
                            xlab="LM time", ylab, ylim, main, logHR=T, ...){
  fm = supermodel$model

  if(CI){
    if (!requireNamespace("msm", quietly = TRUE)) {
      stop("Package \"msm\" must be installed to use this function.", call. = FALSE)}
  }

  if (missing(covars)){ covars <- supermodel$LMcovars }
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
    end_time <- supermodel$end_time

  } else if (end_time > supermodel$end_time & !extend){
    if (!silence) message(paste0("NOTE: arg end_time (=",end_time,") is later than the last LM used in model fitting (=",supermodel$end_time,")",
                                 "\nand has been set back to the last LM used in model fitting. (=",supermodel$end_time,")",
                                 "\nIf you wish to still plot until ",end_time, ", set arg extend=T but note that results after time ",supermodel$end_time," may be unreliable."))
    end_time <- supermodel$end_time

  } else if (end_time > supermodel$end_time & extend){
    if (!silence) warning(paste0("NOTE: arg end_time (=",end_time,") is later than the last LM used in model fitting (=",supermodel$end_time,")",
                                 "\nResults after time ",supermodel$end_time," may be unreliable."))
  }

  if (supermodel$type == "coxph") {
    if (!missing(cause)) {stop("no cause should be input for coxph supermodels.")}
    bet <- fm$coefficients
    func_covars <- supermodel$func_covars
    if(CI){ sig <- stats::vcov(fm) }

  } else if (supermodel$type == "CauseSpecificCox" | supermodel$type == "CSC") {
    if (missing(cause)) { cause <- as.numeric(fm$theCause) }
    if (length(cause) > 1) stop(paste0("Can only predict one cause. Provided are: ", paste(cause, collapse = ", "), sep = ""))

    bet <- fm$models[[cause]]$coefficients
    func_covars <- supermodel$func_covars

    if(CI){ sig <- stats::vcov(fm$models[[cause]]) }
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
    if(CI){
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
    if(CI){
      graphics::lines(t, lower, lty=2)
      graphics::lines(t, upper, lty=2)
    }
  }
}
