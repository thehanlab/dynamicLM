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
plot.dynamicLM <- function(object, covars, conf_int = TRUE, cause, end_time,
                           logHR = TRUE, extend = FALSE, silence = FALSE,
                           xlab = "LM time", ylab, ylim, main, ...) {
  fm = object$model

  if(conf_int) {
    if (!requireNamespace("msm", quietly = TRUE))
      stop("Package \"msm\" must be installed to use this function.",
           call. = FALSE)
  }

  if (missing(covars)) covars <- object$lm_covs
  if (missing(main)) {
    if (is.null(names(covars))) main <- covars
    else main <- names(covars)
  } else {
    if (length(main) != length(covars))
      stop("# of titles given must equal the # of covariates to plot")
  }
  if (missing(ylab)) {
    if (logHR) ylab <- "log HR"
    else ylab <- "HR"
  }
  if (missing(end_time)) {
    end_time <- object$end_time

  } else if (end_time > object$end_time && !extend) {
    if (!silence) message(paste0("NOTE: arg end_time (=",end_time,") is later than the last LM used in model fitting (=",object$end_time,")",
                                 "\nand has been set back to the last LM used in model fitting. (=",object$end_time,")",
                                 "\nIf you wish to still plot until ",end_time, ", set arg extend=T but note that results after time ",object$end_time," may be unreliable."))
    end_time <- object$end_time

  } else if (end_time > object$end_time && extend) {
    if (!silence) warning(paste0("NOTE: arg end_time (=",end_time,") is later than the last LM used in model fitting (=",object$end_time,")",
                                 "\nResults after time ",object$end_time," may be unreliable."))
  }

  if (object$type == "coxph") {
    if (!missing(cause)) stop("no cause should be input for coxph supermodels.")
    bet <- fm$coefficients
    func_covars <- object$func_covars
    if (conf_int) sig <- stats::vcov(fm)

  } else if (object$type == "CauseSpecificCox" | object$type == "CSC") {
    if (missing(cause)) cause <- as.numeric(fm$theCause)
    if (length(cause) > 1)
      stop(paste0("Can only predict one cause. Provided are: ",
                  paste(cause, collapse = ", "), sep = ""))

    bet <- fm$models[[cause]]$coefficients
    func_covars <- object$func_covars

    if(conf_int) sig <- stats::vcov(fm$models[[cause]])
  }
  if (missing(ylim)) set_ylim <- TRUE

  t <- seq(0, end_time, by = 0.1)

  for (i in 1:length(covars)){
    idx <- startsWith(names(bet), covars[i])
    bet_var <- bet[idx]
    HR <- sapply(t, function(x) { # eval HR over times x in t
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
    if (!logHR) {
      HR <- exp(HR)
    }

    if (set_ylim) ylim <- c(min(HR), max(HR))
    if(conf_int) {
      if (logHR) {
        se <- sapply(t, find_se_log, bet_var, sig[idx, idx], func_covars)
      } else {
        se <- sapply(t, find_se, bet_var, sig[idx, idx], func_covars)
      }
      lower <- HR - 1.96 * se
      upper <- HR + 1.96 * se
      if (set_ylim) {
        ylim[1] <- min(min(lower), ylim[1])
        ylim[2] <- max(max(upper), ylim[1])
      }
    }

    plot(t, HR, xlab = xlab, ylab = ylab, main = main[i], type = "l",
         ylim = ylim, ...)
    if (logHR) {
      graphics::lines(t, rep(0, end_time / 0.1 + 1), col = "grey")
    } else {
      graphics::lines(t, rep(1, end_time / 0.1 + 1), col = "grey")
    }
    if(conf_int) {
      graphics::lines(t, lower, lty = 2)
      graphics::lines(t, upper, lty = 2)
    }
  }
}


#' Plot an object output from [score()]: plot the landmark and time-dependent
#' Brier and/or AUC of dynamic landmark supermodels.
#'
#' @param object An object of class "LMScore" output from [score()]
#' @param metrics One or both of "auc" and "brier"
#' @param se Boolean, default TRUE. To include point wise confidence intervals.
#' @param xlab,ylab,y,x,pch,ylim,xlim graphical parameters
#' @param ...
#'
#' @export
#'
plot.LMScore <- function(object, metrics, se = TRUE, xlab, ylab, x, pch, ylim,
                         xlim, ...) {
  if (missing(metrics)) {
    metrics <- c()
    if (!is.null(object$auct)) metrics <- c("auc")
    if (!is.null(object$briert)) metrics <- c(metrics, "brier")
  }

  if (missing(xlab))
    xlab <- "Landmark Time (t)"
  if (missing(pch))
    pch <- 19

  set_ylab <- FALSE
  if (missing(ylab))
    set_ylab <- TRUE
  set_x <- FALSE
  if (missing(x))
    set_x <- TRUE

  plot.metric <- function(df, metric, loc, ylim, xlim) {
    if (set_ylab) ylab <- paste0(metric, "(t, t + ", object$w, ")")

    num_models = length(unique(df$model))
    model_names = df$model[1:num_models]
    tLM = df$tLM
    metric = df[[metric]]
    upper = df[["upper"]]
    lower = df[["lower"]]
    models = df$model

    if (missing(ylim)) {
      ylim = c(min(lower), max(upper))
    }
    if (missing(xlim)) {
      xlim = c(min(tLM), max(tLM))
    }

    plot(tLM, metric, col = models, pch = pch, xlab = xlab, ylab = ylab,
         ylim = ylim, xlim = xlim, ...)

    for (i in 1:num_models){
      idx = models == model_names[i]
      lines(tLM[idx], metric[idx], col = models[idx], pch = pch)
      if (se) {
        lines(tLM[idx], upper[idx], col = models[idx], pch = pch, lty = 2)
        lines(tLM[idx], lower[idx], col = models[idx], pch = pch, lty = 2)
      }
    }
    legend(loc, legend = model_names, col = models, pch = pch, bty = "n")
  }

  if ("auc" %in% metrics) {
    if (is.null(object$auct)) {
      warning("AUC was not set as a metric when calling score() No results to plot. Either call score() again with auc as a metric or do not include it as a metric here.")
    } else {
      if(set_x) x <- "topright"
      plot.metric(object$auct, "AUC", x, ylim, xlim)
    }
  }
  if ("brier" %in% metrics) {
    if (is.null(object$briert)) {
      warning("Brier was not set as a metric when calling score() No results to plot. Either call score() again with auc as a metric or do not include it as a metric here.")
    } else {
      if (set_x) x <- "bottomright"
      plot.metric(object$briert, "Brier", x, ylim, xlim)
    }
  }
}


#' Plot an object output from [calplot()]: plot the calibration plots.
#'
#' @param object An object of class "LMcalibrationPlot" output from [calplot()]
#' @param ...
#'
#' @export
#'
plot.LMcalibrationPlot <- function(object) {
  for (i in 1:length(object)) {
    plot(object[[i]])
  }
}


#' Plot coefficients from an object created by calling `penLM`, analogous to
#' plotting from `glmnet`
#'
#' As in the glmnet package, produces a coefficient profile plot of the
#' coefficient paths.
#'
#' @details
#' If the model is a survival model (i.e., no competing risks), then the output
#' is the same as a call to `glmnet` would produce. For competing risks, the
#' default is only to plot the coefficient profile plot for the cause of
#' interest (first cause) Further events can be examined by setting
#' `all_causes = TRUE`.
#'
#' @param x a fitted `penLM` object
#' @param xvar As in `glmnet`: "What is on the X-axis. "`norm`" plots against the
#'   L1-norm of the coefficients, "`lambda`" against the log-lambda sequence,
#'   and "`dev`" against the percent deviance explained."
#' @param all_causes if penLM fit a cause-specific Cox model, set TRUE to plot
#'   coefficient profile plots for each model.
#' @param silent Set TRUE to hide messages.
#' @param label Set TRUE to label the curves by variable index numbers.
#' @param \dots additional graphical parameters
#' @export
plot.penLM <- function(x, xvar = "norm", all_causes = FALSE, silent = FALSE,
                       label = FALSE, ...) {
  num_causes <- length(x)
  if (all_causes) {
    for (i in 1:num_causes){
      plot(x[[i]], xvar, label, ...)
    }
  } else {
    plot(x[[1]], xvar, label, ...)
    if (length(x) > 1 & !silent)
    message("\n (To print plot paths for remaining cause-specific models, call plot with argument all_causes = TRUE)\n")
  }
}


#' Plot cross-validation curve created by `cv.penLM`, analogous to plotting from
#' `cv.glmnet`
#'
#' The cross-validation curve is plotted as a function of the lambda values
#' used. Upper and lower standard deviation is plotted too.
#'
#' @details
#' If the model is a survival model (i.e., no competing risks), then the output
#' is the same as a call to `cv.glmnet` would produce. For competing risks, the
#' default is only to plot the cross-validation curve for the cause of
#' interest (first cause) Further events can be examined by setting
#' `all_causes = TRUE`.
#'
#' @param x a fitted `cv.penLM` object
#' @param all_causes if penLM fit a cause-specific Cox model, set TRUE to plot
#'   coefficient profile plots for each model.
#' @param silent Set TRUE to hide messages.
#' @param sign.lambda Plot against `log(lambda)` (default) or its negative
#'   if set to -1.
#' @param se.bands Logical. If TRUE, shading is produced to show stand-error
#'   bands. Defaults to TRUE.
#' @param \dots additional graphical parameters
#' @export
plot.cv.penLM <- function(x, all_causes = FALSE, silent = FALSE, label = FALSE,
                          sign.lambda = 1, se.bands = TRUE, ...) {
  num_causes <- length(x)
  if (all_causes) {
    for (i in 1:num_causes){
      plot(x[[i]], sign.lambda, se.bands, ...)
    }
  } else {
    plot(x[[1]], sign.lambda, se.bands, ...)
    if (length(x) > 1 && !silent)
      message("\n (To print plot paths for remaining cause-specific models, call plot with argument all_causes = TRUE)\n")
  }
}



#' Generic function to plot coefficients
#'
#' Can plot positive and negative coefficients in two separate plots or the
#' same. X-axes are the same if separate plots are used.
#'
#' @param coefs (Named) Vector of coefficients
#' @param single.plot Logical, defaults to FALSE. A single plot for both
#'   positive and negative coefficients, or two separate plots.
#' @param max_coefs Default is 10. The maximum number of coefficients to plot.
#'   Can be set to NULL if plotting all coefficients is desired.
#' @param ... Additional arguments to barplot.
#' @export
plot.coefs <- function(coefs, single.plot, max_coefs, ...) {
  pos_coefs <- sort(coefs[coefs>0], decreasing = T)
  neg_coefs <- (-sort(-(coefs[coefs<0]), decreasing = T))

  if (!single.plot){
    xmax <- 1.1 * max(c(pos_coefs, -neg_coefs), na.rm=T)
    ymax <- max(
      min(length(pos_coefs), max_coefs),
      min(length(neg_coefs), max_coefs)
    )
    if (length(neg_coefs)>0){
      neg_coefs <- neg_coefs[1:min(length(neg_coefs),max_coefs)]
      barplot(neg_coefs,
              col="blue",
              names.arg = names(neg_coefs),
              horiz = TRUE, las=1, xpd = F,
              xlim=c(-xmax, 0),
              ylim=c(0,ymax),
              width = 0.8,
              xlab = "Value")
    }
    if (length(pos_coefs)>0) {
      pos_coefs <- pos_coefs[1:min(length(pos_coefs),max_coefs)]
      barplot(pos_coefs,
              col = "blue",
              names.arg = names(pos_coefs),
              horiz = TRUE, las=1, xpd = F,
              xlim=c(0, xmax),
              ylim=c(0,ymax),
              width = 0.8,
              xlab = "Value")
    }
  }

  else {
    if (length(pos_coefs)>0)
      pos_coefs <- pos_coefs[1:min(length(pos_coefs),max_coefs)]
    if (length(neg_coefs)>0)
      neg_coefs <- neg_coefs[1:min(length(neg_coefs),max_coefs)]

    barplot(c(pos_coefs, neg_coefs),
            col="blue",
            names.arg = c(names(pos_coefs),names(neg_coefs)),
            horiz = TRUE, las=1, xpd = F,
            xlab = "Value", ...)
  }
}

#' Plot the coefficients of a penalized Cox supermodel
#'
#' Can plot positive and negative coefficients in two separate plots or the
#' same. X-axes are the same if separate plots are used.
#'
#' @param object a penalized Cox supermodel - created by calling `fitLM` on an
#'   object created from `penLM`/`cv.pen`
#' @param single.plot Logical, defaults to FALSE. A single plot for both
#'   positive and negative coefficients, or two separate plots.
#' @param max_coefs Default is 10. The maximum number of coefficients to plot.
#'   Can be set to NULL if plotting all coefficients is desired.
#' @param ... Additional arguments to barplot.
#' @export
plot.penLMcoxph <- function(
    object,
    single.plot = FALSE,
    max_coefs = 10,
    ...
) {
  coefs <- object$model$coefficients
  plot.coefs(coefs, single.plot, max_coefs, ...)
}

#' Plot the coefficients of a penalized Cause-specific Cox supermodel
#'
#' @param object a penalized cause-specific Cox supermodel - created by calling
#'   `fitLM` on an object created from `penLM`/`cv.pen`
#' @param single.plot Logical, defaults to FALSE. A single plot for both
#'   positive and negative coefficients, or two separate plots.
#' @param max_coefs Default is 10. The maximum number of coefficients to plot.
#'   Can be set to NULL if plotting all coefficients is desired.
#' @param all_causes Logical, default is FALSE. Plot coefficients for all
#'   cause-specific models.
#' @param ... Additional arguments to barplot.
#' @export
plot.penLMCSC <- function(
    object,
    single.plot = FALSE,
    max_coefs = 10,
    all_causes = FALSE,
    ...
) {
  if (!all_causes) {
    coefs <- object$model$models[[1]]$coefficients
    plot.coefs(coefs, single.plot, max_coefs, ...)
  }
  else {
    lapply(1:length(object$model),
           function(i) plot.coefs(object$model$models[[i]]$coefficients, single.plot, max_coefs, ...))
  }
}

#' Plot an object output from `LMScore`: plot the landmark and time-dependent
#' Brier and/or AUC of dynamic landmark supermodels.
#'
#' @param object An object of class "LMScore" output from `LMScore`
#' @param metrics One or both of "auc" and "brier"
#' @param xlab x label
#' @param ylab y label
#' @param x legend location
#' @param pch size of points
#' @param ...
#'
#' @export
#'
plot.LMScore <- function(object, metrics, xlab, ylab, x, pch, ...){
  if (missing(metrics)){
    metrics <- c()
    if (!is.null(object$auct)) metrics <- c("auc")
    if (!is.null(object$briert)) metrics <- c(metrics, "brier")
  }

  if (missing(xlab))
    xlab <- "Landmark Time (tLM)"
  if (missing(pch))
    pch <- 19

  set_ylab <- F
  if (missing(ylab))
    set_ylab <- T
  set_x <- F
  if (missing(x))
    set_x <- T
  plot.metric <- function(df, metric, loc){
    if (set_ylab) ylab <- metric

    num_models = length(unique(df$model))
    model_names = df$model[1:num_models]
    tLM = df$tLM
    metric = df[[metric]]
    models = df$model

    plot(tLM, metric, col = models, pch = pch, xlab=xlab, ylab=ylab, ...)

    for (i in 1:num_models){
      idx = models == model_names[i]
      lines(tLM[idx], metric[idx], col = models[idx], pch = pch)
    }
    legend(loc, legend = model_names, col = models, pch = pch, bty = "n")
  }

  if ("auc" %in% metrics){
    if (is.null(object$auct)){
      warning("AUC was not set as a metric when calling LMScore. No results to plot. Either call LMScore again with auc as a metric or do not include it as a metric here.")
    }
    else {
      if(set_x) x <- "topright"
      plot.metric(object$auct, "AUC", x)
    }
  }
  if ("brier" %in% metrics){
    if (is.null(object$briert)){
      warning("Brier was not set as a metric when calling LMScore. No results to plot. Either call LMScore again with auc as a metric or do not include it as a metric here.")
    }
    else {
      if(set_x) x <- "bottomright"
      plot.metric(object$briert, "Brier", x)
    }
  }
}
