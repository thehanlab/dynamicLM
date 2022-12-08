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
plot.penLM <- function(x, xvar="norm", all_causes = FALSE, silent = FALSE,
                       label=FALSE, ...){
  num_causes <- length(x)
  if (all_causes){
    for (i in 1:num_causes){
      plot(x[[i]], xvar, label, ...)
    }
  } else {
    plot(x[[1]], xvar, label, ...)
    if (length(x) > 1 & !silent) message("\n (To print plot paths for remaining cause-specific models, call plot with argument all_causes = TRUE)\n")
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
plot.cv.penLM <- function(x, all_causes = FALSE, silent = FALSE, label=FALSE,
                          sign.lambda=1, se.bands=TRUE, ...){
  num_causes <- length(x)
  if (all_causes){
    for (i in 1:num_causes){
      plot(x[[i]], sign.lambda, se.bands, ...)
    }
  } else {
    plot(x[[1]], sign.lambda, se.bands, ...)
    if (length(x) > 1 & !silent) message("\n (To print plot paths for remaining cause-specific models, call plot with argument all_causes = TRUE)\n")
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
    ymax <- max(length(pos_coefs), length(neg_coefs))
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
            xlab = "Value")
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
