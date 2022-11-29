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
