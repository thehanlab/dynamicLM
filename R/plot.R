#' Plots the dynamic log-hazard ratio of a cox or CSC supermodel
#'
#' @param x A fitted supermodel
#' @param covars Vector or list of strings indicating the variables to plot
#'   (note these must be given without time interaction).
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
#' @param discrete_grid Defaults to 0.1, how to discretize the grid for plotting
#' @param ... Additional arguments passed to plot
#'
#' @return Plots for each variable in covars showing the dynamic hazard ratio
#' @details See our [GitHub](https://github.com/thehanlab/dynamicLM) for example
#'   code
#' @export
#'
#' @examples
#'
#' data(relapse)
#' outcome <- list(time = "Time", status = "event")
#' covars <- list(fixed = c("male", "stage", "bmi"),
#'                varying = c("treatment"))
#' w <- 60; lms <- c(0, 6, 12, 18)
#' lmdata <- stack_data(relapse, outcome, lms, w, covars, format = "long",
#'                      id = "ID", rtime = "T_txgiven")
#' lmdata <- add_interactions(lmdata, func_covars = c("linear", "quadratic"),
#'                            func_lms = c("linear", "quadratic"))
#' formula <- "Hist(Time, event, LM) ~ male + male_LM1 + male_LM2 +
#'             stage + stage_LM1 + stage_LM2 + bmi + bmi_LM1 + bmi_LM2 +
#'             treatment + treatment_LM1 + treatment_LM2 + LM1 + LM2 + cluster(ID)"
#' supermodel <- dynamic_lm(lmdata, as.formula(formula), "CSC", x = TRUE)
#'
#' par(mfrow = c(2, 3))
#' plot(supermodel)
#'
#' par(mfrow = c(1, 2))
#' plot(supermodel,
#'      covars = c("stage", "bmi"), # subset of covariates to plot
#'      logHR = FALSE,              # plot HR instead of log HR
#'      conf_int = FALSE,           # do not plot confidence intervals
#'      main = c("HR of stage", "HR of BMI"))
#'
plot.dynamicLM <- function(x, covars, conf_int = TRUE, cause, end_time,
                           logHR = TRUE, extend = FALSE, silence = FALSE,
                           xlab = "Landmark time", discrete_grid = 0.1,
                           ylab, ylim, main, ...) {
  fm <- x$model

  if (conf_int) {
    if (!requireNamespace("msm", quietly = TRUE))
      stop("Package \"msm\" must be installed to use this function.",
           call. = FALSE)
  }

  if (missing(covars) || is.null(covars)) {
    covars <- x$lm_covs
  } else {
    if (any(!covars %in% x$lm_covs))
      stop(tidymess(paste(
        covars[!covars %in% x$lm_covs][1], "is not a variable in the model.
        Please check if the argument covars is correct.")))
  }

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
    end_time <- x$end_time

  } else if (end_time > x$end_time && !extend) {
    if (!silence)
      message(paste0(
        "NOTE: arg `end_time` (=", end_time,
        ") is later than the last LM used in model fitting (=", x$end_time, ")",
        "\nand has been set back to the last LM used in model fitting. (=",
        x$end_time, ")", "\nIf you wish to still plot until ", end_time,
        ", set arg extend = TRUE but note that results after time ",
        x$end_time, " may be unreliable."
      ))
    end_time <- x$end_time

  } else if (end_time > x$end_time && extend) {
    if (!silence)
      warning(paste0(
        "NOTE: arg end_time (=", end_time,
        ") is later than the last LM used in model fitting (=", x$end_time, ")",
        "\nResults after time ", x$end_time, " may be unreliable."
      ))
  }

  if (x$type == "coxph") {
    if (!missing(cause)) warning("no `cause` should be input for coxph supermodels.")
    bet <- fm$coefficients
    func_covars <- x$func_covars
    if (conf_int) sig <- stats::vcov(fm)

  } else if (x$type == "CauseSpecificCox" | x$type == "CSC") {
    if (missing(cause)) cause <- as.numeric(fm$theCause)
    if (length(cause) > 1)
      stop(paste0("Can only predict one `cause.` Provided are: ",
                  paste(cause, collapse = ", "), sep = ""))

    bet <- fm$models[[cause]]$coefficients
    func_covars <- x$func_covars

    if (conf_int) sig <- stats::vcov(fm$models[[cause]])
  }
  set_ylim <- FALSE
  if (missing(ylim)) set_ylim <- TRUE

  t <- seq(0, end_time, by = discrete_grid)

  for (i in seq_along(covars)){
    idx <- startsWith(names(bet), covars[i])
    bet_var <- bet[idx]
    HR <- sapply(t, function(x) { # eval HR over times x in t
      sum(sapply(seq_along(bet_var), function(j) {
        var <- bet_var[j]
        name <- names(bet_var)[j]
        if (name == covars[i]) {
          return(var)
        } else {
          idx <- as.numeric(sub(".*_LM(\\d)$", "\\1", name))
          return(func_covars[[idx]](x) * var)
        }
      })) # bet0 + bet1*x + bet2*x^2 + ...
    })

    if (conf_int) {
      se <- sapply(t, find_se_log, bet_var, sig[idx, idx], func_covars)

      lower <- HR - 1.96 * se
      upper <- HR + 1.96 * se
      if (!logHR) {
        HR <- exp(HR)
        lower <- exp(lower)
        upper <- exp(upper)
      }
      if (set_ylim) {
        ylim <- c(min(HR), max(HR))
        ylim[1] <- min(min(lower), ylim[1])
        ylim[2] <- max(max(upper), ylim[1])
      }
    } else {
      if (!logHR) HR <- exp(HR)
      if (set_ylim) ylim <- c(min(HR), max(HR))
    }

    plot(t, HR, xlab = xlab, ylab = ylab, main = main[i], type = "l",
         ylim = ylim, ...)
    if (logHR) {
      graphics::lines(t, rep(0, end_time / discrete_grid + 1), col = "grey")
    } else {
      graphics::lines(t, rep(1, end_time / discrete_grid + 1), col = "grey")
    }
    if(conf_int) {
      graphics::lines(t, lower, lty = 2)
      graphics::lines(t, upper, lty = 2)
    }
  }
}


#' Plot the non-zero coefficients of a penalized Cox landmark supermodel or
#' the dynamic log-hazard ratios
#'
#' Can plot positive and negative coefficients in two separate plots or the
#' same. X-axes are the same if separate plots are used.
#'
#' @param x a penalized Cox supermodel - created by calling [dynamic_lm()] on an
#'   object created from [pen_lm()] or [cv.pen_lm()].
#' @param single_plot Logical, defaults to TRUE. A single plot for both
#'   positive and negative coefficients, or two separate plots.
#' @param max_coefs Default is to plot all coefficients. If specified, gives the
#'   maximum number of coefficients to plot.
#' @param col Fill color for the barplot.
#' @param xlab x-axis Label
#' @param HR Plot the hazard ratio? Default is FALSE. See [plot.dynamicLM()]
#'   for additional arguments.
#' @param covars If HR is TRUE, a vector or list of strings indicating the
#'   variables to plot (note these must be given without time interaction).
#'   Defaults to all non-zero variables.
#' @param ... Additional arguments to barplot or to [plot.dynamicLM()].
#'
#' @details If plotting the log hazard ratios, check [plot.dynamicLM()] to
#'   see further arguments.
#'
#' @export
plot.penLMcoxph <- function(x, single_plot = TRUE, max_coefs = NULL,
                            col = "blue", xlab = "Coefficient value",
                            HR = FALSE, covars = NULL, ...) {
  coefs <- x$model$coefficients
  if (HR) {
    tmp <- x
    class(tmp) <- "dynamicLM"
    if (is.null(covars))
      covars <- intersect(x$lm_covs, names(coefs[coefs > 0]))
    plot(tmp, covars = covars, conf_int = FALSE, ...)
  } else {
    plot.coefs(coefs, single_plot, max_coefs, col, xlab, ...)
  }
}


#' Plot the non-zero coefficients of a penalized cause-specific Cox landmark
#'   supermodel or the dynamic log-hazard ratios
#'
#'Can plot positive and negative coefficients in two separate plots or the
#' same. X-axes are the same if separate plots are used. If plotting the log
#' hazard ratios, check [plot.dynamicLM()] to see further arguments.
#'
#' @param x a penalized cause-specific Cox supermodel - created by calling
#'   [dynamic_lm()] on an object created from [pen_lm()] or [cv.pen_lm()].
#' @param single_plot Logical, defaults to TRUE. A single plot for both
#'   positive and negative coefficients, or two separate plots.
#' @param max_coefs Default is to plot all coefficients. If specified, gives the
#'   maximum number of coefficients to plot.
#' @param all_causes Logical, default is FALSE. Plot coefficients for all
#'   cause-specific models.
#' @param HR Plot the hazard ratio? Default is FALSE. See [plot.dynamicLM()]
#'   for additional arguments.
#' @param covars If HR is TRUE, a vector or list of strings indicating the
#'   variables to plot (note these must be given without time interaction).
#'   Defaults to all variables.
#' @param col Fill color for the barplot.
#' @param xlab x-axis Label

#' @param ... Additional arguments to barplot or to [plot.dynamicLM()].
#'
#' @details If plotting the log hazard ratios, check [plot.dynamicLM()] to
#'   see further arguments.
#' @export
plot.penLMCSC <- function(x,
                          single_plot = TRUE,
                          max_coefs = NULL,
                          all_causes = FALSE,
                          HR = FALSE,
                          covars = NULL,
                          col = "blue",
                          xlab = "Coefficient value",
                           ...) {
  if (HR) {
    tmp <- x
    class(tmp) <- "dynamicLM"
    if (is.null(covars)) {
      coefs <- x$model$models[[1]]$coefficients
      covars <- intersect(x$lm_covs, names(coefs[coefs > 0]))
    }
    plot(tmp, covars = covars, conf_int = FALSE, ...)
  } else {
    if (!all_causes) {
      coefs <- x$model$models[[1]]$coefficients
      plot.coefs(coefs, single_plot, max_coefs, col, xlab, ...)
    } else {
      lapply(seq_along(x$model$models),
             function(i) {
               plot.coefs(x$model$models[[i]]$coefficients,
                          single_plot, max_coefs, col, xlab, ...)
             })
    }
  }
}

# ----------------------------------------------------------------------------
# plot_metric: Helper function for plot.LMScore to plot time-dependent metrics
# ----------------------------------------------------------------------------
plot_metric <- function(df, metric, col, loc, ylim, xlim, main, ylab, xlab, w,
                        font.main, pch, length, legend, legend.title, cex,
                        contrasts = FALSE, ...) {
  model_names <- unique(df$model)
  tLM <- df$tLM
  metric_values <- df[[metric]]
  upper <- df[["upper"]]
  lower <- df[["lower"]]

  cols <- if (is.null(col)) df$model else col[as.numeric(as.factor(df$model))]
  se <- !is.null(lower)

  if (missing(ylim)) {
    ylim <- range(c(metric_values, lower, upper), na.rm = TRUE)
    diff <- abs(0.2 *  (ylim[2] - ylim[1]))
    ylim <- c(ylim[1] - diff, ylim[2] + diff)
  }
  if (missing(xlim)) {
    xlim <- c(min(tLM), max(tLM))
  }
  if (missing(main)) {
    main <- switch(metric,
                   "AUC" = "Time-dependent AUC",
                   "delta.AUC" = "Difference in Time-dependent AUC",
                   "Brier" = "Time-dependent Brier Score",
                   "delta.Brier" = "Difference in Time-dependent Brier Score")
  }
  if (missing(ylab)) {
    ylab <- switch(metric,
                   "AUC" = paste0("AUC(t, t + ", w, ")"),
                   "delta.AUC" = paste0("Difference in AUC(t, t + ", w, ")"),
                   "Brier" = paste0("BS(t, t + ", w, ")"),
                   "delta.Brier" = paste0("Difference in BS(t, t + ", w, ")"))
  }

  plot(tLM, metric_values, col = cols, pch = pch, xlab = xlab, ylab = ylab,
       ylim = ylim, xlim = xlim, main = main, font.main = font.main, ...)
  if (contrasts) graphics::abline(h = 0, col = "black", lwd = 1, lty = 2)

  for (model in model_names){
    idx <- df$model == model
    graphics::lines(tLM[idx], metric_values[idx], col = cols[idx], pch = pch)
    if (se) {
      graphics::arrows(tLM[idx], lower[idx], tLM[idx], upper[idx],
                       col = cols[idx], length = length, angle = 90, code = 3)
    }
  }
  if (legend) {
    graphics::legend(loc, legend = model_names,
                     col = cols[match(model_names, df$model)],
                     pch = pch, bty = "n", cex = cex, title = legend.title)
  }
}

# ----------------------------------------------------------------------------
# plot_summary: Helper function for plot.LMScore to plot time-dependent metrics
# ----------------------------------------------------------------------------
plot_summary <- function(df, col_name, col, ylim, xlab, ylab, main,
                         cex, font.main, length, se, contrasts = FALSE, ...) {

  cols <- if (is.null(col)) df$model else col[as.numeric(as.factor(df$model))]

  if (missing(main)) {
    main <- switch(
      col_name,
      "AUC" = latex2exp::TeX("Summary $\\bar{AUC}_{w}$"),
      "delta.AUC" = latex2exp::TeX("Difference in Summary $\\bar{AUC}_{w}$"),
      "Brier" = latex2exp::TeX("Summary $\\bar{BS}_{w}$"),
      "delta.Brier" = latex2exp::TeX("Difference in Summary $\\bar{BS}_{w}$"))
  }
  if (missing(ylab)) {
    ylab <- switch(
      col_name,
      "AUC" = latex2exp::TeX("$\\bar{AUC}_{w}$"),
      "delta.AUC" = latex2exp::TeX("$\\Delta\\bar{AUC}_{w}$"),
      "Brier" = latex2exp::TeX("$\\bar{BS}_{w}$"),
      "delta.Brier" = latex2exp::TeX("$\\Delta\\bar{BS}_{w}$"))
  }
  if (missing(xlab)) xlab <- ""
  if (xlab == "Landmark Time (t)") xlab <- ""

  if (missing(ylim)) {
    ylim <- range(c(df[[col_name]], df[["lower"]], df[["upper"]]), na.rm = TRUE)
    diff <- abs(0.2 * (ylim[2] - ylim[1]))
    ylim <- c(ylim[1] - diff, ylim[2] + diff)
  }

  if (is.null(df[["lower"]])) se <- FALSE

  # Plot
  mid <- graphics::barplot(df[[col_name]], plot = FALSE)
  graphics::barplot(
    df[[col_name]], col = cols, cex.axis = cex, cex.names = cex,
    names.arg = df$model, ylim = ylim, xlab = xlab, ylab = ylab, main = main,
    font.main = font.main, axis.lty = 1, xpd = ylim[1] == 0, ...
    )

  if (se) {
    graphics::arrows(mid, df[["lower"]], mid, df[["upper"]], col = "black",
                     length = length, angle = 90, code = 3)
  }
  return(mid)
}

#' Plot an object output from [score()]: plot the time-dependent and/or summary
#' Brier and/or AUC of landmark supermodels.
#'
#' @param x An object of class "LMScore" output from [score()]
#' @param metrics One or both of "AUC" and "Brier"
#' @param contrasts Plot the difference between metrics. Default is
#'   FALSE and plots the metrics themselves.
#' @param landmarks Plot time-dependent metrics. Default is TRUE.
#' @param summary Plot the summary metric. Default is TRUE.
#' @param se To include point wise confidence intervals. Default is TRUE.
#' @param add_pairwise_contrasts If plotting summary metrics (`summary = TRUE`,
#'   `landmarks = FALSE`) set this argument TRUE to include the p-values of
#'   significant pairwise contrasts. In this case, arguments `pairwise_heights`
#'   and `width` must be set. The argument `cutoff_contrasts` is optional,
#'   specifying the significance cutoff.
#' @param cutoff_contrasts If `add_pairwise_contrasts`, sets the signifance
#'   level of which tests are considered significant (numeric, default is 0.05).
#' @param pairwise_heights If `add_pairwise_contrasts`, sets the height at which
#'   the p-values are plotted. Given as a vector of heights.
#' @param width If `add_pairwise_contrasts`, the width of the ends of the
#'   contrast bars as a numeric value.
#' @param loc Location for legend.
#' @param xlab,ylab,pch,ylim,xlim,main,font.main,col,cex graphical parameters
#' @param length The width of the ends of the error bars.
#' @param legend Include a legend or not. Default is TRUE.
#' @param legend.title Title of the legend. No title by default.
#' @param auc Plot the AUC or not (if available). Default is TRUE.
#' @param brier Plot the Brier Score or not (if available). Default is TRUE.
#' @param ... Additional arguments to `plot()`
#'
#' @export
plot.LMScore <- function(x,
                         metrics,
                         contrasts = FALSE,
                         landmarks = TRUE,
                         summary = TRUE,
                         se = TRUE,
                         add_pairwise_contrasts = FALSE,
                         cutoff_contrasts = 0.05,
                         pairwise_heights,
                         width,
                         loc, xlab, ylab, pch, ylim, xlim, main, font.main = 1,
                         col = NULL, cex = 1,
                         length = 0.1,
                         legend = TRUE,
                         legend.title = NULL,
                         auc = TRUE,
                         brier = TRUE, ...) {
  if (missing(metrics)) {
    metrics <- c()
    if (!is.null(x$AUC) && auc) metrics <- c("AUC")
    if (!is.null(x$Brier) && brier) metrics <- c(metrics, "Brier")
  }

  if (landmarks) {
    if (missing(xlab)) xlab <- "Landmark Time (t)"
    if (missing(pch)) pch <- 19

    if (missing(loc)) set_x <- TRUE
    else set_x <- FALSE

    for (metric in metrics) {
      if (is.null(x[[metric]])) {
        warning(tidymess(paste(
          metric, "was not set as a metric when calling score(). No results to
          plot. Either call score() again with", metric, "as a metric or do not
          include it as a metric here.")))
      } else {
        if (set_x && (metric == "AUC")) loc <- "topright"
        else if (set_x && (metric == "Brier"))  loc <- "bottomright"
        if (contrasts) {
          df <- x[[metric]]$contrasts
          df$model <- factor(paste(df$model, "-", df$reference))
          label <- paste0("delta.", metric)
        } else {
          df <- x[[metric]]$score
          label <- metric
        }
        plot_metric(df, label, col, loc, ylim, xlim, main, ylab, xlab, x$w,
                    font.main, pch, length, legend, legend.title, cex,
                    contrasts = contrasts, ...)
      }
    }

  }

  if (summary) { ### plot summary metric
    if (!requireNamespace("latex2exp", quietly = TRUE)) {
      stop(tidymess("Package \"latex2exp\" must be installed to use function
                      plot.LMScore on summary metrics"), call. = FALSE)
    }

    for (metric in metrics) {
      metric_summary <- paste0(metric, "_summary")

      if (!is.null(x[[metric_summary]])) {
        if (!contrasts) {
          df <- x[[metric_summary]]$score
          col_name <- metric
        } else {
          df <- x[[metric_summary]]$contrasts
          df$model <- factor(paste(df$model, "-", df$reference))
          col_name <- paste0("delta.", metric)
        }

        mid <- plot_summary(df, col_name, col, ylim, xlab, ylab, main, cex,
                            font.main, length, se, contrasts = contrasts, ...)

        # add significant contrasts
        if (add_pairwise_contrasts && contrasts) {
          warning(tidymess(
            "Plotting pairwise contrasts (add_pairwise_contrasts=TRUE) does not
            make sense when already plotting contrasts (contrasts=TRUE). No
            pairwise contrasts are plotted."))
        }
        else if (add_pairwise_contrasts) {
          if (missing(pairwise_heights) || missing(width)) {
            stop(tidymess("When plotting pairwise contrasts, arguments
                          'pairwise_heights' and 'width' must be given."))
          }
          if (missing(ylim)) {
            warning(tidymess("Using the default ylim may results in pairwise
                             contrasts not being shown, try expanding it."))
          }
          x_axis <- as.numeric(mid)
          names(x_axis) <- df$model

          df <- x[[metric_summary]]$contrasts
          df$pair <- factor(paste(df$model, "-", df$reference))
          col_name <- paste0("delta.", metric)

          num_done <- 0
          for (i in seq_len(nrow(df))) {
            row <- df[i, ]
            p_val <- row$p
            if (p_val < cutoff_contrasts) {
              num_done <- num_done + 1
              x0 <- x_axis[as.character(row$reference)]
              x1 <- x_axis[as.character(row$model)]
              graphics::segments(x0 = x0, x1 = x0,
                                 y0 = pairwise_heights[num_done] - width,
                                 y1 = pairwise_heights[num_done])
              graphics::segments(x0 = x1, x1 = x1,
                                 y0 = pairwise_heights[num_done] - width,
                                 y1 = pairwise_heights[num_done])
              graphics::segments(x0 = x0, x1 = x1,
                                 y0 = pairwise_heights[num_done],
                                 y1 = pairwise_heights[num_done])
              graphics::text(x = (x0 + x1) / 2,
                             y = pairwise_heights[num_done] + width,
                             cex = cex,
                             labels = paste(format.pval(p_val, digits = 2,
                                                        eps = 1e-5)))
            }
          }
        }
      }
    }
  }
}


#' Plot an object output from [calplot()]: plot the calibration plots.
#'
#' @param x An object of class "LMcalibrationPlot" output from [calplot()]
#' @param main Optional title to override default.
#' @param ... Other arguments to pass to pass to plot
#'
#' @export
#'
plot.LMcalibrationPlot <- function(x, main, ...) {
  x$regression_values <- NULL
  add_title <- FALSE
  if (missing(main)) add_title <- TRUE
  for (i in 1:length(x)) {
    plot(x[[i]], ...)

    if (add_title) {
      title <- paste0("Risk calibration")
      graphics::title(main = title)
    } else {
      graphics::title(main = main)
    }
  }
}


#' Plot the coefficient path created by calling [pen_lm()], analogous to
#' plotting from [glmnet()]
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
#' @param x a fitted [pen_lm()] object
#' @param xvar As in [glmnet()]: "What is on the X-axis. "`norm`" plots against
#'   the L1-norm of the coefficients, "`lambda`" against the log-lambda
#'   sequence, and "`dev`" against the percent deviance explained."
#' @param all_causes if pen_lm fit a cause-specific Cox model, set TRUE to plot
#'   coefficient profile plots for each model.
#' @param silent Set TRUE to hide messages.
#' @param label Set TRUE to label the curves by variable index numbers.
#' @param all_causes_title If `all_causes` is set to TRUE, includes a title with
#'   the cause. Defaults to TRUE.
#' @param \dots additional graphical parameters
#' @export
plot.pen_lm <- function(x, xvar = "norm", all_causes = FALSE, silent = FALSE,
                       label = FALSE, all_causes_title = TRUE, ...) {
  if (all_causes) {
    for (i in seq_along(x)){
      plot(x[[i]], xvar, label, ...)
      if (all_causes_title) graphics::title(paste("Cause", i), line = 2.5)
    }
  } else {
    plot(x[[1]], xvar, label, ...)
    if (length(x) > 1 && !silent)
      message(tidymess("\n (To print plot paths for remaining cause-specific
          models, call plot with argument all_causes = TRUE)\n"))
  }
}


#' Plot cross-validation curve created by [cv.pen_lm()], analogous to plotting
#' from [cv.glmnet()]
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
#' @param x a fitted [cv.pen_lm()] object
#' @param all_causes if `pen_lm()` fit a cause-specific Cox model, set TRUE to
#'   plot coefficient profile plots for each model.
#' @param silent Set TRUE to hide messages.
#' @param label Set TRUE to label the curves by variable index numbers.
#' @param sign.lambda Plot against `log(lambda)` (default) or its negative
#'   if set to -1.
#' @param se.bands Logical. If TRUE, shading is produced to show stand-error
#'   bands. Defaults to TRUE.
#' @param all_causes_title If `all_causes` is set to TRUE, includes a title with
#'   the cause. Defaults to TRUE.
#' @param \dots additional graphical parameters
#' @export
plot.cv.pen_lm <- function(x, all_causes = FALSE, silent = FALSE, label = FALSE,
                          sign.lambda = 1, se.bands = TRUE,
                          all_causes_title = TRUE, ...) {
  if (all_causes) {
    for (i in seq_along(x)){
      plot(x[[i]], sign.lambda, se.bands, ...)
      if (all_causes_title) graphics::title(paste("Cause", i), line = 2.5)
    }
  } else {
    plot(x[[1]], sign.lambda, se.bands, ...)
    if (length(x) > 1 && !silent)
      message(tidymess("\n (To print plot paths for remaining cause-specific
          models, call plot with argument all_causes = TRUE)\n"))
  }
}



#' Generic function to plot coefficients
#'
#' Can plot positive and negative coefficients in two separate plots or the
#' same. X-axes are the same if separate plots are used.
#'
#' @param x (Named) Vector of coefficients
#' @param single_plot Logical, defaults to TRUE. A single plot for both
#'   positive and negative coefficients, or two separate plots.
#' @param max_coefs Default is to plot all coefficients. If specified, gives the
#'   maximum number of coefficients to plot.
#' @param col Fill color for the barplot.
#' @param xlab x-axis Label
#' @param ... Additional arguments to barplot.
#' @export
plot.coefs <- function(x, single_plot = TRUE, max_coefs = NULL,
                       col = "blue", xlab = "Coefficient value", ...) {
  coefs <- x
  pos_coefs <- sort(coefs[coefs > 0], decreasing = TRUE)
  neg_coefs <- (-sort(-(coefs[coefs < 0]), decreasing = TRUE))

  if (!single_plot) {
    xmax <- 1.1 * max(c(pos_coefs, -neg_coefs), na.rm = TRUE)
    ymax <- max(
      min(length(pos_coefs), max_coefs),
      min(length(neg_coefs), max_coefs)
    )
    if (length(neg_coefs) > 0) {
      neg_coefs <- neg_coefs[1:min(length(neg_coefs), max_coefs)]
      graphics::barplot(neg_coefs,
                        col = col,
                        names.arg = names(neg_coefs),
                        horiz = TRUE, las = 1, xpd = FALSE,
                        xlim = c(-xmax, 0),
                        ylim = c(0, ymax),
                        width = 0.8,
                        xlab = xlab)
    }
    if (length(pos_coefs) > 0) {
      pos_coefs <- pos_coefs[1:min(length(pos_coefs), max_coefs)]
      graphics::barplot(pos_coefs,
                        col = col,
                        names.arg = names(pos_coefs),
                        horiz = TRUE, las = 1, xpd = FALSE,
                        xlim = c(0, xmax),
                        ylim = c(0, ymax),
                        width = 0.8,
                        xlab = xlab)
    }

  } else {
    all_coefs <- c(pos_coefs, neg_coefs)
    all_coefs <- all_coefs[order(-abs(all_coefs))]
    all_coefs <- all_coefs[1:min(length(all_coefs), max_coefs)]
    all_coefs <- rev(all_coefs)

    graphics::barplot(all_coefs, col = col, names.arg = names(all_coefs),
                      horiz = TRUE, las = 1, xpd = FALSE, xlab = xlab, ...)
  }
}

