#' Print function for object of class LMpred
#'
#' @param x Object of class LMpred
#' @param verbose Boolean, default is FALSE. Print further components.
#' @param ... Arguments passed to print.
#'
#' @return Printed output.
#' @export
#'
print.LMpred <- function(x, verbose = FALSE, ...) {
  cat("$preds\n")
  print(utils::head(x$preds))
  if (nrow(x$preds) > 5) cat(paste(" [ omitted", nrow(x$preds) - 5, "rows ]\n"))

  if (verbose) {
    cat("\n")

    cat("$w\n")
    print(x$w)
    cat("\n")

    cat("$type\n")
    print(x$type)
    cat("\n")

    cat("$LHS\n")
    print(x$LHS)
    cat("\n")

    cat("$data\n")
    print(utils::head(x$data))
    if (nrow(x$data) > 5) cat(paste(" [ omitted", nrow(x$data) - 5, "rows ]\n"))
    cat("\n")

    cat("$cause\n")
    print(x$cause)
    cat("\n")
  }

}

#' Print function for object of class LMdataframe
#'
#' @param x Object of class LMdataframe
#' @param verbose Boolean, default is FALSE. Print further components.
#' @param ... Arguments passed to print.
#'
#' @return Printed output.
#' @export
#'
print.LMdataframe <- function(x, verbose = FALSE, ...) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package \"data.table\" must be installed to use function score.",
         call. = FALSE)
  }

  cat("$data\n")
  print(utils::head(x$data))
  if (nrow(x$data) > 5) cat(paste(" [ omitted", nrow(x$data) - 5, "rows ]\n"))

  if (verbose) {
    cat("\n")
    cat("\n")
    cat("$outcome\n")
    names.outcome <- names(x$outcome)
    for (i in seq_along(x$outcome)) {
      if (is.null(names.outcome[i])) label <- paste0("[[", i, "]]")
      else label <- paste0("$", names.outcome[i])
      cat(paste0("$outcome", label, "\n"))
      print(x$outcome[[i]])
      cat("\n")
    }

    cat("$w\n")
    print(x$w)
    cat("\n")

    cat("$end_time\n")
    print(x$end_time)
    cat("\n")

    names.LMdata <- names(x)

    if ("func_covars" %in% names.LMdata) {
      cat("$func_covars\n")
      names.fc <- names(x$func_covars)
      for (i in seq_along(x$func_covars)) {
        if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
        else paste0("$", names.fc[i])
        cat(paste0("$func_covars$", label, "\n"))
        print(x$func_covars[[i]])
        cat("\n")
      }
    }
    if ("func_lms" %in% names.LMdata) {
      cat("$func_lms\n")
      names.fc <- names(x$func_lms)
      for (i in seq_along(x$func_lms)) {
        if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
        else paste0("$", names.fc[i])
        cat(paste0("$func_lms$", label, "\n"))
        print(x$func_lms[[i]])
        cat("\n")
      }
    }
    if ("lm_covs" %in% names.LMdata) {
      cat("$lm_covs\n")
      print(x$lm_covs)
    }
    if ("all_covs" %in% names.LMdata) {
      cat("\n$all_covs\n")
      print(x$all_covs)
    }
    if ("lm_col" %in% names.LMdata) {
      cat("\n$lm_col\n")
      print(x$lm_col)
    }
    if ("id_col" %in% names.LMdata) {
      cat("\n$id_col\n")
      print(x$id_col)
    }
  }
}

#' Print function for object of class LMScore, i.e., output from [score()]
#'
#' @param x Object of class LMScore
#' @param digits Number of significant digits to include
#' @param landmarks Print the time-dependent metrics at individual landmarks.
#'   Default is TRUE.
#' @param summary Print the summary metrics of the models if they have been
#'   calculated. Default is TRUE.
#' @param ... Arguments passed to print.
#'
#' @importFrom data.table :=
#'
#' @return Printed output.
#' @export
print.LMScore <- function(x, digits = 3, landmarks = TRUE, summary = TRUE, ...) {
  list_indexes <- c("AUC", "Brier", "AUC_summary", "Brier_summary")
  out_class <- c("scoreAUC", "scoreBrier", "scoreAUC", "scoreBrier")
  names <- c("AUC", "Brier Score", "AUC", "Brier Score")

  metrics_to_print <- c()
  if (landmarks) metrics_to_print <- 1:2
  if (summary) metrics_to_print <- c(metrics_to_print, 3:4)

  for (i in metrics_to_print) {
    if (!is.null(x[[list_indexes[i]]])) {
      if (i <=2 )
        cat(paste0("\nMetric: Time-dependent ", names[i], " (w = ", x$w, ")\n"))
      else
        cat(paste0("\nMetric: Averaged time-dependent ", names[i], " (w = ", x$w, ")\n"))

      if (!is.null(x[["B"]]))
        cat("        Bootstrapped over", x[["B"]], "iterations\n")

      obj <- x[[list_indexes[i]]]
      if (!is.null(obj$score$times)) obj$score$times <- NULL
      if (!is.null(obj$contrasts$times)) obj$contrasts$times <- NULL
      class(obj) <- out_class[i]
      print(obj)
      message(paste("NOTE: Predictions are made at time tLM for time tLM +", x$w))
    }
  }
}


#' Print function for object of class LMCSC
#'
#' @param x Object of class LMCSC
#' @param verbose Boolean, default is FALSE. Print further components.
#' @param cause Print the model for a given cause.
#'   If left out, all models are printed.
#' @param ... Arguments passed to print.
#'
#' @return Printed output.
#' @export
#'
print.LMCSC <- function(x, verbose = FALSE, cause, ...) {
  cat(paste0(
    "\nLandmark cause-specific cox super model fit for dynamic prediction of window size ",
    x$w, ":\n\n"))

  cat("$model\n")
  if (missing(cause)) {
    num_causes <- length(x$model$causes)
    for (i in 1:num_causes){
      cat(paste0("----------> Cause: ", i, "\n"))
      cox_model <- x$model$models[[i]]
      cox_model$call <- NULL
      cox_model$loglik <- NA
      print(cox_model)
      cat("\n\n")
    }
  } else {
    cat(paste0("----------> Cause: ", cause, "\n"))
    cox_model <- x$model$models[[cause]]
    cox_model$call <- NULL
    cox_model$loglik <- NA
    print(cox_model)
    cat("\n\n")
  }

  if (verbose) {
    cat("$func_covars\n")
    names.fc <- names(x$func_covars)
    for (i in seq_along(x$func_covars)) {
      if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
      else paste0("$", names.fc[i])
      cat(paste0("$func_covars$", label, "\n"))
      print(x$func_covars[[i]])
      cat("\n")
    }
    cat("$func_lms\n")
    names.fc <- names(x$func_lms)
    for (i in seq_along(x$func_lms)) {
      if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
      else paste0("$", names.fc[i])
      cat(paste0("$func_lms$", label, "\n"))
      print(x$func_lms[[i]])
      cat("\n")
    }

    cat("$w\n")
    print(x$w)
    cat("\n")

    cat("$end_time\n")
    print(x$end_time)
    cat("\n")

    cat("$type\n")
    print(x$type)
    cat("\n")
  }
}


#' Print function for object of class LMcoxph
#'
#' @param x Object of class LMcoxph
#' @param verbose Boolean, default is FALSE. Print further components.
#' @param ... Arguments passed to print.
#'
#' @return Printed output.
#' @export
#'
print.LMcoxph <- function(x, verbose = FALSE, ...) {
  cat(paste0(
    "\nLandmark cox super model fit for dynamic prediction of window size ",
    x$w, ":\n\n"))
  if (verbose) cat("$model\n")

  cox_model <- x$model
  cox_model$call <- NULL
  cox_model$loglik <- NA
  print(cox_model)

  if (verbose) {
    cat("\n\n")
    cat("$func_covars\n")
    names.fc <- names(x$func_covars)
    for (i in seq_along(x$func_covars)){
      if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
      else paste0("$", names.fc[i])
      cat(paste0("$func_covars$", label, "\n"))
      print(x$func_covars[[i]])
      cat("\n")
    }
    cat("$func_lms\n")
    names.fc <- names(x$func_lms)
    for (i in seq_along(x$func_lms)){
      if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
      else paste0("$", names.fc[i])
      cat(paste0("$func_lms$", label, "\n"))
      print(x$func_lms[[i]])
      cat("\n")
    }

    cat("$w\n")
    print(x$w)
    cat("\n")

    cat("$end_time\n")
    print(x$end_time)
    cat("\n")

    cat("$type\n")
    print(x$type)
    cat("\n")
  }
}


#' Print the output from calling `pen_lm()`
#'
#' Similar to the output of printing a`glmnet` object, print a summary of the
#' glmnet path at each step along the path.
#'
#' @details
#' If the model is a survival model (i.e., no competing risks), then the output
#' is the same as a call to glmnet would produce. For competing risks, the
#' default is only to print the output for the cause of interest (first cause).
#' Further events can be examined by setting `all_causes = TRUE`.
#'
#' As in glmnet,
#' "A three-column matrix with columns `Df`, `%Dev` and `Lambda` is printed.
#' The `Df` column is the number of nonzero coefficients (Df is a
#' reasonable name only for lasso fits). `%Dev` is the percent deviance
#' explained (relative to the null deviance)."
#'
#' @param x a pen_lm object
#' @param all_causes if `pen_lm()` fit a cause-specific Cox model, set TRUE to
#'   print a summary of the glmnet path for each model.
#' @param silent Set TRUE to hide messages.
#' @param digits Number of significant digits to include
#' @param \dots additional print arguments
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008).
#'  Regularization Paths for Generalized Linear Models via Coordinate Descent
#' @export
print.pen_lm <- function(x, all_causes = FALSE, silent = FALSE,
                         digits = 3, ...) {
    num_causes <- length(x)
    if (attr(x, "survival.type") == "competing.risk") {
      cat("\n", tidymess("Penalized cause-specific Cox landmark supermodel fit
                         for dynamic prediction:"))
    } else {
      cat("\nPenalized Cox landmark supermodel fit for dynamic prediction:")
    }

    if (all_causes) {
      for (i in 1:num_causes){
        cat(paste0("\n\nCause ", i, ":\n"))
        print(x[[i]])
      }

    } else {
      cat("\n\n")
      if (num_causes > 1) cat(paste0("First Cause:\n"))
      print(x[[1]])
      if (num_causes > 1 && !silent)
        message("\n", tidymess("(To print paths for the remaining
                      cause-specific models, call print with argument
                      all_causes = TRUE)"), "\n")
    }
}


#' Print the output from calling cv.pen_lm()`,

#' Similar to printing the output of printing a `cv.glmnet` object, print a
#' cross-validated penalized cause-specific Cox supermodel
#'
#' @details
#' If the model is a survival model (i.e., no competing risks), then the output
#' is the same as a call to glmnet would produce. For competing risks, the
#' default is only to print the output for the cause of interest (first cause).
#' Further events can be examined by setting `all_causes = TRUE`.
#'
#' @param x a cv.pen_lm object
#' @param all_causes if `cv.pen_lm` fit a cause-specific Cox model, set TRUE to
#'   print a summary of the glmnet path for each model.
#' @param silent Set TRUE to hide messages.
#' @param digits Number of significant digits to include
#' @param \dots additional print arguments
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008).
#'   Regularization Paths for Generalized Linear Models via Coordinate Descent
#' @export
print.cv.pen_lm <- function(x, all_causes = FALSE, silent = FALSE,
                            digits = 3, ...) {
  num_causes <- length(x)

  if (attr(x, "survival.type") == "competing.risk") {
    cat("\n", tidymess("Cross-validated penalized cause-specific Cox landmark
                       supermodel fit for dynamic prediction:"))
  } else {
    cat("\n", tidymess("Cross-validated penalized Cox landmark supermodel
                       fit for dynamic prediction:"))
  }
  if (all_causes){
    for (i in 1:num_causes){
      cat(paste0("\n\nCause ", i, ":\n"))
      print(x[[i]])
    }

  } else {
    cat("\n\n")
    if (num_causes > 1) cat(paste0("First Cause:\n"))
    print(x[[1]])
    if (num_causes > 1 & !silent)
      message("\n", tidymess("(To print for remaining cause-specific models,
                          call print with argument all_causes = TRUE)"), "\n")
  }
}


#' Print function for object of class penLMCSC
#'
#' @param x Object of class penLMCSC
#' @param cause Print the model for a given cause.
#'   If left out, all models are printed.
#' @param verbose Logical, if verbose print func_covars, func_lms, w, end_time
#'   and type.
#' @param ... Arguments passed to print.
#'
#' @return Printed output.
#' @export
#'
print.penLMCSC <- function(x, cause, verbose = FALSE, ...) {
  cat(paste0("\nPenalized landmark cause-specific Cox super model fit for dynamic prediction of window size ",
             x$w, "\n"))
  cat("(Note that zero-valued coefficients are not printed)\n\n")

  if (verbose) cat("$model\n")
  if (missing(cause)) {
    num_causes <- length(x$model$causes)
    for (i in 1:num_causes){
      cat(paste0("----------> Cause: ", i, "\n"))
      cat(paste0("            (lambda = ", x$lambda[[i]], ")\n"))
      coefs = x$model$models[[i]]$coefficients
      non_zero_coefs = coefs[coefs != 0]
      non_zero_coefs = cbind(non_zero_coefs, exp(non_zero_coefs))
      colnames(non_zero_coefs) = c("coef", "exp(coef)")
      print(non_zero_coefs)
      cat("\n\n")
    }
  } else {
    cat(paste0("----------> Cause: ", cause, "\n"))
    cat(paste0("            (lambda = ", x$lambda[[cause]], ")\n"))
    coefs = x$model$models[[cause]]$coefficients
    non_zero_coefs = coefs[coefs != 0]
    non_zero_coefs = cbind(non_zero_coefs, exp(non_zero_coefs))
    colnames(non_zero_coefs) = c("coef", "exp(coef)")
    print(non_zero_coefs)
    cat("\n\n")
  }

  if (verbose) {
    cat("$func_covars\n")
    names.fc = names(x$func_covars)
    for (i in seq_along(x$func_covars)){
      if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
      else paste0("$", names.fc[i])
      cat(paste0("$func_covars$", label, "\n"))
      print(x$func_covars[[i]])
      cat("\n")
    }
    cat("$func_lms\n")
    names.fc = names(x$func_lms)
    for (i in seq_along(x$func_lms)) {
      if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
      else paste0("$", names.fc[i])
      cat(paste0("$func_lms$", label, "\n"))
      print(x$func_lms[[i]])
      cat("\n")
    }

    cat("$w\n")
    print(x$w)
    cat("\n")

    cat("$end_time\n")
    print(x$end_time)
    cat("\n")

    cat("$type\n")
    print(x$type)
    cat("\n")
  }
}


#' Print function for object of class penLMcoxph
#'
#' @param x Object of class penLMcoxph
#' @param verbose Logical, if verbose print func_covars, func_lms, w, end_time
#'   and type.
#' @param ... Arguments passed to print.
#'
#' @return Printed output.
#' @export
#'
print.penLMcoxph <- function(x, verbose = FALSE, ...) {
  cat(paste0("\nPenalized landmark Cox super model fit for dynamic prediction of window size ",
             x$w, " (lambda = ", x$lambda, ")", ":\n\n"))
  if (verbose) cat("$model\n")

  coefs = x$model$coefficients
  non_zero_coefs = coefs[coefs != 0]
  non_zero_coefs = cbind(non_zero_coefs, exp(non_zero_coefs))
  colnames(non_zero_coefs) = c("coef", "exp(coef)")
  print(non_zero_coefs)
  cat("\n\n")

  if (verbose) {
    cat("$func_covars\n")
    names.fc = names(x$func_covars)
    for (i in seq_along(x$func_covars)){
      if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
      else paste0("$", names.fc[i])
      cat(paste0("$func_covars$", label, "\n"))
      print(x$func_covars[[i]])
      cat("\n")
    }
    cat("$func_lms\n")
    names.fc = names(x$func_lms)
    for (i in seq_along(x$func_lms)) {
      if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
      else paste0("$", names.fc[i])
      cat(paste0("$func_lms$", label, "\n"))
      print(x$func_lms[[i]])
      cat("\n")
    }

    cat("$w\n")
    print(x$w)
    cat("\n")

    cat("$end_time\n")
    print(x$end_time)
    cat("\n")

    cat("$type\n")
    print(x$type)
    cat("\n")
  }
}
