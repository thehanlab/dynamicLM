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
    names.outcome = names(x$outcome)
    for (i in 1:length(x$outcome)) {
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
      for (i in 1:length(x$func_covars)) {
        if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
        else paste0("$", names.fc[i])
        cat(paste0("$func_covars$", label, "\n"))
        print(x$func_covars[[i]])
        cat("\n")
      }
    }
    if ("func_lms" %in% names.LMdata) {
      cat("$func_lms\n")
      names.fc = names(x$func_lms)
      for (i in 1:length(x$func_lms)) {
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
      cat("\n")
    }
    if ("all_covs" %in% names.LMdata) {
      cat("$all_covs\n")
      print(x$all_covs)
      cat("\n")
    }
    if ("lm_col" %in% names.LMdata) {
      cat("$lm_col\n")
      print(x$lm_col)
    }
  }
}

#' Print function for object of class LMScore, i.e., output from [score()]
#'
#' @param x Object of class LMScore
#' @param digits Number of significant digits to include
#' @param ... Arguments passed to print.
#'
#' @importFrom data.table :=
#'
#' @return Printed output.
#' @export
#'
print.LMScore <- function(x, digits=3, ...) {

  if (nrow(x$auct)>0) {
    cat(paste0("\nMetric: Time-dependent AUC (window ", x$w, ")\n"))
    cat("\nResults by model:\n")

    AUC = se = times = lower = upper = NULL
    fmt <- paste0("%1.", digits[[1]], "f")
    X <- data.table::copy(x$auct)
    X[, AUC := sprintf(fmt = fmt, 100 * AUC)]
    if (match("se", colnames(X), nomatch = 0)) X[, se := NULL]
    if (match("times", colnames(X), nomatch = 0)) X[, times := NULL]
    if (match("lower", colnames(X), nomatch = 0))
      X[, lower := sprintf(fmt = fmt, 100 * lower)]
    if (match("upper", colnames(X), nomatch = 0))
      X[, upper := sprintf(fmt = fmt, 100 * upper)]

    print(X, digits = digits)

    # cat("\nResults by comparison: \nTODO\n")

    message("NOTE: Values are multiplied by 100 and given in %.")
    message("NOTE: The higher AUC the better.")
    message(paste("NOTE: Predictions are made at time tLM for risk windows of length", x$w))
  }
  if (nrow(x$briert)>0) {
    cat(paste0("\nMetric: Time-dependent Brier Score (window ", x$w, ")\n"))
    cat("\nResults by model:\n")

    Brier = se = times = se.conservative = lower = upper = NULL
    fmt <- paste0("%1.", digits[[1]], "f")
    X <- data.table::copy(x$briert)
    X[, Brier:=sprintf(fmt = fmt, 100 * Brier)]
    if (match("se", colnames(X), nomatch = 0)) X[, se := NULL]
    if (match("times", colnames(X), nomatch = 0)) X[, times := NULL]
    if (match("se.conservative", colnames(X), nomatch = 0))
      X[, se.conservative := NULL]
    if (match("lower", colnames(X), nomatch = 0))
      X[, lower := sprintf(fmt = fmt, 100 * lower)]
    if (match("upper", colnames(X), nomatch = 0))
      X[, upper := sprintf(fmt = fmt, 100 * upper)]
    print(X)

    # cat("\nResults by comparison: \nTODO\n")

    message("NOTE: Values are multiplied by 100 and given in %.")
    message("NOTE: The lower Brier the better.")
    message(paste("NOTE: Predictions are made at time tLM for risk windows of length", x$w))
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
    "\nLandmark cause-specific cox super model fit for dynamic prediction of window size "
    x$w, ":\n\n"))

  cat("$model\n")
  if (missing(cause)) {
    num_causes <- length(x$model$causes)
    for (i in 1:num_causes){
      cat(paste0("----------> Cause: ", i, "\n"))
      cox_model <- x$model$models[[i]]
      cox_model$call <- NULL
      print(cox_model)
      cat("\n\n")
    }
  } else {
    cat(paste0("----------> Cause: ", cause, "\n"))
    cox_model <- x$model$models[[cause]]
    cox_model$call <- NULL
    print(cox_model)
    cat("\n\n")
  }
  if(verbose){
    cat("$func_covars\n")
    names.fc = names(x$func_covars)
    for (i in 1:length(x$func_covars)){
      if (is.null(names.fc[i])) label <- paste0("[[",i,"]]")
      else paste0("$",names.fc[i])
      cat(paste0("$func_covars$",label,"\n"))
      print(x$func_covars[[i]])
      cat("\n")
    }
    cat("$func_LMs\n")
    names.fc = names(x$func_LMs)
    for (i in 1:length(x$func_LMs)){
      if (is.null(names.fc[i])) label <- paste0("[[",i,"]]")
      else paste0("$",names.fc[i])
      cat(paste0("$func_LMs$",label,"\n"))
      print(x$func_LMs[[i]])
      cat("\n")
    }

  if (verbose) {
    cat("$func_covars\n")
    names.fc = names(x$func_covars)
    for (i in 1:length(x$func_covars)){
      if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
      else paste0("$", names.fc[i])
      cat(paste0("$func_covars$", label, "\n"))
      print(x$func_covars[[i]])
      cat("\n")
    }
    cat("$func_lms\n")
    names.fc = names(x$func_lms)
    for (i in 1:length(x$func_lms)){
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
  cat("$model\n")

  cox_model <- x$model
  cox_model$call <- NULL
  print(cox_model)

  if (verbose) {
    cat("\n\n")

    cat("$func_covars\n")
    names.fc = names(x$func_covars)
    for (i in 1:length(x$func_covars)){
      if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
      else paste0("$", names.fc[i])
      cat(paste0("$func_covars$", label, "\n"))
      print(x$func_covars[[i]])
      cat("\n")
    }
    cat("$func_lms\n")
    names.fc = names(x$func_lms)
    for (i in 1:length(x$func_lms)){
      if (is.null(names.fc[i])) label <- paste0("[[", i, "]]")
      else paste0("$",names.fc[i])
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


#' Print the output from calling `LMpen`, similar to the output of printing a `glmnet` object
#'
#' I.e., print a summary of the glmnet path at each step along the path.
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
#' @param x a penLM object
#' @param all_causes if penLM fit a cause-specific Cox model, set TRUE to print
#'   a summary of the glmnet path for each model.
#' @param silent Set TRUE to hide messages.
#' @param digits Number of significant digits to include
#' @param \dots additional print arguments
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008). Regularization Paths for Generalized Linear Models via Coordinate Descent
#' @export
print.penLM <- function (x, all_causes = FALSE, silent = FALSE, digits = 3, ...) {
    num_causes <- length(x)

    cat(paste0("\nPenalized landmark Cox super model fit for dynamic prediction:"))
    if (all_causes){
      for (i in 1:num_causes){
        cat(paste0("\n\nCause ",i,":\n"))
        print(x[[i]])
      }

    } else {
      cat("\n\n")
      if (num_causes > 1) cat(paste0("First Cause:\n"))
      print(x[[1]])
      if (num_causes > 1 & !silent) message("\n (To print paths for remaining cause-specific models, call print with argument all_causes = TRUE)\n")
    }
}

#' Print the output from calling `cv.LMpen`, similar to the output of printing a `cv.glmnet` object
#'
#' I.e., print a cross-validated LMpen object
#' @details
#' If the model is a survival model (i.e., no competing risks), then the output
#' is the same as a call to glmnet would produce. For competing risks, the
#' default is only to print the output for the cause of interest (first cause).
#' Further events can be examined by setting `all_causes = TRUE`.
#'
#' @param x a penLM object
#' @param all_causes if penLM fit a cause-specific Cox model, set TRUE to print
#'   a summary of the glmnet path for each model.
#' @param silent Set TRUE to hide messages.
#' @param digits Number of significant digits to include
#' @param \dots additional print arguments
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008). Regularization Paths for Generalized Linear Models via Coordinate Descent
#' @export
print.cv.penLM <- function (x, all_causes = FALSE, silent = FALSE, digits = 3, ...) {
  num_causes <- length(x)

  cat(paste0("\nCross-validated penalized landmark Cox super model fit for dynamic prediction:"))
  if (all_causes){
    for (i in 1:num_causes){
      cat(paste0("\n\nCause ",i,":\n"))
      print(x[[i]])
    }

  } else {
    cat("\n\n")
    if (num_causes > 1) cat(paste0("First Cause:\n"))
    print(x[[1]])
    if (num_causes > 1 & !silent) message("\n (To print for remaining cause-specific models, call print with argument all_causes = TRUE)\n")
  }
}


#' Print function for object of class penLMCSC
#'
#' @param x Object of class penLMCSC
#' @param verbose Logical, if verbose print func_covars, func_LMs, w, end_time
#'   and type.
#' @param ... Arguments passed to print.
#'
#' @return Printed output.
#' @export
#'
print.penLMCSC <- function(x, verbose=FALSE, ...) {
  cat(paste0("\nPenalized landmark cause-specific Cox super model fit for dynamic prediction of window size ",x$w,"\n"))
  cat("(Note that zero-valued coefficients are not printed)\n\n")

  if(verbose) cat("$model\n")
  num_causes <- length(x$model$causes)
  for (i in 1:num_causes){
    cat(paste0("----------> Cause: ",i,"\n"))
    cat(paste0("            (s = ",x$s[[i]],")\n"))
    coefs = x$model$models[[i]]$coefficients
    non_zero_coefs = coefs[coefs!=0]
    non_zero_coefs = cbind(non_zero_coefs, exp(non_zero_coefs))
    colnames(non_zero_coefs) = c("coef", "exp(coef)")
    print(non_zero_coefs)
    cat("\n\n")
  }

  if (verbose){
    cat("$func_covars\n")
    names.fc = names(x$func_covars)
    for (i in 1:length(x$func_covars)){
      if (is.null(names.fc[i])) label <- paste0("[[",i,"]]")
      else paste0("$",names.fc[i])
      cat(paste0("$func_covars$",label,"\n"))
      print(x$func_covars[[i]])
      cat("\n")
    }
    cat("$func_LMs\n")
    names.fc = names(x$func_LMs)
    for (i in 1:length(x$func_LMs)){
      if (is.null(names.fc[i])) label <- paste0("[[",i,"]]")
      else paste0("$",names.fc[i])
      cat(paste0("$func_LMs$",label,"\n"))
      print(x$func_LMs[[i]])
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
#' @param verbose Logical, if verbose print func_covars, func_LMs, w, end_time
#'   and type.
#' @param ... Arguments passed to print.
#'
#' @return Printed output.
#' @export
#'
print.penLMcoxph <- function(x, verbose=FALSE, ...) {
  cat(paste0("\nPenalized landmark Cox super model fit for dynamic prediction of window size ",x$w," (s = ",x$s,")",":\n\n"))
  if (verbose)cat("$model\n")

  coefs = x$model$coefficients
  non_zero_coefs = coefs[coefs!=0]
  non_zero_coefs = cbind(non_zero_coefs, exp(non_zero_coefs))
  colnames(non_zero_coefs) = c("coef", "exp(coef)")
  print(non_zero_coefs)
  cat("\n\n")

  if(verbose){
    cat("$func_covars\n")
    names.fc = names(x$func_covars)
    for (i in 1:length(x$func_covars)){
      if (is.null(names.fc[i])) label <- paste0("[[",i,"]]")
      else paste0("$",names.fc[i])
      cat(paste0("$func_covars$",label,"\n"))
      print(x$func_covars[[i]])
      cat("\n")
    }
    cat("$func_LMs\n")
    names.fc = names(x$func_LMs)
    for (i in 1:length(x$func_LMs)){
      if (is.null(names.fc[i])) label <- paste0("[[",i,"]]")
      else paste0("$",names.fc[i])
      cat(paste0("$func_LMs$",label,"\n"))
      print(x$func_LMs[[i]])
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
