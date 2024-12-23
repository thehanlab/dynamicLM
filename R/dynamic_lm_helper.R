dynamic_lm_helper <- function(formula, type, data, lmdata, method, cluster,
                              x, w, end_time, func_covars, func_lms, lm_covs,
                              all_covs, outcome, lm_col, original.landmarks,
                              args, ...) {
  # extra LHS of formula
  LHS <- getLHS(formula)

  # check if type and formula match
  if ((type == "coxph" && formula[[2]][[1]] != "Surv") ||
      (type == "CSC" && formula[[2]][[1]] != "Hist")) {
      stop(tidymess("Mismatch between formula's left-hand side (LHS) and type
          of supermodel. For standard survival data (one event and possible
          censoring), use type = \"coxph\" and a LHS of the form Surv(LM, Time,
          event). For competing risks (multiple events and possible censoring),
          use type = \"CSC\" and a LHS of the form Hist(Time, event, LM)"))
  }

  # check if the number of events matches too
  # TODO: do this more rigorously, this assumes censoring is coded with a 0
  not_cox <- (type == "coxph" && length(setdiff(unique(data[[outcome$status]]), 0)) > 1)
  not_csc <- (type == "CSC" && length(setdiff(unique(data[[outcome$status]]), 0)) < 2)
  if (not_cox || not_csc){
    data_type <- ifelse(not_cox, "competing risks data",
                                 "standard survival data")
    input_type <- ifelse(not_cox, "standard survival data",
                                  "competing risks data")
    if (not_cox)
    stop(tidymess(paste0("Mismatch between the number of events and the type of
      model fit. You appear to have ", data_type, ", but have specified ",
      input_type, " (i.e., type = ",    type, "). For standard
      survival data (one event and possible censoring), use type = \"coxph\" and
      a formula with left hand side of the form Surv(LM, Time, event). For
      competing risks (multiple events and possible censoring), use type =
      \"CSC\" and a formula with left hand side of the form
      Hist(Time, event, LM).")))
  }

  # check if we have a cluster term and extract it
  cluster_check <- as.character(stats::as.formula(formula))[3]
  if (!grepl("cluster", cluster_check)) {
    if (missing(cluster)) {
      stop(tidymess("Did you forget to specify a cluster argument or add a '+
                    cluster(..)' term in your formula for your ID variable? No
                    cluster argument was specified in the formula. Standard
                    errors may be estimated incorrectly."))
    } else {
      id_col <- cluster
    }
  } else {
    cluster <- regmatches(cluster_check,
                          gregexpr("(?<=cluster\\().*?(?=\\))",
                                   cluster_check, perl = TRUE))[[1]]
    id_col <- cluster
  }

  # Check the variables are in the data and as expected
  lhs_vars <- all.vars(formula[[3]])
  # lm_covs <- intersect(sub("_[^_]+$", "", lhs_vars), lm_covs)
  lm_covs <- intersect(unique(sub("_LM\\d+$", "", lhs_vars)), lm_covs)

  # check if all columns y are in the dataframe
  missing_cols <- lhs_vars[!(lhs_vars %in% colnames(data))]
  if (length(missing_cols) > 0) {
    idx <- grep("_\\d$", missing_cols)
    if (length(idx) > 0) {
      stop(tidymess(
        "Columns", paste0("`", paste(missing_cols[idx], collapse = ", "), "`"),
        "are not in the dataframe. NOTE: in `dynamicLM` v2,
        landmark-variables interactions are given by variable_name_LMi
        (vs variable_name_i in v1) and landmark transformations are given
        by LMi (vs LM_i in v1)."
      ))
    } else {
      stop(tidymess(
        "Columns", paste0("`", paste(missing_cols, collapse = ", "), "`"),
        "are not in the dataframe."))
    }
  }

  # Fit the model
  if (type == "coxph") {
    superfm <- survival::coxph(formula, data, method = method, ...)
    models <- list(superfm)
    num_causes <- 1
    superfm$call$data <- data
    cl <- "LMcoxph"

  } else if (type == "CauseSpecificCox" || type == "CSC") {
    if (!requireNamespace("riskRegression", quietly = TRUE))
      stop("Package \"riskRegression\" must be installed to use this function.",
           call. = FALSE)
    superfm <- riskRegression::CSC(formula, data, method = method, ...)
    models <- superfm$models
    num_causes <- length(models)
    superfm$call$data <- data
    cl <- "LMCSC"
  }

  # Get linear predictors
  linear.predictors <-
    t(sapply(1:num_causes, function(c) {
      coefs <- models[[c]]$coefficients
      df <- data[, names(coefs)]
      rowSums(data.frame(mapply(`*`, df, coefs)))
    })
    )

  out <- list(model = superfm,
              type = type,
              w = w,
              end_time = end_time,
              func_covars = func_covars,
              func_lms = func_lms,
              lm_covs = lm_covs,
              all_covs = all_covs,
              outcome = outcome,
              LHS = LHS,
              id_col = id_col,
              lm_col = lm_col,
              linear.predictors = linear.predictors,
              original.landmarks = original.landmarks,
              args = args
  )
  if (x == TRUE) out$data <- lmdata
  class(out) <- c(cl, "dynamicLM")

  return(out)
}
