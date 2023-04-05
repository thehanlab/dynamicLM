dynamic_lm_helper <- function(formula, type, data, method, w, end_time,
                              func_covars, func_lms, lm_covs, all_covs,
                              outcome, lm_col, original.landmarks, args, ...) {
  # extra LHS of formula
  LHS = getLHS(formula)

  # check if we have a cluster term and extract it
  cluster_check <- as.character(stats::as.formula(formula))[3]
  if (!grepl("cluster", cluster_check)) {
    if (missing(cluster)) {
      stop("Did you forget to specify a cluster argument or add a '+ cluster(..)' term in your formula for your ID variable? No cluster argument was specified in the formula. Standard errors may be estimated incorrectly.")
    } else {
      id_col <- cluster
    }
  } else {
    cluster <- regmatches(cluster_check,
                          gregexpr("(?<=cluster\\().*?(?=\\))",
                                   cluster_check, perl = TRUE))[[1]]
    id_col <- cluster
  }

  lm_covs <- intersect(sub("_[^_]+$", "", all.vars(formula[[3]])), lm_covs)

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
