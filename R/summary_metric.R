#' Obtain the summary metric from landmark-specific estimates
#'
#' @param metric "AUC" or "Brier"
#' @param df_t time-dependent table of scores (AUC(s,t) or Brier(s,t)), i.e.,
#' a smaller table with entries for each landmark and model used for
#' landmark-estimates and their standard error
#'
#' for example:
#'       tLM      model   times      Brier          se      lower      upper b
#'   1:   0 Null model 59.9999 0.06302245 0.008262163 0.04682891 0.07921599 1
#'   2:   0  dynamicLM 59.9999 0.06140585 0.007931126 0.04586112 0.07695057 1
#'   3:   6 Null model 65.9999 0.08492919 0.011250007 0.06287958 0.10697880 1
#'   4:   6  dynamicLM 65.9999 0.08126512 0.010849048 0.06000138 0.10252887 1
#'   5:  12 Null model 71.9999 0.11151696 0.013936209 0.08420249 0.13883143 1
#'   6:  12  dynamicLM 71.9999 0.10554550 0.013541907 0.07900385 0.13208715 1
#'
#' @param df_c table of contrast of scores.
#'
#' for example:
#'    tLM   times     model  reference  delta.Brier         se        lower        upper         p
#' 1:   0 59.9999 dynamicLM Null model -0.001616602 0.00119091 -0.003950743 0.0007175397 0.1746381
#'
#' @param df_iid df_iid contains the iid decomposition of the score, it is a
#' larger table with entries for each individual, who appear multiple times for
#' different landmarks and different models
#'
#' for example:
#'     tLM ID model   times      IF.Brier b
#'  1:   0  1     0 59.9999 -5.853094e-02 1
#'  2:   0  2     0 59.9999 -8.262274e-05 1
#'  3:   0  3     0 59.9999 -5.852535e-02 1
#'  ...
#'
#' @param conf_int Coverage level of the confidence interval.
#' @param object Either fitted supermodel or risk predictions.
#' @param id_col Column name of the ID (only needed if we incorporate
#'   weighted metrics)
#' @param B Number of bootstrap replicates
#' @param se.fit If FALSE or 0, no standard errors are calculated.
#'
#' @return An list with entries `score` and optionally `contrasts` if
#'   contrasts were calculated.
summary_metric <- function(metric,
                           df_t,
                           df_c,
                           df_iid,
                           conf_int,
                           object,
                           id_col,
                           B,
                           se.fit
                           ) {
  if (!metric %in% c("Brier", "AUC"))
    stop("Only Brier and AUC are handled.")

  if_col <- paste0("IF.", metric)
  delta_col <- paste0("delta.", metric)
  lms <- unique(df_t$tLM)
  num_lms <- length(lms)
  sample_size <- length(unique(df_iid$ID))

  if (any(is.na(df_t[[if_col]]))) {
    stop(tidymess(paste0(
      "There are missing values in the calculated time-dependent ", metric, ",
       the summary metric cannot be calculated.")))
  }

  get_contrasts <- TRUE
  if (is.null(df_c)) get_contrasts <- FALSE
  else if (ncol(df_c) <= 2) get_contrasts <- FALSE

  if (!missing(conf_int)) alpha <- 1 - conf_int

  weight_vector <- rep(1, num_lms)

  ### No bootstrapping: use i.i.d. decomposition
  if (B == 1) {
    #{{{ 1. get the average of the score from df_t by model over times
    #       and of the contrasts from df_c by model over times
    summary_score <- df_t[, lapply(.SD, mean),
                          by = "model", .SDcols = metric]

    if (get_contrasts) {
      contrasts <- merge(unique(df_c[,c("model", "reference")]), summary_score,
                         by.x = "model", by.y = "model")
      contrasts <- merge(contrasts, summary_score,
                         by.x="reference", by.y="model",
                         suffixes = c(".model", ".reference"))
      contrasts[, delta := get(paste0(metric, ".model")) -
                  get(paste0(metric, ".reference"))]
      contrasts <- contrasts[, .(model, reference, delta)]
      data.table::setnames(contrasts, "delta", delta_col)
    }
    #}}}

    #{{{ 2. get the se/CIs
    if (se.fit == TRUE) {
      #{{{{ 2.1. get the covariance of the iid decomp by model across times
      cov_score <- lapply(unique(df_iid$model), function(m) {
        df <- df_iid[df_iid$model == m, c("ID", "tLM", if_col), with = FALSE]

        subsample_sizes <- table(df$tLM)
        if (subsample_sizes[1] != sample_size)
          stop("Subsample size mismatch.")

        adjustment <- subsample_sizes / sample_size
        df <- data.table::dcast(df, ID ~ tLM, value.var = if_col)
        df[, "ID" := NULL]
        df[is.na(df),] <- 0
        df <- sweep(df, MARGIN = 2, STATS = adjustment, FUN = "/")
        cov(df)
      })

      if (get_contrasts) {
        cov_score_contrasts <- lapply(1:nrow(contrasts), function(i) {
          m1 <- contrasts[i, ]$reference
          m2 <- contrasts[i, ]$model
          df1 <- df_iid[df_iid$model == m1,
                        c("ID", "tLM", if_col), with = FALSE]
          df2 <- df_iid[df_iid$model == m2,
                        c("ID", "tLM", if_col), with = FALSE]
          df <- merge(df1, df2, by=c("ID", "tLM"),
                      suffixes = c(".df1", ".df2"))
          df[, (if_col) :=
                  get(paste0(if_col, ".df1")) - get(paste0(if_col, ".df2"))]
          df <- df[, c("ID", "tLM", if_col), with = FALSE]

          subsample_sizes <- table(df$tLM)
          if (subsample_sizes[1] != sample_size)
            stop("Subsample size mismatch.")

          df <- data.table::dcast(df, ID ~ tLM, value.var = if_col)
          df[, ID := NULL]
          df[is.na(df),] <- 0
          adjustment <- subsample_sizes / sample_size
          df <- sweep(df, MARGIN = 2, STATS = adjustment, FUN = "/")
          cov(df)
        })
      }
      #}}}}

      #{{{{ 2.2. apply the delta method to get a confidence interval
      get_se <- function(cov) {
        nabla_g <- weight_vector / sum(weight_vector) # gradient
        var <- nabla_g %*% cov %*% nabla_g            # delta method for var
        sqrt(var / sample_size)                       # variance -> std error
      }
      se_score <- sapply(cov_score, get_se)
      if (get_contrasts) {
        se_score_contrasts <- sapply(cov_score_contrasts, get_se)
      }
      #}}}}

      #{{{{ 2.3. Make a data.table with model, mean, se, lower, upper
      summary_score[, `:=`(
        se = se_score,
        lower = get(metric) - stats::qnorm(1 - alpha / 2) * se_score,
        upper = get(metric) + stats::qnorm(1 - alpha / 2) * se_score
      )]
      if (get_contrasts) {
        contrasts[, se := se_score_contrasts]
        contrasts[, `:=`(
          lower = get(delta_col) - stats::qnorm(1 - alpha / 2) * se,
          upper = get(delta_col) + stats::qnorm(1 - alpha / 2) * se,
          p = 2 * stats::pnorm(abs(get(delta_col) / se), lower.tail=FALSE)
        )]
      }
      #}}}}
    }
    #}}}

  } else { # B > 1 -- bootstrapping: use empirical values
    summary_score1 <- df_t[, lapply(.SD, mean),
                          by = c("model", "bootstrap"), .SDcols = metric]
    summary_score <- clean_bootstraps(summary_score1, metric, alpha,
                                      se.fit = se.fit, summary = TRUE)

    if (get_contrasts) {
      contrasts <- merge(unique(df_c[,c("model", "reference", "bootstrap")]),
                         summary_score1,
                         by.x = c("model", "bootstrap"),
                         by.y = c("model", "bootstrap"))
      contrasts <- merge(contrasts, summary_score1,
                         by.x=c("reference", "bootstrap"),
                         by.y=c("model", "bootstrap"),
                         suffixes = c(".model", ".reference"))
      contrasts[, delta := get(paste0(metric, ".model")) -
                  get(paste0(metric, ".reference"))]
      contrasts <- contrasts[, .(model, reference, delta, bootstrap)]
      data.table::setnames(contrasts, "delta", delta_col)

      contrasts <- clean_bootstraps(contrasts, paste0("delta.", metric), alpha,
                                    se.fit = se.fit, contrasts = TRUE,
                                    summary = TRUE)
    }
  }

  if (get_contrasts)
    return(list(score = summary_score, contrasts = contrasts))
  else
    return(list(score = summary_score))
}


### Code to add alternative weights ###

### * Add argument: weights = FALSE, # (FALSE, TRUE, "survival", "km")
        #' @param weights What kind of weighting to use? Default is NULL which
        #'   gives a simple average. Setting weights to TRUE or "km" or
        #'   "survival" weights estimates by the probability of survival
        #'   P(T > s).
### * The code below would replace the line of `weight_vector <- ...`
### * Add weighted = weights to the outlist
### * Add argument weights = FALSE, to score()
### * Would need to add something like below to checks in score
        # if tolower(weighted) %in% c("km", "survival") || weighted == TRUE
        #   -> then either need preds (to use p1$data with p1$outcome)
        #   -> or need supermodel to be fit with x = TRUE ... (TODO)

# #{{{ 1. get the (weighted) average of the score from df_t by model over times
# #       and of the contrasts from df_c by model over times
# #       First set the weights, then calculate the average
# #       * if (is.null(weights)): use simple weighted average
# #       * if weights is "km" or "survival" or TRUE: weight points by P(T > s)
# #         estimated by Kaplan-Meier
# #       * else user provided function (TODO: maybe later)
#
# if (weights == FALSE) {
#   weight_vector <- rep(1, num_lms)
# } else if (tolower(weights) %in% c("km", "survival") ||
#            weights == TRUE) {
#   # TODO: double check that we've checked that all are the same
#   # i.e., data and outcome of the different entries of object (a list)
#   #       are the same
#
#   outcome <- object[[1]]$outcome
#   type <- object[[1]]$type
#
#   # extract survival data
#   if (inherits(object[[1]], "LMpred")) {
#     # TODO: add id_col to LMpred objects (exists for supermodel)
#
#     # Check id_col
#     if (!(id_col %in% colnames(object[[1]]$data))) {
#       stop(tidymess(paste("An ID column that is in the data must be provided,
#                           you can set this with argument id_col.")))
#     }
#
#     data <- object[[1]]$data[, c(id_col, outcome$time, outcome$status)]
#     data <- data.table::as.data.table(data)
#
#   } else if (FALSE) {
#     # TODO:
#     # check if LMCSC or LMCox (double check classes)
#     # and if it is then check if fit with x = TRUE
#     # otherwise throw an error
#
#   } else {
#     stop("object must be XXXX")
#   }
#
#   # revert from stacked data frame to single times and events
#   data <- data[, lapply(.SD, max), by = id_col,
#                .SDcols = c(outcome$time, outcome$status)]
#   data.table::setnames(data, c(outcome$time, outcome$status),
#                        c("time", "status"))
#
#   # set weights
#   args <- list(formula = stats::as.formula("Hist(time, status) ~ 1"),
#                data = data)
#   fit <- do.call(prodlim::prodlim, args)
#   weight_vector <- predict(fit, times = lms, type = "surv")
#
# } else {
#   # TODO: (later) allow for a more flexible weighting scheme from the user?
#   #       maybe they can input a function get_weights(tLM) which we can
#   #       query for each landmark...
#   #       this really complicates things for the delta method. May not be
#   #       feasible.
#   stop(tidymess("Only simple averages (weights = FALSE) or probability of
#                 survival P(T>s) (weights = TRUE/\"km\"/\"survival\") may be
#                 used."))
# }
# summary_score <- df_t[, lapply(.SD, mean), # weighted.mean, w = weight_vector
#                       by = "model", .SDcols = metric]
