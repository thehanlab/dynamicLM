## TODO: sim.survdata outputs only whole units: user must be aware of this
## TODO: add num.data.frames argument
## TODO: imports (dplyr?, tidyr?, coxed...)
## TODO: importantly!!!!!! handle longitudinal covariates
##       -> we're not updating covariate values at any point here
## TODO: remove multipier arg
## TODO: roxygen

#' Title
#'
#' @param supermodel
#' @param data
#' @param n
#' @param hazard.fun
#' @param replace
#' @param censor
#' @param multiplier
#' @param tmax
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sim_lm_data <- function(supermodel, data, n, hazard.fun = NULL, replace = TRUE,
                        censor = FALSE, multiplier = 1, tmax, ...) {
  # {{{ Get some variables from the supermodel & check for missing variables
  if (inherits(supermodel, "LMCSC")) {
    models <- supermodel$model$models
  } else if  (inherits(supermodel, "LMcoxph")) { # TODO: check
    models <- list(supermodel$model)
  } else {
    stop("The only models supported are Cox and CSC landmark supermodels.")
  }

  w <- supermodel$w
  time_col <- supermodel$outcome$time
  event_col <- supermodel$outcome$status
  lms <- unique(supermodel$data$data[[supermodel$lm_col]])
  func_covars <- supermodel$func_covars
  func_lms <- supermodel$func_lms
  num_causes <- length(models)
  if (missing(tmax)) tmax <- w + 1

  if (missing(data)) stop("Argument data must be specified.")
  if (missing(n)) n <- nrow(data)
  # }}}

  # {{{ Create a function to get the hazard function if needed
  #    of the fitted model
  # TODO: argument checks for hazard.fun if it is user-specified
  if (inherits(hazard.fun, "logical")) if (hazard.fun == TRUE) {
    # instant hazards for each cause-specific Cox model
    inst_haz <- lapply(1:num_causes, function(i) {
      base_data <- 0 * models[[i]]$coefficients
      sf <- survival::survfit(models[[i]], newdata = base_data)
      sf <- data.frame(time = sf$time, surv = sf$surv, Haz = -log(sf$surv))

      lapply(lms, function(s) {
        pred_window <- (s <= sf$time) & (sf$time <= s + w)
        n_times <- sum(pred_window)

        # risk from all 0 data, only effect is LM terms (affects HR)
        # TODO: NB the LM terms don't exist for us as the [s, s+w) intervals
        #       don't overlap, so can just use the OG hazard function...
        #       doesn't really change much but is worth noting

        haz <- sf$Haz[pred_window] #* exp(risk)
        ihaz <- haz[2:n_times] - haz[1:n_times - 1]
        times <- sf$time[pred_window][2:n_times]
        return(list("iH" = ihaz, "time" = times))
      })
    })

    # function to return hazard function for a given cause and landmark
    hazard.fun <- lapply(1:num_causes, function(cause) {
      lapply(seq_along(lms), function(s) {
        ihaz <- inst_haz[[cause]][[s]]
        func <- function(t) {
          sapply(t, function(ti) {
            # add all in one unit before
            idx <- (lms[s] + ti - 1 <= ihaz$time) & (ihaz$time <= lms[s] + ti)
            sum(ihaz$iH[idx]) * multiplier
          })
        }
        return(func)
      })
    })
  }
  # }}}

  # {{{ Simulate data
  # subset of original data to simulate from
  newdata <- data[sample(seq_along(nrow(data)), n, replace = replace), ]
  # this is what we're simulating
  newdata[[time_col]] <- NA
  newdata[[event_col]] <- NA
  idx <- rep(TRUE, n) # inds we need to still simulate survival times for

  i <- 1
  while (i <= length(lms) && sum(idx) != 0) {
    print(paste("-->", sum(idx)))
    # get curr landmark and data to simulate from
    t <- lms[i]
    curr_data <- newdata[idx, ]

    # simulate times for each cause
    times <- matrix(
      sapply(1:num_causes, function(c) {
        risk_scores_c <- sapply(seq_len(nrow(curr_data)), function(row) {
          riskScore(models[[c]], t, curr_data[row, ], func_covars, func_lms)
        })
        risk_scores_c <- data.frame(risk_scores_c)
        coxed::sim.survdata(X = risk_scores_c, T = tmax, beta = 1, censor = 0,
                            hazard.fun = hazard.fun[[c]][[i]], ...)$data$y
      }),
      ncol = num_causes
    )

    # get min time and corresponding event
    time <- apply(times, 1, min)
    event <- apply(times, 1, which.min)

    # find which time and event we update
    keep <- time <= w
    idx_to_change <- idx
    idx_to_change[idx] <- idx[idx] & keep

    # update rel. time and event
    newdata[idx_to_change, time_col] <- time[keep] + t
    newdata[idx_to_change, event_col] <- event[keep]

    # store which still need to be simulated
    idx[idx] <- idx[idx] & !keep
    i <- i + 1
  }

  # censor anyone left
  cutoff <- max(lms) + w
  idx_na <- is.na(newdata[[time_col]]) | is.na(newdata[[event_col]])
  # print(newdata[idx_na,])
  # print(cutoff)
  newdata[idx_na, time_col] <- cutoff
  newdata[idx_na, event_col] <- 0
  # print(newdata[idx_na,])
  # }}}

  # {{{ Censoring

  # TODO: adjust when we have longitudinal data
  # TODO: add covariate dependent censoring (copy sim.survdata?)
  # TODO: add uniform/exponential censoring too (pick one/both?)
  # TODO: then censor = c("unif", "ind", "dep") ? Or maybe for dep
  #       need to specify a formula
  if (is.numeric(censor)) {
    censor_time <- rexp(n, rate = censor)
    # censor_time <- runif(n) * censor
  } else if (censor == TRUE) {
      args <- list(formula = as.formula(paste0("Hist(", time_col, ", ",
                                               event_col, ") ~ 1")),
                   data = data, reverse = TRUE)
      fit <- do.call(prodlim::prodlim, args)
      # print("model fit")
      # print(str(fit))
      censor_time <- sapply(1:n, function(i) {
        # TODO: finding  closest proba could be better
        idx_surv <- which.min(abs(fit$surv - runif(1)))
        idx_time <- fit$surv == fit$surv[idx_surv]

        if (length(idx_time) > 1)
          sample(fit$time[idx_time], 1)
        else
          fit$time[idx_time]
      })
  }
  if (censor != FALSE) {
    # print(censor_time)
    newdata[[time_col]] <- pmin(newdata[[time_col]], censor_time)
    newdata[newdata[[time_col]] == censor_time, event_col] <- 0
  }
  # }}}

  return(newdata)
}
