## TODO: sim.survdata outputs only whole units: user must be aware of this
## TODO: add num.data.frames argument
## TODO: roxygen
## TODO: imports (dplyr?, tidyr?, coxed...)
## TODO: importantly!!!!!! handle longitudinal covariates
##       -> we're not updating covariate values at any point here
## TODO: remove multipier arg

sim_lm_data <- function(supermodel, data, N, hazard.fun = NULL, replace = TRUE,
                        censor = FALSE, multiplier = 1, ...){
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
  Tmax <- w + 1

  if (missing(data)) stop("Argument data must be specified.")
  if (missing(N)) N <- nrow(data)
  # }}}

  # {{{ Create a function to get the hazard function if using the hazard function
  #    of the fitted model
  if (!is.null(hazard.fun)) {
    # instant hazards for each cause-specific Cox model
    instHaz <- lapply(1:num_causes, function(i) {
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

        # risk <- riskScore(models[[i]], s, base_data, func_covars, func_lms)
        haz <- sf$Haz[pred_window] #* exp(risk)
        instHaz <- haz[2:n_times] - haz[1:n_times - 1]
        times <- sf$time[pred_window][2:n_times]
        return(list("iH" = instHaz, "time" = times))
      })
    })

    # function to return hazard function for a given cause and landmark
    hazard.fun <- lapply(1:num_causes, function(cause) {
      lapply(seq_along(lms), function(s) {
        iH <- instHaz[[cause]][[s]]
        func <- function(t) {
          sapply(t, function(ti) {
            idx <- (lms[s] + ti - 1 <= iH$time) & (iH$time <= lms[s] + ti) # add all in one unit before
            sum(iH$iH[idx]) * multiplier
          })
        }
        return (func)
      })
    })
  }
  # }}}

  # {{{ Simulate data
  newdata <- data[sample(1:nrow(data), N, replace = replace), ] # subset of original data
  newdata[[time_col]] <- NA   # this is what we're simulating
  newdata[[event_col]] <- NA  # this is what we're simulating
  idx <- rep(T, N) # inds we need to still simulate survival times for

  i <- 1
  while (i <= length(lms) & sum(idx) != 0) {
    # get curr landmark and data to simulate from
    t <- lms[i]
    curr_data <- newdata[idx, ]

    # simulate times for each cause
    times <- matrix(
      sapply(1:num_causes, function(c) {
        risk_scores_c <- sapply(1:nrow(curr_data), function(row) {
          riskScore(models[[c]], t, curr_data[row, ], func_covars, func_lms)
        })
        risk_scores_c <- data.frame(risk_scores_c)
        coxed::sim.survdata(X = risk_scores_c, T = Tmax, beta = 1, censor = 0,
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
  newdata[idx_na, time_col] <- cutoff
  newdata[idx_na, event_col] <- 0
  # }}}

  # {{{ Censoring

  # TODO: adjust when we have longitudinal data
  # TODO: add covariate dependent censoring (copy sim.survdata?)
  # TODO: add uniform/exponential censoring too
  # TODO: then censor = c("unif", "ind", "dep") ? Or maybe for dep
  #       need to specify a formula

  if (censor == TRUE) {
      args <- list(formula = as.formula(paste0("Hist(", time_col, ", ",
                                               event_col, ") ~ 1")),
                   data = data, reverse = TRUE)
      fit <- do.call(prodlim::prodlim, args)
      print("model fit")
      print(str(fit))
      censoring_times <- sapply(1:N, function(i) {
        # TODO: finding  closest proba could be better
        idx_surv <- which.min(abs(fit$surv - runif(1)))
        idx_time <- fit$surv == fit$surv[idx_surv]

        if (length(idx_time) > 1)
          sample(fit$time[idx_time], 1)
        else
          fit$time[idx_time]
      })
      print(censoring_times)
      newdata[[time_col]] <- pmin(newdata[[time_col]], censoring_times)
      newdata[newdata[[time_col]] == censoring_times, event_col] <- 0
  }
  # }}}

  return(newdata)
}
