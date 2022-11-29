check_evaluation_inputs <- function(
    object,
    times,
    formula,
    data,
    lms,
    id_col="ID",
    split.method = "none",
    B = 1,
    M,
    cores = 1,
    seed,
    cause,
    ...) {

  ### Check inputs ###

  # first checks
  split.method <- tolower(split.method)
  if (split.method == "none") B <- 1

  if (!missing(cores)){
    message("Argument cores is unused. Parallel implementation is not yet implemented.")
  }

  # check naming of objects
  NF <- length(object)
  if (is.null(names(object))) {
    names(object) <- paste("Model", seq(NF))
  } else {
    idx_replace <- (names(object)=="")
    names(object)[idx_replace] <- paste("Model", seq(NF))[idx_replace]
  }
  names(object) <- make.unique(names(object))

  # check no data is input when bootstrapping
  if ((split.method == "bootcv") & (!missing(data))){
    message("Bootstrapping is performed with data from the original model, ",
            "argument data is ignored.")

  }

  # check elements of object are the correct classes
  # check LMpred are given when not bootstrapping
  # check supermodels have data stored if bootstrapping
  for (i in 1:NF){
    name = names(object)[[i]]
    # check compatibility of object type with bootstrapping and predictions
    if (inherits(object[[i]],"LMpred")){
      if (split.method != "none") {
        stop(paste("Cannot bootstrap with deterministic risk predictions:",
                   name))
      }
      if (!missing(data)) {
        stop(paste0("Cannot use deterministic predictions of type LMpred (",
                name,") on new data, data argument should not be specified."))
      }

    }
    # TODO: update to include penalized classes
    else if (!(inherits(object[[i]],c("LMCSC", "LMcoxph")))) {
      stop(paste("all prediction models in object must be of class LMCSC,",
                 "LMcoxph (i.e., output from dynamic_lm) or LMpred (i.e.,",
                 "output from predict.dynamicLM) but", name, "is of class",
                 class(object[[i]])))
    }
    else { # We know it is of class LMCSC or LMcoxph
      if (split.method == "bootcv"){
        LMdata <- object[[i]]$data
        if (is.null(LMdata)) {
          model.class = class(object[[i]])
          stop(paste("Prediction object",name,"does not have data stored.",
                     "In order to bootstrap, objects of class", model.class,
                     "must be fit with argument x=TRUE."))
        }
      }
    }
  }

  # external validation: only need a dataframe
  if (!missing(data)) {
    if (split.method == "bootcv"){
      if (!missing(lms)){
        warning("LMdataframe objects have LM columns, argument lms is ignored.")
      }
    }
    # only need the dataframe
    if (inherits(data,"LMdataframe")) {
      lms <- data$data[[data$lm_col]]
      data <- data$data
    }
    else if (missing(lms)) {
      stop("For external validation with new data, argument lms must be given.")
    }

  } else {
    if (!missing(lms)){
      warning("lms is specified without the argument data.",
              "Did you mean to use argument times instead?",
              "See documentation for more information.")
    }
  }

  # get data if bootstrapping
  if (split.method == "bootcv") {
    LMdata <- object[[1]]$data
    LMdata$data <- LMdata$data[stats::complete.cases(LMdata$data), ]
    data <- LMdata$data
    lm_col <- LMdata$lm_col

    if (!is.null(object[[1]]$id_col)){ id_col <- object[[1]]$id_col }
    if(!(id_col %in% colnames(data))){
      stop(paste("An ID column that is in the data must be provided;", id_col,
                 "is not a column."))
    }
  }

  # more checks
  w <- object[[1]]$w
  if (NF > 1) {
    for (i in 2:NF) {
      if (object[[i]]$w != w)
        stop("prediction window w is not the same for all prediction models.")
    }
  }
  outcome <- object[[1]]$outcome
  if (NF > 1) {
    for (i in 2:NF) {
      for (j in 1:length(outcome)){
        if (object[[i]]$outcome[[j]] != outcome[[j]])
          stop("outcome is not the same for all prediction models.")
      }
    }
  }
  if (missing(cause)){
    if (inherits(object[[1]],"LMpred"))
      cause <- as.numeric(object[[1]]$cause)
    else
      cause <- as.numeric(object[[1]]$model$theCause)
  }
  if (is.null(cause)){
    cause <- 1
    indicator <- 1
  }
  indicator <- cause
  if (missing(formula))
    formula <- object[[1]]$LHS


  ### Split method ###

  perform.boot <- F
  if (split.method == "bootcv"){
    unique.inds = unique(data[[id_col]])
    split.method <- riskRegression::getSplitMethod(split.method, B = B,
                                                   N = length(unique.inds),
                                                   M = M, seed)
    B <- split.method$B
    split.idx <- do.call("cbind", lapply(1:B, split.method$index))
    perform.boot <- !is.null(split.idx)
    if (perform.boot)
      split.idx <- replace(split.idx, seq_along(split.idx),
                           unique.inds[split.idx]) # get IDs of the individuals
  }
  else if (!(split.method %in% c("none", "noplan"))) {
    paste("split.method of type", split.method, "is not supported.")
  }

  ### Create data we need: preds, preds_LM, data ###

  if (inherits(object[[1]],"LMpred")) {
    if (NF > 1) {
      for (i in 2:NF) {
        if (!inherits(object[[i]],"LMpred"))
          stop("All prediction models must be supermodels or of type LMpred")
      }
    }

    # check data
    preds <- data.frame(lapply(object, function(o) o$preds$risk))
    colnames(preds) <- names(object)
    pred_LMs <- object[[1]]$preds$LM
    data <- object[[1]]$data
    num_preds <- nrow(object[[1]]$preds)

    type <- lapply(object,
                   function(o) ifelse(o$type == "coxph", "coxph", "CSC"))

    if (NF > 1) {
      for (i in 2:NF) {
        if (nrow(object[[i]]$preds) != num_preds)
          stop("Number of predictions is not the same for all prediction models.")
        if (!isTRUE(all.equal(object[[i]]$preds$LM, pred_LMs)))
          stop("LM points for individuals across prediction models are not the name.")
      }
    }
  }


  else if (!perform.boot) {
    # TODO: consider including w, extend, silence, complete as args to predict
    if (!missing(data))
      preds = lapply(object, function(o) predict.dynamicLM(o, data, lms, cause))
    else preds = lapply(object, function(o) predict.dynamicLM(o, cause=cause))
    args = match.call()
    args$data = NULL
    args$lms = NULL
    args$object = preds
    return(eval(args))
  }


  else if (perform.boot){
    # pred.list <- parallel::mclapply(1:B,function(b){
    pred.list <- lapply(1:B,function(b){
      id_train_b <- split.idx[,b]
      id_train_b <- data[[id_col]] %in% id_train_b

      data_val_b <- data[!id_train_b, ]
      tLMs_b <- data_val_b[[lm_col]]
      outcomes_val_b <- data_val_b[c(outcome$time, outcome$status, lm_col)]

      data_train_b <- LMdata
      data_train_b$data <- data[id_train_b, ]

      preds.b <- do.call("cbind", lapply(1:NF, function(f) {
        original_model <- object[[f]]
        args <- original_model$args
        args$lmdata <- NULL
        args$lmdata <- data_train_b
        model.b <- eval(args, envir = globalenv())
        base_data <- 0 * model.b$model$coefficients

        pred.b <- try(
          # TODO: consider including w, extend, silence, complete
          #       as args to predict.dynamicLM
          # TODO: use complete = T
          predict.dynamicLM(model.b, newdata = data_val_b, lms = tLMs_b,
                            cause = model.b$cause, complete = FALSE),
          silent = F
        )

        if (inherits(pred.b, "try-error")) {
          return(c())
        } else {
          P <- pred.b$preds$risk
          if (sum(is.na(P)) > 0) {
            message(paste("Dropping bootstrap b =", b, "for model",
                          names(object)[[f]], "due to unreliable predictions."))
            return(c())
          } else {
            return(P)
          }
        }
      }))

      if (length(preds.b) > 0) {
        colnames(preds.b) <- names(object)
        return(cbind(outcomes_val_b, preds.b, b=rep(b, nrow(outcomes_val_b))))
      }
    })#, mc.cores=cores)

    pred.df <- do.call("rbind",pred.list)
    pred.df <- pred.df[stats::complete.cases(pred.df), ] # TODO: remove/fix
    preds <- pred.df[names(object)]
    pred_LMs <- pred.df[[lm_col]]
    data <- pred.df[c(outcome$time, outcome$status, lm_col, "b")]
    num_preds <- nrow(data)
    type <- lapply(object,
                   function(o) ifelse(inherits(o, "LMcoxph"), "coxph", "CSC"))

    rm(pred.df)
  }

  # check LM times
  if (missing(times))
    times <- unique(pred_LMs)
  else {
    if (!all(times %in% unique(pred_LMs))) {
      times = paste(times, collapse = ",")
      pred_LMs = paste(unique(pred_LMs), collapse = ",")
      stop(paste("arg times (= ", times,
                 ") must be a subset of landmark prediction times (= ",
                 pred_LMs, ")"))
    }
  }
  out = list(
    object = object,
    times = times,
    preds = preds,
    pred_LMs = pred_LMs,
    data = data,
    NF = NF,
    w = w,
    formula = formula,
    cause = cause,
    indicator = indicator
  )

  return(out)
}


check_penLM_inputs <- function(x, y, lmdata, xcols, parent_func, ...) {
  parent_func <- eval(parent_func)
  use_lmdata <- TRUE
  # check which data inputs are provided
  if (missing(lmdata)) use_lmdata <- FALSE
  if (!use_lmdata){
    if (missing(x)) 
      stop("argument x is missing with no default, or provide lmdata")
    if (class(x)[1] == "LMdataframe"){
      if (!missing(xcols)) {
        return(parent_func(lmdata = x, xcol = xcol, ...))
      } else if (!missing(y)) {
        if (class(y) == "character") 
          return(parent_func(lmdata = x, xcols = y, ...))
        else 
        `stop("Inputs are mismatched. Arguments (x, y) should be of type (matrix, Surv object) or arguments (lmdata, xcols) should be (LMdataframe, vector of column names)")
      }
      else return(parent_func(lmdata = x, ...))
check_penLM_inputs <- function(x, y, LMdata, xcols, ID_col=NULL, parent_func, CV=FALSE, ...){
  parent_func <- eval(parent_func)
  # check which data inputs are provided
  # i.e., are all provided? & can we replace LMdata and xcols with x and y?
  use_LMdata <- T
  if (missing(LMdata)) use_LMdata <- F
  if (!use_LMdata){
    if (missing(x)) stop("argument x is missing with no default, or provide LMdata")
    if (class(x)[1] == "LMdataframe"){
      if (!missing(xcols)) {
        return(parent_func(LMdata = x, xcol = xcol, ID_col=ID_col, ...))
      } else if (!missing(y)) {
        if (class(y) == "character") return(parent_func(LMdata = x, xcols = y, ID_col=ID_col, ...))
        else stop("Inputs are mismatched. Arguments (x, y) should be of type (matrix, Surv object) or arguments (LMdata, xcols) should be (LMdataframe, vector of column names)")
      }
      else return(parent_func(LMdata = x, ID_col=ID_col, ...))
    }
    if (missing(y)) stop("argument y is missing with no default, or provide lmdata")
    if (!(class(y) %in% c("Surv", "Hist"))) stop("argument y should be a Surv or Hist object")
  }
  if (use_lmdata) {
    if (class(lmdata) == "LMdataframe"){
      if (!missing(x)) message("Argument x was provided but is redundent with argument lmdata. Ignoring x.")
      if (!missing(y)) message("Argument y was provided but is redundent with argument lmdata. Ignoring y.")
    }
    else{
      stop("Argument lmdata must be of class LMdataframe. This can be created using functions stack_data and add_interactions.")
    }
  }

  # get IDs if performing cross-validation
  if (CV) {
    # get ID_col if parent function is cv.penLM
    if (is.null(ID_col) & !use_LMdata) {
      stop("argument ID_col must be provided when using arguments x and y.")
    }
    else if (is.null(ID_col) & use_LMdata){
      ID_col <- LMdata$ID_col
      if (length(LMdata$LMdata[[ID_col]]) == 0)
        stop("The extraction of an ID column from LMdata was unsuccessful, provide argument ID_col")
    }
    # use ID_col to extract IDs
    if (use_LMdata) IDs <- LMdata$LMdata[[ID_col]]
    else {
      IDs <- x[,ID_col]
      if (class(ID_col) == "numeric") x <- x[,-ID_col]
      else x <- x[,colnames(x) != ID_col]
    }
  }

  # if using an LMdataframe, create x and y for glmnet
  if (use_lmdata){
    entry = lmdata$data[[lmdata$LM_col]]
    exit = lmdata$data[[lmdata$outcome$time]]
    status = lmdata$data[[lmdata$outcome$status]]
    y <- Hist(exit, status, entry)

    if (missing(xcols)) {
      if (!is.null(lmdata$all_covs)) xcols <- lmdata$all_covs
      else {
        all_cols = colnames(lmdata$data)
        xcols <- all_cols[!(all_cols %in% c(lmdata$LM_col, lmdata$outcome$time, lmdata$outcome$status, ID_col))]
      }
    }
    x <- as.matrix(lmdata$data[xcols])
  }

  # check if x is competing risks (CR) or not
  # if CR create the correct number of x's: one for each cause
  if (class(y) == "Surv") {
    if (! ("start" %in% attr(y, "dimnames")[[2]])) stop("There is no left-truncated data, which is unusual for a landmark supermodel. Did you forget to include an entry time?")
    y <- list(y)
  }
  else if (class(y) == "Hist") {
    censor_type = attr(y, "cens.type")
    if (censor_type != "rightCensored") stop(paste("Only right-censoring is currently supported, not type", censor_type))

    matrix_cols <- attr(y, "dimnames")[[2]]
    if (sum(c("L", "R") %in% matrix_cols) > 0) stop("Hist object has interval censored times which is not supported.")
    if (sum(c("from", "to") %in% matrix_cols) > 0) stop("Hist object has transition states and multi state models are not supported.")
    if (! ("entry" %in% attr(y, "dimnames")[[2]])) stop("There is no left-truncated data, which is unusual for a landmark supermodel. Did you forget to include an entry time?")

    states <- attr(y, "states")
    if (attr(y, "model") == "survival") {
      y <- list(Surv(y[, 1], y[, 2], y[, 3]))
    } else {
      y <- lapply(1:length(states), function(i) {
        entry = y[, 1]
        exit = y[, 2]
        status = y[, 4] == states[i]
        return(Surv(entry, exit, status))
      })
    }
  }

  if (missing(lmdata)) lmdata <- NULL
  if (missing(xcols)) xcols <- NULL

  out = list(
    x = x,
    y = y,
    lmdata = lmdata,
    xcols = xcols
  )
  if (CV) {
    out$IDs <- IDs
    out$ID_col <- ID_col
  }
  return(out)
}
