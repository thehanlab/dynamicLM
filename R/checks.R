check_evaluation_inputs <- function(
    object,
    times,
    formula,
    data,
    tLM,
    ID_col="ID",
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
    names(object)[(names(object)=="")] <- paste("Model", seq(NF))[(names(object)=="")]
  }
  names(object) <- make.unique(names(object))

  # check no data is input when bootstrapping
  if ((split.method == "bootcv") & (!missing(data))){
    message("Bootstrapping is performed with data from the original model, argument data is ignored.")

  }

  # check elements of object are the correct classes
  # check LMpred are given when not bootstrapping
  # check supermodels have data stored if bootstrapping
  for (i in 1:NF){
    name = names(object)[[i]]
    # check compatibility of object type with bootstrapping and predictions
    if (inherits(object[[i]],"LMpred")){
      if (split.method != "none") {
        stop(paste("Cannot bootstrap with deterministic risk predictions:", name))
      }
      if (!missing(data)) {
        stop(paste0("Cannot use deterministic predictions of type LMpred (",name,") on new data, data argument should not be specified."))
      }

    }
    # TODO: update to include penalized classes
    else if (!(inherits(object[[i]],c("LMCSC", "LMcoxph")))) {
      stop(paste("all prediction models in object must be of class LMCSC, LMcoxph (i.e., output from dynls) or LMpred (i.e., output from predict.dynamicLM) but", name, "is of class",class(object[[i]])))
    }
    else { # We know it is of class LMCSC or LMcoxph
      if (split.method == "bootcv"){
        LMdata <- object[[i]]$data
        if (is.null(LMdata)) {
          model.class = class(object[[i]])
          stop(paste("Prediction object",name,"does not have data stored. In order to bootstrap, objects of class", model.class,"must be fit with argument x=TRUE."))
        }
      }
    }
  }

  # external validation: only need a dataframe
  if (!missing(data)) {
    if (split.method == "bootcv"){
      if (!missing(tLM)){
        warning("LMdataframe objects have LM columns, tLM is ignored.")
      }
    }
    # only need the dataframe
    if (inherits(data,"LMdataframe")) {
      tLM <- data$data[[data$lm_col]]
      data <- data$data
    }
    else if (missing(tLM)) {
      stop("For external validation on new data, argumnet tLM must be specified.")
    }

  } else {
    if (!missing(tLM)){
      warning("tLM is specified without the argument data. Did you mean to use argument times instead? See documentation for more information.")
    }
  }

  # get data if bootstrapping
  if (split.method == "bootcv") {
    LMdata <- object[[1]]$data
    LMdata$data <- LMdata$data[stats::complete.cases(LMdata$data), ]
    data <- LMdata$data
    lm_col <- LMdata$lm_col

    if (!is.null(object[[1]]$ID_col)){ ID_col <- object[[1]]$ID_col }
    if(!(ID_col %in% colnames(data))){
      stop(paste("An ID column that is in the data must be provided;", ID_col, "is not a column."))
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
    unique.inds = unique(data[[ID_col]])
    split.method <- riskRegression::getSplitMethod(split.method, B=B, N=length(unique.inds), M=M, seed)
    B <- split.method$B
    split.idx <- split.method$index
    perform.boot <- !is.null(split.idx)
    if (perform.boot) split.idx <- replace(split.idx, seq_along(split.idx), unique.inds[split.idx])# get IDs of the individuals
  }
  else if (!(split.method %in% c("none", "noplan"))) {
    paste("split.method of type", split.method, "is not supported.")
  }

  ### Create data we need: preds, preds_LM, data ###

  if (inherits(object[[1]],"LMpred")) {
    if (NF > 1) {
      for (i in 2:NF) {
        if (!inherits(object[[i]],"LMpred"))
          stop("All prediction models must either be supermodels or of type LMpred")
      }
    }

    # check data
    preds <- data.frame(lapply(object, function(o) o$preds$risk))
    colnames(preds) <- names(object)
    pred_LMs <- object[[1]]$preds$LM
    data <- object[[1]]$data
    num_preds <- nrow(object[[1]]$preds)

    type <- lapply(object, function(o) ifelse(o$type == "coxph", "coxph", "CSC"))

    if (NF > 1) {
      for (i in 2:NF) {
        if (nrow(object[[i]]$preds) != num_preds)
          stop("number of predictions is not the same for all prediction models.")
        if (!isTRUE(all.equal(object[[i]]$preds$LM, pred_LMs)))
          stop("LM points for individuals across prediction models are not the name.")
      }
    }
  }


  else if (!perform.boot) {
    # TODO: consider including w, extend, silence, complete as args to predict
    if (!missing(data)) preds = lapply(object, function(o) predict.dynamicLM(o, data, tLM, cause))
    else preds = lapply(object, function(o) predict.dynamicLM(o, cause=cause))
    args = match.call()
    args$data = NULL
    args$tLM = NULL
    args$object = preds
    return(eval(args))
  }


  else if (perform.boot){
    # pred.list <- parallel::mclapply(1:B,function(b){
    pred.list <- lapply(1:B,function(b){
      id_train_b <- split.idx[,b]
      id_train_b <- data[[ID_col]] %in% id_train_b

      data_val_b <- data[!id_train_b, ]
      tLMs_b <- data_val_b[[lm_col]]
      outcomes_val_b <- data_val_b[c(outcome$time, outcome$status, lm_col)]

      data_train_b <- LMdata
      data_train_b$data <- data[id_train_b, ]

      preds.b <- do.call("cbind",lapply(1:NF,function(f){
        original_model <- object[[f]]
        args <- original_model$args
        args$data <- NULL
        args$data <- data_train_b
        model.b <- eval(args)
        base_data = 0 * model.b$model$coefficients
        pred.b <- try(
          # TODO: consider including w, extend, silence, complete as args to predict.dynamicLM
          # TODO: use complete=T
          predict.dynamicLM(model.b,newdata=data_val_b,lms=tLMs_b,cause=model.b$cause, complete=F),
          silent = F
        )


        if (inherits(pred.b, "try-error")){
          P <- rep(NA,NROW(data_val_b))
        } else {
          P <- pred.b$preds$risk
        }
        if (sum(is.na(P)) > 0) {
          message(paste("Dropping bootstrap b =",b,"for model",names(object)[[f]],"due to unreliable predictions."))
          return(c())
        }
        else {
          return(P)
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
    type <- lapply(object, function(o) ifelse(inherits(o, "LMcoxph"), "coxph", "CSC"))

    rm(pred.df)
  }

  # check LM times
  if (missing(times))
    times <- unique(pred_LMs)
  else {
    if (!all(times %in% unique(pred_LMs))) {
      times = paste(times, collapse = ",")
      pred_LMs = paste(unique(pred_LMs), collapse = ",")
      stop(paste("arg times (= ",times,") must be a subset of landmark prediction times (= ",pred_LMs,")"))
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
