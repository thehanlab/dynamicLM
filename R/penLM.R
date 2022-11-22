penLM <- function(x, y, LMdata, xcols, ...) {
  use_LMdata <- T
  # check which data inputs are provided
  if (missing(LMdata)) use_LMdata <- F
  if (!use_LMdata){
    if (missing(x)) stop("argument x is missing with no default, or provide LMdata")
    if (class(x)[1] == "LMdataframe"){
        if (!missing(xcols)) {
          return(penLM(LMdata = x, xcol = xcol...))
        } else if (!missing(y)) {
          if (class(y) == "character") return(penLM(LMdata = x, xcols = y, ...))
          else stop("Inputs are mismatched. Arguments (x, y) should be of type (matrix, Surv object) or arguments (LMdata, xcols) should be (LMdataframe, vector of column names)")
        }
        else return(penLM(LMdata = x, ...))
    }
    if (missing(y)) stop("argument y is missing with no default, or provide LMdata")
    if (!(class(y) %in% c("Surv", "Hist"))) stop("argument y should be a Surv or Hist object")
  }
  if (use_LMdata) {
    if (class(LMdata) == "LMdataframe"){
      if (!missing(x)) message("Argument x was provided but is redundent with argument LMdata. Ignoring x.")
      if (!missing(y)) message("Argument y was provided but is redundent with argument LMdata. Ignoring y.")
    }
    else{
      stop("Argument LMdata must be of class LMdataframe. This can be created using functions cutLMsuper and addLMtime.")
    }
  }

  # if using an LMdataframe, create x and y for glmnet
  if (use_LMdata){
    entry = LMdata$LMdata[[LMdata$LM_col]]
    exit = LMdata$LMdata[[LMdata$outcome$time]]
    status = LMdata$LMdata[[LMdata$outcome$status]]
    y <- Hist(exit, status, entry)

    if (missing(xcols)) {
      if (!is.null(LMdata$allLMcovars)) xcols <- LMdata$allLMcovars
      else {
        all_cols = colnames(LMdata$LMdata)
        xcols <- all_cols[!(all_cols %in% c(LMdata$LM_col, LMdata$outcome$time, LMdata$outcome$status))]
      }
    }
    x <- as.matrix(LMdata$LMdata[xcols])
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
  models <- lapply(y, function(yi) {
    glmnet(x = x, y = yi, family = "cox", ...)
  })
  if (length(models) > 1){
    attr(models, "survival.type") = "competing.risk"
  } else {
    attr(models, "survival.type") = "survival"
  }

  class(models) = "LMpen"
  return(models)
}

