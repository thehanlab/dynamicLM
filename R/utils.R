# ----------------------------------------------------------
# find_se_log: Helper function to calculate SE for time varying log HR
# of the form coef[1] + t*coef[2] + t^2*coef[2] + ..
# ----------------------------------------------------------
find_se_log <- function(t, coefs, covar, func_covars){
  if (!requireNamespace("msm", quietly = TRUE)) {
    stop("Package \"msm\" must be installed to use function find_se_log", call. = FALSE)}
  form <- "x1"
  if (length(coefs) > 1){
    for (i in 2:length(coefs)){
      form <- paste(form, sprintf("%s * %f", paste("x",i,sep=""), func_covars[[i-1]](t)),sep=" + ")
    }
  }
  form <- paste("~", form)
  se <- msm::deltamethod(stats::as.formula(form), coefs, covar)
  return(se)
}


# ----------------------------------------------------------
# find_se: Helper function to calculate SE for time varying HR
# of the form coef[1] + t*coef[2] + t^2*coef[2] + ..
# ----------------------------------------------------------
find_se <- function(t, coefs, covar, func_covars){
  if (!requireNamespace("msm", quietly = TRUE)) {
    stop("Package \"msm\" must be installed to use function find_se", call. = FALSE)}
  form <- "x1"
  if (length(coefs) > 1){
    for (i in 2:length(coefs)){
      form <- paste(form, sprintf("%s * %f", paste("x",i,sep=""), func_covars[[i-1]](t)),sep=" + ")
    }
  }
  form <- paste("~ exp(", form,")")
  se <- msm::deltamethod(stats::as.formula(form), coefs, covar)
  return(se)
}


# ----------------------------------------------------------
# replace_na_with_last: Helper function to update a vector by
# replacing any NA value with the previous non-NA value
# ----------------------------------------------------------
replace_na_with_last <- function(x) {
  c(NA, x)[cummax(seq_along(x)*(!is.na(x)))+1]
}

# ----------------------------------------------------------
# getLHS: Function that takes a formula as an input and obtains
# the LHS, and converts it to a Hist object if it is a Surv object
# ----------------------------------------------------------
getLHS <- function(formula){
  LHS = formula[[2]]
  survival_object = LHS[[1]]

  if (length(LHS) == 3){
    exit <- LHS[[2]]
    status <- LHS[[3]]
  } else if (survival_object == "Surv"){
    enter <- LHS[[2]]
    exit <- LHS[[3]]
    status <- LHS[[4]]
  } else if (survival_object == "Hist"){
    exit <- LHS[[2]]
    status <- LHS[[3]]
    enter <- LHS[[4]]
  } else {
    stop("Invalid formula: LHS must be an object of type Hist or Surv")
  }

  exit = deparse(exit)
  status = deparse(status)

  if (length(LHS) == 3) return(paste0("Hist(",exit," ,",status,") ~ 1"))

  enter = deparse(enter)
  return(paste0("Hist(",exit,", ",status,", ",enter,") ~ 1"))
}



# ----------------------------------------------------------
# CSC.extended: cause-specific Cox proportional hazard regression
# extended to accomodate specific arguments to each cause-specific
# Cox model
# TODO: licensing / referencing
# ----------------------------------------------------------
CSC.fixed.coefs <- function(formula,
                            data,
                            cause,
                            cause.specific.coefs,
                            ...){
  fitter <- "coxph"
  surv.type <- "hazard"

  # {{{ formulae & response
  if (inherits(x=formula,what="formula")) formula <- list(formula)
  call <- match.call()
  # get outcome information from formula
  Rform <- update(formula[[1]],".~1")
  response <- eval(Rform[[2]],envir=data)
  if (any(is.na(response)))
    stop("Event history response may not contain missing values")
  time <- response[, "time"]
  status <- response[, "status"]
  event <- prodlim::getEvent(response)
  if ("entry" %in% colnames(response))
    entry <- response[, "entry"]
  else{
    entry <- NULL
  }
  if (any(entry>time)) stop("entry > time detected. Entry time into the study must be strictly greater than outcome time.")
  ## remove event history variables from data
  if(any((this <- match(all.vars(Rform),names(data),nomatch=0))>0)){
    if (data.table::is.data.table(data))
      data <- data[,-this,with=FALSE]
    else
      data <- data[,-this]
  }
  # }}}
  # {{{ sorted unique event times
  eventTimes <- unique(sort(as.numeric(time[status != 0])))
  # }}}
  # {{{ causes
  causes <- prodlim::getStates(response)
  NC <- length(causes)
  if (length(cause.specific.coefs) != NC) stop("There should be one entry for each cause specific argument but there are ", NC, " causes and ",length(cause.specific.args)," elements in cause.specific.args")

  if (length(formula)!=NC[1] && length(formula)>1) stop("Wrong number of formulae. Should be one for each cause ",NC,".")
  if (length(formula)==1) {
    formula <- lapply(1:NC,function(x)formula[[1]])
  }
  # }}}
  # {{{ find the cause of interest
  if (missing(cause)){
    theCause <- causes[1]
  }
  else{
    if ((foundCause <- match(as.character(cause),causes,nomatch=0))==0)
      stop(paste0("Cannot find all requested cause(s) ...\n\n",
                  "Requested cause(s): ", paste0(cause, collapse = ", "),
                  "\n Available causes: ", paste(causes, collapse = ", "),
                  "\n"))
    else{
      theCause <- causes[foundCause]
    }
  }
  otherCauses <- causes[-match(theCause,causes)]
  # }}}
  # {{{ fit Cox models
  CoxModels <- lapply(1:NC,function(x){
    if (x==1)
      causeX <- theCause
    else
      causeX <- otherCauses[x-1]

    statusX <- as.numeric(event==causeX)

    if (is.null(entry))
      workData <- data.frame(time=time,status=statusX)
    else
      workData <- data.frame(time=time,status=statusX,entry=entry)
    if(any(this <- match(names(data),names(workData),nomatch=0)>0)){
      warning(paste("Variables named",paste(names(data)[this],collapse=", "),"in data will be ignored."))
      if (is.data.table(data))
        data <- data[,-this,with=FALSE]
      else
        data <- data[,-this,drop=FALSE]
    }
    workData <- cbind(workData,data)
    if (is.null(entry))
      survresponse <- "survival::Surv(time, status)"
    else
      survresponse <- "survival::Surv(entry, time, status)"
    ## check whether right hand side of formula includes ~.
    allvars <- all.vars(formula[[x]])
    if (any(grepl("^\\.$",allvars))){
      formulaXX <- as.formula(paste0(survresponse,"~."))
    }
    else {
      formulaXX <- update(formula[[x]],paste0(survresponse,"~."))
    }

    args <- list(formulaXX, data = workData)
    extra.args <- list(...)
    fit <- do.call("coxph",c(args, list(x=TRUE, y=TRUE, iter.max=0, init=cause.specific.coefs[[x]]), extra.args))
    fit$call$formula <- formulaXX
    fit$call$data <- workData
    fit
  })
  names(CoxModels) <- paste("Cause",c(theCause,otherCauses))

  # }}}
  out <- list(call=call,
              models=CoxModels,
              response=response,
              eventTimes=eventTimes,
              surv.type=surv.type,
              fitter=fitter,
              theCause=theCause,
              causes=c(theCause,otherCauses))
  class(out) <- "CauseSpecificCox"
  out
}

