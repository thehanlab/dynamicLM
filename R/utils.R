# ----------------------------------------------------------
# tidymess: Helper function to make a string more readable by
# wrapping it
# Reference: https://stackoverflow.com/questions/45693010/how-do-you-format-multiline-r-package-messages
# ----------------------------------------------------------
tidymess <- function(..., prefix = "\n", initial = "") {
  strwrap(..., prefix = prefix, initial = initial)
}

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
