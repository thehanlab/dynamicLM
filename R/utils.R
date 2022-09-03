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
