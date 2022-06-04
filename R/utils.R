# ----------------------------------------------------------
# find_se_log: Helper function to calculate SE for time varying log HR
# of the form coef[1] + t*coef[2] + t^2*coef[2] + ..
# ----------------------------------------------------------
find_se_log <- function(t, coefs, covar, func_covars){
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
# replace_na_with_last: Helper function to update a vector with NA
# replaces a NA value with the previous non-NA value
# ----------------------------------------------------------
replace_na_with_last<-function(x,p=is.na,d=NA) c(d,x)[cummax(seq_along(x)*(!p(x)))+1]
