# ----------------------------------------------------------
# find_se: Helper function to calculate SE for time varying log HR
# of the form coef[1] + t*coef[2] + t^2*coef[2] + ..
# ----------------------------------------------------------
find_se <- function(t, coefs, covar, func_covars){
  form <- "x1"
  if (length(coefs) > 1){
    for (i in 2:length(coefs)){
      form <- paste(form, sprintf("%s * %f", paste("x",i,sep=""), func_covars[[i]](t)),sep=" + ")
    }
  }
  form <- paste("~", form)
  se <- msm::deltamethod(stats::as.formula(form), coefs, covar)
  return(se)
}
