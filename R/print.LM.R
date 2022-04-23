#' Print function for object of class LM.data.frame
#'
#' @param LMdata Object of class LM.data.frame
#'
#' @return Printed output.
#' @export
#'
#' @examples
#'
print.LM.data.frame <- function(LMdata){
  cat("$LMdata\n")
  print(utils::head(LMdata$LMdata))
  if(nrow(LMdata$LMdata) > 5) cat(paste0(" [ omitted ",nrow(LMdata$LMdata)-5," rows ]\n"))
  cat("\n")
  cat("\n")

  cat("$outcome\n")
  names.outcome = names(LMdata$outcome)
  for (i in 1:length(LMdata$outcome)){
    if (is.null(names.outcome[i])) label <- paste0("[[",i,"]]")
    else label <- paste0("$",names.outcome[i])
    cat(paste0("$outcome",label,"\n"))
    print(LMdata$outcome[[i]])
    cat("\n")
  }

  cat("$w\n")
  print(LMdata$w)
  cat("\n")

  cat("$end_time\n")
  print(LMdata$end_time)
  cat("\n")

  names.LMdata = names(LMdata)

  if("func_covars" %in% names.LMdata){
    cat("$func_covars\n")
    names.fc = names(LMdata$func_covars)
    for (i in 1:length(LMdata$outcome)){
      if (is.null(names.fc[i])) label <- paste0("[[",i,"]]")
      else paste0("$",names.fc[i])
      cat(paste0("$func_covars$",label,"\n"))
      print(LMdata$func_covars[[i]])
      cat("\n")
    }
  }
  if("func_LMs" %in% names.LMdata){
    cat("$func_LMs\n")
    names.fc = names(LMdata$func_LMs)
    for (i in 1:length(LMdata$outcome)){
      if (is.null(names.fc[i])) label <- paste0("[[",i,"]]")
      else paste0("$",names.fc[i])
      cat(paste0("$func_LMs$",label,"\n"))
      print(LMdata$func_LMs[[i]])
      cat("\n")
    }
  }
  if("LMcovars" %in% names.LMdata){
    cat("$LMcovars\n")
    print(LMdata$LMcovars)
    cat("\n")
  }
  if("allLMcovars" %in% names.LMdata){
    cat("$allLMcovars\n")
    print(LMdata$allLMcovars)
    cat("\n")
  }
  if("LM_col" %in% names.LMdata){
    cat("$LM_col\n")
    print(LMdata$LM_col)
    cat("\n")
  }
}
