cv.penLM <- function(x, y, LMdata, xcols, measure="deviance",...) {

  checked_input <- match.call()
  checked_input$parent_func = quote(cv.penLM)
  checked_input[[1L]] <- quote(check_penLM_inputs)
  checked_input <- eval(checked_input, parent.frame())

  if(class(checked_input) == "cv.penLM") return(checked_input)

  x = checked_input$x
  y = checked_input$y
  LMdata = checked_input$LMdata
  xcols = checked_input$xcols

  models <- lapply(y, function(yi) {
    cv.glmnet(x = x, y = yi, family = "cox", measure = measure, ...)
  })
  if (length(models) > 1){
    attr(models, "survival.type") = "competing.risk"
  } else {
    attr(models, "survival.type") = "survival"
  }
  if (!is.null(LMdata)) attr(models, "LMdata") = LMdata
  if (!is.null(xcols)) attr(models, "xcols") = xcols
  class(models) = "cv.penLM"
  return(models)
}
