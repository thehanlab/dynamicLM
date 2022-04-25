#' Plots the absolute w-year risk of individuals for different LM points for an event of interest within a window
#'
#' @param superfm Fitted LM super model
#' @param data Data frame of individuals from which to plot risk
#' @param format Character string specifying whether the data are in wide (default) or in long format
#' @param LM_col Character string specifying the column name in data containing the (running) time variable
#' associated with the time-varying covariate(s); only needed if format="long"
#' @param id_col Character string specifying the column name in data containing the subject id; only needed if format="long"
#' @param cause The cause we are looking at, only needed if considering competing risks
#' @param varying Character string specifying column name in the data containing time-varying covariates; only needed if format="wide"
#' @param end_time Final time point to plot risk
#' @param extend Argument to allow for risk to be plot at landmark times that are later than the LMs used in model fitting.
#' Default is FALSE. If set to TRUE, risks may be unreliable.
#' @param silence Silence the message when end_time > LMs used in fitting the model
#' @param pch Passed to points
#' @param lty Vector with line style
#' @param lwd Vector with line widths
#' @param col Vector with colors
#' @param main Title for the plot
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param xlim Limits for the x-axis
#' @param ... Additional arguments passed to plot
#'
#' @return Single plot the absolute w-year risk of individuals
#' @export
#'
plotRisk <- function(superfm, data, format, LM_col, id_col,
                     cause=1, varying,
                     end_time, extend=F, silence=F,
                     pch,lty,lwd,col,main,xlab,ylab,xlim,...){
  # TODO: check that wide format works
  # NOTE: have removed a lot of arguments eg lwd
  # TODO: add w
  if(format=="long"){
    if(missing(id_col)) stop("argument 'id_col' should be specified for long format data")
    if(missing(LM_col)) stop("argument 'LM_col' should be specified for long format data")
    if(! id_col %in% colnames(data)) stop("arg 'id_col' is not a column in data")
    if(! LM_col %in% colnames(data)) stop("arg 'LM_col' is not a column in data")
    unique_ids <- unique(data[[id_col]])
    NF <- length(unique_ids)
  } else if (format=="wide"){
    if(missing(varying)) message("NOTE: wide format but no column with varying covariates is given.")
    NF <- nrow(data)
  } else {stop("format must be wide or long.")}

  if (missing(end_time)) end_time <- superfm$end_time
  else if (end_time > superfm$end_time & !extend){
    if(!silence) message(paste0("NOTE: arg end_time (=",end_time,") is later than the last LM used in model fitting (=",superfm$end_time,")",
                                "\nand has been set back to the last LM used in model fitting. (=",superfm$end_time,")",
                                "\nIf you wish to still plot until ",end_time, ", set arg extend=T but note that results after time ",superfm$end_time," may be unreliable."))
    end_time <- superfm$end_time
  }
  else if (end_time > superfm$end_time & extend){
    if(!silence) warning(paste0("NOTE: arg end_time (=",end_time,") is later than the last LM used in model fitting (=",superfm$end_time,")",
                                "\nResults after time ",superfm$end_time," may be unreliable."))
  }

  ## Set up plot params
  if (missing(lwd))
    lwd <- rep(2, NF)
  if (missing(col)) {col <- 1:NF}
  if (missing(lty))
    lty <- rep(1, NF)
  if (missing(pch))
    pch <- rep(NA_integer_, NF)
  if (length(lwd) < NF)
    lwd <- rep(lwd, NF)
  if (length(lty) < NF)
    lty <- rep(lty, NF)
  if (length(col) < NF)
    col <- rep(col, NF)
  if (length(pch) < NF)
    pch <- rep(pch, NF)
  if(missing(main))
    main <- paste0(superfm$w,"-year dynamic risk prediction")
  if(missing(xlab))
    xlab <- "LM prediction time"
  if(missing(ylab))
    ylab <- "Risk"
  if(missing(xlim))
    xlim <- c(0, end_time)

  ## Create plot
  if (format == "long") {
    for (i in 1:NF){
      id = unique_ids[i]
      data_ind <- data[data[[id_col]] == id,]
      x <- data_ind[[LM_col]]
      idx <- x <= end_time
      x <- x[idx]
      y <- predLMrisk(superfm, data_ind[idx,], x, cause, extend=extend, silence=T)$preds$risk
      if(i==1) plot(stats::stepfun(x,c(y[1],y)),
                    xlab=xlab, ylab=ylab, main=main,
                    pch=pch[i], lty=lty[i], lwd=lwd[i], col=col[i], xlim=xlim, ...)
      else graphics::lines(stats::stepfun(x,c(y[1],y)),pch=pch[i],lty=lty[i],lwd=lwd[i],col=col[i],...)
    }
    graphics::legend(x = "topright",
           legend = unique_ids,
           lty = lty,
           col = col,
           lwd = lwd)

  } else if (format == "wide") {

    ## TODO !!
    ## covs not defined... need to redo

    # if (is.null(end_time)) end_time = max(data[[outcome$time]])
    # time_col=superfm$outcome$time
    #
    # for (row in 1:NF){
    #   ind = data[row,]
    #   t1 = ind[[covs$varying]]
    #   if(t1 != ind[[time_col]] & t1 <= end_time){
    #     x <- c(0, t1)
    #     ind = cutLMsuper(ind, outcome, x, end_time, covs, format="wide")
    #   } else {
    #     x <- c(0)
    #     ind = dynpred::cutLM(ind, outcome, 0, end_time, covs, format="wide")
    #   }
    #
    #   y <- sapply(1:nrow(ind), function(row) {
    #     predLMrisk(superfm,  ind[row,], x[row], extend=extend, silence=T)
    #   })
    #
    #   if (t1!=end_time){
    #     x <- c(x, end_time)
    #     y <- c(y, y[length(y)])
    #   }
    #   if(i==1) plot(stats::stepfun(x,c(y[1],y)), xlab=xlab, ylab=ylab,main=main,pch=pch[i],lty=lty[i],lwd=lwd[i],col=col[i],...)
    #   else graphics::lines(stats::stepfun(x,c(y[1],y)),pch=pch[i],lty=lty[i],lwd=lwd[i],col=col[i],...)
    # }

  }
}
