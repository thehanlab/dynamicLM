#' Plots the absolute w-year risk of individuals for different LM points for an event of interest within a window
#'
#' @param object Fitted LM supermodel
#' @param data Data frame of individuals from which to plot risk
#' @param format Character string specifying whether the data are in wide (default) or in long format
#' @param lm_col Character string specifying the column name in data containing the (running) time variable
#'   associated with the time-varying covariate(s); only needed if format="long"
#' @param id_col Character string specifying the column name in data containing the subject id; only needed if format="long"
#' @param w Prediction window, i.e., predict w-year (/month/..) risk from each of the tLMs.
#'   Defaults to the w used in model fitting.
#'   If w > than that used in model fitting, results are unreliable, but can be produced by setting extend=T.
#' @param cause The cause we are looking at if considering competing risks
#' @param varying Character string specifying column name in the data containing time-varying covariates; only needed if format="wide"
#' @param end_time Final time point to plot risk
#' @param extend Argument to allow for risk to be plot at landmark times that are later than the LMs used in model fitting.
#'   Default is FALSE. If set to TRUE, risks may be unreliable.
#' @param silence Silence the message when end_time > LMs used in fitting the model
#' @param pch Passed to points
#' @param lty Vector with line style
#' @param lwd Vector with line widths
#' @param col Vector with colors
#' @param main Title for the plot
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param xlim Limits for the x-axis
#' @param ylim Limits for the y-axis
#' @param x.legend,y.legend The x and y co-ordinates to be used to position the legend. They can be specified by keyword or in any way which is accepted by xy.coords.
#' @param ... Additional arguments passed to plot
#'
#' @return Single plot the absolute w-year risk of individuals
#' @details See our [GitHub](https://github.com/thehanlab/dynamicLM) for example
#'   code
#' @export
#'
plotrisk <- function(object, data, format, lm_col, id_col, w,
                     cause, varying,
                     end_time, extend=F, silence=F,
                     pch,lty,lwd,col,main,xlab,ylab,xlim,ylim,x.legend,y.legend,...){

  model_w <- object$w
  if(missing(w)){w <- model_w}
  else{
    if(w > model_w && !extend) stop(paste0("Prediction window w (=",w,") is larger than the window used in model fitting (=",model_w,").",
                                           "\nIf you wish to still make predictions at these times, set arg extend=T but note that results may be unreliable."))
    else if (w > model_w & extend) {
      if (!silence) message(paste0("NOTE: Prediction window w (=",w,") is larger than the window used in model fitting (=",model_w,"). ",
                                   "\nPredictions may be unreliable."))
    }
  }

  if(format=="long"){
    if(missing(id_col)) {
      if("ID" %in% colnames(data)) { id_col <- "ID" }
      else if("id" %in% colnames(data)) { id_col <- "id" }
      else { stop("argument 'id_col' should be specified for long format data") }
    }
    if(missing(lm_col)) {
      if("LM" %in% colnames(data)) { lm_col <- "LM" }
      else { stop("argument 'lm_col' should be specified for long format data") }
    }
    if(! id_col %in% colnames(data)) stop("arg 'id_col' is not a column in data")
    if(! lm_col %in% colnames(data)) stop("arg 'lm_col' is not a column in data")
    unique_ids <- unique(data[[id_col]])
    NF <- length(unique_ids)
  } else if (format=="wide"){
    if(missing(varying)) message("NOTE: wide format but no column with varying covariates is given.")
    NF <- nrow(data)
  } else {stop("format must be wide or long.")}

  if (missing(end_time)) end_time <- object$end_time
  else if (end_time > object$end_time & !extend){
    if(!silence) message(paste0("NOTE: arg end_time (=",end_time,") is later than the last LM used in model fitting (=",object$end_time,")",
                                "\nand has been set back to the last LM used in model fitting. (=",object$end_time,")",
                                "\nIf you wish to still plot until ",end_time, ", set arg extend=T but note that results after time ",object$end_time," may be unreliable."))
    end_time <- object$end_time
  }
  else if (end_time > object$end_time & extend){
    if(!silence) warning(paste0("NOTE: arg end_time (=",end_time,") is later than the last LM used in model fitting (=",object$end_time,")",
                                "\nResults after time ",object$end_time," may be unreliable."))
  }
  type <- object$type
  if (type == "coxph") {
    if (!missing(cause) ){ stop("No cause should be specified for a coxph model.")}
    cause <- NULL
  } else if (type == "CauseSpecificCox" | type =="CSC") {
    if (missing(cause)) { cause <- NULL }
  }

  ## Set up plot params
  if (missing(lwd))
    lwd <- rep(2, NF)
  if (missing(col))
    col <- 1:NF
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
    main <- paste0("Dynamic risk prediction (window:", object$w,")")
  if(missing(xlab))
    xlab <- "LM prediction time"
  if(missing(ylab))
    ylab <- "Risk"
  if(missing(xlim))
    xlim <- c(0, end_time)
  if (missing(x.legend)){x.legend="topright"}
  if (missing(y.legend)){y.legend=NULL}

  plotted=F
  encountered_na=F
  ids_na = c()
  ## Create plot
  if (format == "long") {
    for (i in 1:NF){
      id = unique_ids[i]
      data_ind <- data[data[[id_col]] == id,]
      x <- data_ind[[lm_col]]
      idx <- x <= end_time
      x <- x[idx]
      y <- predict.dynamicLM(object, data_ind[idx,], x, cause, extend=extend, silence=T, complete=F)$preds$risk

      ## if some entries have missing values we want to replace them by the most recent score...
      if (sum(is.na(y))!=0){
        if(sum(is.na(y))==length(y)) y = F # cannot be used
        else {
          if (is.na(y[1])){ # find first non-NA index
            firstNonNA <- min(which(!is.na(y)))
            include = firstNonNA:length(y)
            y = y[include]
            x = x[include]
            message(paste0("Individual with ID=",id," had a first entry with missing value. Replaced by 0."))
          }
          encountered_na = T
          ids_na = c(ids_na,id)
          y <- replace_na_with_last(y) # replace later NAs by the most recent score...
        }
      }

      if (inherits(y,"logical")){
        message(paste0("Individual with ID=",id," could not be plotted as they have missing values"))
      } else {
        y <- c(y[1],y)
        if (max(x) < end_time){
          x <- c(x, end_time)
          y <- c(y, y[length(y)])
        }
        if(!plotted){
          plot(stats::stepfun(x,y),
                      xlab=xlab, ylab=ylab, main=main,
                      pch=pch[i], lty=lty[i], lwd=lwd[i], col=col[i], xlim=xlim,ylim=ylim, ...)
          plotted=T
        }
        else{
          graphics::lines(stats::stepfun(x,y),pch=pch[i],lty=lty[i],lwd=lwd[i],col=col[i],...)
        }
      }
    }

    if(plotted) {
      graphics::legend(x = x.legend, y =y.legend,
           legend = unique_ids,
           lty = lty,
           col = col,
           lwd = lwd)
      if(encountered_na) message("Note that individual(s) (",paste(ids_na,collapse=", "),") had entries with missing data.The most recent previous non-missing entry was used instead.")
    }
    else{
      message("No users with non-missing data were provided. No plot could be produced.")
      return(0) # exit
    }

  } else if (format == "wide") {

    idx = data[[varying]] < end_time

    no_change = data[!idx,]
    if (nrow(no_change)>0){
      no_change$LM = 0
      no_change[[varying]] = 0
    }
    change1 = data[idx,]
    if (nrow(change1)>0){
      change1$LM = 0
      change1[[varying]] = 0
    }
    change2 = data[idx,]
    if (nrow(change2)>0){
      change2$LM = change2[[varying]]
      change2[[varying]] = 1
    }

    long_form = rbind(no_change, change1, change2)

    plotrisk(object, long_form, format="long", lm_col="LM", id_col,
                           cause, varying,
                           end_time, extend, silence,
                           pch,lty,lwd,col,main,xlab,ylab,xlim,ylim,x.legend,y.legend,...)


  }
}
