# TODO: add refs
# TODO: add roxygen
get_lm_data <- function (data, outcome, LM, horizon, covs,
                         format = c("wide", "long"), id, rtime, right = TRUE) {
  format <- match.arg(format)
  if (format == "wide") {
    LMdata <- data
    if (!is.null(covs$varying)){
      for (col in covs$covarying)
        LMdata[[col]] <- 1 - as.numeric(LMdata[[col]] > LM)
    }

  }
  else {
    if (missing(id))
      stop("argument 'id' should be specified for long format data")
    if (missing(rtime))
      stop("argument 'rtime' should be specified for long format data")
    ord <- order(data[[id]], data[[rtime]])
    data <- data[ord, ]
    ids <- unique(data[[id]])
    n <- length(ids)
    LMdata <- data[which(!duplicated(data[[id]])), ]
    for (i in 1:n) {
      wh <- which(data[[id]] == ids[i])
      di <- data[wh, ]
      idx <- cut(LM, c(di[[rtime]], Inf), right = right, labels = FALSE)
      if (!is.na(idx)) {
        LMdata[i, ] <- di[idx, ]
      } else {
        LMdata[i, ] <- di[1, ]
        LMdata[i, covs$varying] <- NA
        LMdata[i, rtime] <- NA
      }
    }
  }
  LMdata <- LMdata[LMdata[[outcome$time]] > LM, ]
  if (format == "long")
    LMdata <- LMdata[!is.na(LMdata[[id]]), ]
  LMdata[outcome$status] <- LMdata[[outcome$status]] *
                              as.numeric(LMdata[[outcome$time]] <= horizon)
  LMdata[outcome$time] <- pmin(as.vector(LMdata[[outcome$time]]), horizon)
  LMdata$LM <- LM
  if (format == "long")
    cols <- match(c(id, outcome$time, outcome$status, covs$fixed,
                    covs$varying, rtime, "LM"), names(LMdata))
  else cols <- match(c(outcome$time, outcome$status, covs$fixed,
                       covs$varying, "LM"), names(LMdata))
  return(LMdata[, cols])
}
