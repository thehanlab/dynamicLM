#' Get the PBC data in long data from the survival package
#'
#' This function prepares the PBC data from the `survival` package in long
#' format, merging baseline data with longitudinal data and cleaning it for
#' analysis.
#'
#' @return A data frame with the PBC data in long format
#' @export
#' @import survival
#'
#' @examples
#' pbc_df <- get_pbc_long()
#' head(pbc_df)
#'
get_pbc_long <- function() {
  pbc <- survival::pbc
  pbcseq <- survival::pbcseq

  # only the first 312 patients are in both datasets
  pbc1 <- subset(pbc, id <= 312, select = c(id:sex, stage))
  # merge baseline data with longitudinal data
  pbc_df <- survival::tmerge(pbc1, pbc1, id=id, endpt = event(time, status))
  pbc_df <- survival::tmerge(pbc_df, pbcseq, id=id,
                             # make sure time-dependent covariates (tdc) vary
                             albumin = tdc(day, albumin),
                             alk.phos = tdc(day, alk.phos),
                             ascites = tdc(day, ascites),
                             ast = tdc(day, ast),
                             bili = tdc(day, bili),
                             chol = tdc(day, chol),
                             edema = tdc(day, edema),
                             hepato = tdc(day, hepato),
                             platelet = tdc(day, platelet),
                             protime = tdc(day, protime),
                             spiders = tdc(day, spiders))

  # use complete data
  incomplete_ids <- unique(pbc_df$id[!stats::complete.cases(pbc_df)])
  pbc_df <- pbc_df[!pbc_df$id %in% incomplete_ids, ]

  # convert times to years for easier reading later
  pbc_df$time <- round(pbc_df$time / 365.25, 1)
  pbc_df$tstart <- round(pbc_df$tstart / 365.25, 1)
  pbc_df$tstop <- round(pbc_df$tstop / 365.25, 1)

  # convert factor variables to numeric
  pbc_df$male <- ifelse(pbc_df$sex == "m", 1, 0); pbc_df$sex <- NULL

  return(pbc_df)
}
