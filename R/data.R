#' Time-to-event data of cancer relapse
#'
#' Synthetic dataset containing the time-to-event of cancer relapse (event=1) with the competing risk death (event=2) in long-form with patient information.
#'
#' @format A data frame with 989 rows and 9 columns:
#' \describe{
#'   \item{ID}{Patient ID}
#'   \item{Time}{Time-to-event}
#'   \item{event}{Event of interest (0=censoring, 1=relapse, 2,3=competing risks)}
#'   \item{age.at.time.0}{Patient's age at time of diagnosis}
#'   \item{male}{Sex of patient, 1=male, 0=female}
#'   \item{stage}{Cancer stage at diagnosis}
#'   \item{bmi}{Patient's body mass index at diagnosis}
#'   \item{treatment}{Patient's treatment status, treatment=1=on treatment, treament=0=patient is off treatment}
#'   \item{T_txgiven}{Follow-up time, i.e., time at which updated treatment (tx) information was provided, which is equivalent to the time point at which the patient entry was created.}
#' }
"relapse"
