#' Time-to-event data of cancer relapse
#'
#' Simple synthetic dataset containing the time-to-event of cancer relapse
#' (event=1) with the competing risk in long-form with patient information.
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
#'   \item{treatment}{Patient's treatment status, treatment = 1 = on treatment,
#'                    treament = 0 = patient is off treatment}
#'   \item{T_txgiven}{Follow-up time, i.e., time at which updated treatment (tx)
#'                    information was provided, which is equivalent to the time
#'                    point at which the patient entry was created.}
#' }
"relapse"


#' Time-to-event data of SPLC
#'
#' Synthetic dataset containing the time-to-event of secondary primary lung
#' cancer (SPLC) with competing risks of lung cancer death (cause 2) and
#' other-cause death (cause 3) in long-form with patient information.
#'
#' @format A data frame with 875 rows and 23 columns:
#' \describe{
#'   \item{ID}{Patient ID}
#'   \item{event}{Event of interest (0=censoring, 1=relapse, 2,3=competing risks)}
#'   \item{Time}{Time-to-event}
#'   \item{T.fup}{Follow-up time, i.e., time at which updated covariate
#'                information was provided. This is equivalent to the time point
#'                at which the patient entry was created.}
#'   \item{age.ix}{Patient's age at time of diagnosis}
#'   \item{male}{Sex of patient, 1 = male, 0 = female}
#'   \item{fh}{Family history}
#'   \item{ph}{Prior history}
#'   \item{bmi}{Patient's body mass index at diagnosis}
#'   \item{stage.ix}{Cancer stage at diagnosis (advanced/not)}
#'   \item{surgery.ix}{Surgery (yes/no)}
#'   \item{radiation.ix}{Radiation (yes/no)}
#'   \item{chemo.ix}{Chemotherapy (yes/no)}
#'   \item{smkstatus}{Smoking status. Former = 2, Current = 3}
#'   \item{cigday}{Cigarettes per day.}
#'   \item{packyears}{Number of pack years}
#'   \item{quityears}{Number of quit years}
#'   \item{hist_*}{Histology at diagnosis}
#' }
"splc"

#' Time-to-event data of SPLC (test set)
#'
#' Synthetic dataset containing the time-to-event of secondary primary lung
#' cancer (SPLC) with competing risks of lung cancer death (cause 2) and
#' other-cause death (cause 3) in long-form with patient information.
#'
#' @format A data frame with 607 rows and 24 columns:
#' \describe{
#'   \item{ID}{Patient ID}
#'   \item{event}{Event of interest (0=censoring, 1=relapse, 2,3=competing risks)}
#'   \item{Time}{Time-to-event}
#'   \item{T.fup}{Follow-up time, i.e., time at which updated covariate
#'                information was provided. This is equivalent to the time point
#'                at which the patient entry was created.}
#'   \item{age.ix}{Patient's age at time of diagnosis}
#'   \item{male}{Sex of patient, 1 = male, 0 = female}
#'   \item{fh}{Family history}
#'   \item{ph}{Prior history}
#'   \item{bmi}{Patient's body mass index at diagnosis}
#'   \item{stage.ix}{Cancer stage at diagnosis (advanced/not)}
#'   \item{surgery.ix}{Surgery (yes/no)}
#'   \item{radiation.ix}{Radiation (yes/no)}
#'   \item{chemo.ix}{Chemotherapy (yes/no)}
#'   \item{smkstatus}{Smoking status. Former = 2, Current = 3}
#'   \item{cigday}{Cigarettes per day.}
#'   \item{packyears}{Number of pack years}
#'   \item{quityears}{Number of quit years}
#'   \item{hist_*}{Histology at diagnosis}
#' }
"splc_test"
