#' ADHD Synthetic Dataset
#'
#' Dataset created for CATIE 2023
#'
#' @format A data frame with 150 rows and 13 columns
#' \describe{
#'    \item{ID}{numeric id}
#'    \item{odd}{ODD diagnosis, reflecting whether the child was (coded as 1) or was not (coded as 0)
#'    diagnosed with ODD (oppositional defiant disorder) before the first-stage intervention.}
#'    \item{severity}{ADHD score, reflecting ADHD symptoms at
#'    the end of the previous school year (larger values reflect greater symptoms). Centered and scaled to standard normal.}
#'    \item{priormed}{medication prior to first-stage intervention, reflecting whether the child did
#'    (coded as 1) or did not (coded as 0) receive medication during the previous school year.}
#'    \item{race}{race variable: white (coded 1) versus non- white (coded 0).}
#'    \item{Y0}{baseline school performance (higher values reflect better performance).}
#'    \item{A1}{stage 1 treatment assignment. A1 = -1 if the first stage intervention was Medication;
#'    A2 = 1 if the first stage intervention was Behavioral Intervention.}
#'    \item{R}{response status. R = 0 if participant did not respond to the first stage intervention; R= 1 if he or she responded.}
#'    \item{event_time}{number of months until non-response. range 1-8.}
#'    \item{adherence}{adherence to the stage 1 intervention: whether the child did (coded as 1) or did not (coded as 0) show high adherence to initial treatment.}
#'    \item{Y1}{school performance at the end of stage 1.}
#'    \item{A2}{stage 2 treatment assignment. Non-responders (R=0) to the first stage intervention entered the second stage of intervention and were re-randomized to either ADDING the other type of intervention (coded as -1) or INTENSIFYING the first stage intervention (coded as 1). For responders A2 = NA.}
#'    \item{Y2}{end-of study school performance}
#' }
#' @inheritSection simADHDsmart Operating Characteristics
#' @source Pelham, et. al.
"adhd"

#' Autism Synthetic Dataset
#'
#' Dataset created for CATIE 2023
#'
#' @format A data frame with 200 rows and 11 columns
#'
#' #' \describe{
#'    \item{ID}{numeric id}
#' }
#'
#' @source NULL
"autism"
