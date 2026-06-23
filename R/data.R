#' Prostate Cancer Clinical Trial Data
#'
#' The Byar & Green (1980) prostate cancer dataset,
#' originally from the Veterans Administration Cooperative Urological Research
#' Group (VACURG) randomized trial. Contains 502 subjects randomized to
#' DES (diethylstilbestrol) or placebo.
#'
#' @format A data frame with 502 rows and the following columns:
#' \describe{
#'   \item{id}{Subject identifier.}
#'   \item{A}{Treatment indicator: 1 = DES, 0 = placebo.}
#'   \item{event_time}{Time to event or censoring (months).}
#'   \item{event_type}{Event type: 0 = censored, 1 = prostate cancer death (Y),
#'     2 = other-cause death (D).}
#'   \item{normal_act}{Normal activity: 1 = yes, 0 = no.}
#'   \item{age_cat}{Age category.}
#'   \item{cv_hist}{Cardiovascular disease history.}
#'   \item{hemo_bin}{Hemoglobin level (binarized).}
#' }
#'
#' @source Byar DP, Green SB (1980). "The Choice of Treatment for Cancer
#'   Patients Based on Covariate Information." *Bulletin du Cancer*, 67,
#'   477-490. Data available via `Hmisc::getHdata("prostate")`.
#'
#' @references
#' Rojas-Saunero LP, Young JG, Didelez V, Ikram MA, Swanson SA (2022).
#' "Considering questions before methods in dementia research with competing
#' events and causal goals." *American Journal of Epidemiology*.
#'
#' @examples
#' \dontrun{
#' data(prostate_data)
#' str(prostate_data)
#' table(prostate_data$event_type, prostate_data$A)
#' }
"prostate_data"
