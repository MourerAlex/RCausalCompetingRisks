#' Compute Risk Table Counts by Observed Treatment Arm
#'
#' Returns counts at each `cut_times` value, grouped by observed treatment
#' group (e.g., A = 0 and A = 1). Complements cumulative incidence curves
#' for reporting (adjustedCurves-style risk tables).
#'
#' The user picks **one** count type per call. Counts are computed on
#' observed subjects — NOT on counterfactual arms (`arm_10`, `arm_01`
#' would not correspond to any real subjects).
#'
#' @param fit A `"causal_competing_risks_fit"` object from [causal_competing_risks()].
#' @param count Character. One of:
#'   \describe{
#'     \item{`"at_risk"`}{Number of subjects with a row at k (still under
#'       observation at the start of interval k).}
#'     \item{`"events_y"`}{Number of primary (Y) events occurring in
#'       interval k.}
#'     \item{`"events_d"`}{Number of competing (D) events occurring in
#'       interval k.}
#'     \item{`"censored"`}{Number of subjects censored in interval k.}
#'   }
#'
#' @return A data.frame with:
#'   \describe{
#'     \item{k}{Time point (from `fit$times`).}
#'     \item{`<treatment>_<value>`}{One column per observed treatment level,
#'       named e.g. `A_0` and `A_1`.}
#'   }
#'
#' @examples
#' \dontrun{
#' fit <- causal_competing_risks(pt)
#' risk_table(fit, count = "at_risk")
#' risk_table(fit, count = "events_y")
#' }
#'
#' @seealso [causal_competing_risks()], [plot.causal_competing_risks_risk()]
#' @family accessors
#' @export
risk_table <- function(fit,
                       count = c("at_risk", "events_y", "events_d",
                                 "censored")) {
  stopifnot(inherits(fit, "causal_competing_risks_fit"))
  count <- match.arg(count)
  risk_table_internal(
    pt_data   = fit$person_time,
    id_col    = fit$id_col,
    trt_col   = fit$treatment_col,
    cut_times = fit$times,
    count     = count
  )
}


#' Internal Risk Table Worker (No Class Check)
#'
#' Used by both [risk_table()] and [plot.causal_competing_risks_risk()]. The split lets
#' the plot method compute the table from a `causal_competing_risks_risk` object's stored
#' references without needing a back-reference to the full fit.
#'
#' @param pt_data Person-time data frame.
#' @param id_col,trt_col Character column names.
#' @param cut_times Numeric vector of interval starts.
#' @param count Character; one of `c("at_risk", "events_y", "events_d", "censored")`.
#' @return A data.frame with `k` and one column per observed treatment value.
#' @family internal
#' @keywords internal
risk_table_internal <- function(pt_data, id_col, trt_col, cut_times, count) {
  trt_vals <- sort(unique(pt_data[[trt_col]]))

  result <- data.frame(k = cut_times)              # report clock time
  for (a in trt_vals) {
    col_name <- paste0(trt_col, "_", a)
    vals <- vapply(seq_along(cut_times), function(k_idx) {   # filter by integer index
      rows_a_k <- pt_data[pt_data[[trt_col]] == a & pt_data$k == k_idx, ]
      if (count == "at_risk") {
        length(unique(rows_a_k[[id_col]]))
      } else if (count == "events_y") {
        as.integer(sum(rows_a_k$y_event, na.rm = TRUE))
      } else if (count == "events_d") {
        as.integer(sum(rows_a_k$d_event, na.rm = TRUE))
      } else {
        as.integer(sum(rows_a_k$c_event, na.rm = TRUE))
      }
    }, integer(1))
    result[[col_name]] <- vals
  }
  result
}
