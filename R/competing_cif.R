#' Discrete-Time Cumulative Incidence from Cause-Specific Hazards
#'
#' Shared recursion for the primary-event cumulative incidence \eqn{F_Y} over a
#' discrete grid, from per-interval cause-specific hazards for the primary event
#' \eqn{Y} and the competing event \eqn{D}. This is NOT the Aalen-Johansen
#' estimator: it uses two separate binary cause-specific hazards under the
#' sequential within-interval ordering \eqn{C \to D \to Y} (D screens Y), the
#' discrete-time pooled-hazard convention. Weight-agnostic — hazards may come
#' from a fitted model (g-formula) or a Hajek-weighted pseudo-population (IPW).
#'
#' \deqn{F_Y(k) = \sum_{j \le k} h_Y(j)\,(1 - h_D(j))\,S(j-1), \quad
#'       S(k) = \prod_{j \le k} (1 - h_Y(j))\,(1 - h_D(j))}
#'
#' @param haz_y,haz_d Numeric vectors. Cause-specific hazards by interval for
#'   one unit (a subject, or a marginal arm), ordered \eqn{k = 1, \dots, K}.
#' @return Numeric vector. Cumulative incidence at each interval.
#' @keywords internal
cumulative_incidence <- function(haz_y, haz_d) {
  event_free <- (1 - haz_y) * (1 - haz_d)
  surv <- c(1, cumprod(event_free)[-length(event_free)])
  cumsum(haz_y * (1 - haz_d) * surv)
}


#' Cumulative Incidence from Weighted Person-Time Data
#'
#' Discrete-time cumulative incidence of `Y` in the presence of competing event
#' `D`, computed from per-row event indicators and Hajek-weighted hazards. The
#' caller is responsible for restricting inputs to a single standing-in arm and
#' for building the per-row `weights` (IPTW * IPCW * separable swap); it
#' is weight-agnostic — it only applies whatever weights it is given.
#'
#' Forms weighted hazards by interval via [weighted_hazard_by_k()], aligns them
#' to `cut_times`, and delegates the discrete-time recursion to
#' [cumulative_incidence()].
#'
#' @param y_event,d_event,c_event Numeric vectors of event indicators
#'   (Y, D, censoring), no-NA `{0, 1}` schema.
#' @param k Vector of interval indices.
#' @param weights Numeric vector of per-row weights.
#' @param cut_times Numeric vector of time points to evaluate at.
#' @return Numeric vector of cumulative incidence at each `cut_times[k]`.
#' @keywords internal
cum_inc_from_weighted <- function(y_event, d_event, c_event, k, weights,
                                  cut_times) {
  # At-risk masks (ordering C -> D -> Y): censored rows leave both risk
  # sets; same-interval D exits additionally leave the Y risk set. Subset
  # each hazard's inputs explicitly — mirrors the g-formula fit path
  # (pt_data[y_rows, ]); a boolean "at risk", not NA "missing".
  y_at_risk <- c_event == 0 & d_event == 0
  d_at_risk <- c_event == 0

  haz_y_by_k <- weighted_hazard_by_k(y_event[y_at_risk], k[y_at_risk],
                                     weights[y_at_risk])
  haz_d_by_k <- weighted_hazard_by_k(d_event[d_at_risk], k[d_at_risk],
                                     weights[d_at_risk])

  # Align to integer interval index 1..K_max (missing k -> 0 hazard)
  key <- as.character(seq_along(cut_times))
  haz_y_vec <- unname(haz_y_by_k[key])
  haz_d_vec <- unname(haz_d_by_k[key])
  haz_y_vec[is.na(haz_y_vec)] <- 0
  haz_d_vec[is.na(haz_d_vec)] <- 0

  cumulative_incidence(haz_y_vec, haz_d_vec)
}
