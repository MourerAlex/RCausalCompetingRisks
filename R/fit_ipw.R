#' Separable IPW Estimation and Per-Arm CIF
#'
#' Stensrud 2020 has two algebraically equivalent IPW representations
#' under full isolation:
#'
#' - **Rep 1** stands in `A = a_y` and applies a single D-cumulative-survival
#'   swap weight (W_D from [swap_d_weights()]).
#' - **Rep 2** stands in `A = a_d` and applies a hazard-ratio × lagged-
#'   cumulative-survival swap weight (W_Y from [swap_y_weights()]).
#'
#' Swap-weight construction is rep-specific and lives upstream. The CIF
#' computation given those weights is the shared Hajek estimator
#' [cum_inc_from_weighted()]. Both reps run by default for sensitivity
#' comparison.
#'
#' This file contains the full IPW pipeline: orchestrator
#' [fit_ipw()] (builds all source weights and dispatches), the
#' shared dispatcher [estimate_weighted_cum_inc()], and the per-arm
#' Hajek estimator [weighted_arm_cum_inc()].
#'
#' @keywords internal
NULL


#' IPW Estimation of Separable Effects (Orchestrator)
#'
#' Builds the per-row source weight columns (`w_cens`, `w_a`,
#' `w_d_arm_*`, `w_y_arm_*`) on `pt_data`, applies symmetric
#' truncation (Cole & Hernán 2008), then dispatches to Rep 1 and/or
#' Rep 2 via [estimate_weighted_cum_inc()].
#'
#' @param pt_data Data frame in person-time format.
#' @param models List of fitted glm objects (Y, D, C hazards plus A
#'   propensity, optionally A-numerator for stabilization).
#' @param treatment Character. Treatment column name.
#' @param ipw_reps Character vector. Subset of `c("ipw_rep1", "ipw_rep2")`
#'   indicating which representation(s) to compute.
#' @param cut_times Numeric vector. Time cut points.
#' @param id_col Character. Subject identifier column name.
#' @param truncate Either NULL (no truncation) or a length-2 numeric
#'   vector of percentile bounds. Default `c(0.01, 0.99)`.
#'
#' @return A list with `cumulative_incidence_rep1`,
#'   `cumulative_incidence_rep2`, `pt_data_weighted`, `weight_summary`,
#'   `flagged_ids`, `flagged_log`.
#'
#' @keywords internal
fit_ipw <- function(pt_data,
                         models,
                         treatment,
                         ipw_reps,
                         cut_times,
                         id_col,
                         truncate = c(0.01, 0.99)) {

  # 1. Enrich pt_data with predicted hazards + cumprods + lag (Y only)
  pt_data <- separable_arm_hazards(
    pt_data, models$model_y, models$model_d, id_col
  )

  # 2. Censoring weights (IPCW via the generic wrapper).
  if (!is.null(models$model_c)) {
    pt_data$haz_c <- predict(
      models$model_c, newdata = pt_data, type = "response"
    )
    pt_data$cumprod_one_minus_hazc <- cumprod_survival(
      pt_data$haz_c, pt_data[[id_col]]
    )
    pt_data$w_cens_raw <- ipw_cens(models$model_c, pt_data, id_col)
    pt_data$w_cens     <- pt_data$w_cens_raw
  }

  # 3. Treatment weight (W_A): 1/pi(A|L), or stabilized ratio if numerator
  #    model present. Stored as a column so apply_weight_truncation() and
  #    summarize_weights() can include it in diagnostics.
  if (!is.null(models$model_a)) {
    pt_data$w_a_raw <- ipw_static_trt(
      models$model_a, pt_data, treatment, id_col,
      model_num = models$model_a_num
    )
    pt_data$w_a <- pt_data$w_a_raw
  }

  # 4. Per-target-arm swap weights:
  #    D for Rep 1 (w_d_arm_10/01) and Y for Rep 2 (w_y_arm_10/01).
  pt_data <- swap_d_weights(pt_data, models$model_d, treatment)
  pt_data <- swap_y_weights(pt_data, models$model_y, treatment)

  # 5. Apply symmetric weight truncation to all source weight columns.
  ccr_w_cols <- c("w_cens", "w_a", "w_d_arm_10", "w_d_arm_01", "w_y_arm_10", "w_y_arm_01")
  trunc_result <- apply_weight_truncation(pt_data, id_col, truncate = truncate, w_cols = ccr_w_cols)
  pt_data     <- trunc_result$pt_data
  flagged_ids <- trunc_result$flagged_ids
  flagged_log <- trunc_result$flagged_log

  # 6. Dispatch per representation
  cum_inc_rep1 <- NULL
  cum_inc_rep2 <- NULL
  if ("ipw_rep1" %in% ipw_reps) {
    cum_inc_rep1 <- estimate_weighted_cum_inc(
      pt_data, treatment, cut_times, rep = 1, id_col = id_col
    )
  }
  if ("ipw_rep2" %in% ipw_reps) {
    cum_inc_rep2 <- estimate_weighted_cum_inc(
      pt_data, treatment, cut_times, rep = 2, id_col = id_col
    )
  }

  list(
    cumulative_incidence_rep1 = cum_inc_rep1,
    cumulative_incidence_rep2 = cum_inc_rep2,
    pt_data_weighted          = pt_data,
    weight_summary            = summarize_weights(pt_data, w_base = c("w_cens", "w_a", "w_d_arm_10", "w_d_arm_01", "w_y_arm_10", "w_y_arm_01")),
    flagged_ids               = flagged_ids,
    flagged_log               = flagged_log
  )
}


#' Per-Arm Cumulative Incidence (Rep 1 or Rep 2)
#'
#' Subsets to the standing-in arm, reads the rep-specific swap weight,
#' the censoring weight (IPCW), and the propensity weight (W_A) from
#' pre-built columns on `pt_data`, then computes weighted cumulative
#' incidence. All source weights (`w_cens`, `w_a`, swap weights) are
#' constructed and truncated upstream by [fit_ipw()].
#'
#' @param pt_data Person-time data frame with all weight columns
#'   pre-built by [fit_ipw()].
#' @param treatment Character. Treatment column name.
#' @param a_y,a_d Numeric (0 or 1). Arm configuration.
#' @param cut_times Numeric vector. Time points to evaluate at.
#' @param rep Integer. 1 for Rep 1, 2 for Rep 2.
#' @param id_col Character. Subject id column name (kept for symmetry
#'   with the dispatcher; not used here directly).
#' @return Numeric vector of cumulative incidence at each `cut_times[k]`.
#' @keywords internal
weighted_arm_cum_inc <- function(pt_data, treatment, a_y, a_d, cut_times,
                                  rep, id_col) {

  # Pick standing-in arm and swap-weight columns by rep
  if (rep == 1) {
    stand_in     <- a_y
    w_arm_10_col <- "w_d_arm_10"
    w_arm_01_col <- "w_d_arm_01"
  } else {
    stand_in     <- a_d
    w_arm_10_col <- "w_y_arm_10"
    w_arm_01_col <- "w_y_arm_01"
  }

  # Stand in the chosen arm
  d <- pt_data[pt_data[[treatment]] == stand_in, ]

  # Arm-specific swap weight (1 on diagonal; pre-built columns off-diagonal)
  if (a_y == a_d) {
    d$w_arm <- 1
  } else if (a_y == 1 && a_d == 0) {
    d$w_arm <- d[[w_arm_10_col]]
  } else {
    d$w_arm <- d[[w_arm_01_col]]
  }

  # All source weights are pre-built columns; scalar 1 fallback when absent.
  w_cens <- if ("w_cens" %in% names(d)) d$w_cens else 1
  w_a    <- if ("w_a"    %in% names(d)) d$w_a    else 1

  d$w_total <- w_cens * d$w_arm * w_a

  cum_inc_from_weighted(
    y_event   = d$y_event,
    d_event   = d$d_event,
    c_event   = d$c_event,
    k         = d$k,
    weights   = d$w_total,
    cut_times = cut_times
  )
}


#' Per-Arm CIF Dispatcher (Rep 1 or Rep 2)
#'
#' Iterates the arms defined by [arm_spec()] and returns a wide data
#' frame with one column per arm. `arm_01` enables Decomposition B
#' sensitivity alongside the default Decomposition A.
#'
#' @param pt_data,treatment,cut_times,id_col See [weighted_arm_cum_inc()].
#' @param rep Integer. 1 for Rep 1, 2 for Rep 2.
#' @return Data frame with columns `k` and one column per arm name in
#'   [arm_spec()]`$name`.
#' @keywords internal
estimate_weighted_cum_inc <- function(pt_data, treatment, cut_times, rep,
                                       id_col) {
  spec <- arm_spec()
  cum_inc_list <- lapply(seq_len(nrow(spec)), function(i) {
    weighted_arm_cum_inc(
      pt_data, treatment,
      a_y = spec$a_y[i], a_d = spec$a_d[i],
      cut_times = cut_times, rep = rep,
      id_col = id_col
    )
  })
  names(cum_inc_list) <- spec$name
  do.call(data.frame, c(list(k = cut_times), cum_inc_list))
}
