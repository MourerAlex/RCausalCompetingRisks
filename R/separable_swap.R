#' Separable Swap Weights for IPW (Rep 1's W_D, Rep 2's W_Y)
#'
#' Swap weights are ratios of cumulative survival probabilities under
#' counterfactual arm assignments. They "swap" the observed treatment
#' history into the target counterfactual within the IPW estimator.
#'
#' For Rep 1 (W_D): ratio of D-cumprods under target a_D vs observed a_D.
#' For Rep 2 (W_Y): same construction on Y-cumprods, with the lag-trick
#' from Stensrud 2020 appendix (use lagged Y-cumprods, not the
#' contemporaneous ones).
#'
#' Cumulative survival columns are pre-computed once by
#' [separable_arm_hazards()] and consumed by the constructors below.
#'
#' @keywords internal
NULL


#' Per-Target-Arm D-Swap Weights (Rep 1, pre-truncation)
#'
#' Two columns, each meaningful only on the subset of rows where the
#' observed treatment matches the relevant arm:
#'
#' - `w_d_arm_10`: weight for arm (1, 0). Non-NA on A = 1 rows.
#' - `w_d_arm_01`: weight for arm (0, 1). Non-NA on A = 0 rows.
#'
#' Arms (1, 1) and (0, 0) use `w_arm = 1` implicitly. Raw versions
#' (`*_raw`) are preserved for [reweight()] and diagnostics.
#'
#' Requires [separable_arm_hazards()] to have attached
#' `cumprod_one_minus_hazd_a0` and `cumprod_one_minus_hazd_a1`.
#'
#' @param pt_data Person-time data frame with D cumprods.
#' @param model_d Fitted D-hazard glm, or NULL (skip).
#' @param treatment Character. Observed treatment column name.
#' @return pt_data with `w_d_arm_10_raw`, `w_d_arm_01_raw`, `w_d_arm_10`,
#'   `w_d_arm_01` added.
#' @keywords internal
swap_d_weights <- function(pt_data, model_d, treatment) {
  if (is.null(model_d)) return(pt_data)

  # Arm (1, 0): observed A = 1, target a_D = 0
  pt_data$w_d_arm_10_raw <- ifelse(
    pt_data[[treatment]] == 1,
    pt_data$cumprod_one_minus_hazd_a0 / pt_data$cumprod_one_minus_hazd_a1,
    NA_real_
  )
  # Arm (0, 1): observed A = 0, target a_D = 1
  pt_data$w_d_arm_01_raw <- ifelse(
    pt_data[[treatment]] == 0,
    pt_data$cumprod_one_minus_hazd_a1 / pt_data$cumprod_one_minus_hazd_a0,
    NA_real_
  )
  pt_data$w_d_arm_10 <- pt_data$w_d_arm_10_raw
  pt_data$w_d_arm_01 <- pt_data$w_d_arm_01_raw
  pt_data
}


#' Per-Target-Arm Y-Swap Weights (Rep 2, pre-truncation)
#'
#' Two columns, each meaningful only on the rows Rep 2 stands in
#' (i.e., A = a_d):
#'
#' - `w_y_arm_10`: weight for arm (1, 0). Rep 2 stands in a_d = 0;
#'   non-NA on A = 0 rows.
#' - `w_y_arm_01`: weight for arm (0, 1). Rep 2 stands in a_d = 1;
#'   non-NA on A = 1 rows.
#'
#' Arms (1, 1) and (0, 0) use `w_arm = 1` implicitly. Raw versions
#' (`*_raw`) are preserved for [reweight()] and diagnostics.
#'
#' Two-component construction (Stensrud 2020 appendix):
#' \deqn{W_Y(s) = \frac{h_Y^{a_y}(s)}{h_Y^{a_d}(s)} \cdot
#'                \prod_{j < s} \frac{1 - h_Y^{a_y}(j)}{1 - h_Y^{a_d}(j)}}
#'
#' Requires [separable_arm_hazards()] to have attached `haz_y_a0`,
#' `haz_y_a1`, `lag_cumprod_one_minus_hazy_a0`, and
#' `lag_cumprod_one_minus_hazy_a1`. A small floor avoids 0/0 in extreme
#' cases; positivity violations propagate to `w_total` for upstream
#' diagnostics rather than being silently hidden.
#'
#' @param pt_data Person-time data frame with Y hazard + lagged cumprod.
#' @param model_y Fitted Y-hazard glm, or NULL (skip).
#' @param treatment Character. Observed treatment column name.
#' @return pt_data with `w_y_arm_10_raw`, `w_y_arm_01_raw`, `w_y_arm_10`,
#'   `w_y_arm_01` added.
#' @keywords internal
swap_y_weights <- function(pt_data, model_y, treatment) {
  if (is.null(model_y)) return(pt_data)

  eps <- sqrt(.Machine$double.eps)

  # Arm (1, 0): a_y = 1, a_d = 0; Rep 2 stands in a_d = 0 -> A = 0 rows
  pt_data$w_y_arm_10_raw <- ifelse(
    pt_data[[treatment]] == 0,
    (pt_data$haz_y_a1 / pmax(pt_data$haz_y_a0, eps)) *
      (pt_data$lag_cumprod_one_minus_hazy_a1 /
         pmax(pt_data$lag_cumprod_one_minus_hazy_a0, eps)),
    NA_real_
  )

  # Arm (0, 1): a_y = 0, a_d = 1; Rep 2 stands in a_d = 1 -> A = 1 rows
  pt_data$w_y_arm_01_raw <- ifelse(
    pt_data[[treatment]] == 1,
    (pt_data$haz_y_a0 / pmax(pt_data$haz_y_a1, eps)) *
      (pt_data$lag_cumprod_one_minus_hazy_a0 /
         pmax(pt_data$lag_cumprod_one_minus_hazy_a1, eps)),
    NA_real_
  )

  pt_data$w_y_arm_10 <- pt_data$w_y_arm_10_raw
  pt_data$w_y_arm_01 <- pt_data$w_y_arm_01_raw
  pt_data
}
