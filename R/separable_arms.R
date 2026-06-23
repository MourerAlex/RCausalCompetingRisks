#' Separable Arm Specification
#'
#' Arm dispatch (a_Y, a_D) and the component columns (A_Y, A_D) used by
#' both the g-formula and IPW separable estimators.
#'
#' @keywords internal
NULL


#' Predict Hazards Under Each Separable Counterfactual Arm
#'
#' Attaches per-arm hazard predictions and cumulative-survival columns to
#' `pt_data`, plus Y-only lagged cumulative-survival columns used by Rep 2's
#' W_Y construction (Stensrud 2020 appendix). The treatment-decomposition
#' columns `A_y` and `A_d` are temporarily overwritten via
#' `predict_hazard_under()`; cumulative survival is computed via
#' `cumprod_survival()`.
#'
#' For each non-NULL model, four columns are added:
#' `haz_<event>_a1`, `haz_<event>_a0`, `cumprod_one_minus_haz<event>_a1`,
#' `cumprod_one_minus_haz<event>_a0`. For Y, two additional lagged columns
#' (`lag_cumprod_one_minus_hazy_a1`, `lag_cumprod_one_minus_hazy_a0`) are
#' added because Rep 2's W_Y uses event-free survival up to (not including)
#' the current interval.
#'
#' @param pt_data Person-time data frame.
#' @param model_y Fitted Y-hazard glm, or NULL.
#' @param model_d Fitted D-hazard glm, or NULL.
#' @param id_col Character. Subject id column name.
#' @return `pt_data` enriched with the columns above.
#' @keywords internal
separable_arm_hazards <- function(pt_data, model_y, model_d, id_col) {

  # D-hazards under each counterfactual A_d (Rep 1 swap weights need both)
  if (!is.null(model_d)) {
    pt_data$haz_d_a1 <- predict_counterfactual_hazard(
      model_d, pt_data, "A_d", 1, "D-hazard")
    pt_data$haz_d_a0 <- predict_counterfactual_hazard(
      model_d, pt_data, "A_d", 0, "D-hazard")
    pt_data$cumprod_one_minus_hazd_a1 <- cumprod_survival(
      pt_data$haz_d_a1, pt_data[[id_col]]
    )
    pt_data$cumprod_one_minus_hazd_a0 <- cumprod_survival(
      pt_data$haz_d_a0, pt_data[[id_col]]
    )
  }

  # Y-hazards under each counterfactual A_y (Rep 2 swap weights need both)
  if (!is.null(model_y)) {
    pt_data$haz_y_a1 <- predict_counterfactual_hazard(
      model_y, pt_data, "A_y", 1, "Y-hazard")
    pt_data$haz_y_a0 <- predict_counterfactual_hazard(
      model_y, pt_data, "A_y", 0, "Y-hazard")
    pt_data$cumprod_one_minus_hazy_a1 <- cumprod_survival(
      pt_data$haz_y_a1, pt_data[[id_col]]
    )
    pt_data$cumprod_one_minus_hazy_a0 <- cumprod_survival(
      pt_data$haz_y_a0, pt_data[[id_col]]
    )

    # Y-only lagged cumprod for Rep 2's W_Y lag-trick
    pt_data$lag_cumprod_one_minus_hazy_a1 <- ave(
      pt_data$cumprod_one_minus_hazy_a1, pt_data[[id_col]],
      FUN = function(v) c(1, v[-length(v)])
    )
    pt_data$lag_cumprod_one_minus_hazy_a0 <- ave(
      pt_data$cumprod_one_minus_hazy_a0, pt_data[[id_col]],
      FUN = function(v) c(1, v[-length(v)])
    )
  }

  pt_data
}
