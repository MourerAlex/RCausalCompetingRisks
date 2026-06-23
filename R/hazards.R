#' Fit Hazard Models
#'
#' Fits pooled logistic regression models for the Y-hazard, D-hazard, and
#' (optionally) censoring hazard on person-time data. For each fitted model,
#' collects diagnostics (convergence, min/max fitted probability, positivity
#' flag, captured `glm()` warnings).
#'
#' @param pt_data Data frame in person-time format.
#' @param treatment Character. Treatment column name (used for the censoring
#'   model; Y and D models use `A_y` / `A_d` working copies).
#' @param covariates Character vector. Covariate column names.
#' @param active_methods Character vector. Subset of
#'   `c("gformula", "ipw_rep1", "ipw_rep2")` indicating which methods will run.
#'   Determines which models get fit.
#' @param formulas Named list or NULL. User-specified formulas (names `y`,
#'   `d`, `c`). Any entry absent falls back to the default formula.
#' @param ipcw Logical. When FALSE, the censoring model is not
#'   fit and `model_c` stays NULL.
#'
#' @return A named list with two entries:
#'   \describe{
#'     \item{models}{Named list: `model_y`, `model_d`, `model_c` (glm
#'       objects or NULL).}
#'     \item{checks}{Named list: `y`, `d`, `c` (per-model diagnostics or
#'       NULL). See [check_fitted_positivity()] for the per-model structure.}
#'   }
#'
#' @details
#' ## Default formula
#' `y_event/d_event/c_event ~ A_y/A_d/treatment + k + I(k^2) + I(k^3) + covariates` (additive,
#' no interaction). The time trend is a cubic polynomial in the **integer
#' interval index** `k` (`1..K_max`), inherited from CausalSurvival's grid.
#' Because `k` is an index, the polynomial measures time on the interval-count
#' scale, not the clock scale. With equi-spaced cut points the two are an
#' affine reparameterization (same fit); with **unequally-spaced** cut points
#' they differ — the trend is then smooth in interval rank, not in elapsed
#' time. Supply an explicit `formulas$y/d/c` using a clock-time term if a
#' time-scale trend is required.
#'
#' ## Which models per method
#' - `"gformula"`: model_y + model_d
#' - `"ipw_rep1"`: model_d + model_c
#' - `"ipw_rep2"`: model_y + model_c
#'
#' ## Fit populations (no-NA schema)
#' - Y- and D-hazard (cause-specific): rows with `c_event == 0` (at risk
#'   while uncensored; the competing-event exit row stays in with a 0 on
#'   the other cause).
#' - Censoring (`c_event`) hazard: all person-time rows.
#'
#' ## Treatment handling
#' - Y-hazard model: uses `A_y` (a working-copy column on pt_data)
#' - D-hazard model: uses `A_d` (same idea)
#' - Censoring model: uses the observed `treatment` column directly
#'
#' In observed data all three are identical; `A_y`/`A_d` diverge only in
#' cloned datasets used for cross-arm prediction downstream.
#'
#' @family internal
#' @keywords internal
fit_hazard_models <- function(pt_data,
                              treatment,
                              covariates,
                              active_methods,
                              formulas,
                              ipcw = TRUE) {

  # Build default formula components
  cov_terms <- if (length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }

  time_terms <- "k + I(k^2) + I(k^3)"

  models <- list(model_y = NULL, model_d = NULL, model_c = NULL)
  checks <- list(y = NULL, d = NULL, c = NULL)

  # Fit populations (no-NA schema, ordering C -> D -> Y): the Y-hazard
  # conditions on D_{k+1} = 0, so same-interval D exits leave its risk
  # set (mirrors CS's two-level rule y_rows = !indep_cens & !cond_indep_cens).
  # The D-hazard is at risk while uncensored; the censoring hazard is
  # fit on all person-time rows.
  y_rows <- pt_data$c_event == 0 & pt_data$d_event == 0
  d_rows <- pt_data$c_event == 0

  # --- Y-hazard model (needed for gformula and ipw_rep2) ---
  if (any(c("gformula", "ipw_rep2") %in% active_methods)) {
    fml_y <- formulas$y %||% as.formula(
      paste("y_event ~", paste(c("A_y", time_terms, cov_terms),
                               collapse = " + "))
    )
    fit_result <- fit_logistic(
      fml_y, pt_data[y_rows, , drop = FALSE], "Y-hazard")
    models$model_y <- fit_result$model
    checks$y <- fit_result$check
  }

  # --- D-hazard model (needed for gformula and ipw_rep1) ---
  if (any(c("gformula", "ipw_rep1") %in% active_methods)) {
    fml_d <- formulas$d %||% as.formula(
      paste("d_event ~", paste(c("A_d", time_terms, cov_terms),
                               collapse = " + "))
    )
    fit_result <- fit_logistic(
      fml_d, pt_data[d_rows, , drop = FALSE], "D-hazard")
    models$model_d <- fit_result$model
    checks$d <- fit_result$check
  }

  # --- Censoring model (only if IPW requested AND ipcw = TRUE) ---
  if (any(c("ipw_rep1", "ipw_rep2") %in% active_methods) && ipcw) {
    fml_c <- formulas$c %||% as.formula(
      paste("c_event ~", paste(c(treatment, time_terms, cov_terms),
                               collapse = " + "))
    )
    fit_result <- fit_logistic(fml_c, pt_data, "C-hazard")
    models$model_c <- fit_result$model
    checks$c <- fit_result$check
  }

  list(models = models, checks = checks)
}


# ==============================================================================
# Hazard prediction & survival primitives (framework-agnostic)
# ==============================================================================

