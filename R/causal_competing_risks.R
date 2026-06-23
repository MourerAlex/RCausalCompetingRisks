#' Estimate Separable Effects for Competing Events
#'
#' Takes person-time data (typically from
#' [to_person_time_competing()]), fits discrete-time pooled logistic hazard models,
#' and estimates cumulative incidence under each treatment arm (1,1), (0,0),
#' (1,0), and (0,1) via g-formula and/or IPW (Rep 1, Rep 2).
#'
#' Bootstrap confidence intervals are computed separately via [bootstrap()].
#'
#' @param pt_data A data.frame in person-time format. Typically the output
#'   of [to_person_time_competing()] (class `"person_time"`), in which case validation
#'   and metadata are reused. Users with existing person-time data can pass
#'   a plain data.frame, provided it contains the columns `id`, `treatment`,
#'   `k`, `y_event`, `d_event`, `c_event` (and any covariates).
#' @param id Character. Column name of the subject identifier. Defaults to
#'   the value stored in `pt_data`'s attribute when class is `"person_time"`.
#' @param treatment Character. Column name of the binary treatment.
#' @param covariates Character vector of baseline covariate column names.
#' @param method Character. Estimation method: `"all"` (default — runs
#'   g-formula and IPW), or a character vector subset of
#'   `c("gformula", "ipw")`. The token `"ipw"` runs both
#'   IPW representations (Rep 1 and Rep 2) — they are surfaced under the
#'   output keys `ipw_rep1` and `ipw_rep2` for sensitivity comparison.
#' @param formulas Named list or NULL. Override default model formulas. Names
#'   should be `"y"`, `"d"`, and/or `"c"`.
#' @param truncate Either `NULL` (no truncation) or a length-2 numeric
#'   vector of percentile bounds for symmetric weight truncation
#'   (Cole & Hernán 2008 AJE). Default `NULL` (no truncation). Applied to all
#'   source IPW weight columns (`w_cens`, `w_a`, `w_d_arm_*`, `w_y_arm_*`).
#' @param ipcw Logical. Whether to fit a censoring hazard model
#'   and apply inverse probability of censoring weights in IPW estimators.
#'   Defaults to `TRUE`. Set `FALSE` to assume censoring is fully independent
#'   (e.g., purely administrative) and skip the model fit. Has no effect on
#'   g-formula, which always assumes independent censoring given covariates.
#'
#' @return An S3 object of class `"causal_competing_risks_fit"`.
#'
#' @section Architecture:
#' `causal_competing_risks()` is a thin wrapper around [fit_separable_effects()]. The wrapper
#' handles user-facing concerns (argument normalization, input validation,
#' warning capture, S3 class attachment). The worker [fit_separable_effects()] does
#' the estimation math and is called directly by [bootstrap()] for
#' efficiency — no redundant validation or wrapping inside the replicate
#' loop.
#'
#' @section Identification:
#' The target estimand is the cumulative incidence
#' `P(Y^{a_Y, a_D, c_bar=0}_{K+1} = 1)` for each of the four arms
#' `(a_Y, a_D)`. Identification under the separable-effects framework
#' requires:
#'
#' \itemize{
#'   \item The **generalized decomposition assumption (GDA)**: the
#'     observed treatment `A` decomposes into components `A_Y` and `A_D`
#'     that could in principle be assigned different values in a future
#'     four-arm trial.
#'   \item Standard causal conditions — **treatment exchangeability**
#'     (E1, holds by design in an RCT), **censoring exchangeability**
#'     (E2, ignorable censoring), **consistency**, and **positivity**.
#'   \item The two **dismissible component conditions** —
#'     **Δ1**: the cause-specific hazard of `Y` does not depend on `a_D`
#'     given covariate history; and **Δ2**: the cause-specific hazard
#'     of `D` does not depend on `a_Y` given covariate history.
#' }
#'
#' GDA, Δ1, and Δ2 are untestable from observed data — they are causal
#' / biological claims, not statistical properties. Run [assumptions()]
#' on the fit to inspect the full identification block, including the
#' two decompositions (A and B) of the total effect. A numeric isolation
#' read-out (swap-weight ranges) is planned but currently deferred
#' pending a math check; see `@section Isolation read-out (deferred)`
#' on [assumptions()].
#'
#' @section Censoring handling:
#' Censoring is treated differently by g-formula and IPW:
#'
#' - **G-formula** always assumes **independent censoring given the
#'   adjustment covariates**.
#' - **IPW** can optionally apply inverse probability of censoring weights
#'   (IPCW) via the `ipcw` argument.
#'
#' @section Positivity and zero-hazard safeguards:
#' See [fit_separable_effects()] for detailed notes on how hazards near 0 or 1 are
#' handled by the Rep 2 weight construction and the `check_fitted_positivity()`
#' warning system.
#'
#' @section IPW weight columns:
#' The IPW pipeline writes the following per-row weight columns onto the
#' working person-time data (available as `fit$weights$pt_data_weighted`
#' and summarised in `diagnostic(fit)$weight_summary`). Subscripts use
#' `j` for interval index and `k` for the index of the row.
#'
#' \describe{
#'   \item{`w_cens`}{Inverse probability of remaining uncensored
#'     through interval `k`, from the censoring-hazard model `C ~ A + L`:
#'     \eqn{w_{cens, k} = 1 / \prod_{j \leq k} P(C_{j+1} = 0 \,|\, \text{history}_j, A, L)}.
#'     Populated only when `ipcw = TRUE` and at least one IPW
#'     representation is requested. Used by both IPW Rep 1 and Rep 2.}
#'   \item{`w_a`}{Treatment (propensity) weight from the propensity
#'     model `A ~ L`. Unstabilized:
#'     \eqn{w_{a} = 1 / P(A = a \,|\, L)}; stabilized (when
#'     `formulas$A_num` is provided): the marginal-over-conditional ratio.
#'     Populated whenever any IPW method runs. Used by both Reps.}
#'   \item{`w_d_arm_10`, `w_d_arm_01`}{D-component swap weights for
#'     Rep 1 (D-swap representation). For arm `(a_Y, a_D)` with
#'     `a_Y \neq a_D`, Rep 1 stands in `A = a_Y` and reweights to
#'     swap the D-cumulative-survival ratio:
#'     \eqn{w_{D, k}^{(a_Y, a_D)} =
#'       \prod_{j \leq k}
#'       P(D_{j+1} = 0 \,|\, \text{hist}_j, A = a_D) /
#'       P(D_{j+1} = 0 \,|\, \text{hist}_j, A = a_Y)}.
#'     `arm_10` is `(a_Y = 1, a_D = 0)`; `arm_01` is `(a_Y = 0, a_D = 1)`.
#'     `NA_real_` on rows whose `A` does not equal the stand-in.}
#'   \item{`w_y_arm_10`, `w_y_arm_01`}{Y-component swap weights for
#'     Rep 2 (Y-swap representation). For arm `(a_Y, a_D)` with
#'     `a_Y \neq a_D`, Rep 2 stands in `A = a_D` and reweights using a
#'     hazard-ratio times lagged-cumulative-survival ratio:
#'     \eqn{w_{Y, k}^{(a_Y, a_D)} =
#'       \frac{P(Y_{k+1} = 1 \,|\, \text{hist}_k, A = a_Y)}{P(Y_{k+1} = 1 \,|\, \text{hist}_k, A = a_D)}
#'       \cdot
#'       \frac{\prod_{j < k} P(Y_{j+1} = 0 \,|\, \text{hist}_j, A = a_Y)}{\prod_{j < k} P(Y_{j+1} = 0 \,|\, \text{hist}_j, A = a_D)}}.
#'     The denominator hazards are floored at `sqrt(.Machine$double.eps)`
#'     to avoid `Inf`; rescued rows tend to surface as truncated weights
#'     in `diagnostic(fit)$flagged_log` (their pre-truncation magnitude
#'     puts them past the `truncate` percentile bounds).
#'     `NA_real_` on rows whose `A` does not equal the stand-in.}
#' }
#'
#' All four swap-weight columns and `w_cens` / `w_a` are subjected to
#' the symmetric percentile truncation specified by the `truncate`
#' argument (Cole & Hernán 2008 AJE) before they are multiplied into
#' the per-arm Hajek estimator.
#'
#' @references
#' Stensrud MJ, Young JG, Didelez V, Robins JM, Hernan MA (2020).
#' \doi{10.1080/01621459.2020.1765783}
#'
#' Stensrud MJ, Hernan MA, Tchetgen Tchetgen EJ, Robins JM, Didelez V,
#' Young JG (2021). \doi{10.1007/s10985-021-09530-8}
#'
#' @family main
#' @seealso [to_person_time_competing()], [bootstrap()], [risk()], [contrast()],
#'   [diagnostic()], [assumptions()], [fit_separable_effects()]
#'
#' @examples
#' \dontrun{
#' data(prostate_data)
#' pt <- to_person_time_competing(
#'   prostate_data,
#'   id = "id", time = "event_time", event = "event_type",
#'   treatment = "A",
#'   covariates = c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
#'   event_y = 1, event_d = 2, event_c = 0
#' )
#' fit <- causal_competing_risks(pt)
#' plot(risk(fit))
#'
#' # Add bootstrap CIs
#' boot <- bootstrap(fit, n_boot = 500)
#' plot(risk(fit), ci = boot)
#' contrast(fit, method = "gformula", ci = boot)
#' }
#'
#' @export
causal_competing_risks <- function(pt_data,
                      id = NULL,
                      treatment = NULL,
                      covariates = NULL,
                      method = "all",
                      formulas = NULL,
                      truncate = NULL,
                      ipcw = TRUE) {

  cl <- match.call()
  # Freeze argument VALUES into the stored call so bootstrap()'s
  # eval(fit$call) is self-contained: an argument passed as a variable
  # (e.g. method = chosen_method inside a wrapper function) must not need
  # to still be in scope at bootstrap time. pt_data is left as its original
  # symbol — bootstrap overwrites it per replicate, and freezing it would
  # copy the whole dataset into the call.
  pf <- parent.frame()
  for (nm in setdiff(names(cl)[-1L], "pt_data")) {
    cl[[nm]] <- eval(cl[[nm]], pf)
  }

  # --- Basic input shape checks ---
  validate_input_shape(pt_data, "pt_data")

  # --- Normalize method into active_methods vector ---
  # User-facing surface: c("gformula", "ipw"). The token "ipw" runs both
  # representations: "ipw_rep1" (Rep 1,
  # D-swap) and "ipw_rep2" (Rep 2, Y-swap). The expanded names are what
  # the rest of the pipeline reads on `active_methods` and what the
  # output keys are named on `cumulative_incidence`.
  valid_user_methods <- c("gformula", "ipw")
  if (identical(method, "all")) {
    user_methods <- valid_user_methods
  } else {
    user_methods <- method
  }
  if (!is.character(user_methods) ||
      length(user_methods) < 1 ||
      !all(user_methods %in% valid_user_methods)) {
    stop(
      "Unknown method entries: ",
      paste(setdiff(user_methods, valid_user_methods), collapse = ", "),
      ". Must be 'all' or a character vector of any subset of: ",
      paste(valid_user_methods, collapse = ", "),
      call. = FALSE
    )
  }
  # Expand "ipw" -> c("ipw_rep1", "ipw_rep2") for internal dispatch.
  active_methods <- unlist(lapply(user_methods, function(m) {
    if (m == "ipw") c("ipw_rep1", "ipw_rep2") else m
  }), use.names = FALSE)

  # --- Validate formulas keys (typo defense) ---
  # `formulas` is an optional named list of model-formula overrides.
  # Valid keys: y / d / c (hazard models) and A / A_num (propensity:
  # full and stabilization-numerator). Any other key is almost
  # certainly a typo and would silently fall back to the default
  # formula — error early instead.
  valid_formula_keys <- c("y", "d", "c", "A", "A_num")
  if (!is.null(formulas)) {
    if (!is.list(formulas) || is.null(names(formulas)) ||
        any(names(formulas) == "")) {
      stop(
        "`formulas` must be a named list (keys: ",
        paste(valid_formula_keys, collapse = ", "), ").",
        call. = FALSE
      )
    }
    bad_keys <- setdiff(names(formulas), valid_formula_keys)
    if (length(bad_keys) > 0) {
      stop(
        "Unknown formula key(s): ",
        paste(shQuote(bad_keys), collapse = ", "),
        ". Valid keys: ",
        paste(shQuote(valid_formula_keys), collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  # --- Resolve input origin: class + attributes vs. explicit args ---
  if (inherits(pt_data, "person_time")) {
    # Case 1: data went through to_person_time_competing() — everything is on the
    # attributes. Passing id / treatment / covariates on top is ambiguous
    # and almost always a mistake, so we reject it.
    if (!is.null(id) || !is.null(treatment) || !is.null(covariates)) {
      stop(
        "pt_data is a 'person_time' (from to_person_time_competing()). ",
        "id / treatment / covariates are already stored as attributes. ",
        "Do not pass them again.",
        call. = FALSE
      )
    }
    cut_times      <- attr(pt_data, "cut_times")
    id_col         <- attr(pt_data, "id_col")
    treatment_col  <- attr(pt_data, "treatment_col")
    covariates_vec <- attr(pt_data, "covariates")
  } else {
    # Case 2: user brought their own person-time data without the class
    if (is.null(id) || is.null(treatment)) {
      stop(
        "When pt_data is not from to_person_time_competing(), you must specify ",
        "'id' and 'treatment' column names.",
        call. = FALSE
      )
    }
    covariates_vec <- covariates %||% character()
    validate_person_time(pt_data, id, treatment, covariates_vec)

    id_col        <- id
    treatment_col <- treatment
    cut_times     <- sort(unique(pt_data$k))

    if (!"A_y" %in% names(pt_data)) pt_data$A_y <- pt_data[[treatment_col]]
    if (!"A_d" %in% names(pt_data)) pt_data$A_d <- pt_data[[treatment_col]]
  }

  # --- Warning capture around the worker ---
  collected_warnings <- character()
  out <- withCallingHandlers(
    fit_separable_effects(
      pt_data = pt_data,
      id_col = id_col,
      treatment_col = treatment_col,
      covariates_vec = covariates_vec,
      cut_times = cut_times,
      active_methods = active_methods,
      formulas = formulas,
      ipcw = ipcw,
      truncate = truncate
    ),
    warning = function(w) {
      collected_warnings <<- c(collected_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # One grouped notice (mirrors CS .with_collected_warnings) instead of
  # re-emitting each glm warning individually. The full list lives on
  # fit$warnings and is surfaced by print().
  if (length(collected_warnings) > 0L) {
    warning("causal_competing_risks(): ", length(collected_warnings),
            " warning(s) collected during fit. See fit$warnings for details.",
            call. = FALSE)
  }

  # --- Assemble S3 result ---
  structure(
    list(
      cumulative_incidence = out$cumulative_incidence,
      weights              = out$weights,
      gformula             = NULL,  # v1: reserved for g-formula diagnostics
      models               = out$models,
      model_checks         = out$model_checks,
      person_time          = pt_data,
      warnings             = collected_warnings,
      call                 = cl,
      method               = method,
      active_methods       = active_methods,
      times                = cut_times,
      id_col               = id_col,
      treatment_col        = treatment_col,
      covariates           = covariates_vec,
      formulas             = formulas,
      ipcw    = ipcw,
      truncate             = truncate,
      n                    = length(unique(pt_data[[id_col]]))
    ),
    class = "causal_competing_risks_fit"
  )
}


#' Internal Worker: Fit and Estimate (No Class, No Warnings Wrapping)
#'
#' Internal worker for [causal_competing_risks()]: fits hazard models, dispatches
#' to estimators, returns a plain list. No validation, no S3 class, no
#' warning capture — those are the wrapper's concerns.
#'
#' [bootstrap()] calls this directly inside its replicate loop to avoid the
#' wrapper overhead. Per-replicate warnings are suppressed at the call site.
#'
#' @param pt_data Person-time data.frame (must be fully prepared).
#' @param id_col Character. Subject id column name.
#' @param treatment_col Character. Treatment column name.
#' @param covariates_vec Character vector of covariate column names.
#' @param cut_times Numeric vector of interval cut points.
#' @param active_methods Character vector subset of
#'   `c("gformula", "ipw_rep1", "ipw_rep2")`.
#' @param formulas Named list or NULL. Model formula overrides.
#' @param ipcw Logical. Whether to fit censoring model for IPW.
#' @param truncate Either NULL or length-2 numeric percentile bounds for
#'   symmetric weight truncation.
#'
#' @return A plain list with:
#'   \describe{
#'     \item{cumulative_incidence}{Named list (one entry per method).}
#'     \item{weights}{NULL if no IPW; otherwise list with pt_data_weighted,
#'       weight_summary, flagged_ids, and adjustment settings.}
#'     \item{models}{Named list: model_y, model_d, model_c.}
#'     \item{model_checks}{Named list: y, d, c (per-model diagnostics).}
#'   }
#'
#' @family internal
#' @keywords internal
fit_separable_effects <- function(pt_data,
                          id_col,
                          treatment_col,
                          covariates_vec,
                          cut_times,
                          active_methods,
                          formulas,
                          ipcw,
                          truncate = NULL,
                          stabilize = TRUE) {

  # --- Fit hazard models (Y, D, C) ---
  hazmod <- fit_hazard_models(
    pt_data = pt_data,
    treatment = treatment_col,
    covariates = covariates_vec,
    active_methods = active_methods,
    formulas = formulas,
    ipcw = ipcw
  )
  models <- hazmod$models
  model_checks <- hazmod$checks

  # --- Fit propensity model only when IPW will run ---
  # G-formula doesn't consume model_a; conditional fit avoids wasted work.
  if (any(c("ipw_rep1", "ipw_rep2") %in% active_methods)) {
    prop <- fit_propensity(
      pt_data = pt_data,
      treatment = treatment_col,
      covariates = covariates_vec,
      stabilize = stabilize,
      formula_full = formulas$A,
      formula_num  = formulas$A_num
    )
    models$model_a     <- prop$model_a
    models$model_a_num <- prop$model_a_num
    model_checks$a     <- prop$check_a
    model_checks$a_num <- prop$check_a_num
  }

  # --- Run each requested estimator ---
  ci_list <- list()
  weights_slot <- NULL

  if ("gformula" %in% active_methods) {
    ci_list$gformula <- fit_gformula(
      pt_data = pt_data, models = models,
      cut_times = cut_times, id_col = id_col
    )
  }

  ipw_reps <- intersect(active_methods, c("ipw_rep1", "ipw_rep2"))
  if (length(ipw_reps) > 0) {
    ipw_res <- fit_ipw(
      pt_data = pt_data, models = models,
      treatment = treatment_col,
      ipw_reps = ipw_reps, cut_times = cut_times, id_col = id_col,
      truncate = truncate
    )

    if ("ipw_rep1" %in% ipw_reps &&
        !is.null(ipw_res$cumulative_incidence_rep1)) {
      ci_list$ipw_rep1 <- ipw_res$cumulative_incidence_rep1
    }
    if ("ipw_rep2" %in% ipw_reps &&
        !is.null(ipw_res$cumulative_incidence_rep2)) {
      ci_list$ipw_rep2 <- ipw_res$cumulative_incidence_rep2
    }

    weights_slot <- list(
      pt_data_weighted = ipw_res$pt_data_weighted,
      weight_summary   = ipw_res$weight_summary,
      flagged_ids      = ipw_res$flagged_ids,
      flagged_log      = ipw_res$flagged_log,
      truncate         = truncate
    )
  }

  list(
    cumulative_incidence = ci_list,
    weights              = weights_slot,
    models               = models,
    model_checks         = model_checks
  )
}
