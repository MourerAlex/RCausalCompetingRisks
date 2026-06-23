#' Identifying Assumptions for a Separable-Effects
#'
#' Returns a structured object listing the target estimand and the
#' identifying assumptions under which it is identified from the observed
#' data, plus a numeric isolation read-out (no verdict).
#'
#' The accessor does not run any new computation on the fit beyond
#' summarising the per-row swap weights for the isolation slot. All
#' statements are textual and refer to the estimand and the conditions
#' formalised in Stensrud et al. (2020, 2021).
#'
#' @param fit A `"causal_competing_risks_fit"` object from [causal_competing_risks()].
#'
#' @return An S3 object of class `"causal_competing_risks_assumptions"` with:
#'   \describe{
#'     \item{estimand}{Named list. Target counterfactual notation, the
#'       four arms `(a_Y, a_D)`, and the two decompositions (A and B)
#'       expressed as differences between arm cumulative incidence.}
#'     \item{assumptions}{Named list with one entry per identifying
#'       assumption (GDA, treatment exchangeability E1, censoring
#'       exchangeability E2, consistency, positivity-treatment,
#'       positivity-censoring, no interference, correct model
#'       specification, full isolation). Each entry is a list with
#'       fields `name`, `statement`, `formula` (LaTeX or `NA`),
#'       `status` (`"testable"` / `"untestable"`), `pointer` (diagnostic
#'       accessor or `NA`), and `citation`.}
#'     \item{dismissible_D1}{Named list (same 6-field schema) for the
#'       dismissible component condition Δ1 (cause-specific hazard of Y
#'       does not depend on `a_D`).}
#'     \item{dismissible_D2}{Named list for Δ2 (cause-specific hazard
#'       of D does not depend on `a_Y`).}
#'   }
#'
#' @section Isolation sensitivity (deferred — needs time-varying L):
#' The `full isolation` assumption is recorded but marked `untestable`
#' in v1. Stensrud et al. (2021, Sect. 7) give a falsification: compare
#' the full-isolation (2020) estimate against the generalized
#' Z_k-partition estimator. That comparison requires measured
#' time-varying covariates `L_k, k in {0,...,K}` (equivalently,
#' post-treatment intermediate covariates affected by a component). v1
#' takes baseline confounders only, so the generalized estimator is not
#' computable and the two formulas coincide.
#'
#' @section Isolation read-out (deferred):
#' The locked API spec calls for a numeric isolation read-out — per-row
#' swap-weight ranges `w_d` and `w_y`, no verdict — both in the
#' `$isolation` slot of this object and as a one-line header on
#' `print(fit)`. The empirical mapping `w_d ≈ 1 ↔ Δ2` (and `w_y ≈ 1 ↔
#' Δ1`) needs to be double-checked against the Stensrud Appendix.
#' Until that check is done, the read-out is intentionally omitted to
#' avoid surfacing an unverified causal claim. The textual statements
#' for Δ1 and Δ2 (with diagnostic pointers) remain in place.
#'
#' @references
#' Stensrud MJ, Young JG, Didelez V, Robins JM, Hernán MA (2020).
#' Separable Effects for Causal Inference in the Presence of Competing
#' Events. \doi{10.1080/01621459.2020.1765783}
#'
#' Stensrud MJ, Hernán MA, Tchetgen Tchetgen EJ, Robins JM, Didelez V,
#' Young JG (2021). A generalized theory of separable effects in
#' competing event settings. \doi{10.1007/s10985-021-09530-8}
#'
#' @seealso [causal_competing_risks()], [risk()], [contrast()], [diagnostic()]
#' @family accessors
#' @export
assumptions <- function(fit) {
  stopifnot(inherits(fit, "causal_competing_risks_fit"))

  estimand <- list(
    target = "P(Y^{a_Y, a_D, c_bar = 0}_{K+1} = 1)",
    description = paste0(
      "Cumulative incidence of the primary event Y by interval K+1 ",
      "under hypothetical assignment to the four-arm regime ",
      "(a_Y, a_D), with censoring eliminated."
    ),
    arms = c("(1,1)", "(0,0)", "(1,0)", "(0,1)"),
    decomposition_A = list(
      label = "Decomposition A (vary A_Y first; intermediate arm = arm_10)",
      sde = "SDE-A = F^(1,0) - F^(0,0) = effect of A_Y at a_D = 0",
      sie = "SIE-A = F^(1,1) - F^(1,0) = effect of A_D at a_Y = 1"
    ),
    decomposition_B = list(
      label = "Decomposition B (vary A_D first; intermediate arm = arm_01)",
      sie = "SIE-B = F^(0,1) - F^(0,0) = effect of A_D at a_Y = 0",
      sde = "SDE-B = F^(1,1) - F^(0,1) = effect of A_Y at a_D = 1"
    )
  )

  assumptions_list <- list(

    GDA = list(
      name      = "Generalized decomposition (GDA)",
      statement = paste0(
        "The treatment A can be decomposed into two binary components ",
        "A_Y in {0,1} and A_D in {0,1} such that, in the observed data, ",
        "the determinism A = A_D = A_Y holds, but in a future study A_Y ",
        "and A_D could, in principle, be assigned different values."
      ),
      formula   = "Y^{a_Y = a,\\, a_D = a} = Y^{A = a}",
      status    = "untestable",
      pointer   = NA_character_,
      citation  = "Stensrud et al. (2020), eq. (2); Stensrud et al. (2021)."
    ),

    E1_treatment_exchangeability = list(
      name      = "Treatment exchangeability (E1)",
      statement = paste0(
        "What would happen to everyone under no censoring is independent ",
        "of which treatment was received, given covariates L. Holds by ",
        "randomisation; in observational data it is the ",
        "no-unmeasured-confounding condition."
      ),
      formula   = "(\\bar{Y}^{a,\\bar{c}=0}_{K+1},\\, \\bar{D}^{a,\\bar{c}=0}_{K+1}) \\perp A \\mid L",
      status    = "untestable",
      pointer   = NA_character_,
      citation  = "Stensrud et al. (2020), Section 3.2."
    ),

    E2_censoring_exchangeability = list(
      name      = "Censoring exchangeability (E2)",
      statement = paste0(
        "What would happen at k+1 under no censoring is independent of ",
        "whether censoring occurs at k+1, given observed history. The ",
        "ignorable-censoring condition, which IPCW attempts to address when ipcw = TRUE."
      ),
      formula   = "(Y^{a,\\bar{c}=0}_{k+1},\\, D^{a,\\bar{c}=0}_{k+1}) \\perp C_{k+1} \\mid Y_k = D_k = \\bar{C}_k = 0,\\, L,\\, A",
      status    = "untestable",
      pointer   = "diagnostic(fit)$model_checks$c",
      citation  = "Stensrud et al. (2020), Section 3.3."
    ),

    consistency = list(
      name      = "Consistency",
      statement = paste0(
        "If A = a and C_bar_k = 0, the counterfactual equals the observed ",
        "value. Only applies in the observed arms (a_Y = a_D) and among ",
        "uncensored individuals."
      ),
      formula   = "A = a,\\ \\bar{C}_k = 0 \\;\\Rightarrow\\; Y^{a,\\bar{c}=0}_k = Y_k,\\ D^{a,\\bar{c}=0}_k = D_k",
      status    = "untestable",
      pointer   = NA_character_,
      citation  = "Stensrud et al. (2020), Section 3.4."
    ),

    no_interference = list(
      name      = "No interference",
      statement = paste0(
        "One subject's treatment does not affect another subject's ",
        "potential outcomes. Together with consistency, this is SUTVA ",
        "(the stable-unit-treatment-value assumption)."
      ),
      formula   = NA_character_,
      status    = "untestable",
      pointer   = NA_character_,
      citation  = "Hernan & Robins (2020), What If, Chapter 1."
    ),

    positivity_treatment = list(
      name      = "Positivity (treatment)",
      statement = paste0(
        "For every covariate pattern with positive density, both ",
        "treatments occur with positive probability. Else the identifying ",
        "formulas divide by zero."
      ),
      formula   = "P(L = l) > 0 \\;\\Rightarrow\\; P(A = a \\mid L = l) > 0,\\ a \\in \\{0,1\\}",
      status    = "testable",
      pointer   = "diagnostic(fit)$model_checks$a$min_fitted / $max_fitted; weight_summary (IPW)",
      citation  = "Stensrud et al. (2020), Section 3.5."
    ),

    positivity_censoring = list(
      name      = "Positivity (censoring)",
      statement = paste0(
        "Among event-free uncensored individuals at each interval, some ",
        "remain under observation: P(remain uncensored | history) > 0."
      ),
      formula   = "P(A = a, Y_k = D_k = \\bar{C}_k = 0, L = l) > 0 \\;\\Rightarrow\\; P(C_{k+1} = 0 \\mid \\ldots) > 0",
      status    = "testable",
      pointer   = "diagnostic(fit)$model_checks$c$min_fitted / $max_fitted",
      citation  = "Stensrud et al. (2020), Section 3.5."
    ),

    model_specification = list(
      name      = "Correct model specification",
      statement = paste0(
        "The fitted discrete-time hazard models for Y, D, and (under ",
        "ipcw = TRUE) C, and the propensity model for A, are correctly ",
        "specified."
      ),
      formula   = NA_character_,
      status    = "untestable",
      pointer   = "diagnostic(fit)$model_checks",
      citation  = "Hernan & Robins (2020), What If, Chapter 18."
    ),

    isolation_full = list(
      name      = "Full isolation",
      statement = paste0(
        "No measured covariate is a post-treatment intermediate that is ",
        "affected by one treatment component and also affects the other ",
        "component's outcome (no shared intermediate). v1 targets the 2020 ",
        "full-isolation setting with baseline confounders only. Falsifiable ",
        "only via the 2021 Z_k-partition sensitivity analysis, which ",
        "requires measured time-varying covariates L_k, k in {0,...,K} ",
        "(Stensrud et al. 2021, Sect. 7); with baseline confounders only ",
        "the generalized estimator is not computable, so there is nothing ",
        "to compare and full isolation cannot be tested here."
      ),
      formula   = NA_character_,
      status    = "untestable",
      pointer   = NA_character_,
      citation  = "Stensrud et al. (2021), Section 3, Section 7, Appendix C."
    )
  )

  dismissible_D1 <- list(
    name      = "Dismissible component Delta-1 (for Y)",
    statement = paste0(
      "The cause-specific hazard of Y among event-free individuals does ",
      "not depend on a_D, conditional on the covariate history. A causal ",
      "claim, not a statistical property of the data."
    ),
    formula   = "P(Y^{a_Y, a_D=1, \\bar{c}=0}_{k+1} = 1 \\mid Y_k = 0, D_{k+1} = 0, L) = P(Y^{a_Y, a_D=0, \\bar{c}=0}_{k+1} = 1 \\mid Y_k = 0, D_{k+1} = 0, L)",
    status    = "untestable",
    pointer   = NA_character_,
    citation  = "Stensrud et al. (2020), Section 3.6 (Delta-1); Stensrud et al. (2021), eq. (34)."
  )

  dismissible_D2 <- list(
    name      = "Dismissible component Delta-2 (for D)",
    statement = paste0(
      "The cause-specific hazard of D among event-free individuals does ",
      "not depend on a_Y, conditional on the covariate history. (Same ",
      "structure as Delta-1, roles of Y and D swapped.)"
    ),
    formula   = "P(D^{a_Y=1, a_D, \\bar{c}=0}_{k+1} = 1 \\mid \\ldots) = P(D^{a_Y=0, a_D, \\bar{c}=0}_{k+1} = 1 \\mid \\ldots)",
    status    = "untestable",
    pointer   = NA_character_,
    citation  = "Stensrud et al. (2020), Section 3.6 (Delta-2); Stensrud et al. (2021), eq. (35)."
  )

  # --- Isolation slot intentionally omitted; see @section above. ---
  # The helper isolation_summary() is kept for the future re-enable.

  structure(
    list(
      estimand        = estimand,
      assumptions     = assumptions_list,
      dismissible_D1  = dismissible_D1,
      dismissible_D2  = dismissible_D2
    ),
    class = "causal_competing_risks_assumptions"
  )
}


#' Compute Isolation Summary (No Verdict) — currently unused
#'
#' Pulls the per-row swap weights from the fit and reports their numeric
#' range. Returns NULL fields if no IPW method ran (i.e., gformula-only).
#'
#' Currently NOT consumed by [assumptions()] or [print.causal_competing_risks_fit()];
#' kept for the future re-enable once the empirical link
#' `w_d ≈ 1 ↔ Δ2` (and `w_y ≈ 1 ↔ Δ1`) is math-checked against the
#' Stensrud Appendix. See `@section Isolation read-out (deferred)` on
#' [assumptions()].
#'
#' @param fit A `"causal_competing_risks_fit"` object.
#' @return A named list with `range_w_d`, `range_w_y`, `n_rows_w_d`,
#'   `n_rows_w_y`, and `note`.
#' @family internal
#' @keywords internal
isolation_summary <- function(fit) {

  pt <- fit$weights$pt_data_weighted

  range_or_null <- function(cols) {
    if (is.null(pt)) return(list(range = NULL, n = 0L))
    present <- intersect(cols, names(pt))
    if (length(present) == 0L) return(list(range = NULL, n = 0L))
    vals <- unlist(lapply(present, function(co) pt[[co]]), use.names = FALSE)
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0L) return(list(range = NULL, n = 0L))
    list(range = range(vals), n = length(vals))
  }

  rd <- range_or_null(c("w_d_arm_10", "w_d_arm_01"))
  ry <- range_or_null(c("w_y_arm_10", "w_y_arm_01"))

  note <- paste0(
    "Isolation cannot be falsified from data. Component hazard ratios ",
    "(swap weights) near 1 are consistent with — but do not prove — Δ1 ",
    "(for w_y) and Δ2 (for w_d). Sensitivity analysis under partial ",
    "isolation is the appropriate response when the ranges depart ",
    "meaningfully from 1; see Stensrud & Young (2021)."
  )

  list(
    range_w_d  = rd$range,
    range_w_y  = ry$range,
    n_rows_w_d = rd$n,
    n_rows_w_y = ry$n,
    note       = note
  )
}


#' Print a causal_competing_risks_assumptions Object
#'
#' Renders an identification block: estimand, the standard
#' identifying assumptions, and the two dismissible component conditions
#' Δ1 and Δ2. Numeric isolation read-out is currently deferred; see
#' `@section Isolation read-out (deferred)` on [assumptions()].
#'
#' @param x A `"causal_competing_risks_assumptions"` object.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_competing_risks_assumptions <- function(x, ...) {

  cat("Identifying Assumptions (causal_competing_risks)\n")
  cat("-------------------------------------------\n\n")

  # --- Estimand ---
  cat("Estimand:\n")
  cat("  ", x$estimand$target, "\n", sep = "")
  cat("  ", x$estimand$description, "\n", sep = "")
  cat("\nArms (a_Y, a_D): ",
      paste(x$estimand$arms, collapse = ", "), "\n", sep = "")
  cat("\n", x$estimand$decomposition_A$label, "\n", sep = "")
  cat("  ", x$estimand$decomposition_A$sde, "\n", sep = "")
  cat("  ", x$estimand$decomposition_A$sie, "\n", sep = "")
  cat("\n", x$estimand$decomposition_B$label, "\n", sep = "")
  cat("  ", x$estimand$decomposition_B$sie, "\n", sep = "")
  cat("  ", x$estimand$decomposition_B$sde, "\n", sep = "")

  # --- One-entry renderer (6-field schema) ---
  emit <- function(a) {
    tag <- if (identical(a$status, "testable")) "[testable]" else "[untestable]"
    cat("\n  [", a$name, "]  ", tag, "\n", sep = "")
    cat("    ", a$statement, "\n", sep = "")
    if (!is.null(a$formula) && !is.na(a$formula)) {
      cat("    Formula:  ", a$formula, "\n", sep = "")
    }
    if (!is.null(a$pointer) && !is.na(a$pointer)) {
      cat("    See:      ", a$pointer, "\n", sep = "")
    }
    if (!is.null(a$citation) && !is.na(a$citation)) {
      cat("    Citation: ", a$citation, "\n", sep = "")
    }
  }

  # --- Assumptions ---
  cat("\nIdentifying assumptions:\n")
  for (nm in names(x$assumptions)) emit(x$assumptions[[nm]])

  # --- Dismissible component conditions ---
  for (nm in c("dismissible_D1", "dismissible_D2")) emit(x[[nm]])

  # --- Isolation read-out: deferred pending math check; see @section. ---

  invisible(x)
}


#' Format a causal_competing_risks_assumptions Object as Markdown
#'
#' Renders the same content as `print()` but as a markdown string
#' suitable for vignette / paper reuse. Avoids drift between package
#' text and external documents.
#'
#' @param x A `"causal_competing_risks_assumptions"` object.
#' @param style Character. Currently only `"markdown"` is supported.
#' @param ... Additional arguments (currently unused).
#' @return A single character string (one document) with markdown headers.
#' @export
format.causal_competing_risks_assumptions <- function(x, style = "markdown", ...) {

  if (!identical(style, "markdown")) {
    stop("Only style = \"markdown\" is supported.", call. = FALSE)
  }

  lines <- character()

  push <- function(...) {
    lines <<- c(lines, paste0(..., collapse = ""))
  }

  # --- Estimand ---
  push("## Estimand")
  push("")
  push("`", x$estimand$target, "`")
  push("")
  push(x$estimand$description)
  push("")
  push("**Arms (a_Y, a_D):** ",
       paste(paste0("`", x$estimand$arms, "`"), collapse = ", "))
  push("")
  push("### ", x$estimand$decomposition_A$label)
  push("")
  push("- ", x$estimand$decomposition_A$sde)
  push("- ", x$estimand$decomposition_A$sie)
  push("")
  push("### ", x$estimand$decomposition_B$label)
  push("")
  push("- ", x$estimand$decomposition_B$sie)
  push("- ", x$estimand$decomposition_B$sde)
  push("")

  # --- Assumptions ---
  push("## Identifying assumptions")
  push("")

  emit_assumption <- function(a) {
    tag <- if (identical(a$status, "testable")) "testable" else "untestable"
    push("### ", a$name)
    push("")
    push(a$statement)
    push("")
    if (!is.null(a$formula) && !is.na(a$formula)) {
      push("$$", a$formula, "$$")
      push("")
    }
    push("- **Status:** ", tag)
    if (!is.null(a$pointer) && !is.na(a$pointer)) {
      push("- **Diagnostic:** ", a$pointer)
    }
    if (!is.null(a$citation) && !is.na(a$citation)) {
      push("- **Citation:** ", a$citation)
    }
    push("")
  }

  for (nm in names(x$assumptions)) {
    emit_assumption(x$assumptions[[nm]])
  }
  emit_assumption(x$dismissible_D1)
  emit_assumption(x$dismissible_D2)

  # --- Isolation: deferred pending math check; see @section. ---

  paste(lines, collapse = "\n")
}
