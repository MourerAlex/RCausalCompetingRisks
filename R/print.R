#' Print a causal_competing_risks Object
#'
#' @param x A `"causal_competing_risks_fit"` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns `x`.
#' @export
print.causal_competing_risks_fit <- function(x, ...) {
  title <- "Separable effects fit (causal_competing_risks_fit)"
  cat(title, "\n", sep = "")
  cat(strrep("-", nchar(title)), "\n", sep = "")

  # One-line identification header. Surfaces the estimand framing only;
  # the swap-weight numeric read-out is intentionally suppressed pending
  # a math check against the Stensrud Appendix — see `@section Isolation
  # read-out (deferred)` on [assumptions()] for context. To restore, see
  # commit history.
  cat("Estimand: P(Y^{a_Y, a_D, c_bar=0}_{K+1}=1)",
      "under separable-effects identification\n")
  cat("See `assumptions(fit)` for the full identification block.\n")

  cat("Method(s): ", paste(names(x$cumulative_incidence), collapse = ", "),
      "\n", sep = "")
  cat("N subjects: ", x$n, "\n", sep = "")
  # CCR's k column holds the cut-time value (not a 1..K index), so the
  # count of cut times is length(times) and the final k is max(times).
  cat("Cut times: ", length(x$times), " (T_max = ", max(x$times), ")\n",
      sep = "")

  # Per-method cumulative incidence at the final cut time
  cat("\nCumulative incidence at k = ", max(x$times), ":\n", sep = "")
  for (m in names(x$cumulative_incidence)) {
    df <- x$cumulative_incidence[[m]]
    last <- df[nrow(df), ]
    cat(sprintf(
      "  [%s]  (1,1)=%.4f  (0,0)=%.4f  (1,0)=%.4f  (0,1)=%.4f\n",
      m, last$arm_11, last$arm_00, last$arm_10, last$arm_01
    ))
  }

  # Model-checks summary (no binary positivity flag — just count the
  # continuous signals we surface)
  if (!is.null(x$model_checks)) {
    issues <- 0
    for (chk in x$model_checks) {
      if (is.null(chk)) next
      if (!isTRUE(chk$converged)) issues <- issues + 1
      if (length(chk$glm_warnings) > 0) issues <- issues + 1
    }
    if (issues > 0) {
      cat("\nModel checks:", issues,
          "issue(s) - use `fit$model_checks` to inspect.\n")
    }
  }

  if (length(x$warnings) > 0) {
    cat("\nFit completed with ", length(x$warnings),
        " warning(s) (see fit$warnings).\n", sep = "")
  }

  cat("\nUse risk(), contrast(), diagnostic() to extract components.\n")
  invisible(x)
}


#' Summary of a causal_competing_risks Object
#'
#' Prints per-method cumulative incidence at the selected time, a brief
#' model checks summary, and (when a bootstrap is supplied) the contrast
#' table at that same time.
#'
#' @param object A `"causal_competing_risks_fit"` object.
#' @param ci Optional. A `"causal_competing_risks_bootstrap"` object from
#'   [bootstrap()]. When provided, the contrast table at the selected
#'   `time` is printed alongside the cumulative incidence summary.
#' @param time Numeric scalar or NULL. Time point at which to summarise.
#'   NULL (default) selects the final cut time (`max(object$times)`). A
#'   user-supplied value is snapped to the nearest cut time and a
#'   `message()` is emitted when snapping changes the value.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the per-method list of cumulative incidence
#'   data frames.
#' @export
summary.causal_competing_risks_fit <- function(object, ci = NULL,
                                               time = NULL, ...) {
  if (!is.null(ci)) {
    stopifnot(inherits(ci, "causal_competing_risks_bootstrap"))
  }
  k_at <- snap_time(time, object$times)

  # 1. Banner + 2. method/cohort info box
  .ccr_summary_banner()
  .ccr_summary_info_box(names(object$cumulative_incidence), object$n,
                        length(object$times), max(object$times))

  # 3. Baseline N per observed treatment arm (one row per subject)
  pt   <- object$person_time
  base <- pt[!duplicated(pt[[object$id_col]]), , drop = FALSE]
  arm_levels <- sort(unique(base[[object$treatment_col]]))
  n_per_arm  <- vapply(arm_levels,
                       function(a) sum(base[[object$treatment_col]] == a),
                       integer(1))
  .ccr_summary_baseline(arm_levels, n_per_arm)

  # 4. Counterfactual risk (4 arms, F^a and S^a) per method at k_at
  for (m in names(object$cumulative_incidence)) {
    row <- object$cumulative_incidence[[m]]
    row <- row[row$k == k_at, , drop = FALSE]
    if (nrow(row) == 0L) next
    .ccr_summary_risk(m, row, k_at)
  }

  # 5. Contrasts with proportional CI bars (RD), per method, when ci given
  if (!is.null(ci)) {
    for (m in names(object$cumulative_incidence)) {
      if (!m %in% unique(ci$replicates$method)) next
      ctr <- compute_contrasts(object, method = m, ci = ci)
      ctr <- ctr[ctr$k == k_at & ctr$measure == "rd", , drop = FALSE]
      .ccr_summary_contrasts(ctr, m, k_at, ci)
    }
  }

  # 6. Footer + 7. model-checks tally
  .ccr_summary_footer(ci)
  .ccr_summary_model_checks(object$model_checks)

  invisible(object$cumulative_incidence)
}


# ----------------------------------------------------------------------------
# Internal tile renderers for summary.causal_competing_risks_fit(): an
# ASCII box-drawn console layout mirroring CausalSurvival's summary (banner,
# info box, baseline bars, per-method counterfactual risk, separable
# contrasts with proportional CI bars, footer, model-checks tally).
# ----------------------------------------------------------------------------

#' @keywords internal
.ccr_summary_banner <- function() {
  title   <- "SEPARABLE EFFECTS - SUMMARY"
  inner_w <- 70L
  side    <- (inner_w - nchar(title)) %/% 2L
  middle  <- paste0(strrep(" ", side), title,
                    strrep(" ", inner_w - side - nchar(title)))
  bar     <- strrep("=", inner_w)
  cat("\n+", bar, "+\n", sep = "")
  cat("|", middle, "|\n", sep = "")
  cat("+", bar, "+\n\n", sep = "")
}

#' @keywords internal
.ccr_summary_info_box <- function(methods, N, K, T_max) {
  l1 <- sprintf("Method(s): %s", paste(methods, collapse = ", "))
  l2 <- sprintf("N = %d  |  K = %d  |  T_max = %g", N, K, T_max)
  cw <- max(nchar(l1), nchar(l2)); iw <- cw + 4L
  pad <- function(s) paste0(s, strrep(" ", cw - nchar(s)))
  ind <- strrep(" ", 11L); bar <- strrep("-", iw)
  cat(ind, "+", bar, "+\n", sep = "")
  cat(ind, "|  ", pad(l1), "  |\n", sep = "")
  cat(ind, "|  ", pad(l2), "  |\n", sep = "")
  cat(ind, "+", bar, "+\n\n", sep = "")
}

#' @keywords internal
.ccr_summary_baseline <- function(arm_levels, n_per_arm) {
  cat("  ", strrep("-", 66), "\n", sep = "")
  cat("  BASELINE   at t = 0\n")
  cat("  ", strrep("-", 66), "\n\n", sep = "")
  bar_max <- 20L
  for (i in seq_along(arm_levels)) {
    bar_len <- max(1L, round(n_per_arm[i] / max(n_per_arm) * bar_max))
    cat(sprintf("      arm %s  %s%s   N = %3d\n",
                format(arm_levels[i]),
                strrep("#", bar_len), strrep(" ", bar_max - bar_len),
                n_per_arm[i]))
  }
  cat("\n")
}

#' @keywords internal
.ccr_summary_risk <- function(method, row, k_at) {
  cat("  ", strrep("-", 66), "\n", sep = "")
  cat(sprintf("  COUNTERFACTUAL RISK [%s]   at k = %g\n", method, k_at))
  cat("  ", strrep("-", 66), "\n\n", sep = "")
  cat("      arm (a_Y,a_D)    F^a(t)    S^a(t)\n")
  arms <- c("arm_11", "arm_00", "arm_10", "arm_01")
  labs <- c("(1,1)", "(0,0)", "(1,0)", "(0,1)")
  for (j in seq_along(arms)) {
    if (!arms[j] %in% names(row)) next
    f <- row[[arms[j]]]
    cat(sprintf("          %s        %.3f     %.3f\n", labs[j], f, 1 - f))
  }
  cat("\n")
}

#' Render the separable-contrast section: RD per contrast + proportional bar
#' @keywords internal
.ccr_summary_contrasts <- function(ctr, method, k_at, ci) {
  cat("  ", strrep("-", 66), "\n", sep = "")
  cat(sprintf("  CONTRASTS [%s] on cumulative incidence (RD, k = %g)\n",
              method, k_at))
  cat(sprintf("  (%.0f%% CIs from %d bootstrap replicates)\n",
              (1 - ci$alpha) * 100, ci$n_boot_effective))
  cat("  ", strrep("-", 66), "\n\n", sep = "")

  any_null <- FALSE
  for (i in seq_len(nrow(ctr))) {
    r     <- ctr[i, ]
    decomp <- if (is.na(r$decomp)) "" else sprintf(" (decomp %s)", r$decomp)
    label  <- sprintf("%s%s", r$contrast, decomp)
    inc_null <- r$lower <= 0 && 0 <= r$upper
    if (inc_null) any_null <- TRUE
    cat(sprintf("      %-24s RD = %7.3f   CI [%7.3f, %7.3f]\n",
                label, r$estimate, r$lower, r$upper))
    bar <- .ccr_contrast_bar(r$lower, r$upper, r$estimate, 0)
    cat("           ", bar, if (inc_null) "   *" else "", "\n\n", sep = "")
  }
  if (any_null) {
    cat("      * CI does not exclude the null (RD = 0) at the chosen alpha\n",
        "        (conditional on identifying assumptions - ",
        "see assumptions(fit))\n\n", sep = "")
  }
}

#' Proportional CI bar for an RD contrast (linear scale, null = 0)
#'
#' `[` / `]` mark the CI bounds, `|` the null, `*` the point estimate
#' (`X` when estimate = null). Display range spans
#' `[min(lower, 0), max(upper, 0)]` so the null is always visible.
#' Adapted from CausalSurvival's `.build_contrast_bar` (difference scale).
#' @keywords internal
.ccr_contrast_bar <- function(lower, upper, estimate, null_val, width = 14L) {
  if (!all(is.finite(c(lower, upper, estimate, null_val))) ||
      lower > estimate || estimate > upper) {
    return(sprintf("%7.3f   (non-finite or ill-ordered)   %7.3f",
                   lower, upper))
  }
  display_lo <- min(lower, null_val)
  display_hi <- max(upper, null_val)
  span       <- display_hi - display_lo
  if (span <= 0) {
    return(sprintf("%7.3f   (point CI = null)   %7.3f", lower, upper))
  }
  pos <- function(v) {
    max(1L, min(width, round((v - display_lo) / span * (width - 1L)) + 1L))
  }
  i_lo <- pos(lower); i_hi <- pos(upper)
  i_null <- pos(null_val); i_est <- pos(estimate)
  bar <- rep("-", width)
  bar[i_lo] <- "["; bar[i_hi] <- "]"
  if (i_null != i_lo && i_null != i_hi) bar[i_null] <- "|"
  if (i_est  != i_lo && i_est  != i_hi) bar[i_est] <- if (i_est == i_null) "X" else "*"
  sprintf("%7.3f   %s   %7.3f", display_lo, paste(bar, collapse = ""),
          display_hi)
}

#' @keywords internal
.ccr_summary_footer <- function(ci) {
  cat("  ", strrep("-", 66), "\n", sep = "")
  cat("           identification -> assumptions(fit)\n")
  if (!is.null(ci)) {
    cat(sprintf(
      "           bootstrap detail -> %d replicates @ %.0f%%\n",
      ci$n_boot_effective, (1 - ci$alpha) * 100))
  } else {
    cat("           for CIs: boot <- bootstrap(fit, n_boot = 500); ",
        "summary(fit, ci = boot)\n", sep = "")
  }
  cat("  ", strrep("-", 66), "\n")
}

#' @keywords internal
.ccr_summary_model_checks <- function(model_checks) {
  if (is.null(model_checks)) return(invisible())
  non_converged <- 0L
  for (chk in model_checks) {
    if (is.null(chk)) next
    if (!isTRUE(chk$converged)) non_converged <- non_converged + 1L
  }
  if (non_converged > 0L) {
    cat(sprintf(
      "\n  ! Model checks: %d non-converged model(s). See `fit$model_checks`.\n",
      non_converged))
  }
  invisible()
}


#' Confidence Intervals for causal_competing_risks
#'
#' Deprecated pathway. In the current design, confidence intervals are
#' computed by pairing a fit with a separate [bootstrap()] object.
#'
#' @param object A `"causal_competing_risks_fit"` object.
#' @param parm Not used (included for S3 consistency).
#' @param level Not used.
#' @param ... Additional arguments (currently unused).
#'
#' @return This method always errors — use the `bootstrap()` + `contrast()`
#'   pattern instead.
#' @export
confint.causal_competing_risks_fit <- function(object, parm = NULL, level = 0.95, ...) {
  stop(
    "Confidence intervals are not stored inside `fit` in this package. ",
    "Compute them explicitly:\n",
    "  boot <- bootstrap(fit, n_boot = 500)\n",
    "  contrast(fit, method = '<name>', ci = boot)",
    call. = FALSE
  )
}


#' Print a causal_competing_risks_risk Object
#'
#' Shows per-method per-arm cumulative incidence at the final time point,
#' indicates whether bootstrap CI bands are attached, and points toward
#' `plot()`.
#'
#' @param x A `"causal_competing_risks_risk"` object from [risk()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_competing_risks_risk <- function(x, ...) {
  title <- "Cumulative incidence curves (causal_competing_risks_risk)"
  cat(title, "\n", sep = "")
  cat(strrep("-", nchar(title)), "\n", sep = "")

  methods_avail <- unique(x$risk$method)
  cat("Methods: ", paste(methods_avail, collapse = ", "), "\n", sep = "")
  cat("Bootstrap CIs: ",
      if (any(!is.na(x$risk$lower))) "yes" else "no", "\n\n", sep = "")

  for (m in methods_avail) {
    sub   <- x$risk[x$risk$method == m, , drop = FALSE]
    K_max <- max(sub$k)
    last  <- sub[sub$k == K_max, , drop = FALSE]
    cat(sprintf("[%s] at final time (k = %g):\n", m, K_max))
    for (i in seq_len(nrow(last))) {
      band <- if (is.na(last$lower[i])) "" else
        sprintf("  [%.4f, %.4f]", last$lower[i], last$upper[i])
      cat(sprintf("  (%d,%d) %s = %.4f%s\n",
                  last$a_y[i], last$a_d[i], last$arm[i],
                  last$value[i], band))
    }
    cat("\n")
  }

  cat("Use plot(risk(fit), method = '<name>') to visualize.\n")
  invisible(x)
}


#' Print a causal_competing_risks_contrast Object
#'
#' Shows the contrast table at the selected time point, plus the method
#' and significance level.
#'
#' @param x A `"causal_competing_risks_contrast"` object from [contrast()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_competing_risks_contrast <- function(x, ...) {
  title <- "Causal contrasts (causal_competing_risks_contrast)"
  cat(title, "\n", sep = "")
  cat(strrep("-", nchar(title)), "\n", sep = "")
  cat("Method: ", x$method, "\n", sep = "")
  cat(sprintf("Significance level: %g (%.0f%% CIs)\n\n",
              x$alpha, (1 - x$alpha) * 100))

  k_at <- x$time %||% max(x$contrasts$k)
  cat(sprintf("At k = %g:\n", k_at))
  # Reference value under the null of no causal effect: rd -> 0, rr -> 1.
  out <- x$contrasts[, c("contrast", "decomp", "measure",
                         "estimate", "lower", "upper")]
  out$null_value <- ifelse(out$measure == "rd", 0, 1)
  out <- out[, c("contrast", "decomp", "measure", "null_value",
                 "estimate", "lower", "upper")]
  print(out, row.names = FALSE)

  cat("\nReading: if the [lower, upper] interval includes the `null_value` ",
      "(0 for rd, 1 for rr), the CI does not exclude the null at the chosen ",
      "alpha level (inconclusive). If the interval excludes the `null_value`, ",
      "the CI excludes the null. Both conclusions are conditional on ",
      "identifying assumptions - see assumptions(fit).\n", sep = "")
  cat(sprintf(
    "\n%d rows (1 time point x 10 contrasts). Pass `time = ...` to ",
    nrow(x$contrasts)
  ))
  cat("contrast() for a different cut time.\n")
  invisible(x)
}


#' Print a causal_competing_risks_bootstrap Object
#'
#' Replicate count (requested vs effective), significance level,
#' failed-replicate count, available methods, and pointers to the
#' accessors that pair the bands with a fit.
#'
#' @param x A `"causal_competing_risks_bootstrap"` object from [bootstrap()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_competing_risks_bootstrap <- function(x, ...) {
  title <- "Bootstrap confidence bands (causal_competing_risks_bootstrap)"
  cat(title, "\n", sep = "")
  cat(strrep("-", nchar(title)), "\n", sep = "")
  cat("Replicates requested: ", x$n_boot_requested, "\n", sep = "")
  cat("Replicates effective: ", x$n_boot_effective, "\n", sep = "")
  if (length(x$failed_reps) > 0L) {
    cat("Failed replicates: ", length(x$failed_reps),
        " (see $failed_reps for indices)\n", sep = "")
  }
  cat(sprintf("Significance level: %g (%.0f%% CIs)\n",
              x$alpha, (1 - x$alpha) * 100))
  cat("Methods: ", paste(unique(x$bands$method), collapse = ", "), "\n",
      sep = "")
  cat("\nUse `contrast(fit, method = '<name>', ci = <this>)` for contrast CIs,\n")
  cat("or `plot(risk(fit, ci = <this>), method = '<name>')` for curves with bands.\n")
  invisible(x)
}


#' Print a causal_competing_risks_diagnostic Object
#'
#' Shows per-model fit diagnostics (convergence, fitted probability range,
#' positivity violation flag) and IPW weight summary if available. Uses
#' the `flagged_ids` / `flagged_log` fields (both truncation and trimming
#' cases).
#'
#' @param x A `"causal_competing_risks_diagnostic"` object from [diagnostic()].
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns `x`.
#' @export
print.causal_competing_risks_diagnostic <- function(x, ...) {
  cat("Diagnostics (causal_competing_risks_diagnostic)\n")
  cat("-----------------------------------\n")

  # --- Model checks ---
  if (!is.null(x$model_checks)) {
    cat("\nModel checks:\n")
    any_printed <- FALSE
    for (m in names(x$model_checks)) {
      chk <- x$model_checks[[m]]
      if (is.null(chk)) next
      any_printed <- TRUE
      cat(sprintf(
        "  [%s]  converged=%s  fitted in [%.3g, %.3g]\n",
        chk$label %||% m,
        if (isTRUE(chk$converged)) "TRUE" else "FALSE",
        chk$min_fitted, chk$max_fitted
      ))
      if (length(chk$glm_warnings) > 0) {
        cat("    glm warnings:\n")
        for (w in chk$glm_warnings) {
          cat("      - ", w, "\n", sep = "")
        }
      }
    }
    if (!any_printed) cat("  (no fitted models)\n")
  }

  # --- Weight summary (IPW only) ---
  if (!is.null(x$weight_summary)) {
    cat("\nWeight summary:\n")
    print(x$weight_summary, row.names = FALSE)

    n_flagged <- length(x$flagged_ids)
    n_log_rows <- if (!is.null(x$flagged_log)) nrow(x$flagged_log) else 0
    trunc_label <- if (is.null(x$truncate)) "none" else
      sprintf("[%g, %g]", x$truncate[1], x$truncate[2])
    cat(sprintf(
      "\nFlagged subjects (truncate=%s): %d  (%d row(s) in flagged_log)\n",
      trunc_label, n_flagged, n_log_rows
    ))
  } else {
    cat("\nNo IPW weight diagnostics (no IPW method ran).\n")
  }

  invisible(x)
}
