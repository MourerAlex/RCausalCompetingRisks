#' Bootstrap Confidence Intervals for a causal_competing_risks Fit
#'
#' Performs subject-level bootstrap resampling on a [causal_competing_risks()] fit and
#' constructs percentile confidence intervals for cumulative incidence
#' curves. Returns a separate `"causal_competing_risks_bootstrap"` object that pairs
#' with the fit — pass both to plotting and contrast functions to get
#' confidence bands.
#'
#' @param fit A `"causal_competing_risks_fit"` object from [causal_competing_risks()].
#' @param n_boot Integer. Number of bootstrap replicates (default 500).
#' @param alpha Numeric. Two-sided significance level (default 0.05 for
#'   95% CIs).
#'
#' @param seed Optional integer RNG seed for reproducibility (set once
#'   before the replicate loop).
#' @param verbose Logical. Print progress during the replicate loop
#'   (default `TRUE`).
#'
#' @return An S3 object of class `"causal_competing_risks_bootstrap"`:
#'   \describe{
#'     \item{replicates}{Long data.frame `(boot_id, method, arm, k, value)`
#'       of cumulative incidence estimates per bootstrap draw.}
#'     \item{bands}{Long data.frame `(method, arm, k, lower, upper)`: the
#'       `alpha/2` and `1 - alpha/2` percentile bands per (method, arm, k).}
#'     \item{n_boot_requested, n_boot_effective}{Replicates requested and
#'       the number that survived (failures dropped).}
#'     \item{alpha}{The significance level used.}
#'     \item{failed_reps}{Integer indices of replicates whose refit errored.}
#'     \item{fit_call}{Copy of `fit$call` for provenance.}
#'   }
#'
#' @details
#' ## Resampling scheme
#' Subject-level resampling with replacement. For each replicate, unique
#' IDs are sampled with replacement; the person-time rows for each sampled
#' ID are pulled (via a pre-split lookup). Duplicate draws are given unique
#' synthetic IDs so each resampled subject is treated as distinct.
#'
#' ## Estimation per replicate
#' Calls [fit_separable_effects()] directly (not [causal_competing_risks()]), bypassing the
#' user-facing wrapper's validation and warning capture. Warnings during
#' replicates are suppressed (they would otherwise spam).
#'
#' ## CI construction
#' Percentile method: `alpha/2` and `1 - alpha/2` quantiles of the bootstrap
#' distribution at each time point, computed per arm per method.
#'
#' ## Progress reporting
#' For n_boot > 50: prints every 10 replicates for the first 50, then
#' prints a time estimate for the remaining replicates (based on the first
#' 50), then every 100 after. For n_boot <= 50: prints every 10.
#'
#' @seealso [causal_competing_risks()], [fit_separable_effects()], [contrast()], [risk()]
#'
#' @examples
#' \dontrun{
#' fit <- causal_competing_risks(pt)
#' boot <- bootstrap(fit, n_boot = 500)
#' plot(risk(fit), ci = boot)
#' contrast(fit, method = "gformula", ci = boot)
#' }
#'
#' @param ... Additional arguments (currently unused).
#'
#' @importFrom CausalSurvival bootstrap_init_state bootstrap_progress_reporter
#' @importFrom CausalSurvival bootstrap_resample_pt bootstrap_percentile_bands
#' @importFrom CausalSurvival bootstrap
#' @export
bootstrap.causal_competing_risks_fit <- function(fit, n_boot = 500,
                                                 alpha = 0.05,
                                                 seed = NULL,
                                                 verbose = TRUE, ...) {

  # 1. Validate / canonicalize args (class, n_boot, alpha, seed).
  args <- .validate_bootstrap_args(fit, n_boot, alpha, seed)

  # 2. Resampling state: split person-time by subject id once.
  state <- bootstrap_init_state(args$pt_data, args$id_col)

  # 3. Replicate loop: re-eval the fit's wrapper call on each bootstrap
  #    sample, record per-(method, arm) cumulative incidence in long form,
  #    or note the failure.
  reps_long       <- list()
  failed_reps     <- integer()
  report_progress <- bootstrap_progress_reporter(args$n_boot)
  for (b in seq_len(args$n_boot)) {
    if (verbose) report_progress(b)
    boot_data <- bootstrap_resample_pt(state, args$pt_data, args$id_col)
    call_b <- args$fit$call
    call_b$pt_data <- boot_data
    res <- tryCatch(suppressWarnings(eval(call_b)), error = function(e) NULL)
    if (is.null(res)) {
      failed_reps <- c(failed_reps, b)
      next
    }
    reps_long[[length(reps_long) + 1L]] <-
      .boot_replicate_rows(res, b, args$arm_names)
  }

  # 4. Aggregate per-replicate CIFs into the long replicates table.
  replicates <- if (length(reps_long) == 0L) {
    data.frame(boot_id = integer(), method = character(), arm = character(),
               k = integer(), value = numeric(), stringsAsFactors = FALSE)
  } else do.call(rbind, reps_long)

  # 5. Percentile bands per (method, arm, k).
  bands <- bootstrap_percentile_bands(replicates, args$alpha,
                                      by = c("method", "arm", "k"))

  # 6. Assemble S3 output.
  structure(
    list(
      fit_call         = args$fit$call,
      n_boot_requested = as.integer(args$n_boot),
      n_boot_effective = as.integer(args$n_boot - length(failed_reps)),
      alpha            = args$alpha,
      replicates       = replicates,
      bands            = bands,
      failed_reps      = failed_reps
    ),
    class = "causal_competing_risks_bootstrap"
  )
}


#' Validate and canonicalize bootstrap() arguments
#'
#' Checks `fit` class, integer-positive `n_boot`, `(0, 1)` `alpha`, sets
#' the RNG seed when supplied, and verifies `fit$person_time` is present.
#' Returns a canonical args list.
#'
#' @keywords internal
.validate_bootstrap_args <- function(fit, n_boot, alpha, seed) {
  stopifnot(inherits(fit, "causal_competing_risks_fit"))
  if (!is.numeric(n_boot) || length(n_boot) != 1L ||
      n_boot < 1 || n_boot != round(n_boot)) {
    stop("`n_boot` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be in (0, 1).", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  if (is.null(fit$person_time)) {
    stop("`fit$person_time` is NULL; cannot bootstrap.", call. = FALSE)
  }
  list(
    fit       = fit,
    n_boot    = n_boot,
    alpha     = alpha,
    pt_data   = fit$person_time,
    id_col    = fit$id_col,
    arm_names = arm_spec()$name
  )
}


#' Long rows for one surviving bootstrap replicate
#'
#' Flattens a replicate's per-method cumulative incidence into long rows
#' `(boot_id, method, arm, k, value)`, one block per (method, arm) the
#' method emitted.
#'
#' @keywords internal
.boot_replicate_rows <- function(res, b, arm_names) {
  rows <- list()
  for (m in names(res$cumulative_incidence)) {
    ci_df <- res$cumulative_incidence[[m]]
    if (is.null(ci_df)) next
    for (arm in intersect(arm_names, names(ci_df))) {
      rows[[length(rows) + 1L]] <- data.frame(
        boot_id = b, method = m, arm = arm,
        k = ci_df$k, value = ci_df[[arm]],
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }
  do.call(rbind, rows)
}
