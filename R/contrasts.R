#' Compute Causal Contrasts (Long Format with Confidence Intervals)
#'
#' Builds the long-format contrast table for a single estimation method by
#' combining the method's point estimates (from the fit) with bootstrap
#' quantiles (from a `"causal_competing_risks_bootstrap"` object). Used internally by
#' [contrast()].
#'
#' @param fit A `"causal_competing_risks_fit"` object.
#' @param method Character (length 1). Which method to extract from
#'   `fit$cumulative_incidence` and `ci$replicates`. One of the names in
#'   `fit$active_methods`.
#' @param ci A `"causal_competing_risks_bootstrap"` object from [bootstrap()].
#'
#' @return A long-format data.frame with columns:
#'   \describe{
#'     \item{k}{Time point.}
#'     \item{contrast}{One of `"total"`, `"sep_direct"`, `"sep_indirect"`.}
#'     \item{decomp}{`NA` for totals, `"A"` or `"B"` for direct/indirect.}
#'     \item{measure}{`"rd"` or `"rr"`.}
#'     \item{estimate}{Point estimate.}
#'     \item{lower,upper}{Percentile confidence interval bounds from the
#'       bootstrap distribution (`alpha/2` and `1 - alpha/2`).}
#'   }
#'
#' Rows per time point: 10 (1 total + 2 direct + 2 indirect, each for `rd`
#' and `rr`). For methods that do not emit `arm_01` (IPW Rep 1 and Rep 2),
#' Decomposition B rows have `NA` estimates and CIs.
#'
#' @details
#' Both decompositions share the same reference arm (`arm_00`) and the same
#' total (`arm_11 - arm_00`). They differ only in the **intermediate arm** —
#' the path order from `arm_00` to `arm_11`. Convention: arm `"xy"` means
#' `(a_Y = x, a_D = y)`. `A_Y` acts on the path A → Y not through D (its
#' effect is the *direct* effect, SDE). `A_D` acts on the path A → D → Y
#' (its effect is the *indirect* effect, SIE).
#'
#' ## Decomposition A — vary `A_Y` first, then `A_D` (intermediate arm `arm_10`)
#' - `direct  (SDE-A) = arm_10 - arm_00`  (effect of `A_Y` at `a_D = 0`)
#' - `indirect (SIE-A) = arm_11 - arm_10` (effect of `A_D` at `a_Y = 1`)
#'
#' ## Decomposition B — vary `A_D` first, then `A_Y` (intermediate arm `arm_01`)
#' - `indirect (SIE-B) = arm_01 - arm_00` (effect of `A_D` at `a_Y = 0`)
#' - `direct  (SDE-B) = arm_11 - arm_01`  (effect of `A_Y` at `a_D = 1`)
#'
#' Total = `arm_11 - arm_00` is identical under both decompositions, and the
#' algebraic identity `SDE + SIE = Total` holds for each.
#'
#' @references
#' Stensrud MJ, Young JG, Didelez V, Robins JM, Hernán MA (2020).
#' Separable effects for causal inference in the presence of competing events.
#' *Journal of the American Statistical Association*, §3.
#' \doi{10.1080/01621459.2020.1765783}
#'
#' @family internal
#' @keywords internal
compute_contrasts <- function(fit, method, ci) {

  stopifnot(inherits(fit, "causal_competing_risks_fit"))
  stopifnot(inherits(ci, "causal_competing_risks_bootstrap"))

  if (missing(method) || !is.character(method) || length(method) != 1) {
    stop("'method' must be a single method name.", call. = FALSE)
  }
  if (!method %in% names(fit$cumulative_incidence)) {
    stop(
      "Method '", method, "' not in fit. Available: ",
      paste(names(fit$cumulative_incidence), collapse = ", "),
      call. = FALSE
    )
  }
  if (!method %in% unique(ci$replicates$method)) {
    stop(
      "Method '", method, "' not in bootstrap replicates. Available: ",
      paste(unique(ci$replicates$method), collapse = ", "),
      call. = FALSE
    )
  }

  alpha <- ci$alpha
  lo <- alpha / 2
  hi <- 1 - alpha / 2

  ci_df <- fit$cumulative_incidence[[method]]
  k  <- ci_df$k

  # Point estimates from fit
  a11 <- ci_df$arm_11
  a00 <- ci_df$arm_00
  a10 <- ci_df$arm_10
  a01 <- if ("arm_01" %in% names(ci_df)) ci_df$arm_01 else rep(NA_real_, length(k))

  # Bootstrap replicates -> per-arm [n_boot x n_times] matrices for method
  mats <- boot_arm_matrices(ci$replicates, method, arms = arm_spec()$name)
  b11 <- mats[["arm_11"]]
  b00 <- mats[["arm_00"]]
  b10 <- mats[["arm_10"]]
  b01 <- mats[["arm_01"]]

  boot_total_rd      <- b11 - b00
  boot_total_rr      <- b11 / b00
  boot_direct_A_rd   <- b10 - b00
  boot_direct_A_rr   <- b10 / b00
  boot_indirect_A_rd <- b11 - b10
  boot_indirect_A_rr <- b11 / b10
  boot_direct_B_rd   <- b11 - b01
  boot_direct_B_rr   <- b11 / b01
  boot_indirect_B_rd <- b01 - b00
  boot_indirect_B_rr <- b01 / b00

  # Quantile helper over n_boot for each time point
  q_at_k <- function(boot_mat, p) {
    apply(boot_mat, 2, quantile, probs = p, na.rm = TRUE)
  }

  # Build rows for one contrast entry
  make_rows <- function(contrast_name, decomp_label, measure_label,
                        point_est, boot_mat) {
    data.frame(
      k        = k,
      contrast = contrast_name,
      decomp   = decomp_label,
      measure  = measure_label,
      estimate = point_est,
      lower    = unname(q_at_k(boot_mat, lo)),
      upper    = unname(q_at_k(boot_mat, hi)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  rbind(
    # --- RD super-block ---
    make_rows("total",    NA_character_, "rd", a11 - a00, boot_total_rd),
    make_rows("sep_direct",   "A",           "rd", a10 - a00, boot_direct_A_rd),
    make_rows("sep_indirect", "A",           "rd", a11 - a10, boot_indirect_A_rd),
    make_rows("sep_direct",   "B",           "rd", a11 - a01, boot_direct_B_rd),
    make_rows("sep_indirect", "B",           "rd", a01 - a00, boot_indirect_B_rd),

    # --- RR super-block ---
    make_rows("total",    NA_character_, "rr", a11 / a00, boot_total_rr),
    make_rows("sep_direct",   "A",           "rr", a10 / a00, boot_direct_A_rr),
    make_rows("sep_indirect", "A",           "rr", a11 / a10, boot_indirect_A_rr),
    make_rows("sep_direct",   "B",           "rr", a11 / a01, boot_direct_B_rr),
    make_rows("sep_indirect", "B",           "rr", a01 / a00, boot_indirect_B_rr)
  )
}
