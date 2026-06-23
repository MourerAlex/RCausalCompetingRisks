#' Plot Cumulative Incidence Curves
#'
#' Plots cumulative incidence for each arm of one estimation method. Pairs
#' with [risk()] — pass a `"causal_competing_risks_risk"` object. Optional bootstrap CI
#' ribbons when the paired `risk()` was built with `ci = boot`. Optional
#' contrast annotations (dotted vertical bridges + numeric labels) at
#' user-specified time points, and an optional risk table stacked below.
#'
#' @param x A `"causal_competing_risks_risk"` object from [risk()].
#' @param method Character (length 1). Required. One of
#'   `names(x$cumulative_incidence)`. No default.
#' @param arms Character vector. Which arms to draw. Defaults to the 3 core
#'   arms (`arm_11`, `arm_00`, `arm_10`); include `"arm_01"` explicitly for
#'   Decomposition B sensitivity.
#' @param eval_times Numeric vector or NULL. Time points at which to draw
#'   annotations. Only has an effect when `contrast_annotations` is also
#'   provided.
#' @param contrast_annotations Character vector. Max 2. Which contrasts to
#'   annotate with dotted vertical bridges + numeric labels at `eval_times`.
#'   One or two of:
#'   `"total"`, `"sep_direct_A"`, `"sep_indirect_A"`, `"sep_direct_B"`,
#'   `"sep_indirect_B"`.
#' @param curves Logical. `TRUE` (default) draws the cumulative-incidence
#'   curves panel; `FALSE` suppresses it so only the requested `risk_table`
#'   panel(s) render (table-only mode). Requires a non-NULL `risk_table`.
#' @param risk_table NULL (default, no table) or a character vector with
#'   entries in `"at_risk"`, `"events_y"`, `"events_d"`, `"censored"`. One
#'   counts panel per entry is stacked below the curves via [/]. The
#'   `risk()` accessor stores the needed references — no need to re-pass
#'   `fit`.
#' @param risk_table_height Numeric. Height of the risk table relative to
#'   the main plot (which is 1). Default `0.23` (the table is ~23% of the
#'   main plot's height). Increase to give the table more room.
#' @param arm_colors Named character vector overriding arm hex colors.
#'   Defaults to Okabe-Ito. Only the arms you pass are overridden; others
#'   keep the default.
#' @param arm_labels Named character vector overriding arm legend labels.
#' @param title Character or NULL. Plot title. Default:
#'   `"Cumulative incidence (<method>)"`.
#' @param subtitle Character or NULL. Plot subtitle. Default: none.
#' @param x_label Character. X-axis label. Default `"Time (interval k)"`.
#' @param y_label Character. Y-axis label. Default `"Cumulative incidence"`.
#' @param base_size Numeric. Base font size for `theme_minimal()`. Scales
#'   all text (title, axes, legend). Default 11.
#' @param annotation_size Numeric. Size of contrast annotation label text
#'   (passed to `geom_label(size = ...)`, which uses ggplot's mm units).
#'   Default 2.8.
#' @param linewidth Numeric. Width of the cumulative incidence step lines.
#'   Default 0.8.
#' @param ribbon_alpha Numeric in [0, 1]. Transparency of the CI ribbons.
#'   Default 0.15.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object (or a patchwork object when `risk_table` is set).
#'
#' @details
#' Default colors follow the Okabe-Ito colorblind-safe palette:
#' - `arm_11` (treated)              black    `#000000`
#' - `arm_00` (control)              blue     `#0072B2`
#' - `arm_10` (separable, via A_Y=1) green    `#009E73`
#' - `arm_01` (separable, via A_D=1) vermillion `#D55E00`
#'
#' Override individual entries via `arm_colors`. Non-overridden arms keep
#' their defaults.
#'
#' @family plot
#' @export
plot.causal_competing_risks_risk <- function(x,
                                method,
                                arms = NULL,
                                eval_times = NULL,
                                contrast_annotations = NULL,
                                curves = TRUE,
                                risk_table = NULL,
                                risk_table_height = 0.23,
                                arm_colors = NULL,
                                arm_labels = NULL,
                                title = NULL,
                                subtitle = NULL,
                                x_label = "Time (interval k)",
                                y_label = "Cumulative incidence",
                                base_size = 11,
                                annotation_size = 2.8,
                                linewidth = 0.8,
                                ribbon_alpha = 0.15,
                                ...) {

  # --- Validate method (required, no default) ---
  avail_methods <- unique(x$risk$method)
  if (missing(method)) {
    stop(
      "'method' is required. Available: ",
      paste(avail_methods, collapse = ", "),
      call. = FALSE
    )
  }
  if (!method %in% avail_methods) {
    stop(
      "'method' must be one of: ", paste(avail_methods, collapse = ", "),
      call. = FALSE
    )
  }

  # Long-format slice for this method: (arm, a_y, a_d, k, value, lower, upper)
  risk_method <- x$risk[x$risk$method == method, ]
  have_bands  <- any(!is.na(risk_method$lower)) && any(!is.na(risk_method$upper))

  # Per-arm [n_boot x n_times] replicate matrices for this method, if
  # bootstrap was supplied. Used to compute PROPER contrast CIs at eval_times.
  boot_method <- if (!is.null(x$replicates) && method %in% unique(x$replicates$method)) {
    boot_arm_matrices(x$replicates, method, arms = arm_spec()$name)
  } else NULL
  alpha_val <- x$alpha  # NULL if no bootstrap

  # --- Validate arms ---
  spec       <- arm_spec()
  all_arms   <- spec$name
  avail_arms <- unique(risk_method$arm)
  if (is.null(arms)) {
    # Default: 3 core arms. arm_01 is a Decomposition B sensitivity — cluttery
    # to show by default alongside the main three. User includes it
    # explicitly when they want the 4th curve.
    arms <- intersect(setdiff(all_arms, "arm_01"), avail_arms)
  } else {
    bad <- setdiff(arms, all_arms)
    if (length(bad) > 0) {
      stop("Unknown arm(s): ", paste(bad, collapse = ", "),
           ". Must be subset of: ", paste(all_arms, collapse = ", "),
           call. = FALSE)
    }
    missing_in_data <- setdiff(arms, avail_arms)
    if (length(missing_in_data) > 0) {
      stop("Arm(s) not present in this method's output: ",
           paste(missing_in_data, collapse = ", "),
           call. = FALSE)
    }
  }

  # --- Validate contrast_annotations ---
  valid_annotations <- c("total", "sep_direct_A", "sep_indirect_A",
                         "sep_direct_B", "sep_indirect_B")
  if (!is.null(contrast_annotations)) {
    bad <- setdiff(contrast_annotations, valid_annotations)
    if (length(bad) > 0) {
      stop("Unknown contrast_annotations: ", paste(bad, collapse = ", "),
           ". Must be subset of: ", paste(valid_annotations, collapse = ", "),
           call. = FALSE)
    }
    if (length(contrast_annotations) > 2) {
      stop("contrast_annotations is limited to 2 entries (more gets cluttered).",
           call. = FALSE)
    }
    if (is.null(eval_times)) {
      stop("contrast_annotations requires eval_times to be specified.",
           call. = FALSE)
    }
  }

  # --- Arm palette (Okabe-Ito default, user-overridable per-arm) ---
  default_arm_colors <- setNames(spec$color, spec$name)
  default_arm_labels <- setNames(spec$label, spec$name)
  if (!is.null(arm_colors)) {
    bad_c <- setdiff(names(arm_colors), names(default_arm_colors))
    if (length(bad_c) > 0) {
      stop("arm_colors has unknown entries: ", paste(bad_c, collapse = ", "),
           call. = FALSE)
    }
    # Overlay user values onto defaults
    default_arm_colors[names(arm_colors)] <- arm_colors
  }
  if (!is.null(arm_labels)) {
    bad_l <- setdiff(names(arm_labels), names(default_arm_labels))
    if (length(bad_l) > 0) {
      stop("arm_labels has unknown entries: ", paste(bad_l, collapse = ", "),
           call. = FALSE)
    }
    default_arm_labels[names(arm_labels)] <- arm_labels
  }
  arm_colors <- default_arm_colors
  arm_labels <- default_arm_labels

  # --- Shared x-axis ticks for main plot AND risk table ---
  # Computed once here so both panels use IDENTICAL breaks + limits, which
  # patchwork then aligns to the same column positions. Endpoints (0 and
  # max(cut_times)) are always included as ticks, with pretty() filling in
  # between.
  k_grid <- sort(unique(risk_method$k))
  cut_times_full <- c(0, k_grid)
  shared_x_limits <- c(0, max(cut_times_full))
  pretty_inner <- pretty(shared_x_limits, n = 5)
  pretty_inner <- pretty_inner[
    pretty_inner > 0 & pretty_inner < shared_x_limits[2]
  ]
  shared_x_ticks <- sort(unique(c(shared_x_limits, pretty_inner)))

  # --- Build plot data from long-format $risk slice ---
  # Prepend a (k = 0, cum_inc = 0) row per arm so curves start at the
  # origin and the x-axis can extend to 0, matching the risk-table x-axis
  # (which always starts at 0).
  body_rows <- risk_method[risk_method$arm %in% arms, c("k", "arm", "value")]
  names(body_rows)[names(body_rows) == "value"] <- "cum_inc"
  origin_rows <- data.frame(
    k       = 0,
    arm     = arms,
    cum_inc = 0,
    stringsAsFactors = FALSE
  )
  plot_data <- rbind(origin_rows, body_rows)
  plot_data$arm_label <- arm_labels[plot_data$arm]

  # --- Base plot ---
  p <- ggplot(plot_data, aes(
    x = .data$k, y = .data$cum_inc,
    color = .data$arm_label, group = .data$arm_label
  )) +
    geom_step(linewidth = linewidth) +
    scale_color_manual(
      values = setNames(arm_colors[arms], arm_labels[arms])
    ) +
    # Shared x-scale (matches the risk table when present)
    scale_x_continuous(
      breaks = shared_x_ticks,
      limits = shared_x_limits,
      expand = expansion(mult = 0.02)
    ) +
    labs(
      x     = x_label,
      y     = y_label,
      color = "Arm",
      title = title %||% paste0("Cumulative incidence (", method, ")"),
      subtitle = subtitle
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = "bottom",
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(face = "plain"),
      axis.title.x  = element_text(face = "bold"),
      axis.title.y  = element_text(face = "bold"),
      legend.title  = element_text(face = "bold")
    )

  # --- CI ribbons when bootstrap available (step-like to match geom_step) ---
  if (have_bands) {
    # Prepend (k = 0, lower = 0, upper = 0) so ribbons start at origin
    # alongside the curves.
    body_bands <- risk_method[risk_method$arm %in% arms,
                              c("k", "arm", "lower", "upper")]
    origin_bands <- data.frame(
      k     = 0,
      arm   = arms,
      lower = 0,
      upper = 0,
      stringsAsFactors = FALSE
    )
    ribbon_data <- rbind(origin_bands, body_bands)

    # Step-transform: geom_ribbon draws linear polygons between (x, ymin, ymax)
    # points. To match the stepped cumulative incidence curves, we duplicate
    # points so the ribbon stays flat from k[i] to k[i+1] then "jumps" at
    # k[i+1]. For each arm with n rows:
    #   k_step     = c(k[1], k[2], k[2], k[3], k[3], ..., k[n])
    #   lower_step = c(l[1], l[1], l[2], l[2], ...,            l[n])
    #   upper_step = c(u[1], u[1], u[2], u[2], ...,            u[n])
    ribbon_data <- do.call(rbind, lapply(arms, function(a) {
      d <- ribbon_data[ribbon_data$arm == a, ]
      d <- d[order(d$k), ]
      n <- nrow(d)
      if (n < 2) return(d)
      k_step     <- c(d$k[1],     rep(d$k[2:n], each = 2))
      lower_step <- c(rep(d$lower[1:(n - 1)], each = 2), d$lower[n])
      upper_step <- c(rep(d$upper[1:(n - 1)], each = 2), d$upper[n])
      data.frame(
        k     = k_step,
        lower = lower_step,
        upper = upper_step,
        arm   = a,
        stringsAsFactors = FALSE
      )
    }))
    ribbon_data$arm_label <- arm_labels[ribbon_data$arm]

    p <- p + geom_ribbon(
      data = ribbon_data,
      aes(
        x = .data$k, ymin = .data$lower, ymax = .data$upper,
        fill = .data$arm_label
      ),
      alpha = ribbon_alpha, inherit.aes = FALSE
    ) +
      scale_fill_manual(
        values = setNames(arm_colors[arms], arm_labels[arms]),
        guide = "none"
      )
  }

  # --- Contrast annotations at eval_times ---
  if (!is.null(contrast_annotations) && length(eval_times) > 0) {
    annotation_df <- build_contrast_annotations(
      risk_method, have_bands, boot_method,
      contrast_annotations, eval_times, alpha_val
    )

    # Dotted vertical bridges between each pair at each eval time
    p <- p + geom_segment(
      data = annotation_df,
      aes(
        x = .data$k, xend = .data$k,
        y = .data$y_low, yend = .data$y_high
      ),
      linetype = "dotted", color = "grey30", inherit.aes = FALSE
    ) +
      # Label sits below the lowest point of interest at that eval_time (which
      # is either the lower CI band of the lower arm, or the lower arm curve
      # itself if no CI). Keeps the bridge visually clean.
      geom_label(
        data = annotation_df,
        aes(
          x = .data$k,
          y = .data$y_label,
          label = .data$label
        ),
        size = annotation_size, inherit.aes = FALSE, label.size = 0.15,
        vjust = 1
      )
  }

  # --- Risk table panel(s) (one per requested count, stacked via patchwork) ---
  if (!is.null(risk_table)) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop(
        "risk_table requires the 'patchwork' package. Install via ",
        "install.packages('patchwork').",
        call. = FALSE
      )
    }
    if (is.null(x$person_time) || is.null(x$id_col) ||
        is.null(x$treatment_col) || is.null(x$times)) {
      stop(
        "risk_table needs person-time data on the risk() object. ",
        "Did you build it via risk(fit) from a 'causal_competing_risks' fit?",
        call. = FALSE
      )
    }
    valid_counts <- c("at_risk", "events_y", "events_d", "censored")
    if (!is.character(risk_table)) {
      stop("`risk_table` must be NULL or a character vector with entries in: ",
           paste(shQuote(valid_counts), collapse = ", "), ".",
           call. = FALSE)
    }
    bad_rt <- setdiff(risk_table, valid_counts)
    if (length(bad_rt) > 0L) {
      stop("`risk_table` has unknown entries: ",
           paste(shQuote(bad_rt), collapse = ", "),
           ". Valid: ", paste(shQuote(valid_counts), collapse = ", "), ".",
           call. = FALSE)
    }

    n_tbl <- length(risk_table)
    tbl_plots <- lapply(seq_along(risk_table), function(i) {
      build_risk_table_plot(
        pt_data   = x$person_time,
        id_col    = x$id_col,
        trt_col   = x$treatment_col,
        cut_times = x$times,
        count     = risk_table[[i]],
        base_size = base_size,
        x_breaks  = shared_x_ticks,
        x_limits  = shared_x_limits,
        # Bottom-most panel gets the cut-time tick labels; the ones above
        # stay clean (curves panel / upper tables carry no x-axis text).
        show_x_axis = (i == n_tbl)
      )
    })

    if (isTRUE(curves)) {
      # curves on top + risk-table panel(s) below. `risk_table_height` is
      # the per-panel ratio (curves panel = 1); total bottom = N * h.
      p <- wrap_plots(
        c(list(p), tbl_plots), ncol = 1,
        heights = c(1, rep(risk_table_height, length(tbl_plots)))
      )
    } else {
      # Table-only: drop curves panel, stack the table panels equally.
      p <- wrap_plots(
        tbl_plots, ncol = 1,
        heights = rep(1, length(tbl_plots))
      )
    }
  } else if (!isTRUE(curves)) {
    stop("`curves = FALSE` requires a non-NULL `risk_table` ",
         "(nothing else would render).", call. = FALSE)
  }

  p
}


#' Build Risk Table Plot for Stacking Below Cumulative Incidence Curves
#'
#' Renders the output of [risk_table()] as a minimalist ggplot (tile + text)
#' suitable for patchwork stacking below a curves plot.
#'
#' @param fit A "causal_competing_risks_fit" object.
#' @param count Character, one of the accepted values in [risk_table()].
#' @param base_size Base font size (inherited from the parent plot call).
#' @param show_x_axis Logical. Draw the x-axis tick labels/ticks. Set
#'   `FALSE` for upper panels when stacking multiple tables so only the
#'   bottom-most panel carries the shared axis.
#' @return A ggplot2 object.
#' @family internal
#' @keywords internal
build_risk_table_plot <- function(pt_data, id_col, trt_col,
                                  cut_times, count, base_size = 11,
                                  x_breaks = NULL, x_limits = NULL,
                                  show_x_axis = TRUE) {
  # Include k = 0 (baseline) as an explicit time point so the table covers
  # the full range from origin to max(cut_times).
  table_times <- c(0, cut_times)

  tbl <- risk_table_internal(pt_data, id_col, trt_col, table_times, count)
  arm_cols <- setdiff(names(tbl), "k")

  # --- Tick positions: use breaks supplied by the parent plot when given,
  # otherwise fall back to a pretty() default. The "supplied" path is what
  # plot.causal_competing_risks_risk uses to guarantee perfect alignment with the main
  # plot's x-axis.
  pretty_ticks <- if (!is.null(x_breaks)) {
    x_breaks
  } else {
    inner <- pretty(c(0, max(table_times)), n = 5)
    inner <- inner[inner > 0 & inner < max(table_times)]
    sort(unique(c(0, max(table_times), inner)))
  }
  scale_limits <- x_limits %||% c(0, max(table_times))

  # For each pretty tick, find the nearest table_times index and pull the
  # count there. We display the number AT THE PRETTY TICK position.
  snapped_idx <- vapply(pretty_ticks, function(b) {
    which.min(abs(table_times - b))
  }, integer(1))

  long <- do.call(rbind, lapply(arm_cols, function(ac) {
    data.frame(
      k     = pretty_ticks,
      arm   = ac,
      n     = tbl[[ac]][snapped_idx],
      stringsAsFactors = FALSE
    )
  }))
  long$arm <- factor(long$arm, levels = rev(arm_cols))

  title_map <- c(
    at_risk  = "Number at risk",
    events_y = "Number of events (Y)",
    events_d = "Number of events (D)",
    censored = "Number censored"
  )
  title_text <- title_map[count] %||% count

  ggplot(long, aes(x = .data$k, y = .data$arm)) +
    geom_text(
      aes(label = .data$n),
      size  = base_size * 0.32,
      hjust = 0.5
    ) +
    scale_x_continuous(
      breaks = pretty_ticks,
      labels = pretty_ticks,
      limits = scale_limits,
      expand = expansion(mult = 0.02)
    ) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(
      x     = NULL,
      y     = trt_col,
      title = title_text
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      # Strip everything but axis lines on bottom + left for the
      # adjustedCurves-style frame.
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.line.x.bottom = element_line(color = "black",
                                                 linewidth = 0.4),
      axis.line.y.left   = element_line(color = "black",
                                                 linewidth = 0.4),
      # Only the bottom-most stacked panel shows the x-axis tick labels;
      # upper panels stay clean so the shared axis reads once.
      axis.text.x      = if (show_x_axis) element_text() else element_blank(),
      axis.ticks.x     = if (show_x_axis) {
                           element_line(color = "black", linewidth = 0.3)
                         } else element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_text(face = "bold"),
      axis.title.y     = element_text(face = "bold", angle = 90),
      plot.title       = element_text(
                           face = "bold",
                           size = base_size * 0.95
                         ),
      plot.title.position = "plot"
    )
}


#' Build Contrast Annotation Rows for plot.causal_competing_risks_risk
#'
#' For each (annotation, eval_time) pair, compute the bridge endpoints
#' (y_low, y_high) from the two arms in the contrast, plus a numeric label
#' showing the point estimate. If bootstrap replicates (`boot_method`) are
#' supplied, computes a PROPER per-contrast percentile CI by taking
#' per-replicate differences then quantiles (NOT the naive difference of
#' per-arm CIs, which is too wide because it ignores cross-arm covariance
#' within each bootstrap replicate).
#'
#' Also computes `y_label` — the y-coordinate where the label should be
#' drawn. This sits below the lower CI band at that time point (or below
#' the lower arm curve if no CI bands), slightly offset for legibility.
#'
#' Emits a one-time `message()` if any `eval_times` had to be snapped to the
#' nearest `cut_times` (when they didn't match exactly).
#'
#' @param ci_df Data.frame with `k` and the arm columns (point estimates).
#' @param ci_bands Data.frame with `k` and `{arm}_lower` / `{arm}_upper`
#'   columns, or NULL if no bootstrap was supplied.
#' @param boot_method Named list of per-arm `[n_boot x n_times]` replicate
#'   matrices for the chosen method (from [boot_arm_matrices()]), or NULL
#'   if no bootstrap.
#' @param contrast_annotations Character vector (validated upstream).
#' @param eval_times Numeric vector of time points (user-specified, may
#'   not match `cut_times` exactly).
#' @param alpha_val Numeric in (0, 1) or NULL. Significance level for the
#'   contrast CI (from the bootstrap object). NULL means no CI in label.
#' @return Data.frame with columns `k`, `annotation`, `y_low`, `y_high`,
#'   `y_label`, `label`.
#' @family internal
#' @keywords internal
build_contrast_annotations <- function(risk_method, have_bands, boot_method,
                                       contrast_annotations, eval_times,
                                       alpha_val = NULL) {

  # Mapping from annotation name to (arm_target, arm_reference).
  # For the bridge: two curves of arms[1] and arms[2]. The contrast is
  # arms[1] - arms[2] (i.e. target minus reference).
  pair_map <- list(
    total          = c("arm_11", "arm_00"),
    sep_direct_A   = c("arm_10", "arm_00"),
    sep_indirect_A = c("arm_11", "arm_10"),
    sep_direct_B   = c("arm_11", "arm_01"),
    sep_indirect_B = c("arm_01", "arm_00")
  )

  k_grid     <- sort(unique(risk_method$k))
  avail_arms <- unique(risk_method$arm)

  # Per-arm vectorised look-ups indexed by k_grid order.
  arm_value <- function(arm_name) {
    sub <- risk_method[risk_method$arm == arm_name, ]
    sub$value[match(k_grid, sub$k)]
  }
  arm_lower <- function(arm_name) {
    sub <- risk_method[risk_method$arm == arm_name, ]
    sub$lower[match(k_grid, sub$k)]
  }

  # --- Snap eval_times to nearest cut_time and report if snapped ---
  snap_info <- lapply(eval_times, function(et) {
    idx <- which.min(abs(k_grid - et))
    list(idx = idx, k_at = k_grid[idx], requested = et)
  })
  snapped <- vapply(snap_info, function(s) {
    !isTRUE(all.equal(s$k_at, s$requested))
  }, logical(1))
  if (any(snapped)) {
    msg_pieces <- vapply(snap_info[snapped], function(s) {
      sprintf("%g -> %g", s$requested, s$k_at)
    }, character(1))
    message(
      "eval_times snapped to nearest cut_times: ",
      paste(msg_pieces, collapse = ", ")
    )
  }

  # CI quantile bounds (if bootstrap available)
  have_ci <- !is.null(boot_method) && !is.null(alpha_val)
  if (have_ci) {
    lo <- alpha_val / 2
    hi <- 1 - alpha_val / 2
  }

  rows <- list()
  for (ann in contrast_annotations) {
    arms_pair <- pair_map[[ann]]
    a1 <- arms_pair[1]
    a2 <- arms_pair[2]

    if (!all(arms_pair %in% avail_arms)) next  # missing arm — skip

    a1_vals <- arm_value(a1)
    a2_vals <- arm_value(a2)
    a1_lows <- if (have_bands) arm_lower(a1) else NULL
    a2_lows <- if (have_bands) arm_lower(a2) else NULL

    for (si in seq_along(snap_info)) {
      s <- snap_info[[si]]
      idx <- s$idx
      k_at <- s$k_at

      y1 <- a1_vals[idx]
      y2 <- a2_vals[idx]
      est <- y1 - y2

      # Proper per-contrast bootstrap CI: per-replicate difference, then
      # quantile. NOT y1_lower - y2_upper (too wide, ignores covariance).
      if (have_ci &&
          a1 %in% names(boot_method) &&
          a2 %in% names(boot_method)) {
        reps_c <- boot_method[[a1]][, idx] - boot_method[[a2]][, idx]
        c_lo <- unname(quantile(reps_c, probs = lo, na.rm = TRUE))
        c_hi <- unname(quantile(reps_c, probs = hi, na.rm = TRUE))
        label <- sprintf("RD=%.3f\n[%.3f, %.3f]", est, c_lo, c_hi)
      } else {
        label <- sprintf("RD=%.3f", est)
      }

      # Label y-position: below the lower CI band of the lower arm (if we
      # have bands), else below the lower arm curve itself. Offset downward
      # by a small fraction of the plot range so the label doesn't overlap.
      low_arm <- if (y1 < y2) a1 else a2
      if (have_bands) {
        low_lows <- if (low_arm == a1) a1_lows else a2_lows
        y_label <- low_lows[idx] - 0.03
      } else {
        y_label <- min(y1, y2) - 0.03
      }

      rows[[length(rows) + 1]] <- data.frame(
        k = k_at,
        annotation = ann,
        y_low  = min(y1, y2),
        y_high = max(y1, y2),
        y_label = y_label,
        label  = label,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}


#' Plot Effect-Over-Time Curves  (PLACEHOLDER — to be rewritten)
#' @keywords internal
#' @export
plot.causal_competing_risks_contrast <- function(x, ...) {
  message("plot.causal_competing_risks_contrast: rewrite pending. x$contrasts is in long format.")
  invisible(NULL)
}


#' Plot Weight Diagnostics  (PLACEHOLDER — deferred)
#' @keywords internal
#' @export
plot.causal_competing_risks_diagnostic <- function(x, ...) {
  message("plot.causal_competing_risks_diagnostic: not yet implemented.")
  invisible(NULL)
}
