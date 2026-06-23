#' Prepare Subject-Level Competing-Events Data for [causal_competing_risks()]
#'
#' Discretizes subject-level competing-events data (one row per subject,
#' a single multi-valued `event` column) into person-time format. This is
#' a thin wrapper over [CausalSurvival::to_person_time()]: the generic
#' single-event machinery (the `(k-1, k]` grid, `T_max` administrative
#' truncation, three-way censoring split, treatment standardization) is
#' inherited, and the competing event D is layered on top.
#'
#' The competing event is passed through the single-event engine as an
#' event (not as censoring) so that a D landing exactly at `T_max` is not
#' lost to admin-truncation, then recovered into its own `d_event` flag.
#'
#' @param data Subject-level data.frame (one row per subject).
#' @param id,time,treatment,covariates See [CausalSurvival::to_person_time()].
#' @param event Character. Name of the multi-valued event-type column.
#' @param event_y,event_d,event_c Values in `event` marking the primary
#'   event (Y), the competing event (D), and censoring, respectively.
#' @param ipcw,T_max,cut_points,time_varying See
#'   [CausalSurvival::to_person_time()]. `cut_points` follows the CS
#'   convention: `NULL` = 12 equi-spaced intervals; a positive integer =
#'   that many intervals; a vector = interior cut points within `(0, T_max)`.
#' @param ... Unused. Reserved to catch the renamed `n_intervals` argument
#'   with a clear error.
#'
#' @return A `c("person_time", "data.frame")` carrying CS's columns plus
#'   `d_event`, a single collapsed `c_event` (CS's `cond_indep_cens` /
#'   `indep_cens` merged), `A_y`, `A_d`, and an `event_labels` attribute
#'   `list(y, d, c)`.
#'   Per row at most one of `y_event, d_event, c_event` is 1.
#'
#' @seealso [CausalSurvival::to_person_time()], [causal_competing_risks()]
#' @export
to_person_time_competing <- function(data,
                           id = "id",
                           time = "event_time",
                           event = "event_type",
                           treatment = "A",
                           covariates = character(),
                           event_y,
                           event_d,
                           event_c,
                           ipcw = TRUE,
                           T_max = NULL,
                           cut_points = NULL,
                           time_varying = NULL,
                           ...) {

  # 0a. Catch the renamed n_intervals argument (and any other stray dots)
  #     with a clear error instead of R's bare "unused argument".
  if (...length() > 0) {
    dots <- names(list(...))
    if ("n_intervals" %in% dots) {
      stop("'n_intervals' was renamed: pass cut_points = <integer> ",
           "(number of intervals).", call. = FALSE)
    }
    stop("unused argument(s): ",
         paste(if (is.null(dots)) "<unnamed>" else dots, collapse = ", "),
         call. = FALSE)
  }

  # 0b. Basic input shape (mirrors CS's checks, which would otherwise only
  #     fire after the event guard below has already read names(data)).
  if (is.null(data)) {
    stop("data is NULL.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("data must be a data.frame. Got class: ",
         paste(class(data), collapse = ", "), call. = FALSE)
  }
  if (nrow(data) == 0) {
    stop("data has 0 rows.", call. = FALSE)
  }

  # 1. Competing-events guard (CCR-specific): event column present, the
  #    three declared labels distinct, and no event value outside them.
  if (!event %in% names(data)) {
    stop("event column '", event, "' not found in data.", call. = FALSE)
  }
  declared <- c(event_y, event_d, event_c)
  if (anyDuplicated(declared) > 0) {
    stop("event_y, event_d, event_c must be distinct.", call. = FALSE)
  }
  extra <- setdiff(unique(data[[event]]), declared)
  if (length(extra) > 0) {
    stop("event column '", event, "' has value(s) not in ",
         "{event_y, event_d, event_c}: ",
         paste(extra, collapse = ", "), ".", call. = FALSE)
  }

  # 2. Way B: one binary status = 1 for EITHER event (Y or D),
  #    0 for censoring. D rides the "had-the-event" path so it cannot be
  #    eaten by CS's admin-truncation rule when it lands at t == T_max.
  data[[".status_ccr"]] <- as.integer(data[[event]] %in% c(event_y, event_d))

  # 3. Delegate discretization to CausalSurvival: inherits the (k-1, k]
  #    grid, no time shift, k = findInterval, T_max admin-truncation,
  #    three-way censoring split, and treatment standardization.
  #    Explicit :: keeps the delegation to the single-event engine
  #    unambiguous regardless of package attach order.
  pt <- CausalSurvival::to_person_time(
    data         = data,
    id           = id,
    time         = time,
    status       = ".status_ccr",
    ipcw         = ipcw,
    T_max        = T_max,
    treatment    = treatment,
    covariates   = covariates,
    cut_points   = cut_points,
    time_varying = time_varying
  )

  # 4. Recover the competing event: CS flagged every event row as
  #    y_event == 1. The exit rows of D-subjects become d_event; their
  #    y_event is cleared. ($<- preserves CS attrs + person_time class.)
  id_d        <- data[[id]][data[[event]] == event_d]
  d_exit_rows <- pt$y_event == 1L & pt[[id]] %in% id_d
  pt$d_event              <- as.integer(d_exit_rows)
  pt$y_event[d_exit_rows] <- 0L

  # 4b. Collapse CS's two-way censoring split (cond_indep_cens / indep_cens)
  #     into CCR's single censoring indicator. CCR does not expose per-subject
  #     independent censoring; under the default scalar ipcw all censoring is
  #     conditionally-independent anyway, so the split carries no extra
  #     information here.
  pt$c_event         <- as.integer(pt$cond_indep_cens | pt$indep_cens)
  pt$cond_indep_cens <- NULL
  pt$indep_cens      <- NULL

  # 5. Cross-arm working copies of the standardized treatment, used by
  #    the g-formula / IPW workers to predict under the (A_y, A_d) arms.
  pt$A_y <- pt[[treatment]]
  pt$A_d <- pt[[treatment]]

  # 6. Stamp the three-slot event-label map (CS only knows Y + censoring).
  attr(pt, "event_labels") <- list(y = event_y, d = event_d, c = event_c)
  pt
}
