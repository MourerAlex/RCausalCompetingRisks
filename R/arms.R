#' Specification of the Four Treatment Arms
#'
#' The separable-effects framework decomposes the binary treatment `A` into
#' two binary components `A_Y` (path `A → Y`, the direct channel) and
#' `A_D` (path `A → D → Y`, the indirect channel), so the counterfactual
#' estimand is indexed by the pair `(a_Y, a_D) ∈ {0, 1}^2` — exactly four
#' arms.
#'
#' Each row of [arm_spec()] gives one arm's name, its `(a_y, a_d)` pair
#' (the semantic key), and a default plotting color and label. All
#' modules that iterate the arm space — IPW dispatcher, bootstrap array
#' construction, plotting, long-format `$risk` pivot — read from this
#' helper instead of hardcoding the four names.
#'
#' Default colors follow the Okabe-Ito colorblind-safe palette.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{name}{Character. Arm name (e.g. `"arm_11"`).}
#'     \item{a_y, a_d}{Integer (0 or 1). The two component-treatment
#'       values that define the arm.}
#'     \item{color}{Character. Default hex color for plotting.}
#'     \item{label}{Character. Default human-readable legend label.}
#'   }
#' @family internal
#' @keywords internal
arm_spec <- function() {
  data.frame(
    name  = c("arm_11", "arm_00", "arm_10", "arm_01"),
    a_y   = c(1L, 0L, 1L, 0L),
    a_d   = c(1L, 0L, 0L, 1L),
    color = c("#000000", "#0072B2", "#009E73", "#D55E00"),
    label = c("(1,1) Treated",
              "(0,0) Control",
              "(1,0) Separable direct (A_Y)",
              "(0,1) Separable direct (A_D)"),
    stringsAsFactors = FALSE
  )
}
