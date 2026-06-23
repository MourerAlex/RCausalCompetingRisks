#' Null-Coalescing Operator
#'
#' Returns `x` if not NULL, otherwise `y`. Used internally.
#'
#' @param x,y Values to coalesce.
#' @return `x` if not NULL, else `y`.
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#' Per-arm bootstrap matrices for one method
#'
#' Pivots the long `replicates` table (filtered to `method`) into a named
#' list of `[n_boot x n_times]` matrices, one per arm — rows are bootstrap
#' replicates (sorted by `boot_id`), columns are cut times (sorted by `k`).
#' Mirrors the old squeezed array slice `boot_array[, method, arm, ]`. Arms
#' absent from the method get an all-`NA` matrix so downstream cross-arm
#' algebra keeps its shape (and yields `NA` CIs, matching the `NA` point
#' estimates for unsupported decompositions).
#'
#' @keywords internal
boot_arm_matrices <- function(replicates, method, arms = NULL) {
  reps_m   <- replicates[replicates$method == method, , drop = FALSE]
  boot_ids <- sort(unique(reps_m$boot_id))
  ks       <- sort(unique(reps_m$k))
  present  <- unique(reps_m$arm)
  if (is.null(arms)) arms <- present
  na_mat <- matrix(NA_real_, nrow = length(boot_ids), ncol = length(ks))
  mats <- lapply(arms, function(a) {
    if (!a %in% present) return(na_mat)
    t(tapply(reps_m$value[reps_m$arm == a],
             list(reps_m$k[reps_m$arm == a], reps_m$boot_id[reps_m$arm == a]),
             identity))
  })
  names(mats) <- arms
  mats
}
