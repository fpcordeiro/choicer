# Internal utility functions for choicer package

#' Convert elapsed time to formatted string
#' @param time A proc_time object from system.time()
#' @returns Formatted string like "0h:1m:23s"
#' @noRd
convertTime <- function(time) {
  et <- time["elapsed"]
  if (et < 1) {
    s <- round(et, 2)
  } else {
    s <- round(et, 0)
  }
  h <- s %/% 3600
  s <- s - 3600 * h
  m <- s %/% 60
  s <- s - 60 * m
  return(paste(h, "h:", m, "m:", s, "s", sep = ""))
}

#' Remove columns in the null space of a matrix
#' @param mat A matrix
#' @param tol Tolerance for rank determination
#' @returns Matrix with linearly dependent columns removed
#' @noRd
remove_nullspace_cols <- function(mat, tol = 1e-7) {
  if (is.null(mat)) return(mat)
  if (ncol(mat) == 1) return(mat)
  qrdecomp <- qr(mat, tol = tol)
  rank <- qrdecomp$rank
  if (rank == ncol(mat)) return(mat)
  bad_cols_idx <- qrdecomp$pivot[(rank + 1):ncol(mat)]
  mat <- mat[, setdiff(1:ncol(mat), bad_cols_idx), drop = FALSE]
  return(mat)
}

#' Check for collinearity and remove dependent columns
#' @param X A matrix
#' @returns List with `mat` (cleaned matrix) and `dropped` (names of dropped columns)
#' @noRd
check_collinearity <- function(X) {
  colnames_before <- colnames(X)
  X <- remove_nullspace_cols(X)
  colnames_after <- colnames(X)
  colnames_diff <- setdiff(colnames_before, colnames_after)
  if (length(colnames_diff) > 0) {
    cat("The following variables were dropped due to collinearity:\n")
    cat(colnames_diff, "\n")
  }
  return(list(mat = X, dropped = colnames_diff))
}

#' Extract lower triangular elements (vech)
#' @param M A square matrix
#' @returns Vector of lower triangular elements including diagonal
#' @noRd
vech <- function(M) M[lower.tri(M, diag = TRUE)]
