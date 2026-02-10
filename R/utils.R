# Internal utility functions for choicer package

#' Null-coalescing operator
#'
#' Returns `y` when `x` is NULL, otherwise returns `x`.
#' @param x Value to check.
#' @param y Default value if `x` is NULL.
#' @returns `x` if not NULL, otherwise `y`.
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

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
    message("The following variables were dropped due to collinearity: ",
            paste(colnames_diff, collapse = ", "))
  }
  return(list(mat = X, dropped = colnames_diff))
}

#' Extract lower triangular elements (vech)
#' @param M A square matrix
#' @returns Vector of lower triangular elements including diagonal
#' @noRd
vech <- function(M) M[lower.tri(M, diag = TRUE)]

#' Resolve a variable name or index to a 1-based integer index
#'
#' Accepts a character variable name or a 1-based integer index. Validates
#' against the column names of the relevant matrix.
#'
#' @param var Character name or integer index.
#' @param col_names Character vector of valid column names.
#' @returns Integer index (1-based).
#' @noRd
resolve_var_index <- function(var, col_names) {
  if (is.character(var)) {
    idx <- match(var, col_names)
    if (is.na(idx)) {
      stop("Variable '", var, "' not found. Available: ",
           paste(col_names, collapse = ", "))
    }
    return(idx)
  }
  if (is.numeric(var)) {
    var <- as.integer(var)
    if (var < 1L || var > length(col_names)) {
      stop("Variable index ", var, " out of range [1, ", length(col_names), "].")
    }
    return(var)
  }
  stop("'elast_var' must be a character variable name or integer index.")
}

#' Label a J x J matrix with alternative names
#'
#' Adds row and column names from \code{alt_mapping} to a square matrix.
#'
#' @param mat A J x J matrix.
#' @param alt_mapping data.table with alternative labels in column 2.
#' @returns The matrix with row/column names set.
#' @noRd
label_matrix <- function(mat, alt_mapping) {
  alt_labels <- alt_mapping[[2]]
  if (length(alt_labels) == nrow(mat)) {
    rownames(mat) <- alt_labels
    colnames(mat) <- alt_labels
  }
  mat
}
