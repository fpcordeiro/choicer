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

#' Extract lower triangular elements (column-major vech)
#'
#' Column-major lower-triangular vectorization. For a K x K matrix M,
#' returns a length K(K+1)/2 vector ordered by column:
#' c(M_11, M_21, ..., M_K1, M_22, M_32, ..., M_K2, ..., M_KK).
#'
#' Note: this is the conventional `vech` ordering. The choicer C++
#' engine uses the row-major variant; see `vech_row()` below.
#' @param M A square matrix
#' @returns Vector of lower triangular elements including diagonal
#' @noRd
vech_col <- function(M) M[lower.tri(M, diag = TRUE)]

# Row-major vech: lower-triangular vectorization in row-major order.
# For a K x K matrix M, returns a length K(K+1)/2 vector
# c(M_11, M_21, M_22, M_31, M_32, M_33, ...).
# Matches the convention used by build_L_mat() and jacobian_vech_Sigma()
# in src/mxlogit.cpp. Contrast with the column-major [vech_col()].
vech_row <- function(M) {
  K <- nrow(M)
  out <- numeric(K * (K + 1) / 2)
  idx <- 1L
  for (i in seq_len(K)) {
    for (j in seq_len(i)) {
      out[idx] <- M[i, j]
      idx <- idx + 1L
    }
  }
  out
}

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

#' Per-column scale vector for a design matrix
#'
#' Returns the per-column scale (sample SD or a robust SD-equivalent) used to
#' standardize a design matrix before optimization. Column names are preserved.
#' @param M A numeric matrix.
#' @param method One of "sd", "mad", or "iqr".
#' @returns Named numeric vector of column scales.
#' @noRd
.column_scales <- function(M, method) {
  scale_fn <- switch(
    method,
    sd  = stats::sd,
    mad = stats::mad,
    iqr = function(x) stats::IQR(x) / 1.349
  )
  apply(M, 2, scale_fn)
}

#' Validate that column scales are not near-zero
#'
#' Raises an informative error for columns in `idx` whose scale is below `eps`
#' (near-constant columns that cannot be standardized).
#' @param s Named numeric vector of column scales.
#' @param method Scaling method (for the error message).
#' @param label Block label (e.g., "fixed-coefficient").
#' @param idx Integer indices of columns to check (default: all).
#' @param eps Numeric threshold below which a scale is "too small".
#' @returns Invisibly, `s`.
#' @noRd
.assert_scales_ok <- function(s, method, label, idx = seq_along(s), eps = 1e-8) {
  bad <- s[idx] < eps
  if (any(bad)) {
    off <- idx[bad]
    stop("scale_vars='", method, "': ", label,
         " column(s) with scale < ", eps, ": ",
         paste0(names(s)[off], "=", signif(s[off], 3), collapse = ", "))
  }
  invisible(s)
}

#' Back-transform scaled-space estimates to natural units
#'
#' Applies the delta-method back-transform
#' `theta_natural = bt_mult * theta_scaled + bt_shift` and
#' `vcov_natural = (bt_mult bt_mult') o vcov_scaled`, re-deriving SEs while
#' guarding against NA/negative variances, and restores parameter names.
#' @param theta_hat Numeric vector of scaled-space estimates.
#' @param vcov_result List with `vcov` (matrix or NULL) and `se`.
#' @param bt_mult Numeric multiplier vector (length n_params).
#' @param bt_shift Numeric shift vector (length n_params).
#' @param param_names Character vector of parameter names.
#' @returns List with `theta` and `vcov_result`.
#' @noRd
.backtransform_estimates <- function(theta_hat, vcov_result, bt_mult, bt_shift,
                                     param_names) {
  theta_hat <- theta_hat * bt_mult + bt_shift
  names(theta_hat) <- param_names
  if (!is.null(vcov_result$vcov)) {
    vcov_result$vcov <- vcov_result$vcov * tcrossprod(bt_mult)
    rownames(vcov_result$vcov) <- param_names
    colnames(vcov_result$vcov) <- param_names
    diag_v <- diag(vcov_result$vcov)
    se <- rep(NA_real_, length(theta_hat))
    ok <- !is.na(diag_v) & diag_v >= 0
    se[ok] <- sqrt(diag_v[ok])
    names(se) <- param_names
    vcov_result$se <- se
  }
  list(theta = theta_hat, vcov_result = vcov_result)
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
