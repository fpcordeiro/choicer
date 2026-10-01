# Internal utility functions for choicer package

#' Null-coalescing operator
#'
#' Returns `y` when `x` is NULL, otherwise returns `x`.
#' @param x Value to check.
#' @param y Default value if `x` is NULL.
#' @returns `x` if not NULL, otherwise `y`.
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

# Measured Rcpp::wrap()+R-list-overhead multiplier for the beta_i draw cube
# (profiling report `_benchmarks/hb_profiling_report.md` §4: true RSS is
# 1.7-1.9x the naive 8*K*N*R_keep byte estimate). Using the conservative
# (higher) end for a fail-fast guard.
.HB_BETA_I_WRAP_FACTOR <- 1.9

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
#'
#' Keeps the columns that `qr(mat, tol)` retains. Designs with a million or
#' more elements are factored by row chunks (`.qr_rank_pivot()`), without a
#' design-sized copy.
#' @param mat A matrix
#' @param tol Tolerance for rank determination
#' @returns Matrix with linearly dependent columns removed
#' @noRd
remove_nullspace_cols <- function(mat, tol = 1e-7) {
  if (is.null(mat)) return(mat)
  if (ncol(mat) == 1) return(mat)
  qrdecomp <- if (as.numeric(nrow(mat)) * ncol(mat) < 1e6) {
    qr(mat, tol = tol)
  } else {
    .qr_rank_pivot(mat, tol)
  }
  rank <- qrdecomp$rank
  if (rank == ncol(mat)) return(mat)
  bad_cols_idx <- qrdecomp$pivot[(rank + 1):ncol(mat)]
  mat <- mat[, setdiff(1:ncol(mat), bad_cols_idx), drop = FALSE]
  return(mat)
}

#' Rank and column pivot of qr(mat, tol) for a tall matrix, by row chunks
#'
#' Tall-skinny QR: the R factor is updated chunk by chunk,
#' R = qr.R(qr(rbind(R, chunk))), so it keeps the Gram matrix of the rows seen
#' so far and never grows past ncol(mat) rows. qr()'s LINPACK rule (limited
#' column pivoting against `tol` times each column's original norm) depends
#' on `mat` only through that Gram matrix, so applied to the final R it makes
#' the same decisions in exact arithmetic; they can differ only where
#' rounding decides, and there qr() itself changes with row order or BLAS.
#' The chunks go through LAPACK (dgeqp3: about 2p matrix-vector BLAS calls
#' per chunk, against LINPACK's p^2 vector calls, whose per-call overhead
#' can dominate at this size); its column pivoting is undone, as only the
#' Gram matrix matters. Beyond `mat`, memory holds a few chunk-sized buffers
#' and a p x p factor, and there is no limit of 2^31 - 1 elements.
#'
#' @param mat A numeric matrix.
#' @param tol Tolerance for rank determination, as in qr().
#' @param rows Rows per chunk (about 2^20 elements by default).
#' @returns The qr() of the final R factor (its rank and pivot are used).
#' @noRd
.qr_rank_pivot <- function(mat, tol = 1e-7,
                           rows = max(ncol(mat), floor(2^20 / ncol(mat)))) {
  n <- nrow(mat)
  R <- NULL
  for (s in seq(1, n, by = rows)) {
    chunk <- mat[s:min(n, s + rows - 1), , drop = FALSE]
    q <- qr(rbind(R, chunk), LAPACK = TRUE)
    R <- qr.R(q)[, order(q$pivot), drop = FALSE]
  }
  qr(R, tol = tol)
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

#' Stop when a fit's parameter names repeat
#'
#' Coefficients, `vcov()`, `summary()`, `wtp()` and named bounds look
#' parameters up by name, so a covariate named like a parameter the model
#' generates would leave two parameters under one name. Each element of
#' `labels` names the whole parameter vector (the stored names, or the
#' labels `summary()` prints) and is checked on its own.
#'
#' @param labels List of character vectors of parameter names, the stored
#'   names first.
#' @param generated The model's generated names, for the message.
#' @returns `NULL`, invisibly; stops listing the repeated names in parameter
#'   order.
#' @noRd
.check_param_names <- function(labels, generated) {
  dups <- unique(unlist(lapply(labels, function(nm) nm[duplicated(nm)])))
  if (length(dups) > 0L) {
    dups <- dups[order(match(dups, labels[[1L]]))]
    stop("Parameter names must be unique; repeated: ",
         paste0("'", dups, "'", collapse = ", "),
         ". Covariates may not take the names the model gives its own ",
         "parameters (", paste(generated, collapse = ", "), ").",
         call. = FALSE)
  }
  invisible(NULL)
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
  bad <- !is.finite(s[idx]) | s[idx] < eps
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

# Columns of the `alt_mapping` summary that the prepare_*_data() functions
# return beside the user's alternative column.
.ALT_MAPPING_COLS <- c("alt_int", "N_OBS", "N_CHOICES", "TAKE_RATE", "MKT_SHARE")

#' Reject user column names that collide with choicer's own
#'
#' The prepare_*_data() functions copy the user's index columns (ids,
#' alternatives, choices, weights, clusters, decision makers) into a private
#' working table. Every working column they add to it, and every result they
#' read back from a table grouped by a user column, is named with the prefix
#' ".choicer_" or read by position, so user columns may not carry that
#' prefix. The alternative column reappears in the returned `alt_mapping`
#' beside fixed columns, so it may not take one of their names.
#'
#' @param needed Names of the user columns the preparation uses.
#' @param alt_col Name of the alternative column.
#' @returns `NULL`, invisibly; errors on a collision.
#' @noRd
.check_col_names <- function(needed, alt_col) {
  internal <- unique(needed[startsWith(needed, ".choicer_")])
  if (length(internal) > 0) {
    stop("Column names beginning with '.choicer_' are reserved for ",
         "choicer's working columns; rename ",
         paste0("'", internal, "'", collapse = ", "), ".", call. = FALSE)
  }
  if (alt_col %in% .ALT_MAPPING_COLS) {
    stop("The alternative column cannot be named '", alt_col, "': the ",
         "returned `alt_mapping` reserves ",
         paste0("'", .ALT_MAPPING_COLS, "'", collapse = ", "),
         " for its own columns. Rename the alternative column.",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Stop when a covariate is also an index column
#'
#' The MNL, MXL and NL preparations dropped a covariate that is also their
#' id, alternative, choice, weight, cluster or decision-maker column from
#' their working table after gathering the design. When they next read the
#' column they stopped with an unrelated error, or, for the choice column,
#' counted any object of its name visible to the call instead (a vector in
#' the workspace, say). Such a covariate now stops them, at the same point,
#' with an error that says why.
#'
#' @param covariates Names of the covariate columns.
#' @param index_cols Names of the preparation's index columns.
#' @returns `NULL`, invisibly; errors on an overlap.
#' @noRd
.check_index_covariates <- function(covariates, index_cols) {
  both <- unique(covariates[covariates %in% index_cols])
  if (length(both) > 0) {
    stop("A covariate cannot also be the id, alternative, choice, weight, ",
         "cluster or decision-maker column: ",
         paste0("'", both, "'", collapse = ", "), ". To use one as a ",
         "covariate, copy it under another name.", call. = FALSE)
  }
  invisible(NULL)
}

#' Code the alternatives 1..J in a preparation's working table
#'
#' Adds `.choicer_alt_int`, the position of each row's alternative in
#' `levels`. The codes are computed outside `dt[...]`, where a user column
#' named `levels` would mask the argument, and enter it under a reserved
#' name, which no user column can take. `:=` grows a table that has no spare
#' column slot (`options(datatable.alloccol = 0)`) by reallocating it, so the
#' caller must keep the returned table.
#'
#' @param dt The preparation's working data.table.
#' @param alt_col Name of the alternative column.
#' @param levels Alternative labels in code order.
#' @returns `dt` with the column added, invisibly.
#' @noRd
.code_alternatives <- function(dt, alt_col, levels) {
  .choicer_codes <- as.integer(factor(dt[[alt_col]], levels = levels))
  dt[, .choicer_alt_int := .choicer_codes]
}

#' Collapse a row-level column to one value per choice situation
#'
#' Used by the prepare_*_data() functions for columns that must be constant
#' within a choice situation (e.g. cluster labels). Called AFTER filtering
#' and ordering, so alignment is by id, never by position.
#'
#' @param dt data.table after filtering/ordering (long format).
#' @param col Name of the column to collapse.
#' @param id_col Name of the choice-situation id column.
#' @param ids Vector of situation ids in prepared order.
#' @returns Vector of length \code{length(ids)}, one value per situation.
#' @noRd
.collapse_situation_col <- function(dt, col, id_col, ids) {
  nuniq <- .n_distinct_by(dt, col, id_col)
  if (any(nuniq != 1L)) {
    stop("`", col, "` must be constant within each '", id_col,
         "' (one value per choice situation).", call. = FALSE)
  }
  .first_by(dt, col, id_col, ids)
}

#' Number of distinct values of a column within each choice situation
#'
#' `dt[, uniqueN(get(col)), by = id_col]$V1` without sorting every situation
#' separately (a per-group uniqueN() allocates sort buffers for each of
#' millions of groups): counts the distinct (id, value) pairs of each id.
#' The pair table has fixed column names, so any user column name works
#' (including V1, N, or names with commas).
#'
#' @param dt A data.table.
#' @param col Name of the column whose values are counted.
#' @param id_col Name of the choice-situation id column.
#' @returns Integer vector with one count per choice situation, in order of
#'   first appearance of each id. NA is a value, as in uniqueN(); so is NaN,
#'   apart from NA; -0 equals 0.
#' @noRd
.n_distinct_by <- function(dt, col, id_col) {
  # Read-only: pairs shares dt's column vectors (setDT() copies nothing).
  pairs <- data.table::setDT(list(id = dt[[id_col]], value = dt[[col]]))
  pairs[, .N, by = c("id", "value")][, .N, by = "id"][["N"]]
}

#' Value of a column at the first row of each choice situation
#'
#' `dt[, get(col)[1L], by = id_col]` aligned to `ids`, without evaluating R
#' code for every situation and without reading the result by the name "V1",
#' which is also the id column's name when that column is called V1.
#'
#' @param dt A data.table.
#' @param col Name of the column to read.
#' @param id_col Name of the choice-situation id column.
#' @param ids Situation ids in prepared order.
#' @returns Vector of length `length(ids)`.
#' @noRd
.first_by <- function(dt, col, id_col, ids) {
  first <- !duplicated(dt[[id_col]])
  dt[[col]][first][match(ids, dt[[id_col]][first])]
}

#' Whether a choice column holds only 0 and 1
#'
#' As numbers or logicals: a character or factor column of "0" and "1"
#' passes `%in% c(0, 1)` by coercion, but its choices cannot be summed.
#'
#' @param x The choice column.
#' @returns `TRUE` or `FALSE`.
#' @noRd
.is_zero_one <- function(x) {
  !is.character(x) && !is.factor(x) && all(x %in% c(0, 1))
}

#' A working table's choice-situation keys under fixed names
#'
#' A private table holding the columns `by` of `dt` as `k1`, `k2`, ..., and
#' the further columns given in `...`. setDT() copies nothing: the table
#' shares `dt`'s vectors. Grouped by those names, it meets no user column
#' name, so a column named like a local (`id_col`, `choice_col`) cannot take
#' the local's place, as it would inside `dt[...]`.
#'
#' @param dt A data.table.
#' @param by Names of the columns identifying a choice situation.
#' @param ... Further named columns, as vectors.
#' @returns A data.table.
#' @noRd
.key_table <- function(dt, by, ...) {
  keys <- lapply(by, function(cc) dt[[cc]])
  names(keys) <- paste0("k", seq_along(by))
  data.table::setDT(c(keys, list(...)))
}

#' Number of rows of each choice situation
#'
#' `dt[, .N, by = by]` counted on a private table (`.key_table()`).
#'
#' @param dt A data.table.
#' @param by Names of the columns identifying a choice situation.
#' @returns Integer vector, one count per situation in order of first
#'   appearance.
#' @noRd
.n_rows_by <- function(dt, by) {
  keys <- paste0("k", seq_along(by))
  .key_table(dt, by)[, .N, by = keys][["N"]]
}

#' Number of chosen rows of each choice situation
#'
#' `dt[, sum(get(choice_col)), by = by]` summed on a private table
#' (`.key_table()`), by data.table in one pass in C (GForce) rather than by
#' an R call per situation.
#'
#' @param dt A data.table.
#' @param choice_col Name of the 0/1 choice column (numeric or logical).
#' @param by Names of the columns identifying a choice situation.
#' @returns One count per situation in order of first appearance, for
#'   comparison with 1: unlike `sum()`, GForce keeps the choice column's
#'   attributes.
#' @noRd
.n_chosen_by <- function(dt, choice_col, by) {
  keys <- paste0("k", seq_along(by))
  .key_table(dt, by, chosen = dt[[choice_col]])[
    , list(n = sum(chosen)), by = keys][["n"]]
}

#' Rows and choices of each alternative
#'
#' The first columns of a preparation's `alt_mapping`: `alt_int`, the
#' alternative's code, the alternative column, and the alternative's numbers
#' of rows (`N_OBS`) and of choices (`N_CHOICES`), keyed by the first two.
#' Counted on a private table with fixed names, as `.n_chosen_by()` counts
#' choices, but summed by `base::sum()` for each alternative: GForce's sum
#' would keep the choice column's attributes (those of a pdata.frame's
#' columns, say), which `sum()` drops, and there are only J groups.
#'
#' @param dt A preparation's working table, holding `.choicer_alt_int`.
#' @param alt_col Name of the alternative column.
#' @param choice_col Name of the 0/1 choice column.
#' @returns A data.table keyed by `alt_int` and `alt_col`.
#' @noRd
.alt_counts <- function(dt, alt_col, choice_col) {
  tab <- data.table::setDT(list(alt_int = dt$.choicer_alt_int,
                                alt = dt[[alt_col]],
                                chosen = dt[[choice_col]]))
  counts <- tab[, list(N_OBS = .N, N_CHOICES = base::sum(chosen)),
                keyby = c("alt_int", "alt")]
  data.table::setnames(counts, "alt", alt_col)
  counts
}

#' Copy the named columns of the user's data into a new data.table
#'
#' Deep-copies only the columns of `data` named in `cols`, in input order, so
#' the preparation that follows (`:=`, `setorderv()` by reference) never
#' touches the caller's data and never copies columns it does not use.
#' `data.table::as.data.table()` of a list returns a deep copy (documented in
#' ?as.data.table; checked in the 1.10.4-3 and 1.18.6.1 sources).
#'
#' Other inputs keep the preps' previous route, a full `as.data.table()` copy
#' with the unused columns dropped, because `as.data.table()` itself changes
#' them: it splits matrix and data-frame columns into renamed columns that
#' `cols` may refer to, converts POSIXlt columns, and (data.table >= 1.17.0)
#' passes data frame subclasses such as plm's pdata.frame through their own
#' `as.data.frame()` method; with duplicated names its `:=` removes only the
#' first match. The direct route is therefore limited to data.frame,
#' data.table and tibble inputs with unique names, no matrix or data-frame
#' columns and atomic `cols`. Missing columns are left for the caller to
#' report.
#'
#' @param data User data.
#' @param cols Names of the columns to copy; duplicates allowed.
#' @returns A data.table.
#' @noRd
.copy_cols <- function(data, cols) {
  if (.direct_route(data, cols)) {
    keep <- which(names(data) %in% cols)
    return(data.table::as.data.table(.subset(data, keep)))
  }
  dt <- data.table::as.data.table(data)[]
  vars_to_drop <- setdiff(names(dt), cols)
  if (length(vars_to_drop) > 0) {
    dt[, (vars_to_drop) := NULL]
  }
  dt
}

#' Whether the preps may read `cols` of the user's data directly
#'
#' See `.copy_cols()` for why other inputs take the full-copy route.
#'
#' @param data User data.
#' @param cols Names of the columns to be read.
#' @returns `TRUE` or `FALSE`.
#' @noRd
.direct_route <- function(data, cols) {
  nm <- names(data)
  (identical(class(data), "data.frame") || data.table::is.data.table(data) ||
     inherits(data, "tbl_df")) &&
    !anyDuplicated(nm) &&
    all(vapply(unclass(data), function(v) is.null(dim(v)), logical(1L))) &&
    all(vapply(which(nm %in% cols), function(j) is.atomic(.subset2(data, j)),
               logical(1L)))
}

#' The data frame the preps read covariates from, and the columns to scan
#'
#' The user's data itself, read in place and never modified, when
#' `.direct_route()` allows; then only the needed columns are scanned for
#' missing values. Otherwise the table the preps have always built: a full
#' `as.data.table()` copy with the unused columns dropped by `:=`, all of
#' whose columns are scanned, as `.SD` was (with duplicated names, `:=` drops
#' only the first match, and the survivor's missing values still count).
#' When a scanned column is integer64, bit64 is loaded first (see
#' `.load_bit64_for()`), so that the checks reading `src` in place see its
#' values.
#'
#' @param data User data.
#' @param needed Names of the columns the model uses.
#' @returns A list: `src`, `data` or a data.table copy of it, and `scan`,
#'   positions of the columns of `src` to scan for missing values.
#' @noRd
.prep_source <- function(data, needed) {
  if (.direct_route(data, needed)) {
    src <- data
    scan <- which(names(data) %in% needed)
  } else {
    src <- .copy_cols(data, needed)
    scan <- seq_along(src)
  }
  .load_bit64_for(src, scan)
  list(src = src, scan = scan)
}

#' Rows with a missing value in any of the given columns
#'
#' Column-by-column `rowSums(is.na(x[cols])) > 0`, without its rows-by-columns
#' logical matrix; columns without missing values cost one pass of anyNA().
#' integer64 columns need bit64's methods loaded (`.load_bit64_for()`).
#'
#' @param x A data frame or data.table.
#' @param cols Column positions to scan (default: all).
#' @returns Logical vector with one element per row.
#' @noRd
.rows_with_na <- function(x, cols = seq_along(x)) {
  has_na <- logical(nrow(x))
  for (j in cols) {
    col <- .subset2(x, j)
    if (anyNA(col)) has_na <- has_na | is.na(col)
  }
  has_na
}

#' Stop before building a design matrix the kernels cannot address
#'
#' choicer builds with RcppArmadillo's default 32-bit word (ARMA_32BIT_WORD),
#' and the kernels view X and W without copying them: past 2^32 - 1 elements
#' the element count wraps and rows are misread without an error.
#'
#' @param n Number of rows.
#' @param cols Column names.
#' @param what Name of the matrix, for the message.
#' @returns Invisibly, `NULL`.
#' @noRd
.check_design_size <- function(n, cols, what) {
  size <- as.numeric(n) * length(cols)
  if (size > 2^32 - 1) {
    stop(what, " would have ", format(n, big.mark = ",", scientific = FALSE),
         " rows and ", length(cols), " columns, ",
         format(size, big.mark = ",", scientific = FALSE), " values, more ",
         "than 2^32 - 1, the most the estimation kernels can address.",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Design matrix gathered from the source columns
#'
#' `X[i, ] = src[rows[i], cols]` in double storage: the values and layout of
#' `as.matrix()` on the sorted rows, filled in C++ straight from `src`
#' without copying the covariates into the prep's working table, reading
#' integer columns' raw storage as `as.matrix()` does. integer64 columns are
#' read as their values instead, as `.int64_to_double()` converts them and
#' with its warning for values of magnitude 2^53 or more among the rows read.
#' The estimation kernels take double matrices and would otherwise convert
#' an all-integer design on every call. Columns are looked up by position,
#' first match, as `dt[, ..cols]` did. With `base`, the rows are differenced
#' in the same pass, `X[i, ] = src[rows[i], cols] - src[base[i], cols]`, as
#' the multinomial probit's design is.
#'
#' @param src Data frame holding the columns (see `.prep_source()`).
#' @param cols Names of numeric columns of `src`.
#' @param rows Row indices into `src`, in prepared order.
#' @param what Name of the matrix, for the size check's message.
#' @param base `NULL`, or row indices into `src`, one per element of `rows`,
#'   of the values to subtract.
#' @returns A `length(rows)` x `length(cols)` double matrix with column names
#'   `as.character(cols)` (a named `cols` leaves no names behind, as with
#'   `as.matrix()`); for zero columns, `as.matrix()`'s 0 x 0 logical matrix,
#'   which the callers' final checks reject.
#' @noRd
.gather_matrix <- function(src, cols, rows, what, base = NULL) {
  if (!length(cols)) return(as.matrix(data.table::data.table()))
  .check_design_size(length(rows), cols, what)
  X <- prep_gather_design(lapply(match(cols, names(src)),
                                 function(j) .subset2(src, j)), rows, base)
  big <- attr(X, "int64_big")
  if (!is.null(big)) {
    attr(X, "int64_big") <- NULL
    for (cc in unique(as.character(cols)[big])) {
      .warn_int64_big(paste0("Column '", cc, "'"))
    }
  }
  dimnames(X) <- list(NULL, as.character(cols))
  X
}

#' Whether prepared choice situations are in ascending-id order
#'
#' The cross-sectional prepared order is ascending id; a panel mixed logit
#' orders situations by decision maker first. Positional (unnamed)
#' per-situation inputs are unambiguous only when the two coincide. Radix
#' ordering matches data.table's C-locale sort.
#'
#' @param ids Choice-situation ids in prepared order.
#' @returns `TRUE` if `ids` is in ascending order.
#' @noRd
.ids_sorted <- function(ids) {
  identical(order(ids, method = "radix"), seq_along(ids))
}

#' First choice situation of each likelihood unit
#'
#' Likelihood units are decision makers in a panel mixed logit (`d$Ti`
#' non-NULL, situations sorted by decision maker) and choice situations
#' otherwise, so MNL, NL and cross-sectional MXL data get `seq_along(d$M)`.
#'
#' @param d Prepared or stored data: a list with `M` and, for a panel fit, `Ti`.
#' @returns Integer index of each unit's first choice situation.
#' @noRd
.unit_first <- function(d) {
  Ti <- d[["Ti"]]
  if (is.null(Ti)) return(seq_along(d[["M"]]))
  cumsum(c(1L, Ti[-length(Ti)]))
}

#' Collapse a per-situation vector to one value per likelihood unit
#'
#' The identity in the cross-section (`d$Ti` NULL). In a panel mixed logit,
#' errors unless `x` is constant within each decision maker (weights and
#' cluster labels are decision-maker attributes there) and returns each
#' unit's value. NA-safe; works for numeric, character and factor vectors.
#'
#' @param x Vector with one entry per choice situation (prepared order), or
#'   NULL.
#' @param d Prepared or stored data: a list with `M` and, for a panel fit, `Ti`.
#' @param what Label of `x` for the error messages.
#' @returns `x` (cross-section or NULL `x`), else `x[.unit_first(d)]`.
#' @noRd
.to_units <- function(x, d, what) {
  Ti <- d[["Ti"]]
  if (is.null(x) || is.null(Ti)) return(x)
  if (length(x) != sum(Ti)) {
    stop(what, ": got ", length(x), " values for ", sum(Ti),
         " choice situations.", call. = FALSE)
  }
  first <- .unit_first(d)
  code <- match(x, unique(x))  # integer labels; match() pairs NA with NA
  if (any(code != rep(code[first], Ti))) {
    stop(what, " must be constant within each decision maker: with ",
         "`person_col` the likelihood has one term per decision maker, so ",
         "weights are decision-maker weights and clusters must nest ",
         "decision makers.", call. = FALSE)
  }
  x[first]
}

#' Whether bit64 can be loaded
#'
#' Its own function so that the tests can mock bit64's absence.
#'
#' @returns `TRUE` or `FALSE`.
#' @noRd
.bit64_available <- function() requireNamespace("bit64", quietly = TRUE)

#' Stop unless bit64 can be loaded to read an integer64 object
#'
#' @param what Name of the object for the message, e.g. "Column 'x1'".
#' @returns Invisibly, `NULL`.
#' @noRd
.need_bit64 <- function(what) {
  if (!.bit64_available()) {
    stop(what, " is of class integer64; install the bit64 package to read ",
         "it.", call. = FALSE)
  }
  invisible(NULL)
}

#' Warn that integer64 values were rounded to doubles
#'
#' The warning has class `choicer_int64_big`, so that `.int64_warn_once()`
#' can drop repeats.
#'
#' @param what Name of the converted object, e.g. "Column 'x1'".
#' @returns Invisibly, the warning message.
#' @noRd
.warn_int64_big <- function(what) {
  warning(warningCondition(
    paste0(what, " has integer64 values of magnitude 2^53 or more; they ",
           "were rounded to the nearest double."),
    class = "choicer_int64_big"))
}

#' One rounding warning per column across a preparation's conversions
#'
#' A preparation can convert a column more than once: a fixed covariate that
#' is also a random coefficient, a structural covariate that is also an
#' alternative-level one. The returned function evaluates its argument and
#' lets through only the first `.warn_int64_big()` warning about each
#' column, a record its calls share.
#'
#' @returns A function of one argument.
#' @noRd
.int64_warn_once <- function() {
  seen <- character(0)
  function(expr) {
    withCallingHandlers(expr, choicer_int64_big = function(w) {
      if (conditionMessage(w) %in% seen) invokeRestart("muffleWarning")
      seen <<- c(seen, conditionMessage(w))
    })
  }
}

#' Load bit64 when a column read in place is integer64
#'
#' bit64's integer64 keeps 64-bit integers in the storage of a double vector,
#' which base R reads as raw bits: `anyNA()` and `is.na()` miss
#' `NA_integer64_` (its bits are those of -0) and flag small negative values
#' (their bits are NaNs), and `is.finite()`, `[` and `==` misread them alike.
#' bit64 registers methods that read the values when its namespace loads,
#' which holding an integer64 column does not guarantee (one read from an
#' `.rds` file arrives without it).
#'
#' @param x A data frame.
#' @param cols Positions of the columns that will be read.
#' @returns Invisibly, `NULL`; stops when a column is integer64 and bit64
#'   cannot be loaded.
#' @noRd
.load_bit64_for <- function(x, cols = seq_along(x)) {
  for (j in cols) {
    if (inherits(.subset2(x, j), "integer64")) {
      return(.need_bit64(paste0("Column '", names(x)[j], "'")))
    }
  }
  invisible(NULL)
}

#' The values of an integer64 vector as doubles
#'
#' `is.numeric()` accepts integer64 (bit64 defines no method), but code that
#' reads its storage directly, `as.matrix()` and the C++ kernels among it,
#' reads the raw bit patterns as doubles: 1 as 4.9e-324, 2^40 as 5.4e-312,
#' small negative values as NaN. `as.double()` dispatches to bit64's method,
#' which is exact below 2^53 in magnitude and maps `NA_integer64_` to `NA`;
#' larger values round to the nearest double, and the warning bit64 gives for
#' them, which does not say what was converted, is replaced by one that does.
#' The preps' gathers read integer64 columns in place with the same values
#' and warning (`prep_gather_design()`).
#'
#' @param x A vector or matrix, returned unchanged unless it is integer64.
#' @param what Name of `x` for the messages, e.g. "Column 'x1'".
#' @returns `x` as a double vector, or a double matrix with `x`'s `dim` and
#'   `dimnames`.
#' @noRd
.int64_to_double <- function(x, what) {
  if (!inherits(x, "integer64")) return(x)
  .need_bit64(what)
  out <- suppressWarnings(as.double(x))
  if (any(abs(out) >= 2^53, na.rm = TRUE)) .warn_int64_big(what)
  dim(out) <- dim(x)             # as.double() drops them
  dimnames(out) <- dimnames(x)
  out
}

#' Convert the integer64 columns of a private table to double, in place
#'
#' For tables the caller owns and reads with R code: a prep's working table
#' of index columns (for its weight column), and the copies the prediction
#' helpers make of `newdata`. Called before their missing-value checks, so
#' every later step sees doubles. See `.int64_to_double()`.
#'
#' @param dt A data.table owned by the caller; modified in place.
#' @param cols Names of the columns read as numbers (covariates, weights).
#' @returns `dt`, invisibly.
#' @noRd
.int64_cols_to_double <- function(dt, cols) {
  for (j in which(names(dt) %in% cols)) {
    v <- .subset2(dt, j)
    if (inherits(v, "integer64")) {
      data.table::set(dt, j = j, value = .int64_to_double(
        v, paste0("Column '", names(dt)[j], "'")))
    }
  }
  invisible(dt)
}
