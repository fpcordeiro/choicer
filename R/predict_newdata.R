# Counterfactual prediction inputs
#
# Helpers that turn user-supplied `newdata` (a long-format data.frame or a
# pre-built "modified design matrix" list) into the kernel inputs used by the
# predict() methods. Stored coefficients and design matrices are natural-scale,
# so no scaling logic is needed here.

#' Validate prediction weights
#'
#' Returns a length-`N` numeric weight vector, defaulting to equal weights.
#'
#' @param weights Optional numeric vector (or NULL).
#' @param N Number of choice situations.
#' @returns Numeric vector of length `N`.
#' @noRd
.validate_pred_weights <- function(weights, N) {
  if (is.null(weights)) return(rep(1, N))
  if (!is.numeric(weights) || length(weights) != N) {
    stop("'weights' must be a numeric vector with one entry per choice ",
         "situation (expected length ", N, ", got ", length(weights), ").")
  }
  if (any(!is.finite(weights))) {
    stop("'weights' must contain only finite values.")
  }
  as.numeric(weights)
}

#' Extract a numeric design matrix from a data.table
#'
#' Subsets `dt` to `cols` (in that order) and returns a double matrix, so the
#' column order always matches the fitted coefficient order.
#'
#' @param dt A data.table.
#' @param cols Character vector of column names.
#' @returns A numeric matrix with `length(cols)` columns.
#' @noRd
.newdata_matrix <- function(dt, cols) {
  if (length(cols) == 0) {
    return(matrix(numeric(0), nrow = nrow(dt), ncol = 0))
  }
  mat <- as.matrix(dt[, cols, with = FALSE])
  storage.mode(mat) <- "double"
  mat
}

#' Build prediction inputs from a newdata data.frame
#'
#' Converts long-format out-of-sample data into the kernel inputs
#' (`X`, `W`, `alt_idx`, `M`, `N`, `weights`) for counterfactual prediction.
#' Unlike `prepare_mnl_data()` / `prepare_nl_data()`, no choice column is
#' required, alternative codes come from the *fit-time* `alt_mapping` (so ASCs
#' stay aligned even when some alternatives are absent from `newdata`), and no
#' columns are dropped for collinearity. Per-id subsets of alternatives
#' (varying choice sets) are allowed.
#'
#' @param object A fitted `choicer_fit` object.
#' @param newdata A data.frame in the same long format used at fit time:
#'   one row per (id, alternative) pair with the fit-time id, alternative,
#'   and covariate columns. A choice column is not required.
#' @param weights Optional numeric vector with one weight per choice situation
#'   in `newdata` (after dropping outside-option rows). Defaults to 1.
#' @returns List with `X`, `W` (NULL unless MXL), `alt_idx` (1-based integer),
#'   `M` (alternatives per id), `N`, and `weights` (length `N`). Rows are
#'   ordered by id, then by fit-time alternative code.
#' @noRd
prepare_newdata <- function(object, newdata, weights = NULL) {
  spec <- object$data_spec
  if (is.null(spec) || is.null(spec$id_col) || is.null(spec$alt_col)) {
    stop("Cannot resolve newdata columns: the fitted object carries no ",
         "'data_spec' with id/alt column names. Use the list form of ",
         "'newdata' (X, alt_idx, M) instead.")
  }
  id_col <- spec$id_col
  alt_col <- spec$alt_col

  # Required covariates: X columns from the beta block of the coefficient
  # vector; W columns (MXL only) from the stored column-scale names.
  x_cols <- names(object$coefficients)[object$param_map$beta]
  w_cols <- if (identical(object$model, "mxl")) names(object$sW)

  dt <- data.table::as.data.table(newdata)

  needed <- unique(c(id_col, alt_col, x_cols, w_cols))
  missing_cols <- setdiff(needed, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing columns in newdata: ",
         paste(missing_cols, collapse = ", "))
  }

  # Subset to used columns; this always copies, so the user's object is
  # never mutated by the reordering below.
  dt <- dt[, needed, with = FALSE]

  # Remove outside-option rows when modelling the outside option implicitly
  # (mirrors prepare_mnl_data(); the outside option is handled in C++).
  if (isTRUE(object$include_outside_option) &&
      !is.null(spec$outside_opt_label)) {
    dt <- dt[get(alt_col) != spec$outside_opt_label]
    if (nrow(dt) == 0) {
      stop("No inside alternatives remain after removing outside option rows.")
    }
  }

  # No missing values allowed in any used column
  na_cols <- names(dt)[vapply(dt, anyNA, logical(1L))]
  if (length(na_cols) > 0) {
    stop("newdata contains missing values in column(s): ",
         paste(na_cols, collapse = ", "))
  }

  # Covariates must be numeric
  cov_cols <- unique(c(x_cols, w_cols))
  is_num <- vapply(dt[, cov_cols, with = FALSE], is.numeric, logical(1L))
  if (!all(is_num)) {
    stop("All covariates must be numeric. Non-numeric column(s): ",
         paste(cov_cols[!is_num], collapse = ", "))
  }

  # Map alternative labels through the fit-time mapping (authoritative for
  # ASC alignment); labels unseen at fit time are an error.
  am <- object$alt_mapping
  pos <- match(dt[[alt_col]], am[[alt_col]])
  if (anyNA(pos)) {
    unseen <- unique(dt[[alt_col]][is.na(pos)])
    stop("newdata contains alternatives not seen at fit time: ",
         paste(unseen, collapse = ", "))
  }
  dt[, alt_int := am$alt_int[pos]]

  # No duplicate (id, alternative) pairs
  if (anyDuplicated(dt, by = c(id_col, "alt_int")) > 0) {
    stop("newdata contains duplicated (", id_col, ", ", alt_col, ") pairs.")
  }

  # Order rows: ascending id, ascending alternative code within id
  # (same convention as the prepare_*_data() functions).
  data.table::setorderv(dt, c(id_col, "alt_int"))

  X <- .newdata_matrix(dt, x_cols)
  W <- if (!is.null(w_cols)) .newdata_matrix(dt, w_cols)
  alt_idx <- as.integer(dt$alt_int)
  M <- dt[, .N, by = id_col][["N"]]
  N <- length(M)

  list(
    X = X,
    W = W,
    alt_idx = alt_idx,
    M = M,
    N = N,
    weights = .validate_pred_weights(weights, N)
  )
}

#' Validate a pre-built newdata list ("modified design matrix" path)
#'
#' Advanced pathway for policy simulation: the caller supplies kernel inputs
#' directly as `list(X, alt_idx, M, [W], [weights])` (e.g., `object$data` with
#' a perturbed `X`). Dimensions and column names are validated against the
#' fitted model; `alt_idx` must use the fit-time integer codes from
#' `object$alt_mapping`.
#'
#' @param object A fitted `choicer_fit` object.
#' @param newdata List with `X`, `alt_idx`, `M`, plus `W` for MXL and an
#'   optional `weights` element.
#' @param weights Optional numeric vector of length `length(M)`; takes
#'   precedence over `newdata$weights`.
#' @returns List with the same shape as `prepare_newdata()`.
#' @noRd
validate_newdata_list <- function(object, newdata, weights = NULL) {
  missing_fields <- setdiff(c("X", "alt_idx", "M"), names(newdata))
  if (length(missing_fields) > 0) {
    stop("newdata list is missing element(s): ",
         paste(missing_fields, collapse = ", "))
  }

  x_cols <- names(object$coefficients)[object$param_map$beta]
  X <- newdata$X
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("newdata$X must be a numeric matrix.")
  }
  if (ncol(X) != length(x_cols)) {
    stop("newdata$X has ", ncol(X), " column(s); the fitted model expects ",
         length(x_cols), " (", paste(x_cols, collapse = ", "), ").")
  }
  if (!is.null(colnames(X)) && !identical(colnames(X), x_cols)) {
    stop("newdata$X column names must match the fitted covariates: ",
         paste(x_cols, collapse = ", "))
  }
  if (any(!is.finite(X))) {
    stop("newdata$X must contain only finite values.")
  }
  storage.mode(X) <- "double"

  M <- newdata$M
  if (!is.numeric(M) || anyNA(M) || any(M < 1)) {
    stop("newdata$M must contain positive integers ",
         "(alternatives per choice situation).")
  }
  M <- as.integer(M)
  if (sum(M) != nrow(X)) {
    stop("sum(newdata$M) must equal nrow(newdata$X) (got ",
         sum(M), " vs ", nrow(X), ").")
  }
  N <- length(M)

  alt_idx <- newdata$alt_idx
  J <- max(object$alt_mapping$alt_int)
  if (!is.numeric(alt_idx) || length(alt_idx) != nrow(X)) {
    stop("newdata$alt_idx must have one entry per row of newdata$X.")
  }
  alt_idx <- as.integer(alt_idx)
  if (anyNA(alt_idx) || any(alt_idx < 1L) || any(alt_idx > J)) {
    stop("newdata$alt_idx must contain integers in 1..", J,
         " (fit-time alternative codes; see object$alt_mapping).")
  }

  W <- NULL
  if (identical(object$model, "mxl")) {
    w_cols <- names(object$sW)
    W <- newdata$W
    if (is.null(W)) {
      stop("newdata list must include 'W' (random-coefficient design matrix) ",
           "for mixed logit prediction.")
    }
    if (!is.matrix(W) || !is.numeric(W)) {
      stop("newdata$W must be a numeric matrix.")
    }
    if (ncol(W) != length(w_cols)) {
      stop("newdata$W has ", ncol(W), " column(s); the fitted model expects ",
           length(w_cols), " (", paste(w_cols, collapse = ", "), ").")
    }
    if (!is.null(colnames(W)) && !identical(colnames(W), w_cols)) {
      stop("newdata$W column names must match the fitted random-coefficient ",
           "variables: ", paste(w_cols, collapse = ", "))
    }
    if (nrow(W) != nrow(X)) {
      stop("newdata$W must have the same number of rows as newdata$X.")
    }
    if (any(!is.finite(W))) {
      stop("newdata$W must contain only finite values.")
    }
    storage.mode(W) <- "double"
  }

  list(
    X = X,
    W = W,
    alt_idx = alt_idx,
    M = M,
    N = N,
    weights = .validate_pred_weights(weights %||% newdata$weights, N)
  )
}

#' Dispatch newdata to the data.frame or list pathway
#'
#' @param object A fitted `choicer_fit` object.
#' @param newdata A data.frame (long format) or list (`X`, `alt_idx`, `M`, ...).
#' @param weights Optional numeric vector of prediction weights.
#' @returns List with `X`, `W`, `alt_idx`, `M`, `N`, `weights`.
#' @noRd
resolve_predict_newdata <- function(object, newdata, weights = NULL) {
  if (is.data.frame(newdata)) {
    prepare_newdata(object, newdata, weights = weights)
  } else if (is.list(newdata)) {
    validate_newdata_list(object, newdata, weights = weights)
  } else {
    stop("'newdata' must be a data.frame in long format or a list with ",
         "X, alt_idx, and M.")
  }
}
