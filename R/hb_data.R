# Shared data-preparation machinery for the hierarchical Bayes models
# (prepare_hmnl_data / prepare_hmnp_data in R/hmnlogit_utils.R and
# R/hmnprobit_utils.R). The two preps share ~90% of their logic, factored
# into .prepare_hb_panel() below; the wrappers only add model-specific
# extras (rc_dist alignment for the HMNL) and the class tag.
#
# Design:
#   * Two-level (person, task) indexing: person_col groups choice situations
#     into respondents sharing one beta_i; person_col = NULL makes each task
#     its own respondent (Ti all 1 — the cross-sectional mode).
#   * X carries STRUCTURAL covariates only — no ASC dummy columns. The
#     alternative-level effect delta_j is indexed by alt_of_row (1..J), not
#     carried as a dense design block: the memory fix at large J.
#   * The outside option is implicit, reusing the prepare_mnl_data
#     convention (R/mnlogit_utils.R:432): physical outside rows (identified
#     by outside_opt_label) are removed, the kernels add the outside term,
#     and choice_pos = 0 encodes "outside chosen".
#   * Z (J x P) is the alternative-level mean-function design for
#     delta_j = z_j' theta + xi_j, one deduplicated row per inside
#     alternative, intercept always present.
#   * cf_residual_col (Petrin-Train control function) is appended to X as an
#     ordinary covariate; provenance is recorded in data_spec.

#' Stable internal keys for hierarchical-Bayes choice situations
#'
#' Length-prefixing makes the composite key unambiguous even when identifiers
#' contain the separator. The keys are deliberately internal: user-facing
#' output continues to use the original person/task columns.
#' @noRd
.hb_task_keys <- function(id, person = NULL) {
  id_chr <- enc2utf8(as.character(id))
  id_key <- paste0(nchar(id_chr, type = "bytes"), ":", id_chr)
  if (is.null(person)) return(id_key)

  person_chr <- enc2utf8(as.character(person))
  paste0(nchar(person_chr, type = "bytes"), ":", person_chr, "|", id_key)
}

#' Rows with a non-finite value in any of the given columns
#'
#' Column-by-column `!is.finite()`, the counterpart of `.rows_with_na()`:
#' no rows-by-columns logical matrix, and a column holding only finite values
#' is not ORed into the result.
#'
#' @param x A data frame or data.table.
#' @param cols Column positions to scan.
#' @returns Logical vector with one element per row.
#' @noRd
.rows_not_finite <- function(x, cols) {
  bad <- logical(nrow(x))
  for (j in cols) {
    finite <- is.finite(.subset2(x, j))
    if (!all(finite)) bad <- bad | !finite
  }
  bad
}

#' Drop every choice situation that has a flagged row
#'
#' An anti-join on keys with fixed names: data.table parses the strings of
#' `on =` as join conditions, so it cannot join on a task column named, say,
#' "t>0" or " task". NA keys match each other, as grouping by task does.
#'
#' @param dt The preparation's working table.
#' @param flagged Logical vector, one element per row of `dt`.
#' @param task_by Names of the person and task columns of `dt`.
#' @returns List with `dt`, the rows of the other situations in their
#'   order, and `n_tasks`, the number of situations dropped.
#' @noRd
.drop_flagged_tasks <- function(dt, flagged, task_by) {
  keys <- data.table::setDT(list(person = dt[[task_by[1L]]],
                                 task = dt[[task_by[2L]]]))
  bad <- unique(keys[flagged])
  keep <- rep(TRUE, nrow(dt))
  keep[keys[bad, on = c("person", "task"), which = TRUE]] <- FALSE
  list(dt = dt[keep], n_tasks = nrow(bad))
}

#' Shared panel preparation for the hierarchical Bayes preps
#'
#' Internal workhorse behind [prepare_hmnl_data()] and [prepare_hmnp_data()].
#' Sorts by (person, task, alternative), builds the structural design matrix
#' `X`, the alternative index `alt_of_row`, the alternative-level design `Z`,
#' the task/person index vectors, and the implicit-outside-option choice
#' encoding. See the file header for the design contract.
#'
#' @param data Data frame containing choice data.
#' @param id_col Name of the column identifying choice situations (tasks).
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the 0/1 chosen-alternative column.
#' @param covariate_cols Names of the structural covariate columns.
#' @param person_col Name of the respondent column; `NULL` makes every choice
#'   situation its own respondent.
#' @param alt_covariate_cols Names of alternative-level covariate columns
#'   (constant within each alternative) for the delta mean function.
#' @param outside_opt_label Label of physical outside-option rows to remove
#'   when the outside good is modelled implicitly.
#' @param cf_residual_col Name of a user-supplied first-stage residual column
#'   (control function), appended to `X` as an ordinary covariate.
#' @param include_outside_option Logical; model an implicit outside option
#'   with systematic utility 0.
#' @returns Unclassed list with the fields documented in
#'   [prepare_hmnl_data()].
#' @noRd
.prepare_hb_panel <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    person_col = NULL,
    alt_covariate_cols = NULL,
    outside_opt_label = NULL,
    cf_residual_col = NULL,
    include_outside_option = TRUE
) {
  ## Preliminary housekeeping --------------------------------------------------
  needed <- unique(c(person_col, id_col, alt_col, choice_col, covariate_cols,
                     alt_covariate_cols, cf_residual_col))
  # Names first: they need no data, and reading the data can stop for a
  # missing bit64.
  .check_col_names(needed, alt_col)
  # The covariates are read in place from `src` (the caller's data, never
  # modified, when possible); only the index columns are copied into `dt`,
  # and X and Z are gathered by source row.
  prep_src <- .prep_source(data, needed)
  src <- prep_src$src

  if (!is.null(cf_residual_col) && cf_residual_col %in% covariate_cols) {
    stop("`cf_residual_col` must not also appear in `covariate_cols`; ",
         "it is appended to the design matrix automatically.")
  }

  if (!all(needed %in% names(src)))
    stop("Missing columns: ",
         paste(setdiff(needed, names(src)), collapse = ", "))
  dt <- .copy_cols(src, c(person_col, id_col, alt_col, choice_col))
  dt[, .choicer_row := seq_len(.N)]

  # Endogeneity reminder: price-like covariates without a control-function
  # residual mean delta_j = z_j'theta + xi_j is exogenous only conditional on
  # Z. Informational (message, not warning) — supplying cf_residual_col is
  # the user's call. Alternative-level covariates are scanned too: a price in
  # Z is exactly the BLP case where correlation with xi_j bites.
  scan_cols <- c(covariate_cols, alt_covariate_cols)
  price_like <- grepl("price|cost|fee|tuition", scan_cols, ignore.case = TRUE)
  if (is.null(cf_residual_col) && any(price_like)) {
    message("Covariate(s) ", paste(scan_cols[price_like], collapse = ", "),
            " look like price/cost variables but no `cf_residual_col` was ",
            "supplied. If they are endogenous, consider a control-function ",
            "residual (Petrin & Train 2010).")
  }

  ## Remove outside-option rows when modelling it implicitly ------------------
  ## (mirrors prepare_mnl_data, R/mnlogit_utils.R:432)
  if (include_outside_option && !is.null(outside_opt_label)) {
    # computed outside dt[...], where a column named outside_opt_label would
    # mask the argument
    keep <- dt[[alt_col]] != outside_opt_label
    dt <- dt[keep]
    if (nrow(dt) == 0) {
      stop("No inside alternatives remain after removing outside option rows.")
    }
  }

  ## Two-level indexing: person over task --------------------------------------
  ## person_col = NULL: each choice situation is its own respondent (Ti = 1).
  ## Tasks are keyed by (person, id) so task ids only need to be unique
  ## within a respondent. (The column is read outside dt[...], where a column
  ## named id_col or person_col would take the local's place; `:=` copies it.)
  .choicer_resp <- dt[[if (is.null(person_col)) id_col else person_col]]
  dt[, .choicer_person := .choicer_resp]
  rm(.choicer_resp)
  task_by <- c(".choicer_person", id_col)

  ## Drop tasks with missing observations --------------------------------------
  ## A flagged row takes its whole task with it; the anti-join matches NA
  ## keys as grouping by task does.
  has_na <- .rows_with_na(src, prep_src$scan)[dt$.choicer_row]
  if (any(has_na)) {
    dropped <- .drop_flagged_tasks(dt, has_na, task_by)
    dt <- dropped$dt
    warning("Removed ", dropped$n_tasks,
            " choice situations containing missing values.")
  }
  rm(has_na)
  if (nrow(dt) == 0) {
    stop("All choice situations removed due to missing values.")
  }

  ## Sanity checks -------------------------------------------------------------

  ## Covariates (incl. cf residual and alt-level covariates) must be numeric
  x_cols <- c(covariate_cols, cf_residual_col)
  num_cols <- unique(c(x_cols, alt_covariate_cols))
  num_pos <- match(num_cols, names(src))
  if (!all(vapply(num_pos, function(j) is.numeric(.subset2(src, j)), NA)))
    stop("All covariates must be numeric.")

  ## Non-finite covariate values (Inf/-Inf/NaN) are as fatal as NAs: same
  ## graceful task-drop path, instead of failing the terminal
  ## stopifnot(all(is.finite(X))) with an unactionable assertion.
  has_bad <- .rows_not_finite(src, num_pos)[dt$.choicer_row]
  if (any(has_bad)) {
    dropped <- .drop_flagged_tasks(dt, has_bad, task_by)
    dt <- dropped$dt
    warning("Removed ", dropped$n_tasks,
            " choice situations containing non-finite covariate values.")
  }
  rm(has_bad)
  if (nrow(dt) == 0) {
    stop("All choice situations removed due to non-finite covariate values.")
  }

  ## choice column must be 0/1 with the outside-option convention of
  ## prepare_mnl_data (R/mnlogit_utils.R:459-472): exactly one '1' per task,
  ## or at most one when an all-zeros task means "outside chosen".
  if (!.is_zero_one(dt[[choice_col]]))
    stop("`", choice_col, "` must contain only 0 and 1.")

  n_chosen <- .n_chosen_by(dt, choice_col, task_by)
  if (include_outside_option == FALSE && any(n_chosen != 1)) {
    stop("Each ", id_col, " must have exactly one chosen alternative (one '1' in ",
         choice_col, ").")
  }
  if (include_outside_option && any(n_chosen > 1)) {
    stop("Each ", id_col, " must have at most one chosen alternative (one '1' in ",
         choice_col, "). A choice situation with no explicit choice is ",
         "assumed to be outside option.")
  }

  ## Create integer alternative codes (inside alternatives, 1..J) --------------
  levels <- sort(unique(dt[[alt_col]]))
  dt <- .code_alternatives(dt, alt_col, levels)
  J <- length(levels)

  ## An alternative may appear at most once per choice situation: the kernels'
  ## incremental delta-phase denominator updates assume each (task, j) pair is
  ## a single row, and a duplicate would silently corrupt them. One pass over
  ## the (person, task, alternative) keys; the situations are counted only
  ## when a key repeats.
  keys <- data.table::setDT(list(person = dt$.choicer_person,
                                 task = dt[[id_col]],
                                 alt = dt$.choicer_alt_int))
  if (anyDuplicated(keys)) {
    n_dup <- nrow(unique(keys[duplicated(keys), c("person", "task")]))
    stop("Each alternative may appear at most once per choice situation; ",
         n_dup, " choice situation(s) contain duplicated alternatives.")
  }
  rm(keys)

  ## Order rows ----------------------------------------------------------------
  ##   between persons          : ascending person
  ##   within person, between tasks: ascending task id
  ##   within task              : ascending alternative code
  ## This sort is the single source of truth for every downstream index
  ## (alt_of_row, choice_pos, the kernel CSR offsets).
  data.table::setorderv(dt, c(".choicer_person", id_col, ".choicer_alt_int"))

  ## Position within the task, and the task's number (1..n_tasks in sorted
  ## order, where each task's rows are contiguous), computed outside dt[...],
  ## where a column named task_by would take the local's place.
  .choicer_idx <- data.table::rowidv(dt, cols = task_by)
  .choicer_task <- cumsum(.choicer_idx == 1L)
  dt[, `:=`(.choicer_idx_in_group = .choicer_idx,
            .choicer_task_idx = .choicer_task)]
  rm(.choicer_idx, .choicer_task)

  # Retain sorted task identities so welfare counterfactuals can match the
  # baseline and policy states by identity rather than silently by position.
  task_identity <- unique(dt[, ..task_by])
  task_keys <- .hb_task_keys(
    task_identity[[id_col]],
    if (!is.null(person_col)) task_identity[[".choicer_person"]]
  )

  ## Task-constant covariates ---------------------------------------------------
  ## A covariate with no within-task variation is unidentified WITHOUT an
  ## outside good (it cancels from every softmax/utility contrast and
  ## flattens the pooled MLE) — dropped with a warning. WITH a first-class
  ## outside good it shifts all inside utilities relative to the outside and
  ## is genuinely identified — kept, with an informational message.
  ## A column at a time from the source: a covariate is constant within every
  ## task when each row equals its task's first row (the values are finite
  ## here). Rows are sorted by task and .choicer_task_idx numbers the tasks in
  ## that order, so first_src is the source row of each row's task's first row.
  first_src <-
    dt$.choicer_row[dt$.choicer_idx_in_group == 1L][dt$.choicer_task_idx]
  task_const <- vapply(x_cols, function(cc) {
    v <- .subset2(src, match(cc, names(src)))
    all(v[dt$.choicer_row] == v[first_src])
  }, logical(1L))
  rm(first_src)
  dropped_task_const <- character(0)
  if (any(task_const)) {
    const_cols <- x_cols[task_const]
    if (include_outside_option) {
      message("Covariate(s) constant within every choice situation kept: ",
              paste(const_cols, collapse = ", "),
              " (identified relative to the outside option).")
    } else {
      warning("Covariate(s) constant within every choice situation are not ",
              "identified without an outside option and were dropped: ",
              paste(const_cols, collapse = ", "), call. = FALSE)
      if (!is.null(cf_residual_col) && cf_residual_col %in% const_cols) {
        # Losing the control function is a substantive modelling change, not
        # just a design-matrix cleanup — call it out by name.
        warning("The control-function residual `", cf_residual_col, "` was ",
                "among the dropped task-constant columns: the endogeneity ",
                "correction is NOT active in this fit.", call. = FALSE)
      }
      dropped_task_const <- const_cols
      x_cols <- setdiff(x_cols, const_cols)
      if (length(x_cols) == 0) {
        stop("No covariates remain after dropping columns constant within ",
             "every choice situation.")
      }
    }
  }

  ## Build objects -------------------------------------------------------------
  ## Structural design matrix: covariates only, cf residual (if any) last.
  ## NO ASC dummies — delta_j is indexed by alt_of_row, never carried in X.
  warn_once <- .int64_warn_once()  # a column can be in both X and Z
  X <- warn_once(.gather_matrix(src, x_cols, dt$.choicer_row,  # rows x K_struct
                                "The design matrix X"))
  X_res <- check_collinearity(X)
  X <- X_res$mat
  dropped_vars <- c(dropped_task_const, X_res$dropped)
  K_struct <- ncol(X)

  ## Alternative index per row (1..J); doubles as alt_idx for the pooled-MLE
  ## init, which reuses the identical X/M/choice_pos through the existing
  ## frequentist kernels.
  alt_of_row <- as.integer(dt$.choicer_alt_int)

  ## M[t] - # inside alternatives per task (with the implicit outside the
  ## effective choice set is M + 1)
  M <- .n_rows_by(dt, task_by)
  n_tasks <- length(M)
  if (!include_outside_option && any(M < 2)) {
    stop("Each choice situation must contain at least 2 alternatives when ",
         "include_outside_option = FALSE.")
  }

  ## choice_pos[t] - 1-based index of the chosen row *within* its task;
  ## 0 = outside option chosen (only with include_outside_option = TRUE)
  choice_pos <- integer(n_tasks)
  chosen <- dt[[choice_col]] == 1
  choice_pos[dt$.choicer_task_idx[chosen]] <- dt$.choicer_idx_in_group[chosen]
  rm(chosen)

  ## Person-level indexing: Ti tasks per person, in sorted person order
  person_task <- unique(dt[, .(.choicer_person, .choicer_task_idx)])
  Ti <- person_task[, .N, by = .choicer_person][["N"]]
  person_ids <- unique(person_task$.choicer_person)
  N_persons <- length(person_ids)

  ## Alternative-level design Z (J x P) ----------------------------------------
  z_res <- warn_once(.resolve_alt_covariates(
    src, dt$.choicer_row, dt$.choicer_alt_int, alt_covariate_cols, levels))
  Z <- z_res$Z
  P <- ncol(Z)
  dt[, .choicer_row := NULL]

  ## Alternatives summary (mirrors prepare_mnl_data) ---------------------------
  ## One inside-alternative aggregation; the outside branch only prepends its
  ## synthetic alt_int = 0 row.
  alt_mapping <- .alt_counts(dt, alt_col, choice_col)
  if (include_outside_option) {
    outside_alt_mapping <- data.table::data.table(
      alt_int = 0L, N_OBS = n_tasks, N_CHOICES = sum(choice_pos == 0L)
    )
    outside_alt_mapping[[alt_col]] <- outside_opt_label %||% NA
    alt_mapping <- list(outside_alt_mapping, alt_mapping) |>
      data.table::rbindlist(use.names = TRUE, fill = TRUE)
    data.table::setcolorder(alt_mapping,
                            c("alt_int", alt_col, "N_OBS", "N_CHOICES"))
  }
  alt_mapping[, `:=`(
    TAKE_RATE = N_CHOICES / N_OBS,
    MKT_SHARE = N_CHOICES / sum(N_CHOICES)
  )]

  ## Parameter index map (robust to collinearity/task-constant drops)
  param_map <- list(
    beta  = stats::setNames(seq_len(K_struct), colnames(X)),
    theta = stats::setNames(seq_len(P), colnames(Z))
  )

  ## Final validity checks -----------------------------------------------------
  stopifnot(
    length(alt_of_row) == nrow(X),
    length(choice_pos) == n_tasks,
    length(M)          == n_tasks,
    sum(Ti)            == n_tasks,
    length(person_ids) == length(Ti),
    all(choice_pos >= 0L & choice_pos <= M),
    nrow(Z)            == J,
    all(is.finite(X)),
    all(is.finite(Z))
  )

  ## return output -------------------------------------------------------------
  list(
    X           = X,
    alt_of_row  = alt_of_row,
    alt_idx     = alt_of_row,          # alias for the pooled-MLE init kernels
    Z           = Z,
    M           = M,
    choice_pos  = choice_pos,
    Ti          = Ti,
    person_ids  = person_ids,
    task_keys   = task_keys,
    N_persons   = N_persons,
    n_tasks     = n_tasks,
    J           = as.integer(J),
    K_struct    = K_struct,
    P           = P,
    include_outside_option = include_outside_option,
    alt_mapping = alt_mapping,
    param_map   = param_map,
    dropped_cols   = if (length(dropped_vars) > 0) dropped_vars else NULL,
    dropped_z_cols = if (length(z_res$dropped) > 0) z_res$dropped else NULL,
    data_spec = list(
      id_col = id_col,
      alt_col = alt_col,
      choice_col = choice_col,
      covariate_cols = covariate_cols,
      person_col = person_col,
      alt_covariate_cols = alt_covariate_cols,
      outside_opt_label = outside_opt_label,
      cf_residual_col = cf_residual_col,
      include_outside_option = include_outside_option
    )
  )
}

#' Build the alternative-level mean-function design Z
#'
#' Deduplicates `alt_covariate_cols` to one row per inside alternative
#' (validating that each column is constant within its alternative), prepends
#' an always-present intercept column, drops non-intercept columns that are
#' constant across alternatives (identified only through the intercept), and
#' removes any remaining collinear columns. With `alt_covariate_cols = NULL`
#' the design is intercept-only (P = 1), so theta_0 is the common inside-good
#' level relative to the outside option.
#'
#' @param src Data frame holding the alternative covariate columns (see
#'   `.prep_source()`).
#' @param rows Source row of each prepared row, in prepared (sorted) order.
#' @param alt_int Integer alternative code (`1..J`) of each prepared row.
#' @param alt_covariate_cols Names of alternative-level covariate columns, or
#'   `NULL` for an intercept-only design.
#' @param levels Sorted vector of inside-alternative labels (length J).
#' @returns List with `Z` (J x P matrix, intercept first) and `dropped`
#'   (names of dropped Z columns, possibly empty).
#' @noRd
.resolve_alt_covariates <- function(src, rows, alt_int, alt_covariate_cols,
                                    levels) {
  J <- length(levels)
  if (is.null(alt_covariate_cols)) {
    Z <- matrix(1, nrow = J, ncol = 1,
                dimnames = list(NULL, "(Intercept)"))
    return(list(Z = Z, dropped = character(0)))
  }

  ## Constant-within-alternative validation: z_j is a property of the
  ## alternative, so any within-alternative variation is a data error. Each
  ## column is read on the prepared rows, one at a time, beside the
  ## alternative codes; its value at an alternative's first row is z_j.
  pos <- match(alt_covariate_cols, names(src))
  constant <- logical(length(pos))
  z_first <- vector("list", length(pos))
  for (k in seq_along(pos)) {
    pairs <- data.table::setDT(list(alt = alt_int,
                                    value = .subset2(src, pos[k])[rows]))
    constant[k] <- all(.n_distinct_by(pairs, "value", "alt") == 1L)
    # as.matrix() below would read integer64 values as raw bits
    z_first[[k]] <- .int64_to_double(
      .first_by(pairs, "value", "alt", seq_len(J)),
      paste0("Column '", alt_covariate_cols[k], "'"))
  }
  bad <- alt_covariate_cols[!constant]
  if (length(bad) > 0) {
    stop("`alt_covariate_cols` must be constant within each alternative: ",
         paste(bad, collapse = ", "))
  }

  ## One row per alternative, in alt_int (= sorted label) order.
  Zmat <- as.matrix(data.table::setDT(stats::setNames(z_first,
                                                      alt_covariate_cols)))

  ## Non-intercept columns constant ACROSS alternatives carry no information
  ## beyond the intercept — dropped with a message (the intercept itself is
  ## always kept: theta_0 is identified against the outside good).
  const_across <- vapply(
    seq_len(ncol(Zmat)),
    function(k) max(Zmat[, k]) - min(Zmat[, k]) == 0,
    logical(1L)
  )
  dropped <- character(0)
  if (any(const_across)) {
    dropped <- colnames(Zmat)[const_across]
    message("Alternative-level covariate(s) constant across alternatives ",
            "dropped from Z (only the intercept identifies a common level): ",
            paste(dropped, collapse = ", "))
    Zmat <- Zmat[, !const_across, drop = FALSE]
  }

  Z <- cbind("(Intercept)" = rep(1, J), Zmat)
  Z_res <- check_collinearity(Z)
  Z <- Z_res$mat
  if (!("(Intercept)" %in% colnames(Z))) {
    stop("Internal error: the Z intercept column was dropped as collinear.")
  }
  dropped <- c(dropped, Z_res$dropped)

  list(Z = Z, dropped = dropped)
}
