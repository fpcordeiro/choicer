# Shared data-preparation machinery for the hierarchical Bayes models
# (prepare_hmnl_data / prepare_hmnp_data in R/hmnlogit_utils.R and
# R/hmnprobit_utils.R). The two preps share ~90% of their logic, factored
# into .prepare_hb_panel() below; the wrappers only add model-specific
# extras (rc_dist alignment for the HMNL) and the class tag.
#
# Design (see _plans/hierarchical_bayes_plan.md, "Data prep"):
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
  dt <- data.table::as.data.table(data)[]

  if (!is.null(cf_residual_col) && cf_residual_col %in% covariate_cols) {
    stop("`cf_residual_col` must not also appear in `covariate_cols`; ",
         "it is appended to the design matrix automatically.")
  }

  needed <- unique(c(person_col, id_col, alt_col, choice_col, covariate_cols,
                     alt_covariate_cols, cf_residual_col))
  if (!all(needed %in% names(dt)))
    stop("Missing columns: ",
         paste(setdiff(needed, names(dt)), collapse = ", "))

  # Drop non-relevant variables
  vars_to_drop <- setdiff(names(dt), needed)
  if (length(vars_to_drop) > 0) {
    dt[, (vars_to_drop) := NULL]
  }

  # Endogeneity reminder: price-like covariates without a control-function
  # residual mean delta_j = z_j'theta + xi_j is exogenous only conditional on
  # Z. Informational (message, not warning) — supplying cf_residual_col is
  # the user's call.
  price_like <- grepl("price|cost|fee|tuition", covariate_cols,
                      ignore.case = TRUE)
  if (is.null(cf_residual_col) && any(price_like)) {
    message("Covariate(s) ", paste(covariate_cols[price_like], collapse = ", "),
            " look like price/cost variables but no `cf_residual_col` was ",
            "supplied. If they are endogenous, consider a control-function ",
            "residual (Petrin & Train 2010).")
  }

  ## Remove outside-option rows when modelling it implicitly ------------------
  ## (mirrors prepare_mnl_data, R/mnlogit_utils.R:432)
  if (include_outside_option && !is.null(outside_opt_label)) {
    dt <- dt[get(alt_col) != outside_opt_label]
    if (nrow(dt) == 0) {
      stop("No inside alternatives remain after removing outside option rows.")
    }
  }

  ## Two-level indexing: person over task --------------------------------------
  ## person_col = NULL: each choice situation is its own respondent (Ti = 1).
  ## Tasks are keyed by (person, id) so task ids only need to be unique
  ## within a respondent.
  if (is.null(person_col)) {
    dt[, HB_PERSON := get(id_col)]
  } else {
    dt[, HB_PERSON := get(person_col)]
  }
  task_by <- c("HB_PERSON", id_col)

  ## Drop tasks with missing observations --------------------------------------
  dt[, HAS_NA := rowSums(is.na(.SD)) > 0]
  dt[, TASK_HAS_NA := any(HAS_NA), by = task_by]
  n_bad_tasks <- nrow(unique(dt[TASK_HAS_NA == TRUE, ..task_by]))
  if (n_bad_tasks > 0) {
    dt <- dt[TASK_HAS_NA == FALSE]
    warning("Removed ", n_bad_tasks,
            " choice situations containing missing values.")
  }
  if (nrow(dt) == 0) {
    stop("All choice situations removed due to missing values.")
  }
  dt[, c("HAS_NA", "TASK_HAS_NA") := NULL]

  ## Sanity checks -------------------------------------------------------------

  ## Covariates (incl. cf residual and alt-level covariates) must be numeric
  x_cols <- c(covariate_cols, cf_residual_col)
  num_cols <- unique(c(x_cols, alt_covariate_cols))
  if (!all(vapply(dt[, ..num_cols], is.numeric, logical(1L))))
    stop("All covariates must be numeric.")

  ## choice column must be 0/1 with the outside-option convention of
  ## prepare_mnl_data (R/mnlogit_utils.R:459-472): exactly one '1' per task,
  ## or at most one when an all-zeros task means "outside chosen".
  bad_choice <- dt[[choice_col]] %in% c(0, 1) == FALSE
  if (any(bad_choice))
    stop("`", choice_col, "` must contain only 0 and 1.")

  by_task <- dt[, .(chosen = sum(get(choice_col))), by = task_by]
  if (include_outside_option == FALSE && any(by_task$chosen != 1)) {
    stop("Each ", id_col, " must have exactly one chosen alternative (one '1' in ",
         choice_col, ").")
  }
  if (include_outside_option && any(by_task$chosen > 1)) {
    stop("Each ", id_col, " must have at most one chosen alternative (one '1' in ",
         choice_col, "). A choice situation with no explicit choice is ",
         "assumed to be outside option.")
  }

  ## Create integer alternative codes (inside alternatives, 1..J) --------------
  levels <- sort(unique(dt[[alt_col]]))
  dt[, alt_int := as.integer(factor(get(alt_col), levels = levels))]
  J <- length(levels)

  ## Order rows ----------------------------------------------------------------
  ##   between persons          : ascending person
  ##   within person, between tasks: ascending task id
  ##   within task              : ascending alternative code
  ## This sort is the single source of truth for every downstream index
  ## (alt_of_row, choice_pos, the kernel CSR offsets).
  data.table::setorderv(dt, c("HB_PERSON", id_col, "alt_int"))

  dt[, idx_in_group := seq_len(.N), by = task_by]
  dt[, task_idx := .GRP, by = task_by]      # 1..n_tasks in sorted order

  ## Task-constant covariates ---------------------------------------------------
  ## A covariate with no within-task variation is unidentified WITHOUT an
  ## outside good (it cancels from every softmax/utility contrast and
  ## flattens the pooled MLE) — dropped with a warning. WITH a first-class
  ## outside good it shifts all inside utilities relative to the outside and
  ## is genuinely identified — kept, with an informational message.
  rng_by_task <- dt[, lapply(.SD, function(v) max(v) - min(v)),
                    by = task_by, .SDcols = x_cols]
  task_const <- vapply(x_cols, function(cc) all(rng_by_task[[cc]] == 0),
                       logical(1L))
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
  X <- as.matrix(dt[, ..x_cols])                       # total_rows x K_struct
  X_res <- check_collinearity(X)
  X <- X_res$mat
  dropped_vars <- c(dropped_task_const, X_res$dropped)
  K_struct <- ncol(X)

  ## Alternative index per row (1..J); doubles as alt_idx for the pooled-MLE
  ## init, which reuses the identical X/M/choice_pos through the existing
  ## frequentist kernels.
  alt_of_row <- as.integer(dt$alt_int)

  ## M[t] - # inside alternatives per task (with the implicit outside the
  ## effective choice set is M + 1)
  M <- dt[, .N, by = task_by][["N"]]
  n_tasks <- length(M)
  if (!include_outside_option && any(M < 2)) {
    stop("Each choice situation must contain at least 2 alternatives when ",
         "include_outside_option = FALSE.")
  }

  ## choice_pos[t] - 1-based index of the chosen row *within* its task;
  ## 0 = outside option chosen (only with include_outside_option = TRUE)
  choice_pos <- integer(n_tasks)
  chosen_dt <- dt[get(choice_col) == 1, .(task_idx, pos = idx_in_group)]
  choice_pos[chosen_dt$task_idx] <- chosen_dt$pos

  ## Person-level indexing: Ti tasks per person, in sorted person order
  person_task <- unique(dt[, .(HB_PERSON, task_idx)])
  Ti <- person_task[, .N, by = HB_PERSON][["N"]]
  person_ids <- unique(person_task$HB_PERSON)
  N_persons <- length(person_ids)

  ## Alternative-level design Z (J x P) ----------------------------------------
  z_res <- .resolve_alt_covariates(dt, alt_covariate_cols, levels)
  Z <- z_res$Z
  P <- ncol(Z)

  ## Alternatives summary (mirrors prepare_mnl_data) ---------------------------
  if (include_outside_option) {
    inside_alt_mapping <- dt[
      , .(N_OBS = .N, N_CHOICES = sum(get(choice_col))),
      keyby = c("alt_int", alt_col)
    ]
    outside_alt_mapping <- data.table::data.table(
      alt_int = 0L, N_OBS = n_tasks, N_CHOICES = sum(choice_pos == 0L)
    )
    outside_alt_mapping[[alt_col]] <- outside_opt_label %||% NA
    alt_mapping <- list(outside_alt_mapping, inside_alt_mapping) |>
      data.table::rbindlist(use.names = TRUE, fill = TRUE)
    data.table::setcolorder(alt_mapping,
                            c("alt_int", alt_col, "N_OBS", "N_CHOICES"))
  } else {
    alt_mapping <- dt[
      , .(N_OBS = .N, N_CHOICES = sum(get(choice_col))),
      keyby = c("alt_int", alt_col)
    ]
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
#' @param dt Sorted prep data.table carrying `alt_int` and the alternative
#'   covariate columns.
#' @param alt_covariate_cols Names of alternative-level covariate columns, or
#'   `NULL` for an intercept-only design.
#' @param levels Sorted vector of inside-alternative labels (length J).
#' @returns List with `Z` (J x P matrix, intercept first) and `dropped`
#'   (names of dropped Z columns, possibly empty).
#' @noRd
.resolve_alt_covariates <- function(dt, alt_covariate_cols, levels) {
  J <- length(levels)
  if (is.null(alt_covariate_cols)) {
    Z <- matrix(1, nrow = J, ncol = 1,
                dimnames = list(NULL, "(Intercept)"))
    return(list(Z = Z, dropped = character(0)))
  }

  ## Constant-within-alternative validation: z_j is a property of the
  ## alternative, so any within-alternative variation is a data error.
  nuniq <- dt[, lapply(.SD, data.table::uniqueN),
              by = alt_int, .SDcols = alt_covariate_cols]
  bad <- alt_covariate_cols[
    vapply(alt_covariate_cols, function(cc) any(nuniq[[cc]] != 1L),
           logical(1L))
  ]
  if (length(bad) > 0) {
    stop("`alt_covariate_cols` must be constant within each alternative: ",
         paste(bad, collapse = ", "))
  }

  ## One row per alternative, in alt_int (= sorted label) order.
  zdt <- dt[, lapply(.SD, function(v) v[1L]),
            keyby = alt_int, .SDcols = alt_covariate_cols]
  Zmat <- as.matrix(zdt[, ..alt_covariate_cols])

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
