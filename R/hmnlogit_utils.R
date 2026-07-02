# Hierarchical multinomial logit (HMNL) utilities. Phase 0 ships the data
# prep; run_hmnlogit() and the Gibbs kernel arrive with the HMNL phase of
# _plans/hierarchical_bayes_plan.md.

#' Prepare inputs for hierarchical multinomial logit estimation
#'
#' Prepares and validates panel (or cross-sectional) choice data for the
#' hierarchical Bayesian multinomial logit. The model has two random-effect
#' levels: respondent-level structural tastes \eqn{\beta_i \sim N(b, W)}
#' over the `covariate_cols`, and a global alternative-level effect
#' \eqn{\delta_j = z_j'\theta + \xi_j}, \eqn{\xi_j \sim N(0, \sigma_d^2)},
#' with mean-function design \eqn{z_j} built from `alt_covariate_cols`.
#'
#' \strong{Structure.} The design matrix `X` carries structural covariates
#' only — no alternative-specific-constant dummies. The alternative effect
#' \eqn{\delta_j} is indexed by `alt_of_row` (integer codes `1..J`), so
#' memory and compute scale with the number of rows, not with `J` extra
#' design columns.
#'
#' \strong{Outside option.} With `include_outside_option = TRUE` (the
#' default) the outside good is modelled implicitly, following the
#' [prepare_mnl_data()] convention: physical outside rows (identified by
#' `outside_opt_label`) are removed, the estimation kernels add the outside
#' term (systematic utility 0), and a choice situation whose inside rows are
#' all `0` in `choice_col` is coded as "outside chosen" (`choice_pos = 0`).
#' The outside option anchors the location of \eqn{\delta} (mean utility
#' relative to the outside good).
#'
#' \strong{Cross-section vs panel.} `person_col` groups choice situations
#' into respondents sharing one \eqn{\beta_i}. With `person_col = NULL`
#' (default) every choice situation is its own respondent (`Ti` all 1) —
#' the cross-sectional random-coefficients mode.
#'
#' \strong{Control function.} `cf_residual_col` (a user-supplied first-stage
#' residual, Petrin & Train 2010) is appended to `X` as an ordinary
#' covariate; its provenance is recorded in `data_spec`. The first stage is
#' NOT run here — supplying a valid residual is the user's responsibility.
#'
#' @param data Data frame containing choice data.
#' @param id_col Name of the column identifying choice situations (tasks).
#'   Task ids only need to be unique within a respondent.
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating the chosen alternative
#'   (1 = chosen, 0 = not chosen).
#' @param covariate_cols Vector of names of structural covariate columns
#'   (the random-coefficient dimensions).
#' @param person_col Name of the respondent column grouping choice
#'   situations. `NULL` (default) makes each choice situation its own
#'   respondent.
#' @param alt_covariate_cols Names of alternative-level covariate columns
#'   (constant within each alternative) forming the \eqn{\delta} mean
#'   function. `NULL` (default) gives an intercept-only design (P = 1).
#' @param outside_opt_label Label of physical outside-option rows, removed
#'   when `include_outside_option = TRUE` (the outside good is implicit).
#' @param cf_residual_col Name of a first-stage residual column (control
#'   function for an endogenous covariate), appended to `X`. Default `NULL`.
#' @param include_outside_option Logical; if `TRUE` (default) an implicit
#'   outside option with systematic utility 0 is part of every choice set.
#' @param rc_dist Integer vector, one entry per column of `covariate_cols`:
#'   `0` for a normal random coefficient, `1` for log-normal (the
#'   coefficient enters utility as `exp(beta_ik)`; hierarchy normal on the
#'   log scale). Default `NULL` is all-normal. Automatically aligned through
#'   dropped columns; a `cf_residual_col` coordinate is always normal.
#' @returns A list of class `c("choicer_data_hmnl", "list")` containing:
#'   \itemize{
#'     \item `X`: Structural design matrix (total_rows x K_struct), no ASC
#'       columns; `cf_residual_col` last when supplied.
#'     \item `alt_of_row`: Integer alternative code per row (`1..J`).
#'     \item `alt_idx`: Alias of `alt_of_row` for the pooled-MLE init.
#'     \item `Z`: Alternative-level design (J x P), intercept first.
#'     \item `M`: Inside alternatives per choice situation.
#'     \item `choice_pos`: 1-based within-task position of the chosen row;
#'       `0` = outside option chosen.
#'     \item `Ti`: Choice situations per respondent.
#'     \item `person_ids`, `N_persons`, `n_tasks`, `J`, `K_struct`, `P`.
#'     \item `include_outside_option`: Logical flag.
#'     \item `alt_mapping`: Data.table mapping alternatives to summary
#'       statistics (outside option is `alt_int = 0`).
#'     \item `param_map`: Named list of index vectors (`beta`, `theta`),
#'       robust to collinearity drops.
#'     \item `rc_dist`: Integer vector aligned with the columns of `X`.
#'     \item `dropped_cols`, `dropped_z_cols`: Dropped column names, if any.
#'     \item `data_spec`: Column-name metadata (incl. `person_col`,
#'       `outside_opt_label`, `cf_residual_col`, `alt_covariate_cols`).
#'   }
#' @seealso [prepare_hmnp_data()] for the hierarchical probit counterpart.
#' @examples
#' library(data.table)
#' set.seed(42)
#' N <- 20; T <- 3; J <- 4
#' dt <- data.table(
#'   pid  = rep(1:N, each = T * J),
#'   task = rep(seq_len(N * T), each = J),
#'   alt  = rep(1:J, N * T)
#' )
#' dt[, `:=`(x1 = rnorm(.N), x2 = runif(.N, -1, 1))]
#' dt[, quality := alt / J, by = alt]        # alternative-level covariate
#' dt[, choice := 0L]
#' # leave some tasks all-zero: outside option chosen
#' dt[, choice := if (runif(1) < 0.8) sample(c(1L, rep(0L, J - 1))) else 0L,
#'    by = task]
#' input <- prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "x2"),
#'                            person_col = "pid",
#'                            alt_covariate_cols = "quality")
#' str(input$Z)
#' input$alt_mapping
#' @export
prepare_hmnl_data <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    person_col = NULL,
    alt_covariate_cols = NULL,
    outside_opt_label = NULL,
    cf_residual_col = NULL,
    include_outside_option = TRUE,
    rc_dist = NULL
) {
  ## rc_dist is validated against the user-supplied covariates BEFORE the
  ## shared prep so error messages refer to the input, then aligned to the
  ## surviving X columns afterwards (collinearity / task-constant drops).
  K_user <- length(covariate_cols)
  if (is.null(rc_dist)) rc_dist <- rep(0L, K_user)
  rc_dist <- as.integer(rc_dist)
  if (length(rc_dist) != K_user) {
    stop("`rc_dist` must have length ", K_user,
         " (one entry per column of `covariate_cols`).")
  }
  if (!all(rc_dist %in% c(0L, 1L))) {
    stop("`rc_dist` entries must be 0 (normal) or 1 (log-normal).")
  }

  panel <- .prepare_hb_panel(
    data = data, id_col = id_col, alt_col = alt_col,
    choice_col = choice_col, covariate_cols = covariate_cols,
    person_col = person_col, alt_covariate_cols = alt_covariate_cols,
    outside_opt_label = outside_opt_label, cf_residual_col = cf_residual_col,
    include_outside_option = include_outside_option
  )

  ## Align rc_dist through drops: X columns are the kept structural
  ## covariates (in covariate_cols order) plus, last, the cf residual —
  ## whose coefficient is always a normal coordinate.
  rc_map <- stats::setNames(rc_dist, covariate_cols)
  kept <- colnames(panel$X)
  rc_aligned <- integer(length(kept))
  names(rc_aligned) <- kept
  is_user_cov <- kept %in% covariate_cols
  rc_aligned[is_user_cov] <- rc_map[kept[is_user_cov]]
  panel$rc_dist <- rc_aligned

  structure(panel, class = c("choicer_data_hmnl", "list"))
}
