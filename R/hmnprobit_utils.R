# Hierarchical multinomial probit (HMNP) utilities. Phase 0 ships the data
# prep; run_hmnprobit() and the Gibbs kernel arrive with the HMNP phase of
# _plans/hierarchical_bayes_plan.md.

#' Prepare inputs for hierarchical multinomial probit estimation
#'
#' Prepares and validates panel (or cross-sectional) choice data for the
#' hierarchical Bayesian multinomial probit with iid \eqn{N(0, \sigma^2)}
#' utility shocks. The model shares its two-level random-effect structure
#' with [prepare_hmnl_data()]: respondent-level structural tastes
#' \eqn{\beta_i \sim N(b, W)} over the `covariate_cols` (normal only — the
#' probit keeps full conjugacy), and a global alternative-level effect
#' \eqn{\delta_j = z_j'\theta + \xi_j}, \eqn{\xi_j \sim N(0, \sigma_d^2)}.
#'
#' The returned structure is identical to [prepare_hmnl_data()] (both preps
#' share one internal engine), except there is no `rc_dist` field. Unlike
#' [prepare_mnp_data()], utilities are NOT differenced against a base
#' alternative: the iid-shock model works in un-differenced utility space,
#' so unbalanced choice sets are supported and the outside option is
#' implicit (its latent utility is a stochastic \eqn{N(0, \sigma^2)} draw in
#' the kernel, systematic utility 0).
#'
#' @inheritParams prepare_hmnl_data
#' @returns A list of class `c("choicer_data_hmnp", "list")` with the same
#'   components as [prepare_hmnl_data()] (minus `rc_dist`).
#' @seealso [prepare_hmnl_data()] for the component-by-component description.
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
#' dt[, choice := 0L]
#' dt[, choice := if (runif(1) < 0.8) sample(c(1L, rep(0L, J - 1))) else 0L,
#'    by = task]
#' input <- prepare_hmnp_data(dt, "task", "alt", "choice", c("x1", "x2"),
#'                            person_col = "pid")
#' input$Ti[1:5]
#' input$alt_mapping
#' @export
prepare_hmnp_data <- function(
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
  panel <- .prepare_hb_panel(
    data = data, id_col = id_col, alt_col = alt_col,
    choice_col = choice_col, covariate_cols = covariate_cols,
    person_col = person_col, alt_covariate_cols = alt_covariate_cols,
    outside_opt_label = outside_opt_label, cf_residual_col = cf_residual_col,
    include_outside_option = include_outside_option
  )
  structure(panel, class = c("choicer_data_hmnp", "list"))
}
