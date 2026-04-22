# Parameter recovery diagnostics.
# Consumes `choicer_sim` from R/simulation.R and fitted `choicer_fit` objects.

# Internal helpers =============================================================

.safe_div <- function(num, denom) {
  ifelse(denom == 0 | !is.finite(denom), NA_real_, num / denom)
}

# Convergence dispatch ---------------------------------------------------------
# Different optimizers use incompatible convergence codes. nloptr returns
# positive integers for success (1-4 = tolerance-based stop), stats::optim
# returns 0 on success. Custom optimizers are treated permissively.
.is_converged <- function(fit) {
  code <- fit$convergence
  if (is.null(code) || length(code) != 1L || !is.finite(code)) return(FALSE)
  nm <- fit$optimizer$name %||% "nloptr"
  switch(nm,
    nloptr = isTRUE(as.integer(code) %in% c(1L, 2L, 3L, 4L)),
    optim  = isTRUE(as.integer(code) == 0L),
    # custom: be permissive (both conventions are seen in the wild)
    isTRUE(as.integer(code) %in% c(0L, 1L, 2L, 3L, 4L))
  )
}

.build_recovery_rows <- function(group, idx, truth_vec, coefs, se, z) {
  if (length(idx) != length(truth_vec)) {
    stop(sprintf(
      "recovery_table: mismatched length in '%s' block (fit has %d params, truth has %d). Check the DGP and fit settings.",
      group, length(idx), length(truth_vec)
    ))
  }
  est    <- unname(coefs[idx])
  se_i   <- unname(se[idx])
  pnames <- names(coefs)[idx]
  if (is.null(pnames)) pnames <- paste0(group, "_", seq_along(idx))
  data.table::data.table(
    parameter    = pnames,
    group        = group,
    true         = as.numeric(truth_vec),
    estimate     = est,
    se           = se_i,
    bias         = est - as.numeric(truth_vec),
    rel_bias_pct = 100 * .safe_div(est - as.numeric(truth_vec), as.numeric(truth_vec)),
    z_vs_true    = .safe_div(est - as.numeric(truth_vec), se_i),
    lower_ci     = est - z * se_i,
    upper_ci     = est + z * se_i,
    covers       = (as.numeric(truth_vec) >= est - z * se_i) &
                   (as.numeric(truth_vec) <= est + z * se_i)
  )
}

# Drop the normalized ASC from truth$delta when the estimator has fixed the
# first inside alternative's ASC to zero for identification. Detection is
# purely structural: the ASC block is shorter than truth$delta by exactly one,
# there is no outside option baked into the fit, and the data contains no
# "outside" alt (alt == 0 / j == 0 / etc.) that would have absorbed the drop.
.maybe_drop_normalized_asc <- function(fit, truth_delta) {
  pm <- fit$param_map
  n_asc <- length(pm$asc %||% integer(0))
  n_delta <- length(truth_delta)
  if (n_asc == n_delta) return(truth_delta)
  if (n_asc == n_delta - 1L && !isTRUE(fit$include_outside_option)) {
    # `run_*logit()` fixes the first inside alternative's ASC to zero when
    # there is no outside option. Truth's delta therefore needs its first
    # entry dropped to match what the estimator returns.
    return(truth_delta[-1L])
  }
  truth_delta
}

.align_truth_to_coefs <- function(truth, coefs, param_map, fit = NULL) {
  # Returns list of recovery-row blocks (one per group). Blocks present in both
  # `param_map` and `truth` are included; others are silently skipped (for
  # example a fit without ASCs, or a DGP without mu).
  mapping <- list(
    beta   = list(pm = "beta",   truth = "beta"),
    mu     = list(pm = "mu",     truth = "mu"),
    sigma  = list(pm = "sigma",  truth = "L_params"),
    lambda = list(pm = "lambda", truth = "lambdas"),
    asc    = list(pm = "asc",    truth = "delta")
  )
  lapply(names(mapping), function(group) {
    m <- mapping[[group]]
    idx <- param_map[[m$pm]]
    truth_vec <- truth[[m$truth]]
    if (is.null(idx) || is.null(truth_vec)) return(NULL)
    if (group == "asc" && !is.null(fit)) {
      truth_vec <- .maybe_drop_normalized_asc(fit, truth_vec)
    }
    list(group = group, idx = idx, truth_vec = as.numeric(truth_vec))
  })
}

# recovery_table() S3 generic ==================================================

#' Parameter recovery table
#'
#' Compares fitted coefficients to a set of true parameter values on the same
#' scale as the estimator's internal parameterization. Returns one row per
#' estimated parameter with true value, estimate, standard error, bias,
#' relative bias (%), z-score against the truth, Wald CI, and a coverage
#' indicator.
#'
#' For MXL fits the `sigma` block compares the raw Cholesky parameters
#' (`L_params`), not the reconstructed covariance matrix. For log-normal
#' random-coefficient means the raw `mu` estimate is compared directly; callers
#' who want recovery on the DGP scale (`exp(mu)`) should transform both sides
#' before calling.
#'
#' When the estimator has normalized the first inside alternative's ASC to
#' zero (which happens for MNL/MXL with `include_outside_option = FALSE` and
#' no outside option baked into the fit), the first entry of `truth$delta` is
#' dropped before the comparison so lengths match.
#'
#' @param object A `choicer_fit` object (MNL, MXL, or NL) or a `choicer_mc`
#'   result.
#' @param truth Either a `choicer_sim` object (whose `$true_params` will be
#'   used) or a named list of true parameter values.
#' @param level Confidence level for the Wald CI and coverage indicator.
#'   Default `0.95`.
#' @param ... Unused.
#' @return See class-specific methods.
#' @examples
#' \donttest{
#' sim <- simulate_mnl_data(N = 2000, J = 4, seed = 123)
#' fit <- run_mnlogit(
#'   data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = c("x1", "x2"),
#'   outside_opt_label = 0L, include_outside_option = FALSE, use_asc = TRUE
#' )
#' recovery_table(fit, sim)
#' }
#' @export
recovery_table <- function(object, truth = NULL, level = 0.95, ...) {
  UseMethod("recovery_table")
}

#' @describeIn recovery_table Returns a `choicer_recovery` object (a
#'   `data.table`) with columns `parameter`, `group`, `true`, `estimate`,
#'   `se`, `bias`, `rel_bias_pct`, `z_vs_true`, `lower_ci`, `upper_ci`,
#'   `covers`.
#' @export
recovery_table.choicer_fit <- function(object, truth = NULL, level = 0.95, ...) {
  if (inherits(truth, "choicer_sim")) truth <- truth$true_params
  if (!is.list(truth)) stop("`truth` must be a named list or a `choicer_sim`.")
  z <- stats::qnorm(1 - (1 - level) / 2)

  coefs <- stats::coef(object)
  se    <- object$se %||% rep(NA_real_, length(coefs))
  if (length(se) != length(coefs)) {
    stop(sprintf(
      "recovery_table: length(se)=%d != length(coefficients)=%d for this fit.",
      length(se), length(coefs)
    ))
  }
  pm    <- object$param_map

  blocks <- .align_truth_to_coefs(truth, coefs, pm, fit = object)
  rows <- lapply(
    Filter(Negate(is.null), blocks),
    function(b) .build_recovery_rows(b$group, b$idx, b$truth_vec, coefs, se, z)
  )
  out <- data.table::rbindlist(rows)
  class(out) <- c("choicer_recovery", class(out))
  attr(out, "level") <- level
  attr(out, "model") <- setdiff(class(object), c("choicer_fit", "list"))[1]
  out
}

#' @describeIn recovery_table For a `choicer_mc` object, delegates to
#'   `summary(object, level)` and returns a `choicer_mc_summary`. Inspect
#'   `object$replications` directly for per-rep detail.
#' @export
recovery_table.choicer_mc <- function(object, truth = NULL, level = 0.95, ...) {
  summary(object, level = level)
}

#' @export
print.choicer_recovery <- function(x, ...) {
  cat("<choicer_recovery>",
      if (!is.null(attr(x, "model"))) paste0(" model=", attr(x, "model")) else "",
      if (!is.null(attr(x, "level"))) paste0(" level=", attr(x, "level")) else "",
      "\n", sep = "")
  show <- data.table::copy(x)
  class(show) <- class(show)[!class(show) %in% "choicer_recovery"]
  numeric_cols <- setdiff(names(show), c("parameter", "group", "covers"))
  for (col in numeric_cols) show[, (col) := round(get(col), 4)]
  print(show)
  invisible(x)
}

# monte_carlo() ================================================================

#' Monte Carlo parameter recovery
#'
#' Replicates a `(DGP -> fit)` cycle `R` times with independent seeds and
#' collects per-parameter estimates, standard errors, bias, and coverage.
#' Returns a `choicer_mc` object; call `summary()` for aggregated statistics
#' (mean estimate, bias, RMSE, coverage rate, convergence rate).
#'
#' Each iteration calls `sim_fun(seed = seed + r - 1L)`, then `fit_fun(sim)`.
#' Write `sim_fun` as a closure that captures `N`, `J`, and other DGP settings
#' and forwards `seed`. Write `fit_fun` as a closure that takes a
#' `choicer_sim` and returns a fitted `choicer_fit` object, wrapping any
#' data-preparation, draws, or optimizer-control setup.
#'
#' @param sim_fun Function of `seed` returning a `choicer_sim`.
#' @param fit_fun Function of a `choicer_sim` returning a `choicer_fit`.
#' @param R Number of replications.
#' @param seed Base integer seed. Replication `r` uses `seed + r - 1L`.
#' @param parallel Logical; if `TRUE` and `future.apply` is available, run
#'   replications in parallel using the user's active `future::plan()`.
#' @param progress Logical; print a one-line progress update per iteration in
#'   serial mode. Ignored when `parallel = TRUE`.
#' @param ... Unused.
#' @return A `choicer_mc` object: a list with elements `replications` (a long
#'   `data.table` with one row per estimated parameter per replication) and
#'   `meta` (run metadata).
#' @examples
#' \donttest{
#' sim_fun <- function(seed) simulate_mnl_data(N = 1000, J = 4, seed = seed)
#' fit_fun <- function(sim) run_mnlogit(
#'   data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = c("x1", "x2"), outside_opt_label = 0L,
#'   include_outside_option = FALSE, use_asc = TRUE,
#'   control = list(print_level = 0L)
#' )
#' mc <- monte_carlo(sim_fun, fit_fun, R = 5, seed = 1L, progress = FALSE)
#' summary(mc)
#' }
#' @export
monte_carlo <- function(sim_fun, fit_fun, R = 100, seed = 1L,
                        parallel = FALSE, progress = TRUE, ...) {
  seeds <- as.integer(seed) + seq_len(R) - 1L

  na_row <- function(r, elapsed, err_msg) {
    data.table::data.table(
      rep_id = r, seed = seeds[r], parameter = NA_character_, group = NA_character_,
      true = NA_real_, estimate = NA_real_, se = NA_real_,
      bias = NA_real_, rel_bias_pct = NA_real_, z_vs_true = NA_real_,
      lower_ci = NA_real_, upper_ci = NA_real_, covers = NA,
      loglik = NA_real_, converged = FALSE, time_sec = elapsed,
      error = err_msg
    )
  }

  run_rep <- function(r) {
    sim <- sim_fun(seed = seeds[r])
    tic <- Sys.time()
    fit <- tryCatch(fit_fun(sim), error = function(e) e)
    toc <- Sys.time()
    elapsed <- as.numeric(difftime(toc, tic, units = "secs"))

    if (inherits(fit, "error")) {
      return(na_row(r, elapsed, conditionMessage(fit)))
    }

    rt <- tryCatch(recovery_table(fit, sim), error = function(e) e)
    if (inherits(rt, "error")) {
      return(na_row(r, elapsed, conditionMessage(rt)))
    }

    rt[, `:=`(
      rep_id    = r,
      seed      = seeds[r],
      loglik    = as.numeric(stats::logLik(fit)),
      converged = .is_converged(fit),
      time_sec  = elapsed,
      error     = NA_character_
    )]
    data.table::setcolorder(rt, c("rep_id", "seed", "parameter", "group", "true",
                                  "estimate", "se", "bias", "rel_bias_pct",
                                  "z_vs_true", "lower_ci", "upper_ci", "covers",
                                  "loglik", "converged", "time_sec", "error"))
    rt
  }

  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      warning("`future.apply` not installed; falling back to serial execution.")
      parallel <- FALSE
    }
  }

  if (parallel) {
    reps <- future.apply::future_lapply(seq_len(R), run_rep, future.seed = TRUE)
  } else {
    reps <- vector("list", R)
    for (r in seq_len(R)) {
      if (progress) message(sprintf("[monte_carlo] rep %d/%d (seed=%d)", r, R, seeds[r]))
      reps[[r]] <- run_rep(r)
    }
  }

  replications <- data.table::rbindlist(reps, use.names = TRUE, fill = TRUE)

  structure(
    list(
      replications = replications,
      meta = list(
        R = R, seed = seed, seeds = seeds,
        parallel = parallel, timestamp = Sys.time()
      )
    ),
    class = "choicer_mc"
  )
}

#' @export
summary.choicer_mc <- function(object, level = 0.95, ...) {
  reps <- object$replications
  keep <- reps[!is.na(parameter) & converged == TRUE]

  agg <- keep[, .(
    R_success  = .N,
    mean_est   = mean(estimate, na.rm = TRUE),
    median_est = stats::median(estimate, na.rm = TRUE),
    sd_est     = stats::sd(estimate, na.rm = TRUE),
    mean_se    = mean(se, na.rm = TRUE),
    bias       = mean(estimate - true, na.rm = TRUE),
    rmse       = sqrt(mean((estimate - true)^2, na.rm = TRUE)),
    coverage   = mean(covers, na.rm = TRUE)
  ), by = .(parameter, group, true)]
  data.table::setcolorder(agg, c("parameter", "group", "true", "R_success",
                                 "mean_est", "median_est", "sd_est", "mean_se",
                                 "bias", "rmse", "coverage"))

  conv_rate <- mean(reps[, any(converged, na.rm = TRUE), by = rep_id]$V1)

  structure(
    list(
      summary = agg,
      meta = object$meta,
      conv_rate = conv_rate,
      n_reps = object$meta$R,
      level = level
    ),
    class = "choicer_mc_summary"
  )
}

#' @export
print.choicer_mc <- function(x, ...) {
  cat("<choicer_mc> R=", x$meta$R, " seed=", x$meta$seed,
      " parallel=", x$meta$parallel, "\n", sep = "")
  cat("  replications: ", nrow(x$replications), " rows\n", sep = "")
  cat("  columns: ", paste(names(x$replications), collapse = ", "), "\n", sep = "")
  cat("Call summary() for aggregated recovery statistics.\n")
  invisible(x)
}

#' @export
print.choicer_mc_summary <- function(x, ...) {
  cat("<choicer_mc_summary> R=", x$n_reps,
      " convergence_rate=", round(x$conv_rate, 3),
      " coverage_level=", x$level, "\n", sep = "")
  show <- data.table::copy(x$summary)
  numeric_cols <- setdiff(names(show), c("parameter", "group", "R_success"))
  for (col in numeric_cols) show[, (col) := round(get(col), 4)]
  print(show)
  invisible(x)
}
