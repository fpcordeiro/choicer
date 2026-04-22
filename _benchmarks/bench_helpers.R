# Cross-package benchmark helpers. Sourced by _benchmarks/*.R scripts.
#
# These helpers fit a `choicer_sim` object (from choicer::simulate_*_data())
# with one or more reference packages (mlogit, logitr, mixl, gmnl, apollo) and
# return a long data.table suitable for comparing coefficients, SEs, log-
# likelihoods, and wall time across implementations.
#
# Not part of the installed package: reference packages live in library()
# calls here instead of in DESCRIPTION Suggests. Install them manually before
# running a benchmark script (e.g. install.packages(c("mlogit", "logitr"))).

suppressPackageStartupMessages({
  library(data.table)
  library(choicer)
})

# Output schema ================================================================
#
# A `choicer_benchmark` is a long data.table with one row per (package, parameter):
#
#   package     character  name of the reference package ("choicer", "mlogit", ...)
#   parameter   character  parameter name in that package's output
#   group       character  one of beta, mu, sigma, lambda, asc, other
#   estimate    numeric    point estimate
#   se          numeric    standard error
#   loglik      numeric    final log-likelihood (repeated per parameter)
#   time_sec    numeric    wall time for fit() (repeated per parameter)
#   converged   logical    convergence flag (repeated per parameter)
#   n_params    integer    total number of estimated parameters (repeated)

.bench_row <- function(package, parameter, group, estimate, se, loglik, time_sec,
                       converged, n_params) {
  data.table(
    package   = package,
    parameter = parameter,
    group     = group,
    estimate  = as.numeric(estimate),
    se        = as.numeric(se),
    loglik    = as.numeric(loglik),
    time_sec  = as.numeric(time_sec),
    converged = as.logical(converged),
    n_params  = as.integer(n_params)
  )
}

.require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "Package '%s' is required for this dispatcher. Install it with install.packages('%s').",
      pkg, pkg
    ))
  }
}

# Top-level dispatcher =========================================================

#' Benchmark a choicer_sim across multiple packages
#'
#' @param sim A `choicer_sim` object from simulate_mnl_data / simulate_mxl_data
#'   / simulate_nl_data.
#' @param packages Character vector of reference packages. Recognized names:
#'   "choicer", "mlogit", "logitr", "mixl", "gmnl", "apollo". Not all packages
#'   support every model; invalid combinations raise a clear error.
#' @param ... Extra arguments forwarded to the per-package dispatcher.
#' @return A `choicer_benchmark` (a data.table).
benchmark_fit <- function(sim, packages = c("choicer"), ...) {
  stopifnot(inherits(sim, "choicer_sim"))
  model <- sim$model
  rows <- list()
  for (pkg in packages) {
    disp <- .dispatch(model, pkg)
    message(sprintf("[benchmark_fit] %s x %s ...", model, pkg))
    rows[[pkg]] <- tryCatch(disp(sim, ...), error = function(e) {
      warning(sprintf("%s/%s failed: %s", model, pkg, conditionMessage(e)))
      NULL
    })
  }
  out <- rbindlist(Filter(Negate(is.null), rows), use.names = TRUE, fill = TRUE)
  class(out) <- c("choicer_benchmark", class(out))
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Duplicate of .is_converged() from R/recovery.R. Kept here because this file
# is source()'d (not loaded as part of the package), so namespace-internal
# helpers aren't reachable.
.is_converged <- function(fit) {
  code <- fit$convergence
  if (is.null(code) || length(code) != 1L || !is.finite(code)) return(FALSE)
  nm <- fit$optimizer$name %||% "nloptr"
  switch(nm,
    nloptr = isTRUE(as.integer(code) %in% c(1L, 2L, 3L, 4L)),
    optim  = isTRUE(as.integer(code) == 0L),
    isTRUE(as.integer(code) %in% c(0L, 1L, 2L, 3L, 4L))
  )
}

.dispatch <- function(model, pkg) {
  key <- paste0("bench_", model, "_", pkg)
  # globalenv() covers scripts that define their own bench_* dispatchers
  # before calling benchmark_fit(); the second lookup finds the built-in
  # dispatchers defined later in this file.
  fn <- get0(key, envir = globalenv(), mode = "function") %||%
        get0(key, envir = environment(.dispatch), mode = "function")
  if (is.null(fn)) {
    stop(sprintf("No dispatcher for model='%s' package='%s' (expected `%s`).",
                 model, pkg, key))
  }
  fn
}

# choicer dispatchers ==========================================================

bench_mnl_choicer <- function(sim, ...) {
  tic <- Sys.time()
  fit <- run_mnlogit(
    data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = paste0("x", seq_len(sim$settings$K_x)),
    outside_opt_label = 0L, include_outside_option = FALSE, use_asc = TRUE,
    control = list(print_level = 0L)
  )
  elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
  pm <- fit$param_map
  cf <- coef(fit); se <- fit$se
  groups <- rep("other", length(cf))
  groups[pm$beta] <- "beta"
  if (!is.null(pm$asc)) groups[pm$asc] <- "asc"
  .bench_row("choicer", names(cf), groups, cf, se,
             as.numeric(logLik(fit)), elapsed,
             .is_converged(fit), length(cf))
}

bench_mxl_choicer <- function(sim, S = NULL, rc_dist = NULL, rc_mean = TRUE, ...) {
  input <- prepare_mxl_data(
    data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = paste0("x", seq_len(sim$settings$K_x)),
    random_var_cols = paste0("w", seq_len(sim$settings$K_w)),
    outside_opt_label = 0L, include_outside_option = FALSE,
    rc_correlation = sim$true_params$rc_correlation
  )
  if (is.null(S)) S <- 50 * ceiling((1.5 * sqrt(input$N)) / 50)
  eta <- get_halton_normals(S, input$N, sim$settings$K_w)
  tic <- Sys.time()
  fit <- run_mxlogit(
    input_data = input, eta_draws = eta,
    rc_dist = rc_dist %||% sim$true_params$rc_dist,
    rc_mean = rc_mean, use_asc = TRUE,
    control = list(print_level = 0L)
  )
  elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
  pm <- fit$param_map
  cf <- coef(fit); se <- fit$se
  groups <- rep("other", length(cf))
  groups[pm$beta] <- "beta"
  if (!is.null(pm$mu))    groups[pm$mu]    <- "mu"
  if (!is.null(pm$sigma)) groups[pm$sigma] <- "sigma"
  if (!is.null(pm$asc))   groups[pm$asc]   <- "asc"
  .bench_row("choicer", names(cf), groups, cf, se,
             as.numeric(logLik(fit)), elapsed,
             .is_converged(fit), length(cf))
}

bench_nl_choicer <- function(sim, ...) {
  # `sim$data` already contains a `nest` column (integer; 0L = outside).
  tic <- Sys.time()
  fit <- run_nestlogit(
    data = sim$data, id_col = "id", alt_col = "j", choice_col = "choice",
    covariate_cols = c("X", "W"), nest_col = "nest",
    use_asc = TRUE, include_outside_option = TRUE, outside_opt_label = 0L,
    control = list(print_level = 0L)
  )
  elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
  pm <- fit$param_map
  cf <- coef(fit); se <- fit$se
  groups <- rep("other", length(cf))
  groups[pm$beta] <- "beta"
  if (!is.null(pm$lambda)) groups[pm$lambda] <- "lambda"
  if (!is.null(pm$asc))    groups[pm$asc]    <- "asc"
  .bench_row("choicer", names(cf), groups, cf, se,
             as.numeric(logLik(fit)), elapsed,
             .is_converged(fit), length(cf))
}

# mlogit dispatchers ===========================================================

bench_mnl_mlogit <- function(sim, ...) {
  .require("mlogit"); .require("dfidx")
  dt <- copy(sim$data)
  dt[, alt_factor := factor(alt)]
  x_cols <- paste0("x", seq_len(sim$settings$K_x))
  form <- stats::as.formula(paste("choice ~", paste(x_cols, collapse = " + ")))

  tic <- Sys.time()
  dt_idx <- dfidx::dfidx(as.data.frame(dt), idx = c("id", "alt_factor"))
  m <- mlogit::mlogit(form, data = dt_idx, print.level = 0)
  elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))

  s <- summary(m)
  est <- stats::coef(m)
  se  <- s$CoefTable[, "Std. Error"]
  # mlogit names alt-specific intercepts "(intercept):<level>"
  groups <- ifelse(names(est) %in% x_cols, "beta", "asc")
  .bench_row("mlogit", names(est), groups, est, se,
             as.numeric(stats::logLik(m)), elapsed,
             TRUE, length(est))
}

bench_nl_mlogit <- function(sim, ...) {
  .require("mlogit"); .require("dfidx")
  # mlogit cannot fit a singleton outside-option nest; restrict to ids whose
  # chosen alternative is inside the nesting structure, then drop j == 0.
  dt <- copy(sim$data)
  ids_inside <- unique(dt[choice == 1L & j > 0L, id])
  dt <- dt[id %in% ids_inside & j > 0L]
  dt[, alt_factor := factor(j)]
  nest_spec <- list()
  for (g in seq_along(sim$settings$nest_structure)) {
    nest_spec[[paste0("nest_", g)]] <- as.character(sim$settings$nest_structure[[g]])
  }

  tic <- Sys.time()
  dt_idx <- dfidx::dfidx(as.data.frame(dt), idx = c("id", "alt_factor"))
  m <- mlogit::mlogit(choice ~ X + W, data = dt_idx, nests = nest_spec, print.level = 0)
  elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))

  s <- summary(m)
  est <- stats::coef(m)
  se  <- s$CoefTable[, "Std. Error"]
  groups <- vapply(names(est), function(nm) {
    if (nm %in% c("X", "W")) "beta"
    else if (grepl("^iv", nm) || grepl("nest", nm)) "lambda"
    else "asc"
  }, character(1))
  .bench_row("mlogit", names(est), groups, est, se,
             as.numeric(stats::logLik(m)), elapsed,
             TRUE, length(est))
}

# logitr dispatchers ===========================================================

bench_mnl_logitr <- function(sim, ...) {
  .require("logitr")
  dt <- copy(sim$data)
  dt[, alt_factor := factor(alt)]
  x_cols <- paste0("x", seq_len(sim$settings$K_x))

  tic <- Sys.time()
  fit <- logitr::logitr(
    data = dt, outcome = "choice", obsID = "id",
    pars = c(x_cols, "alt_factor"),
    options = list(
      algorithm = "NLOPT_LD_LBFGS",
      xtol_rel = 1e-8, maxeval = 1000L, print_level = 0L
    )
  )
  elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))

  est <- fit$coefficients
  V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  se <- if (!is.null(V)) sqrt(pmax(diag(V), 0)) else rep(NA_real_, length(est))
  if (length(se) != length(est)) se <- rep(NA_real_, length(est))
  groups <- ifelse(names(est) %in% x_cols, "beta", "asc")
  .bench_row("logitr", names(est), groups, est, se,
             as.numeric(fit$logLik), elapsed,
             isTRUE(fit$status >= 0), length(est))
}

bench_mxl_logitr <- function(sim, S = NULL, rc_dist = NULL, ...) {
  .require("logitr"); .require("randtoolbox")
  dt <- copy(sim$data)
  dt[, alt_factor := factor(alt)]
  K_x <- sim$settings$K_x
  K_w <- sim$settings$K_w
  x_cols <- paste0("x", seq_len(K_x))
  w_cols <- paste0("w", seq_len(K_w))
  if (is.null(S)) S <- 50 * ceiling((1.5 * sqrt(sim$settings$N)) / 50)
  if (is.null(rc_dist)) rc_dist <- sim$true_params$rc_dist

  # logitr randPars: "n" = normal, "ln" = log-normal
  rand_spec <- stats::setNames(
    ifelse(rc_dist == 1L, "ln", "n"),
    w_cols
  )
  # logitr's standardDraws expects one column per random parameter block.
  # Dimension grows with correlation; let logitr generate draws internally
  # via numDraws to avoid version-specific layout mismatches.
  tic <- Sys.time()
  fit <- logitr::logitr(
    data = dt, outcome = "choice", obsID = "id",
    pars = c(x_cols, w_cols, "alt_factor"),
    randPars = rand_spec,
    correlation = sim$true_params$rc_correlation,
    numDraws = S,
    options = list(
      algorithm = "NLOPT_LD_LBFGS",
      xtol_rel = 1e-8, maxeval = 2000L, print_level = 0L
    )
  )
  elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))

  est <- fit$coefficients
  se  <- tryCatch(sqrt(diag(stats::vcov(fit))), error = function(e) rep(NA_real_, length(est)))
  groups <- vapply(names(est), function(nm) {
    if (nm %in% x_cols) "beta"
    else if (grepl("^sd_", nm) || grepl("^sigma", nm)) "sigma"
    else if (nm %in% w_cols) "mu"
    else "asc"
  }, character(1))
  .bench_row("logitr", names(est), groups, est, se,
             as.numeric(fit$logLik), elapsed,
             isTRUE(fit$status >= 0), length(est))
}

# Stubs for packages without working dispatchers yet ===========================
# These raise a clear error so benchmark_fit() surfaces the gap rather than
# silently skipping. Implement them by adapting the reference package's own
# documentation/vignette and mirroring the structure of bench_mxl_logitr().

bench_mxl_mixl   <- function(sim, ...) stop("Not yet implemented: MXL dispatcher for 'mixl'. See _benchmarks/bench_helpers.R.")
bench_mxl_gmnl   <- function(sim, ...) stop("Not yet implemented: MXL dispatcher for 'gmnl'. See _benchmarks/bench_helpers.R.")
bench_mxl_apollo <- function(sim, ...) stop("Not yet implemented: MXL dispatcher for 'apollo'. See _benchmarks/bench_helpers.R.")

# Print method =================================================================

print.choicer_benchmark <- function(x, ...) {
  cat("<choicer_benchmark> packages=",
      paste(unique(x$package), collapse = ", "), "\n", sep = "")
  show <- copy(x)
  class(show) <- class(show)[!class(show) %in% "choicer_benchmark"]
  for (col in c("estimate", "se", "loglik", "time_sec")) {
    show[, (col) := round(get(col), 4)]
  }
  print(show)
  invisible(x)
}
