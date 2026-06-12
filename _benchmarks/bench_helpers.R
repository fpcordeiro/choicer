# Cross-package benchmark helpers. Sourced by _benchmarks/*.R scripts.
#
# These helpers fit a `choicer_sim` object (from choicer::simulate_*_data())
# with one or more reference packages (mlogit, logitr, mixl, gmnl, apollo for
# the MLE models; bayesm, MNP for the Bayesian multinomial probit) and return
# a long data.table suitable for comparing coefficients, SEs, log-
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
#
# Bayesian (MNP) benchmarks reinterpret the columns: estimate / se are the
# posterior mean / SD of the identified draws, loglik and converged are NA,
# and an extra `ess` column carries the effective sample size of each
# parameter's kept draws.

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
#'   "choicer", "mlogit", "logitr", "mixl", "gmnl", "apollo", "bayesm",
#'   "MNP". Not all packages support every model; invalid combinations raise
#'   a clear error.
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

# Bayesian MNP dispatchers =====================================================
#
# All three samplers target the same model (latent utility differences vs a
# base alternative, multivariate normal errors). Comparisons are made on the
# identified scale: each kept draw is normalized by sigma_11, so estimate =
# E[beta / sqrt(sigma_11)] and Sigma rows are E[Sigma / sigma_11] -- directly
# comparable to `sim$true_params` from simulate_mnp_data(). Parameter names
# are canonicalized to choicer's (x*, ASC_*, Sigma_ij with row-wise lower
# triangle), so the long table lines up across packages.
#
# Dispatchers accept a shared `mcmc = list(R, burn, thin)` so every package
# runs the same chain length on the same dataset.

# Same initial-positive-sequence ESS estimator as
# inst/simulations/mnp_simulation.R (duplicated: this file is source()'d,
# not part of the package).
.ess <- function(x) {
  n <- length(x)
  rho <- as.numeric(stats::acf(x, lag.max = min(200L, n - 1L), plot = FALSE)$acf)[-1]
  cut <- which(rho < 0)[1]
  if (!is.na(cut)) rho <- rho[seq_len(cut - 1L)]
  n / (1 + 2 * sum(rho))
}

.mnp_mcmc_defaults <- function(mcmc) {
  list(
    R    = as.integer(mcmc$R %||% 20000L),
    burn = as.integer(mcmc$burn %||% 5000L),
    thin = as.integer(mcmc$thin %||% 5L)
  )
}

# choicer's Sigma draw naming: row-wise lower triangle (11, 21, 22, 31, ...)
.sigma_names <- function(p) {
  nm <- character(p * (p + 1L) / 2L)
  k <- 1L
  for (i in seq_len(p)) {
    for (j in seq_len(i)) {
      nm[k] <- sprintf("Sigma_%d%d", i, j)
      k <- k + 1L
    }
  }
  nm
}

.mnp_bench_rows <- function(package, beta_id, sigma_id, groups_beta, time_sec) {
  draws <- cbind(beta_id, sigma_id)
  out <- .bench_row(
    package, colnames(draws),
    c(groups_beta, rep("sigma", ncol(sigma_id))),
    colMeans(draws), apply(draws, 2, stats::sd),
    NA_real_, time_sec, NA, ncol(draws)
  )
  out[, ess := apply(draws, 2, .ess)]
  # Sigma_11 == 1 by construction on the identified scale; ESS is undefined
  out[parameter == "Sigma_11", ess := NA_real_]
  out[]
}

bench_mnp_choicer <- function(sim, mcmc = list(), ...) {
  m <- .mnp_mcmc_defaults(mcmc)
  tic <- Sys.time()
  fit <- run_mnprobit(
    data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = paste0("x", seq_len(sim$settings$K_x)),
    use_asc = TRUE,
    mcmc = list(R = m$R, burn = m$burn, thin = m$thin)
  )
  elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
  pm <- fit$param_map
  groups <- rep("other", fit$n_params)
  groups[pm$beta] <- "beta"
  if (!is.null(pm$asc)) groups[pm$asc] <- "asc"
  .mnp_bench_rows("choicer", fit$draws$beta, fit$draws$sigma, groups, elapsed)
}

bench_mnp_bayesm <- function(sim, mcmc = list(), ...) {
  .require("bayesm")
  m <- .mnp_mcmc_defaults(mcmc)
  if (m$burn %% m$thin != 0L) {
    stop("bayesm dispatcher requires mcmc$burn divisible by mcmc$thin.")
  }
  # Closest apples-to-apples reference: rmnpGibbs implements the same
  # McCulloch-Rossi (1994) Gibbs sampler, and runs here on the *identical*
  # differenced design matrix and priors as choicer (via prepare_mnp_data()).
  # bayesm codes the base alternative as y = J; choicer codes it as 0.
  input <- prepare_mnp_data(
    sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = paste0("x", seq_len(sim$settings$K_x)), use_asc = TRUE
  )
  y <- ifelse(input$y == 0L, input$J, input$y)
  nu <- input$p + 3

  tic <- Sys.time()
  out <- bayesm::rmnpGibbs(
    Data = list(y = y, X = input$X, p = input$J),
    Prior = list(betabar = rep(0, input$K), A = 0.01 * diag(input$K),
                 nu = nu, V = nu * diag(input$p)),
    Mcmc = list(R = m$R, keep = m$thin, nprint = 0)
  )
  elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))

  # rmnpGibbs has no burn-in argument: drop the first burn/thin kept rows,
  # then normalize per draw by sigma_11.
  drop_rows <- seq_len(m$burn %/% m$thin)
  betadraw <- out$betadraw[-drop_rows, , drop = FALSE]
  sigmadraw <- out$sigmadraw[-drop_rows, , drop = FALSE]
  s11 <- sigmadraw[, 1L]
  beta_id <- betadraw / sqrt(s11)
  colnames(beta_id) <- colnames(input$X)

  # sigmadraw rows are vec(Sigma) (column-major p x p); reorder to choicer's
  # row-wise lower triangle
  p <- input$p
  lower_idx <- unlist(lapply(seq_len(p), function(i) (seq_len(i) - 1L) * p + i))
  sigma_id <- (sigmadraw / s11)[, lower_idx, drop = FALSE]
  colnames(sigma_id) <- .sigma_names(p)

  groups <- rep("other", ncol(beta_id))
  groups[input$param_map$beta] <- "beta"
  if (!is.null(input$param_map$asc)) groups[input$param_map$asc] <- "asc"
  .mnp_bench_rows("bayesm", beta_id, sigma_id, groups, elapsed)
}

bench_mnp_MNP <- function(sim, mcmc = list(), ...) {
  .require("MNP")
  m <- .mnp_mcmc_defaults(mcmc)
  x_cols <- paste0("x", seq_len(sim$settings$K_x))
  base_alt <- sim$settings$base_alt

  # MNP::mnp() wants wide data (one row per choice situation) with a factor
  # response, plus a list of per-alternative covariate matrices (choiceX).
  dt <- copy(sim$data)
  wide <- dcast(dt, id ~ alt, value.var = x_cols)
  wide <- merge(wide, dt[choice == 1L, .(id, choice = factor(alt))], by = "id")
  setorder(wide, id)
  alts <- sort(unique(dt$alt))
  choiceX <- lapply(alts, function(a) {
    as.matrix(wide[, paste0(x_cols, "_", a), with = FALSE])
  })
  names(choiceX) <- as.character(alts)

  # Imai-van Dyk (2005) marginal data augmentation sampler. Priors are
  # matched where the parameterizations align: p.var = 100 is choicer's
  # prior precision A = 0.01 * I, p.df its nu = p + 3. MNP identifies the
  # scale in-sampler (trace restriction by default), so its implied prior on
  # Sigma is not identical to the inverse-Wishart used by choicer/bayesm;
  # draws are renormalized to the sigma_11 = 1 scale below either way.
  # MNP's `thin` counts skipped draws, hence thin - 1 to keep every thin-th.
  # mnp() forwards `choiceX` unevaluated and evaluates it against `data` with
  # fallback to the global environment, so a local variable is not visible
  # here; do.call() inlines the evaluated list into the call.
  p <- length(alts) - 1L
  tic <- Sys.time()
  fit <- do.call(MNP::mnp, list(
    choice ~ 1, data = as.data.frame(wide),
    choiceX = choiceX, cXnames = x_cols,
    base = as.character(base_alt),
    n.draws = m$R, burnin = m$burn, thin = m$thin - 1L,
    p.var = 100, p.df = p + 3,
    verbose = FALSE
  ))
  elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))

  # fit$param columns: coefficients, then the upper triangle of Sigma named
  # by alternative labels ("2:2", "2:3", ...). Note "(Intercept):2" also
  # contains ":", so covariance columns are matched by exact name.
  draws <- fit$param
  non_base <- setdiff(alts, base_alt)
  sig_cols <- unlist(lapply(seq_len(p), function(i) {
    vapply(seq_len(i), function(j) paste0(non_base[j], ":", non_base[i]), character(1))
  }))
  stopifnot(all(sig_cols %in% colnames(draws)))

  s11 <- draws[, sig_cols[1L]]
  sigma_id <- draws[, sig_cols, drop = FALSE] / s11
  colnames(sigma_id) <- .sigma_names(p)

  beta_id <- draws[, setdiff(colnames(draws), sig_cols), drop = FALSE] / sqrt(s11)
  colnames(beta_id) <- sub("^\\(Intercept\\):", "ASC_", colnames(beta_id))
  # match choicer's column order: covariates first, then ASCs
  beta_id <- beta_id[, c(x_cols, setdiff(colnames(beta_id), x_cols)), drop = FALSE]

  groups <- ifelse(colnames(beta_id) %in% x_cols, "beta", "asc")
  .mnp_bench_rows("MNP", beta_id, sigma_id, groups, elapsed)
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
  if ("ess" %in% names(show)) show[, ess := round(ess)]
  print(show)
  invisible(x)
}
