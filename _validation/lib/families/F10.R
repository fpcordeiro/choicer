# ============================================================================
# F10: the share-matching constants, run_mxlogit(concentrate_out_asc =
# "blp"), against maximum likelihood (plan §5; claim C1): the port of
# _validation/blp_concentration_mc.R, whose four designs and DGP it keeps
# statement for statement, so that the mixed logit article's §5.9 table and
# NEWS stay reproducible from the runner (plan §3). Findings only.
# Designs (cells), each 2,000 choice situations among J = 10 alternatives,
# every situation a random subset of at least half of its market's:
#   cross    a cross-section with an outside option
#   panel    500 decision makers x 4 situations, an outside option
#   base     a cross-section without an outside option (alternative 1 the
#            reference)
#   markets  two markets of five alternatives (each decision maker in one),
#            an outside option, targets within each market
# DGP (base R and data.table: neither choicer's simulator nor its kernels
# draw the data or the choices; choicer enters only through the fit's
# inputs below, prepare_mxl_data() and the estimation draws of
# get_halton_normals()):
# V = x1 + (-1 + e1_n) w1 + (0.5 + 0.7 e2_n) w2 + delta_j + Gumbel, x1 ~
# N(0, 1), w ~ N(0, 1.5^2), e ~ N(0, 1) per decision maker, delta evenly
# spaced from -1 to 1; the choice by inverse CDF of the exact logit
# probabilities given e. Fitted with X = (x1, w1, w2), W = (w1, w2), a
# diagonal Cholesky (rc_mean = FALSE), as the old script did: through
# prepare_mxl_data() and the store-mode cube get_halton_normals(S = 100, U,
# 2) passed as eta_draws (unscrambled Halton, the same for every rep), not
# the suite's generate-mode draws: these are what the article's numbers
# were computed with.
# Seeds. Reps 1..1000 run on the old script's seeds, 1000 * (the design's
# position in cross, panel, base, markets) + rep, so that the port gives
# the old script's data and fits rep by rep; later reps on the runner's own
# seed. The population's within-market shares (the exact targets) are the
# old script's pop_shares(): the logit probabilities averaged over 10^6
# simulated situations (seed 424242), frozen below as literals
# (f10_pop_shares() recomputes them, about 25 s per design;
# checks/wp7_sim.R compares).
# Methods (every fit through fit_two_starts(), lib/fit.R: from choicer's
# default start, the old script's, and from theta0, the better certified
# one kept, its record in diag$start; multi-start on reps 10, 20, ...;
# routes hessian and robust):
#   ml          polish = "newton" (max_iter 50, tol_decrement 1e-14);
#   blp_sample  "blp" with the sample's shares;
#   blp_exact   "blp" with the population's shares, as exact;
#   blp_n1, blp_n10  "blp" with shares from n = N and n = 10 N multinomial
#               choices per market (N = 2,000 / markets), target_vcov their
#               multinomial covariance; drawn on the old seed + 77 k.
# Findings (f10_findings(records), the §5.9 table with MCSEs): per design,
# estimator and group (theta_c, the constants) the SD against ML's on the
# reps both fitted, mean SE over SD, 95% coverage of the Hessian and
# sandwich intervals, the largest |bias| / SD, and the work per fit; then
# the article's table, ranges over the designs. Sourced by load_family()
# into a child of the library's environment.
# ============================================================================

F10_CELLS <- c("cross", "panel", "base", "markets")
F10_DELTA <- seq(-1, 1, length.out = 10L)

# The population's within-market shares of the inside alternatives (the
# exact targets), by design: f10_pop_shares() of 2026-11-02, the old
# script's pop_shares() to the last bit.
# Hexadecimal constants (sprintf("%a")): R's parser reads some 17-digit
# decimals an ulp off (one share each of cross and markets), and the old
# script's fits use these exact doubles. In decimals, cross 0.0503 to
# 0.1556, panel 0.0504 to 0.1554, base 0.0531 to 0.1620, markets 0.1291 to
# 0.2044 (market 1) and 0.1457 to 0.2252 (market 2), rising with delta.
F10_SHARES <- list(
  cross = c(0x1.9be903526ad95p-5, 0x1.d728f57be7ea1p-5, 0x1.0cea20a3615dap-4,
            0x1.319d82e43aad6p-4, 0x1.5be34116af00ap-4, 0x1.891ce1747b9f5p-4,
            0x1.bdc0a24f9aad5p-4, 0x1.f7acb6e55d696p-4, 0x1.1b4d74d1909f7p-3,
            0x1.3e95a62399167p-3),
  panel = c(0x1.9ce8482a0c4b6p-5, 0x1.d468fd1b10dbep-5, 0x1.0cf37ba126c9bp-4,
            0x1.31c832e0be2bp-4, 0x1.5cd3d6e7cd7a8p-4, 0x1.8974683a5cf66p-4,
            0x1.be53058ca6568p-4, 0x1.f6f4efdd44fffp-4, 0x1.1b4a4b54197f5p-3,
            0x1.3e46dbfeff794p-3),
  base = c(0x1.b34ea8aaf6857p-5, 0x1.f1267282ede19p-5, 0x1.1b5b8468cd131p-4,
           0x1.416d21e403119p-4, 0x1.6d528e069f13ap-4, 0x1.9c4832fc9c188p-4,
           0x1.d2922b517515bp-4, 0x1.07214d81c710ep-3, 0x1.278d05b4a63e6p-3,
           0x1.4bb99cad59852p-3),
  markets = c(0x1.086dfb1361bbdp-3, 0x1.28dcfbac5b8a2p-3, 0x1.4dd405bbb8eep-3,
              0x1.752bf3b4c9839p-3, 0x1.a2ab994919226p-3, 0x1.2a4ee40909964p-3,
              0x1.4cb8de58f629dp-3, 0x1.75a0dee901235p-3, 0x1.9ee7edd78e65bp-3,
              0x1.cd3beb51a2b29p-3))

# The old script's sim_raw(): N decision makers x T situations; each
# decision maker's market; each situation's choice set; covariates, tastes
# e and utilities V (without the Gumbel term). Its random streams are the
# old script's, call for call.
f10_raw <- function(seed, N, J, T, markets, beta = c(1, -1, 0.5), sig = c(1, 0.7),
                    delta = seq(-1, 1, length.out = J), wsd = 1.5) {
  vs_set_seed(seed)
  n_sit <- N * T
  pid_sit <- rep(seq_len(N), each = T)
  mk_alt <- rep(seq_len(markets), each = J %/% markets)
  # each decision maker's market (one market draws nothing)
  mk_sit <- if (markets > 1L) {
    rep(sample.int(markets, N, replace = TRUE), each = T)
  } else {
    rep(1L, n_sit)
  }
  alts <- lapply(seq_len(n_sit), function(t) {
    pool <- which(mk_alt == mk_sit[t])
    k <- ceiling(length(pool) / 2):length(pool)
    k <- k[sample.int(length(k), 1L)]
    sort(pool[sample.int(length(pool), k)])
  })
  id <- rep(seq_len(n_sit), lengths(alts))
  dt <- data.table::data.table(id = id, alt = unlist(alts), pid = pid_sit[id], mk = mk_sit[id])
  dt[, x1 := stats::rnorm(.N)]
  dt[, w1 := stats::rnorm(.N, sd = wsd)]
  dt[, w2 := stats::rnorm(.N, sd = wsd)]
  e <- matrix(stats::rnorm(2L * N), 2L)
  dt[, V := beta[1] * x1 + (beta[2] + sig[1] * e[1, pid]) * w1 +
            (beta[3] + sig[2] * e[2, pid]) * w2 + delta[alt]]
  dt[]
}

# The old script's seed of a rep (reps 1..1000; later reps keep the
# runner's seed).
f10_old_seed <- function(cell, seed) {
  rep <- vs_seed_rep("F10", cell$cell, seed)
  if (rep <= 1000L) 1000L * match(cell$cell, F10_CELLS) + rep else seed
}

# A rep's data: the long table with the choices (the old sim_data()), its
# prepared form and the estimation draws (meta$input, meta$eta), theta0 in
# choicer's names, the exact targets and each alternative's market.
f10_simulate <- function(cell, seed) {
  s0 <- f10_old_seed(cell, seed)
  dt <- f10_raw(s0, cell$N, cell$J, cell$T, cell$markets)
  dt[, p := exp(V) / (sum(exp(V)) + cell$outside), by = id]
  r <- stats::runif(max(dt$id))
  dt[, choice := { cp <- cumsum(p); as.integer(r[id[1L]] <= cp & r[id[1L]] > cp - p) },
     by = id]
  input <- suppressMessages(choicer::prepare_mxl_data(
    dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "w1", "w2"), random_var_cols = c("w1", "w2"),
    include_outside_option = cell$outside, outside_opt_label = if (cell$outside) 0L,
    person_col = if (cell$T > 1L) "pid"))
  U <- if (is.null(input$Ti)) length(input$M) else length(input$Ti)
  delta <- if (cell$outside) F10_DELTA else F10_DELTA[-1L] - F10_DELTA[1L]
  names(delta) <- paste0("ASC_", if (cell$outside) 1:10 else 2:10)
  truth <- c(x1 = 1, w1 = -1, w2 = 0.5, L_11 = 0, L_22 = log(0.7), delta)
  list(data = dt, truth = truth,
       meta = list(input = input, eta = choicer::get_halton_normals(cell$S, U, 2L),
                   old_seed = s0, shares = F10_SHARES[[cell$cell]],
                   market = rep(seq_len(cell$markets), each = cell$J %/% cell$markets)))
}

# The old script's pop_shares(): the logit probabilities averaged over
# n_pop simulated situations, within each market, in alternative order.
f10_pop_shares <- function(cell, n_pop = 1e6) {
  dt <- f10_raw(424242L, as.integer(n_pop / cell$T), cell$J, cell$T, cell$markets)
  dt[, p := exp(V) / (sum(exp(V)) + cell$outside), by = id]
  n_m <- dt[, .(n = data.table::uniqueN(id)), by = mk]
  s <- dt[, .(s = sum(p)), by = .(alt, mk)][n_m, on = "mk"][, s := s / n][order(alt)]
  list(within = s$s, market = s$mk)
}

# --- Targets and fits -------------------------------------------------------

# The targets of a "blp" method: NULL (the sample's), the exact shares, or
# shares from n = k N / markets choices per market (k = 1, 10) with their
# multinomial covariance, a draw made positive by redrawing (the old
# script's, on its seed + 77 k).
f10_targets <- function(sim, cell, kind) {
  m <- sim$meta
  if (kind == "sample") return(list(shares = NULL, vcov = NULL))
  if (kind == "exact") return(list(shares = m$shares, vcov = NULL))
  k <- c(n1 = 1L, n10 = 10L)[[kind]]
  J <- length(m$shares)
  n <- k * length(m$input$M) / cell$markets
  vs_set_seed(m$old_seed + 77L * k)
  s_hat <- numeric(J)
  Vs <- matrix(0, J, J)
  for (mk in seq_len(cell$markets)) {
    im <- which(m$market == mk)
    pin <- m$shares[im]
    probs <- if (cell$outside) c(1 - sum(pin), pin) else pin
    repeat {
      sh <- as.numeric(stats::rmultinom(1L, n, probs)) / n
      if (all(sh > 0)) break
    }
    shin <- if (cell$outside) sh[-1L] else sh
    s_hat[im] <- shin
    Vs[im, im] <- (diag(shin, length(im)) - tcrossprod(shin)) / n
  }
  list(shares = s_hat, vcov = Vs)
}

# One estimator's fit on the rep's prepared data and draws, from the
# default start and from theta0 (fit_two_starts()), and on reps 10, 20, ...
# the multi-start check of the global maximum (attribute "multistart").
f10_fit <- function(sim, cell, seed, what, ...) {
  m <- sim$meta
  run <- function(theta_init) {
    choicer::run_mxlogit(input_data = m$input, eta_draws = m$eta, theta_init = theta_init, ...)
  }
  fit <- fit_two_starts(run, sim$truth, what)
  if (vs_seed_rep("F10", cell$cell, seed) %% 10 == 0) {
    attr(fit, "multistart") <- multistart(fit, refit = function(theta_init) {
      require_converged(run(theta_init), what)
    }, seed = seed)
  }
  fit
}

# A method's result: the est table (hessian and robust SEs), the joint Wald
# test on theta_c, and the work: the start record, the multi-start check
# and, for "blp", the share-matching record.
f10_result <- function(fit, sim) {
  V <- fit_vcovs(fit, c("hessian", "robust"))
  common <- names(coef(fit))[-as.integer(fit$param_map$asc)]
  conc <- fit$optimizer$concentration
  list(est = est_table(fit, sim$truth, V = V, gated = FALSE),
       vcov_c = vcov_common(fit, V = V$hessian),
       joint_wald = joint_wald(coef(fit)[common], sim$truth, V$hessian),
       diag = c(fit_diag(fit), list(
         start = attr(fit, "start"), multistart = attr(fit, "multistart"),
         old_seed = sim$meta$old_seed,
         concentration = if (!is.null(conc)) {
           unlist(conc[c("evaluations", "inner_iterations", "solves", "failed", "decrement",
                         "max_residual")])
         })))
}

f10_ml <- function(sim, cell, seed) {
  fit <- f10_fit(sim, cell, seed, "the ML fit", polish = "newton",
                 polish_control = list(max_iter = 50L, tol_decrement = 1e-14))
  f10_result(fit, sim)
}

# A "blp" method with targets of `kind` (f10_targets()).
f10_blp <- function(kind) {
  force(kind)
  function(sim, cell, seed) {
    tg <- f10_targets(sim, cell, kind)
    fit <- f10_fit(sim, cell, seed, sprintf("the blp fit (%s shares)", kind),
                   concentrate_out_asc = "blp", target_shares = tg$shares,
                   target_vcov = tg$vcov)
    f10_result(fit, sim)
  }
}

f10_estimands <- function(fit, sim) NULL

# --- Findings: the article's table --------------------------------------------

F10_LABELS <- c(ml = "ML", blp_sample = "blp, sample shares", blp_exact = "blp, exact shares",
                blp_n10 = "blp, shares from n = 10 N", blp_n1 = "blp, shares from n = N")

# The §5.9 quantities of one design, estimator and group of parameters
# (columns `j`), on the reps both it and ML fitted: errors E (n x p) and
# ML's M, SEs se (hessian) and ser (robust). Each value with its MCSE by
# the delta method on per-rep influence terms (the SDs' and mean SEs'
# sampling error included; parameters within a rep correlated):
#   sd_vs_ml    mean_j SD_j / SD_ml,j
#   se_over_sd  mean_j mean(se_j) / SD_j
#   cover95     mean over reps and j of |err| <= 1.96 se (the old script's
#               1.96), and cover95_robust with the sandwich's
#   max_bias_sd max_j |mean err_j| / SD_j, the MCSE of the maximizing j's.
.f10_group <- function(E, M, se, ser) {
  n <- nrow(E)
  dev <- function(X) sweep(X, 2L, colMeans(X))
  a <- dev(E)^2
  b <- dev(M)^2
  abar <- colMeans(a)
  bbar <- colMeans(b)
  sd_j <- apply(E, 2L, stats::sd)
  mcse <- function(psi) stats::sd(psi) / sqrt(n)
  r <- sqrt(abar / bbar)
  psi_r <- drop(sweep(0.5 * (sweep(a, 2L, abar, "/") - sweep(b, 2L, bbar, "/")), 2L, r, "*") %*%
                  rep(1 / ncol(E), ncol(E)))
  sbar <- colMeans(se)
  q <- sbar / sd_j
  psi_q <- drop(sweep(sweep(se, 2L, sbar, "/") - 0.5 * sweep(a, 2L, abar, "/") - 0.5, 2L, q,
                      "*") %*% rep(1 / ncol(E), ncol(E)))
  cov_h <- rowMeans(abs(E) <= 1.96 * se)
  cov_r <- rowMeans(abs(E) <= 1.96 * ser)
  bias <- colMeans(E) / sd_j
  k <- which.max(abs(bias))
  psi_b <- sign(bias[k]) * ((E[, k] - mean(E[, k])) / sd_j[k] -
                              bias[k] * 0.5 * (a[, k] / abar[k] - 1))
  c(sd_vs_ml = mean(r), sd_vs_ml_mcse = mcse(psi_r),
    se_over_sd = mean(q), se_over_sd_mcse = mcse(psi_q),
    cover95 = mean(cov_h), cover95_mcse = mcse(cov_h),
    cover95_robust = mean(cov_r), cover95_robust_mcse = mcse(cov_r),
    max_bias_sd = abs(bias[[k]]), max_bias_sd_mcse = mcse(psi_b))
}

# The findings of F10's records (every cell and method): a list of
#   by_design  per design, estimator and group (theta_c, constants): fits
#              (paired with ML), failed reps, the §5.9 quantities with MCSEs,
#              and the work: median outer evaluations, inner evaluations per
#              solve, failed inner solves, reps that warned (a warning of
#              either start's fit: run_unit() records the method's
#              warnings, the kept fit's and the other start's alike; the
#              multi-start check muffles its refits'), the largest outer
#              decrement, the kept fits whose optimizer code is outside
#              nloptr's 1..4 (stops require_converged() accepted: ML's by
#              the polish, "blp"'s by the outer decrement), the median
#              seconds per rep (both starts),
#              and the reps where the theta0 start won;
#   article    per estimator, the ranges over the designs: the constants'
#              and theta_c's SD over ML's, the coverage (both groups), the
#              largest bias / SD of theta_c (the article's table).
f10_findings <- function(records) {
  recs <- Filter(function(r) is.list(r) && identical(r$family, "F10"), records)
  ok <- Filter(function(r) isTRUE(r$ok) && !is.null(r$est), recs)
  rows <- list()
  for (cl in intersect(F10_CELLS, unique(vapply(recs, `[[`, "", "cell")))) {
    ml <- Filter(function(r) r$cell == cl && r$method == "ml", ok)
    ml_reps <- vapply(ml, `[[`, 1L, "rep")
    for (m in intersect(names(F10_LABELS), unique(vapply(recs, `[[`, "", "method")))) {
      all_m <- Filter(function(r) r$cell == cl && r$method == m, recs)
      rm <- Filter(function(r) r$cell == cl && r$method == m && r$rep %in% ml_reps, ok)
      if (length(rm) < 2L) next
      reps <- vapply(rm, `[[`, 1L, "rep")
      mm <- ml[match(reps, ml_reps)]
      mat <- function(rs, col) do.call(rbind, lapply(rs, function(r) r$est[[col]]))
      E <- mat(rm, "est") - mat(rm, "truth")
      M <- mat(mm, "est") - mat(mm, "truth")
      se <- mat(rm, "se_hessian")
      ser <- mat(rm, "se_robust")
      asc <- rm[[1]]$est$block == "asc"
      dg <- lapply(rm, `[[`, "diag")
      conc <- function(f) vapply(dg, function(d) as.numeric(d$concentration[f] %||% NA), 1)
      start <- vapply(dg, function(d) d$start$winner %||% NA_character_, "")
      work <- data.frame(
        evals = stats::median(conc("evaluations")),
        inner_per_solve = stats::median(conc("inner_iterations") / conc("solves")),
        failed_solves = sum(conc("failed")),
        warned = sum(vapply(all_m, function(r) any(r$warnings$stage == "method"), TRUE)),
        max_decrement = suppressWarnings(max(vapply(dg, function(d) d$decrement %||% NA, 1),
                                             conc("decrement"), na.rm = TRUE)),
        stops_not_converged = sum(!vapply(dg, function(d) isTRUE(d$code %in% 1:4), TRUE)),
        seconds = stats::median(vapply(rm, function(r) r$diag$unit_sec %||% NA_real_, 1)),
        theta0_won = sum(start == "theta0", na.rm = TRUE))
      for (g in c("theta_c", "constants")) {
        j <- if (g == "constants") asc else !asc
        v <- .f10_group(E[, j, drop = FALSE], M[, j, drop = FALSE], se[, j, drop = FALSE],
                        ser[, j, drop = FALSE])
        rows[[length(rows) + 1L]] <- data.frame(
          design = cl, estimator = F10_LABELS[[m]], params = g, fits = length(rm),
          failed = length(all_m) - length(Filter(function(r) isTRUE(r$ok), all_m)),
          as.list(v), work, stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) return(list())
  by_design <- do.call(rbind, rows)
  rng <- function(x) sprintf("%.3f to %.3f", min(x), max(x))
  article <- do.call(rbind, lapply(intersect(F10_LABELS, by_design$estimator), function(e) {
    d <- by_design[by_design$estimator == e, , drop = FALSE]
    th <- d$params == "theta_c"
    data.frame(targets = e, constants_sd_vs_ml = rng(d$sd_vs_ml[!th]),
               theta_c_sd_vs_ml = rng(d$sd_vs_ml[th]), cover95 = rng(d$cover95),
               max_bias_sd_theta_c = rng(d$max_bias_sd[th]), designs = sum(th),
               stringsAsFactors = FALSE)
  }))
  list(by_design = by_design, article = article)
}

family <- list(
  id = "F10",
  claims = "C1",
  dgp = "inline (blp_concentration_mc.R)",
  cells = data.frame(cell = F10_CELLS, N = c(2000L, 500L, 2000L, 2000L), T = c(1L, 4L, 1L, 1L),
                     outside = c(TRUE, TRUE, FALSE, TRUE), markets = c(1L, 1L, 1L, 2L),
                     J = 10L, S = 100L, stringsAsFactors = FALSE),
  methods = list(ml = f10_ml, blp_sample = f10_blp("sample"), blp_exact = f10_blp("exact"),
                 blp_n10 = f10_blp("n10"), blp_n1 = f10_blp("n1")),
  methods_by_cell = NULL,
  R = list(mc = c(cross = 1000, panel = 1000, base = 1000, markets = 1000)),
  simulate = f10_simulate,
  estimands = f10_estimands,
  gates = NULL,
  findings = f10_findings
)
