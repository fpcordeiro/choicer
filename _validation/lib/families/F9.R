# ============================================================================
# F9: the Bayesian samplers (plan §1 C7, §2 "Bayesian", §5 F9): MNP
# (run_mnprobit()) on D2's design, HMNL (run_hmnlogit()) and HMNP
# (run_hmnprobit()) on D1-L; data from DGP "DB" (lib/dgp_bayes.R). Sourced by
# load_family() into a child of the library's environment.
#
# SBC cells (gated; C7's exactness at any N): each rep draws the truth from
# DB_PRIOR$sbc, simulates, runs the sampler with the same prior (one chain:
# `burn` iterations, then `post` kept every `thin`-th), and ranks every
# identified quantity and the log-likelihood among 99 equally spaced kept
# draws (sbc_diag(): diag$sbc). Gate: sbc_ecdf on every quantity, Holm over
# the three cells (plan §2). Chains: 5,000 + 30,000 kept every third (10,000
# draws), not the plan's 5,000 + 10,000: at 10,000 straight draws the
# smallest bulk ESS was 73-167 (MNP, N = 300), 96-117 (HMNL) and 134-173
# (HMNP), at plan §2's floor of 99; with the thinning, 204-544 for the SBC
# cells (2026-11-02, 1-4 reps each).
# BvM cells (findings; C7's Bernstein-von Mises): theta0 fixed (DB_THETA0),
# prior DB_PRIOR$bvm. est holds each N-indexed quantity (MNP beta and Sigma,
# HMNL b, W and delta) with est = posterior mean and se_hessian = posterior
# SD, so that se_ratio is posterior SD / ESE and the tails are Wald-type's;
# derived holds the equal-tailed 95% credible interval of each as
# "cri:<quantity>", est its midpoint and se its half-width over
# qnorm(0.975), so that metrics_cell()'s tails (truth < est - z se, truth >
# est + z se) are exactly the credible interval's. joint_wald: the common
# quantities (no constants) against the posterior covariance. theta and
# sigma_d2 gain information only with J: SBC only (C7). diag$mcmc holds
# each quantity's bulk and tail ESS (choicer::ess(): ess_bulk, ess_tail)
# and their minima (ess_min, ess_tail_min).
# BvM chains are lengthened until the smallest tail ESS n is about 400. An
# empirical 2.5% quantile errs by about 2.67 / sqrt(n) posterior SDs, which
# widens each credible-interval tail from 2.5% to Phi(-1.96 / sqrt(1 +
# 2.67^2 / n)): 2.60% at n = 400, 2.70% at 200, 2.91% at 100 (the tolerance
# is 0.75 pp). Smallest tail ESS over the quantities, reps 1-3 (2026-11-02):
#   bvm_mnp_500   5k + 30k / 3:  576-738
#   bvm_mnp_2000  5k + 48k / 4:  562-601  (318-382 at 30k / 3)
#   bvm_hmnl_200  5k + 60k / 6:  488-560  (209-406 at 30k / 3)
#   bvm_hmnl_800  5k + 90k / 9:  466, 201, 677  (47-266 at 30k / 3)
# always at a Sigma or W entry, W[feature,feature] in the HMNL (the other
# HMNL quantities above 1,600 at 800 persons). At 800 persons a rep's
# W[feature,feature] tails can thus run about 0.2 pp wide, and the other
# quantities' less than 0.05 pp; 400 in every rep would take about 180k
# iterations (about 250 s a rep, 35 CPU-h for the cell's 500 reps).
# Quick cells (q_*; plan Q4): the SBC cells with fewer persons or situations
# and short chains, R = 100; q4_check() (lib/dgp_bayes.R) applies Holm at
# 0.001 to their sbc_ecdf rows.
# The HB samplers' fit-time convergence warning (rank-R-hat > 1.01 or bulk
# ESS < 400) is recorded by run_unit(), and its rate is diag$converged's.
# ============================================================================

f9_cells <- function() {
  data.frame(
    cell  = c("sbc_mnp", "sbc_hmnl", "sbc_hmnp", "bvm_mnp_500", "bvm_mnp_2000", "bvm_hmnl_200",
              "bvm_hmnl_800", "q_sbc_mnp", "q_sbc_hmnl", "q_sbc_hmnp"),
    kind  = c("sbc", "sbc", "sbc", "bvm", "bvm", "bvm", "bvm", "sbc", "sbc", "sbc"),
    model = c("mnp", "hmnl", "hmnp", "mnp", "mnp", "hmnl", "hmnl", "mnp", "hmnl", "hmnp"),
    N     = c(300L, 100L, 100L, 500L, 2000L, 200L, 800L, 200L, 50L, 50L),
    T     = c(1L, 8L, 8L, 1L, 1L, 8L, 8L, 1L, 8L, 8L),
    burn  = c(rep(5000L, 7L), 1000L, 2000L, 2000L),
    post  = c(30000L, 30000L, 30000L, 30000L, 48000L, 60000L, 90000L, 8000L, 8000L, 8000L),
    thin  = c(3L, 3L, 3L, 3L, 4L, 6L, 9L, 2L, 2L, 2L),
    stringsAsFactors = FALSE)
}

F9_L <- 99L
F9_VARIANT <- c(mnp = "default", hmnl = "hmnl", hmnp = "hmnp")

f9_simulate <- function(cell, seed) {
  prior <- if (cell$kind == "sbc") DB_PRIOR$sbc[[cell$model]]
  db_simulate(cell$N, seed, F9_VARIANT[[cell$model]], prior = prior, T = cell$T)
}

# The quantities BvM records: every identified one for MNP; b, W and delta
# for HMNL.
f9_bvm_quantities <- function(nm) nm[db_block(nm) != "other"]

# Derived estimands: each quantity's equal-tailed 95% credible interval as
# "cri:<quantity>", est its midpoint and se its half-width / qnorm(0.975).
f9_cri <- function(draws, truth) {
  ci <- apply(draws, 2L, stats::quantile, c(0.025, 0.975), names = FALSE)
  data.frame(name = paste0("cri:", colnames(draws)), est = (ci[1, ] + ci[2, ]) / 2,
             se = (ci[2, ] - ci[1, ]) / (2 * stats::qnorm(0.975)),
             truth = unname(truth[colnames(draws)]), stringsAsFactors = FALSE)
}

f9_estimands <- function(fit, sim) {
  D <- db_draws(fit, sim$meta$model)
  f9_cri(D[, f9_bvm_quantities(colnames(D)), drop = FALSE], sim$truth)
}

# The log-likelihood at the truth and at the thinned draws `idx` (SBC's test
# quantity): MNP by GHK with uniforms fixed per rep (subseed 4); HMNL's at
# every kept draw from the sampler's trace; HMNP's from the kept taste draws.
f9_loglik <- function(model, fit, D, sim, seed, idx) {
  if (model == "mnp") {
    m <- sim$meta
    U <- vs_with_seed(db_subseed(seed, 4L),
                      array(stats::runif(m$N * DB_GHK_DRAWS * DB_P), c(m$N, DB_GHK_DRAWS, DB_P)))
    ll0 <- db_mnp_loglik(sim$truth[DB_BETA], .db_unvech(sim$truth[.DB_SIGMA_NAMES]), m$Xd, m$y, U)
    ll <- vapply(idx, function(k) {
      db_mnp_loglik(fit$draws$beta[k, ], .db_unvech(fit$draws$sigma[k, ]), m$Xd, m$y, U)
    }, numeric(1))
    return(list(truth = ll0, draws = ll))
  }
  lay <- db_hb_layout(sim$data)
  delta0 <- sim$meta$delta[DB_HB_ALTS]
  if (model == "hmnl") {
    return(list(truth = db_hmnl_loglik(lay, sim$meta$tastes, delta0),
                draws = fit$chains[[1L]]$loglik))
  }
  B <- db_beta_i_draws(fit, lay, idx)
  dl <- D[idx, paste0("delta_", DB_HB_ALTS), drop = FALSE]
  list(truth = db_hmnp_loglik(lay, sim$meta$tastes, delta0),
       draws = vapply(seq_along(idx), function(k) db_hmnp_loglik(lay, B[, , k], dl[k, ]), 1))
}

# A method: the model's sampler on the rep's data; SBC cells return diag$sbc,
# BvM cells the est, derived, joint_wald and vcov_c of the posterior.
f9_method <- function(model) {
  force(model)
  function(sim, cell, seed) {
    sbc <- cell$kind == "sbc"
    mcmc <- list(R = cell$burn + cell$post, burn = cell$burn, thin = cell$thin,
                 seed = db_subseed(seed, 3L))
    warned <- FALSE
    t0 <- proc.time()[["elapsed"]]
    fit <- withCallingHandlers(
      db_fit(model, sim$data, DB_PRIOR[[cell$kind]][[model]], mcmc,
             keep_beta_i = if (sbc && model == "hmnp") "draws" else "none"),
      warning = function(w) {
        if (startsWith(conditionMessage(w), "Convergence check failed")) warned <<- TRUE
      })
    sec <- proc.time()[["elapsed"]] - t0
    D <- db_draws(fit, model)
    diag <- list(code = NA_integer_, converged = if (model == "mnp") NA else !warned,
                 iterations = mcmc$R, time_sec = sec,
                 mcmc = list(burn = cell$burn, thin = cell$thin, kept = nrow(D),
                             accept_beta = fit$accept$mean_beta %||% NA_real_,
                             accept_delta = fit$accept$mean_delta %||% NA_real_))
    if (sbc) {
      truth <- sim$truth[intersect(names(sim$truth), colnames(D))]
      idx <- sbc_thin_index(nrow(D), F9_L)
      ll <- f9_loglik(model, fit, D, sim, seed, idx)
      diag$sbc <- sbc_diag(truth, D, L = F9_L, loglik = ll$truth, loglik_draws = ll$draws)
      return(list(diag = diag))
    }
    q <- f9_bvm_quantities(colnames(D))
    D <- D[, q, drop = FALSE]
    m <- colMeans(D)
    V <- stats::cov(D)
    blk <- db_block(q)
    com <- q[blk != "asc"]
    e <- choicer::ess(D)
    diag$mcmc[c("ess_bulk", "ess_tail")] <- list(e[, "bulk"], e[, "tail"])
    diag$mcmc$ess_min <- min(e[, "bulk"], na.rm = TRUE)
    diag$mcmc$ess_tail_min <- min(e[, "tail"], na.rm = TRUE)
    list(est = data.frame(name = q, block = blk, est = unname(m), truth = unname(sim$truth[q]),
                          se_hessian = sqrt(diag(V)), se_bhhh = NA_real_, se_robust = NA_real_,
                          se_cluster = NA_real_, gated = FALSE, stringsAsFactors = FALSE),
         vcov_c = list(route = "posterior", common = V[com, com, drop = FALSE],
                       asc_var = diag(V)[q[blk == "asc"]]),
         joint_wald = joint_wald(m[com], sim$truth, V[com, com, drop = FALSE]),
         derived = f9_cri(D, sim$truth),
         diag = diag)
  }
}

f9_sbc_cells <- c("sbc_mnp", "sbc_hmnl", "sbc_hmnp")
f9_quick_cells <- c("q_sbc_mnp", "q_sbc_hmnl", "q_sbc_hmnp")
f9_cell_tab <- f9_cells()

family <- list(
  id = "F9",
  claims = "C7",
  dgp = "DB",
  cells = f9_cell_tab,
  methods = list(mnp = f9_method("mnp"), hmnl = f9_method("hmnl"), hmnp = f9_method("hmnp")),
  methods_by_cell = stats::setNames(as.list(f9_cell_tab$model), f9_cell_tab$cell),
  R = list(quick = stats::setNames(rep(100, 3L), f9_quick_cells),
           mc = c(sbc_mnp = 1000, sbc_hmnl = 1000, sbc_hmnp = 1000, bvm_mnp_500 = 500,
                  bvm_mnp_2000 = 500, bvm_hmnl_200 = 500, bvm_hmnl_800 = 500)),
  simulate = f9_simulate,
  estimands = f9_estimands,
  gates = list(cell = f9_sbc_cells, metrics = "sbc_ecdf"),
  gates_quick = list(cell = f9_quick_cells, metrics = "sbc_ecdf")
)
