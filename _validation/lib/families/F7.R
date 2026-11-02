# ============================================================================
# F7: population scale (plan §5; claims C1, C6) on the two population designs
# of lib/dgp_geo.R, one size each, data from the port (engine "cpp"):
#   claims  D3-claims: 625,000 patients, about 10^6 admissions (a panel,
#           three in four admitted once), sets of 3 to 60 of J = 2,000
#           hospitals within 40 km (3.9 x 10^7 rows), no outside option;
#   census  D4-census: 10^6 students in 150 markets, each student's 15
#           nearest of the market's 20 schools (J = 3,000; 1.5 x 10^7 rows)
#           and an outside option.
# Each design is simulated per model: MNL data (cells "claims.mnl",
# "census.mnl"; variant "mnl") and MXL data (cells "claims", "census";
# variant "default"). Alternatives never chosen are dropped (plan §4;
# diag$drop). Methods:
#   mnl  run_mnlogit() with every constant, polish = "newton" (both designs);
#        route hessian.
#   mxl  run_mxlogit() with every constant, S = 100 generate-mode draws seeded
#        from the rep, polish = "newton", two starts (fit_two_starts()): the
#        claims panel (person_col = "patient": the patient's tastes shared
#        by their admissions; log-normal -distance and normal teaching x
#        severity with means), the census cross-section (D4's rc_mean =
#        FALSE design); route hessian; multi-start on reps 10, 20, ...
#   blp  census only: concentrate_out_asc = "blp" with the sample's shares
#        (R = 25, the same reps' data as mxl).
# Threads come from the runner (mc.R --threads; 16 in the mc tier). Peak
# RSS at one thread (2026-11-02, macOS): claims 8.5 GB (MNL rep) and 8.1 GB
# (an MXL fit's preparation, Hessian and variance), census 4.2 and 4.0 GB.
# Derived, every method (truth at_truth()'s, S x 4 for the MXL; no SE):
#   claims  closing hospital F7_CLOSE_HOSP (the 10th largest): the diversion
#           to the hospital taking most of its patients at theta0, and the
#           change in consumer surplus per admission (utils; logsum());
#   census  closing school f7_schools()$close (market 1's largest
#           constant): the diversion to its first recipient at theta0 and
#           to the outside option, and the change in consumer surplus per
#           student; raising the quality of f7_schools()$quality (market
#           1's second largest) by one SD of quality across schools: the
#           change in its share of all students and in consumer surplus.
# Recipients at theta0 come from geo_recipients() on 40,000 units, once per
# process; diag$asc_expected holds each constant's expected count in the rep
# (pooled_tails() by expected-count tercile: the entry's pooled_groups()).
# Gates: mnl as F1 (C1, with the constants' incidental effects absent from
# the slopes, C6): the slopes' bias, Hessian SE ratio and tails, and the
# joint Wald size on the slopes, both designs. Findings: the MXL's bias/ESE,
# SE ratio and coverage; pooled constant coverage; the counterfactuals'
# bias/ESE; time and peak memory per fit (diag).
# ============================================================================

F7_CLOSE_HOSP <- 10L

.f7 <- new.env()

# The census's target schools: market 1's largest and second largest
# constants at theta0 (the MXL variant's), and one SD of quality across
# schools; computed once per process from the fixed world.
f7_schools <- function() {
  if (is.null(.f7$schools)) {
    dgp <- get_dgp("D4-census")
    w <- geo_schools(dgp$design, geo_draw(dgp$design, 0L, 1L)$world)
    th <- dgp_theta0(dgp, "default")
    m1 <- w$alt[w$market == 1L]
    top <- m1[order(th[paste0("ASC_", m1)], decreasing = TRUE)][1:2]
    .f7$schools <- list(close = top[1], quality = top[2], sd_q = stats::sd(w$quality))
  }
  .f7$schools
}

# The first recipient (and, with an outside option, the outside option) of
# the closed alternative at theta0, per design and variant, once per process.
f7_recipients <- function(design, variant) {
  key <- paste(design, variant)
  if (is.null(.f7[[key]])) {
    if (design == "claims") {
      sim <- get_dgp("D3-claims")$simulate(N = 40000L, seed = 1L, variant = variant,
                                           engine = "cpp")
      d <- geo_recipients(sim, F7_CLOSE_HOSP)
    } else {
      sim <- get_dgp("D4-census")$simulate(N = 40000L, seed = 1L, variant = variant,
                                           engine = "cpp")
      d <- geo_recipients(sim, f7_schools()$close)
    }
    .f7[[key]] <- c(names(d)[names(d) != "outside"][1L], if ("outside" %in% names(d)) "outside")
  }
  .f7[[key]]
}

f7_simulate <- function(cell, seed) {
  geo_drop(get_dgp(cell$dgp)$simulate(N = cell$N, seed = seed, variant = cell$variant,
                                      engine = "cpp"))
}

# The counterfactuals (header) at the fit and at theta0.
f7_estimands <- function(fit, sim) {
  claims <- sim$meta$design$name == "D3-claims"
  fit0 <- at_truth(fit, sim$truth)
  n <- if (claims) sim$meta$A else sim$meta$N
  if (claims) {
    h <- F7_CLOSE_HOSP
    to <- f7_recipients("claims", sim$meta$variant)
  } else {
    sc <- f7_schools()
    h <- sc$close
    to <- f7_recipients("census", sim$meta$variant)
  }
  cl <- geo_closure(fit, sim$data, h, to, n)
  cl0 <- geo_closure(fit0, sim$data, h, to, n)
  out <- data.frame(name = c(paste0("div_to_", to), "dcs_close"),
                    est = c(unname(cl$div), cl$dcs), se = NA_real_,
                    truth = c(unname(cl0$div), cl0$dcs), stringsAsFactors = FALSE)
  if (!claims) {
    ses <- sim$meta$students$ses
    q <- geo_quality(fit, sim$data, sc$quality, sc$sd_q, ses, n)
    q0 <- geo_quality(fit0, sim$data, sc$quality, sc$sd_q, ses, n)
    out <- rbind(out, data.frame(name = c("dshare_q", "dcs_q"), est = c(q$dshare, q$dcs),
                                 se = NA_real_, truth = c(q0$dshare, q0$dcs),
                                 stringsAsFactors = FALSE))
  }
  out
}

.f7_diag <- function(fit, sim) {
  diag <- fit_diag(fit)
  diag$drop <- sim$meta$drop
  diag$asc_expected <- geo_asc_expected(fit, sim)
  diag$start <- attr(fit, "start")
  diag$rows <- nrow(sim$data)
  diag
}

.f7_out <- function(fit, sim, gated, diag) {
  V <- fit_vcovs(fit, "hessian")
  cf <- coef(fit)
  common <- if (gated) names(cf)[fit$param_map$beta] else names(cf)[-as.integer(fit$param_map$asc)]
  list(est = est_table(fit, sim$truth, V = V, gated = if (gated) "beta" else FALSE),
       vcov_c = vcov_common(fit, V = V$hessian),
       joint_wald = joint_wald(cf[common], sim$truth, V$hessian),
       derived = f7_estimands(fit, sim),
       diag = diag)
}

f7_mnl <- function(sim, cell, seed) {
  x <- if (cell$dgp == "D3-claims") D3_X$mnl else D4_X
  fit <- require_converged(choicer::run_mnlogit(
    sim$data, "id", "alt", "choice", x, include_outside_option = sim$meta$spec$outside,
    polish = "newton"), "the MNL fit")
  .f7_out(fit, sim, TRUE, .f7_diag(fit, sim))
}

# One MXL fit: the claims panel or the census cross-section.
f7_run <- function(sim, cell, seed, ...) {
  if (cell$dgp == "D3-claims") {
    choicer::run_mxlogit(sim$data, "id", "alt", "choice", D3_X$default, D3_W,
                         rc_dist = c(1L, 0L), rc_mean = TRUE, person_col = "patient",
                         S = cell$S, draws = "generate", seed = seed + 1L, ...)
  } else {
    choicer::run_mxlogit(sim$data, "id", "alt", "choice", D4_X, D4_W, rc_dist = c(1L, 0L),
                         rc_mean = FALSE, include_outside_option = TRUE, S = cell$S,
                         draws = "generate", seed = seed + 1L, ...)
  }
}

f7_mxl <- function(sim, cell, seed) {
  fit <- fit_two_starts(function(theta_init) {
    f7_run(sim, cell, seed, polish = "newton", theta_init = theta_init)
  }, sim$truth, "the MXL fit")
  diag <- .f7_diag(fit, sim)
  if (vs_seed_rep("F7", cell$cell, seed) %% 10 == 0) {
    diag$multistart <- multistart(fit, refit = function(theta_init) {
      require_converged(f7_run(sim, cell, seed, polish = "newton", theta_init = theta_init),
                        "the MXL fit")
    }, seed = seed)
  }
  .f7_out(fit, sim, FALSE, diag)
}

f7_blp <- function(sim, cell, seed) {
  fit <- fit_two_starts(function(theta_init) {
    f7_run(sim, cell, seed, concentrate_out_asc = "blp", theta_init = theta_init)
  }, sim$truth, "the blp fit")
  diag <- .f7_diag(fit, sim)
  diag$concentration <- unlist(fit$optimizer$concentration[c("evaluations", "solves", "failed",
                                                             "max_residual")])
  .f7_out(fit, sim, FALSE, diag)
}

family <- list(
  id = "F7",
  claims = c("C1", "C6"),
  dgp = "D3-claims, D4-census",
  cells = data.frame(cell = c("claims.mnl", "claims", "census.mnl", "census"),
                     N = c(625000L, 625000L, 1000000L, 1000000L),
                     J = c(2000L, 2000L, 3000L, 3000L), S = c(NA, 100L, NA, 100L),
                     dgp = rep(c("D3-claims", "D4-census"), each = 2L),
                     variant = rep(c("mnl", "default"), 2L), stringsAsFactors = FALSE),
  methods = list(mnl = f7_mnl, mxl = f7_mxl, blp = f7_blp),
  methods_by_cell = list(claims.mnl = "mnl", claims = "mxl", census.mnl = "mnl",
                         census = c("mxl", "blp")),
  R = list(mc = c(claims.mnl = 1000, claims = 100, census.mnl = 1000, census = 100)),
  R_by_method = list(census = c(blp = 25)),
  simulate = f7_simulate,
  estimands = f7_estimands,
  pooled_groups = geo_count_groups,
  gates = list(cell = list(mnl = c("claims.mnl", "census.mnl")),
               metrics = list(mnl = c("bias_std@beta", "se_ratio:hessian@beta",
                                      "tail_lower:hessian@beta", "tail_upper:hessian@beta",
                                      "joint_wald_size")))
)
