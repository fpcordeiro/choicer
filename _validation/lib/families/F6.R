# ============================================================================
# F6: many alternatives (plan §5; claim C6) on D3 (lib/dgp_geo.R): hospital
# choice among the hospitals within 30 km of each patient (3 to 30 of J =
# 200, mean 21.5), no outside option, the largest hospital the reference, a
# cross-section. Cells N in {5k, 20k, 80k} patients (N/J = 25, 100, 400;
# J/sqrt(N) = 2.83, 1.41, 0.71), each simulated twice: MNL data (cells
# "5k.mnl", ...; variant "mnl") for the MNL, whose claim needs a correct
# model, and MXL data (cells "5k", ...; variant "default") for the MXL and
# "blp", paired. Data from the port (engine "cpp", identical to the R
# reference); hospitals never chosen in a rep are dropped (plan §4,
# drop_unbounded_alts(); counts in diag$drop: about 10 of 200 expected at
# 5k, 1.5 at 20k, 0.1 at 80k; checks/wp6_geo.R --only=calib).
# Methods:
#   mnl  run_mnlogit() with every constant, polish = "newton"; routes hessian
#        and robust (the sandwich, ungated: biased down in finite samples
#        with many constants, Pfaffermayr 2019).
#   mxl  run_mxlogit() with every constant: a log-normal taste on -distance
#        and a normal one on teaching x severity, both with means (rc_mean =
#        TRUE: both vary by patient within a hospital), uncorrelated; S = 200
#        generate-mode draws seeded from the rep; polish = "newton"; from
#        choicer's default start and from theta0 (fit_two_starts(), record in
#        diag$start); route hessian; multi-start on reps 10, 20, ...
#   blp  the same with concentrate_out_asc = "blp" and the sample's shares
#        (polish "none"; one market: the sets connect every hospital), at
#        20k only (R = 500, the same reps' data as mxl).
# Derived, every method: the diversion from hospital F6_CLOSE (the 10th
# largest), closed, to the three hospitals taking most of its patients at
# theta0 (geo_recipients() on a simulation of 40,000 patients, per variant,
# once per process): D_hk = (S_k without h - S_k) / S_h over the patients
# offered h, predict(newdata, type = "shares") (geo_closure()), and the
# closure's change in consumer surplus per patient (dcs_close, utils: the
# change in logsum() summed over those patients, over N), truths
# at_truth()'s (S x 4 for the MXL); choicer gives no SE. diag$asc_expected
# holds each constant's expected count at theta0 in the rep, for
# pooled_tails() by expected-count tercile (the entry's pooled_groups()).
# Gates (C6: the MNL is a two-way Poisson model, Baker 1994; its slopes
# carry no leading-order incidental bias): mnl at every rung, the slopes'
# bias and Hessian SE ratio. Findings: the MNL's robust (sandwich) SE ratio
# below 1 at small N/J (C6, Pfaffermayr 2019); the MXL's N bias / J and
# bias/ESE against J/sqrt(N) (S1's probe on D4: about -0.16 J/sqrt(N));
# whether "blp" reduces it; pooled constant coverage by expected-count
# tercile; the closure's diversion and consumer-surplus bias (O(J/N) in the
# MNL, C6).
# ============================================================================

F6_CLOSE <- 10L

.f6 <- new.env()

# The three largest recipients of F6_CLOSE's patients at theta0, per
# variant, computed once per process.
f6_recipients <- function(variant) {
  if (is.null(.f6[[variant]])) {
    sim <- get_dgp("D3")$simulate(N = 40000L, seed = 1L, variant = variant, engine = "cpp")
    .f6[[variant]] <- names(geo_recipients(sim, F6_CLOSE))[1:3]
  }
  .f6[[variant]]
}

f6_simulate <- function(cell, seed) {
  geo_drop(get_dgp("D3")$simulate(N = cell$N, seed = seed, variant = cell$variant,
                                  engine = "cpp"))
}

# The closure's diversion to the three recipients and its change in
# consumer surplus, at the fit and at theta0.
f6_estimands <- function(fit, sim) {
  to <- f6_recipients(sim$meta$variant)
  est <- geo_closure(fit, sim$data, F6_CLOSE, to, sim$meta$A)
  tru <- geo_closure(at_truth(fit, sim$truth), sim$data, F6_CLOSE, to, sim$meta$A)
  data.frame(name = c(paste0("div_to_", to), "dcs_close"), est = c(unname(est$div), est$dcs),
             se = NA_real_, truth = c(unname(tru$div), tru$dcs), stringsAsFactors = FALSE)
}

# The record's diagnostics: fit_diag(), the counts dropped, the constants'
# expected counts, and the two starts' record (MXL).
.f6_diag <- function(fit, sim) {
  diag <- fit_diag(fit)
  diag$drop <- sim$meta$drop
  diag$asc_expected <- geo_asc_expected(fit, sim)
  diag$start <- attr(fit, "start")
  diag
}

.f6_out <- function(fit, sim, gated, diag, routes = "hessian") {
  V <- fit_vcovs(fit, routes)
  cf <- coef(fit)
  list(est = est_table(fit, sim$truth, V = V, gated = gated),
       vcov_c = vcov_common(fit, V = V$hessian),
       joint_wald = joint_wald(cf[-as.integer(fit$param_map$asc)], sim$truth, V$hessian),
       derived = f6_estimands(fit, sim),
       diag = diag)
}

f6_mnl <- function(sim, cell, seed) {
  fit <- require_converged(choicer::run_mnlogit(sim$data, "id", "alt", "choice", D3_X$mnl,
                                                polish = "newton"), "the MNL fit")
  .f6_out(fit, sim, "beta", .f6_diag(fit, sim), c("hessian", "robust"))
}

f6_run <- function(sim, cell, seed, ...) {
  choicer::run_mxlogit(sim$data, "id", "alt", "choice", D3_X$default, D3_W,
                       rc_dist = c(1L, 0L), rc_mean = TRUE, S = cell$S, draws = "generate",
                       seed = seed + 1L, ...)
}

f6_mxl <- function(sim, cell, seed) {
  fit <- fit_two_starts(function(theta_init) {
    f6_run(sim, cell, seed, polish = "newton", theta_init = theta_init)
  }, sim$truth, "the MXL fit")
  diag <- .f6_diag(fit, sim)
  if (vs_seed_rep("F6", cell$cell, seed) %% 10 == 0) {
    diag$multistart <- multistart(fit, refit = function(theta_init) {
      require_converged(f6_run(sim, cell, seed, polish = "newton", theta_init = theta_init),
                        "the MXL fit")
    }, seed = seed)
  }
  .f6_out(fit, sim, FALSE, diag)
}

f6_blp <- function(sim, cell, seed) {
  fit <- fit_two_starts(function(theta_init) {
    f6_run(sim, cell, seed, concentrate_out_asc = "blp", theta_init = theta_init)
  }, sim$truth, "the blp fit")
  diag <- .f6_diag(fit, sim)
  diag$concentration <- unlist(fit$optimizer$concentration[c("evaluations", "solves", "failed",
                                                             "max_residual")])
  .f6_out(fit, sim, FALSE, diag)
}

.f6_cells <- function() {
  N <- c(5000L, 20000L, 80000L)
  d <- data.frame(cell = c("5k", "20k", "80k", "5k.mnl", "20k.mnl", "80k.mnl"), N = rep(N, 2L),
                  J = 200L, S = rep(c(200L, NA_integer_), each = 3L),
                  variant = rep(c("default", "mnl"), each = 3L), stringsAsFactors = FALSE)
  d$N_over_J <- d$N / d$J
  d$J_over_sqrtN <- d$J / sqrt(d$N)
  d
}

family <- list(
  id = "F6",
  claims = "C6",
  dgp = "D3",
  cells = .f6_cells(),
  methods = list(mnl = f6_mnl, mxl = f6_mxl, blp = f6_blp),
  methods_by_cell = list("5k" = "mxl", "20k" = c("mxl", "blp"), "80k" = "mxl",
                         "5k.mnl" = "mnl", "20k.mnl" = "mnl", "80k.mnl" = "mnl"),
  R = list(mc = c("5k" = 1000, "20k" = 1000, "80k" = 500,
                  "5k.mnl" = 1000, "20k.mnl" = 1000, "80k.mnl" = 500)),
  R_by_method = list("20k" = c(blp = 500)),
  simulate = f6_simulate,
  estimands = f6_estimands,
  pooled_groups = geo_count_groups,
  gates = list(cell = list(mnl = c("5k.mnl", "20k.mnl", "80k.mnl")),
               metrics = list(mnl = c("bias_std@beta", "se_ratio:hessian@beta")))
)
