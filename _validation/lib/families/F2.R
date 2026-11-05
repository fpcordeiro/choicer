# ============================================================================
# F2: nested logit (plan §5; claims C1, C5) on D4-40's NL variant: school
# choice among each student's 15 nearest of one market's 40 schools (40
# constants; S1, 2026-11-02: 40 schools, not D4-small's 160, so that the
# gates test C1 and not the incidental constants, C6; J/sqrt(N) = 1.26,
# 0.63, 0.32 on the rungs), nests by sector (lambda: public .5, charter .7,
# private .85), the outside option a degenerate nest. Cells N in {1k, 4k,
# 16k} students. Method `ml`: run_nestlogit() with polish = "newton" (the
# optimizer's stop certified, as Q2 asks), routes hessian, bhhh, robust;
# derived: wtp_q_ses = -b_qs / b_d, wtp()'s ratio of the quality x SES and
# distance coefficients: the extra km a student one SD above the mean SES
# travels for one SD of school quality, in the mean-SES student's distance
# disutility (quality's main effect is in the constants, so the SES gradient
# is the quality WTP the model identifies). Multi-start on reps 10, 20, ...
# Never-chosen alternatives are dropped in every rep (drop_unbounded_alts(),
# plan §4's policy; counts in diag$drop). Every fit, multi-start refits
# included, must converge (require_converged(), lib/fit.R).
# Gates at 16k as F1's: bias, Hessian SE ratio and tails for beta, lambda and
# the constants whose expected count at the cell is at least
# F2_ASC_GATE_COUNT (a constant's bias and skewness are O(1/n_j), so those of
# rarely chosen schools are findings; 38 of 40 at 16k); joint Wald size on
# beta and lambda; the WTP's tails.
# ============================================================================

F2_ASC_GATE_COUNT <- 50

.f2 <- new.env()

# Each school's expected share at theta0 (mean probability over 40,000
# students of one seed), computed once per process.
f2_shares <- function() {
  if (is.null(.f2$shares)) {
    .f2$shares <- get_dgp("D4-40")$simulate(N = 40000L, seed = 1L, variant = "nl")$meta$shares
  }
  .f2$shares
}

f2_simulate <- function(cell, seed) {
  geo_drop(get_dgp("D4-40")$simulate(N = cell$N, seed = seed, variant = "nl"))
}

f2_fit <- function(sim, theta_init = NULL) {
  require_converged(choicer::run_nestlogit(sim$data, "id", "alt", "choice", D4_X,
                                           nest_col = "sector", include_outside_option = TRUE,
                                           theta_init = theta_init, polish = "newton"),
                    "the NL fit")
}

# The WTP for quality x SES in distance units; its truth is wtp() at theta0.
f2_estimands <- function(fit, sim) {
  w <- choicer::wtp(fit, price_var = "dist", attr_vars = "q_ses")
  w0 <- choicer::wtp(at_truth(fit, sim$truth), price_var = "dist", attr_vars = "q_ses")
  data.frame(name = "wtp_q_ses", est = w$Estimate, se = w$Std_Error, truth = w0$Estimate,
             stringsAsFactors = FALSE)
}

f2_ml <- function(sim, cell, seed) {
  fit <- f2_fit(sim)
  V <- fit_vcovs(fit, c("hessian", "bhhh", "robust"))
  cf <- coef(fit)
  common <- names(cf)[-as.integer(fit$param_map$asc)]
  big <- names(which(cell$N * f2_shares() >= F2_ASC_GATE_COUNT))
  gated <- c("beta", "lambda", intersect(paste0("ASC_", big), names(cf)))
  diag <- c(fit_diag(fit), list(drop = sim$meta$drop))
  if (vs_seed_rep("F2", cell$cell, seed) %% 10 == 0) {
    diag$multistart <- multistart(fit, refit = function(theta_init) f2_fit(sim, theta_init),
                                  seed = seed)
  }
  list(est = est_table(fit, sim$truth, V = V, gated = gated),
       vcov_c = vcov_common(fit, V = V$hessian),
       joint_wald = joint_wald(cf[common], sim$truth, V$hessian),
       derived = f2_estimands(fit, sim),
       diag = diag)
}

family <- list(
  id = "F2",
  claims = c("C1", "C5"),
  dgp = "D4-40",
  cells = data.frame(cell = c("1k", "4k", "16k"), N = c(1000L, 4000L, 16000L),
                     stringsAsFactors = FALSE),
  methods = list(ml = f2_ml),
  methods_by_cell = NULL,
  R = list(quick = c("1k" = 100), mc = c("1k" = 5000, "4k" = 5000, "16k" = 5000)),
  simulate = f2_simulate,
  estimands = f2_estimands,
  gates = list(cell = "16k",
               metrics = list(ml = c("bias_std", "se_ratio:hessian", "tail_lower:hessian",
                                     "tail_upper:hessian", "joint_wald_size",
                                     "tail_lower:hessian@derived",
                                     "tail_upper:hessian@derived"))),
  # The quick tier at 1k (S1; J/sqrt(N) = 1.26), on the family's seeds, 2026-11-05: R = 300 and, the
  # cell being cheap, R = 3000 (below; MCSE 0.018 for bias, 0.013 for an SE ratio, 0.3 pp for a
  # tail); reps 1-100, Q3's own, pass every row. Kept: the bias of dist_ses and q_ses (0.027, 0.060),
  # the slopes' Hessian SE ratio (0.992-0.994) and tails (1.7-3.5%), the WTP's SE ratio (0.999).
  # Dropped: dist's bias, -0.111 (R = 300: -0.125 +- 0.058; C6 predicted about -0.2 from the mixed
  # logit's -0.16 J/sqrt(N)), which alone fails about 1% of R = 100 runs; lambda (C5, C6): Lambda_1's
  # bias -0.101, upper tails 8.7%, 6.2%, 4.4%, SE ratios down to 0.953, fails about 2% of runs; the
  # joint Wald size on beta and lambda, 9.7% +- 0.5; the constants (about 20 choices each at 1k,
  # O(1/n_j): the rarest one's bias 0.64 and SE ratio 1.47 at R = 300). Q1.2 and Q2 check lambda.
  # False failures at R = 100, 1,000 bootstrap samples of the 3000: 0.6%, recentered 1.1%.
  gates_quick = list(cell = "1k",
                     metrics = list(ml = c("bias_std@dist_ses", "bias_std@q_ses",
                                           "se_ratio:hessian@beta", "tail_lower:hessian@beta",
                                           "tail_upper:hessian@beta", "se_ratio:hessian@derived")))
)
