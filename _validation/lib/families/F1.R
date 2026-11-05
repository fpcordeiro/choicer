# ============================================================================
# F1 MNL (plan §5): first-order inference of the multinomial logit on D2's
# MNL population (lib/dgp_mode.R), N in {500, 2000, 8000} travellers.
# Method `ml`: run_mnlogit(), polished by Newton; SEs by the Hessian, BHHH and
# robust routes; derived estimand the value of travel time, wtp() of travel
# in vcost units (delta SE), its truth wtp() at theta0; and the LR test of
# IIA, an NL with the ground nest {bus, car, train} fitted to the same data
# (lambda = 1 is interior: df = 1), from the MNL's estimate at lambda = 1 and
# polished, so that LR >= 0 measures the statistic, not the optimizer.
# Gates at 8k (C1): bias, Hessian SE ratio and tails of every coefficient,
# the joint Wald size on the slopes, the LR size, and the VTTS interval's
# two-sided coverage, not its tails: b_vcost is weakly determined
# (|b_vcost| / SE = 2.0, 4.0, 7.9 at 500, 2k, 8k), so the ratio is skewed
# and its one-sided errors are O(N^-1/2) (C5): under the normal
# approximation of the coefficients the tails are 5.1% / 0.1% at 8k
# (coverage 94.8%), 10.9% / 0.0% at 500 (checks/wp3_mode.R).
# Findings (C5): N bias, sqrt(N) ESE, tail imbalance, BHHH and robust SE
# ratios, mean LR / q, the VTTS tails.
# Sourced by load_family() into a child of the library's environment.
# ============================================================================

f1_simulate <- function(cell, seed) {
  get_dgp("D2")$simulate(N = cell$N, seed = seed, variant = "default")
}

# The value of travel time (vcost units per minute), wtp()'s delta SE; truth
# wtp() at theta0.
f1_estimands <- function(fit, sim) {
  w <- choicer::wtp(fit, price_var = "vcost", attr_vars = "travel")
  w0 <- choicer::wtp(at_truth(fit, sim$truth), price_var = "vcost", attr_vars = "travel")
  data.frame(name = "vtts", est = w$Estimate, se = w$Std_Error, truth = w0$Estimate,
             stringsAsFactors = FALSE)
}

# Both fits are polished by Newton and must converge (require_converged(),
# lib/fit.R). The unrestricted fit of the LR test starts at the MNL's
# estimate with Lambda_1 = 1, in the "bhhh" coordinates (D2's covariates
# span units to hundreds: from its default start, scale_vars = "none" can
# fail); its lambda and decrement go to diag$lr_nl.
f1_ml <- function(sim, cell, seed) {
  fit <- require_converged(choicer::run_mnlogit(sim$data, "id", "alt", "choice", D2_X,
                                                polish = "newton"))
  V <- fit_vcovs(fit, c("hessian", "bhhh", "robust"))
  cf <- coef(fit)
  slopes <- names(cf)[fit$param_map$beta]
  asc <- names(cf)[fit$param_map$asc]
  nl <- require_converged(choicer::run_nestlogit(
    sim$data, "id", "alt", "choice", D2_X, nest_col = "nest",
    theta_init = c(cf[slopes], Lambda_1 = 1, cf[asc]), scale_vars = "bhhh",
    polish = "newton"), "the NL of the LR test")
  diag <- fit_diag(fit)
  diag$lr_nl <- c(lambda = coef(nl)[["Lambda_1"]], decrement = fit_diag(nl)$decrement)
  list(est = est_table(fit, sim$truth, V = V),
       vcov_c = vcov_common(fit, V = V$hessian),
       joint_wald = joint_wald(cf[slopes], sim$truth, V$hessian),
       lr = lr_test(nl$loglik, fit$loglik, 1),
       derived = f1_estimands(fit, sim),
       diag = diag)
}

family <- list(
  id = "F1",
  claims = c("C1", "C5"),
  dgp = "D2",
  cells = data.frame(cell = c("500", "2k", "8k"), N = c(500L, 2000L, 8000L),
                     stringsAsFactors = FALSE),
  methods = list(ml = f1_ml),
  methods_by_cell = NULL,
  R = list(quick = c("500" = 100), mc = c("500" = 5000, "2k" = 5000, "8k" = 5000)),
  simulate = f1_simulate,
  estimands = f1_estimands,
  gates = list(cell = "8k",
               metrics = list(ml = c("bias_std", "se_ratio:hessian", "tail_lower:hessian",
                                     "tail_upper:hessian", "joint_wald_size",
                                     "coverage:hessian@vtts", "lr_size"))),
  # The quick tier at 500 (S1), on the family's seeds, 2026-11-05: R = 300 and, the cell being
  # cheap, R = 3000 (below; MCSE 0.018 for bias, 0.013 for an SE ratio, 0.3 pp for a tail); reps
  # 1-100, Q3's own, pass every row. Kept: bias of every coefficient but wait (|value| at most
  # 0.083), the Hessian SE ratio (0.980-1.000) and tails (1.9-3.0%) of every coefficient, the joint
  # Wald size (4.3%) and the LR size (5.5%). Dropped: wait's bias, -0.165 (R = 300: -0.133 +- 0.058):
  # C5's O(N^-1/2) is real for the strongest coefficient at N = 500, and it alone fails 2% of R = 100
  # runs; the VTTS, its ratio skewed by a weak b_vcost (|b| / SE = 2.0; C5): lower tail 11.2% +- 0.6,
  # coverage 88.8%.
  # False failures at R = 100, 1,000 bootstrap samples of the 3000: 0.6%, recentered 0.6%.
  gates_quick = list(cell = "500",
                     metrics = list(ml = c("bias_std@travel", "bias_std@vcost",
                                           "bias_std@income_air", "bias_std@asc",
                                           "se_ratio:hessian", "tail_lower:hessian",
                                           "tail_upper:hessian", "joint_wald_size", "lr_size")))
)
