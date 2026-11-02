# ============================================================================
# F2: nested logit (plan §5; claims C1, C5, C6) on D4-small's NL variant:
# school choice among each student's 15 nearest schools of 4 markets x 40
# (160 constants), nests by sector (lambda: public .5, charter .7, private
# .85), the outside option a degenerate nest. Cells N in {1k, 4k, 16k}
# students. Method `ml`: run_nestlogit() with polish = "newton" (the
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
# F2_ASC_GATE_COUNT (incidental constants, C6: a constant's bias and
# skewness are O(1/n_j), so those of rarely chosen schools are findings);
# joint Wald size on beta and lambda; the WTP's tails.
# ============================================================================

F2_ASC_GATE_COUNT <- 50

.f2 <- new.env()

# Each school's expected share at theta0 (mean probability over 40,000
# students of one seed), computed once per process.
f2_shares <- function() {
  if (is.null(.f2$shares)) {
    .f2$shares <- get_dgp("D4")$simulate(N = 40000L, seed = 1L, variant = "nl")$meta$shares
  }
  .f2$shares
}

f2_simulate <- function(cell, seed) {
  geo_drop(get_dgp("D4")$simulate(N = cell$N, seed = seed, variant = "nl"))
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
  claims = c("C1", "C5", "C6"),
  dgp = "D4",
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
                                     "tail_upper:hessian@derived")))
)
