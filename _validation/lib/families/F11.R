# ============================================================================
# F11 misspecification (plan §5, C3): D2's MXL population (a normal
# coefficient on travel, coefficient of variation 1; lib/dgp_mode.R) fitted
# by the MNL. The MNL is then a quasi-likelihood: theta_hat -> theta_star,
# the pseudo-true value (D2_THETA_STAR, exact over the 210 traveller
# types), and only the sandwich H^-1 B H^-1 estimates its variance (White
# 1982). N in {500, 8000}. Method `ml`: run_mnlogit(), polished by Newton;
# SEs by the Hessian and robust routes; the joint Wald test on the slopes
# under each (named tests "hessian" and "robust"); the common block of the
# covariance by the robust route.
# Gates at 8k: bias against theta_star, the robust SE ratio and tails, and
# the size of the robust joint Wald test (only the robust one is valid, C3).
# Finding: the Hessian SE ratio differs from 1 (population values at
# theta_star, Hessian SE / robust SE: travel 0.79, income_air 0.84, the
# constants 0.85-0.88, wait 0.97, vcost 0.99).
# Sourced by load_family() into a child of the library's environment.
# ============================================================================

# D2's "mxl" data; the truth is the MNL's pseudo-true value (the MXL's
# theta0 stays in meta$theta0).
f11_simulate <- function(cell, seed) {
  sim <- get_dgp("D2")$simulate(N = cell$N, seed = seed, variant = "mxl")
  sim$meta$theta0 <- sim$truth
  sim$truth <- sim$meta$theta_star
  sim
}

f11_estimands <- function(fit, sim) {
  data.frame(name = character(), est = numeric(), se = numeric(), truth = numeric(),
             stringsAsFactors = FALSE)
}

f11_ml <- function(sim, cell, seed) {
  fit <- require_converged(choicer::run_mnlogit(sim$data, "id", "alt", "choice", D2_X,
                                                polish = "newton"))
  V <- fit_vcovs(fit, c("hessian", "robust"))
  cf <- coef(fit)
  slopes <- names(cf)[fit$param_map$beta]
  list(est = est_table(fit, sim$truth, V = V),
       vcov_c = vcov_common(fit, route = "robust", V = V$robust),
       joint_wald = list(hessian = joint_wald(cf[slopes], sim$truth, V$hessian),
                         robust = joint_wald(cf[slopes], sim$truth, V$robust)),
       diag = fit_diag(fit))
}

family <- list(
  id = "F11",
  claims = "C3",
  dgp = "D2",
  cells = data.frame(cell = c("500", "8k"), N = c(500L, 8000L), stringsAsFactors = FALSE),
  methods = list(ml = f11_ml),
  methods_by_cell = NULL,
  R = list(quick = NULL, mc = c("500" = 2000, "8k" = 2000)),
  simulate = f11_simulate,
  estimands = f11_estimands,
  gates = list(cell = "8k",
               metrics = list(ml = c("bias_std", "se_ratio:robust", "tail_lower:robust",
                                     "tail_upper:robust", "joint_wald_size:robust")))
)
