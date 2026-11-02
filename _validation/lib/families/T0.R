# ============================================================================
# T0: a plumbing fixture, not a validation family. A multinomial logit with
# three alternatives, two covariates and constants (alternative 1 the
# reference), no outside option, simulated inline in base R, for the
# runner's and the report's tests. Method `ml` fits it; method `fail` stops
# on even reps and warns on odd ones (then fits as `ml`), so that failure
# and warning counts can be tested. Sourced by load_family() into a child of
# the library's environment.
# ============================================================================

t0_theta0 <- c(x1 = 0.8, x2 = -1, ASC_2 = 0.5, ASC_3 = -0.5)

# N situations x 3 alternatives; utilities x1 * 0.8 - x2 + delta_j plus
# Gumbel errors, the chosen alternative the utility's maximum.
t0_simulate <- function(cell, seed) {
  vs_set_seed(seed)
  N <- cell$N
  J <- 3L
  x1 <- stats::rnorm(N * J)
  x2 <- stats::rnorm(N * J)
  v <- t0_theta0[["x1"]] * x1 + t0_theta0[["x2"]] * x2 +
    rep(c(0, t0_theta0[["ASC_2"]], t0_theta0[["ASC_3"]]), N)
  u <- matrix(v - log(-log(stats::runif(N * J))), N, J, byrow = TRUE)
  chosen <- max.col(u, ties.method = "first")
  data <- data.frame(id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N), x1 = x1, x2 = x2,
                     choice = as.integer(rep(seq_len(J), N) == rep(chosen, each = J)))
  list(data = data, truth = t0_theta0, meta = list(N = N, J = J))
}

# Derived estimand: the WTP for x1 in units of x2 (x2 is the price), with
# wtp()'s delta-method SE; its truth is wtp() at theta0.
t0_estimands <- function(fit, sim) {
  w <- choicer::wtp(fit, price_var = "x2", attr_vars = "x1")
  w0 <- choicer::wtp(at_truth(fit, sim$truth), price_var = "x2", attr_vars = "x1")
  data.frame(name = "wtp_x1", est = w$Estimate, se = w$Std_Error, truth = w0$Estimate,
             stringsAsFactors = FALSE)
}

t0_ml <- function(sim, cell, seed) {
  fit <- choicer::run_mnlogit(sim$data, "id", "alt", "choice", c("x1", "x2"))
  V <- fit_vcovs(fit, c("hessian", "bhhh", "robust"))
  b <- fit$param_map$beta
  list(est = est_table(fit, sim$truth, V = V, gated = "beta"),
       vcov_c = vcov_common(fit, V = V$hessian),
       joint_wald = joint_wald(coef(fit)[b], sim$truth, V$hessian),
       derived = t0_estimands(fit, sim),
       diag = fit_diag(fit))
}

# The rep is recovered from the seed (vs_seed() offsets it by the rep).
t0_fail <- function(sim, cell, seed) {
  rep <- vs_seed_rep("T0", cell$cell, seed)
  if (rep %% 2 == 0) stop("T0 fixture: deliberate failure on even rep ", rep)
  warning(warningCondition(paste("T0 fixture: deliberate warning on odd rep", rep),
                           class = "vs_fixture_warning"))
  t0_ml(sim, cell, seed)
}

family <- list(
  id = "T0",
  claims = "C1",
  dgp = "inline",
  cells = data.frame(cell = c("200", "800"), N = c(200L, 800L), stringsAsFactors = FALSE),
  methods = list(ml = t0_ml, fail = t0_fail),
  methods_by_cell = NULL,
  R = list(quick = c("200" = 4, "800" = 4), mc = c("200" = 4, "800" = 4)),
  simulate = t0_simulate,
  estimands = t0_estimands,
  gates = list(cell = "800",
               metrics = list(ml = c("bias_std", "se_ratio:hessian", "tail_lower:hessian",
                                     "tail_upper:hessian")))
)
