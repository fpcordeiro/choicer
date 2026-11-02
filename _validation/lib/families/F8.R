# ============================================================================
# F8 choice-based sampling (plan §5, C2). Per rep, a frame of M = 3 N /
# min_j Q_j i.i.d. travellers from D2's MNL population (lib/dgp_mode.R),
# then a sample stratified by the chosen mode with design shares H uniform
# (1/4 each) against Q, N = 2000:
#   quota   fixed strata, 500 each: sample_by_choice(n_per_alt = 500);
#   random  Bernoulli sampling, each frame member kept with probability
#           f_j = N H_j / (M Q_j) of its chosen mode j, so n_j is random with
#           E n_j = N / 4. Not sample_by_choice(frac_per_alt =), which keeps
#           round(f_j x available) per stratum: fixed sizes given the frame,
#           nearly quotas.
# WESML weights are wesml_weights(Q = the true Q, H = the design's 1/4):
# known weights Q_j / H_j. Under random sizes the sample is i.i.d. over the
# frame's members (kept or not), so the plain robust meat is exact; under
# quotas the WESML variance subtracts the stratum mean scores (Manski &
# Lerman 1977) and the plain meat is conservative (C2). Weights from the
# realized shares would give the quota estimator under random sizes too.
# Here the population predicts (checks/wp3_mode.R) the plain meat's excess
# under quotas in the constants only: robust SE / SD 1.008-1.018, slopes
# 1.000; and the unweighted MNL's Hessian SE / SD 1.011-1.024 for the
# constants, 1.000 for the slopes.
# Methods: mnl_unweighted (weights and provenance stripped), mnl_wesml and
# nl_wesml (se_method "sandwich"), nl_unweighted. One population per rep,
# the MNL's: there the NL with the ground nest is correct at lambda = 1, so
# nl_wesml estimates theta0 with Lambda_1 = 1, and the unweighted fits the
# sampled model exactly: the MNL with constants shifted by
#   log(H_j / Q_j) - log(H_air / Q_air)    (air is choicer's reference),
# slopes and lambda unchanged (Manski & Lerman 1977; Prentice & Pyke 1979).
# So nl_unweighted's predicted bias is 0 here; its bias under an NL truth
# would need the NL population.
# Gates (random sizes): mnl_unweighted bias, Hessian SE ratio and tails of
# every coefficient and the size of its Hessian joint Wald test on the
# slopes; WESML robust tails. Under quotas everything is a
# finding; there the unweighted constants' Hessian variance is not the
# sampling variance (Prentice & Pyke: the slopes' is), so est$gated is the
# slopes only (quick.R's Q3 gates `quota` with mnl_unweighted).
# Sourced by load_family() into a child of the library's environment.
# ============================================================================

F8_N <- 2000L
F8_H <- c(air = 0.25, bus = 0.25, car = 0.25, train = 0.25)
F8_FRAME <- as.integer(ceiling(3 * F8_N / min(D2_Q)))

# The truths: theta0 (WESML), theta0 with the constants shifted
# (unweighted), each with Lambda_1 = 1 for the NL.
f8_truths <- function() {
  th <- D2_THETA0$default
  asc <- paste0("ASC_", D2_ALTS[-1L])
  shift <- log(F8_H[D2_ALTS[-1L]] / D2_Q[D2_ALTS[-1L]]) - log(F8_H[["air"]] / D2_Q[["air"]])
  shifted <- th
  shifted[asc] <- th[asc] + unname(shift)
  nl <- function(t) c(t[D2_X], Lambda_1 = 1, t[asc])
  list(mnl_wesml = th, mnl_unweighted = shifted, nl_wesml = nl(th), nl_unweighted = nl(shifted))
}

# The frame, then the cell's sample with WESML weights (data, with
# provenance) and the same rows without them (meta$unweighted). Two seeds
# drawn from `seed`: the frame's and the sampling's.
f8_simulate <- function(cell, seed) {
  s <- vs_with_seed(seed, sample.int(2147483646L, 2L))
  pop <- get_dgp("D2")$simulate(N = cell$frame, seed = s[1], variant = "default")
  frame <- pop$data
  smp <- vs_with_seed(s[2], if (cell$cell == "quota") {
    out <- choicer::sample_by_choice(frame, "id", "alt", "choice", n_per_alt = cell$n_per_alt)
    out[, .wesml_weight := NULL]
    out
  } else {
    y <- frame$alt[frame$choice == 1L]                      # frame members in id order
    f <- F8_N * F8_H[y] / (cell$frame * D2_Q[y])
    keep <- which(stats::runif(cell$frame) < f)
    frame[id %in% keep]
  })
  for (a in c("Q", "H", "choice_sampling")) data.table::setattr(smp, a, NULL)
  data.table::setorder(smp, id, alt)
  unweighted <- data.table::copy(smp)
  data <- choicer::wesml_weights(smp, "id", "alt", "choice", Q = D2_Q, H = F8_H, attach = TRUE)
  n_j <- table(factor(smp$alt[smp$choice == 1L], levels = D2_ALTS))
  list(data = data, truth = D2_THETA0$default,
       meta = list(unweighted = unweighted, truths = f8_truths(), n_j = n_j,
                   N = sum(n_j), frame = cell$frame, design = cell$cell))
}

f8_estimands <- function(fit, sim) {
  data.frame(name = character(), est = numeric(), se = numeric(), truth = numeric(),
             stringsAsFactors = FALSE)
}

# One method's fit and record. `weighted` fits the WESML weights with the
# sandwich; `nested` fits the NL with the ground nest in the "bhhh"
# coordinates (with run_nestlogit()'s default "none", 18 and 15 of 30 WESML
# fits (quota, random) stopped at the start with NLOPT_FAILURE; with
# "bhhh", none of 120 NL fits did, and the log-likelihoods agreed to 2e-9
# where both converged), multi-start on every 10th rep. Every fit is
# polished by Newton and must converge (require_converged(), lib/fit.R).
# est$gated: every coefficient under random sizes; the slopes of the
# unweighted MNL under quotas; nothing else.
f8_fit <- function(sim, cell, seed, method, weighted, nested) {
  data <- if (weighted) sim$data else sim$meta$unweighted
  call_fit <- function(theta_init = NULL) require_converged(if (nested && weighted) {
    choicer::run_nestlogit(data, "id", "alt", "choice", D2_X, nest_col = "nest",
                           theta_init = theta_init, scale_vars = "bhhh", polish = "newton",
                           weights_col = ".wesml_weight", se_method = "sandwich")
  } else if (nested) {
    choicer::run_nestlogit(data, "id", "alt", "choice", D2_X, nest_col = "nest",
                           theta_init = theta_init, scale_vars = "bhhh", polish = "newton")
  } else if (weighted) {
    choicer::run_mnlogit(data, "id", "alt", "choice", D2_X, polish = "newton",
                         weights_col = ".wesml_weight", se_method = "sandwich")
  } else {
    choicer::run_mnlogit(data, "id", "alt", "choice", D2_X, polish = "newton")
  })
  fit <- call_fit()
  truth <- sim$meta$truths[[method]]
  V <- fit_vcovs(fit, c("hessian", "robust"))
  cf <- coef(fit)
  common <- names(cf)[setdiff(seq_along(cf), fit$param_map$asc)]
  gated <- if (cell$cell == "random") TRUE else if (method == "mnl_unweighted") "beta" else FALSE
  diag <- fit_diag(fit)
  if (nested && vs_seed_rep("F8", cell$cell, seed) %% 10 == 0) {
    diag$multistart <- multistart(fit, refit = call_fit, seed = seed)
  }
  list(est = est_table(fit, truth, V = V, gated = gated),
       vcov_c = vcov_common(fit, route = if (weighted) "robust" else "hessian",
                            V = if (weighted) V$robust else V$hessian),
       joint_wald = list(hessian = joint_wald(cf[common], truth, V$hessian),
                         robust = joint_wald(cf[common], truth, V$robust)),
       diag = diag)
}

# A method of the family: f8_fit() with its choices fixed.
f8_method <- function(method, weighted, nested) {
  force(method)
  force(weighted)
  force(nested)
  function(sim, cell, seed) f8_fit(sim, cell, seed, method, weighted, nested)
}

family <- list(
  id = "F8",
  claims = "C2",
  dgp = "D2",
  cells = data.frame(cell = c("quota", "random"), N = F8_N, frame = F8_FRAME,
                     n_per_alt = c(F8_N %/% 4L, NA), stringsAsFactors = FALSE),
  methods = list(mnl_unweighted = f8_method("mnl_unweighted", FALSE, FALSE),
                 mnl_wesml = f8_method("mnl_wesml", TRUE, FALSE),
                 nl_wesml = f8_method("nl_wesml", TRUE, TRUE),
                 nl_unweighted = f8_method("nl_unweighted", FALSE, TRUE)),
  methods_by_cell = NULL,
  R = list(quick = c(quota = 100), mc = c(quota = 2000, random = 2000)),
  simulate = f8_simulate,
  estimands = f8_estimands,
  gates = list(cell = "random",
               metrics = list(mnl_unweighted = c("bias_std", "se_ratio:hessian",
                                                 "tail_lower:hessian", "tail_upper:hessian",
                                                 "joint_wald_size:hessian"),
                              mnl_wesml = c("tail_lower:robust", "tail_upper:robust"),
                              nl_wesml = c("tail_lower:robust", "tail_upper:robust")))
)
