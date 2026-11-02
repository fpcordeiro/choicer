# ============================================================================
# F8 choice-based sampling (plan §5, C2). Per rep, a frame of M = 3 N /
# min_j Q_j i.i.d. travellers from one of D2's populations (lib/dgp_mode.R),
# then a sample stratified by the chosen mode with design shares H uniform
# (1/4 each) against the population's shares Q, N = 2000:
#   quota      the MNL population; fixed strata, 500 each:
#              sample_by_choice(n_per_alt = 500);
#   random     the MNL population; Bernoulli sampling, each frame member kept
#              with probability f_j = N H_j / (M Q_j) of its chosen mode j, so
#              n_j is random with E n_j = N / 4. Not sample_by_choice(
#              frac_per_alt =), which keeps round(f_j x available) per
#              stratum: fixed sizes given the frame, nearly quotas;
#   random_nl  the NL population (lambda_ground = .545), sampled as `random`.
# WESML weights are wesml_weights(Q = the population's exact shares, H = the
# design's 1/4): known weights Q_j / H_j. Under random sizes the sample is
# i.i.d. over the frame's members (kept or not), so the plain robust meat is
# exact; under quotas the WESML variance subtracts the stratum mean scores
# (Manski & Lerman 1977) and the plain meat is conservative (C2). Weights
# from the realized shares would give the quota estimator under random sizes
# too. The population predicts (checks/wp3_mode.R) the plain meat's excess
# under quotas in the constants only: robust SE / SD 1.008-1.018, slopes
# 1.000 (a finding); and the unweighted MNL's Hessian SE / SD 1.011-1.024 for
# the constants, 1.000 for the slopes.
# Methods and what each estimates:
#   MNL population (quota, random): mnl_unweighted (weights and provenance
#   stripped), mnl_wesml and nl_wesml (se_method "sandwich"), nl_unweighted.
#   The NL with the ground nest is correct at lambda = 1, so the WESML fits
#   estimate theta0 (Lambda_1 = 1), and the unweighted fits the sampled model
#   exactly: the MNL with constants shifted by
#     log(H_j / Q_j) - log(H_air / Q_air)    (air is choicer's reference),
#   slopes and lambda unchanged (Manski & Lerman 1977; Prentice & Pyke 1979).
#   NL population (random_nl): nl_wesml estimates the NL's theta0;
#   nl_unweighted and mnl_misspec (the unweighted MNL, misspecified here)
#   estimate their pseudo-true values on the sample's (type, mode)
#   distribution, pi_i P_ij f_j (exact, F8_NL_PSEUDO_TRUE). The sampled model
#   is no NL with shifted constants (f differs within the ground nest), so
#   the unweighted NL misses the Manski-Lerman shift: its pseudo-true value
#   differs from theta0_nl with the constants so shifted by -2.19, 0.29 and
#   -1.71 robust SE at N = 2000 (bus, car, train), -0.73 (wait), 0.71
#   (income_air), 0.52 (Lambda_1). The records report against the
#   pseudo-true values; the derived estimand `vtts` (-b_travel / b_vcost)
#   reports every method against the population's own VTTS (-0.3033 for
#   the NL's), where the MNL's misspecification shows: its pseudo-true VTTS
#   is -0.1905, 4.2 SE at N = 2000 (the unweighted NL's -0.2947, 0.15 SE).
# Gates: under random sizes (MNL population) mnl_unweighted's bias, Hessian
# SE ratio and tails of every coefficient and its Hessian joint Wald size on
# the slopes, and the WESML robust tails; nl_wesml's robust tails in
# random_nl too (each against its own truth). Everything else is a finding.
# Under quotas the unweighted constants' Hessian variance is not the
# sampling variance (Prentice & Pyke: the slopes' is), so est$gated is the
# slopes only there (quick.R's Q3 gates `quota` with mnl_unweighted). The
# misspecified MNL is its own method, so that mnl_unweighted's gates do not
# reach the NL population.
# Sourced by load_family() into a child of the library's environment.
# ============================================================================

F8_N <- 2000L
F8_H <- c(air = 0.25, bus = 0.25, car = 0.25, train = 0.25)
F8_FRAME <- as.integer(ceiling(3 * F8_N / min(D2_Q)))

# The unweighted fits' pseudo-true values on a Bernoulli sample (H = F8_H)
# from the NL population: f8_pseudo_true() on 2026-11-02, choicer
# 0.2.1.9000, one thread (decrements below 1e-13).
F8_NL_PSEUDO_TRUE <- list(
  nl_unweighted = c(wait = -0.090335283063099014, travel = -0.0029670623899442838,
                    vcost = -0.010069397947579318, income_air = 0.005387679283693751,
                    Lambda_1 = 0.56628613067110789, ASC_bus = -0.46536790721587606,
                    ASC_car = -4.1215873145604913, ASC_train = -0.48147943110085967),
  mnl_misspec = c(wait = -0.1204312769115345, travel = -0.0038962222051065706,
                  vcost = -0.020456960014024621, income_air = 0.00066270885766312566,
                  ASC_bus = -1.7912829609964853, ASC_car = -6.7118307659514587,
                  ASC_train = -1.7162545999292558))

# Reproduces F8_NL_PSEUDO_TRUE: the NL and the MNL fitted, unweighted, to the
# sample's (type, mode) distribution (d2_cb_weights()).
f8_pseudo_true <- function() {
  W <- d2_cb_weights("nl", unname(F8_H[D2_ALTS]))
  list(nl_unweighted = coef(d2_pseudo_fit(W, "nl")), mnl_misspec = coef(d2_pseudo_fit(W)))
}

# Each method's truth under the cell's population, and the population's VTTS.
f8_truths <- function(population) {
  if (population == "nl") {
    th <- D2_THETA0$nl
    return(c(list(nl_wesml = th, vtts = -th[["travel"]] / th[["vcost"]]), F8_NL_PSEUDO_TRUE))
  }
  th <- D2_THETA0$default
  Q <- d2_pop_shares("default")
  asc <- paste0("ASC_", D2_ALTS[-1L])
  shifted <- th
  shifted[asc] <- th[asc] + unname(log(F8_H[D2_ALTS[-1L]] / Q[-1L]) - log(F8_H[["air"]] / Q[1L]))
  nl <- function(t) c(t[D2_X], Lambda_1 = 1, t[asc])
  list(mnl_wesml = th, mnl_unweighted = shifted, nl_wesml = nl(th), nl_unweighted = nl(shifted),
       vtts = -th[["travel"]] / th[["vcost"]])
}

# The frame, then the cell's sample with WESML weights (data, with
# provenance) and the same rows without them (meta$unweighted). Two seeds
# drawn from `seed`: the frame's and the sampling's.
f8_simulate <- function(cell, seed) {
  s <- vs_with_seed(seed, sample.int(2147483646L, 2L))
  frame <- get_dgp("D2")$simulate(N = cell$frame, seed = s[1], variant = cell$population)$data
  Q <- d2_pop_shares(cell$population)
  smp <- vs_with_seed(s[2], if (cell$sampling == "quota") {
    out <- choicer::sample_by_choice(frame, "id", "alt", "choice", n_per_alt = cell$n_per_alt)
    out[, .wesml_weight := NULL]
    out
  } else {
    y <- frame$alt[frame$choice == 1L]                      # frame members in id order
    f <- F8_N * F8_H[y] / (cell$frame * Q[y])
    keep <- which(stats::runif(cell$frame) < f)
    frame[id %in% keep]
  })
  for (a in c("Q", "H", "choice_sampling")) data.table::setattr(smp, a, NULL)
  data.table::setorder(smp, id, alt)
  unweighted <- data.table::copy(smp)
  data <- choicer::wesml_weights(smp, "id", "alt", "choice", Q = Q, H = F8_H, attach = TRUE)
  n_j <- table(factor(smp$alt[smp$choice == 1L], levels = D2_ALTS))
  list(data = data, truth = D2_THETA0[[cell$population]],
       meta = list(unweighted = unweighted, truths = f8_truths(cell$population), n_j = n_j,
                   N = sum(n_j), frame = cell$frame, design = cell$cell))
}

# The value of travel time, wtp()'s delta SE (the fit's own variance
# route); truth the population's VTTS.
f8_estimands <- function(fit, sim) {
  w <- choicer::wtp(fit, price_var = "vcost", attr_vars = "travel")
  data.frame(name = "vtts", est = w$Estimate, se = w$Std_Error, truth = sim$meta$truths$vtts,
             stringsAsFactors = FALSE)
}

# One method's fit and record. `weighted` fits the WESML weights with the
# sandwich; `nested` fits the NL with the ground nest in the "bhhh"
# coordinates (with run_nestlogit()'s default "none", 18 and 15 of 30 WESML
# fits (quota, random) stopped at the start with NLOPT_FAILURE; with
# "bhhh", none of 120 NL fits did, and the log-likelihoods agreed to 2e-9
# where both converged), multi-start on every 10th rep. Every fit is
# polished by Newton and must converge (require_converged(), lib/fit.R).
# est$gated: every coefficient under random sizes from the MNL population
# and nl_wesml's in random_nl; the slopes of the unweighted MNL under
# quotas; nothing else.
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
  gated <- switch(cell$cell, random = TRUE, random_nl = method == "nl_wesml",
                  quota = if (method == "mnl_unweighted") "beta" else FALSE)
  diag <- fit_diag(fit)
  if (nested && vs_seed_rep("F8", cell$cell, seed) %% 10 == 0) {
    diag$multistart <- multistart(fit, refit = call_fit, seed = seed)
  }
  list(est = est_table(fit, truth, V = V, gated = gated),
       vcov_c = vcov_common(fit, route = if (weighted) "robust" else "hessian",
                            V = if (weighted) V$robust else V$hessian),
       joint_wald = list(hessian = joint_wald(cf[common], truth, V$hessian),
                         robust = joint_wald(cf[common], truth, V$robust)),
       derived = f8_estimands(fit, sim),
       diag = diag)
}

# A method of the family: f8_fit() with its choices fixed.
f8_method <- function(method, weighted, nested) {
  force(method)
  force(weighted)
  force(nested)
  function(sim, cell, seed) f8_fit(sim, cell, seed, method, weighted, nested)
}

F8_MNL_POP <- c("mnl_unweighted", "mnl_wesml", "nl_wesml", "nl_unweighted")

family <- list(
  id = "F8",
  claims = "C2",
  dgp = "D2",
  cells = data.frame(cell = c("quota", "random", "random_nl"), N = F8_N, frame = F8_FRAME,
                     population = c("default", "default", "nl"),
                     sampling = c("quota", "bernoulli", "bernoulli"),
                     n_per_alt = c(F8_N %/% 4L, NA, NA), stringsAsFactors = FALSE),
  methods = list(mnl_unweighted = f8_method("mnl_unweighted", FALSE, FALSE),
                 mnl_wesml = f8_method("mnl_wesml", TRUE, FALSE),
                 nl_wesml = f8_method("nl_wesml", TRUE, TRUE),
                 nl_unweighted = f8_method("nl_unweighted", FALSE, TRUE),
                 mnl_misspec = f8_method("mnl_misspec", FALSE, FALSE)),
  methods_by_cell = list(quota = F8_MNL_POP, random = F8_MNL_POP,
                         random_nl = c("nl_wesml", "nl_unweighted", "mnl_misspec")),
  R = list(quick = c(quota = 100), mc = c(quota = 2000, random = 2000, random_nl = 2000)),
  simulate = f8_simulate,
  estimands = f8_estimands,
  gates = list(cell = c("random", "random_nl"),
               metrics = list(mnl_unweighted = c("bias_std", "se_ratio:hessian",
                                                 "tail_lower:hessian", "tail_upper:hessian",
                                                 "joint_wald_size:hessian"),
                              mnl_wesml = c("tail_lower:robust", "tail_upper:robust"),
                              nl_wesml = c("tail_lower:robust", "tail_upper:robust")))
)
