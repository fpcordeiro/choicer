# ============================================================================
# The quick tier's own checks (plan §5), beside the references of
# lib/references.R, lib/references_quad.R and lib/dgp_bayes.R:
#   Q0  every DGP against choicer at theta0: the probabilities its choices are
#       drawn from against predict() (MNL, NL; max abs diff < 1e-12), and the
#       suite's simulated log-likelihood against choicer's on a shared draws
#       cube, every parameter fixed by lower = upper (MXL; relative < 1e-12);
#   Q2  one large-N fit per family (C1): max |z| < 4, the BHHH and robust SEs
#       within 10% of the Hessian's, the polish certified.
# Each returns check rows data.frame(test, check, value, threshold, pass). A
# comparison is .ref_row() (lib/references.R), whose pass is value < threshold
# and FALSE for a NaN; q_row() builds the rest: a verdict without a threshold,
# and, with pass NA and no threshold, an information row. Definitions only.
# ============================================================================

# A row of the quick tier with its verdict given; pass NA and no threshold mark
# an information row. quick.R fails a row with a threshold and pass NA.
q_row <- function(test, check, value = NA_real_, threshold = NA_real_, pass = NA) {
  data.frame(test = test, check = check, value = value, threshold = threshold, pass = pass,
             stringsAsFactors = FALSE)
}

# What a run of quick.R is: "tier" for the whole quick tier as the plan sets it,
# else the options that shortened or changed it. Stamped on every row of its
# CSVs and on its last line, so that a shortened run never reads as the tier.
q_run_kind <- function(only = NULL, reps = 0L, perturb = 0, fixtures = FALSE) {
  k <- c(if (!is.null(only)) paste0("only=", only), if (reps > 0L) paste0("reps=", reps),
         if (perturb != 0) paste0("perturb=", perturb), if (fixtures) "fixtures")
  if (length(k)) paste(k, collapse = " ") else "tier"
}

# --- Q0: the DGPs against choicer --------------------------------------------

Q0_TOL <- c(prob = 1e-12, loglik = 1e-12)

.q0_quiet <- function(expr) suppressMessages(suppressWarnings(expr))

# MNL or NL: the probabilities a dataset's choices were drawn from (p per row,
# p0 per situation with an outside option) against predict() at theta. The
# fit, on the same data, only carries the prepared data and names: at_truth()
# puts theta in place, the NL's lambda included. `perturb` is added on
# choicer's side only (the failure injection of checks/wp9_quick.R).
.q0_predict <- function(label, data, theta, fit_fun, p, p0 = NULL, perturb = 0) {
  fit <- .q0_quiet(fit_fun(data))
  named <- identical(names(coef(fit)), names(theta))
  pr <- stats::predict(at_truth(fit, theta + perturb), type = "probabilities")
  diff <- max(abs(pr$choice_prob[, 1L] - p),
              if (!is.null(p0)) abs(pr$choice_prob_outside[, 1L] - p0))
  rbind(.ref_row("Q0", paste(label, "names(coef(fit)) = names(theta0)"), as.numeric(!named), 0.5),
        .ref_row("Q0", sprintf("%s (%d rows): max |p - predict()| at theta0", label, nrow(data)),
                 diff, Q0_TOL[["prob"]]))
}

# A spec's mixed logit in choicer's arguments; `person` makes it a panel.
.q0_mxl_model <- function(spec, person = NULL) {
  list(x = spec$beta, w = spec$random, dist = spec$rc_dist, mean = isTRUE(spec$rc_mean),
       corr = isTRUE(spec$correlated), outside = isTRUE(spec$outside), person = person)
}

# choicer's log-likelihood at theta on the cube eta, every parameter fixed.
.q0_fixed_fit <- function(data, m, theta, eta) {
  idat <- choicer::prepare_mxl_data(data, "id", "alt", "choice", m$x, m$w,
                                    include_outside_option = m$outside, rc_correlation = m$corr,
                                    person_col = m$person)
  .q0_quiet(choicer::run_mxlogit(input_data = idat, eta_draws = eta, rc_dist = m$dist,
                                 rc_mean = m$mean, rc_correlation = m$corr,
                                 include_outside_option = m$outside, theta_init = theta,
                                 lower = theta, upper = theta))
}

.q0_rel_row <- function(label, ll, cl) {
  .ref_row("Q0", label, abs(ll - cl) / abs(ll), Q0_TOL[["loglik"]])
}

# MXL: the suite's simulated log-likelihood on a shared normal cube (S = 50,
# one block per likelihood unit) against choicer's at theta.
.q0_cube <- function(label, data, spec, theta, seed, person = NULL, perturb = 0) {
  d <- design(data, spec, person = person)
  U <- if (is.null(person)) max(d$sit) else max(d$unit)
  K <- length(spec$random)
  eta <- vs_with_seed(seed, array(stats::rnorm(K * 50L * U), c(K, 50L, U)))
  ll <- mxl_loglik(d, unpack_theta("mxl", spec, theta), eta, chosen_rows(data$choice, d$sit))
  fit <- .q0_fixed_fit(data, .q0_mxl_model(spec, person), theta + perturb, eta)
  rbind(.ref_row("Q0", paste(label, "names(coef(fit)) = names(theta0)"),
                 as.numeric(!identical(names(coef(fit)), names(theta))), 0.5),
        .q0_rel_row(sprintf("%s (%d rows): simulated loglik on a shared cube, rel", label,
                            nrow(data)), ll, fit$loglik))
}

# MXL: the log-likelihood of the probabilities the choices were drawn from
# (given the DGP's own tastes) against choicer's on those tastes' normals as a
# one-draw cube: it ties the DGP's probabilities, not only the suite's engine,
# to choicer's mapping from eta to the coefficients.
.q0_own_tastes <- function(label, data, m, theta, eta, probs, sit, U) {
  lp <- list(logp = log(probs$p), logp0 = if (!is.null(probs$p0)) log(probs$p0))
  ll <- sum(.chosen_logp(lp, chosen_rows(data$choice, sit)))
  cl <- .q0_fixed_fit(data, m, theta, array(eta, c(dim(eta)[1L], 1L, U)))$loglik
  .q0_rel_row(paste(label, "the DGP's own tastes as a one-draw cube, rel"), ll, cl)
}

# Q0's rows: every DGP and variant the families use (plan §4), on small samples.
q0_checks <- function(perturb = 0) {
  rows <- list()
  add <- function(r) rows[[length(rows) + 1L]] <<- r

  # D2 mode_rp (F11's taste: its normal recovered from the DGP's coefficient).
  mnl <- function(x, outside = FALSE) {
    function(d) choicer::run_mnlogit(d, "id", "alt", "choice", x, include_outside_option = outside)
  }
  nl <- function(x, nest, outside = FALSE) {
    function(d) choicer::run_nestlogit(d, "id", "alt", "choice", x, nest_col = nest,
                                       include_outside_option = outside)
  }
  d2 <- get_dgp("D2")
  s <- d2$simulate(N = 1000L, seed = 11L, variant = "default")
  add(.q0_predict("D2 default (MNL)", s$data, s$truth, mnl(D2_X), s$meta$p, perturb = perturb))
  s <- d2$simulate(N = 1000L, seed = 11L, variant = "nl")
  add(.q0_predict("D2 nl", s$data, s$truth, nl(D2_X, "nest"), s$meta$p, perturb = perturb))
  s <- d2$simulate(N = 1000L, seed = 12L, variant = "mxl")
  sp <- d2_spec("mxl")
  add(.q0_cube("D2 mxl", s$data, sp, s$truth, 13L, perturb = perturb))
  par <- unpack_theta("mxl", sp, s$truth)
  add(.q0_own_tastes("D2 mxl:", s$data, .q0_mxl_model(sp), s$truth + perturb,
                     t((s$meta$tastes - par$mu) / par$L[1L, 1L]), list(p = s$meta$p),
                     data.table::rleid(s$data$id), 1000L))

  # D1 sp_panel: the panel and the composite likelihood; D1-L's logit, written
  # as choicer's MXL without means (its tastes' means in the fixed
  # coefficients), on its own tastes.
  d1 <- get_dgp("D1")
  s <- d1$simulate(N = 300L, seed = 21L, keep_probs = TRUE)
  sp <- s$meta$spec
  add(.q0_cube("D1 panel", s$data, sp, s$truth, 101L, person = "person", perturb = perturb))
  add(.q0_cube("D1 composite", s$data, sp, s$truth, 201L, perturb = perturb))
  add(.q0_own_tastes("D1 panel:", s$data, .q0_mxl_model(sp, "person"), s$truth + perturb,
                     s$meta$eta, s$meta$probs, data.table::rleid(s$data$id), 300L))
  s <- d1$simulate(N = 300L, seed = 31L, variant = "labeled", keep_probs = TRUE)
  th <- s$meta$theta
  names(th)[seq_along(SPL_X)] <- SPL_X
  dw <- data.table::copy(s$data)[, `:=`(price_w = price, feature_w = feature)]
  m <- list(x = SPL_X, w = c("price_w", "feature_w"), dist = c(0L, 0L), mean = FALSE, corr = TRUE,
            outside = TRUE, person = "person")
  add(.q0_own_tastes("D1-L:", dw, m, th + perturb, s$meta$eta, s$meta$probs,
                     data.table::rleid(s$data$id), 300L))

  # D4 school: the families' 40-school design (R engine, as F2 and F3), the
  # census (the port, as F7).
  for (v in c("mnl", "nl")) {
    s <- get_dgp("D4-40")$simulate(N = 1000L, seed = 41L, variant = v, keep_probs = TRUE)
    f <- if (v == "mnl") mnl(D4_X, TRUE) else nl(D4_X, "sector", TRUE)
    add(.q0_predict(paste("D4-40", v), s$data, geo_truth(s$meta, unique(s$data$alt)), f,
                    s$meta$probs$p, s$meta$probs$p0, perturb))
  }
  s <- get_dgp("D4-census")$simulate(N = 1500L, seed = 42L, variant = "mnl", keep_probs = TRUE,
                                     engine = "cpp")
  add(.q0_predict("D4-census mnl", s$data, geo_truth(s$meta, unique(s$data$alt)), mnl(D4_X, TRUE),
                  s$meta$probs$p, s$meta$probs$p0, perturb))

  # D3 hospital: the many-alternatives design and the claims panel (the port,
  # as F6 and F7), no outside option.
  for (id in c("D3", "D3-claims")) {
    s <- get_dgp(id)$simulate(N = 1500L, seed = 51L, variant = "mnl", keep_probs = TRUE,
                              engine = "cpp")
    add(.q0_predict(paste(id, "mnl"), s$data, geo_truth(s$meta, unique(s$data$alt)), mnl(D3_X$mnl),
                    s$meta$probs$p, perturb = perturb))
  }

  # MXL on the geographic designs: a shared cube and the DGP's own tastes.
  for (id in c("D4-40", "D4-census", "D3", "D3-claims")) {
    cpp <- id != "D4-40"
    s <- get_dgp(id)$simulate(N = 800L, seed = 61L, variant = "default", keep_probs = TRUE,
                              engine = if (cpp) "cpp" else "r")
    person <- if (isTRUE(s$meta$design$panel)) "patient"
    th <- geo_truth(s$meta, unique(s$data$alt))
    sp <- s$meta$spec
    sp$alts <- sp$alts[sp$alts %in% unique(s$data$alt)]
    add(.q0_cube(paste(id, "mxl"), s$data, sp, th, 71L, person = person, perturb = perturb))
    add(.q0_own_tastes(paste(id, "mxl:"), s$data, .q0_mxl_model(sp, person), th + perturb,
                       s$meta$eta, s$meta$probs, data.table::rleid(s$data$id), s$meta$N))
  }
  do.call(rbind, rows)
}

# --- Q2: one large-N fit per family --------------------------------------------

Q2_Z_MAX <- 4
Q2_SE_TOL <- 0.10

# One fit per family (plan §5, no "blp"), through the family's own simulate()
# and fit functions (fe, the environment of the family's file), at the family's
# top rung (`base`) unless `N` says otherwise; each under about 3 minutes on a
# worker (2026-11-05, one busy worker: F3 150 s, F6 133 s, F4 59 s, F2 19 s, F1
# and F6's MNL 6 and 4 s). F1's and F2's are cheap, so their N is far above the
# top rung: the larger N, the smaller the gap from theta0 a z can see. The mixed
# logits stay at the rungs where C4 and C6 leave C1 standing: F3 16k (sqrt(N)/S =
# 0.24, J/sqrt(N) = 0.32), F4 3.2k (sqrt(U)/S = 0.11); F6 20k is the most the
# budget holds (N/J = 100, J/sqrt(N) = 1.41, sqrt(N)/S = 0.71: C6 and C4 predict
# a bias of a fraction of an SD, which max |z| < 4 tolerates), and F6's MNL, at
# seconds, runs at the top rung beside it. Measured there: max |z| 2.6, the
# routes' SEs within 5.6% for the common parameters.
Q2_UNITS <- list(
  F1 = list(family = "F1", base = "8k", N = 1000000L, cost = 1, fit = function(fe, sim, cell, seed) {
    require_converged(choicer::run_mnlogit(sim$data, "id", "alt", "choice", D2_X, polish = "newton"))
  }),
  F2 = list(family = "F2", base = "16k", N = 200000L, cost = 2, fit = function(fe, sim, cell, seed) {
    fe$f2_fit(sim)
  }),
  F3 = list(family = "F3", base = "16k", cost = 5, fit = function(fe, sim, cell, seed) {
    fit_two_starts(function(theta_init) {
      fe$f3_run(sim, cell, seed, polish = "newton", theta_init = theta_init)
    }, sim$truth, "the ML fit")
  }),
  F4 = list(family = "F4", base = "3.2k", cost = 3, fit = function(fe, sim, cell, seed) {
    fit_two_starts(function(theta_init) {
      fe$f4_run(sim, seed, cell$S, person_col = "person", theta_init = theta_init)
    }, sim$truth, "the panel fit")
  }),
  F6 = list(family = "F6", base = "20k", cost = 4, fit = function(fe, sim, cell, seed) {
    fit_two_starts(function(theta_init) {
      fe$f6_run(sim, cell, seed, polish = "newton", theta_init = theta_init)
    }, sim$truth, "the MXL fit")
  }),
  F6.mnl = list(family = "F6", base = "80k.mnl", cost = 1, fit = function(fe, sim, cell, seed) {
    require_converged(choicer::run_mnlogit(sim$data, "id", "alt", "choice", D3_X$mnl,
                                           polish = "newton"), "the MNL fit")
  })
)

# Each parameter's tolerance on |SE_route / SE_hessian - 1|, the plan's 10%,
# wider where the sampling error is: for the constant of an alternative chosen n
# times, B_jj is about n and H_jj its expectation, so the routes' SEs differ by
# about 1 / (2 sqrt(n)), and 2 / sqrt(n) is 4 SD at every n.
.q2_se_tol <- function(fit, data) {
  nm <- names(coef(fit))
  tol <- stats::setNames(rep(Q2_SE_TOL, length(nm)), nm)
  asc <- nm[as.integer(fit$param_map$asc)]
  count <- table(data$alt[data$choice == 1L])
  n <- as.vector(count)[match(sub("^ASC_", "", asc), names(count))]
  n[is.na(n)] <- 0
  tol[asc] <- pmax(Q2_SE_TOL, 2 / sqrt(n))
  tol
}

# Q2's rows for a fit at the truth. Under C1 z = (theta_hat - theta0) / SE is
# about N(0, 1) per parameter: P(max |z| >= 4) is at most p * 6.3e-5 over p
# parameters (0.3% at p = 50, 1.3% at p = 210), so at a fixed seed a failure
# says that the estimator or its SE moved, not that the draw was unlucky. H = B
# (C1) has the three variance routes agree within sampling error (.q2_se_tol()).
# A route that fails, or a NaN SE of any parameter, gives a NaN and fails.
.q2_rows <- function(k, fit, sim, label) {
  th <- coef(fit)
  se <- lapply(fit_vcovs(fit, c("hessian", "bhhh", "robust")), function(V) {
    .vs_se(V) %||% stats::setNames(rep(NA_real_, length(th)), names(th))
  })
  z <- abs(th - .vs_align(sim$truth, names(th), "truth")) / se$hessian
  tol <- .q2_se_tol(fit, sim$data)
  chk <- function(check) sprintf("%s %s: %s", k, label, check)
  worst <- function(x) if (anyNA(x)) names(x)[is.na(x)][1L] else names(x)[which.max(x)]
  pol <- fit$optimizer$polish
  rbind(
    q_row("Q2", chk(sprintf("polish certified (decrement over the mean weight), code %s",
                            fit$convergence)),
          pol$decrement %||% NA_real_, pass = isTRUE(pol$converged)),
    .ref_row("Q2", chk(sprintf("max |z| over %d parameters (worst %s; P(>= %g | C1) <= %.2f%%)",
                               length(th), worst(z), Q2_Z_MAX,
                               100 * length(th) * 2 * stats::pnorm(-Q2_Z_MAX))),
             max(z), Q2_Z_MAX),
    do.call(rbind, lapply(c("bhhh", "robust"), function(r) {
      dev <- abs(se[[r]] / se$hessian - 1)
      w <- worst(dev / tol)
      .ref_row("Q2", chk(sprintf(paste("max |SE_%s / SE_hessian - 1| / tol over %d parameters",
                                       "(worst %s: %.3f, tol %.3f; tol 10%%, or 2 / sqrt(n)",
                                       "for a constant chosen n times)"),
                                 r, length(th), w, dev[[w]], tol[[w]])),
               max(dev / tol), 1)
    })))
}

# One unit of Q2 on a worker: simulate, fit, check. An error is one failed row
# (quick.R's q2_task()). Warnings are counted in an information row.
q2_unit <- function(k, L) {
  u <- Q2_UNITS[[k]]
  fam <- load_family(u$family, L)
  cell <- family_cell(fam, u$base)
  if (!is.null(u$N)) cell$N <- u$N
  seed <- vs_seed("Q2", k, 1)
  w <- character()
  watched <- function(expr) {
    withCallingHandlers(expr, warning = function(cnd) {
      w <<- c(w, class(cnd)[1])
      invokeRestart("muffleWarning")
    }, message = function(cnd) invokeRestart("muffleMessage"))
  }
  now <- function() proc.time()[["elapsed"]]
  t0 <- now()
  vs_set_seed(seed)
  sim <- watched(fam$simulate(cell, seed))
  t1 <- now()
  vs_set_seed(seed)
  fit <- watched(u$fit(environment(fam$simulate), sim, cell, seed))
  t2 <- now()
  label <- sprintf("(N = %s%s)", format(cell$N, big.mark = ",", trim = TRUE),
                   if (!is.null(cell$S) && !is.na(cell$S)) sprintf(", S = %d", cell$S) else "")
  rows <- watched(.q2_rows(k, fit, sim, label))
  st <- attr(fit, "start")
  rbind(rows, q_row("Q2", sprintf(
    "%s %s: simulate %.1f s, fit %.1f s, variances %.1f s; %d rows%s; warnings: %s", k, label,
    t1 - t0, t2 - t1, now() - t2, nrow(sim$data),
    if (!is.null(st)) sprintf("; start kept: %s, dlogL %.2e", st$winner, st$dlogL) else "",
    if (length(w)) paste(names(table(w)), table(w), sep = " x", collapse = ", ") else "none"),
    now() - t0))
}
