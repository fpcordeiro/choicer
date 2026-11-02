# ============================================================================
# The quick tier's cross-package references on the shipped mode_choice
# (plan §5): Q1.1, choicer's MNL against survival::clogit and mlogit; Q1.2,
# choicer's NL against mlogit's nested logit. Each returns check rows
# data.frame(test, check, value, threshold, pass) (quick.R's format); the
# check text carries each optimizer's convergence measure, so that a row
# shows the comparison measures the estimators, not the stopping rules.
# The data and specification are D2's (lib/dgp_mode.R: d2_mode_data(),
# D2_X): wait, travel, vcost, income x air and constants, air the
# reference, unweighted. The quadrature reference (Q1.3) is
# lib/references_quad.R. Definitions only.
# ============================================================================

# One check row; pass is value < threshold.
.ref_row <- function(test, check, value, threshold) {
  data.frame(test = test, check = check, value = value, threshold = threshold,
             pass = isTRUE(value < threshold), stringsAsFactors = FALSE)
}

# mode_choice for the references: D2's columns plus 0/1 dummies of bus, car
# and train (clogit's constants).
.ref_mode_data <- function() {
  dt <- d2_mode_data()
  for (a in D2_ALTS[-1L]) data.table::set(dt, j = a, value = as.numeric(dt$alt == a))
  dt
}

# mlogit's long data (dfidx), the alternative index named `alt`.
.ref_mlogit_data <- function(dt) {
  dfidx::dfidx(as.data.frame(dt[, c("id", "alt", "choice", D2_X), with = FALSE]),
               idx = c("id", "alt"))
}

# Coefficients and SEs of a reference renamed to choicer's names: `map` is
# reference name -> choicer name.
.ref_align <- function(est, se, map, nm) {
  est <- est[names(map)]
  se <- se[names(map)]
  names(est) <- names(se) <- unname(map)
  list(est = est[nm], se = se[nm])
}

# Q1.1: choicer's MNL (Newton polish) against clogit (exact partial
# likelihood, strata = travellers; one choice per stratum, so it is the MNL
# likelihood) and mlogit (Newton-Raphson), both at tightened tolerances:
# max |dtheta| / SE < 1e-4, max relative SE difference < 1e-4 (Hessian SEs
# on every side), |dlogL| < 1e-6.
ref_q11 <- function() {
  vs_require(c("survival", "mlogit", "dfidx"))
  dt <- .ref_mode_data()
  fit <- suppressMessages(choicer::run_mnlogit(dt, "id", "alt", "choice", D2_X,
                                               polish = "newton"))
  nm <- names(coef(fit))
  se <- fit$se
  ll <- fit$loglik
  dec <- fit$optimizer$polish$decrement
  asc <- stats::setNames(paste0("ASC_", D2_ALTS[-1L]), D2_ALTS[-1L])

  # clogit() evaluates its formula's strata() and its call to coxph() in the
  # caller's frame, so both are made to resolve in survival's namespace
  # (survival need not be attached).
  form <- stats::as.formula(paste("choice ~", paste(c(D2_X, D2_ALTS[-1L]), collapse = " + "),
                                  "+ strata(id)"), env = new.env(parent = asNamespace("survival")))
  clogit_ns <- function(...) clogit(...)
  environment(clogit_ns) <- asNamespace("survival")
  cl <- clogit_ns(form, data = as.data.frame(dt),
                  control = survival::coxph.control(eps = 1e-11, iter.max = 100))
  a_cl <- .ref_align(stats::coef(cl), sqrt(diag(stats::vcov(cl))),
                     c(stats::setNames(D2_X, D2_X), asc), nm)

  ml <- mlogit::mlogit(stats::as.formula(paste("choice ~", paste(D2_X, collapse = " + "))),
                       data = .ref_mlogit_data(dt), reflevel = "air", method = "nr",
                       tol = 1e-14, ftol = 1e-14, iterlim = 200)
  a_ml <- .ref_align(stats::coef(ml), sqrt(diag(stats::vcov(ml))),
                     c(stats::setNames(D2_X, D2_X),
                       stats::setNames(asc, paste0("(Intercept):", D2_ALTS[-1L]))), nm)

  conv <- sprintf("choicer decrement %.1e; clogit %d iterations; mlogit code %s",
                  dec, cl$iter, format(ml$est.stat$code))
  refs <- list(clogit = list(a = a_cl, ll = cl$loglik[2]),
               mlogit = list(a = a_ml, ll = as.numeric(stats::logLik(ml))))
  rows <- lapply(names(refs), function(k) {
    a <- refs[[k]]$a
    rbind(.ref_row("Q1.1", sprintf("MNL vs %s: max |dtheta|/SE (%s)", k, conv),
                   max(abs(coef(fit) - a$est) / se), 1e-4),
          .ref_row("Q1.1", sprintf("MNL vs %s: max |dSE|/SE", k), max(abs(se - a$se) / a$se),
                   1e-4),
          .ref_row("Q1.1", sprintf("MNL vs %s: |dlogL| (logL %.6f)", k, ll),
                   abs(ll - refs[[k]]$ll), 1e-6))
  })
  do.call(rbind, rows)
}

# Q1.2: choicer's NL with the ground nest {bus, car, train} (air its own
# nest) against mlogit's nested logit with un.nest.el = TRUE: one
# dissimilarity `iv` for both nests, which in RU2 (mlogit's default,
# unscaled = FALSE) leaves the singleton fly nest unaffected, so iv is
# choicer's Lambda_1. max |dtheta| / SE < 1e-3 (choicer's Hessian SEs),
# |dlogL| < 1e-5.
ref_q12 <- function() {
  vs_require(c("mlogit", "dfidx"))
  dt <- .ref_mode_data()
  fit <- suppressMessages(choicer::run_nestlogit(dt, "id", "alt", "choice", D2_X,
                                                 nest_col = "nest", polish = "newton"))
  nm <- names(coef(fit))
  asc <- stats::setNames(paste0("ASC_", D2_ALTS[-1L]), paste0("(Intercept):", D2_ALTS[-1L]))
  nests <- list(fly = "air", ground = c("bus", "car", "train"))
  ml <- mlogit::mlogit(stats::as.formula(paste("choice ~", paste(D2_X, collapse = " + "))),
                       data = .ref_mlogit_data(dt), reflevel = "air", nests = nests,
                       un.nest.el = TRUE, method = "bfgs", tol = 1e-14, ftol = 1e-14,
                       iterlim = 5000)
  a <- .ref_align(stats::coef(ml), sqrt(diag(stats::vcov(ml))),
                  c(stats::setNames(D2_X, D2_X), iv = "Lambda_1", asc), nm)
  ll_ml <- as.numeric(stats::logLik(ml))
  conv <- sprintf("choicer decrement %.1e; mlogit code %s, %d iterations",
                  fit$optimizer$polish$decrement, format(ml$est.stat$code),
                  as.integer(ml$est.stat$nb.iter))
  rbind(
    .ref_row("Q1.2", sprintf("NL vs mlogit (un.nest.el): max |dtheta|/SE (%s)", conv),
             max(abs(coef(fit) - a$est) / fit$se), 1e-3),
    .ref_row("Q1.2", sprintf("NL vs mlogit: |dlogL| (logL %.6f; Lambda_1 %.6f)", fit$loglik,
                             coef(fit)[["Lambda_1"]]), abs(fit$loglik - ll_ml), 1e-5))
}
