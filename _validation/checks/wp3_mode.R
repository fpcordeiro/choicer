# ============================================================================
# Self-check of WP3: D2 mode_rp (lib/dgp_mode.R), the families F1, F8, F11
# and the references Q1.1, Q1.2 (lib/references.R).
#   1. DGP = model: D2's choice probabilities at theta0 equal predict() at
#      theta0 (MNL, NL; max abs diff < 1e-12); MXL: the suite's simulated
#      log-likelihood on a shared eta_draws cube equals choicer's fit with
#      every parameter fixed (< 1e-10 relative).
#   2. Calibration: the population's expected shares at theta0 equal Q
#      (MNL: the WESML first-order conditions; MXL: the matched constants;
#      < 1e-8); the Gauss-Hermite rule of the MXL against a trapezoid rule;
#      theta_star refitted from the frozen theta0 equals D2_THETA_STAR; the
#      MNL and NL populations' weighted fits return their theta0; F8's
#      frozen pseudo-true values under the NL population are reproduced.
#   3. theta_star exact against a simulated N = 10^7 fit (the MNL
#      log-likelihood of 10^7 travellers is that of their (type, mode)
#      counts, so the fit is a weighted MNL on the counts, checked against
#      an unaggregated fit of the first 10^6): |difference| / SD < 4.
#   4. The exact constant shift on a choice-based sample of 10^5 (quotas of
#      25,000, unweighted MNL): constants within 3 SE of theta0's shifted by
#      log(H_j / Q_j) - log(H_air / Q_air), slopes within 3 SE of theta0's.
#   5. ref_q11() and ref_q12() pass.
#   6. Smoke: F1, F8 and F11 at every cell and method, reps 1 and 2,
#      through run_rep(): every record ok, finite estimates and SEs.
#   INFO lines (no verdict): identification strength (NL lambda, the free
#   MXL fit's sigma), population SEs and |b_vcost| / SE per F1 cell, the
#   VTTS interval's tails under the normal approximation of the
#   coefficients, F8's and F11's predicted SE ratios from the population,
#   F8 random_nl's gaps from the Manski-Lerman shift and in the VTTS.
#
# Usage (about 30 seconds; one OpenMP thread, the dev build first on R_LIBS):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 Rscript _validation/checks/wp3_mode.R
# Exits with status 1 on any FAIL.
# ============================================================================

root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                               value = TRUE)[1]))))
sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
b <- boot$vs_boot(root)
suppressPackageStartupMessages(library(data.table))
vs_threads(1L)
options(width = 160L)
t_start <- proc.time()[["elapsed"]]
cat("choicer", b$build$version, "from", b$build$path, "\n")

verdicts <- logical(0)
report <- function(check, label, value, tol, pass = isTRUE(value < tol)) {
  cat(sprintf("%s  %-3s %-66s %10.3e  (tol %.0e)\n", if (pass) "PASS" else "FAIL", check, label,
              value, tol))
  verdicts[length(verdicts) + 1L] <<- pass
}
info <- function(...) cat("      ", ..., "\n", sep = "")
quiet <- function(expr) suppressMessages(.d2_quiet(expr))
dgp <- get_dgp("D2")

# --- 1. DGP = model ------------------------------------------------------------------

for (v in c("default", "nl")) {
  sim <- dgp$simulate(N = 2000L, seed = 11L, variant = v)
  fit <- quiet(if (v == "nl") {
    choicer::run_nestlogit(sim$data, "id", "alt", "choice", D2_X, nest_col = "nest")
  } else {
    choicer::run_mnlogit(sim$data, "id", "alt", "choice", D2_X)
  })
  report("1", sprintf("%s: names(coef(fit)) = names(theta0)", v),
         as.numeric(!identical(names(coef(fit)), names(sim$truth))), 0.5)
  fit$coefficients[names(sim$truth)] <- sim$truth
  pr <- predict(fit, type = "probabilities")$choice_prob[, 1]
  report("1", sprintf("%s: max |p_DGP - predict()| at theta0", v), max(abs(pr - sim$meta$p)), 1e-12)
}
sim <- dgp$simulate(N = 1000L, seed = 12L, variant = "mxl")
sp <- d2_spec("mxl")
idat <- choicer::prepare_mxl_data(sim$data, "id", "alt", "choice", sp$beta, sp$random)
eta <- vs_with_seed(13L, array(stats::rnorm(50L * 1000L), c(1L, 50L, 1000L)))
th <- sim$truth
fit <- quiet(choicer::run_mxlogit(input_data = idat, eta_draws = eta, rc_dist = 0L, rc_mean = TRUE,
                                  theta_init = th, lower = th, upper = th))
d <- design(sim$data, sp)
ll <- mxl_loglik(d, unpack_theta("mxl", sp, th), eta, chosen_rows(sim$data$choice, d$sit))
report("1", "mxl: names(coef(fit)) = names(theta0)",
       as.numeric(!identical(names(coef(fit)), names(th))), 0.5)
report("1", "mxl: simulated loglik, suite vs choicer at theta0 (rel)",
       abs(ll - fit$loglik) / abs(ll), 1e-10)

# --- 2. Calibration --------------------------------------------------------------------

ty <- d2_types()
info(sprintf("types %d, choosers per mode (air bus car train) %s, sum(pi) - 1 = %.1e", ty$n,
             paste(as.integer(ty$n_y), collapse = " "), sum(ty$pi) - 1))
for (v in D2_VARIANTS) {
  sh <- d2_pop_shares(v)
  info(sprintf("%-7s population shares (air bus car train) %s; max |share - Q| %.2e", v,
               paste(sprintf("%.6f", sh), collapse = " "), max(abs(sh - D2_Q[D2_ALTS]))))
}
report("2", "default: population shares at theta0 = Q (WESML FOCs)",
       max(abs(d2_pop_shares("default") - D2_Q[D2_ALTS])), 1e-8)
report("2", "mxl: population shares at theta0 = Q (matched constants)",
       max(abs(d2_pop_shares("mxl") - D2_Q[D2_ALTS])), 1e-8)
# The trapezoid rule on [-12, 12] converges geometrically for this analytic
# integrand: an independent reference for the Gauss-Hermite rule.
P_tr <- local({
  z <- seq(-12, 12, by = 0.02)
  rep_long <- ty$long[rep(seq_len(nrow(ty$long)), length(z))]
  rep_long[, rep_id := rep(seq_len(ty$n * length(z)), each = 4L)]
  d_tr <- design(rep_long, sp, id = "rep_id")
  par <- unpack_theta("mxl", sp, D2_THETA0$mxl)
  bz <- rc_coefs(par$mu, par$L, matrix(z, 1L), par$rc_dist, par$rc_mean)
  p <- mxl_cond_probs(d_tr, par, matrix(rep(bz, each = ty$n), ncol = 1L))$p
  matrix(matrix(p, 4L * ty$n, length(z)) %*% (stats::dnorm(z) * 0.02), ty$n, 4L, byrow = TRUE)
})
report("2", sprintf("mxl: P_true, Gauss-Hermite (%d nodes) vs trapezoid (h = .02)", D2_GH_NODES),
       max(abs(d2_type_probs("mxl") - P_tr)), 1e-12)
star_fit <- d2_theta_star_fit()
report("2", "theta_star refitted from the frozen theta0 vs D2_THETA_STAR (rel)",
       max(abs(coef(star_fit)[names(D2_THETA_STAR)] - D2_THETA_STAR) / abs(D2_THETA_STAR)), 1e-6)

# --- INFO: identification strength and population predictions ------------------------

mc <- d2_wesml_data()
wfit <- function(f, ...) quiet(f(mc, "id", "alt", "choice", weights_col = ".wesml_weight",
                                 se_method = "sandwich", polish = "newton", ...))
nl <- wfit(choicer::run_nestlogit, covariate_cols = D2_X, nest_col = "nest")
info(sprintf(paste("NL (WESML, mode_choice): Lambda_1 = %.4f, robust SE %.4f (Hessian %.4f);",
                   "(1 - lambda) / SE = %.2f"), coef(nl)[["Lambda_1"]], nl$se[["Lambda_1"]],
             sqrt(diag(stats::vcov(nl, type = "hessian")))[["Lambda_1"]],
             (1 - coef(nl)[["Lambda_1"]]) / nl$se[["Lambda_1"]]))
mx <- wfit(choicer::run_mxlogit, covariate_cols = sp$beta, random_var_cols = sp$random,
           rc_dist = 0L, rc_mean = TRUE, S = 5000L, draws = "generate", seed = 1L)
info(sprintf("MXL (WESML, free sigma): sigma = %.3g, sigma / SE(sigma) = 1 / SE(L_11) = %.3f",
             exp(coef(mx)[["L_11"]]), 1 / mx$se[["L_11"]]),
     sprintf(", Mu_travel %.5f", coef(mx)[["Mu_travel"]]))
info(sprintf("MXL truth: Mu_travel %.6f, sigma %.6f (CV %.3f); P(b_travel > 0) = %.3f",
             D2_THETA0$mxl[["Mu_travel"]], exp(D2_THETA0$mxl[["L_11"]]),
             exp(D2_THETA0$mxl[["L_11"]]) / abs(D2_THETA0$mxl[["Mu_travel"]]),
             stats::pnorm(-abs(D2_THETA0$mxl[["Mu_travel"]]) / exp(D2_THETA0$mxl[["L_11"]]))))

# Population information of the MNL at theta0 (pseudo-situations weighted by
# pi_i P_ij(theta0)): its fit returns theta0 (Fisher consistency) and H^-1.
pop <- d2_pseudo_fit(ty$pi * d2_type_probs("default"))
report("2", "MNL population fit at theta0 returns theta0 (Fisher consistency)",
       max(abs(coef(pop) - D2_THETA0$default) / sqrt(diag(stats::vcov(pop, type = "hessian")))),
       1e-6)
Hinv <- stats::vcov(pop, type = "hessian")                 # per traveller (total weight 1)
th0 <- D2_THETA0$default
for (N in c(500, 2000, 8000)) {
  se <- sqrt(diag(Hinv) / N)
  info(sprintf("F1 N = %5d: SE %s; |b_vcost| / SE %.2f", N,
               paste(sprintf("%s %.3g", names(se), se), collapse = ", "),
               abs(th0[["vcost"]]) / se[["vcost"]]))
}
# VTTS = -b_travel / b_vcost (wtp()'s sign) and its delta-method interval
# under b ~ N(theta0, H^-1 / N) (the asymptotic normal of the coefficients).
v0 <- -th0[["travel"]] / th0[["vcost"]]
Vtc <- Hinv[c("travel", "vcost"), c("travel", "vcost")]
B <- vs_with_seed(21L, matrix(stats::rnorm(4e5), ncol = 2L) %*% chol(Vtc))
for (N in c(500, 2000, 8000, 32000, 128000)) {
  bt <- th0[["travel"]] + B[, 1] / sqrt(N)
  bc <- th0[["vcost"]] + B[, 2] / sqrt(N)
  r <- -bt / bc
  g1 <- -1 / bc
  g2 <- bt / bc^2
  se <- sqrt((g1^2 * Vtc[1, 1] + 2 * g1 * g2 * Vtc[1, 2] + g2^2 * Vtc[2, 2]) / N)
  lo <- mean(v0 < r - 1.959964 * se)
  up <- mean(v0 > r + 1.959964 * se)
  info(sprintf("F1 N = %6d: VTTS %.4f, predicted delta interval: tails %.4f / %.4f, coverage %.4f",
               N, v0, lo, up, 1 - lo - up))
}

# F11: the population's Hessian and robust SEs at theta_star.
Vh <- stats::vcov(star_fit, type = "hessian")
Vs <- Vh %*% solve(stats::vcov(star_fit, type = "bhhh")) %*% Vh
rat <- sqrt(diag(Vh) / diag(Vs))
info("F11 predicted Hessian SE / robust SE at theta_star: ",
     paste(sprintf("%s %.3f", names(rat), rat), collapse = ", "))

# F8: per-type MNL scores and Hessians, then the WESML and unweighted
# variances under quotas (H_j = 1/4) and Bernoulli sampling.
f8_pred <- function() {
  K <- length(th0)
  Z <- array(0, c(ty$n, 4L, K))                            # type x mode x parameter
  Zx <- as.matrix(ty$long[, D2_X, with = FALSE])
  for (k in 1:4) {
    Z[, k, seq_along(D2_X)] <- Zx[seq(k, nrow(Zx), by = 4L), ]
    if (k > 1L) Z[, k, length(D2_X) + k - 1L] <- 1
  }
  probs <- function(theta) {
    V <- sapply(1:4, function(k) Z[, k, ] %*% theta)
    E <- exp(V - apply(V, 1, max))
    E / rowSums(E)
  }
  parts <- function(theta) {                               # scores s[i, j, ] and Hessians
    P <- probs(theta)
    zbar <- apply(Z * array(P, dim(Z)), c(1, 3), sum)
    S <- Z - aperm(array(zbar, c(ty$n, K, 4L)), c(1, 3, 2))
    Hm <- lapply(seq_len(ty$n), function(i) crossprod(S[i, , ] * sqrt(P[i, ])))
    list(P = P, S = S, Hm = Hm)
  }
  Q <- D2_Q[D2_ALTS]
  h <- rep(0.25, 4)
  a <- parts(th0)
  Pj <- a$P * ty$pi                                        # pi_i P_ij: the joint of (type, mode)
  Ej <- function(f) {                                      # E[f | y = j] over the types
    lapply(1:4, function(j) {
      Reduce(`+`, lapply(seq_len(ty$n), function(i) Pj[i, j] * f(i, j))) / Q[j]
    })
  }
  mu <- Ej(function(i, j) a$S[i, j, ])
  Ess <- Ej(function(i, j) tcrossprod(a$S[i, j, ]))
  A <- Reduce(`+`, lapply(seq_len(ty$n), function(i) ty$pi[i] * a$Hm[[i]]))
  M_plain <- Reduce(`+`, lapply(1:4, function(j) Q[j]^2 / h[j] * Ess[[j]]))
  M_quota <- Reduce(`+`, lapply(1:4, function(j) Q[j]^2 / h[j] * (Ess[[j]] - tcrossprod(mu[[j]]))))
  Ai <- solve(A)
  v_plain <- diag(Ai %*% M_plain %*% Ai)
  v_quota <- diag(Ai %*% M_quota %*% Ai)
  # Unweighted: the sampled model at the shifted constants; types by mode as
  # in the population (quota) with h_j per mode.
  ths <- th0
  asc <- paste0("ASC_", D2_ALTS[-1L])
  ths[asc] <- th0[asc] + log(h[-1] / Q[-1]) - log(h[1] / Q[1])
  s <- parts(ths)
  Wq <- sweep(a$P * ty$pi, 2, h / Q, `*`)                  # (type, mode) shares in a quota sample
  Hs <- Reduce(`+`, lapply(seq_len(ty$n), function(i) sum(Wq[i, ]) * s$Hm[[i]]))
  mus <- lapply(1:4, function(j) colSums(Wq[, j] * s$S[, j, ]) / h[j])
  Ess_s <- Reduce(`+`, lapply(1:4, function(j) crossprod(s$S[, j, ] * sqrt(Wq[, j]))))
  M_us <- Ess_s - Reduce(`+`, lapply(1:4, function(j) h[j] * tcrossprod(mus[[j]])))
  Hsi <- solve(Hs)
  list(names = names(th0), wesml_robust_quota = sqrt(v_plain / v_quota),
       wesml_hessian_quota = sqrt(diag(Ai) / v_quota),
       wesml_hessian_random = sqrt(diag(Ai) / v_plain),
       unweighted_hessian_quota = sqrt(diag(Hsi) / diag(Hsi %*% M_us %*% Hsi)))
}
pr8 <- f8_pred()
for (k in setdiff(names(pr8), "names")) {
  info(sprintf("F8 predicted SE ratio, %-26s %s", k,
               paste(sprintf("%s %.3f", pr8$names, pr8[[k]]), collapse = ", ")))
}

# F8 random_nl: the NL population's WESML target is theta0_nl (its weighted
# population fit returns it); the unweighted fits' frozen pseudo-true values
# are reproduced; their gaps from the Manski-Lerman shift and the VTTS, in
# robust SE at N = 2000.
f8_env <- environment(load_family("F8", b$L)$simulate)
th_nl <- D2_THETA0$nl
wpop <- d2_pseudo_fit(ty$pi * d2_type_probs("nl"), "nl")
report("2", "NL population fit at theta0_nl returns it (Fisher consistency)",
       max(abs(coef(wpop) - th_nl) / sqrt(diag(stats::vcov(wpop, type = "hessian")))), 1e-6)
pt <- f8_env$f8_pseudo_true()
for (m in names(pt)) {
  fz <- f8_env$F8_NL_PSEUDO_TRUE[[m]]
  report("2", sprintf("F8 random_nl: %s pseudo-true refitted vs frozen (rel)", m),
         max(abs(pt[[m]][names(fz)] - fz) / abs(fz)), 1e-6)
}
W8 <- d2_cb_weights("nl", unname(f8_env$F8_H[D2_ALTS]))
Qn <- d2_pop_shares("nl")
asc <- paste0("ASC_", D2_ALTS[-1L])
shift <- log(0.25 / Qn[-1L]) - log(0.25 / Qn[1L])
vtts_nl <- -th_nl[["travel"]] / th_nl[["vcost"]]
for (m in names(pt)) {
  f8 <- d2_pseudo_fit(W8, if (m == "nl_unweighted") "nl" else "mnl")
  Vh8 <- stats::vcov(f8, type = "hessian")
  V8 <- Vh8 %*% solve(stats::vcov(f8, type = "bhhh")) %*% Vh8 / 2000
  naive <- if (m == "nl_unweighted") th_nl else coef(d2_pseudo_fit(ty$pi * d2_type_probs("nl")))
  naive[asc] <- naive[asc] + shift
  t8 <- coef(f8)
  gap <- (t8 - naive[names(t8)]) / sqrt(diag(V8))
  g <- c(-1 / t8[["vcost"]], t8[["travel"]] / t8[["vcost"]]^2)
  se_v <- sqrt(drop(g %*% V8[c("travel", "vcost"), c("travel", "vcost")] %*% g))
  info(sprintf("F8 random_nl %s: (pseudo-true - Manski-Lerman shifted %s) / SE: %s", m,
               if (m == "nl_unweighted") "theta0_nl" else "population MNL",
               paste(sprintf("%s %.2f", names(gap), gap), collapse = ", ")))
  info(sprintf("    its VTTS %.4f vs the population's %.4f: %.2f SE; Hessian / robust SE %s",
               -t8[["travel"]] / t8[["vcost"]], vtts_nl,
               (-t8[["travel"]] / t8[["vcost"]] - vtts_nl) / se_v,
               paste(sprintf("%.3f", sqrt(diag(Vh8) / 2000 / diag(V8))), collapse = " ")))
}

# --- 3. theta_star exact vs simulated ----------------------------------------------------

t0 <- proc.time()[["elapsed"]]
counts <- matrix(0, ty$n, 4L)
first <- NULL
for (k in 1:10) {
  s <- dgp$simulate(N = 1e6, seed = 100L + k, variant = "mxl")
  if (k == 1L) first <- s$data
  ch <- s$data[choice == 1L]
  counts <- counts + unclass(table(factor(ch$type, levels = seq_len(ty$n)),
                                   factor(ch$alt, levels = D2_ALTS)))
}
# The polish's decrement is over the mean weight (about 1200 on the counts
# of 10^6), so both fits compared below are tightened.
tight <- list(tol_decrement = 1e-15)
sim_fit <- d2_pseudo_fit(counts, polish_control = tight)
dz <- coef(sim_fit)[names(D2_THETA_STAR)] - D2_THETA_STAR
sd_pop <- sqrt(diag(Vs) / 1e7)                             # the exact sampling SD at N = 10^7
W <- drop(dz %*% solve(Vs / 1e7, dz))
info(sprintf("10^7 travellers simulated and fitted in %.1f s; (theta_hat - theta_star) / SD: %s",
             proc.time()[["elapsed"]] - t0,
             paste(sprintf("%s %.2f", names(dz), dz / sd_pop), collapse = ", ")))
report("3", "theta_star vs the simulated N = 10^7 fit: max |diff| / SD", max(abs(dz / sd_pop)), 4)
report("3", sprintf("same, joint: W = %.2f on 7 df, its p-value", W),
       stats::pchisq(W, 7, lower.tail = FALSE), 1e-3,
       pass = stats::pchisq(W, 7, lower.tail = FALSE) > 1e-3)
c1 <- unclass(table(factor(first[choice == 1L]$type, levels = seq_len(ty$n)),
                    factor(first[choice == 1L]$alt, levels = D2_ALTS)))
agg <- d2_pseudo_fit(c1, polish_control = tight)
t0 <- proc.time()[["elapsed"]]
full <- quiet(choicer::run_mnlogit(first, "id", "alt", "choice", D2_X, polish = "newton",
                                   polish_control = tight))
info(sprintf("unaggregated fit of the first 10^6 (4e6 rows): %.1f s; decrements %.1e, %.1e",
             proc.time()[["elapsed"]] - t0, fit_diag(agg)$decrement, fit_diag(full)$decrement))
report("3", "first 10^6: fit on counts = unaggregated fit (max |diff| / SE)",
       max(abs(coef(agg) - coef(full)) / full$se), 1e-6)
rm(first, full)

# --- 4. The exact constant shift under choice-based sampling ----------------------------

frame <- dgp$simulate(N = 400000L, seed = 31L, variant = "default")$data
smp <- vs_with_seed(32L, choicer::sample_by_choice(frame, "id", "alt", "choice",
                                                   n_per_alt = 25000L))
for (a in c("Q", "H", "choice_sampling")) data.table::setattr(smp, a, NULL)
smp[, .wesml_weight := NULL]
cb <- quiet(choicer::run_mnlogit(smp, "id", "alt", "choice", D2_X, polish = "newton"))
asc <- paste0("ASC_", D2_ALTS[-1L])
shifted <- th0
shifted[asc] <- th0[asc] + log(0.25 / D2_Q[D2_ALTS[-1L]]) - log(0.25 / D2_Q[["air"]])
zc <- (coef(cb) - shifted[names(coef(cb))]) / cb$se
info(sprintf("choice-based N = %d: (theta_hat - shifted truth) / SE: %s", sum(smp$choice),
             paste(sprintf("%s %.2f", names(zc), zc), collapse = ", ")))
report("4", "unweighted MNL slopes within 3 SE of theta0", max(abs(zc[D2_X])), 3)
report("4", "unweighted MNL constants within 3 SE of the shifted truth", max(abs(zc[asc])), 3)
rm(frame, smp)

# --- 5. References -----------------------------------------------------------------------

for (r in list(ref_q11(), ref_q12())) {
  for (i in seq_len(nrow(r))) report("5", paste(r$test[i], r$check[i]), r$value[i], r$threshold[i])
}

# --- 6. Smoke: F1, F8, F11 at reps 1 and 2 ------------------------------------------------

tab <- list()
for (id in c("F1", "F8", "F11")) {
  fam <- load_family(id, b$L)
  for (cl in fam$cells$cell) for (r in 1:2) {
    recs <- run_rep(fam, cl, r, cell_methods(fam, cl))
    for (rec in recs) {
      fin <- isTRUE(rec$ok) && all(is.finite(rec$est$est)) && all(is.finite(rec$est$se_hessian))
      tab[[length(tab) + 1L]] <- data.frame(
        family = id, cell = cl, method = rec$method, rep = r, ok = rec$ok, finite = fin,
        sec = round(rec$diag$unit_sec, 2), sim_sec = round(rec$diag$sim_sec, 2),
        iterations = rec$diag$iterations, code = rec$diag$code,
        decrement = signif(rec$diag$decrement, 2), warnings = nrow(rec$warnings),
        error = if (is.null(rec$error)) "" else substr(rec$error, 1, 60), stringsAsFactors = FALSE)
    }
  }
}
tab <- do.call(rbind, tab)
print(tab, row.names = FALSE)
report("6", "smoke: every record ok with finite estimates and Hessian SEs",
       sum(!tab$finite), 0.5)

cat(sprintf("\n%d checks, %d FAIL; %.1f s\n", length(verdicts), sum(!verdicts),
            proc.time()[["elapsed"]] - t_start))
if (any(!verdicts)) quit(status = 1L)
