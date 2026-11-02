# ============================================================================
# Null and power studies of the suite's metrics (WP1). Fake rep records
# (lib/fit.R's new_rep_record(); no model is fitted) with known sampling
# distributions test lib/metrics.R and lib/sbc.R:
#   0  the formulas, gating (auto-gated failure rates included) and both
#      gate rules against direct computation; route, joint Wald and global
#      maximum failure rates on constructed reps;
#   1  MCSE calibration: sd and mean of (value - true value) / mcse;
#   2  gate size, rules "plan" and "band": an exact null family of 33 gated
#      rows fails in at most 2% of runs (a noisy-SE null reported);
#   3  gate power against the deviations the tolerances name: the measured
#      rate equals the design power (closed form), under both rules;
#   4  joint Wald and LR size and mean under the null;
#   5  coverage pooled over 2,000 constants with a shared shift per rep;
#   6  the exact SBC ECDF test: its p against a simulated null, uniform p,
#      the band, power against misfits; families of 30 SBC quantities;
#   7  rate_fit() recovers the slope -1/2 of ESE in N;
#   8  edge cases: R = 1 and 2, every rep failed, NA routes, truths and tests.
# A meta-replication is one simulated set of R reps. Prints PASS, FAIL or
# INFO (reported, not judged) per claim with the measured numbers, and exits
# 1 on any FAIL. Seeded: every run prints the same numbers.
#
# Usage (full under 8 minutes; --fast, fewer meta-replications and bands
# widened to 3.5 of their SEs, under 1 minute; --only= runs some sections):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 \
#     Rscript _validation/checks/wp1_metrics_null.R [--fast] [--only=0,3] [--seed=N]
# ============================================================================

root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                               value = TRUE)[1]))))
source(file.path(root, "lib", "families.R"), local = (boot <- new.env()))
attach(boot$vs_source_lib(file.path(root, "lib")), name = "vs_lib", warn.conflicts = FALSE)
opt <- vs_args(defaults = list(fast = FALSE, only = "", seed = 20261102L))
only <- if (nzchar(opt$only)) strsplit(opt$only, ",", fixed = TRUE)[[1]] else as.character(0:8)

# Meta-replications per study.
M <- if (opt$fast) {
  list(null100 = 150, null1000 = 40, pooled = 50, sbc = 200, sbc_power = 30, rate = 150,
       sbc_fam = c("100" = 40, "1000" = 20))
} else {
  list(null100 = 1000, null1000 = 500, pooled = 300, sbc = 1000, sbc_power = 200, rate = 500,
       sbc_fam = c("100" = 500, "1000" = 300))
}
LEVEL <- 0.95
ZQ <- qnorm(1 - (1 - LEVEL) / 2)
A2 <- (1 - LEVEL) / 2                    # each tail's target
NU <- 50                                 # se_hessian ~ sigma sqrt(chi^2_NU / NU)
NU_ROB <- 30                             # se_robust, a noisier route
L_SBC <- 99L

# --- Reporting --------------------------------------------------------------

results <- character()
# One verdict line, then the details (one line per element of `...` pasted).
report <- function(name, ok, ...) {
  verdict <- if (is.na(ok)) "INFO" else if (isTRUE(ok)) "PASS" else "FAIL"
  cat(sprintf("%-4s  %s\n", verdict, name))
  det <- unlist(strsplit(paste0(...), "\n", fixed = TRUE))
  if (length(det)) cat(paste0("        ", det), sep = "\n")
  results[name] <<- verdict
}
section <- function(id, title, expr) {
  if (!id %in% only) return(invisible())
  cat(sprintf("\n== %s  %s\n", id, title))
  vs_set_seed(opt$seed + as.integer(id))
  t0 <- proc.time()[["elapsed"]]
  tryCatch(expr, error = function(e) report(paste0(id, " stopped by an error"), FALSE,
                                            conditionMessage(e)))
  cat(sprintf("        [%.1f s]\n", proc.time()[["elapsed"]] - t0))
}
f3 <- function(x) formatC(x, digits = 3, format = "f")
f4 <- function(x) formatC(x, digits = 4, format = "f")
# A proportion out of n meta-replications, with its binomial SE.
rate_se <- function(x) sprintf("%s (SE %s, n %d)", f4(mean(x)), f4(sqrt(mean(x) * (1 - mean(x)) /
                                                                       length(x))), length(x))
# Acceptance bands for the sd and mean of n z-scores: the brief's 0.15, or
# 3.5 SEs where --fast leaves too few meta-replications for it.
band <- function(n) {
  if (!opt$fast) return(c(sd = 0.15, mean = 0.15))
  c(sd = max(0.15, 3.5 / sqrt(2 * n)), mean = max(0.15, 3.5 / sqrt(n)))
}

# --- Fake records -------------------------------------------------------------

# A record skeleton of one cell: estimands `name` in blocks `block` with
# truths `truth`, built by new_rep_record() once; fill_records() copies it.
skeleton <- function(cell, name, block, truth, family = "N", method = "m") {
  est <- data.frame(name = name, block = block, est = truth, truth = truth, se_hessian = NA_real_,
                    se_bhhh = NA_real_, se_robust = NA_real_, se_cluster = NA_real_, gated = TRUE,
                    stringsAsFactors = FALSE)
  new_rep_record(family, cell, method, rep = 1, seed = NA_real_, est = est)
}

# R records of one cell from a skeleton: rep r takes row r of `est` (R x k),
# of each `se` matrix (route -> R x k) and of `truth`, and the r-th joint
# Wald and LR tests (lists of equal-length vectors W or LR, df, p).
fill_records <- function(sk, est, se = list(), truth = NULL, wald = NULL, lr = NULL) {
  base <- unclass(sk$est)
  lapply(seq_len(nrow(est)), function(r) {
    x <- base
    x$est <- est[r, ]
    for (k in names(se)) x[[paste0("se_", k)]] <- se[[k]][r, ]
    if (!is.null(truth)) x$truth <- truth[r, ]
    sk$rep <- r
    sk$est <- structure(x, class = "data.frame")
    if (!is.null(wald)) sk$joint_wald <- list(W = wald$W[r], df = wald$df[r], p = wald$p[r])
    if (!is.null(lr)) sk$lr <- list(LR = lr$LR[r], df = lr$df[r], p = lr$p[r])
    sk
  })
}

# Chi-square tests of statistic `stat` (named W or LR) with `df` degrees.
chisq_tests <- function(stat, df, name = "W") {
  stats::setNames(list(stat, rep(df, length(stat)), pchisq(stat, df, lower.tail = FALSE)),
                  c(name, "df", "p"))
}

# Rows of a metrics table by metric, estimand (NA: cell-level) and route.
rows_of <- function(tab, metric, estimand = NA, route = NA) {
  m <- tab$metric == metric
  m <- m & (if (is.na(estimand)) is.na(tab$estimand) else tab$estimand %in% estimand)
  m & (if (is.na(route)) is.na(tab$route) else tab$route %in% route)
}
row_key <- function(metric, estimand, route) paste(metric, estimand, route, sep = "|")
# pooled_tails()'s rows of one tail ("tail_lower"), named with or without a
# "pooled_" prefix; a group, when given, is read from `group` or `estimand`.
pooled_rows <- function(tab, metric, group = NULL) {
  m <- sub("^pooled_", "", tab$metric) == metric
  if (is.null(group)) return(m)
  m & (if (!is.null(tab$group)) tab$group %in% group else tab$estimand %in% group)
}

# --- Exact moments under the null ----------------------------------------------

# E[s] / sigma for the sd of n normal values; E[sigma / s]; E sqrt(chi^2_k / k).
c4 <- function(n) sqrt(2 / (n - 1)) * exp(lgamma(n / 2) - lgamma((n - 1) / 2))
e_inv_s <- function(n) sqrt((n - 1) / 2) * exp(lgamma((n - 2) / 2) - lgamma((n - 1) / 2))
e_sqrt_chi <- function(k) sqrt(2 / k) * exp(lgamma((k + 1) / 2) - lgamma(k / 2))

# The error of estimand `tail04`: X = (Z + a (Z^2 - 1)) / sqrt(1 + 2 a^2),
# mean 0, variance 1, right-skewed; P(X > c) and P(X < -c) from the roots of
# a Z^2 + Z - a -/+ c sqrt(1 + 2 a^2). With se = sigma, a sets tail_lower
# P(X > z) to 0.04.
skew_right <- function(a, c) {
  cc <- c * sqrt(1 + 2 * a^2)
  d <- sqrt(1 + 4 * a * (a + cc))
  pnorm((-1 + d) / (2 * a), lower.tail = FALSE) + pnorm((-1 - d) / (2 * a))
}
skew_left <- function(a, c) {
  cc <- c * sqrt(1 + 2 * a^2)
  d <- sqrt(pmax(1 + 4 * a * (a - cc), 0))
  pnorm((-1 + d) / (2 * a)) - pnorm((-1 - d) / (2 * a))
}
A_SKEW <- uniroot(function(a) skew_right(a, ZQ) - 0.04, c(1e-3, 0.5), tol = 1e-12)$root

# --- The null family ----------------------------------------------------------
# Cell A: six estimands at several scales, b1-b3 equicorrelated (0.6), b4-b5
# correlated -0.4; joint Wald on the true V (W ~ chi^2_6), LR ~ chi^2_2. Two
# nulls share these draws and tests: "exact", se_hessian = sigma (every
# metric's true value is its target), and "noisy", se_hessian with chi^2_50
# noise and se_robust with chi^2_30 noise (realistic SEs: the Wald tails are
# P(t_nu > z), 0.0278 for the Hessian, not 0.025).
# Cell B (exact SEs): estimand `null` and five deviations, all with sigma 2:
# se90 and se97 (se scaled 0.90, 0.97), bias20 and bias05 (bias 0.20 and
# 0.05 sigma), tail04 (skewed: tail_lower 0.04, bias 0, se ratio 1); its
# joint Wald has size 0.08 (W = c chi^2_6). Its `null` estimand with noisy
# SEs completes the noisy family.
# Gated: bias_std, se_ratio, both tails on the Hessian route, joint Wald and
# LR size (and, auto-gated, each cell's Hessian route_failure_rate and A's
# joint Wald failure rate), so that A plus one B estimand is a family of
# N_GATED = 33 gated rows.

A_NAMES <- paste0("b", 1:6)
A_BLOCK <- c("beta", "beta", "beta", "mu", "sigma", "lambda")
A_THETA <- c(1, -0.5, 2, 0, 0.3, -1)
A_SD <- c(0.1, 0.1, 1, 10, 0.01, 1)
A_CORR <- diag(6)
A_CORR[1:3, 1:3] <- 0.6
A_CORR[4, 5] <- A_CORR[5, 4] <- -0.4
diag(A_CORR) <- 1
A_CHOL <- chol(A_CORR)
B_NAMES <- c("null", "se90", "se97", "bias20", "bias05", "tail04")
B_SE_SCALE <- c(1, 0.90, 0.97, 1, 1, 1)
B_SHIFT <- c(0, 0, 0, 0.20, 0.05, 0)
B_THETA <- 0.5
B_SD <- 2
WALD08 <- qchisq(0.95, 6) / qchisq(0.92, 6)
GATES <- list(cell = c("A", "B"),
              metrics = c("bias_std", "se_ratio:hessian", "tail_lower:hessian",
                          "tail_upper:hessian", "joint_wald_size", "lr_size"))
SK_A <- skeleton("A", A_NAMES, A_BLOCK, A_THETA)
SK_B <- skeleton("B", B_NAMES, rep("beta", 6), rep(B_THETA, 6))
SK_B0 <- skeleton("B", "null", "beta", B_THETA)
TESTS <- c("joint_wald_size", "joint_wald_failure_rate", "lr_size", "lr_mean_q")
# Gated rows of the null families: A's 24 estimand rows, its Hessian
# route_failure_rate, joint Wald size and failure rate and LR size (the
# failure rates auto-gated), and B's 4 estimand rows and route_failure_rate.
N_GATED <- 33L

# Cell A's records under both nulls; the exact ones carry no tests (they
# are the noisy ones', and metrics_cell() spends most of its time on tests).
draw_A <- function(R) {
  E <- matrix(rnorm(R * 6), R)
  sd <- rep(A_SD, each = R)
  est <- rep(A_THETA, each = R) + (E %*% A_CHOL) * sd
  seh <- sqrt(matrix(rchisq(R * 6, NU), R) / NU) * sd
  ser <- sqrt(matrix(rchisq(R * 6, NU_ROB), R) / NU_ROB) * sd
  list(noisy = fill_records(SK_A, est, list(hessian = seh, robust = ser),
                            wald = chisq_tests(rowSums(E^2), 6),
                            lr = chisq_tests(rchisq(R, 2), 2, "LR")),
       exact = fill_records(SK_A, est, list(hessian = matrix(sd, R))))
}
# Cell B's records: the deviations with exact SEs, and `null` with noisy SEs.
draw_B <- function(R) {
  Z <- matrix(rnorm(R * 6), R)
  Z[, 6] <- (Z[, 6] + A_SKEW * (Z[, 6]^2 - 1)) / sqrt(1 + 2 * A_SKEW^2)
  est <- B_THETA + B_SD * (Z + rep(B_SHIFT, each = R))
  se_noisy <- B_SD * sqrt(rchisq(R, NU) / NU)
  list(exact = fill_records(SK_B, est, list(hessian = matrix(B_SD * rep(B_SE_SCALE, each = R), R)),
                            wald = chisq_tests(WALD08 * rchisq(R, 6), 6)),
       noisy = fill_records(SK_B0, est[, 1, drop = FALSE], list(hessian = matrix(se_noisy))))
}

# The families of one meta-replication: under the exact null, "null" (A, B's
# null estimand and B's route_failure_rate, N_GATED gated rows), one per B
# deviation (A, that estimand and B's route_failure_rate), and "wald08" (the
# null family and B's joint Wald size and failure rate, N_GATED + 2);
# "noisy" (the null family under the noisy null).
FAMILIES <- c(null = "null", se90 = "se90", se97 = "se97", bias20 = "bias20",
              bias05 = "bias05", tail04 = "tail04", wald08 = "null")
stack_families <- function(tA, tB, tA_noisy, tB_noisy) {
  ex <- lapply(names(FAMILIES), function(f) {
    keep <- tB$estimand %in% FAMILIES[[f]] | tB$metric == "route_failure_rate" |
      (f == "wald08" & tB$metric %in% c("joint_wald_size", "joint_wald_failure_rate"))
    t <- rbind(tA, tB[keep, , drop = FALSE])
    t$family <- f
    t
  })
  noisy <- rbind(tA_noisy, tB_noisy)
  noisy$family <- "noisy"
  do.call(rbind, c(ex, list(noisy)))
}
# The row each deviation is meant to move.
DEV_ROW <- list(se90 = c("se_ratio", "se90", "hessian"), se97 = c("se_ratio", "se97", "hessian"),
                bias20 = c("bias_std", "bias20", NA), bias05 = c("bias_std", "bias05", NA),
                tail04 = c("tail_lower", "tail04", "hessian"),
                wald08 = c("joint_wald_size", NA, NA))
RULES <- c("plan", "band")

CAL_METRICS <- c("bias_std", "ese", "se_ratio", "tail_lower", "tail_upper", "coverage",
                 "joint_wald_size", "lr_size", "lr_mean_q")

# One null study: M meta-replications of cells A and B at R reps. Returns
# the calibration rows (the noisy null's estimands and tests), and per gate
# rule each family's failure, each deviant row's failure and the rows that
# failed in the exact "null" family; and the gated count of each family.
null_study <- function(R, M, tol_set) {
  cal <- vector("list", M)
  fam <- dev <- fails <- stats::setNames(rep(list(vector("list", M)), 2), RULES)
  n_gated <- NULL
  for (m in seq_len(M)) {
    rA <- draw_A(R)
    rB <- draw_B(R)
    tA_noisy <- metrics_cell(rA$noisy, GATES, tol_set, LEVEL)
    tA <- metrics_cell(rA$exact, GATES, tol_set, LEVEL)   # its tests: "na" placeholders
    tA <- rbind(tA[!tA$metric %in% TESTS, , drop = FALSE],
                tA_noisy[tA_noisy$metric %in% TESTS, , drop = FALSE])
    tB <- metrics_cell(rB$exact, GATES, tol_set, LEVEL)
    tB_noisy <- metrics_cell(rB$noisy, GATES, tol_set, LEVEL)
    tB_noisy <- tB_noisy[!is.na(tB_noisy$estimand) | tB_noisy$metric == "route_failure_rate", ,
                         drop = FALSE]                  # without its tests' placeholders
    cal[[m]] <- rbind(tA_noisy[tA_noisy$metric %in% CAL_METRICS, CAL_COLS],
                      tB_noisy[tB_noisy$metric %in% CAL_METRICS, CAL_COLS])
    st <- stack_families(tA, tB, tA_noisy, tB_noisy)
    if (m == 1L) n_gated <- tapply(st$gated, st$family, sum)
    for (rule in RULES) {
      g <- gate_family(st, rule = rule)
      fail <- g$verdict %in% "fail"
      fam[[rule]][[m]] <- vapply(c(names(FAMILIES), "noisy"),
                                 function(f) any(fail[g$family == f]), NA)
      dev[[rule]][[m]] <- vapply(names(DEV_ROW), function(f) {
        d <- DEV_ROW[[f]]
        i <- g$family == f & rows_of(g, d[1], d[2], d[3]) & (f != "wald08" | g$cell == "B")
        if (sum(i) != 1L) stop("family ", f, ": ", sum(i), " rows match ", toString(d))
        fail[i]
      }, NA)
      i <- fail & g$family == "null"
      fails[[rule]][[m]] <- row_key(g$metric[i], g$estimand[i], g$route[i])
    }
  }
  bind <- function(x) lapply(x, function(l) do.call(rbind, l))
  list(cal = do.call(rbind, cal), fam = bind(fam), dev = bind(dev),
       null_fails = lapply(fails, unlist), n_gated = n_gated, R = R, M = M)
}
CAL_COLS <- c("metric", "estimand", "route", "value", "mcse")

# The finite-R expectation E[value] of each calibration row under the noisy
# null (exact: chi-square and gamma moments, the t tail), and the population
# value.
true_values <- function(cal, R) {
  nu <- ifelse(cal$route %in% "robust", NU_ROB, NU)
  sd <- c(stats::setNames(A_SD, A_NAMES), null = B_SD)[cal$estimand]
  tail <- pt(ZQ, nu, lower.tail = FALSE)
  pop <- with(cal, ifelse(metric == "bias_std", 0, ifelse(metric == "ese", sd,
                   ifelse(metric == "se_ratio", 1, ifelse(metric %in% c("tail_lower", "tail_upper"),
                   tail, ifelse(metric == "coverage", 1 - 2 * tail,
                   ifelse(metric %in% c("joint_wald_size", "lr_size"), 0.05, 1)))))))
  exp_r <- ifelse(cal$metric == "ese", sd * c4(R),
                  ifelse(cal$metric == "se_ratio", e_sqrt_chi(R * nu) * e_inv_s(R), pop))
  list(expected = unname(exp_r), population = unname(pop))
}

# --- Design power: P(a gated row fails) in closed form ----------------------------
# The row is Holm's smallest of k gated p-values (the other rows null), so it
# fails when its p < alpha / k and, under the plan rule, |value - target| >
# tol. Each p is metrics.R's (accepted by the lead): exact binomial for
# proportions, t with R - 1 df for bias_std, normal otherwise; under the
# band rule, one-sided against the nearer edge of [target - tol, target + tol].

# A proportion whose count is Binomial(R, v).
power_prop <- function(v, target, tol, R, rule, k = 30, alpha = 0.01) {
  x <- 0:R
  dev <- x / R - target
  out <- abs(dev) > tol * (1 + 1e-9)
  p <- if (rule == "plan") {
    pmin(1, 2 * pmin(pbinom(x, R, target), pbinom(x - 1, R, target, lower.tail = FALSE)))
  } else {
    ifelse(dev > 0, pbinom(x - 1, R, min(target + tol, 1), lower.tail = FALSE),
           pbinom(x, R, max(target - tol, 0)))
  }
  sum(dbinom(x[out & p < alpha / k], R, v))
}
# bias_std = b with sqrt(R) b ~ noncentral t_{R-1}(delta sqrt(R)) (bias delta
# sigma, exact SEs): it fails when |t| exceeds a critical c.
power_bias <- function(delta, tol, R, rule, k = 30, alpha = 0.01) {
  c <- if (rule == "plan") max(qt(1 - alpha / (2 * k), R - 1), tol * sqrt(R)) else
    tol * sqrt(R) + qt(1 - alpha / k, R - 1)
  ncp <- delta * sqrt(R)
  pt(c, R - 1, ncp, lower.tail = FALSE) + pt(-c, R - 1, ncp)
}
# se_ratio = r = kse sigma / s (exact SEs scaled by kse), MCSE r / sqrt(2 (R
# - 1)); (R - 1) s^2 / sigma^2 ~ chi^2_{R-1}, so r below r_lo or above r_hi
# fails.
power_ratio <- function(kse, tol, R, rule, k = 30, alpha = 0.01) {
  q <- sqrt(2 * (R - 1))
  if (rule == "plan") {
    z <- qnorm(1 - alpha / (2 * k))
    lo <- min(1 - tol, 1 / (1 + z / q))
    hi <- if (q > z) max(1 + tol, 1 / (1 - z / q)) else Inf
  } else {
    z <- qnorm(1 - alpha / k)
    lo <- (1 - tol) / (1 + z / q)
    hi <- if (q > z) (1 + tol) * q / (q - z) else Inf
  }
  pchisq((R - 1) * (kse / lo)^2, R - 1, lower.tail = FALSE) +
    pchisq((R - 1) * (kse / hi)^2, R - 1)
}

# --- 0. Formulas, gating and the gate rule ----------------------------------------

# The fixed record set: 44 reps, 4 failed; b1, b2 (beta, correlated), s1
# (sigma; est$gated FALSE in rep 3); derived `wtp` with a truth per rep;
# routes hessian and robust (bhhh, cluster all NA); one joint Wald and one
# LR test per rep; a multi-start check in 8 reps, 1 of them a missed
# global maximum (dlogL > 0.01 and max |dtheta| / SE > 0.1).
fixed_records <- function(method = "m1", cell = "A", family = "F0", wald_list = FALSE) {
  vs_with_seed(11, {
    fails <- c(5L, 17L, 30L, 44L)
    ms <- c(0, 0, 2e-6, 0, 1e-7, 5, 0, 0.02)
    mt <- c(0, 0, 1, 0, 0, 2, 0, 0.05)
    lapply(1:44, function(r) {
      if (r %in% fails) {
        return(new_rep_record(family, cell, method, r, NA_real_, ok = FALSE, error = "injected"))
      }
      th <- c(1, -1, 0.3)
      sd <- c(0.5, 0.2, 0.1)
      z <- rnorm(3)
      z[2] <- 0.7 * z[1] + sqrt(1 - 0.49) * z[2]
      est <- data.frame(name = c("b1", "b2", "s1"), block = c("beta", "beta", "sigma"),
                        est = th + sd * z, truth = th,
                        se_hessian = sd * sqrt(rchisq(3, 5) / 5), se_bhhh = NA_real_,
                        se_robust = sd * sqrt(rchisq(3, 3) / 3), se_cluster = NA_real_,
                        gated = c(TRUE, TRUE, r != 3L), stringsAsFactors = FALSE)
      wtp_truth <- 2 + rnorm(1, sd = 0.1)
      derived <- data.frame(name = "wtp", est = wtp_truth + rnorm(1, sd = 0.4),
                            se = 0.4 * sqrt(rchisq(1, 4) / 4), truth = wtp_truth)
      W <- 1.5 * rchisq(1, 3)
      jw <- list(W = W, df = 3, p = pchisq(W, 3, lower.tail = FALSE))
      if (wald_list) {
        W2 <- rchisq(1, 2)
        jw <- list(common = jw, beta = list(W = W2, df = 2, p = pchisq(W2, 2, lower.tail = FALSE)))
      }
      LR <- rchisq(1, 2)
      i <- match(r, c(1:4, 6:9))
      new_rep_record(family, cell, method, r, NA_real_, est = est, derived = derived,
                     joint_wald = jw,
                     lr = list(LR = LR, df = 2, p = pchisq(LR, 2, lower.tail = FALSE)),
                     diag = list(multistart = if (!is.na(i)) {
                       list(dlogL = ms[i], max_dtheta_se = mt[i], n_failed = 0L, n_warned = 0L)
                     }))
    })
  })
}

# Constructed reps of cell A (family "F0", method "m1"): b1, b2 with
# se_hessian 0.1 and se_robust 0.12 (bhhh, cluster never computed), a joint
# Wald test; `se_h` / `se_r` name, per estimand, the reps missing that SE;
# W is NA in reps `w_na` and infinite in reps `w_inf`.
fail_records <- function(R, se_h = list(), se_r = list(), w_na = integer(), w_inf = integer()) {
  vs_with_seed(13, lapply(seq_len(R), function(r) {
    est <- data.frame(name = c("b1", "b2"), block = "beta", est = c(1, -1) + rnorm(2, sd = 0.1),
                      truth = c(1, -1), se_hessian = 0.1, se_bhhh = NA_real_, se_robust = 0.12,
                      se_cluster = NA_real_, gated = TRUE, stringsAsFactors = FALSE)
    for (nm in names(se_h)) if (r %in% se_h[[nm]]) est$se_hessian[est$name == nm] <- NA
    for (nm in names(se_r)) if (r %in% se_r[[nm]]) est$se_robust[est$name == nm] <- NA
    W <- if (r %in% w_na) NA_real_ else if (r %in% w_inf) Inf else rchisq(1, 2)
    new_rep_record("F0", "A", "m1", r, NA_real_, est = est,
                   joint_wald = list(W = W, df = 2, p = pchisq(W, 2, lower.tail = FALSE)))
  }))
}

# The metrics of a cell computed directly from the API's definitions: value,
# the proportion MCSEs at the estimate, at the target and at the smoothed
# estimate (x + 1) / (R + 2), the target, tol.
oracle_cell <- function(records, tol_set = "mc") {
  tol <- tolerances(tol_set)
  ok <- vapply(records, `[[`, NA, "ok")
  recs <- records[ok]
  R <- length(recs)
  out <- list()
  add <- function(metric, estimand, route, value, mcse, target, tol, mcse_t = mcse,
                  mcse_s = mcse) {
    out[[length(out) + 1L]] <<- data.frame(metric, estimand, route, value, mcse, mcse_t, mcse_s,
                                           target, tol, stringsAsFactors = FALSE)
  }
  prop <- function(metric, estimand, route, hit, target, tol) {
    v <- mean(hit)
    s <- (sum(hit) + 1) / (R + 2)
    add(metric, estimand, route, v, sqrt(v * (1 - v) / R), target, tol,
        sqrt(target * (1 - target) / R), sqrt(s * (1 - s) / R))
  }
  one <- function(name, est, truth, ses) {
    d <- est - truth
    ese <- sd(d)
    add("bias_std", name, NA, mean(d) / ese, 1 / sqrt(R), 0, tol[["bias_std"]])
    add("ese", name, NA, ese, ese / sqrt(2 * (R - 1)), NA, NA)
    for (r in names(ses)) {
      s <- ses[[r]]
      if (all(is.na(s))) next
      mod <- sqrt(mean(s^2))
      ratio <- mod / ese      # Morris Table 6, relative error in ModSE, by the delta method
      add("se_ratio", name, r, ratio,
          ratio * sqrt(var(s^2) / (4 * R * mod^4) + 1 / (2 * (R - 1))), 1, tol[["se_ratio"]])
      lo <- truth < est - ZQ * s
      up <- truth > est + ZQ * s
      prop("tail_lower", name, r, lo, A2, tol[["tail"]])
      prop("tail_upper", name, r, up, A2, tol[["tail"]])
      prop("coverage", name, r, !lo & !up, LEVEL, tol[["coverage"]])
    }
  }
  E <- do.call(rbind, lapply(recs, `[[`, "est"))
  for (nm in unique(E$name)) {
    e <- E[E$name == nm, ]
    one(nm, e$est, e$truth, list(hessian = e$se_hessian, bhhh = e$se_bhhh,
                                  robust = e$se_robust, cluster = e$se_cluster))
  }
  D <- do.call(rbind, lapply(recs, `[[`, "derived"))
  for (nm in unique(D$name)) {
    e <- D[D$name == nm, ]
    one(nm, e$est, e$truth, list(hessian = e$se))
  }
  # A route fails in a rep whose est table lacks a finite SE on it for any
  # estimand; a row per route with an SE in some rep.
  for (r in c("hessian", "bhhh", "robust", "cluster")) {
    fin <- lapply(recs, function(x) is.finite(x$est[[paste0("se_", r)]]))
    if (any(unlist(fin))) {
      prop("route_failure_rate", NA, r, !vapply(fin, all, NA), 0, tol[["failure"]])
    }
  }
  W <- vapply(recs, function(r) r$joint_wald$W, 0)
  wp <- vapply(recs, function(r) r$joint_wald$p, 0)
  prop("joint_wald_size", NA, NA, wp < 0.05, 0.05, tol[["size"]])
  prop("joint_wald_failure_rate", NA, NA, !is.finite(W), 0, tol[["failure"]])
  lp <- vapply(recs, function(r) r$lr$p, 0)
  lq <- vapply(recs, function(r) r$lr$LR / r$lr$df, 0)
  prop("lr_size", NA, NA, lp < 0.05, 0.05, tol[["size"]])
  add("lr_mean_q", NA, NA, mean(lq), sd(lq) / sqrt(R), 1, NA)
  add("failure_rate", NA, NA, mean(!ok), NA, NA, NA)
  ms <- Filter(Negate(is.null), lapply(recs, function(r) r$diag$multistart))
  add("global_max_rate", NA, NA, missed_rate(ms), NA, 0, NA)
  do.call(rbind, out)
}
# The share of multi-start checks (with a dlogL) that found a missed global
# maximum: dlogL > 0.01 and max |dtheta| / SE > 0.1.
missed_rate <- function(ms) {
  dl <- vapply(ms, function(m) as.numeric(m$dlogL %||% NA)[1], 0)
  dt <- vapply(ms, function(m) as.numeric(m$max_dtheta_se %||% NA)[1], 0)
  k <- !is.na(dl)
  mean(dl[k] > 0.01 & !is.na(dt[k]) & dt[k] > 0.1)
}

PROP_METRICS <- c("tail_lower", "tail_upper", "coverage", "joint_wald_size", "lr_size",
                  "joint_wald_failure_rate", "route_failure_rate")

# The gate rule applied directly (plan §2, api_metrics.md): p = 2 pnorm(-|value
# - target| / mcse) (a zero MCSE: 1 without a deviation, 0 with one); with
# exact = TRUE, small-R references instead: bias_std's p from t with R - 1
# df, a proportion's the exact two-sided binomial p (twice the smaller tail
# of its count out of R).
oracle_gate <- function(tab, alpha = 0.01, exact = FALSE) {
  dev <- abs(tab$value - tab$target)
  zz <- ifelse(tab$mcse > 0, dev / tab$mcse, ifelse(dev == 0, 0, Inf))
  p <- 2 * pnorm(-zz)
  if (exact) {
    s <- tab$metric == "bias_std" & !is.na(zz) & tab$R >= 2
    p[s] <- 2 * pt(-zz[s], tab$R[s] - 1)
    i <- tab$metric %in% PROP_METRICS & !is.na(tab$value)
    x <- round(tab$value[i] * tab$R[i])
    p[i] <- pmin(1, 2 * pmin(pbinom(x, tab$R[i], tab$target[i]),
                             pbinom(x - 1, tab$R[i], tab$target[i], lower.tail = FALSE)))
  }
  oracle_verdict(tab, p, alpha, dev > tab$tol * (1 + 1e-9))
}
# The band rule (the lead's, 2026-11-02): p one-sided that the true value
# lies in [target - tol, target + tol], against the edge nearer the
# estimate (1 inside the band): normal, t with R - 1 df for bias_std, the
# exact binomial tail at the edge for proportions; fail where p_holm < alpha.
oracle_band <- function(tab, alpha = 0.01) {
  dev <- tab$value - tab$target
  out <- !is.na(dev) & !is.na(tab$tol) & abs(dev) > tab$tol * (1 + 1e-9)
  p <- ifelse(is.na(dev) | is.na(tab$tol), NA_real_, 1)
  zz <- ifelse(tab$mcse > 0, (abs(dev) - tab$tol) / tab$mcse, Inf)
  i <- out & !tab$metric %in% c("bias_std", PROP_METRICS)
  p[i] <- pnorm(-zz[i])
  i <- out & tab$metric == "bias_std" & tab$R >= 2
  p[i] <- pt(-zz[i], tab$R[i] - 1)
  i <- which(out & tab$metric %in% PROP_METRICS)
  x <- round(tab$value[i] * tab$R[i])
  up <- dev[i] > 0
  p[i] <- ifelse(up, pbinom(x - 1, tab$R[i], pmin(tab$target[i] + tab$tol[i], 1),
                            lower.tail = FALSE),
                 pbinom(x, tab$R[i], pmax(tab$target[i] - tab$tol[i], 0)))
  oracle_verdict(tab, p, alpha, TRUE)
}
# Holm over each family's gated rows with a p, then the verdicts.
# An sbc_ecdf row keeps its own p (the exact ECDF test's) and fails on
# p_holm < alpha alone, under both rules.
oracle_verdict <- function(tab, p, alpha, beyond_tol) {
  sbc <- tab$metric == "sbc_ecdf"
  p[sbc] <- tab$p[sbc]
  ph <- rep(NA_real_, nrow(tab))
  for (f in unique(tab$family)) {
    i <- which(tab$family == f & tab$gated & !is.na(p))
    ph[i] <- p.adjust(p[i], "holm")
  }
  na <- is.na(ph) | (!sbc & (is.na(tab$value) | is.na(tab$mcse)))
  verdict <- ifelse(!tab$gated, "finding", ifelse(na, "na",
                    ifelse(ph < alpha & (sbc | beyond_tol), "fail", "pass")))
  list(p = p, p_holm = ph, verdict = verdict)
}
same <- function(a, b, tol = 1e-9) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) &
                                                           abs(a - b) <= tol * pmax(1, abs(b)))

section("0", "formulas, gating and the gate rule on a fixed record set", {
  recs <- fixed_records()
  rv <- do.call(new_rep_record, recs[[1]][names(formals(new_rep_record))])
  sk <- SK_A
  re <- fill_records(sk, matrix(A_THETA, 1), list(hessian = matrix(A_SD, 1)))[[1]]
  rv2 <- do.call(new_rep_record, re[names(formals(new_rep_record))])
  report("0.1 test records validate: new_rep_record() rebuilds them unchanged",
         identical(rv, recs[[1]]) && identical(rv2, re))

  report("0.2 tolerances(): the API's mc and quick sets, with `failure` (the lead's)",
         identical(tolerances("mc"), c(bias_std = 0.10, se_ratio = 0.05, tail = 0.0075,
                                       coverage = 0.015, size = 0.015, failure = 0.005)) &&
           identical(tolerances("quick"), c(bias_std = 0.30, se_ratio = 0.15, tail = 0.025,
                                            coverage = 0.05, size = 0.025, failure = 0.02)))

  gates <- list(cell = "A", metrics = list(
    m1 = c("bias_std@beta", "se_ratio:hessian", "tail_lower:robust", "tail_upper:hessian@wtp",
           "joint_wald_size", "lr_size"),
    m2 = "bias_std", m3 = "bias_std@derived"))
  tab <- metrics_cell(recs, gates, "mc", LEVEL)
  cols <- c("family", "cell", "method", "metric", "estimand", "route", "value", "mcse", "target",
            "tol", "p", "gated", "verdict", "R", "n_failed")
  miss <- setdiff(cols, names(tab))
  orc <- oracle_cell(recs)
  ko <- row_key(orc$metric, orc$estimand, orc$route)
  kt <- row_key(tab$metric, tab$estimand, tab$route)
  i <- match(ko, kt)
  is_prop <- orc$metric %in% PROP_METRICS
  at_est <- same(tab$mcse[i], orc$mcse)
  at_tgt <- same(tab$mcse[i], orc$mcse_t)
  at_smo <- same(tab$mcse[i], orc$mcse_s)
  conv <- if (all(at_est[is_prop])) "the estimate (the API's reading)" else
    if (all(at_tgt[is_prop])) "the target" else
      if (all(at_smo[is_prop])) "the smoothed estimate (x + 1) / (R + 2)" else "no one rule"
  mcse_ok <- ifelse(orc$metric %in% c("failure_rate", "global_max_rate"), TRUE,
                    ifelse(is_prop, at_est | at_tgt | at_smo, at_est))
  bad <- is.na(i) | !same(tab$value[i], orc$value) | !mcse_ok | !same(tab$target[i], orc$target) |
    !same(tab$tol[i], orc$tol)
  extra <- setdiff(kt, ko)
  cell_rows <- tab$metric %in% c("failure_rate", "global_max_rate")
  report("0.3 metrics_cell(): columns, rows, value, MCSE, target, tol, R, n_failed",
         !length(miss) && !any(bad) && !length(extra) && all(tab$n_failed == 4) &&
           all(tab$R[!cell_rows] == 40),
         sprintf("%d rows expected, %d mismatched, %d unexpected; missing columns: %s\n",
                 length(ko), sum(bad), length(extra), if (length(miss)) toString(miss) else "none"),
         if (any(bad)) paste0("  mismatch ", ko[bad], ": got ", signif(tab$value[i[bad]], 6), " / ",
                              signif(tab$mcse[i[bad]], 6), ", expected ", signif(orc$value[bad], 6),
                              " / ", signif(orc$mcse[bad], 6), " (MCSE at the target ",
                              signif(orc$mcse_t[bad], 6), ", smoothed ",
                              signif(orc$mcse_s[bad], 6), ")\n", collapse = ""),
         if (length(extra)) paste0("  unexpected ", extra, "\n", collapse = ""),
         sprintf("R: %s for metrics (claim 40 = ok reps), %s for failure_rate, %s for %s",
                 toString(unique(tab$R[!cell_rows])), toString(tab$R[tab$metric == "failure_rate"]),
                 toString(tab$R[tab$metric == "global_max_rate"]), "global_max_rate"),
         sprintf("; n_failed %s; proportion MCSEs evaluated at %s", toString(unique(tab$n_failed)),
                 conv))

  want <- c(row_key("bias_std", c("b1", "b2"), NA), row_key("se_ratio", c("b1", "b2"), "hessian"),
            row_key("tail_lower", c("b1", "b2"), "robust"), row_key("tail_upper", "wtp", "hessian"),
            row_key("joint_wald_size", NA, NA), row_key("lr_size", NA, NA),
            # auto-gated: the routes of gated SE metrics, the gated Wald's failures
            row_key("route_failure_rate", NA, c("hessian", "robust")),
            row_key("joint_wald_failure_rate", NA, NA))
  got <- kt[tab$gated]
  relabel <- function(rs, method = NULL, cell = NULL) lapply(rs, function(r) {
    if (!is.null(method)) r$method <- method
    if (!is.null(cell)) r$cell <- cell
    r
  })
  g2 <- metrics_cell(relabel(recs, method = "m2"), gates)
  g3 <- metrics_cell(relabel(recs, method = "m3"), gates)
  gz <- metrics_cell(relabel(recs, cell = "Z"), gates)
  gn <- metrics_cell(recs, NULL)
  k2 <- row_key(g2$metric, g2$estimand, g2$route)[g2$gated]
  k3 <- row_key(g3$metric, g3$estimand, g3$route)[g3$gated]
  report("0.4 gating: metric[:route][@block|@name], per method, est$gated, cells",
         setequal(got, want) && setequal(k2, row_key("bias_std", c("b1", "b2"), NA)) &&
           setequal(k3, row_key("bias_std", "wtp", NA)) && !any(gz$gated) && !any(gn$gated) &&
           all(is.na(tab$verdict)),
         "m1 gated: ", toString(got), "\nm2 gated: ", toString(k2), "\nm3 gated: ", toString(k3),
         sprintf("\ncell not gated: %d gated rows; gates NULL: %d", sum(gz$gated), sum(gn$gated)))

  wl <- metrics_cell(fixed_records(wald_list = TRUE), gates)
  wr <- wl[grepl("joint_wald_size", wl$metric), ]
  wl_txt <- paste(wr$metric, wr$estimand, sep = "/")
  w_ok <- nrow(wr) == 2L && all(grepl("common", wl_txt) | grepl("beta", wl_txt)) &&
    any(grepl("common", wl_txt)) && any(grepl("beta", wl_txt))
  report("0.5 a named list of joint Wald tests gives one row per test", w_ok,
         "rows: ", toString(paste0(wl_txt, " = ", f4(wr$value), if (nrow(wr)) " gated " else "",
                                   wr$gated)))

  # Two families, each with a fail, an NA, a small deviation with a tiny
  # MCSE (pass by tolerance) and a large deviation within noise.
  tx <- rbind(transform(tab, family = "Fa"), transform(tab, family = "Fb"))
  ga <- which(tx$family == "Fa" & tx$gated)
  gb <- which(tx$family == "Fb" & tx$gated)
  tx$value[ga[1]] <- tx$target[ga[1]] + 0.5            # bias_std 0.5: a fail
  tx$mcse[ga[1]] <- 0.05
  tx$value[ga[2]] <- NA                                # na
  tx$mcse[ga[3]] <- 1e-5                               # deviation within tol: pass
  tx$value[ga[3]] <- tx$target[ga[3]] + 0.5 * tx$tol[ga[3]]
  tx$mcse[gb[1]] <- NA                                 # na
  tx$value[gb[2]] <- tx$target[gb[2]] + 3 * tx$tol[gb[2]]   # beyond tol, p within noise
  tx$mcse[gb[2]] <- 10 * tx$tol[gb[2]]
  # The plan rule's p: t_{R-1} for bias_std and exact binomial for
  # proportions (the lead accepted these small-R references).
  for (rule in RULES) for (alpha in c(0.01, 0.2)) {
    g <- gate_family(tx, alpha, rule = rule)
    o <- if (rule == "plan") oracle_gate(tx, alpha, exact = TRUE) else oracle_band(tx, alpha)
    gi <- tx$gated & !is.na(o$p_holm)
    bad <- which(!same(g$p, o$p) | !(g$verdict %in% o$verdict) |
                   (gi & !same(g$p_holm, o$p_holm)))
    report(sprintf("0.6 gate_family(rule = \"%s\", alpha = %s): p, Holm per family, verdicts",
                   rule, alpha), !length(bad),
           "verdicts: ", paste(names(table(g$verdict)), table(g$verdict), collapse = ", "),
           if (length(bad)) paste0("\nrows unlike the rule: ",
                                   toString(row_key(tx$metric, tx$estimand, tx$route)[bad])))
  }
  ga1 <- gate_family(tx[tx$family == "Fa", ], 0.01)
  report("0.7 Holm is within a family: Fa alone gives Fa's p_holm in the stack",
         all(same(ga1$p_holm, gate_family(tx, 0.01)$p_holm[tx$family == "Fa"])))

  # pooled_tails() on 30 reps of 40 constants, with and without groups.
  J <- 40L
  sd <- seq(0.1, 2, length.out = J)
  th <- seq(-1, 1, length.out = J)
  nm <- sprintf("ASC_%02d", seq_len(J))
  sk <- skeleton("P", c("x1", nm), c("beta", rep("asc", J)), c(0, th))
  D <- matrix(rnorm(30 * J, sd = 1.3), 30) * rep(sd, each = 30) + rep(th, each = 30)
  pr <- fill_records(sk, cbind(rnorm(30), D), list(hessian = matrix(c(1, sd), 30, J + 1,
                                                                    byrow = TRUE)))
  grp <- stats::setNames(rep(c("g1", "g2", "g3", "g4"), each = 10), nm)
  lo <- rep(th, each = 30) < D - ZQ * rep(sd, each = 30)
  up <- rep(th, each = 30) > D + ZQ * rep(sd, each = 30)
  direct <- function(hit, cols) {
    s <- rowMeans(hit[, cols, drop = FALSE])
    c(value = mean(s), mcse = sd(s) / sqrt(30))
  }
  pt <- pooled_tails(pr, block = "asc", level = LEVEL)
  pg <- pooled_tails(pr, block = "asc", groups = grp, level = LEVEL)
  chk <- function(tab, metric, hit, cols, group = NULL) {
    r <- tab[pooled_rows(tab, metric, group), ]
    d <- direct(hit, cols)
    nrow(r) == 1L && same(r$value, d[["value"]]) && same(r$mcse, d[["mcse"]])
  }
  ok <- chk(pt, "tail_lower", lo, seq_len(J)) && chk(pt, "tail_upper", up, seq_len(J)) &&
    all(vapply(unique(grp), function(gg) {
      cols <- which(grp == gg)
      chk(pg, "tail_lower", lo, cols, gg) && chk(pg, "tail_upper", up, cols, gg)
    }, NA))
  report("0.8 pooled_tails(): per-rep shares averaged, between-rep MCSE, per group", ok,
         "metrics: ", toString(unique(pt$metric)), "; estimand: ", toString(unique(pt$estimand)),
         "; rows: ", nrow(pt), " (", nrow(pg), " with 4 groups, estimand ",
         toString(unique(pg$estimand)), ")")

  tl <- skew_right(A_SKEW, ZQ)
  tu <- skew_left(A_SKEW, ZQ)
  x <- vs_with_seed(5, {
    z <- rnorm(2e6)
    (z + A_SKEW * (z^2 - 1)) / sqrt(1 + 2 * A_SKEW^2)
  })
  report("0.9 the tail04 deviant: the skewed error's tails (exact) and moments",
         abs(tl - 0.04) < 1e-8 && abs(mean(x > ZQ) - tl) < 4 * sqrt(tl / 2e6) &&
           abs(mean(x < -ZQ) - tu) < 4 * sqrt(tu / 2e6) && abs(mean(x)) < 0.003 &&
           abs(sd(x) - 1) < 0.003,
         sprintf(paste0("a = %.5f; tail_lower %.5f (simulated %.5f), tail_upper %.5f (simulated",
                        " %.5f); mean %.4f, sd %.4f (2e6 draws)"), A_SKEW, tl, mean(x > ZQ), tu,
                 mean(x < -ZQ), mean(x), sd(x)))

  report("0.10 cell A's joint Wald (rowSums(E^2)) equals fit.R's joint_wald() on the true V", {
    E <- matrix(rnorm(18), 3)
    est <- rep(A_THETA, each = 3) + (E %*% A_CHOL) * rep(A_SD, each = 3)
    V <- diag(A_SD) %*% A_CORR %*% diag(A_SD)
    dimnames(V) <- list(A_NAMES, A_NAMES)
    w <- vapply(1:3, function(r) joint_wald(stats::setNames(est[r, ], A_NAMES),
                                            stats::setNames(A_THETA, A_NAMES), V)$W, 0)
    all(same(w, rowSums(E^2), 1e-10))
  })

  # Route and joint Wald failures on 20 constructed reps (fail_records()):
  # se_hessian missing for b2 in reps 3 and 7 and for b1 in 7 and 12 (3
  # reps), se_robust for b1 in rep 5 (1 rep), bhhh and cluster never
  # computed; W NA in reps 2 and 9 and infinite in rep 15 (3 reps).
  fr <- fail_records(20L, se_h = list(b1 = c(7, 12), b2 = c(3, 7)), se_r = list(b1 = 5),
                     w_na = c(2, 9), w_inf = 15)
  ft <- metrics_cell(fr, NULL)
  smo <- function(x, n) sqrt((x + 1) / (n + 2) * (1 - (x + 1) / (n + 2)) / n)
  want <- data.frame(key = c(row_key("route_failure_rate", NA, c("hessian", "robust")),
                             row_key("joint_wald_failure_rate", NA, NA)),
                     value = c(3, 1, 3) / 20, mcse = smo(c(3, 1, 3), 20))
  i <- match(want$key, row_key(ft$metric, ft$estimand, ft$route))
  wsz <- ft[ft$metric == "joint_wald_size", ]
  ok <- !anyNA(i) && all(same(ft$value[i], want$value)) && all(same(ft$mcse[i], want$mcse)) &&
    all(ft$target[i] == 0) && all(ft$tol[i] == tolerances("mc")[["failure"]]) &&
    all(ft$R[i] == 20) && !any(ft$route %in% c("bhhh", "cluster")) && wsz$R == 17
  report("0.11 route_failure_rate and joint_wald_failure_rate: exact values", ok,
         "got ", toString(paste0(want$key, " = ", ft$value[i], " (mcse ", signif(ft$mcse[i], 4),
                                 ")")),
         sprintf("\njoint_wald_size over R = %d reps (claim: the 17 with a finite W; the ", wsz$R),
         "non-finite ones count only in joint_wald_failure_rate)")

  # Auto-gating: gating a route's SE metric gates its route_failure_rate (a
  # row added, value 1, where the route has no SE in any rep), gating the
  # joint Wald's size its failure rate.
  gk <- function(g) row_key(g$metric, g$estimand, g$route)[g$gated &
                    g$metric %in% c("route_failure_rate", "joint_wald_failure_rate")]
  g1 <- metrics_cell(fr, list(cell = "A", metrics = c("se_ratio:hessian", "joint_wald_size")))
  g2 <- metrics_cell(fr, list(cell = "A", metrics = "tail_lower:robust"))
  g3 <- metrics_cell(fr, list(cell = "A", metrics = "se_ratio:bhhh"))
  b3 <- g3[g3$metric == "route_failure_rate" & g3$route %in% "bhhh", ]
  v1 <- vapply(RULES, function(r) {
    g <- gate_family(g1, rule = r)
    toString(g$verdict[g$metric %in% c("route_failure_rate", "joint_wald_failure_rate") &
                         g$gated])
  }, "")
  report("0.12 auto-gating: an SE metric's route_failure_rate, the joint Wald's failure rate",
         setequal(gk(g1), c(row_key("route_failure_rate", NA, "hessian"),
                            row_key("joint_wald_failure_rate", NA, NA))) &&
           setequal(gk(g2), row_key("route_failure_rate", NA, "robust")) &&
           nrow(b3) == 1L && isTRUE(b3$gated) && isTRUE(b3$value == 1) &&
           all(v1 == "fail, fail"),
         "se_ratio:hessian + joint_wald_size gate: ", toString(gk(g1)),
         "\ntail_lower:robust gates: ", toString(gk(g2)),
         sprintf("\nse_ratio:bhhh (never computed) gates route_failure_rate|bhhh = %s",
                 toString(b3$value)),
         sprintf("\nverdicts of the 3/20 failure rates: plan %s; band %s", v1[["plan"]],
                 v1[["band"]]))

  g0 <- vapply(c(20L, 1000L), function(R) vapply(RULES, function(r) {
    g <- gate_family(metrics_cell(fail_records(R), list(cell = "A", metrics = c(
      "se_ratio:hessian", "tail_lower:robust", "joint_wald_size"))), rule = r)
    i <- g$gated & g$metric %in% c("route_failure_rate", "joint_wald_failure_rate")
    paste0(sum(i), ":", toString(unique(g$verdict[i])), ":", toString(unique(g$value[i])))
  }, ""), c(plan = "", band = ""))
  report("0.13 no route or Wald failures: the failure-rate rows pass under both rules",
         all(g0 == "3:pass:0"), "R = 20, 1000 (rows:verdicts:values): ", toString(g0))

  # Missed global maxima: dlogL > 0.01 and max |dtheta| / SE > 0.1, both
  # strict; a check without dlogL is left out; reps 10-12 ran none.
  gm_dl <- c(0.02, 0.02, 0.005, 1e-5, 5, 0.011, NA, 0.01, 0.5)
  gm_dt <- c(0.2, 0.05, 0.5, 1, 3, NA, 1, 0.5, 0.1)
  gmr <- lapply(seq_len(12), function(r) {
    x <- fail_records(1L)[[1]]
    x$rep <- r
    if (r <= 9) x$diag$multistart <- list(dlogL = gm_dl[r], max_dtheta_se = gm_dt[r],
                                          n_failed = 0L, n_warned = 0L)
    x
  })
  gt <- metrics_cell(gmr, NULL)
  gv <- gt[gt$metric == "global_max_rate", ]
  report("0.14 global_max_rate: dlogL > 0.01 and max |dtheta| / SE > 0.1 among checks run",
         nrow(gv) == 1L && same(gv$value, 2 / 8) && gv$R == 8 && same(gv$target, 0) &&
           !isTRUE(gv$gated),
         sprintf("value %s over R = %s (expected 2/8: (0.02, 0.2) and (5, 3); not the boundary ",
                 toString(gv$value), toString(gv$R)),
         "pairs (0.01, 0.5) and (0.5, 0.1), nor (0.02, 0.05), (0.011, NA); (NA, 1) left out)")
})

# --- 1-4. The null family at R = 100 and 1000 -------------------------------------

STUDY <- list()
run_null <- function() {
  if (!is.null(STUDY$r1000)) return(invisible())
  t0 <- proc.time()[["elapsed"]]
  vs_set_seed(opt$seed + 100L)
  STUDY$r100 <<- null_study(100, M$null100, "quick")
  t1 <- proc.time()[["elapsed"]]
  STUDY$r1000 <<- null_study(1000, M$null1000, "mc")
  cat(sprintf("        null studies: R = 100 x %d in %.1f s, R = 1000 x %d in %.1f s\n",
              M$null100, t1 - t0, M$null1000, proc.time()[["elapsed"]] - t1))
}

section("1", "MCSE calibration under the noisy null", {
  run_null()
  for (s in STUDY[c("r100", "r1000")]) {
    tv <- true_values(s$cal, s$R)
    z <- (s$cal$value - tv$expected) / s$cal$mcse
    zp <- (s$cal$value - tv$population) / s$cal$mcse
    for (met in CAL_METRICS) {
      i <- s$cal$metric == met
      zi <- z[i]
      fin <- is.finite(zi)
      b <- band(s$M)
      sdz <- sd(zi[fin])
      mz <- mean(zi[fin])
      # The brief's seven metrics, and coverage and lr_size, which a family
      # may gate. At R = 100 the binomial rows' and lr_mean_q's mean(z) is
      # reported, not judged (the lead, 2026-11-02).
      why <- if (s$R < 1000 && met %in% PROP_METRICS) {
        paste("a plug-in MCSE about the estimate is biased at expected counts of 2.5 to 5;",
              "the gate tests counts exactly and does not read it")
      } else if (s$R < 1000 && met == "lr_mean_q") {
        "LR / q is skewed (chi^2_2 / 2, skewness 2); lr_mean_q is ungated"
      }
      ok <- all(fin) && abs(sdz - 1) <= b[["sd"]] && (!is.null(why) || abs(mz) <= b[["mean"]])
      ev <- unique(signif(tv$expected[i], 6))
      extra <- ""
      if (met %in% PROP_METRICS) {
        p0 <- tv$expected[i]
        zs <- (s$cal$value[i] - p0) / sqrt(p0 * (1 - p0) / s$R)
        extra <- sprintf("\nwith the MCSE at the true p: sd %s, mean %s; MCSE 0 in %d rows",
                         f3(sd(zs)), f3(mean(zs)), sum(s$cal$mcse[i] == 0, na.rm = TRUE))
      }
      if (met %in% c("ese", "se_ratio")) {
        extra <- sprintf("\nmean z against the population value (%s): %s",
                         if (met == "ese") "sigma" else "1", f3(mean(zp[i][is.finite(zp[i])])))
      }
      report(sprintf("1 R=%-4d %-15s sd(z) %s  mean(z) %s%s", s$R, met, f3(sdz), f3(mz),
                     if (is.null(why)) "" else " (INFO)"), ok,
             sprintf("n %d z-scores (%d meta-reps), %d not finite; bands: |sd - 1| <= %s, ",
                     length(zi), s$M, sum(!fin), f3(b[["sd"]])),
             if (is.null(why)) sprintf("|mean| <= %s", f3(b[["mean"]])) else
               paste0("mean(z) not judged: ", why),
             sprintf("\ntrue value E[value] at R (exact): %s", toString(head(ev, 4))),
             if (length(ev) > 4) " ..." else "", extra)
    }
  }
  cat("        true values: bias_std 0 (symmetry); ese sigma c4(R); se_ratio",
      "E sqrt(chi2_{R nu}/(R nu)) E[sigma/s] (independent);\n        tails",
      "P(t_nu > z) (Z / S ~ t_nu; nu = 50 hessian, 30 robust); coverage 1 - 2 tails; sizes 0.05;",
      "lr_mean_q 1.\n")
})

section("2", "gate size: null families of 33 gated rows, both rules", {
  run_null()
  for (s in STUDY[c("r1000", "r100")]) {
    ts <- if (s$R == 1000) "mc" else "quick"
    tol <- tolerances(ts)
    for (rule in RULES) {
      f <- s$fam[[rule]][, "null"]
      nf <- sort(table(s$null_fails[[rule]]), decreasing = TRUE)
      report(sprintf("2 R=%d, %s tolerances, rule %s: exact null families failing %s", s$R, ts,
                     rule, rate_se(f)),
             s$n_gated[["null"]] == N_GATED && mean(f) <= 0.02,
             sprintf("gated rows: %d; SEs exact, every true value at its target; ",
                     s$n_gated[["null"]]), "claim: share <= 0.02",
             if (length(nf)) paste0("\nfailed rows (count): ", toString(head(paste0(
               names(nf), " (", nf, ")"), 8))))
    }
    # The noisy null: each Hessian tail is P(t_50 > z) = 0.0278, inside the band.
    tail <- pt(ZQ, NU, lower.tail = FALSE)
    q <- vapply(RULES, function(r) 1 - (1 - power_prop(tail, A2, tol[["tail"]], s$R, r))^14, 0)
    report(sprintf("2 R=%d, %s tolerances: noisy null (nu = 50) families failing: plan %s, band %s",
                   s$R, ts, f4(mean(s$fam$plan[, "noisy"])), f4(mean(s$fam$band[, "noisy"]))), NA,
           sprintf("n %d, gated rows %d; expected from the 14 Hessian tail rows alone (true tail ",
                   s$M, s$n_gated[["noisy"]]),
           sprintf("%s): plan %s, band %s", f4(tail), f4(q[["plan"]]), f4(q[["band"]])))
  }
})

section("3", "gate power against the named deviations (R = 1000, mc tolerances)", {
  run_null()
  s <- STUDY$r1000
  R <- s$R
  tol <- tolerances("mc")
  # The design's power (exact SEs; see power_prop() and friends), Holm over
  # each family's gated rows (s$n_gated: N_GATED, wald08's N_GATED + 2); in
  # tail04's, its tail_upper (0.0035 against 0.025) is rejected first in
  # nearly every run, which leaves one fewer.
  kf <- s$n_gated
  kf[["tail04"]] <- kf[["tail04"]] - 1L
  fn <- list(
    se90 = function(R_, rule, k) power_ratio(0.90, tol[["se_ratio"]], R_, rule, k),
    se97 = function(R_, rule, k) power_ratio(0.97, tol[["se_ratio"]], R_, rule, k),
    bias20 = function(R_, rule, k) power_bias(0.20, tol[["bias_std"]], R_, rule, k),
    bias05 = function(R_, rule, k) power_bias(0.05, tol[["bias_std"]], R_, rule, k),
    tail04 = function(R_, rule, k) power_prop(0.04, A2, tol[["tail"]], R_, rule, k),
    wald08 = function(R_, rule, k) power_prop(0.08, 0.05, tol[["size"]], R_, rule, k))
  outside <- c("se90", "bias20", "tail04", "wald08")
  grid <- c(seq(20, 1000, by = 10), seq(1050, 20000, by = 50))
  for (d in names(fn)) {
    pw <- list()
    pw[[d]] <- function(R_, rule, k = kf[[d]]) fn[[d]](R_, rule, k)
    for (rule in RULES) {
      caught <- s$dev[[rule]][, d]
      des <- pw[[d]](R, rule)
      se <- max(sqrt(des * (1 - des) / s$M), 1 / s$M)
      report(sprintf("3 %-6s rule %s: its row fails in %s; design %s", d, rule, rate_se(caught),
                     f4(des)), abs(mean(caught) - des) <= 3 * se,
             sprintf("claim: |measured - design| <= 3 MCSE (%s); family fails in %s",
                     f4(3 * se), f4(mean(s$fam[[rule]][, d]))))
    }
    need <- vapply(RULES, function(rule) {
      if (!d %in% outside) return("-")
      hit <- which(vapply(grid, function(r) pw[[d]](r, rule) >= 0.9, NA))
      if (length(hit)) format(grid[hit[1]], big.mark = ",") else "> 20,000"
    }, "")
    report(sprintf("3 %-6s design power at R = 1000: plan %s (alone %s), band %s (alone %s)", d,
                   f3(pw[[d]](R, "plan")), f3(pw[[d]](R, "plan", 1)), f3(pw[[d]](R, "band")),
                   f3(pw[[d]](R, "band", 1))), NA,
           if (d %in% outside) sprintf("first R with 90%% in its family: plan %s, band %s",
                                       need[["plan"]], need[["band"]]) else
             "inside the tolerance: power falls to 0 with R under both rules")
  }
  cat(sprintf(paste0("        deviant truths (exact SEs): se90 tails %s each, se97 %s; tail04",
                     " tail_lower 0.04, tail_upper %s (a = %.4f)\n"), f4(pnorm(-0.90 * ZQ)),
              f4(pnorm(-0.97 * ZQ)), f4(skew_left(A_SKEW, ZQ)), A_SKEW))
})

section("4", "joint Wald and LR under the null (V-hat = V, LR ~ chi^2_2)", {
  run_null()
  for (s in STUDY[c("r100", "r1000")]) {
    for (met in c("joint_wald_size", "lr_size", "lr_mean_q")) {
      v <- s$cal$value[s$cal$metric == met]
      target <- if (met == "lr_mean_q") 1 else 0.05
      se <- sd(v) / sqrt(length(v))
      report(sprintf("4 R=%-4d %-15s mean %s (SE %s) over %d meta-reps", s$R, met, f4(mean(v)),
                     f4(se), length(v)), abs(mean(v) - target) <= 3 * se,
             sprintf("claim: |mean - %s| <= 3 SE", target))
    }
  }
})

# --- 5. Pooled constant coverage -------------------------------------------------

section("5", "coverage pooled over 2,000 constants with a shared shift per rep", {
  J <- 2000L
  R <- 100L
  rho <- 0.1
  sd <- exp(rnorm(J, log(0.5), 0.5))
  th <- rnorm(J)
  nm <- sprintf("ASC_%04d", seq_len(J))
  sk <- skeleton("P", c("x1", nm), c("beta", rep("asc", J)), c(1, th))
  se <- matrix(c(0.1, sd), R, J + 1, byrow = TRUE)
  out <- matrix(NA_real_, M$pooled, 4, dimnames = list(NULL, c("lo", "lo_mcse", "up", "up_mcse")))
  for (m in seq_len(M$pooled)) {
    D <- (sqrt(rho) * rnorm(R) + sqrt(1 - rho) * matrix(rnorm(R * J), R)) * rep(sd, each = R)
    recs <- fill_records(sk, cbind(1 + 0.1 * rnorm(R), D + rep(th, each = R)), list(hessian = se))
    pt <- pooled_tails(recs, block = "asc", level = LEVEL)
    lo <- pt[pooled_rows(pt, "tail_lower"), ]
    up <- pt[pooled_rows(pt, "tail_upper"), ]
    out[m, ] <- c(lo$value, lo$mcse, up$value, up$mcse)
  }
  for (side in c("lo", "up")) {
    v <- out[, side]
    z <- (v - A2) / out[, paste0(side, "_mcse")]
    zn <- (v - A2) / sqrt(v * (1 - v) / (R * J))
    b <- band(M$pooled)
    report(sprintf("5 tail_%s pooled: sd(z) %s, mean(z) %s; naive binomial MCSE: sd(z) %s",
                   if (side == "lo") "lower" else "upper", f3(sd(z)), f3(mean(z)), f3(sd(zn))),
           abs(sd(z) - 1) <= b[["sd"]] && abs(mean(z) - 0) <= b[["mean"]] && sd(zn) > 1.5,
           sprintf(paste0("%d meta-reps of R = %d reps x %d constants, intra-rep correlation %s; ",
                          "true tail %s (exact: each z ~ N(0, 1));\nmean value %s, mean MCSE %s,",
                          " naive MCSE %s; bands |sd - 1| <= %s, |mean| <= %s, naive sd > 1.5"),
                   M$pooled, R, J, rho, A2, f4(mean(v)), f4(mean(out[, paste0(side, "_mcse")])),
                   f4(mean(sqrt(v * (1 - v) / (R * J)))), f3(b[["sd"]]), f3(b[["mean"]])))
  }
})

# --- 6. SBC -------------------------------------------------------------------

# The ECDF test's statistic from rank counts ((L + 1) x n; counts[k, ] =
# #{rank = k - 1}), computed here: the smallest over z_i = i / (L + 1) of the
# two-sided binomial tail of #{rank < i} under Binomial(R, z_i), capped at 1.
ecdf_stat <- function(counts, R, L) {
  cum <- 0
  stat <- rep(Inf, ncol(counts))
  for (i in seq_len(L)) {
    cum <- cum + counts[i, ]
    z <- i / (L + 1)
    stat <- pmin(stat, pmin(1, 2 * pmin(pbinom(cum, R, z),
                                        pbinom(cum - 1, R, z, lower.tail = FALSE))))
  }
  stat
}
# Whether each row of `ranks` leaves its test's band somewhere (F(z_i) =
# #{rank < i} / R; a point without an acceptance set counts as outside).
out_of_band <- function(ranks, res, R) {
  vapply(seq_len(nrow(ranks)), function(m) {
    Fh <- cumsum(tabulate(ranks[m, ] + 1L, L_SBC + 1L))[seq_len(L_SBC)] / R
    lo <- res[[m]]$lower
    hi <- res[[m]]$upper
    any(is.na(lo) | Fh < lo | Fh > hi)
  }, NA)
}
# n SBC ranks among L_SBC draws from a posterior whose sd is k times the
# truth's and whose mean is shifted by `shift` of its sd: Binomial(L,
# Phi((z - shift) / k)), z ~ N(0, 1).
misfit_ranks <- function(n, k, shift) {
  rbinom(n, L_SBC, pnorm((rnorm(n) - shift) / k))
}
# One run of the SBC family: R records whose diag$sbc holds 30 null ranks
# (q01-q30) and, with `narrow`, a too-narrow quantity; metrics_cell() gated
# on sbc_ecdf, then gate_family() over the families "null" (q01-q30) and
# "narrow" (q02-q30 and narrow).
SK_S <- new_rep_record("S", "S", "m", rep = 1, seed = NA_real_)
sbc_family <- function(R, narrow) {
  rk <- matrix(sample.int(L_SBC + 1L, R * 30, replace = TRUE) - 1L, R,
               dimnames = list(NULL, sprintf("q%02d", 1:30)))
  if (narrow) rk <- cbind(rk, narrow = misfit_ranks(R, 0.8, 0))
  recs <- lapply(seq_len(R), function(r) {
    x <- SK_S
    x$rep <- r
    x$diag$sbc <- list(L = L_SBC, ranks = rk[r, ])
    x
  })
  tab <- metrics_cell(recs, list(cell = "S", metrics = "sbc_ecdf"))
  fam <- transform(tab[!tab$estimand %in% "narrow", ], family = "null")
  if (narrow) fam <- rbind(fam, transform(tab[!tab$estimand %in% "q01", ], family = "narrow"))
  gate_family(fam)
}

section("6", "SBC: sbc_rank(), the exact ECDF test (L = 99) and SBC families", {
  r1 <- sbc_rank(0.5, c(0.1, 0.6, 0.4))
  tie <- vs_with_seed(3, vapply(1:4000, function(i) sbc_rank(1, c(0, 1, 1, 1)), 0))
  tt <- tabulate(tie + 1L, 5L)
  report("6.0 sbc_rank(): count of draws below; ties split uniformly",
         r1 == 2 && is.numeric(r1) && sbc_rank(-Inf, 1:5) == 0 && sbc_rank(Inf, 1:5) == 5 &&
           tt[1] == 0 && sum(tt[2:5]) == 4000 && all(abs(tt[2:5] / 4000 - 0.25) < 0.03),
         "ranks of 1 among (0, 1, 1, 1): ", toString(tt[2:5]), " at 1..4")

  probe <- sample.int(L_SBC + 1L, 1000, replace = TRUE) - 1L
  seed_before <- .Random.seed
  t0 <- proc.time()[["elapsed"]]
  tst <- sbc_ecdf_test(probe, L_SBC)
  t1 <- proc.time()[["elapsed"]] - t0
  same_rng <- identical(seed_before, .Random.seed)
  report("6.1 sbc_ecdf_test() fields, timing and the caller's RNG state (R = 1000)",
         all(c("stat", "p", "gamma", "z", "lower", "upper", "R") %in% names(tst)) &&
           length(tst$z) == L_SBC && length(tst$lower) == L_SBC &&
           length(tst$upper) == L_SBC && same_rng,
         sprintf("%.3f s for the first call (the band included); caller's RNG state %s", t1,
                 if (same_rng) "restored" else "CHANGED by the call"))

  # The exact p against a simulation of the discrete null at small (R, L):
  # the statistic recomputed here from the rank counts.
  n_sim <- if (opt$fast) 5e4 else 2e5
  for (RL in list(c(20, 9), c(60, 19))) {
    R <- RL[1]
    L <- RL[2]
    counts <- rmultinom(n_sim, R, rep(1, L + 1))
    st <- ecdf_stat(counts, R, L)
    pts <- vapply(c(0.002, 0.03, 0.3, 0.8), function(q) which.min(abs(st - quantile(st, q))), 1L)
    res <- vapply(pts, function(j) {
      t <- sbc_ecdf_test(rep(0:L, counts[, j]), L)
      ps <- mean(st <= st[j] * (1 + 1e-12))
      c(p = t$p, sim = ps, z = (t$p - ps) / sqrt(ps * (1 - ps) / n_sim),
        same_stat = same(t$stat, st[j], 1e-12))
    }, numeric(4))
    report(sprintf("6.2 exact p against %s simulated nulls, R = %d, L = %d: max |z| %s", n_sim, R,
                   L, f3(max(abs(res["z", ])))),
           all(abs(res["z", ]) <= 3) && all(res["same_stat", ] == 1),
           "p exact / simulated: ", toString(sprintf("%.5f / %.5f", res["p", ], res["sim", ])),
           "; statistic recomputed here equal: ", all(res["same_stat", ] == 1),
           "\nthe simulated p counts statistics within 1e-12 of the observed one as ties")
  }
  # The statistic is symmetric under ranks -> L - ranks (the tail at z_i of
  # count c equals the tail at z_{L+1-i} of R - c), so the p-values of a rank
  # vector and its mirror image are equal.
  for (RL in list(c(20, 9), c(100, 99), c(1000, 99))) {
    d <- vapply(seq_len(if (opt$fast) 40 else 100), function(k) {
      r <- sample.int(RL[2] + 1L, RL[1], replace = TRUE) - 1L
      c(sbc_ecdf_test(r, RL[2])$p, sbc_ecdf_test(RL[2] - r, RL[2])$p)
    }, numeric(2))
    dp <- abs(d[1, ] - d[2, ])
    report(sprintf("6.2 mirror symmetry, R = %d, L = %d: p differs in %d of %d pairs, max %s",
                   RL[1], RL[2], sum(dp > 1e-9), ncol(d), formatC(max(dp), digits = 3,
                                                                  format = "g")),
           all(dp <= 1e-9), "claim: p(ranks) = p(L - ranks) within 1e-9")
  }

  # The null p uniform; the band left exactly where p <= 0.05. Ranks drawn
  # before any call.
  agree <- logical()
  for (R in c(100L, 1000L)) {
    ranks <- matrix(sample.int(L_SBC + 1L, M$sbc * R, replace = TRUE) - 1L, M$sbc)
    res <- lapply(seq_len(M$sbc), function(m) sbc_ecdf_test(ranks[m, ], L_SBC))
    p <- vapply(res, `[[`, 0, "p")
    agree <- c(agree, out_of_band(ranks, res, R) == (p <= 0.05))
    a <- c(0.05, 0.01, 0.001)
    rate <- vapply(a, function(x) mean(p <= x), 0)
    se <- sqrt(a * (1 - a) / M$sbc)
    report(sprintf("6.3 null, R=%d: rejection at 0.05 %s, 0.01 %s, 0.001 %s; mean p %s", R,
                   f4(rate[1]), f4(rate[2]), f4(rate[3]), f3(mean(p))),
           all(abs(rate - a) <= 3 * se),
           sprintf("%d meta-reps; claim: within 3 binomial SEs (%s); ", M$sbc,
                   toString(f4(3 * se))),
           sprintf("p deciles: %s", toString(f3(quantile(p, 1:9 / 10)))))
  }
  alts <- list(narrow = c(k = 0.8, shift = 0), wide = c(k = 1.25, shift = 0),
               shifted = c(k = 1, shift = 0.2))
  for (a in names(alts)) {
    pw <- vapply(c(100L, 1000L), function(R) {
      ranks <- misfit_ranks(M$sbc_power * R, alts[[a]][["k"]], alts[[a]][["shift"]])
      ranks <- matrix(ranks, M$sbc_power)
      res <- lapply(seq_len(M$sbc_power), function(m) sbc_ecdf_test(ranks[m, ], L_SBC))
      p <- vapply(res, `[[`, 0, "p")
      agree <<- c(agree, out_of_band(ranks, res, R) == (p <= 0.05))
      mean(p <= 0.05)
    }, 0)
    report(sprintf("6.4 power at 0.05 against a %s posterior: R=100 %s, R=1000 %s", a, f3(pw[1]),
                   f3(pw[2])), NA,
           sprintf("posterior sd %s x the truth's, shifted by %s posterior sd; %d meta-reps each",
                   alts[[a]][["k"]], alts[[a]][["shift"]], M$sbc_power))
  }
  report(sprintf("6.5 the ECDF leaves the 0.05 band exactly when p <= 0.05: %d of %d runs agree",
                 sum(agree), length(agree)), all(agree), "null and misfit runs above")

  sev <- sbc_ecdf_test(misfit_ranks(1000, 0.5, 0), L_SBC)
  report(sprintf("6.6 a severe misfit (posterior sd 0.5 x, R = 1000): p = %.3g", sev$p),
         is.finite(sev$p) && sev$p >= 0 && sev$p < 1e-6,
         "claim: p < 1e-6 (no Monte Carlo floor)")

  # Families of 30 SBC quantities under the null, gated on sbc_ecdf; at
  # R = 1000 a 31st quantity, too narrow (posterior sd 0.8 x), replaces one
  # of them in a second family.
  for (R in c(100L, 1000L)) {
    Ms <- M$sbc_fam[[as.character(R)]]
    fails <- matrix(NA, Ms, 2, dimnames = list(NULL, c("null", "narrow")))
    for (m in seq_len(Ms)) {
      g <- sbc_family(R, narrow = R == 1000L)
      if (m == 1L) {
        o <- oracle_gate(g, exact = TRUE)
        gb <- gate_family(g, rule = "band")
        ok1 <- identical(as.character(g$verdict), o$verdict) && all(same(g$p_holm, o$p_holm)) &&
          identical(g$verdict, gb$verdict)
        ng <- tapply(g$gated, g$family, sum)
      }
      f <- g$verdict %in% "fail"
      fails[m, ] <- c(any(f[g$family == "null"]),
                      any(f[g$family == "narrow" & g$estimand %in% "narrow"]))
    }
    report(sprintf("6.7 SBC family of %d quantities under the null, R=%d: failing %s", ng[["null"]],
                   R, rate_se(fails[, "null"])),
           ok1 && ng[["null"]] == 30 && mean(fails[, "null"]) <= 0.02,
           "Holm at alpha 0.01; claim: share <= 0.02; first run's verdicts and p_holm equal the ",
           "rule applied here (sbc_ecdf on p_holm < alpha alone) and the band rule's: ", ok1)
    if (R == 1000L) {
      report(sprintf("6.8 a too-narrow quantity (posterior sd 0.8 x) among 29 null, R=1000: %s",
                     paste("caught in", rate_se(fails[, "narrow"]))),
             mean(fails[, "narrow"]) >= 0.9, "claim: caught >= 0.90")
    }
  }
})

# --- 7. rate_fit() --------------------------------------------------------------

section("7", "rate_fit(): the slope of ESE in N", {
  N <- c(500, 2000, 8000)
  R <- c(400L, 200L, 100L)
  cells <- as.character(N)
  # A real metrics table: one cell per N, ESE = 2 / sqrt(N).
  tab <- do.call(rbind, lapply(seq_along(N), function(j) {
    sk <- skeleton(cells[j], "b1", "beta", 1, family = "RL")
    s <- 2 / sqrt(N[j])
    metrics_cell(fill_records(sk, matrix(1 + s * rnorm(R[j]), R[j]),
                              list(hessian = matrix(s, R[j], 1))), NULL)
  }))
  x <- stats::setNames(N, cells)
  rf0 <- rate_fit(tab, x, "ese", estimand = "b1")
  idx <- vapply(cells, function(cl) which(tab$cell == cl & rows_of(tab, "ese", "b1")), 1L)
  kap <- 1 / sqrt(2 * (R - 1))
  sl <- se <- numeric(M$rate)
  for (m in seq_len(M$rate)) {
    v <- 2 / sqrt(N) * (1 + kap * rnorm(3))
    tab$value[idx] <- v
    tab$mcse[idx] <- v * kap
    rf <- rate_fit(tab, x, "ese", estimand = "b1")
    sl[m] <- rf[["slope"]]
    se[m] <- rf[["se"]]
  }
  cover <- abs(sl + 0.5) <= 2 * se
  report(sprintf("7 slope within 2 SE of -0.5 in %s", rate_se(cover)), mean(cover) >= 0.90,
         sprintf("mean slope %s, sd %s, mean se %s (sd / se %s); ", f4(mean(sl)), f4(sd(sl)),
                 f4(mean(se)), f3(sd(sl) / mean(se))),
         sprintf("one fit on records: %s (se %s, n %s)", f4(rf0[["slope"]]), f4(rf0[["se"]]),
                 rf0[["n"]]),
         "\ncells N = 500, 2000, 8000 at R = 400, 200, 100; value noise at the MCSE ",
         "ESE / sqrt(2 (R - 1)); claim: >= 0.90")
})

# --- 8. Edge cases ----------------------------------------------------------------

G_EDGE <- list(cell = "E", metrics = c("bias_std", "se_ratio:hessian", "tail_lower:hessian",
                                       "tail_upper:hessian", "joint_wald_size", "lr_size"))
# A small cell E: b1, b2 with se_hessian; failed reps `fail`; NA truth of b2
# in reps `truth_na`; NA se_hessian of `se_na` estimands; NA est of b1 in
# reps `est_na`; NA Wald in reps `wald_na`; estimates fixed at the truth.
edge_records <- function(R, fail = integer(), truth_na = integer(), se_na = character(),
                         est_na = integer(), wald_na = integer(), est_const = FALSE) {
  lapply(seq_len(R), function(r) {
    if (r %in% fail) return(new_rep_record("W", "E", "m", r, NA_real_, ok = FALSE, error = "x"))
    th <- c(1, -1)
    sd <- c(0.3, 0.2)
    est <- data.frame(name = c("b1", "b2"), block = "beta",
                      est = if (est_const) th else th + rnorm(2, sd = sd), truth = th,
                      se_hessian = sd * sqrt(rchisq(2, 50) / 50), se_bhhh = NA_real_,
                      se_robust = NA_real_, se_cluster = NA_real_, gated = TRUE,
                      stringsAsFactors = FALSE)
    if (r %in% truth_na) est$truth[2] <- NA
    if (r %in% est_na) est$est[1] <- NA
    est$se_hessian[est$name %in% se_na] <- NA
    W <- if (r %in% wald_na) NA_real_ else rchisq(1, 2)
    LR <- rchisq(1, 1)
    new_rep_record("W", "E", "m", r, NA_real_, est = est,
                   joint_wald = list(W = W, df = 2, p = pchisq(W, 2, lower.tail = FALSE)),
                   lr = list(LR = LR, df = 1, p = pchisq(LR, 1, lower.tail = FALSE)))
  })
}
# The table and its violations of the verdict rule: a gated row with an NA
# value or MCSE not "na", an ungated row not "finding", any NA verdict.
edge_run <- function(recs, gates = G_EDGE) {
  g <- gate_family(metrics_cell(recs, gates))
  na <- is.na(g$value) | is.na(g$mcse)
  bad <- is.na(g$verdict) | (g$gated & na & g$verdict != "na") | (!g$gated & g$verdict != "finding")
  list(g = g, bad = which(bad %in% TRUE))
}
verdicts <- function(g) {
  v <- table(factor(g$verdict[g$gated], c("pass", "fail", "na")))
  paste(names(v), v, collapse = ", ")
}

section("8", "edge cases", {
  cases <- list(
    "R = 1" = list(recs = edge_records(1), all_na = TRUE),
    "R = 2" = list(recs = edge_records(2)),
    "every rep failed (5)" = list(recs = edge_records(5, fail = 1:5), all_na = TRUE),
    "se_hessian of b2 all NA (route gated)" = list(recs = edge_records(20, se_na = "b2"),
                                                   gated_na = c("b2", "hessian"),
                                                   fails = list("route_failure_rate", "hessian",
                                                                1)),
    "se_cluster all NA in the cell, se_ratio:cluster gated" = list(
      recs = edge_records(20), gated_na = c("b1", "cluster"),
      gates = list(cell = "E", metrics = c("se_ratio:cluster", "bias_std")),
      fails = list("route_failure_rate", "cluster", 1)),
    "se_cluster all NA in the cell, not gated" = list(recs = edge_records(20),
                                                      no_route = "cluster"),
    "truth of b2 NA in every rep" = list(recs = edge_records(20, truth_na = 1:20), b2_na = TRUE),
    "truth of b2 NA in 3 of 20 reps" = list(recs = edge_records(20, truth_na = c(2, 9, 14))),
    "est of b1 NA in 2 of 20 reps" = list(recs = edge_records(20, est_na = c(4, 11))),
    "joint Wald NA in 2 of 20 reps" = list(recs = edge_records(20, wald_na = c(3, 7)),
                                           fails = list("joint_wald_failure_rate", NA, 0.1)),
    "estimates equal to the truth (ESE 0)" = list(recs = edge_records(20, est_const = TRUE)))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    res <- tryCatch(edge_run(cs$recs, cs$gates %||% G_EDGE),
                    error = function(e) conditionMessage(e))
    if (is.character(res)) {
      report(paste("8", nm, ": metrics_cell() + gate_family()"), FALSE, "error: ", res)
      next
    }
    g <- res$g
    ok <- !length(res$bad)
    det <- sprintf("%d rows, %d gated (%s); rule violations: %d", nrow(g), sum(g$gated),
                   verdicts(g), length(res$bad))
    if (isTRUE(cs$all_na)) {
      ok <- ok && all(g$verdict[g$gated] == "na")
      det <- paste0(det, "; claim: every gated row \"na\"")
    }
    if (isTRUE(cs$b2_na)) {
      i <- g$gated & g$estimand %in% "b2"
      ok <- ok && all(g$verdict[i] == "na") && all(is.na(g$value[i]))
      det <- paste0(det, sprintf("; b2's %d gated rows NA and \"na\"", sum(i)))
    }
    if (nm == "every rep failed (5)") {
      fr <- g$value[rows_of(g, "failure_rate")]
      ok <- ok && length(fr) == 1L && fr == 1
      det <- paste0(det, sprintf("; failure_rate %s; gated rows: %s", toString(fr),
                                 toString(row_key(g$metric, g$estimand, g$route)[g$gated])))
    }
    if (!is.null(cs$gated_na)) {
      # A gated route whose SEs are all NA cannot be judged: the family must
      # show it, by the estimand's rows or a placeholder row (estimand NA),
      # NA and "na".
      i <- (g$estimand %in% cs$gated_na[1] | is.na(g$estimand)) & g$route %in% cs$gated_na[2] &
        g$gated & g$metric != "route_failure_rate"
      ok <- ok && sum(i) >= 1L && all(g$verdict[i] == "na")
      det <- paste0(det, sprintf("; gated %s rows of %s or a placeholder: %d%s; claim: \"na\"",
                                 cs$gated_na[2], cs$gated_na[1], sum(i),
                                 if (sum(i)) paste0(" (", toString(paste(g$metric[i], g$estimand[i],
                                                                         g$verdict[i])), ")") else
                                   " (dropped: the gate cannot see them)"))
    }
    if (!is.null(cs$no_route)) {
      i <- g$route %in% cs$no_route
      ok <- ok && !any(i)
      det <- paste0(det, sprintf("; %s rows: %d (claim: none)", cs$no_route, sum(i)))
    }
    if (!is.null(cs$fails)) {
      # The failures count against the gate: the auto-gated failure rate fails.
      fr <- cs$fails
      i <- g$metric == fr[[1]] & (if (is.na(fr[[2]])) is.na(g$route) else g$route %in% fr[[2]])
      ok <- ok && sum(i) == 1L && isTRUE(g$gated[i]) && same(g$value[i], fr[[3]]) &&
        g$verdict[i] == "fail"
      det <- paste0(det, sprintf("\n%s%s = %s, gated %s, %s (claim: value %s, \"fail\")", fr[[1]],
                                 if (is.na(fr[[2]])) "" else paste0(":", fr[[2]]),
                                 toString(g$value[i]), toString(g$gated[i]),
                                 toString(g$verdict[i]), fr[[3]]))
    }
    zero <- g$gated & g$mcse %in% 0 & g$verdict %in% "fail"
    if (any(zero)) {
      det <- paste0(det, sprintf("\nfailed with MCSE 0: %s",
                                 toString(row_key(g$metric, g$estimand, g$route)[zero])))
    }
    report(paste("8", nm), ok, det)
  }
})

n_fail <- sum(results == "FAIL")
cat(sprintf("\n%d claims: %d PASS, %d FAIL, %d INFO\n", length(results), sum(results == "PASS"),
            n_fail, sum(results == "INFO")))
quit(save = "no", status = if (n_fail) 1L else 0L)
