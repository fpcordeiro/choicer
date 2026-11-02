# Self-check of lib/simulate.R (WP0b): the suite's probability engines agree
# with choicer, its choice draws reproduce its probabilities, and its
# parameter names are choicer's.
#   1. MNL (constants with and without an outside option, use_asc = FALSE)
#      and NL (two nests with a parameter and a singleton, with and without
#      an outside option): probabilities at theta0 (not the estimate) equal
#      predict() with the fit's coefficients replaced, max abs diff < 1e-12;
#      the weighted log-likelihood at the estimate equals logLik(fit).
#   2. MXL (cross-section and panel, normal and shifted log-normal, correlated
#      and uncorrelated, K_w = 2): the simulated log-likelihood on a shared
#      eta_draws cube equals choicer's (every parameter fixed by lower =
#      upper; and the internal kernel) to < 1e-10 relative; probabilities
#      equal predict() on choicer's store-mode draws.
#   3. Choice frequencies over 10^6 situations of 5 types (MNL, NL, MXL with
#      tastes from draw_tastes()) against the suite's probabilities (MXL: the
#      Gauss-Hermite integral): chi-square p-value per type, FAIL below 1e-4.
#   4. theta_names() and theta_blocks() equal names(coef(fit)) and
#      fit$param_map for every fit; pack_theta() and unpack_theta() invert.
#
# Usage (one OpenMP thread, the dev build first on R_LIBS):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 \
#     Rscript _validation/checks/wp0_simulate.R
# Exits with status 1 on any FAIL.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else "."
source(file.path(here, "..", "lib", "simulate.R"))
suppressPackageStartupMessages({
  library(choicer)
  library(data.table)
})
if (packageVersion("choicer") < "0.2.1.9000") stop("choicer >= 0.2.1.9000 needed.")
cat("choicer", format(packageVersion("choicer")), "from",
    getNamespaceInfo("choicer", "path"), "\n")
set_num_threads(1L)
t_start <- proc.time()[["elapsed"]]

verdicts <- logical(0)
report <- function(check, label, value, tol, pass = isTRUE(value < tol)) {
  cat(sprintf("%s  %-3s %-58s %10.3e  (tol %.0e)\n", if (pass) "PASS" else "FAIL",
              check, label, value, tol))
  verdicts[length(verdicts) + 1L] <<- pass
}
info <- function(...) cat("      ", ..., "\n", sep = "")

# Fits record choicer's warnings (counted, printed at the end) and drop its
# messages.
fit_warnings <- character(0)
quiet_fit <- function(expr) {
  withCallingHandlers(suppressMessages(expr), warning = function(cnd) {
    fit_warnings <<- c(fit_warnings, substr(conditionMessage(cnd), 1, 100))
    invokeRestart("muffleWarning")
  })
}

# Check 4 on one fit.
check_names <- function(model, spec, fit, label) {
  ok <- identical(theta_names(model, spec), names(coef(fit)))
  pm <- lapply(fit$param_map, as.integer)
  tb <- lapply(theta_blocks(model, spec), as.integer)
  ok_pm <- identical(pm[order(names(pm))], tb[order(names(tb))])
  report("4", paste(label, "names / param_map"), as.numeric(!(ok && ok_pm)), 0.5)
  if (!ok) info("suite: ", paste(theta_names(model, spec), collapse = " "),
                " | choicer: ", paste(names(coef(fit)), collapse = " "))
}

# A random long design: N situations, each offering a random 2..J of `alts`,
# rows sorted by situation and alternative code (choicer's order); x1, x2,
# w1, w2 ~ N(0, 1); decision makers of T situations each; weights per
# decision maker.
toy_data <- function(N, alts, seed, T = 1L) {
  set.seed(seed)
  J <- length(alts)
  m <- sample(2:J, N, replace = TRUE)
  offered <- apply(matrix(runif(N * J), N, J), 1, rank) <= rep(m, each = J)
  idx <- which(offered) - 1L                         # J x N, column-major
  sit <- idx %/% J + 1L
  n <- length(sit)
  dt <- data.table(id = sit, alt = sort(alts)[idx %% J + 1L],
                   x1 = rnorm(n), x2 = rnorm(n), w1 = rnorm(n), w2 = rnorm(n),
                   person = (sit - 1L) %/% T + 1L)
  w_person <- runif(max(dt$person), 0.5, 2)
  dt[, w := w_person[person]]
  dt
}

# --- 1. MNL and NL probabilities against predict() ------------------------------

labels <- c("car", "bus", "air", "train", "walk")
asc_all <- c(air = 0.2, bus = -0.3, car = 0.4, train = 0.1, walk = -0.5)
nest_of <- c(air = "fast", train = "fast", bus = "road", car = "road", walk = "foot")
lambda0 <- c(fast = 0.5, road = 0.75)

for (model in c("mnl", "nl")) {
  cfgs <- if (model == "mnl") {
    list(list(outside = FALSE, use_asc = TRUE), list(outside = TRUE, use_asc = TRUE),
         list(outside = FALSE, use_asc = FALSE))
  } else {
    list(list(outside = FALSE, use_asc = TRUE), list(outside = TRUE, use_asc = TRUE))
  }
  for (cfg in cfgs) for (seed in 1:3) {
    label <- sprintf("%s outside=%s asc=%s seed=%d", toupper(model), cfg$outside,
                     cfg$use_asc, seed)
    dt <- toy_data(300L, labels, seed)
    if (model == "nl") dt[, nest := nest_of[alt]]
    spec <- list(beta = c("x1", "x2"), alts = labels, outside = cfg$outside,
                 use_asc = cfg$use_asc, nests = unname(nest_of[labels]))
    blocks <- list(beta = c(x2 = -0.5, x1 = 0.8),
                   asc = if (cfg$use_asc) asc_all,
                   lambda = if (model == "nl") lambda0)
    theta0 <- pack_theta(model, spec, blocks[lengths(blocks) > 0])
    par0 <- unpack_theta(model, spec, theta0)
    d <- design(dt, spec)
    p0 <- if (model == "mnl") mnl_probs(d, par0) else nl_probs(d, par0)
    dt[, choice := draw_choice(p0$p, p0$p0, d$sit, seed = 100 + seed)$choice]
    fit <- quiet_fit(if (model == "mnl") {
      run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"), use_asc = cfg$use_asc,
                  include_outside_option = cfg$outside, weights_col = "w")
    } else {
      run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), nest_col = "nest",
                    use_asc = cfg$use_asc, include_outside_option = cfg$outside,
                    weights_col = "w")
    })
    check_names(model, spec, fit, label)
    report("1", paste(label, "alt codes = choicer's"),
           as.numeric(!identical(as.integer(fit$data$alt_idx), d$alt)), 0.5)
    # Log-likelihood at the estimate, weighted
    w_sit <- dt[, w[1L], by = id][[2L]]
    par_hat <- unpack_theta(model, spec, coef(fit))
    chosen <- chosen_rows(dt$choice, d$sit)
    ll <- if (model == "mnl") mnl_loglik(d, par_hat, chosen, w_sit) else
      nl_loglik(d, par_hat, chosen, w_sit)
    report("1", paste(label, "loglik vs logLik(fit), rel"),
           abs(ll - as.numeric(logLik(fit))) / abs(ll), 1e-10)
    # predict() at theta0, which is not the estimate
    gap <- max(abs(coef(fit)[names(theta0)] - theta0))
    pr_hat <- predict(fit, type = "probabilities")$choice_prob[, 1]
    fit$coefficients[names(theta0)] <- theta0
    pr <- predict(fit, type = "probabilities")$choice_prob[, 1]
    diff <- max(abs(pr - p0$p))
    if (cfg$outside) diff <- max(diff, abs(1 - rowsum(pr, d$sit)[, 1] - p0$p0))
    report("1", paste(label, "max|p - predict| at theta0"), diff, 1e-12)
    info(sprintf("max|theta0 - estimate| = %.3f; predict moved by %.3e", gap,
                 max(abs(pr - pr_hat))))
  }
}

# --- 2. MXL simulated log-likelihood and probabilities ----------------------------

mxl_alts <- c(1L, 2L, 10L, 3L)                       # codes sort numerically
mxl_cfgs <- list()
for (panel in c(FALSE, TRUE)) for (corr in c(FALSE, TRUE)) for (dist in list(c(0L, 0L), c(1L, 0L))) {
  mxl_cfgs[[length(mxl_cfgs) + 1L]] <- list(panel = panel, corr = corr, dist = dist,
                                            rc_mean = TRUE, outside = corr)
}
mxl_cfgs <- c(mxl_cfgs, list(
  list(panel = FALSE, corr = TRUE, dist = c(1L, 0L), rc_mean = FALSE, outside = FALSE),
  list(panel = TRUE, corr = FALSE, dist = c(0L, 1L), rc_mean = FALSE, outside = TRUE)))
L_corr <- matrix(c(0.8, 0.3, 0, 0.6), 2)
S <- 25L

for (i in seq_along(mxl_cfgs)) {
  cfg <- mxl_cfgs[[i]]
  label <- sprintf("MXL %s corr=%s dist=%s mean=%s out=%s", if (cfg$panel) "panel" else "cross",
                   cfg$corr, paste(cfg$dist, collapse = ""), cfg$rc_mean, cfg$outside)
  T <- if (cfg$panel) 4L else 1L
  dt <- toy_data(300L, mxl_alts, 20 + i, T)
  spec <- list(beta = "x1", alts = mxl_alts, outside = cfg$outside, random = c("w1", "w2"),
               rc_dist = cfg$dist, rc_mean = cfg$rc_mean, correlated = cfg$corr)
  blocks <- list(beta = 0.5, mu = if (cfg$rc_mean) c(w1 = 0.4, w2 = -0.3),
                 L = if (cfg$corr) L_corr else diag(diag(L_corr)),
                 asc = c(`1` = 0.1, `2` = 0.3, `3` = -0.2, `10` = 0.5))
  theta0 <- pack_theta("mxl", spec, blocks[lengths(blocks) > 0])
  par0 <- unpack_theta("mxl", spec, theta0)
  d <- design(dt, spec, person = if (cfg$panel) "person")
  U <- if (cfg$panel) max(d$unit) else max(d$sit)
  tastes <- draw_tastes(U, par0$mu, par0$L, par0$rc_dist, par0$rc_mean, seed = 200 + i)
  pc <- mxl_cond_probs(d, par0, tastes)
  dt[, choice := draw_choice(pc$p, pc$p0, d$sit, seed = 300 + i)$choice]
  idat <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("w1", "w2"), weights_col = "w",
                           include_outside_option = cfg$outside, rc_correlation = cfg$corr,
                           person_col = if (cfg$panel) "person")
  set.seed(400 + i)
  eta <- array(rnorm(2 * S * U), c(2L, S, U))
  # Every parameter fixed at theta0: fit$loglik is choicer's at theta0
  fit <- quiet_fit(run_mxlogit(input_data = idat, eta_draws = eta, rc_dist = cfg$dist,
                               rc_mean = cfg$rc_mean, rc_correlation = cfg$corr,
                               include_outside_option = cfg$outside, theta_init = theta0,
                               lower = theta0, upper = theta0))
  check_names("mxl", spec, fit, label)
  w_unit <- if (cfg$panel) dt[, w[1L], by = person][[2L]] else dt[, w[1L], by = id][[2L]]
  chosen <- chosen_rows(dt$choice, d$sit)
  ll <- mxl_loglik(d, par0, eta, chosen, w_unit)
  report("2", paste(label, "loglik vs fixed fit, rel"), abs(ll - fit$loglik) / abs(ll), 1e-10)
  k <- choicer:::mxl_loglik_gradient_parallel(
    theta0, idat$X, idat$W, idat$alt_idx, idat$choice_idx, idat$M, idat$weights, eta,
    cfg$dist, cfg$corr, cfg$rc_mean, TRUE, cfg$outside, Ti = idat$Ti)
  report("2", paste(label, "loglik vs kernel, rel"), abs(ll + k$objective) / abs(ll), 1e-10)
  # predict() on choicer's store-mode draws (one block per situation), at a
  # theta other than the one the fit holds
  theta1 <- theta0 + 0.05
  fit$coefficients[names(theta1)] <- theta1
  pr <- predict(fit, type = "probabilities")
  d_cs <- d
  d_cs$unit <- NULL
  mine <- mxl_probs(d_cs, unpack_theta("mxl", spec, theta1),
                    get_halton_normals(S, max(d$sit), 2L))
  diff <- max(abs(pr$choice_prob[, 1] - mine$p))
  if (cfg$outside) diff <- max(diff, abs(pr$choice_prob_outside[, 1] - mine$p0))
  report("2", paste(label, "max|p - predict| (Halton)"), diff, 1e-12)
}

# --- 3. Choice frequencies ---------------------------------------------------------

# 5 situation types (sets of 2..5 alternatives), each repeated R times.
R <- 200000L
types <- toy_data(5L, labels, 7L)
spec <- list(beta = c("x1", "x2"), alts = labels, outside = TRUE,
             nests = unname(nest_of[labels]), random = c("w1", "w2"),
             rc_dist = c(1L, 0L), rc_mean = TRUE, correlated = TRUE)
d5 <- design(types, spec)
n5 <- length(d5$sit)
repl <- function(d, R) {
  r <- rep(seq_len(n5), R)
  list(sit = d$sit[r] + 5L * rep(0:(R - 1L), each = n5), alt = d$alt[r],
       X = d$X[r, , drop = FALSE], W = d$W[r, , drop = FALSE])
}
big <- repl(d5, R)
cat(sprintf("      frequency design: %d situations, %d rows\n", 5L * R, length(big$sit)))

# Chi-square goodness of fit per type: counts by position (and outside)
# against R times the probabilities.
gof <- function(model, chosen, p5) {
  pos <- chosen - (rep(0:(R - 1L), each = 5L) * n5)   # row within the 5-type block
  type <- rep(1:5, R)
  for (t in 1:5) {
    rows <- which(d5$sit == t)
    obs <- c(sum(chosen[type == t] == 0L), tabulate(pos[type == t & chosen > 0L] - rows[1L] + 1L,
                                                    length(rows)))
    expct <- R * c(p5$p0[t], p5$p[rows])
    stat <- sum((obs - expct)^2 / expct)
    pv <- pchisq(stat, length(obs) - 1L, lower.tail = FALSE)
    report("3", sprintf("%s type %d (%d alts + outside): chi2 %.2f, df %d, p-value", model, t,
                        length(rows), stat, length(obs) - 1L), pv, 1e-4, pass = pv >= 1e-4)
  }
}

par_mnl <- unpack_theta("mnl", spec, pack_theta("mnl", spec, list(beta = c(0.8, -0.5),
                                                                  asc = asc_all)))
tm <- system.time(pb <- mnl_probs(big, par_mnl))[["elapsed"]]
tc <- system.time(ch <- draw_choice(pb$p, pb$p0, big$sit, seed = 11))[["elapsed"]]
info(sprintf("mnl_probs %.2f s, draw_choice %.2f s", tm, tc))
gof("MNL", ch$chosen, mnl_probs(d5, par_mnl))

par_nl <- unpack_theta("nl", spec, pack_theta("nl", spec, list(beta = c(0.8, -0.5),
                                                               asc = asc_all, lambda = lambda0)))
tm <- system.time(pb <- nl_probs(big, par_nl))[["elapsed"]]
info(sprintf("nl_probs %.2f s", tm))
gof("NL", draw_choice(pb$p, pb$p0, big$sit, seed = 12)$chosen, nl_probs(d5, par_nl))

par_mxl <- unpack_theta("mxl", spec, pack_theta("mxl", spec, list(
  beta = c(0.8, -0.5), mu = c(-0.5, 0.3), L = L_corr, asc = asc_all)))
tm <- system.time({
  tastes <- draw_tastes(5L * R, par_mxl$mu, par_mxl$L, par_mxl$rc_dist, par_mxl$rc_mean, seed = 13)
  pb <- mxl_cond_probs(big, par_mxl, tastes)
})[["elapsed"]]
info(sprintf("draw_tastes + mxl_cond_probs %.2f s", tm))
# The integral over eta ~ N(0, I_2) by a 40 x 40 Gauss-Hermite product rule
Q <- 40L
jac <- matrix(0, Q, Q)
jac[cbind(1:(Q - 1L), 2:Q)] <- jac[cbind(2:Q, 1:(Q - 1L))] <- sqrt(1:(Q - 1L) / 2)
eg <- eigen(jac, symmetric = TRUE)
z1 <- sqrt(2) * eg$values
w1 <- eg$vectors[1L, ]^2
nodes <- rbind(rep(z1, Q), rep(z1, each = Q))
wq <- rep(w1, Q) * rep(w1, each = Q)
dq <- repl(d5, Q * Q)
B <- rc_coefs(par_mxl$mu, par_mxl$L, nodes, par_mxl$rc_dist, par_mxl$rc_mean)
pq <- mxl_cond_probs(dq, par_mxl, t(B)[rep(seq_len(Q * Q), each = 5L), ])
p5 <- list(p = drop(matrix(pq$p, n5) %*% wq), p0 = drop(matrix(pq$p0, 5L) %*% wq))
info(sprintf("Gauss-Hermite: total probability per type - 1: %.1e",
             max(abs(rowsum(p5$p, d5$sit)[, 1] + p5$p0 - 1))))
gof("MXL", draw_choice(pb$p, pb$p0, big$sit, seed = 14)$chosen, p5)

# --- 4. pack_theta() and unpack_theta() invert each other --------------------------

spec_r <- list(beta = c("x1", "x2"), alts = labels, outside = FALSE,
               nests = unname(nest_of[labels]), random = c("w1", "w2"),
               rc_dist = c(1L, 0L), rc_mean = TRUE, correlated = TRUE)
th <- pack_theta("mxl", spec_r, list(beta = c(0.8, -0.5), mu = c(-0.5, 0.3), L = L_corr,
                                     asc = asc_all))
pr_ <- unpack_theta("mxl", spec_r, rev(th))          # by name, any order
back <- pack_theta("mxl", spec_r, list(beta = pr_$beta, mu = pr_$mu, sigma = pr_$sigma,
                                       asc = pr_$delta))
err <- max(abs(back - th), abs(pr_$L - L_corr), abs(pr_$delta - (asc_all - asc_all[["air"]])))
report("4", "MXL pack/unpack round trip (asc normalized to air)", err, 1e-15)
th_nl <- pack_theta("nl", spec_r, list(beta = c(0.8, -0.5), lambda = lambda0, asc = asc_all))
pn <- unpack_theta("nl", spec_r, th_nl)
err <- max(abs(pn$lambda - c(fast = 0.5, foot = 1, road = 0.75)),
           abs(pack_theta("nl", spec_r, list(beta = pn$beta, lambda = unname(lambda0),
                                             asc = unname(th_nl[grep("^ASC_", names(th_nl))])))
               - th_nl))
report("4", "NL pack/unpack: lambda per nest code, singleton at 1", err, 1e-15)

# --- Summary ------------------------------------------------------------------------

tab <- table(fit_warnings)
for (k in seq_along(tab)) info(sprintf("fit warning x%d: %s", tab[[k]], names(tab)[k]))
cat(sprintf("\n%d checks, %d FAIL; %.1f s\n", length(verdicts), sum(!verdicts),
            proc.time()[["elapsed"]] - t_start))
if (any(!verdicts)) quit(status = 1L)
