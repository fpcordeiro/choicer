# _validation/blp_concentration_mc.R
#
# Monte Carlo of run_mxlogit(concentrate_out_asc = "blp"), the share-matching
# constants, against maximum likelihood on the same data and draws: the
# numbers behind §5.9 of the mixed logit's math article and NEWS.
#
# Usage, from the package root, with the build under test first in R_LIBS:
#   R_LIBS=<lib> Rscript _validation/blp_concentration_mc.R [QUICK=TRUE] \
#     [R=1000] [designs=cross,panel,base,markets] [workers=<cores - 2>] [out=<dir>]
#
# Designs (2,000 choice situations, J = 10 alternatives, S = 100 draws):
#   cross    a cross-section with an outside option
#   panel    500 decision makers with 4 situations each
#   base     a cross-section without an outside option
#   markets  two markets of five alternatives, each decision maker in one
# DGP: V = x1 + (-1 + e1_n) w1 + (0.5 + 0.7 e2_n) w2 + delta_j + Gumbel, with
# x1 ~ N(0, 1), w ~ N(0, 1.5^2), e ~ N(0, 1) per decision maker, delta from -1
# to 1, choice sets random subsets of at least half of the market's
# alternatives; fitted with X = (x1, w1, w2), W = (w1, w2), a diagonal
# Cholesky (rc_mean = FALSE).
# Estimators: ML (Newton-certified); "blp" with the sample's shares, with the
# DGP's exact within-market shares (from 10^6 simulated situations), and
# with noisy targets from n = N and n = 10 N choices per market (target_vcov).
# Output in <out> (default _validation/output/blp_mc): mc_<design>.rds, one
# record per replication, and summary.csv. The printed table gives, per
# design and estimator, the SD of the estimates against ML's, mean SE / SD,
# and 95% coverage with the "hessian" and the sandwich standard errors, for
# theta_c and the constants (averaged), the largest bias over SD, the work per
# fit, and the convergence record: failed inner solves, fits that warned, the
# largest outer decrement and the optimizer stops with a code other than
# nloptr's convergence codes 1 to 4.
# QUICK=TRUE: 50 replications of the cross design.

suppressPackageStartupMessages({
  library(choicer)
  library(data.table)
  library(parallel)
})
opt <- list(QUICK = "FALSE", R = "1000", designs = "cross,panel,base,markets",
            workers = as.character(max(1L, detectCores() - 2L)),
            out = file.path("_validation", "output", "blp_mc"))
for (a in commandArgs(TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1L]]
  opt[[kv[1L]]] <- kv[2L]
}
quick <- isTRUE(as.logical(opt$QUICK))
R <- if (quick) 50L else as.integer(opt$R)
designs <- if (quick) "cross" else strsplit(opt$designs, ",", fixed = TRUE)[[1L]]
workers <- as.integer(opt$workers)
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
set_num_threads(1L)  # the replications run in parallel instead

# --- The DGP -----------------------------------------------------------------

sim_raw <- function(seed, N, J, T, outside, markets, beta = c(1, -1, 0.5),
                    sig = c(1, 0.7), delta = seq(-1, 1, length.out = J), wsd = 1.5) {
  set.seed(seed)
  n_sit <- N * T
  pid_sit <- rep(seq_len(N), each = T)
  mk_alt <- rep(seq_len(markets), each = J %/% markets)
  # each decision maker's market (one market draws nothing: the random
  # streams of the runs behind the docs' numbers)
  mk_sit <- if (markets > 1L) {
    rep(sample.int(markets, N, replace = TRUE), each = T)
  } else {
    rep(1L, n_sit)
  }
  alts <- lapply(seq_len(n_sit), function(t) {
    pool <- which(mk_alt == mk_sit[t])
    k <- ceiling(length(pool) / 2):length(pool)
    k <- k[sample.int(length(k), 1L)]
    sort(pool[sample.int(length(pool), k)])
  })
  id <- rep(seq_len(n_sit), lengths(alts))
  dt <- data.table(id = id, alt = unlist(alts), pid = pid_sit[id], mk = mk_sit[id])
  dt[, x1 := stats::rnorm(.N)]
  dt[, w1 := stats::rnorm(.N, sd = wsd)]
  dt[, w2 := stats::rnorm(.N, sd = wsd)]
  e <- matrix(stats::rnorm(2L * N), 2L)
  dt[, V := beta[1] * x1 + (beta[2] + sig[1] * e[1, pid]) * w1 +
            (beta[3] + sig[2] * e[2, pid]) * w2 + delta[alt]]
  dt[]
}

design_args <- function(design) switch(design,
  cross = list(N = 2000L, T = 1L, outside = TRUE, markets = 1L),
  panel = list(N = 500L, T = 4L, outside = TRUE, markets = 1L),
  base = list(N = 2000L, T = 1L, outside = FALSE, markets = 1L),
  markets = list(N = 2000L, T = 1L, outside = TRUE, markets = 2L),
  stop("unknown design ", design))

sim_data <- function(seed, a, J = 10L, S = 100L) {
  dt <- sim_raw(seed, a$N, J, a$T, a$outside, a$markets)
  dt[, p := exp(V) / (sum(exp(V)) + a$outside), by = id]
  r <- stats::runif(max(dt$id))
  dt[, choice := { cp <- cumsum(p); as.integer(r[id[1L]] <= cp & r[id[1L]] > cp - p) },
     by = id]
  d <- suppressMessages(prepare_mxl_data(
    dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "w1", "w2"), random_var_cols = c("w1", "w2"),
    include_outside_option = a$outside, outside_opt_label = if (a$outside) 0L,
    person_col = if (a$T > 1L) "pid"))
  U <- if (is.null(d$Ti)) length(d$M) else length(d$Ti)
  delta <- seq(-1, 1, length.out = J)
  list(d = d, eta = get_halton_normals(S, U, 2L),
       truth = c(1, -1, 0.5, log(c(1, 0.7)), if (a$outside) delta else delta[-1L] - delta[1L]))
}

# The DGP's within-market shares of the inside alternatives: the exact
# probabilities averaged over 10^6 simulated situations
pop_shares <- function(a, J = 10L, n_pop = 1e6) {
  dt <- sim_raw(424242L, as.integer(n_pop / a$T), J, a$T, a$outside, a$markets)
  dt[, p := exp(V) / (sum(exp(V)) + a$outside), by = id]
  n_m <- dt[, .(n = uniqueN(id)), by = mk]
  s <- dt[, .(s = sum(p)), by = .(alt, mk)][n_m, on = "mk"][, s := s / n][order(alt)]
  list(within = s$s, market = s$mk)
}

# --- One replication ---------------------------------------------------------

quiet <- function(expr) {
  warns <- character()
  val <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(val = val, warnings = warns)
}

fit_blp <- function(x, ts = NULL, tv = NULL) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(quiet(run_mxlogit(input_data = x$d, eta_draws = x$eta,
                                  concentrate_out_asc = "blp", target_shares = ts,
                                  target_vcov = tv)),
                error = function(e) list(error = conditionMessage(e)))
  if (!is.null(r$error)) return(list(error = r$error))
  f <- r$val
  cr <- f$optimizer$concentration
  ser <- tryCatch(sqrt(pmax(diag(unname(suppressMessages(vcov(f, type = "robust")))), 0)),
                  error = function(e) rep(NA_real_, length(coef(f))))
  list(est = unname(coef(f)), se = unname(f$se), se_r = ser, evals = cr$evaluations,
       inner = cr$inner_iterations, solves = cr$solves, failed = cr$failed,
       decrement = cr$decrement, convergence = f$convergence,
       warnings = r$warnings, time = proc.time()[["elapsed"]] - t0)
}

one_rep <- function(r, design, a, pop) {
  seed <- 1000L * match(design, c("cross", "panel", "base", "markets")) + r
  x <- sim_data(seed, a)
  rec <- list(rep = r, seed = seed, truth = x$truth)
  t0 <- proc.time()[["elapsed"]]
  ml <- tryCatch(suppressWarnings(suppressMessages(run_mxlogit(
    input_data = x$d, eta_draws = x$eta, polish = "newton",
    polish_control = list(max_iter = 50L, tol_decrement = 1e-14)))),
    error = function(e) NULL)
  if (is.null(ml) || is.null(ml$vcov)) return(c(rec, list(error = "ML failed")))
  rec$ml <- list(est = unname(coef(ml)), se = unname(ml$se),
                 se_r = sqrt(diag(unname(suppressMessages(vcov(ml, type = "robust"))))),
                 time = proc.time()[["elapsed"]] - t0)
  rec$blp_sample <- fit_blp(x)
  rec$blp_exact <- fit_blp(x, pop$within)
  J <- length(pop$within)
  n_sit <- length(x$d$M)
  for (k in c(1L, 10L)) {
    n <- k * n_sit / a$markets
    set.seed(seed + 77L * k)
    s_hat <- numeric(J)
    Vs <- matrix(0, J, J)
    for (m in seq_len(a$markets)) {
      im <- which(pop$market == m)
      pin <- pop$within[im]
      probs <- if (a$outside) c(1 - sum(pin), pin) else pin
      repeat {
        sh <- as.numeric(stats::rmultinom(1L, n, probs)) / n
        if (all(sh > 0)) break
      }
      shin <- if (a$outside) sh[-1L] else sh
      s_hat[im] <- shin
      Vs[im, im] <- (diag(shin, length(im)) - tcrossprod(shin)) / n
    }
    rec[[paste0("blp_n", k)]] <- fit_blp(x, s_hat, Vs)
  }
  rec
}

# --- The runs and the summary ------------------------------------------------

ests <- c(ml = "ML", blp_sample = "blp, sample shares", blp_exact = "blp, exact shares",
          blp_n1 = "blp, shares from n = N", blp_n10 = "blp, shares from n = 10 N")
rows <- list()
for (design in designs) {
  a <- design_args(design)
  pop <- pop_shares(a)
  t0 <- Sys.time()
  res <- mclapply(seq_len(R), one_rep, design = design, a = a, pop = pop,
                  mc.cores = workers, mc.preschedule = FALSE)
  saveRDS(res, file.path(opt$out, sprintf("mc_%s.rds", design)))
  res <- Filter(function(r) is.list(r) && is.null(r$error), res)
  cat(sprintf("%s: %d replications, %.1f min\n", design, length(res),
              as.numeric(Sys.time() - t0, units = "mins")))
  truth <- res[[1L]]$truth
  p <- 5L
  dj <- (p + 1L):length(truth)
  sd_ml <- NULL
  for (e in names(ests)) {
    ok <- Filter(function(r) !is.null(r[[e]]) && is.null(r[[e]]$error), res)
    if (!length(ok)) next
    est <- do.call(rbind, lapply(ok, function(r) r[[e]]$est))
    se <- do.call(rbind, lapply(ok, function(r) r[[e]]$se))
    ser <- do.call(rbind, lapply(ok, function(r) r[[e]]$se_r))
    err <- sweep(est, 2L, truth)
    sdv <- apply(est, 2L, stats::sd)
    if (e == "ml") sd_ml <- sdv
    grp <- list(theta_c = seq_len(p), constants = dj)
    for (g in names(grp)) {
      j <- grp[[g]]
      rows[[length(rows) + 1L]] <- data.frame(
        design = design, estimator = ests[[e]], params = g, fits = nrow(est),
        failed = length(res) - nrow(est),
        sd_vs_ml = mean(sdv[j] / sd_ml[j]),
        se_over_sd = mean(colMeans(se[, j, drop = FALSE]) / sdv[j]),
        cover95 = mean(abs(err[, j]) <= 1.96 * se[, j]),
        cover95_sandwich = mean(abs(err[, j]) <= 1.96 * ser[, j]),
        max_bias_sd = max(abs(colMeans(err[, j, drop = FALSE])) / sdv[j]),
        evals = if (e == "ml") NA else stats::median(sapply(ok, function(r) r[[e]]$evals)),
        inner_per_solve = if (e == "ml") NA else
          stats::median(sapply(ok, function(r) r[[e]]$inner / r[[e]]$solves)),
        failed_solves = if (e == "ml") NA else sum(sapply(ok, function(r) r[[e]]$failed)),
        warned = if (e == "ml") NA else sum(sapply(ok, function(r) length(r[[e]]$warnings) > 0L)),
        max_decrement = if (e == "ml") NA else
          max(sapply(ok, function(r) r[[e]]$decrement), na.rm = TRUE),
        stops_not_converged = if (e == "ml") NA else
          sum(sapply(ok, function(r) !isTRUE(r[[e]]$convergence %in% 1:4))),
        seconds = stats::median(sapply(ok, function(r) r[[e]]$time)))
    }
  }
}
tab <- do.call(rbind, rows)
utils::write.csv(tab, file.path(opt$out, "summary.csv"), row.names = FALSE)
num <- vapply(tab, is.numeric, TRUE) & names(tab) != "max_decrement"
tab[num] <- lapply(tab[num], round, 3)
tab$max_decrement <- signif(tab$max_decrement, 2)
op <- options(width = 200)
print(tab, row.names = FALSE)
options(op)
