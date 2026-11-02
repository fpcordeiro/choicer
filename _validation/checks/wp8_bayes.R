# ============================================================================
# WP8's checks: the Bayesian samplers' data (lib/dgp_bayes.R), the SBC loop
# (lib/sbc.R, loop half) and family F9. Sections (--only, comma-separated):
#   prior   prior draws (--n_prior each) against the priors' moments: the
#           normals' means and variances, the inverse Wisharts' and inverse
#           gammas' means; choicer's own IW and IG draws (riwishart_cpp(),
#           hb_test_sigma_d2_gibbs()) against the same means, so that the
#           parameterizations match the samplers'.
#   data    data from prior draws: structure (balanced, one choice per MNP
#           situation, at most one per HB task) and the prior predictive
#           spread of the shares (--n_pp datasets per model); theta0's
#           coherence (D1-L's truths; the MNP's constants recomputed by
#           db_mnp_calibrate()); the DGPs' choice frequencies against the
#           suite's probabilities (GHK for MNP, Gauss-Hermite for HMNP, the
#           logit for HMNL) at theta0, 4 binomial SEs.
#   loglik  db_hmnl_loglik() against run_hmnlogit()'s loglik trace at kept
#           draws (beta_i draws kept), < 1e-10 relative.
#   smoke   F9 through run_rep(): every cell once, rep 2 too at the bottom
#           rungs (bvm_mnp_500, bvm_hmnl_200); records re-validated; the BvM
#           records' tail ESS per quantity; the credible-interval encoding;
#           metrics_cell() on the records; a table of times (and smallest
#           bulk and tail ESS) and the mc tier's estimate (R$mc x seconds).
#   q14     ref_q14(): choicer's MNP against bayesm::rmnpGibbs.
#   q4      (not in the default) Q4's mini-SBC: F9's quick cells, R$quick
#           reps each, on --workers PSOCK workers at one thread, q4_check();
#           the wall time and the sum of the reps' times (one thread's).
# Usage (the dev build first on R_LIBS; one thread; each section < 10 min):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 \
#     Rscript _validation/checks/wp8_bayes.R [--only=prior,data] [--workers=8]
# Exits with status 1 on any FAIL.
# ============================================================================

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg))), ".."))
sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
b <- boot$vs_boot(root)
suppressPackageStartupMessages(library(data.table))
args <- vs_args(defaults = list(only = "prior,data,loglik,smoke,q14", n_prior = 4000L,
                                n_pp = 200L, workers = 8L))
vs_threads(1L)
options(width = 160L)
only <- trimws(strsplit(args$only, ",", fixed = TRUE)[[1]])
cat(sprintf("wp8_bayes: choicer %s from %s; sections %s\n", b$build$version, b$build$path,
            paste(only, collapse = ", ")))
t_start <- proc.time()[["elapsed"]]

verdicts <- logical(0)
report <- function(label, value, tol, pass = isTRUE(value < tol)) {
  cat(sprintf("%s  %-74s %11.4g  (threshold %.3g)\n", if (pass) "PASS" else "FAIL", label,
              value, tol))
  verdicts[length(verdicts) + 1L] <<- pass
}
now <- function() proc.time()[["elapsed"]]
fam <- load_family("F9", b$L)

# --- prior ------------------------------------------------------------------------

# |sample mean - mean| / its SE, over a matrix of draws (one column each).
z_mean <- function(x, mu) abs(colMeans(x) - mu) / (apply(x, 2L, stats::sd) / sqrt(nrow(x)))

if ("prior" %in% only) {
  cat("\n== prior: draws against the priors' moments ==\n")
  n <- args$n_prior
  for (kind in names(DB_PRIOR)) for (model in names(DB_PRIOR[[kind]])) {
    pr <- DB_PRIOR[[kind]][[model]]
    d <- vs_with_seed(1L, replicate(n, db_prior_draw(model, pr), simplify = FALSE))
    lab <- sprintf("%s/%s", kind, model)
    if (model == "mnp") {
      B <- t(vapply(d, `[[`, numeric(7L), "beta"))
      S <- t(vapply(d, function(x) .db_vech(x$Sigma), numeric(6L)))
      ES <- .db_vech(pr$V / (pr$nu - DB_P - 1))
    } else {
      B <- cbind(t(vapply(d, `[[`, numeric(2L), "b")),
                 t(vapply(d, `[[`, numeric(2L), "theta")))
      pr_mu <- c(pr$b_bar, pr$theta_bar)
      S <- t(vapply(d, function(x) c(.db_vech(x$W), x$sigma_d2, sigma2 = x$sigma2 %||% NA,
                                     xi2 = mean(x$xi^2)), numeric(6L)))
      ES <- c(.db_vech(pr$V / (pr$nu - 3)), pr$sd_prior$d0 / (pr$sd_prior$c0 - 1),
              if (model == "hmnp") pr$s0 / (pr$a0 - 1) else NA,
              pr$sd_prior$d0 / (pr$sd_prior$c0 - 1))
    }
    mu <- if (model == "mnp") pr$beta_bar else pr_mu
    vv <- if (model == "mnp") diag(solve(pr$A)) else c(diag(solve(pr$A)), diag(solve(pr$A_theta)))
    report(sprintf("%s: normals' means, max |z|", lab), max(z_mean(B, mu)), 4)
    rv <- apply(B, 2L, stats::var) / vv - 1
    report(sprintf("%s: normals' variances, max |var / prior - 1| / (4 sqrt(2/n))", lab),
           max(abs(rv)) / (4 * sqrt(2 / (n - 1))), 1)
    k <- !is.na(ES)
    report(sprintf("%s: IW / IG means (W or Sigma, sigma_d2, sigma2, xi^2), max |z|", lab),
           max(z_mean(S[, k, drop = FALSE], ES[k])), 4)
  }
  # choicer's own draws with the same arguments: IW(nu, V) and IG(c0, d0).
  pr <- DB_PRIOR$sbc$mnp
  ci <- t(vapply(seq_len(n), function(s) .db_vech(choicer:::riwishart_cpp(pr$nu, pr$V, s)),
                 numeric(6L)))
  report("choicer's riwishart_cpp(nu, V): mean against V / (nu - p - 1), max |z|",
         max(z_mean(ci, .db_vech(pr$V / (pr$nu - DB_P - 1)))), 4)
  sdp <- DB_PRIOR$sbc$hmnl$sd_prior
  ig <- choicer:::hb_test_sigma_d2_gibbs(0, n, 1, FALSE, 1, sdp$c0, sdp$d0)
  mine <- vs_with_seed(2L, replicate(n, .db_rinvgamma(sdp$c0 + 0.5, sdp$d0)))
  report("choicer's sigma_d2 | xi = 0 against IG(c0 + 1/2, d0) here: |z| of the means",
         abs(mean(ig) - mean(mine)) / sqrt(stats::var(ig) / n + stats::var(mine) / n), 4)
}

# --- data -------------------------------------------------------------------------

if ("data" %in% only) {
  cat("\n== data: datasets from prior draws, theta0, and the choice probabilities ==\n")
  for (v in DB_VARIANTS) {
    m <- DB_MODEL[[v]]
    s <- db_simulate(if (m == "mnp") 300L else 100L, 11L, v, prior = DB_PRIOR$sbc[[m]])
    dt <- s$data
    per <- dt[, .N, by = id]$N
    nch <- dt[, sum(choice), by = id]$V1
    ok <- all(per == if (m == "mnp") 4L else 6L) &&
      (if (m == "mnp") all(nch == 1L) else all(nch <= 1L)) && all(is.finite(s$truth))
    report(sprintf("%s: a prior draw's data are balanced with valid choices (1 = yes)", v),
           as.numeric(ok), 0.5, pass = ok)
    cat(sprintf("     shares: %s\n", paste(sprintf("%s %.3f", names(s$meta$shares),
                                                  s$meta$shares), collapse = ", ")))
    t0 <- now()
    pp <- db_prior_predictive(v, n = args$n_pp)
    cat(sprintf("     prior predictive, %d datasets (%.0f s): share quantiles 5/50/95%%\n",
                args$n_pp, now() - t0))
    print(round(pp$quantiles, 3))
    cat(sprintf("     datasets with a share < 0.02: %.3f; > 0.9: %.3f\n",
                pp$extreme_rate[[1]], pp$extreme_rate[[2]]))
    report(sprintf("%s: prior predictive datasets with a share < 0.02 or > 0.9, rate", v),
           mean(apply(pp$shares, 1L, min) < 0.02 | apply(pp$shares, 1L, max) > 0.9), 0.20)
  }
  # theta0's coherence.
  d1 <- get_dgp("D1")$theta0$labeled
  report("hmnl theta0 = D1-L's (WP4), max |diff|", max(abs(DB_THETA0$hmnl - d1[names(DB_THETA0$hmnl)])),
         1e-15)
  sh <- db_simulate(50L, 1L, "hmnp")$truth
  report("hmnp theta0 = db_simulate()'s truth at theta0, max |diff|",
         max(abs(DB_THETA0$hmnp - sh[names(DB_THETA0$hmnp)])), 1e-14)
  t0 <- now()
  cal <- db_mnp_calibrate()
  report(sprintf("mnp theta0 = db_mnp_calibrate()'s (%.0f s), max |diff|", now() - t0),
         max(abs(cal$theta - DB_THETA0$default[names(cal$theta)])), 1e-9)
  ps <- db_mnp_pop_shares(DB_THETA0$default, seed = 2L)
  report("mnp theta0: population shares against Q (other GHK draws), max |diff|",
         max(abs(ps - DB_Q[names(ps)])), 1e-3)
  # The DGPs' choice frequencies against the suite's probabilities at theta0.
  N <- 50000L
  s <- db_simulate(N, 3L, "default")
  U <- vs_with_seed(4L, array(stats::runif(N * DB_GHK_DRAWS * DB_P), c(N, DB_GHK_DRAWS, DB_P)))
  th <- DB_THETA0$default
  P <- vapply(0:DB_P, function(c) {
    db_mnp_prob(th[DB_BETA], .db_unvech(th[.DB_SIGMA_NAMES]), s$meta$Xd, rep(c, N), U)
  }, numeric(N))
  cnt <- tabulate(s$meta$y + 1L, DB_P + 1L)
  z <- abs(cnt - colSums(P)) / sqrt(colSums(P * (1 - P)))
  report(sprintf("mnp: choice counts (N = %d) against GHK's expected counts, max |z|", N),
         max(z), 4)
  for (v in c("hmnl", "hmnp")) {
    s <- db_simulate(5000L, 4L, v)
    lay <- db_hb_layout(s$data)
    obs <- tabulate(lay$chosen + 1L, lay$J + 1L)
    P <- vapply(0:lay$J, function(k) {
      l2 <- lay
      l2$chosen <- rep(k, lay$n_tasks)
      exp(if (v == "hmnl") db_hmnl_logp(l2, s$meta$tastes, s$meta$delta[DB_HB_ALTS])
          else db_hmnp_logp(l2, s$meta$tastes, s$meta$delta[DB_HB_ALTS]))
    }, numeric(lay$n_tasks))
    report(sprintf("%s: probabilities sum to 1 over the 7 options (24-node GH), max |sum - 1|", v),
           max(abs(rowSums(P) - 1)), 1e-5)
    z <- abs(obs - colSums(P)) / sqrt(colSums(P * (1 - P)))
    report(sprintf("%s: choice counts (40,000 tasks) against the expected, max |z|", v),
           max(z), 4)
  }
}

# --- loglik -----------------------------------------------------------------------

if ("loglik" %in% only) {
  cat("\n== loglik: db_hmnl_loglik() against run_hmnlogit()'s trace ==\n")
  s <- db_simulate(50L, 5L, "hmnl", prior = DB_PRIOR$sbc$hmnl)
  fit <- suppressWarnings(suppressMessages(db_fit(
    "hmnl", s$data, DB_PRIOR$sbc$hmnl, list(R = 700, burn = 200, seed = 9), keep_beta_i = "draws")))
  lay <- db_hb_layout(s$data)
  D <- db_draws(fit, "hmnl")
  idx <- c(1L, 100L, 250L, 400L, 500L)
  B <- db_beta_i_draws(fit, lay, idx)
  mine <- vapply(seq_along(idx), function(k) {
    db_hmnl_loglik(lay, B[, , k], D[idx[k], paste0("delta_", DB_HB_ALTS)])
  }, 1)
  tr <- fit$chains[[1L]]$loglik[idx]
  report("HMNL: the suite's log-likelihood at 5 kept draws vs the trace, max rel. diff",
         max(abs(mine - tr) / abs(tr)), 1e-10)
}

# --- smoke ------------------------------------------------------------------------

if ("smoke" %in% only) {
  cat("\n== smoke: run_rep() at every F9 cell ==\n")
  tab <- list()
  recs <- list()
  for (cl in fam$cells$cell) {
    reps <- if (cl %in% c("bvm_mnp_500", "bvm_hmnl_200")) 1:2 else 1L
    for (r in reps) {
      rec <- run_rep(fam, cl, r, cell_methods(fam, cl))[[1L]]
      recs[[length(recs) + 1L]] <- rec
      valid <- tryCatch({
        do.call(new_rep_record, rec[c("family", "cell", "method", "rep", "seed", "ok", "error",
                                      "warnings", "est", "vcov_c", "joint_wald", "lr", "derived",
                                      "diag", "meta")])
        TRUE
      }, error = function(e) FALSE)
      s <- rec$diag$sbc
      tab[[length(tab) + 1L]] <- data.frame(
        family = rec$family, cell = rec$cell, method = rec$method, rep = rec$rep, ok = rec$ok,
        valid = valid, sim_sec = round(rec$diag$sim_sec, 2), unit_sec = round(rec$diag$unit_sec, 1),
        fit_sec = round(rec$diag$time_sec %||% NA_real_, 1),
        iterations = rec$diag$iterations %||% NA, code = rec$diag$code %||% NA,
        ess_min = round(s$ess_min %||% rec$diag$mcmc$ess_min %||% NA_real_),
        ess_tail_min = round(rec$diag$mcmc$ess_tail_min %||% NA_real_),
        flagged = s$flagged %||% NA, n_ranks = length(s$ranks),
        conv_warning = !isTRUE(rec$diag$converged) && !is.na(rec$diag$converged),
        warnings = paste(unique(rec$warnings$class), collapse = " "),
        error = substr(rec$error %||% "", 1, 60), stringsAsFactors = FALSE)
      report(sprintf("F9/%s/%s rep %d: ok and valid", rec$cell, rec$method, rec$rep),
             as.numeric(rec$ok && valid), 0.5, pass = rec$ok && valid)
      if (rec$ok && is.null(s)) {
        et <- rec$diag$mcmc$ess_tail
        okt <- identical(names(et), rec$est$name) && all(is.finite(et))
        report(sprintf("F9/%s rep %d: tail ESS per quantity in diag$mcmc (min %.0f)", rec$cell,
                       rec$rep, rec$diag$mcmc$ess_tail_min), as.numeric(okt), 0.5, pass = okt)
      }
      if (rec$ok && !is.null(s)) {
        need <- c(names(rec$diag$sbc$ranks))
        okr <- all(c("loglik") %in% need) && all(s$ranks >= 0 & s$ranks <= s$L)
        report(sprintf("F9/%s: %d SBC ranks in 0..99, the log-likelihood's among them",
                       rec$cell, length(s$ranks)), as.numeric(okr), 0.5, pass = okr)
      }
    }
  }
  tab <- do.call(rbind, tab)
  print(tab, row.names = FALSE)
  # The credible interval's encoding: est - z se and est + z se are the
  # draws' 2.5% and 97.5% quantiles.
  x <- vs_with_seed(7L, matrix(stats::rgamma(30000L, 2), ncol = 3L,
                               dimnames = list(NULL, c("a", "b", "c"))))
  cr <- get("f9_cri", envir = environment(fam$simulate))(x, c(a = 1, b = 1, c = 1))
  z <- stats::qnorm(0.975)
  q <- apply(x, 2L, stats::quantile, c(0.025, 0.975), names = FALSE)
  report("cri: est -/+ z se against the 2.5% and 97.5% quantiles, max rel. diff",
         max(abs(c(cr$est - z * cr$se, cr$est + z * cr$se) - c(q[1, ], q[2, ])) / abs(c(q))), 1e-14)
  # metrics_cell() on the smoke records.
  key <- vapply(recs, function(r) paste(r$cell, r$method), "")
  mt <- do.call(rbind, lapply(unname(split(recs, key)), metrics_cell, gates = fam$gates))
  n_sbc <- sum(mt$metric == "sbc_ecdf")
  report(sprintf("metrics_cell(): %d sbc_ecdf rows (13 + 15 + 15 + quick 43)", n_sbc),
         abs(n_sbc - 86), 0.5)
  nb <- sum(mt$metric == "tail_lower" & grepl("^bvm", mt$cell) & mt$block == "derived")
  report(sprintf("metrics_cell(): %d credible-interval tail rows at the BvM cells (12+12+11+11)",
                 nb), abs(nb - 46), 0.5)
  # Times and the mc tier's estimate.
  sec <- tapply(tab$unit_sec + tab$sim_sec, tab$cell, mean)
  R_mc <- fam$R$mc
  est <- data.frame(cell = names(R_mc), R_mc = unname(R_mc), sec_per_rep = round(sec[names(R_mc)], 1),
                    cpu_h = round(R_mc * sec[names(R_mc)] / 3600, 2))
  cat("\nmc tier estimate (one thread):\n")
  print(est, row.names = FALSE)
  cat(sprintf("F9 mc total: %.1f CPU-h\n", sum(est$cpu_h)))
  rq <- fam$R$quick
  cat(sprintf("Q4 (quick cells, R = %s): %.0f s at one thread\n", paste(rq, collapse = "/"),
              sum(rq * sec[names(rq)])))
}

# --- q14 --------------------------------------------------------------------------

if ("q14" %in% only) {
  cat("\n== q14: ref_q14(), choicer's MNP against bayesm::rmnpGibbs ==\n")
  t0 <- now()
  rows <- ref_q14()
  print(rows, row.names = FALSE, digits = 4)
  g <- rows[!is.na(rows$pass), ]
  report(sprintf("Q1.4: %d check rows pass (%.0f s)", nrow(g), now() - t0), sum(!g$pass), 0.5)
}

# --- q4 ---------------------------------------------------------------------------

if ("q4" %in% only) {
  cat("\n== q4: the mini-SBC on PSOCK workers ==\n")
  rq <- fam$R$quick
  tasks <- unlist(lapply(names(rq), function(cl) lapply(seq_len(rq[[cl]]), function(r) {
    list(cell = cl, rep = r)
  })), recursive = FALSE)
  t0 <- now()
  cl <- parallel::makePSOCKcluster(args$workers)
  parallel::clusterExport(cl, "root", envir = globalenv())
  parallel::clusterEvalQ(cl, {
    sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
    b <- boot$vs_boot(root)
    vs_threads(1L)
    fam <- load_family("F9", b$L)
    NULL
  })
  t_boot <- now() - t0
  recs <- parallel::parLapplyLB(cl, tasks, function(tk) {
    run_rep(fam, tk$cell, tk$rep, cell_methods(fam, tk$cell))[[1L]]
  }, chunk.size = 1L)
  wall <- now() - t0
  parallel::stopCluster(cl)
  q4 <- q4_check(fam, recs)
  print(q4$rows, row.names = FALSE, digits = 3)
  cpu <- sum(vapply(recs, function(r) r$diag$unit_sec + r$diag$sim_sec, 1))
  fl <- tapply(vapply(recs, function(r) isTRUE(r$diag$sbc$flagged), TRUE),
               vapply(recs, `[[`, "", "cell"), mean)
  cat(sprintf("Q4: %d reps; wall %.0f s on %d workers (%.1f s starting them); one thread %.0f s\n",
              length(recs), wall, args$workers, t_boot, cpu))
  g <- q4$tab[q4$tab$metric == "sbc_ecdf", ]
  pm <- tapply(g$p, g$cell, min)
  cat("smallest raw ECDF-test p by cell:", paste(names(pm), signif(pm, 2), collapse = ", "), "\n")
  cat("flagged (ESS < 99) rate by cell:", paste(names(fl), round(fl, 2), collapse = ", "), "\n")
  report(sprintf("Q4: %d quantities, Holm p > 0.001 (failures)", nrow(q4$rows)),
         sum(!q4$rows$pass), 0.5)
  report("Q4: every rep ok", sum(!vapply(recs, `[[`, TRUE, "ok")), 0.5)
}

cat(sprintf("\nwp8_bayes: %d PASS, %d FAIL (%.0f s)\n", sum(verdicts), sum(!verdicts),
            proc.time()[["elapsed"]] - t_start))
quit(save = "no", status = if (all(verdicts)) 0L else 1L)
