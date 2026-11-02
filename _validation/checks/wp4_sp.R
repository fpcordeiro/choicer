# ============================================================================
# WP4's checks: D1 and D1-L (lib/dgp_sp.R) and family F4. Sections (--only,
# comma-separated):
#   q0     the DGP against choicer: reproducibility and the level balance;
#          D1's simulated log-likelihood on a shared cube against choicer's
#          fixed-parameter fit, panel (person_col) and cross-section, and on
#          the DGP's own tastes (a one-draw cube per person: the
#          probabilities the choices were drawn from), < 1e-10 relative;
#          D1-L's probabilities against an independent computation from
#          (b, W, delta) (< 1e-12) and its own-tastes likelihood against
#          choicer's MXL with the same model (fixed means in X, the tastes'
#          deviations in W, rc_mean = FALSE); D1-L's probit choices for one
#          fixed task over 10^6 persons against the probit probabilities by
#          Gauss-Hermite quadrature (tastes and the chosen row's shock),
#          chi-square p-value, FAIL below 1e-4.
#   calib  the moments at --calib_N persons (seed 1) against plan §4's
#          targets, D1's coefficient moments analytic beside simulated, and
#          sqrt(U)/S at every F4 cell.
#   ident  identification strength: one panel fit at 200 and at 800
#          persons (S = 500): each L parameter over its SE, and the implied
#          SDs (of g = L eta and of the coefficients) over their delta SEs;
#          the median over five more fits at 200.
#   hb     D1-L in prepare_hmnl_data()'s layout: one short run_hmnlogit()
#          (logit data) and one run_hmnprobit() (probit data), 100 persons.
#   smoke  F4 through run_rep() (--smoke=F4 or F4:<cell>[:<method>]): R = 1,
#          R = 2 at the bottom rung plus rep 10 there (the multi-start); each
#          record re-validated; a table of times.
#   plateau  F4's reps 1..K per cell (--plateau=200:40,800:10, ...; not in
#          the default sections): failures and fits with an L diagonal
#          driven to -Inf (diag$plateau: L_11, the time coefficient's SD to
#          0; L_22, Sigma singular), each refit from theta0 and with
#          scale_vars = "none" to see whether it is the maximum.
#   xtime  the composite fit at 3.2k with S = 200 and S = 800 (sqrt(U)/S 0.8
#          and 0.2): times and the estimates' gap over the Hessian SE.
# Usage (the dev build first on R_LIBS; one thread; each section < 10 min):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 \
#     Rscript _validation/checks/wp4_sp.R --only=q0,calib
# Exits with status 1 on any FAIL.
# ============================================================================

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg))), ".."))
sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
b <- boot$vs_boot(root)
suppressPackageStartupMessages(library(data.table))
args <- vs_args(defaults = list(only = "q0,calib,ident,hb,smoke,xtime", smoke = "F4",
                                calib_N = 100000L, plateau = "200:40"))
vs_threads(1L)
options(width = 160L)
only <- trimws(strsplit(args$only, ",", fixed = TRUE)[[1]])
cat(sprintf("wp4_sp: choicer %s from %s; sections %s\n", b$build$version, b$build$path,
            paste(only, collapse = ", ")))
t_start <- proc.time()[["elapsed"]]

verdicts <- logical(0)
report <- function(label, value, tol, pass = isTRUE(value < tol)) {
  cat(sprintf("%s  %-70s %11.4g  (threshold %.3g)\n", if (pass) "PASS" else "FAIL", label,
              value, tol))
  verdicts[length(verdicts) + 1L] <<- pass
}
quiet <- function(expr) suppressMessages(expr)
now <- function() proc.time()[["elapsed"]]
dgp <- get_dgp("D1")

# choicer's log-likelihood of D1 data at theta, every parameter fixed, on
# the draws cube eta (K_w x S x U; U persons with person_col, tasks without).
d1_fixed_ll <- function(data, theta, eta, person_col = NULL) {
  idat <- prepare_mxl_data(data, "id", "alt", "choice", SP_X, SP_W, include_outside_option = TRUE,
                           rc_correlation = TRUE, person_col = person_col)
  quiet(run_mxlogit(input_data = idat, eta_draws = eta, rc_dist = c(1L, 0L), rc_mean = TRUE,
                    rc_correlation = TRUE, include_outside_option = TRUE, theta_init = theta,
                    lower = theta, upper = theta))
}

# D1-L's probit choice probabilities for one task (rows alt, price, feature;
# the opt-out first), independently of the DGP's code: tastes b + Q
# sqrt(Lambda) eta (W = Q Lambda Q') on a 40 x 40 Gauss-Hermite rule, and
# given them P_j = E_e prod_{k != j} Phi(V_j + e - V_k) over the chosen
# row's shock (64 nodes), the opt-out's V = 0, every shock N(0, 1).
probit_ref <- function(hb, task, n_b = 40L, n_e = 64L) {
  gh <- function(n) {
    jac <- matrix(0, n, n)
    jac[cbind(1:(n - 1L), 2:n)] <- jac[cbind(2:n, 1:(n - 1L))] <- sqrt(1:(n - 1L) / 2)
    eg <- eigen(jac, symmetric = TRUE)
    list(z = sqrt(2) * eg$values, w = eg$vectors[1L, ]^2)
  }
  g <- gh(n_b)
  e <- gh(n_e)
  ev <- eigen(hb$W, symmetric = TRUE)
  B <- hb$b + ev$vectors %*% diag(sqrt(ev$values)) %*%
    rbind(rep(g$z, n_b), rep(g$z, each = n_b))
  wb <- rep(g$w, n_b) * rep(g$w, each = n_b)
  V <- rbind(0, hb$delta[task$alt] + task$price %o% B[1L, ] + task$feature %o% B[2L, ])
  J <- nrow(V)
  P <- matrix(0, J, ncol(V))
  for (j in seq_len(J)) for (q in seq_along(e$z)) {
    lp <- stats::pnorm(matrix(rep(V[j, ] + e$z[q], each = J - 1L), J - 1L) - V[-j, , drop = FALSE],
                       log.p = TRUE)
    P[j, ] <- P[j, ] + e$w[q] * exp(colSums(lp))
  }
  drop(P %*% wb)
}

# The log-likelihood of the probabilities a dataset was drawn from.
own_ll <- function(sim, d) {
  ch <- chosen_rows(sim$data$choice, d$sit)
  p <- sim$meta$probs
  sum(log(ifelse(ch > 0L, p$p[pmax(ch, 1L)], p$p0)))
}

# --- q0 ------------------------------------------------------------------------------

if ("q0" %in% only) {
  cat("\n== q0: D1 and D1-L against choicer ==\n")
  for (v in dgp$variants) {
    s1 <- dgp$simulate(N = 300L, seed = 11L, variant = v)
    s2 <- dgp$simulate(N = 300L, seed = 11L, variant = v)
    s3 <- dgp$simulate(N = 300L, seed = 12L, variant = v)
    report(sprintf("%s: simulate() twice at one seed identical, another seed not", v),
           as.numeric(!(identical(s1$data, s2$data) && identical(s1$truth, s2$truth) &&
                          !identical(s1$data$choice, s3$data$choice))), 0.5)
    lev <- SP_LEVELS[[v]]
    J <- length(SP_ALTS[[v]])
    off <- 0
    for (a in names(lev)) {
      tab <- s1$data[, .N, by = c("person", "alt", a)]
      off <- off + sum(tab$N != s1$meta$T / length(lev[[a]])) +
        abs(nrow(tab) - 300L * J * length(lev[[a]]))
    }
    report(sprintf("%s: every level T / levels times per person and position", v), off, 0.5)
    report(sprintf("%s: rows in choicer's order (person, id, alt) and one choice at most", v),
           as.numeric(is.unsorted(s1$data$id) || any(s1$data[, sum(choice), by = id][[2L]] > 1L)),
           0.5)
  }

  # D1: the panel and the cross-section on shared cubes, and the DGP's tastes.
  for (sd in 1:2) {
    sim <- dgp$simulate(N = 300L, seed = 20L + sd, keep_probs = TRUE)
    th <- sim$truth
    spec <- sim$meta$spec
    par <- unpack_theta("mxl", spec, th)
    d_p <- design(sim$data, spec, person = "person")
    d_x <- design(sim$data, spec)
    chosen <- chosen_rows(sim$data$choice, d_p$sit)
    eta_p <- vs_with_seed(100L + sd, array(stats::rnorm(2L * 50L * 300L), c(2L, 50L, 300L)))
    fp <- d1_fixed_ll(sim$data, th, eta_p, "person")
    report(sprintf("D1 seed %d: choicer's names = theta0's", sd),
           as.numeric(!identical(names(coef(fp)), names(th))), 0.5)
    ll <- mxl_loglik(d_p, par, eta_p, chosen)
    report(sprintf("D1 seed %d panel: simulated ll on a shared cube (S = 50), rel", sd),
           abs(ll - fp$loglik) / abs(ll), 1e-10)
    eta_x <- vs_with_seed(200L + sd, array(stats::rnorm(2L * 50L * 2400L), c(2L, 50L, 2400L)))
    ll <- mxl_loglik(d_x, par, eta_x, chosen)
    report(sprintf("D1 seed %d cross-section: simulated ll on a shared cube (S = 50), rel", sd),
           abs(ll - d1_fixed_ll(sim$data, th, eta_x)$loglik) / abs(ll), 1e-10)
    ll1 <- own_ll(sim, d_p)
    cl1 <- d1_fixed_ll(sim$data, th, array(sim$meta$eta, c(2L, 1L, 300L)), "person")$loglik
    report(sprintf("D1 seed %d: the DGP's own tastes as a one-draw cube per person, rel", sd),
           abs(ll1 - cl1) / abs(ll1), 1e-10)
  }

  # D1-L: an independent computation, then choicer's MXL with the same model.
  for (sd in 1:2) {
    sl <- dgp$simulate(N = 300L, seed = 30L + sd, variant = "labeled", keep_probs = TRUE)
    hb <- sl$meta$hb
    dt <- sl$data
    B <- t(hb$b + t(chol(hb$W)) %*% sl$meta$eta)
    V <- hb$delta[dt$alt] + dt$price * B[dt$person, 1L] + dt$feature * B[dt$person, 2L]
    den <- 1 + rowsum(exp(V), dt$id)[, 1L]
    diff <- max(abs(B - sl$meta$tastes), abs(exp(V) / den[as.character(dt$id)] - sl$meta$probs$p),
                abs(1 / den - sl$meta$probs$p0))
    report(sprintf("D1-L seed %d: probabilities against (b, W, delta) directly", sd), diff, 1e-12)
    th <- sl$meta$theta
    names(th)[1:2] <- SPL_X
    dw <- copy(dt)[, `:=`(price_w = price, feature_w = feature)]
    idat <- prepare_mxl_data(dw, "id", "alt", "choice", SPL_X, c("price_w", "feature_w"),
                             include_outside_option = TRUE, rc_correlation = TRUE,
                             person_col = "person")
    f <- quiet(run_mxlogit(input_data = idat, eta_draws = array(sl$meta$eta, c(2L, 1L, 300L)),
                           rc_dist = c(0L, 0L), rc_mean = FALSE, rc_correlation = TRUE,
                           include_outside_option = TRUE, theta_init = th, lower = th, upper = th))
    ll1 <- own_ll(sl, design(dt, sl$meta$spec, person = "person"))
    report(sprintf("D1-L seed %d: own tastes against choicer's MXL (one-draw cube), rel", sd),
           abs(ll1 - f$loglik) / abs(ll1), 1e-10)
  }
  sp <- dgp$simulate(N = 2000L, seed = 33L, variant = "labeled", shock = "probit")
  report("D1-L probit: at most one choice per task",
         as.numeric(any(sp$data[, sum(choice), by = id][[2L]] > 1L)), 0.5)
  cat(sprintf("      D1-L probit at theta0 (calibrated for the logit): opt-out share %.3f\n",
              dgp$moments(sp)[["outside_share"]]))
  # One fixed task for 10^6 persons, choices drawn the DGP's way (tastes by
  # rc_coefs(), then .sp_probit_choice()), against the quadrature.
  tr <- sp_truth("labeled")
  spec <- sp_spec("labeled")
  par <- unpack_theta("mxl", spec, tr$theta)
  task <- data.table(alt = SP_ALTS$labeled, price = c(1, 2, 3, 4, 2, 3),
                     feature = c(1, 0, 1, 0, 1, 0))
  R <- 1000000L
  big <- task[rep(seq_len(6L), R)]
  big[, person := rep(seq_len(R), each = 6L)]
  big[, id := person]
  t0 <- now()
  ch <- vs_with_seed(61L, {
    eta <- matrix(stats::rnorm(2L * R), 2L)
    .sp_probit_choice(design(big, spec, person = "person"), par,
                      t(rc_coefs(par$mu, par$L, eta, par$rc_dist, par$rc_mean)))
  })
  obs <- tabulate(ifelse(ch$chosen > 0L, (ch$chosen - 1L) %% 6L + 2L, 1L), 7L)
  ref <- probit_ref(tr$hb, task)
  stat <- sum((obs - R * ref)^2 / (R * ref))
  pv <- stats::pchisq(stat, 6, lower.tail = FALSE)
  cat(sprintf("      quadrature: sum - 1 = %.1e; opt-out, P1..P6: %s; %.1f s\n", sum(ref) - 1,
              paste(sprintf("%.4f", ref), collapse = " "), now() - t0))
  report(sprintf("D1-L probit, one task x 10^6 persons: chi2 %.2f, df 6, p-value", stat), pv,
         1e-4, pv >= 1e-4)
}

# --- calib -----------------------------------------------------------------------------

if ("calib" %in% only) {
  cat(sprintf("\n== calib: moments at N = %s persons x 8 tasks (seed 1) against the targets ==\n",
              format(args$calib_N, big.mark = ",")))
  tg <- dgp$targets
  for (v in dgp$variants) {
    t0 <- now()
    sim <- dgp$simulate(N = args$calib_N, seed = 1L, variant = v)
    m <- dgp$moments(sim)
    cat(sprintf("-- variant %s, simulated in %.1f s\n", v, now() - t0))
    i <- match(names(m), tg$moment)
    tab <- data.frame(moment = names(m), value = signif(unname(m), 4), lo = tg$lo[i], hi = tg$hi[i],
                      within = ifelse(is.na(i), NA, m >= tg$lo[i] & m <= tg$hi[i]))
    if (v == "default") {
      a <- sp_coef_moments()
      tab$analytic <- signif(unname(a[tab$moment]), 4)
    }
    print(tab, row.names = FALSE)
    for (k in intersect(tg$moment, names(m))) {
      j <- tg$moment == k
      report(sprintf("%s: %s in [%.2f, %.2f]", v, k, tg$lo[j], tg$hi[j]), m[[k]], tg$hi[j],
             m[[k]] >= tg$lo[j] && m[[k]] <= tg$hi[j])
    }
  }
  cat("\ntheta0 (default):\n")
  print(dgp$theta0$default, digits = 6)
  cat("theta0 (labeled):\n")
  print(dgp$theta0$labeled, digits = 6)
  cat("\nsqrt(U)/S per F4 cell (U: persons for panel, tasks for xsec):\n")
  fam <- load_family("F4", b$L)
  cl <- fam$cells
  print(data.frame(cell = cl$cell, persons = cl$N, tasks = cl$N * cl$T, S = cl$S,
                   panel = round(sqrt(cl$N) / cl$S, 3), S_xsec = cl$S_xsec,
                   xsec = round(sqrt(cl$N * cl$T) / cl$S_xsec, 3),
                   xsec_S800 = round(sqrt(cl$N * cl$T) / 800, 3),
                   gated = cl$cell %in% fam$gates$cell), row.names = FALSE)
}

# --- ident -----------------------------------------------------------------------------

# The implied SDs of a panel fit from its sigma block, through Sigma's vech
# (sigma_vech(): Sigma_11, Sigma_21, Sigma_22): sd(g_1), sd(g_2),
# corr(g_1, g_2), and the coefficients' SDs sd(b_t) (shifted log-normal) and
# sd(b_q) = sd(g_2).
implied_sd <- function(s) {
  v <- sigma_vech(s, 2L, TRUE)$est
  c(sd_g1 = sqrt(v[1]), sd_g2 = sqrt(v[3]), corr_g = v[2] / sqrt(v[1] * v[3]),
    sd_b_time = sqrt((exp(v[1]) - 1) * exp(v[1])), sd_b_quality = sqrt(v[3]))
}

if ("ident" %in% only) {
  cat("\n== ident: sigma / SE(sigma) from panel fits (S = 500) ==\n")
  sg <- c("L_11", "L_21", "L_22")
  truth_sd <- implied_sd(dgp$theta0$default[sg])
  fit_at <- function(N, seed) {
    sim <- dgp$simulate(N = N, seed = seed)
    t0 <- now()
    fit <- quiet(run_mxlogit(sim$data, "id", "alt", "choice", SP_X, SP_W, rc_dist = c(1L, 0L),
                             rc_mean = TRUE, rc_correlation = TRUE, include_outside_option = TRUE,
                             person_col = "person", S = 500L, draws = "generate",
                             seed = seed + 1L, polish = "newton"))
    s <- coef(fit)[sg]
    G <- numDeriv::jacobian(function(x) implied_sd(x), unname(s))
    imp <- implied_sd(s)
    list(N = N, seed = seed, sec = now() - t0, L = s, L_se = fit$se[sg], imp = imp,
         imp_se = stats::setNames(sqrt(diag(G %*% fit$vcov[sg, sg] %*% t(G))), names(imp)))
  }
  for (r in list(fit_at(200L, 501L), fit_at(800L, 502L))) {
    cat(sprintf("-- %d persons (seed %d), %.1f s\n", r$N, r$seed, r$sec))
    est <- c(r$L, r$imp)
    se <- c(r$L_se, r$imp_se)
    print(data.frame(param = names(est), truth = round(c(dgp$theta0$default[sg], truth_sd), 4),
                     est = round(est, 4), se = round(se, 4), est_over_se = round(est / se, 2)),
          row.names = FALSE)
    if (r$N == 200L) {
      for (k in c("sd_g1", "sd_g2")) {
        report(sprintf("200 persons: %s / SE (the plan requires >= 5 for D4 only)", k),
               r$imp[[k]] / r$imp_se[[k]], 5, r$imp[[k]] / r$imp_se[[k]] >= 5)
      }
    }
  }
  more <- lapply(503:507, function(s) fit_at(200L, s))
  ratio <- t(vapply(more, function(r) c(seed = r$seed, r$L, r$imp / r$imp_se), numeric(9)))
  cat("      five more fits at 200 persons: L, and sigma / SE(sigma)\n")
  print(round(ratio, 2))
  pl <- ratio[, "L_11"] < -5 | ratio[, "L_22"] < -5
  cat(sprintf("      median over the %d fits off the plateau (no L diagonal below -5):\n",
              sum(!pl)))
  print(round(apply(ratio[!pl, -(1:4), drop = FALSE], 2, stats::median), 2))
}

# --- hb --------------------------------------------------------------------------------

if ("hb" %in% only) {
  cat("\n== hb: D1-L through run_hmnlogit() and run_hmnprobit(), 100 persons ==\n")
  run_hb <- function(fn, prep, sim) {
    w <- character(0)
    t0 <- now()
    out <- tryCatch(withCallingHandlers({
      inp <- quiet(prep(sim$data, "id", "alt", "choice", SPL_X, person_col = "person",
                        alt_covariate_cols = "z"))
      quiet(fn(input_data = inp, mcmc = list(R = 600L, burn = 300L, seed = 1L), chains = 1))
    }, warning = function(c) {
      w <<- c(w, substr(conditionMessage(c), 1, 90))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    list(fit = out, warnings = w, sec = now() - t0)
  }
  hb_table <- function(fit, sim) {
    tr <- sim$truth
    dr <- fit$draws
    est <- c(colMeans(dr$b), colMeans(dr$w_vech), colMeans(dr$theta),
             sigma_d2 = mean(dr$sigma_d2),
             stats::setNames(colMeans(dr$delta), paste0("delta_", colnames(dr$delta))))
    data.frame(param = names(tr), truth = round(unname(tr), 3),
               post_mean = round(unname(est[names(tr)]), 3))
  }
  sl <- dgp$simulate(N = 100L, seed = 41L, variant = "labeled")
  sp <- dgp$simulate(N = 100L, seed = 42L, variant = "labeled", shock = "probit")
  for (r in list(list("HMNL", run_hmnlogit, prepare_hmnl_data, sl),
                 list("HMNP", run_hmnprobit, prepare_hmnp_data, sp))) {
    res <- run_hb(r[[2]], r[[3]], r[[4]])
    ok <- !inherits(res$fit, "error")
    report(sprintf("%s on D1-L (%s data), R = 600, burn 300: runs (%.1f s)", r[[1]],
                   r[[4]]$meta$shock, res$sec), as.numeric(!ok), 0.5)
    if (!ok) {
      cat("      error:", conditionMessage(res$fit), "\n")
      next
    }
    for (x in unique(res$warnings)) cat("      warning:", x, "\n")
    print(hb_table(res$fit, r[[4]]), row.names = FALSE)
  }
}

# --- smoke -----------------------------------------------------------------------------

if ("smoke" %in% only) {
  cat("\n== smoke: run_rep() at every cell x method ==\n")
  tok <- trimws(strsplit(args$smoke, ",", fixed = TRUE)[[1]])
  tab <- list()
  for (tk in tok) {
    p <- strsplit(tk, ":", fixed = TRUE)[[1]]
    fam <- load_family(p[1], b$L)
    cells <- if (length(p) >= 2L) p[2] else fam$cells$cell
    for (cl in cells) {
      ms <- cell_methods(fam, cl)
      if (length(p) >= 3L) ms <- intersect(ms, p[3])
      reps <- if (cl == "200") c(1L, 2L, 10L) else 1L
      for (r in reps) {
        recs <- run_rep(fam, cl, r, ms)
        for (rec in recs) {
          valid <- tryCatch({
            do.call(new_rep_record, rec[c("family", "cell", "method", "rep", "seed", "ok", "error",
                                          "warnings", "est", "vcov_c", "joint_wald", "lr",
                                          "derived", "diag", "meta")])
            TRUE
          }, error = function(e) FALSE)
          w <- unique(rec$warnings$class)
          jw <- rec$joint_wald
          if (!is.null(jw) && is.null(jw$p)) jw <- jw[[1L]]
          tab[[length(tab) + 1L]] <- data.frame(
            family = rec$family, cell = rec$cell, method = rec$method, rep = rec$rep, ok = rec$ok,
            valid = valid, sim_sec = round(rec$diag$sim_sec, 2),
            unit_sec = round(rec$diag$unit_sec, 1),
            fit_sec = round(rec$diag$time_sec %||% NA_real_, 1),
            iterations = rec$diag$iterations %||% NA, code = rec$diag$code %||% NA,
            converged = rec$diag$converged %||% NA,
            decrement = signif(rec$diag$decrement %||% NA_real_, 2),
            lr = signif(rec$lr$LR %||% NA_real_, 3),
            lr_fit_sec = round(rec$diag$lr_restricted[["time_sec"]] %||% NA_real_, 1),
            ms_dlogL = signif(rec$diag$multistart$dlogL %||% NA_real_, 2),
            joint_p = signif(jw$p %||% NA_real_, 2),
            warnings = if (length(w)) paste(w, collapse = " ") else "",
            error = substr(rec$error %||% "", 1, 60), stringsAsFactors = FALSE)
          report(sprintf("%s/%s/%s rep %d: ok and valid", rec$family, rec$cell, rec$method,
                         rec$rep), as.numeric(!(rec$ok && valid)), 0.5)
          if (rec$ok && rec$rep == 1L) {
            e <- rec$est
            se <- if (rec$method == "xsec") e$se_cluster else e$se_hessian
            cat(sprintf("      z (%s SE): %s\n", if (rec$method == "xsec") "cluster" else "Hessian",
                        paste(sprintf("%s %.2f", e$name, (e$est - e$truth) / se), collapse = ", ")))
            cat(sprintf("      derived: %s\n", paste(sprintf(
              "%s %.4g (se %.2g, truth %.4g)", rec$derived$name, rec$derived$est, rec$derived$se,
              rec$derived$truth), collapse = "; ")))
          }
        }
      }
    }
  }
  print(do.call(rbind, tab), row.names = FALSE)
}

# --- plateau ---------------------------------------------------------------------------

# Reps 1..K of F4's methods per cell (--plateau=<cell>:<K>,...): failures,
# fits on a plateau (an L diagonal below -5: L_11, the time coefficient's SD
# to 0; L_22, Sigma singular, |corr(g_1, g_2)| -> 1), and, for each plateau
# fit, the default fit's log-likelihood against refits started at theta0
# and with scale_vars = "none" (dlogL > 0: the plateau is not the maximum).
if ("plateau" %in% only) {
  cat("\n== plateau: F4's reps, failures and plateaus (an L diagonal below -5) ==\n")
  fam <- load_family("F4", b$L)
  fenv <- environment(fam$methods$panel)
  rows <- list()
  rescue <- list()
  for (tk in trimws(strsplit(args$plateau, ",", fixed = TRUE)[[1]])) {
    p <- strsplit(tk, ":", fixed = TRUE)[[1]]
    cl <- p[1]
    cell <- family_cell(fam, cl)
    t0 <- now()
    recs <- unlist(lapply(seq_len(as.integer(p[2])),
                          function(r) run_rep(fam, cl, r, c("panel", "xsec"))), recursive = FALSE)
    for (m in c("panel", "xsec")) {
      rr <- Filter(function(x) x$method == m, recs)
      ok <- vapply(rr, function(x) x$ok, logical(1))
      pl <- vapply(rr, function(x) isTRUE(x$diag$plateau), logical(1))
      rows[[length(rows) + 1L]] <- data.frame(
        cell = cl, method = m, reps = length(rr), ok = sum(ok), plateau = sum(pl),
        median_sec = round(stats::median(vapply(rr, function(x) x$diag$unit_sec, numeric(1))), 1),
        mean_lr = if (m == "panel") round(mean(unlist(lapply(rr, function(x) x$lr$LR))), 3) else NA,
        errors = paste(unique(substr(unlist(lapply(rr[!ok], function(x) x$error)), 1, 50)),
                       collapse = " | "), stringsAsFactors = FALSE)
      for (x in rr[pl]) {
        sim <- fam$simulate(cell, x$seed)
        S <- if (m == "panel") cell$S else cell$S_xsec
        extra <- if (m == "panel") list(person_col = "person") else list(cluster_col = "person")
        fit_m <- function(...) {
          tryCatch(quiet(do.call(fenv$f4_fit, c(list(sim, x$seed, S, "rescue"), extra, list(...)))),
                   error = function(e) NULL)
        }
        f0 <- fit_m()
        f_t <- fit_m(theta_init = sim$truth)
        f_n <- fit_m(scale_vars = "none")
        ll <- function(f) if (is.null(f)) NA_real_ else f$loglik
        rescue[[length(rescue) + 1L]] <- data.frame(
          cell = cl, method = m, rep = x$rep, min_L_diag = round(x$diag$min_L_diag, 1),
          dll_theta0_start = round(ll(f_t) - ll(f0), 3),
          dll_scale_none = round(ll(f_n) - ll(f0), 3),
          L11_theta0_start = round(if (is.null(f_t)) NA else coef(f_t)[["L_11"]], 3),
          L22_theta0_start = round(if (is.null(f_t)) NA else coef(f_t)[["L_22"]], 3))
      }
    }
    cat(sprintf("      cell %s: %.0f s\n", cl, now() - t0))
  }
  print(do.call(rbind, rows), row.names = FALSE)
  if (length(rescue)) {
    cat("plateau fits against refits (dlogL > 0: the plateau is not the maximum):\n")
    print(do.call(rbind, rescue), row.names = FALSE)
  }
}

# --- xtime -----------------------------------------------------------------------------

if ("xtime" %in% only) {
  cat("\n== xtime: the composite fit at 3.2k persons (25,600 tasks), S = 200 and 800 ==\n")
  fam <- load_family("F4", b$L)
  cell <- family_cell(fam, "3.2k")
  seed <- vs_seed("F4", "3.2k", 1L)
  sim <- fam$simulate(cell, seed)
  fits <- list()
  for (S in c(200L, 800L)) {
    t0 <- now()
    fits[[as.character(S)]] <- quiet(f <- run_mxlogit(
      sim$data, "id", "alt", "choice", SP_X, SP_W, rc_dist = c(1L, 0L), rc_mean = TRUE,
      rc_correlation = TRUE, include_outside_option = TRUE, S = S, draws = "generate",
      seed = seed + 1L, polish = "newton", cluster_col = "person"))
    cat(sprintf(paste("      S = %d (sqrt(U)/S = %.2f): %.1f s wall, fit %.1f s, %d iterations,",
                      "certified %s, decrement %.1e, min L diagonal %.2f\n"),
                S, sqrt(25600) / S, now() - t0, f$optimizer$elapsed_time, f$optimizer$iterations,
                f$optimizer$polish$converged, f$optimizer$polish$decrement,
                min(coef(f)[c("L_11", "L_22")])))
  }
  V200 <- stats::vcov(fits[["200"]], type = "hessian")
  gap <- (coef(fits[["200"]]) - coef(fits[["800"]])) / sqrt(diag(V200))
  cat("      (theta_S200 - theta_S800) / Hessian SE at S = 200:\n")
  print(round(gap, 3))
}

cat(sprintf("\n%d checks, %d FAIL; %.1f s\n", length(verdicts), sum(!verdicts),
            proc.time()[["elapsed"]] - t_start))
if (any(!verdicts)) quit(status = 1L)
