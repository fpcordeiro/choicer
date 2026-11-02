# ============================================================================
# WP5's checks: D4's designs (lib/dgp_geo.R; --dgp=D4-40, the default and
# F2/F3's, or D4, D4-small), families F2 and F3, and the Q1.3 reference
# (lib/references_quad.R). Sections (--only, comma-separated):
#   q0     the probabilities D4 draws choices from at theta0 against
#          choicer: MNL and NL against predict() (max abs diff < 1e-12);
#          MXL, the suite's simulated log-likelihood on a shared cube and on
#          the DGP's own tastes (a one-draw cube) against choicer's
#          fixed-parameter fit (< 1e-10 relative); the draws' order and
#          reproducibility.
#   calib  the calibration moments at N = 200,000 against plan §4's targets
#          (per variant), and per rung of F2/F3: sqrt(N)/S, J/sqrt(N), the
#          expected counts per school and the expected number never chosen.
#   ident  identification strength sigma_q / SE(sigma_q) at N = 1000: the
#          16k fit's SE scaled by 4, and the median over five 1k fits
#          (target >= 5, plan §4).
#   drop   drop_unbounded_alts() on a constructed case.
#   q13    ref_q13() on D4-40 (--q13_N, --q13_S).
#   smoke  F2 and F3, every cell x method (--smoke=F2,F3 or F3:16k:blp, ...)
#          through run_rep(): R = 1, R = 2 at the bottom rung plus rep 10
#          there (the multi-start); each record re-validated.
# Usage (the dev build first on R_LIBS; one thread; each section < 10 min):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 \
#     Rscript _validation/checks/wp5_geo.R --only=q0,calib,drop
# Exits with status 1 on any FAIL.
# ============================================================================

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg))), ".."))
sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
b <- boot$vs_boot(root)
suppressPackageStartupMessages(library(data.table))
args <- vs_args(defaults = list(only = "q0,calib,ident,drop,q13,smoke", smoke = "F2,F3",
                                dgp = "D4-40", q13_N = 2000L, q13_S = 4000L))
vs_threads(1L)
options(width = 160L)
only <- trimws(strsplit(args$only, ",", fixed = TRUE)[[1]])
cat(sprintf("wp5_geo: choicer %s from %s; sections %s\n", b$build$version, b$build$path,
            paste(only, collapse = ", ")))
t_start <- proc.time()[["elapsed"]]

verdicts <- logical(0)
report <- function(label, value, tol, pass = isTRUE(value < tol)) {
  cat(sprintf("%s  %-66s %11.4g  (threshold %.3g)\n", if (pass) "PASS" else "FAIL", label,
              value, tol))
  verdicts[length(verdicts) + 1L] <<- pass
}
quiet <- function(expr) {
  withCallingHandlers(suppressMessages(expr),
                      choicer_unbounded_constants = function(w) invokeRestart("muffleWarning"))
}
dgp <- get_dgp(args$dgp)
cat(sprintf("DGP %s: %s\n", dgp$id, dgp$name))

# --- q0 ------------------------------------------------------------------------------

if ("q0" %in% only) {
  cat("\n== q0: D4's probabilities at theta0 against choicer ==\n")
  # Reproducibility, and the choice uniforms draw_choice() is handed.
  s1 <- dgp$simulate(N = 500L, seed = 11L, variant = "nl")
  s2 <- dgp$simulate(N = 500L, seed = 11L, variant = "nl")
  report("simulate() twice at one seed: identical data and truth",
         as.numeric(!(identical(s1$data, s2$data) && identical(s1$truth, s2$truth))), 0.5)
  dr <- geo_draw(dgp$design, 500L, 11L)
  last <- vs_with_seed(11L, stats::runif(3L * 500L + 500L + 2L * 500L + 500L))
  report("the choice uniforms are the sample stream's last N",
         as.numeric(!identical(utils::tail(last, 500L), dr$sample$choice)), 0.5)
  report("the world does not depend on N or the seed",
         as.numeric(!identical(geo_draw(dgp$design, 10L, 1L)$world, dr$world)), 0.5)
  for (v in c("mnl", "nl")) {
    for (sd in 1:2) {
      sim <- dgp$simulate(N = 3000L, seed = sd, variant = v, keep_probs = TRUE)
      fit <- quiet(if (v == "mnl") {
        run_mnlogit(sim$data, "id", "alt", "choice", D4_X, include_outside_option = TRUE)
      } else {
        run_nestlogit(sim$data, "id", "alt", "choice", D4_X, nest_col = "sector",
                      include_outside_option = TRUE)
      })
      ok_names <- identical(names(coef(fit)), names(sim$truth))
      fit$coefficients[names(sim$truth)] <- sim$truth
      if (v == "nl") fit$lambda <- fit$coefficients[fit$param_map$lambda]
      pr <- predict(fit, type = "probabilities")
      diff <- max(abs(pr$choice_prob[, 1] - sim$meta$probs$p),
                  abs(pr$choice_prob_outside[, 1] - sim$meta$probs$p0))
      report(sprintf("%s seed %d: names = choicer's", toupper(v), sd), as.numeric(!ok_names), 0.5)
      report(sprintf("%s seed %d (%d rows): max|p - predict()| at theta0", toupper(v), sd,
                     nrow(sim$data)), diff, 1e-12)
    }
  }
  for (sd in 1:2) {
    sim <- dgp$simulate(N = 1000L, seed = sd, variant = "default", keep_probs = TRUE)
    th <- sim$truth
    spec <- sim$meta$spec
    d <- design(sim$data, spec)
    chosen <- chosen_rows(sim$data$choice, d$sit)
    idat <- prepare_mxl_data(sim$data, "id", "alt", "choice", D4_X, D4_W,
                             include_outside_option = TRUE, rc_correlation = FALSE)
    fixed_ll <- function(eta) {
      quiet(run_mxlogit(input_data = idat, eta_draws = eta, rc_dist = c(1L, 0L), rc_mean = FALSE,
                        rc_correlation = FALSE, include_outside_option = TRUE, theta_init = th,
                        lower = th, upper = th))$loglik
    }
    eta <- vs_with_seed(100L + sd, array(stats::rnorm(2L * 50L * 1000L), c(2L, 50L, 1000L)))
    ll <- mxl_loglik(d, unpack_theta("mxl", spec, th), eta, chosen)
    cl <- fixed_ll(eta)
    report(sprintf("MXL seed %d: simulated ll on a shared cube (S = 50), rel", sd),
           abs(ll - cl) / abs(ll), 1e-10)
    eta1 <- array(sim$meta$eta, c(2L, 1L, 1000L))
    p <- sim$meta$probs
    ll1 <- sum(log(ifelse(chosen > 0L, p$p[pmax(chosen, 1L)], p$p0)))
    report(sprintf("MXL seed %d: the DGP's own tastes as a one-draw cube, rel", sd),
           abs(ll1 - fixed_ll(eta1)) / abs(ll1), 1e-10)
  }
}

# --- calib -----------------------------------------------------------------------------

if ("calib" %in% only) {
  cat("\n== calib: moments at N = 200,000 (seed 1) against the targets ==\n")
  tg <- dgp$targets
  for (v in dgp$variants) {
    t0 <- proc.time()[["elapsed"]]
    sim <- dgp$simulate(N = 200000L, seed = 1L, variant = v)
    m <- dgp$moments(sim)
    cat(sprintf("-- variant %s (%s), simulated in %.1f s\n", v, D4_VARIANTS[[v]],
                proc.time()[["elapsed"]] - t0))
    i <- match(names(m), tg$moment)
    print(data.frame(moment = names(m), value = signif(unname(m), 4), lo = tg$lo[i], hi = tg$hi[i],
                     within = ifelse(is.na(i), NA, m >= tg$lo[i] & m <= tg$hi[i])),
          row.names = FALSE)
    for (k in tg$moment) {
      report(sprintf("%s: %s in [%.2f, %.2f]", v, k, tg$lo[tg$moment == k], tg$hi[tg$moment == k]),
             m[[k]], tg$hi[tg$moment == k], m[[k]] >= tg$lo[tg$moment == k] &&
               m[[k]] <= tg$hi[tg$moment == k])
    }
    # Per rung: expected counts per school, expected number never chosen
    # (sum_j exp(-N s_j)), F2's gated constants, sqrt(N)/S and J/sqrt(N).
    sh <- sim$meta$shares
    S_rung <- c("1000" = 500, "4000" = 500, "16000" = 520)       # F3's cells
    rung <- do.call(rbind, lapply(c(1000, 4000, 16000), function(N) {
      e <- N * sh
      data.frame(variant = v, N = N, J = length(sh), J_over_sqrtN = round(length(sh) / sqrt(N), 3),
                 S = S_rung[[as.character(N)]],
                 sqrtN_over_S = round(sqrt(N) / S_rung[[as.character(N)]], 3),
                 count_min = round(min(e), 1),
                 count_p10 = round(unname(stats::quantile(e, 0.1)), 1),
                 count_med = round(stats::median(e), 1),
                 exp_never = round(sum(exp(-e)), 2), n_count_ge50 = sum(e >= 50))
    }))
    print(rung, row.names = FALSE)
  }
  print(dgp$theta0$default[1:5])
  print(dgp$theta0$nl[1:6])
}

# --- ident -----------------------------------------------------------------------------

if ("ident" %in% only) {
  cat("\n== ident: sigma_q / SE(sigma_q) at N = 1000 (target >= 5) ==\n")
  fit_at <- function(N, seed, S) {
    sim <- geo_drop(dgp$simulate(N = N, seed = seed, variant = "default"))
    t0 <- proc.time()[["elapsed"]]
    fit <- suppressMessages(run_mxlogit(sim$data, "id", "alt", "choice", D4_X, D4_W,
                                        rc_dist = c(1L, 0L), rc_mean = FALSE,
                                        include_outside_option = TRUE, S = S, draws = "generate",
                                        seed = seed + 1L))
    c(N = N, seed = seed, S = S, sec = proc.time()[["elapsed"]] - t0,
      L_22 = coef(fit)[["L_22"]], se_22 = fit$se[["L_22"]], L_11 = coef(fit)[["L_11"]],
      se_11 = fit$se[["L_11"]], dropped = sim$meta$drop[["alts"]])
  }
  r1 <- as.data.frame(do.call(rbind, lapply(1:5, function(s) fit_at(1000L, 900L + s, 500L))))
  rb <- fit_at(16000L, 907L, 200L)
  print(rbind(round(r1, 3), round(rb, 3)), row.names = FALSE)
  th0 <- dgp$theta0$default
  cat(sprintf("truth: L_11 = log(sigma_d) = %.4f, L_22 = log(sigma_q) = %.4f\n",
              th0[["L_11"]], th0[["L_22"]]))
  # sigma / SE(sigma) = 1 / SE(log sigma) (delta method)
  s16 <- 1 / (rb[["se_22"]] * 4)
  report("sigma_q/SE(sigma_q) at N = 1000, from the 16k fit (S = 200) scaled by 4", s16, 5,
         s16 >= 5)
  s1k <- stats::median(1 / r1$se_22)
  report("sigma_q/SE(sigma_q) at N = 1000, median of five 1k fits (S = 500)", s1k, 5, s1k >= 5)
  cat(sprintf("      sigma_d/SE(sigma_d) at N = 1000: %.2f (16k scaled), %.2f (median of 1k)\n",
              1 / (rb[["se_11"]] * 4), stats::median(1 / r1$se_11)))
}

# --- drop ------------------------------------------------------------------------------

if ("drop" %in% only) {
  cat("\n== drop: drop_unbounded_alts() on a constructed case ==\n")
  # Round 1: A is never chosen (its rows go, and situation 8, which offered
  # only A, is left without rows: its student chose the outside option); C
  # is chosen wherever offered (situations 2 and 3 go). Round 2: D, offered
  # in 2 and 4, is now offered only in 4, where it is chosen: 4 goes. Round
  # 3: B (chosen in 1 and 6) and F (chosen in 5) are stable.
  toy <- data.table(
    id     = c(1, 1, 1,  2, 2, 2,  3, 3,  4, 4, 4,  5, 5, 5,  6, 6,  7, 7,  8),
    alt    = c("A", "B", "F",  "B", "C", "D",  "B", "C",  "B", "D", "F",  "A", "B", "F",
               "B", "F",  "B", "F",  "A"),
    choice = c(0, 1, 0,  0, 1, 0,  0, 1,  0, 1, 0,  0, 0, 1,  1, 0,  0, 0,  0))
  r <- drop_unbounded_alts(toy)
  print(r$counts)
  ok <- identical(sort(r$dropped_alts), c("A", "C", "D")) &&
    identical(sort(r$dropped_ids), c(2, 3, 4, 8)) && identical(unique(r$data$id), c(1, 5, 6, 7)) &&
    identical(r$data$alt, rep(c("B", "F"), 4L)) && r$counts[["rounds"]] == 2L &&
    r$counts[["rows"]] == 11L
  report("dropped A, C, D and situations 2, 3, 4, 8 in two rounds", as.numeric(!ok), 0.5)
  sim <- geo_drop(dgp$simulate(N = 1000L, seed = 3L, variant = "nl"))
  st <- sim$data[, .(n = .N, ch = sum(choice)), by = alt]
  report(sprintf("D4-nl N = 1000 seed 3: no alternative never or always chosen (%d dropped)",
                 sim$meta$drop[["alts"]]), as.numeric(any(st$ch == 0 | st$ch == st$n)), 0.5)
  fit <- withCallingHandlers(suppressMessages(run_nestlogit(
    sim$data, "id", "alt", "choice", D4_X, nest_col = "sector", include_outside_option = TRUE)),
    choicer_unbounded_constants = function(w) stop("choicer still warns"))
  report("geo_truth(): theta0's names equal the fit's after the drop",
         as.numeric(!identical(names(coef(fit)), names(sim$truth))), 0.5)
}

# --- q13 -------------------------------------------------------------------------------

if ("q13" %in% only) {
  cat("\n== q13: MXL MSL against the quadrature ML ==\n")
  t0 <- proc.time()[["elapsed"]]
  rows <- ref_q13(N = args$q13_N, S = args$q13_S, verbose = TRUE)
  print(rows, row.names = FALSE)
  for (i in which(!is.na(rows$pass))) {
    report(rows$check[i], rows$value[i], rows$threshold[i] %||% NA, rows$pass[i])
  }
  cat(sprintf("      ref_q13() took %.1f s\n", proc.time()[["elapsed"]] - t0))
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
      bottom <- cl == fam$cells$cell[which.min(fam$cells$N)]
      reps <- if (bottom) c(1L, 2L, 10L) else 1L
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
          tab[[length(tab) + 1L]] <- data.frame(
            family = rec$family, cell = rec$cell, method = rec$method, rep = rec$rep, ok = rec$ok,
            valid = valid, sim_sec = round(rec$diag$sim_sec, 1),
            seconds = round(rec$diag$unit_sec, 1),
            fit_sec = round(rec$diag$time_sec %||% NA_real_, 1),
            iterations = rec$diag$iterations %||% NA, code = rec$diag$code %||% NA,
            converged = rec$diag$converged %||% NA,
            decrement = signif(rec$diag$decrement %||% NA_real_, 2),
            alts_dropped = rec$diag$drop[["alts"]] %||% NA,
            ms_dlogL = signif(rec$diag$multistart$dlogL %||% NA_real_, 2),
            joint_p = signif(rec$joint_wald$p %||% NA_real_, 2),
            warnings = if (length(w)) paste(w, collapse = " ") else "",
            error = substr(rec$error %||% "", 1, 60), stringsAsFactors = FALSE)
          report(sprintf("%s/%s/%s rep %d: ok and valid", rec$family, rec$cell, rec$method,
                         rec$rep), as.numeric(!(rec$ok && valid)), 0.5)
          if (rec$ok && rec$rep == 1L) {
            e <- rec$est
            z <- (e$est - e$truth) / e$se_hessian
            cat(sprintf("      common z: %s\n", paste(sprintf("%s %.2f", e$name[e$block != "asc"],
                                                              z[e$block != "asc"]),
                                                      collapse = ", ")))
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

cat(sprintf("\n%d checks, %d FAIL; %.1f s\n", length(verdicts), sum(!verdicts),
            proc.time()[["elapsed"]] - t_start))
if (any(!verdicts)) quit(status = 1L)
