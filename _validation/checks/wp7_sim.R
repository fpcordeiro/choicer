# ============================================================================
# WP7's checks: families F5 (simulation error) and F10 (the share-matching
# constants, the port of blp_concentration_mc.R). Sections (--only,
# comma-separated):
#   load    both families load and validate; F5's regimes and their sources.
#   rates   F5's rate uncertainty on synthetic records, no fits (under a minute):
#           b_se and slope_se over complete units on a correlated ladder, against
#           the diagonal formula on independent S points, with one dataset, with
#           source rows on the ladder, and calibrated over ensembles.
#   shares  F10's frozen population shares against the old script's
#           pop_shares() and the port's f10_pop_shares() (about 25 s a
#           design; --designs=).
#   repro   the port against the old script rep by rep (--repro=cross:1-3,
#           panel:1, ...): the old one_rep() (its functions parsed from the
#           script in git history, commit ae38e12, nothing else run) and the
#           port's run_rep() on the same seeds; where the default start won
#           fit_two_starts(), the estimates and SEs must be the old ones to
#           1e-12 (a tie within fit_two_starts()'s tol keeps the default
#           start's fit); where the theta0 start won, the log-likelihood gain
#           and the shift in SE units are printed (the rule changed the rep).
#   smoke   F10 through run_rep() at R = 2 (--smoke_reps=1:2), every cell and
#           method, and rep 10 of cross (the multi-start check); each record
#           re-validated; metrics_cell() on each unit; f10_findings() on the
#           records; a table of times and the mc budget.
#   spec    F5's models against the source families' own fits, every
#           regime (--spec=F3,F4,...): the source's method on a small cell of
#           the same model
#           (F3 1k.q, F4 200.q, F6 5k; F7's claims and census at N /
#           --spec_div, default 100), rep 1, through run_unit(), and F5's
#           model with the same S, draws seed and start (F7's: fit_newton()
#           from theta0); estimates equal.
#   f5      one dataset per regime (--f5=F4,F3,F6,F7.claims,F7.census)
#           through run_unit() at a reduced ladder (--f5_ref, --f5_gen,
#           --f5_store, --f5_seeds; the source's fit at its own S is full
#           size), F7's regimes at N / --f5_div (default 100, as F7's own
#           smoke divides; checks/wp6_geo.R); the record re-validated,
#           f5_findings() on it, the per-fit times and their extrapolation
#           to the regime's full ladder (linear in draws, and in rows for
#           F7: the warm fits' seconds per draw times the full S, times the
#           divisor).
# F6 and F7 simulate through the port (lib/geo_gen.cpp): set VS_GEO_CACHE to
# a build cache of the caller's own (geo_gen_compile()).
# Usage (the dev build first on R_LIBS; one thread; each section < 10 min):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 \
#     Rscript _validation/checks/wp7_sim.R --only=load,repro [--out=DIR]
# --out (default: no files) saves each section's records and tables as .rds.
# Exits with status 1 on any FAIL.
# ============================================================================

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg))), ".."))
sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
b <- boot$vs_boot(root)
suppressPackageStartupMessages(library(data.table))
args <- vs_args(defaults = list(only = "load,rates,repro", designs = "cross,panel,base,markets",
                                repro = "cross:1-3,panel:1,base:1,markets:1", smoke_reps = "1:2",
                                f5 = "F4", f5_ref = 400L, f5_gen = "50,100", f5_store = "100",
                                f5_seeds = 2L, f5_div = 100, spec_div = 100,
                                spec = "F3,F4,F6,F7.claims,F7.census", out = ""))
vs_threads(1L)
options(width = 180L)
only <- trimws(strsplit(args$only, ",", fixed = TRUE)[[1]])
cat(sprintf("wp7_sim: choicer %s from %s; sections %s\n", b$build$version, b$build$path,
            paste(only, collapse = ", ")))

verdicts <- logical(0)
report <- function(label, value, tol, pass = isTRUE(value < tol)) {
  cat(sprintf("%s  %-78s %11.4g  (threshold %.3g)\n", if (pass) "PASS" else "FAIL", label,
              value, tol))
  verdicts[length(verdicts) + 1L] <<- pass
}
now <- function() proc.time()[["elapsed"]]
save_out <- function(obj, name) {
  if (nzchar(args$out)) {
    dir.create(args$out, showWarnings = FALSE, recursive = TRUE)
    saveRDS(obj, file.path(args$out, paste0(name, ".rds")))
  }
}
split_arg <- function(x) trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
valid_record <- function(rec) {
  tryCatch({
    do.call(new_rep_record, rec[c("family", "cell", "method", "rep", "seed", "ok", "error",
                                  "warnings", "est", "vcov_c", "joint_wald", "lr", "derived",
                                  "diag", "meta")])
    TRUE
  }, error = function(e) FALSE)
}
f10 <- load_family("F10", b$L)
f10_env <- environment(f10$simulate)

# The old script's functions, parsed from its text (no top-level code runs).
# F10 replaced the script and it left the tree; it lives in git history at
# commit ae38e12 (on main), read here with `git show`, so sections shares and
# repro need a clone with that commit.
OLD_SCRIPT <- "ae38e12:_validation/blp_concentration_mc.R"
old_script <- function() {
  src <- suppressWarnings(system2("git", c("-C", shQuote(root), "show", OLD_SCRIPT), stdout = TRUE,
                                  stderr = FALSE))
  if (!is.null(attr(src, "status"))) {
    stop("`git show ", OLD_SCRIPT, "` failed: sections shares and repro need a clone of the ",
         "repository with that commit", call. = FALSE)
  }
  ex <- parse(text = src, keep.source = FALSE)
  env <- new.env(parent = globalenv())
  for (e in ex) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
        identical(e[[3]][[1]], as.name("function"))) eval(e, env)
  }
  env
}

# --- load ------------------------------------------------------------------------

if ("load" %in% only) {
  cat("\n== load: F5 and F10 validate ==\n")
  f5 <- tryCatch(load_family("F5", b$L), error = function(e) conditionMessage(e))
  report("F5 loads and validates", as.numeric(!is.list(f5)), 0.5)
  report("F10 loads and validates", as.numeric(!is.list(f10)), 0.5)
  if (is.list(f5)) {
    print(f5$cells, row.names = FALSE)
    for (k in seq_len(nrow(f5$cells))) {
      cl <- f5$cells[k, ]
      src <- load_family(cl$src, b$L)
      row <- tryCatch(family_cell(src, cl$src_cell), error = function(e) NULL)
      report(sprintf("F5 regime %s: %s has cell %s with N = %d", cl$cell, cl$src, cl$src_cell,
                     cl$N), as.numeric(is.null(row) || row$N != cl$N), 0.5)
    }
  }
  print(f10$cells, row.names = FALSE)
}

# --- rates -----------------------------------------------------------------------

# F5's b_se and slope_se on synthetic records, no fits. A record holds z (S by
# seed), with store and source rows where given, as named vectors (S -> z). The
# reference's SE is 1, so z is the estimate.
if ("rates" %in% only) {
  cat("\n== rates: F5's rate uncertainty, from complete units ==\n")
  f5 <- load_family("F5", b$L)
  S <- c(50L, 100L, 200L, 400L, 800L)
  S_ref <- 3200L
  x <- 1 / S - 1 / S_ref
  rec5 <- function(d, z, store = NULL, src = NULL) {
    nm <- c(as.vector(outer(paste0("price/S", S), paste0("g", seq_len(ncol(z))), paste, sep = "/")),
            sprintf("price/S%s/store", names(store)), sprintf("price/S%s/src", names(src)))
    e <- data.frame(name = c("price/ref", nm), block = "beta", est = c(0, z, store, src), truth = 0,
                    se_hessian = NA_real_, se_bhhh = NA_real_, se_robust = NA_real_,
                    se_cluster = NA_real_, gated = FALSE, stringsAsFactors = FALSE)
    new_rep_record("F5", "F4", "msl", d, seed = d, est = e,
                   diag = list(se_ref = c(price = 1), S_ref = S_ref, units = 3200L,
                               src = c(S_src = 500L)))
  }
  # The diagonal formula's b_se from the S-level MCSEs of one ladder: right only
  # when the S points are independent.
  diag_se <- function(by_S) {
    d <- by_S[by_S$kind == "generate", ]
    1 / sqrt(sum((1 / d$S - 1 / S_ref)^2 / d$bias_se_mcse^2))
  }

  # z_d(S) = B_d (1/S - 1/3200), the same at both seeds: b_d = B_d, so b = mean(B)
  # and its MCSE is sd(B) / sqrt(30); B_d scales every S alike, so the slope does
  # not depend on it. The diagonal formula reports 2.033325, too small by sqrt(5).
  set.seed(104)
  B <- stats::rnorm(30, 100, 30)
  lad <- f5$findings(lapply(1:30, function(d) rec5(d, cbind(B[d] * x, B[d] * x))))$rates
  se_B <- stats::sd(B) / sqrt(30)
  report("correlated ladder: sd(B) / sqrt(30) is the review's 4.546653", abs(se_B - 4.546653),
         1e-6)
  report("correlated ladder: one rates row, b = mean(B)", abs(lad$b - mean(B)), 1e-9,
         pass = nrow(lad) == 1L && abs(lad$b - mean(B)) < 1e-9)
  report("correlated ladder: b_se = sd(B) / sqrt(30)", abs(lad$b_se - se_B), 1e-9)
  report("correlated ladder: slope_se is 0", lad$slope_se, 1e-10)

  # Independent noise at each S and seed, no part shared across S (200 datasets,
  # so that the ratio is not itself noisy): the diagonal formula is right here, and
  # the units must agree with it.
  set.seed(105)
  ind <- f5$findings(lapply(1:200, function(d) {
    rec5(d, 100 * cbind(x, x) + sqrt(30 / S) * matrix(stats::rnorm(10), 5L, 2L))
  }))
  cat(sprintf("      independent S points: b_se %.4f, diagonal %.4f, slope_se %.4f\n",
              ind$rates$b_se, diag_se(ind$by_S), ind$rates$slope_se))
  report("independent S points: b_se within 15% of the diagonal formula",
         abs(ind$rates$b_se / diag_se(ind$by_S) - 1), 0.15)

  # One dataset (the F7 regimes): the unit is the seed, and store mode, one fit per
  # S, has no units to spread over.
  set.seed(106)
  z1 <- 100 * cbind(x, x) + matrix(stats::rnorm(10, 0, 0.1), 5L, 2L)
  one <- tryCatch(f5$findings(list(rec5(1L, z1, store = c("100" = 0.9, "200" = 0.4)))),
                  error = function(e) conditionMessage(e))
  report("one dataset, two seeds, store rows: findings() returns its tables",
         as.numeric(!is.list(one)), 0.5)
  if (is.list(one)) {
    g <- one$rates[one$rates$kind == "generate", ]
    s <- one$rates[one$rates$kind == "store", ]
    d <- one$by_S[one$by_S$kind == "generate", ]
    cb <- (x / d$bias_se_mcse^2) / sum(x^2 / d$bias_se_mcse^2)
    report("one dataset: b_se is the two seeds' spread / sqrt(2)",
           abs(g$b_se - stats::sd(drop(cb %*% z1)) / sqrt(2)), 1e-9)
    report("one dataset: store's b_se and slope_se are NA, its b is not",
           as.numeric(!(nrow(s) == 1L && is.na(s$b_se) && is.na(s$slope_se) && is.finite(s$b))),
           0.5)
  }

  # Source rows at an S on the ladder: step 1's fit, off the ladder's seeds. They
  # must reach by_S and stay out of the rates.
  recs <- lapply(1:30, function(d) rec5(d, cbind(B[d] * x, B[d] * x)))
  src <- f5$findings(lapply(1:30, function(d) {
    rec5(d, cbind(B[d] * x, B[d] * x), src = c("200" = 50 + d))
  }))
  report("source rows at S = 200 reach by_S", as.numeric(!any(src$by_S$kind == "source" &
                                                                src$by_S$S == 200L)), 0.5)
  report("source rows at S = 200 leave the rates as they were",
         as.numeric(!identical(src$rates, f5$findings(recs)$rates)), 0.5)

  # Calibration: B_d a shared part per dataset, and each seed's noise correlated
  # 0.9 across S (nested draws). The mean reported SE over the SD of b; the
  # diagonal formula's ratio beside it shows the design tells the two apart (0.50).
  set.seed(107)
  rho <- 0.9
  t0 <- now()
  cal <- vapply(seq_len(300L), function(i) {
    B <- stats::rnorm(30, 100, 30)
    f <- f5$findings(lapply(1:30, function(d) {
      seeds <- function() B[d] * x + sqrt(30 / S) * (rho * stats::rnorm(1) +
                                                       sqrt(1 - rho^2) * stats::rnorm(5))
      rec5(d, cbind(seeds(), seeds()))
    }))
    c(b = f$rates$b, b_se = f$rates$b_se, diag = diag_se(f$by_S))
  }, numeric(3))
  sd_b <- stats::sd(cal["b", ])
  cat(sprintf(paste("      300 ensembles (%.0f s): SD of b %.3f; mean b_se %.3f, ratio %.3f;",
                    "diagonal %.3f, ratio %.3f\n"),
              now() - t0, sd_b, mean(cal["b_se", ]), mean(cal["b_se", ]) / sd_b,
              mean(cal["diag", ]), mean(cal["diag", ]) / sd_b))
  report("calibration: mean b_se / SD(b) within 0.85-1.15", abs(mean(cal["b_se", ]) / sd_b - 1),
         0.15)
  report("calibration: the design tells them apart, diagonal formula / SD(b) below 0.7",
         mean(cal["diag", ]) / sd_b, 0.7)
}

# --- shares ----------------------------------------------------------------------

if ("shares" %in% only) {
  cat("\n== shares: F10's frozen exact targets ==\n")
  old <- old_script()
  for (d in split_arg(args$designs)) {
    t0 <- now()
    po <- old$pop_shares(old$design_args(d))
    pn <- f10_env$f10_pop_shares(family_cell(f10, d))
    cl <- family_cell(f10, d)
    report(sprintf("%s: old pop_shares() = F10_SHARES, bit for bit (%.0f s)", d, now() - t0),
           max(abs(po$within - f10_env$F10_SHARES[[d]])), 1e-300,
           pass = identical(po$within, f10_env$F10_SHARES[[d]]))
    report(sprintf("%s: f10_pop_shares() = old pop_shares(), markets too", d),
           max(abs(pn$within - po$within)), 1e-300,
           pass = identical(pn, po) &&
             identical(as.integer(po$market), rep(seq_len(cl$markets), each = 10L / cl$markets)))
  }
}

# --- repro -----------------------------------------------------------------------

# The old record's estimator `e` as est / se / se_r vectors, and the port's.
if ("repro" %in% only) {
  cat("\n== repro: the port against blp_concentration_mc.R, rep by rep ==\n")
  old <- old_script()
  map <- c(ml = "ml", blp_sample = "blp_sample", blp_exact = "blp_exact", blp_n1 = "blp_n1",
           blp_n10 = "blp_n10")
  rows <- list()
  recs_all <- list()
  for (tk in split_arg(args$repro)) {
    p <- strsplit(tk, ":", fixed = TRUE)[[1]]
    d <- p[1]
    ab <- as.integer(strsplit(p[2], "-", fixed = TRUE)[[1]])
    reps <- seq.int(ab[1], ab[length(ab)])
    cl <- family_cell(f10, d)
    pop <- list(within = f10_env$F10_SHARES[[d]],
                market = rep(seq_len(cl$markets), each = 10L / cl$markets))
    a <- old$design_args(d)
    for (r in reps) {
      t0 <- now()
      o <- old$one_rep(r, d, a, pop)
      t_old <- now() - t0
      t0 <- now()
      recs <- run_rep(f10, d, r, names(map))
      t_new <- now() - t0
      recs_all <- c(recs_all, recs)
      report(sprintf("%s rep %d: the port's data are the old seed's (old seed %d)", d, r, o$seed),
             as.numeric(!identical(recs[[1]]$diag$old_seed, as.integer(o$seed))), 0.5)
      for (m in names(map)) {
        rec <- recs[[m]]
        oe <- o[[map[[m]]]]
        if (is.null(oe) || !is.null(oe$error) || !isTRUE(rec$ok)) {
          report(sprintf("%s rep %d %s: both fitted", d, r, m), 1, 0.5)
          cat("      old:", oe$error %||% "ok", "| port:", rec$error %||% "ok", "\n")
          next
        }
        e <- rec$est
        st <- rec$diag$start
        dd <- c(est = max(abs(e$est - oe$est)), se = max(abs(e$se_hessian - oe$se)),
                se_r = max(abs(e$se_robust - oe$se_r)))
        shift <- max(abs(e$est - oe$est) / e$se_hessian)
        rows[[length(rows) + 1L]] <- data.frame(
          design = d, rep = r, method = m, winner = st$winner,
          dlogL_default_minus_theta0 = signif(st$dlogL, 3), max_dest = signif(dd[["est"]], 3),
          max_dse = signif(dd[["se"]], 3), max_dse_robust = signif(dd[["se_r"]], 3),
          shift_se = signif(shift, 3), code = rec$diag$code,
          old_code = if (m == "ml") NA else oe$convergence, stringsAsFactors = FALSE)
        if (identical(st$winner, "default")) {
          report(sprintf("%s rep %d %s: default start won; est and SEs = the old script's", d, r,
                         m), max(dd), 1e-12)
        } else {
          cat(sprintf(paste("INFO  %s rep %d %s: the theta0 start won (loglik gain %.3g); the",
                            "rep moved by up to %.3g SE from the old script's\n"), d, r, m,
                      -st$dlogL, shift))
        }
      }
      cat(sprintf("      %s rep %d: old script %.1f s, port %.1f s (two starts)\n", d, r, t_old,
                  t_new))
    }
  }
  tab <- do.call(rbind, rows)
  print(tab, row.names = FALSE)
  save_out(list(table = tab, records = recs_all), "repro")
}

# --- smoke -----------------------------------------------------------------------

if ("smoke" %in% only) {
  cat("\n== smoke: F10 through run_rep(), every cell and method ==\n")
  ab <- as.integer(strsplit(args$smoke_reps, ":", fixed = TRUE)[[1]])
  reps <- seq.int(ab[1], ab[length(ab)])
  tab <- list()
  recs_all <- list()
  for (cl in f10$cells$cell) {
    for (r in c(reps, if (cl == "cross") 10L)) {
      recs <- run_rep(f10, cl, r, cell_methods(f10, cl))
      recs_all <- c(recs_all, recs)
      for (rec in recs) {
        valid <- valid_record(rec)
        w <- unique(rec$warnings$class)
        tab[[length(tab) + 1L]] <- data.frame(
          family = rec$family, cell = rec$cell, method = rec$method, rep = rec$rep, ok = rec$ok,
          valid = valid, sim_sec = round(rec$diag$sim_sec, 2),
          unit_sec = round(rec$diag$unit_sec, 2),
          fit_sec = round(rec$diag$time_sec %||% NA_real_, 2),
          iterations = rec$diag$iterations %||% NA, code = rec$diag$code %||% NA,
          decrement = signif(c(rec$diag$decrement, rec$diag$concentration[["decrement"]])[1] %||%
                               NA_real_, 2),
          evals = rec$diag$concentration[["evaluations"]] %||% NA,
          start = rec$diag$start$winner %||% NA_character_,
          start_dll = signif(rec$diag$start$dlogL %||% NA_real_, 3),
          ms_dlogL = signif(rec$diag$multistart$dlogL %||% NA_real_, 2),
          joint_p = signif(rec$joint_wald$p %||% NA_real_, 2),
          warnings = if (length(w)) paste(w, collapse = " ") else "",
          error = substr(rec$error %||% "", 1, 60), stringsAsFactors = FALSE)
        report(sprintf("F10/%s/%s rep %d: ok and valid%s", rec$cell, rec$method, rec$rep,
                       if (rec$rep %% 10 == 0) ", with the multi-start check" else ""),
               as.numeric(!(rec$ok && valid && (rec$rep %% 10 != 0 ||
                                                  isTRUE(is.finite(rec$diag$multistart$dlogL))))),
               0.5)
      }
    }
  }
  tab <- do.call(rbind, tab)
  print(tab, row.names = FALSE)
  units <- unique(tab[c("cell", "method")])
  for (k in seq_len(nrow(units))) {
    rs <- Filter(function(r) r$cell == units$cell[k] && r$method == units$method[k], recs_all)
    mt <- tryCatch(gate_family(metrics_cell(rs, f10$gates)),
                   error = function(e) conditionMessage(e))
    report(sprintf("metrics_cell() on F10/%s/%s: %s rows, none gated", units$cell[k],
                   units$method[k], if (is.data.frame(mt)) nrow(mt) else "error"),
           as.numeric(!is.data.frame(mt) || any(mt$gated)), 0.5)
  }
  fnd <- tryCatch(f10$findings(recs_all), error = function(e) conditionMessage(e))
  report("f10_findings() on the smoke records gives the two tables",
         as.numeric(!is.list(fnd) || !all(c("by_design", "article") %in% names(fnd))), 0.5)
  if (is.list(fnd)) {
    num <- vapply(fnd$by_design, is.numeric, TRUE)
    bd <- fnd$by_design
    bd[num] <- lapply(bd[num], signif, 3)
    print(bd, row.names = FALSE)
    print(fnd$article, row.names = FALSE)
  }
  # The budget: the median seconds per method over the reps without the
  # multi-start, plus its cost on one rep in ten (rep 10's excess).
  plain <- tab$rep %% 10 != 0
  med <- tapply(tab$unit_sec[plain], tab$method[plain], stats::median)
  ms <- tab[!plain, , drop = FALSE]
  extra <- if (nrow(ms)) sum(pmax(ms$unit_sec - med[ms$method], 0)) / 10 else 0
  cat(sprintf("\nPer rep (median over cells and reps, seconds, both starts): %s\n",
              paste(sprintf("%s %.1f", names(med), med), collapse = ", ")))
  per_rep <- sum(med) + extra
  cat(sprintf(paste("mc budget: %.1f s per rep (5 methods, the multi-start's %.1f s included) x 4",
                    "designs x 1000 reps = %.1f CPU-h (plan ~22)\n"), per_rep, extra,
              per_rep * 4000 / 3600))
  save_out(list(table = tab, records = recs_all, findings = fnd), "smoke")
}

# --- spec ------------------------------------------------------------------------

# Each regime's spec cell and method in its source family (a small cell of
# the same model; F7's at N / --spec_div).
SPEC <- list(F3 = c("1k.q", "ml"), F4 = c("200.q", "panel"), F6 = c("5k", "mxl"),
             F7.claims = c("claims", "mxl"), F7.census = c("census", "mxl"))

if ("spec" %in% only) {
  cat("\n== spec: F5's models against the source families' fits ==\n")
  f5 <- load_family("F5", b$L)
  for (rg in intersect(split_arg(args$spec), f5$cells$cell)) {
    src_id <- f5$cells$src[f5$cells$cell == rg]
    src <- load_family(src_id, b$L)
    cl <- SPEC[[rg]][1]
    m <- SPEC[[rg]][2]
    if (src_id == "F7") {
      k <- src$cells$cell == cl
      src$cells$N[k] <- as.integer(src$cells$N[k] / as.numeric(args$spec_div))
    }
    row <- family_cell(src, cl)
    t0 <- now()
    rec <- run_unit(src, cl, m, 1L)
    if (!isTRUE(rec$ok)) {
      report(sprintf("%s: %s/%s rep 1 fitted", rg, src_id, cl), 1, 0.5)
      cat("      error:", rec$error, "\n")
      next
    }
    seed <- vs_seed(src_id, cl, 1L)
    vs_set_seed(seed)
    sim <- src$simulate(row, seed)
    model <- function(...) {
      suppressMessages(suppressWarnings(environment(f5$simulate)$F5_MODELS[[rg]](
        sim, S = row$S, draws = "generate", seed = seed + 1L, ...)))
    }
    # F7's fit of rep 1 is Newton steps from theta0 (f7_estimate()); the
    # others' the two starts' winner.
    if (src_id == "F7") {
      winner <- if (isTRUE(rec$diag$newton_only)) "Newton, theta0" else "fallback, theta0"
      fit <- fit_newton(model, sim$truth)
    } else {
      winner <- rec$diag$start$winner %||% "default"
      fit <- model(theta_init = if (winner == "theta0") sim$truth)
    }
    d <- max(abs(coef(fit)[rec$est$name] - rec$est$est))
    report(sprintf("%s: F5's model = %s/%s/%s's fit, rep 1, N = %d, S = %d (%s start; %.0f s)",
                   rg, src_id, cl, m, row$N, row$S, winner, now() - t0), d, 1e-12)
  }
}

# --- f5 --------------------------------------------------------------------------

if ("f5" %in% only) {
  cat("\n== f5: one dataset per regime at a reduced ladder ==\n")
  f5 <- load_family("F5", b$L)
  full <- f5$cells
  small <- full
  small$S_ref <- as.integer(args$f5_ref)
  small$S_gen <- as.character(args$f5_gen)
  small$S_store <- as.character(args$f5_store)
  small$seeds <- as.integer(args$f5_seeds)
  f5s <- f5
  f5s$cells <- small
  # F7's regimes at N / --f5_div: F7 with those cells divided, put in F5's
  # cache of loaded sources (f5s shares f5's environment).
  div <- stats::setNames(ifelse(small$src == "F7", as.numeric(args$f5_div), 1), small$cell)
  if (any(small$src == "F7")) {
    f7 <- load_family("F7", b$L)
    for (k in which(small$src == "F7")) {
      i <- f7$cells$cell == small$src_cell[k]
      f7$cells$N[i] <- as.integer(f7$cells$N[i] / div[[small$cell[k]]])
    }
    assign("F7", f7, envir = environment(f5$simulate)$.f5_loaded)
  }
  recs <- list()
  for (cl in intersect(split_arg(args$f5), small$cell)) {
    t0 <- now()
    rec <- run_unit(f5s, cl, "msl", 1L)
    recs[[cl]] <- rec
    report(sprintf("F5/%s rep 1 (N / %g, S_ref %d, ladder %s, store %s): ok and valid (%.0f s)",
                   cl, div[[cl]], small$S_ref[1], small$S_gen[1], small$S_store[1], now() - t0),
           as.numeric(!(rec$ok && valid_record(rec))), 0.5)
    if (!rec$ok) {
      cat("      error:", rec$error, "\n")
      next
    }
    cat(sprintf("      warnings: %s\n", paste(unique(rec$warnings$class), collapse = " ")))
    fits <- rec$diag$fits
    two <- !is.null(rec$diag$start)
    if (two) {
      cat(sprintf("      the source's fit (two starts): winner %s, dlogL (default - theta0) %.3g\n",
                  rec$diag$start$winner, rec$diag$start$dlogL))
    } else {
      src <- fits[fits$variant == "src", ]
      cat(sprintf("      the source's fit (Newton steps from theta0): newton_only %s, %d steps\n",
                  src$newton_only, src$newton_steps))
    }
    fits$sec_per_draw <- fits$time_sec / fits$S
    print(fits, row.names = FALSE)
    # Extrapolation to the regime's full ladder, linear in draws and, at N /
    # div, in rows (times div): step 1 (the source's S, its route) as
    # measured (the unit's time less the other fits'); the warm reference at
    # its seconds per draw times the full S_ref; the warm ladder at the
    # ladder's median seconds per draw times each S.
    fr <- full[full$cell == cl, ]
    is_src <- fits$variant == "src"
    is_ref <- fits$variant == "ref"
    t_src <- (rec$diag$unit_sec - sum(fits$time_sec[!is_src])) * div[[cl]]
    ref_spd <- fits$sec_per_draw[is_ref] * div[[cl]]
    lad_spd <- stats::median(fits$sec_per_draw[!is_src & !is_ref]) * div[[cl]]
    S_lad <- c(rep(environment(f5$simulate)$.f5_S(fr$S_gen), fr$seeds),
               environment(f5$simulate)$.f5_S(fr$S_store))
    est_ref <- ref_spd * fr$S_ref
    est_lad <- lad_spd * sum(S_lad)
    tot <- t_src + est_ref + est_lad
    cat(sprintf(paste("      full ladder (extrapolated, linear in draws%s): the source's fit",
                      "%.0f s (S = %d, %s), reference %.0f s (S = %d, %.4f s per draw),",
                      "ladder %.0f s (%d fits, %.4f s per draw); %.2f h per dataset; x %d",
                      "datasets = %.1f CPU-h\n"),
                if (div[[cl]] > 1) sprintf(" and in rows, x %g", div[[cl]]) else "",
                t_src, fits$S[is_src], if (two) "two starts" else "Newton from theta0",
                est_ref, fr$S_ref, ref_spd, est_lad, length(S_lad),
                lad_spd, tot / 3600, f5$R$mc[[cl]], tot * f5$R$mc[[cl]] / 3600))
  }
  fnd <- tryCatch(f5$findings(recs), error = function(e) conditionMessage(e))
  report("f5_findings() on the records gives by_S, rates and regimes",
         as.numeric(!is.list(fnd) || !all(c("by_S", "rates", "regimes") %in% names(fnd))), 0.5)
  if (is.list(fnd)) {
    for (k in names(fnd)) {
      x <- fnd[[k]]
      num <- vapply(x, is.numeric, TRUE)
      x[num] <- lapply(x[num], signif, 3)
      cat(sprintf("\n%s:\n", k))
      print(x, row.names = FALSE)
    }
  }
  mt <- lapply(recs, function(r) tryCatch(gate_family(metrics_cell(list(r), f5$gates)),
                                         error = function(e) conditionMessage(e)))
  report("metrics_cell() on each F5 record, none gated",
         as.numeric(!all(vapply(mt, function(m) is.data.frame(m) && !any(m$gated), TRUE))), 0.5)
  save_out(list(records = recs, findings = fnd), "f5")
}

cat(sprintf("\nwp7_sim: %d PASS, %d FAIL\n", sum(verdicts), sum(!verdicts)))
if (any(!verdicts)) quit(status = 1L)
