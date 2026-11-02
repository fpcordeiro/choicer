# ============================================================================
# WP6's checks: the geographic generator's C++ port (lib/geo_gen.cpp) against
# its R reference (lib/dgp_geo.R), the designs D3, D3-claims and D4-census,
# and families F6 and F7. Always first: the port's build (geo_gen_compile():
# the distances' bitwise self-test, and a held cache lock reported). Sections
# (--only, comma-separated):
#   q0     D3's probabilities at theta0 against choicer: MNL predict() (max
#          abs diff < 1e-12) on D3, D3-claims and D4-census; MXL the suite's
#          simulated log-likelihood on a shared cube and on the DGP's own
#          tastes (a one-draw cube) against choicer's fixed-parameter fit (<
#          1e-10 relative), cross-section (D3) and panel (D3-claims,
#          person_col); the draws' order and reproducibility.
#   port   the port against the reference: every design and variant, three
#          seeds: identical data (sets, distances, covariates, choices),
#          probabilities within 1e-14, expected counts within 1e-12
#          relative; the radius search on edge cases (K_min binding for
#          everyone, K_max binding, fewer hospitals than K_min, K_min = K_max,
#          lattice points with ties at the radius and cell edges rounded).
#   calib  the moments at N = 200,000 (seed 1) against the targets (D3 and
#          D3-claims per variant, D4-census per variant), and per rung of F6
#          and at F7's sizes: expected counts per alternative and the
#          expected number never chosen.
#   full   one generation at full size with the port (--full=claims or
#          census): seconds, rows, R's heap; run it under /usr/bin/time -l
#          for the process's peak RSS.
#   smoke  F6 and F7 through run_rep() (--smoke=F6,F7, or F6:20k:blp, ...):
#          R = 1 per cell x method, plus reps 2 and 10 at F6's bottom rungs
#          (the multi-start); F7's MXL cells ("claims", "census") at N /
#          --f7_div and S = --f7_S, its MNL cells at full size unless
#          --f7_mnl_div > 1; --S=<cell>:<S>,... overrides a cell's S (timing at
#          a smaller S); each record re-validated.
# Usage (the dev build first on R_LIBS; one thread; each section < 10 min):
#   R_LIBS=<lib> OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 VS_GEO_CACHE=<dir> \
#     Rscript _validation/checks/wp6_geo.R --only=q0,port,calib
# VS_GEO_CACHE is the port's build cache (geo_gen_compile()); default
# _validation/output/geo_gen_cache. Exits with status 1 on any FAIL.
# ============================================================================

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(normalizePath(sub("^--file=", "", file_arg))), ".."))
sys.source(file.path(root, "lib", "common.R"), envir = (boot <- new.env()), keep.source = FALSE)
b <- boot$vs_boot(root)
suppressPackageStartupMessages(library(data.table))
args <- vs_args(defaults = list(only = "q0,port,calib", smoke = "F6,F7", full = "claims",
                                f7_div = 25L, f7_S = 100L, f7_mnl_div = 1L, S = NULL))
vs_threads(1L)
options(width = 160L)
only <- trimws(strsplit(args$only, ",", fixed = TRUE)[[1]])
cat(sprintf("wp6_geo: choicer %s from %s; sections %s\n", b$build$version, b$build$path,
            paste(only, collapse = ", ")))
t_start <- proc.time()[["elapsed"]]

verdicts <- logical(0)
report <- function(label, value, tol, pass = isTRUE(value < tol)) {
  cat(sprintf("%s  %-74s %11.4g  (threshold %.3g)\n", if (pass) "PASS" else "FAIL", label,
              value, tol))
  verdicts[length(verdicts) + 1L] <<- pass
}

# The port's build: compiled or loaded under the cache's lock, its
# distances checked against R's bit for bit; a held lock makes
# geo_gen_compile() wait, then stop naming it.
t0 <- proc.time()[["elapsed"]]
cache <- geo_gen_compile()
cat(sprintf("geo_gen_compile(): %s, %.1f s\n", cache, proc.time()[["elapsed"]] - t0))
report("the port's distances equal R's bit for bit (.geo_selftest())",
       as.numeric(!isTRUE(tryCatch(.geo_selftest(cache), error = function(e) FALSE))), 0.5)
lock <- .geo_lock(cache, 60)
err <- tryCatch({ geo_gen_compile(wait = 1); "" }, error = function(e) conditionMessage(e))
unlink(lock, recursive = TRUE)
report("a held lock: geo_gen_compile(wait = 1) stops naming it",
       as.numeric(!grepl("locked by another process", err, fixed = TRUE)), 0.5)
quiet <- function(expr) {
  withCallingHandlers(suppressMessages(expr),
                      choicer_unbounded_constants = function(w) invokeRestart("muffleWarning"))
}
now <- function() proc.time()[["elapsed"]]

# --- q0 ------------------------------------------------------------------------------

if ("q0" %in% only) {
  cat("\n== q0: D3's probabilities at theta0 against choicer ==\n")
  dgp <- get_dgp("D3-claims")
  s1 <- dgp$simulate(N = 400L, seed = 11L, variant = "default")
  s2 <- dgp$simulate(N = 400L, seed = 11L, variant = "default")
  report("D3-claims simulate() twice at one seed: identical data and truth",
         as.numeric(!(identical(s1$data, s2$data) && identical(s1$truth, s2$truth))), 0.5)
  dr <- geo_draw(dgp$design, 400L, 11L)
  A <- sum(dr$T)
  stream <- vs_with_seed(11L, stats::runif(3L * 400L + 400L + 400L + A + 2L * 400L + A))
  report("D3 panel: the choice uniforms are the sample stream's last A",
         as.numeric(!identical(utils::tail(stream, A), dr$sample$choice)), 0.5)
  report("D3 panel: the admissions are min(T_max, floor(u^(-1/alpha))) of the third block",
         as.numeric(!identical(.d3_admissions(stream[1600L + seq_len(400L)], dgp$design),
                               dr$T)), 0.5)
  report("D3: the world does not depend on N or the seed",
         as.numeric(!identical(geo_draw(dgp$design, 10L, 1L)$world, dr$world)), 0.5)
  # MNL: predict() at theta0 on the alternatives present.
  for (id in c("D3", "D3-claims", "D4-census")) for (sd in 1:2) {
    N <- if (id == "D3") 3000L else 2500L
    sim <- get_dgp(id)$simulate(N = N, seed = sd, variant = "mnl", keep_probs = TRUE,
                                engine = "cpp")
    th <- geo_truth(sim$meta, unique(sim$data$alt))
    x <- if (id == "D4-census") D4_X else D3_X$mnl
    fit <- quiet(run_mnlogit(sim$data, "id", "alt", "choice", x,
                             include_outside_option = sim$meta$spec$outside))
    ok_names <- identical(names(coef(fit)), names(th))
    fit$coefficients[names(th)] <- th
    pr <- predict(fit, type = "probabilities")
    diff <- max(abs(pr$choice_prob[, 1] - sim$meta$probs$p),
                if (sim$meta$spec$outside) abs(pr$choice_prob_outside[, 1] - sim$meta$probs$p0))
    report(sprintf("%s MNL seed %d: names = choicer's", id, sd), as.numeric(!ok_names), 0.5)
    report(sprintf("%s MNL seed %d (%d rows): max|p - predict()| at theta0", id, sd,
                   nrow(sim$data)), diff, 1e-12)
  }
  # MXL: the simulated log-likelihood on a shared cube, and the DGP's tastes.
  for (id in c("D3", "D3-claims")) for (sd in 1:2) {
    N <- if (id == "D3") 1500L else 800L
    sim <- get_dgp(id)$simulate(N = N, seed = sd, variant = "default", keep_probs = TRUE,
                                engine = "cpp")
    panel <- isTRUE(sim$meta$design$panel)
    th <- geo_truth(sim$meta, unique(sim$data$alt))
    spec <- sim$meta$spec
    spec$alts <- sort(unique(sim$data$alt))
    d <- design(sim$data, spec, person = if (panel) "patient")
    chosen <- chosen_rows(sim$data$choice, d$sit)
    idat <- prepare_mxl_data(sim$data, "id", "alt", "choice", D3_X$default, D3_W,
                             person_col = if (panel) "patient", rc_correlation = FALSE)
    fixed_ll <- function(eta) {
      quiet(run_mxlogit(input_data = idat, eta_draws = eta, rc_dist = c(1L, 0L), rc_mean = TRUE,
                        rc_correlation = FALSE, theta_init = th, lower = th, upper = th))$loglik
    }
    U <- if (panel) sim$meta$N else sim$meta$A
    eta <- vs_with_seed(100L + sd, array(stats::rnorm(2L * 40L * U), c(2L, 40L, U)))
    ll <- mxl_loglik(d, unpack_theta("mxl", spec, th), eta, chosen)
    cl <- fixed_ll(eta)
    report(sprintf("%s MXL seed %d (%s): simulated ll on a shared cube (S = 40), rel", id, sd,
                   if (panel) "panel" else "cross-section"), abs(ll - cl) / abs(ll), 1e-10)
    eta1 <- array(sim$meta$eta, c(2L, 1L, sim$meta$N))
    p <- sim$meta$probs$p
    ll1 <- sum(log(p[chosen]))
    report(sprintf("%s MXL seed %d: the DGP's own tastes as a one-draw cube, rel", id, sd),
           abs(ll1 - fixed_ll(eta1)) / abs(ll1), 1e-10)
  }
}

# --- port ----------------------------------------------------------------------------

if ("port" %in% only) {
  cat("\n== port: lib/geo_gen.cpp against the R reference ==\n")
  for (id in c("D3", "D3-claims", "D4", "D4-40", "D4-census")) {
    dgp <- get_dgp(id)
    N <- c(D3 = 3000L, `D3-claims` = 2000L, D4 = 2000L, `D4-40` = 2000L, `D4-census` = 3000L)[[id]]
    for (v in dgp$variants) for (sd in 1:3) {
      t0 <- now()
      r <- dgp$simulate(N = N, seed = sd, variant = v, keep_probs = TRUE)
      t1 <- now()
      cp <- dgp$simulate(N = N, seed = sd, variant = v, keep_probs = TRUE, engine = "cpp")
      t2 <- now()
      same <- identical(r$data, cp$data) && identical(r$truth, cp$truth)
      dp <- max(abs(r$meta$probs$p - cp$meta$probs$p), abs(r$meta$probs$p0 - cp$meta$probs$p0))
      de <- max(abs(r$meta$expected - cp$meta$expected) / pmax(abs(r$meta$expected), 1e-300))
      report(sprintf("%s %s seed %d (%d rows; R %.1fs, C++ %.2fs): identical data", id, v, sd,
                     nrow(r$data), t1 - t0, t2 - t1), as.numeric(!same), 0.5)
      report(sprintf("%s %s seed %d: max|p_R - p_C++|", id, v, sd), dp, 1e-14)
      report(sprintf("%s %s seed %d: expected counts, max rel diff", id, v, sd), de, 1e-12)
    }
  }
  # The radius search on edge cases, against the brute force.
  base <- GEO_D3
  cases <- list(
    "radius 0.01 (K_min binding everywhere)" = list(radius = 0.01),
    "radius 40 (K_max binding everywhere)" = list(radius = 40),
    "two hospitals, K_min = 3" = list(hospitals = 2L),
    "K_min = K_max = 7" = list(K_min = 7L, K_max = 7L),
    "one town, radius 1" = list(clusters = 1L, radius = 1))
  # Lattices 1.1 apart, hospitals at i * 11 / 10 and patients at seq(0, by =
  # 1.1) or i * 11 / 10, radius 2.2 or 4.4: floor(v / cell) rounds points on
  # a cell edge into the next cell and many hospitals tie at the radius (the
  # search without its ring of slack missed 42 and 84 rows here).
  lattice <- function(r, pat) {
    i <- 0:20
    hx <- i * 11 / 10
    px <- if (pat == "seq") seq(0, by = 1.1, length.out = 21L) else hx
    list(h = data.table(j = seq_len(441L), x = rep(hx, 21L), y = rep(hx, each = 21L)),
         p = data.table(patient = seq_len(441L), x = rep(px, 21L), y = rep(px, each = 21L)),
         des = utils::modifyList(base, list(side = 22 + 1e-9, radius = r, K_min = 1L,
                                            K_max = 10000L)))
  }
  cases[["lattice, radius 2.2, ties at the radius and on cell edges"]] <- lattice(2.2, "seq")
  cases[["lattice, radius 4.4, ties at the radius and on cell edges"]] <- lattice(4.4, "div")
  for (k in names(cases)) {
    if (startsWith(k, "lattice")) {
      des <- cases[[k]]$des
      h <- cases[[k]]$h
      p <- cases[[k]]$p
    } else {
      des <- utils::modifyList(base, cases[[k]])
      dr <- geo_draw(des, 1500L, 3L)
      h <- geo_hospitals(des, dr$world)
      p <- d3_patients(des, dr)
    }
    ref <- geo_sets_radius(p, h, des)
    cpp <- setDT(.geo_gen()$geo_radius(p$x, p$y, rep(1L, nrow(p)), h$x, h$y, des$side, des$radius,
                                        des$K_min, des$K_max))
    same <- identical(ref$patient, cpp$patient) && identical(ref$j, cpp$j) &&
      identical(ref$dist, cpp$dist)
    report(sprintf("radius search, %s (%d rows): = brute force", k, nrow(ref)),
           as.numeric(!same), 0.5)
  }
}

# --- calib ---------------------------------------------------------------------------

if ("calib" %in% only) {
  cat("\n== calib: moments at N = 200,000 (seed 1) against the targets ==\n")
  for (id in c("D3", "D3-claims", "D4-census")) {
    dgp <- get_dgp(id)
    tg <- dgp$targets
    for (v in dgp$variants) {
      t0 <- now()
      sim <- dgp$simulate(N = 200000L, seed = 1L, variant = v, engine = "cpp")
      m <- dgp$moments(sim)
      cat(sprintf("-- %s variant %s: %d rows, simulated and summarized in %.1f s\n", id, v,
                  nrow(sim$data), now() - t0))
      i <- match(names(m), tg$moment)
      print(data.frame(moment = names(m), value = signif(unname(m), 4), lo = tg$lo[i],
                       hi = tg$hi[i],
                       within = ifelse(is.na(i), NA, m >= tg$lo[i] & m <= tg$hi[i])),
            row.names = FALSE)
      for (k in tg$moment) {
        lo <- tg$lo[tg$moment == k]
        hi <- tg$hi[tg$moment == k]
        report(sprintf("%s %s: %s in [%.2f, %.2f]", id, v, k, lo, hi), m[[k]], hi,
               m[[k]] >= lo && m[[k]] <= hi)
      }
      if (startsWith(id, "D3")) {
        report(sprintf("%s %s: the sets connect every hospital (components)", id, v),
               m[["components"]], 1.5)
      }
      # Expected counts per alternative and the expected number never chosen
      # (sum_j exp(-n s_j)) at the families' sizes.
      sh <- sim$meta$expected / if (startsWith(id, "D3")) sim$meta$A else sim$meta$N
      sizes <- switch(id, D3 = c(5000, 20000, 80000), `D3-claims` = 625000 * sim$meta$A / 2e5,
                      `D4-census` = 1e6)
      print(do.call(rbind, lapply(sizes, function(n) {
        e <- n * sh
        data.frame(design = id, variant = v, n = round(n), J = length(sh),
                   J_over_sqrtN = round(length(sh) / sqrt(n), 3), count_min = round(min(e), 2),
                   count_p10 = round(unname(stats::quantile(e, 0.1)), 1),
                   count_med = round(stats::median(e), 1), exp_never = round(sum(exp(-e)), 2))
      })), row.names = FALSE)
    }
  }
}

# --- full ----------------------------------------------------------------------------

if ("full" %in% only) {
  for (des in trimws(strsplit(args$full, ",", fixed = TRUE)[[1]])) {
    id <- c(claims = "D3-claims", census = "D4-census")[[des]]
    N <- c(claims = 625000L, census = 1000000L)[[des]]
    for (v in c("mnl", "default")) {
      gc(reset = TRUE)
      t0 <- now()
      sim <- get_dgp(id)$simulate(N = N, seed = 1L, variant = v, engine = "cpp")
      t1 <- now()
      sim <- geo_drop(sim)
      t2 <- now()
      g <- gc()
      cat(sprintf(paste("full %s %s: N = %d, %d situations, %d rows, %d columns; generated in",
                        "%.1f s, drop %.1f s (%d dropped); R heap peak %.0f MB, data %.0f MB\n"),
                  id, v, N, if (startsWith(id, "D3")) sim$meta$A else N, nrow(sim$data),
                  ncol(sim$data), t1 - t0, t2 - t1, sim$meta$drop[["alts"]],
                  sum(g[, which(colnames(g) == "max used") + 1L]),
                  as.numeric(utils::object.size(sim$data)) / 2^20))
      report(sprintf("full %s %s: generated within 10 minutes", id, v), t1 - t0, 600)
      rm(sim)
    }
  }
}

# --- smoke ---------------------------------------------------------------------------

if ("smoke" %in% only) {
  cat("\n== smoke: run_rep() at every cell x method ==\n")
  S_over <- if (is.null(args$S)) list() else {
    kv <- strsplit(trimws(strsplit(args$S, ",", fixed = TRUE)[[1]]), ":", fixed = TRUE)
    stats::setNames(lapply(kv, function(x) as.integer(x[2])), vapply(kv, `[`, "", 1L))
  }
  tok <- trimws(strsplit(args$smoke, ",", fixed = TRUE)[[1]])
  tab <- list()
  for (tk in tok) {
    p <- strsplit(tk, ":", fixed = TRUE)[[1]]
    fam <- load_family(p[1], b$L)
    if (fam$id == "F7") {
      mxl <- fam$cells$cell %in% c("claims", "census")
      fam$cells$N[mxl] <- as.integer(fam$cells$N[mxl] / args$f7_div)
      fam$cells$S[mxl] <- as.integer(args$f7_S)
      fam$cells$N[!mxl] <- as.integer(fam$cells$N[!mxl] / args$f7_mnl_div)
    }
    for (k in names(S_over)) fam$cells$S[fam$cells$cell == k] <- S_over[[k]]
    cells <- if (length(p) >= 2L) p[2] else fam$cells$cell
    for (cl in cells) {
      ms <- cell_methods(fam, cl)
      if (length(p) >= 3L) ms <- intersect(ms, p[3])
      bottom <- fam$id == "F6" && fam$cells$N[fam$cells$cell == cl] == 5000L
      reps <- if (bottom) c(1L, 2L, 10L) else 1L
      for (r in reps) {
        t0 <- now()
        recs <- run_rep(fam, cl, r, ms)
        for (rec in recs) {
          valid <- tryCatch({
            do.call(new_rep_record, rec[c("family", "cell", "method", "rep", "seed", "ok", "error",
                                          "warnings", "est", "vcov_c", "joint_wald", "lr",
                                          "derived", "diag", "meta")])
            TRUE
          }, error = function(e) FALSE)
          w <- unique(rec$warnings$class)
          st <- rec$diag$start
          row <- fam$cells[fam$cells$cell == cl, ]
          tab[[length(tab) + 1L]] <- data.frame(
            family = rec$family, cell = rec$cell, method = rec$method, rep = rec$rep,
            N = row$N, S = row$S, rows = rec$diag$rows %||% NA, ok = rec$ok, valid = valid,
            sim_sec = round(rec$diag$sim_sec, 1), seconds = round(rec$diag$unit_sec, 1),
            fit_sec = round(rec$diag$time_sec %||% NA_real_, 1),
            iterations = rec$diag$iterations %||% NA, code = rec$diag$code %||% NA,
            decrement = signif(rec$diag$decrement %||% NA_real_, 2),
            start = if (is.null(st)) "" else sprintf("%s %.2g", st$winner, st$dlogL %||% NA),
            dropped = rec$diag$drop[["alts"]] %||% NA,
            heap_mb = round(rec$diag$peak_mem_mb[["heap_mb"]]),
            rss_mb = round(rec$diag$peak_mem_mb[["rss_mb"]]),
            ms_dlogL = signif(rec$diag$multistart$dlogL %||% NA_real_, 2),
            joint_p = signif(rec$joint_wald$p %||% NA_real_, 2),
            warnings = if (length(w)) paste(w, collapse = " ") else "",
            error = substr(rec$error %||% "", 1, 80), stringsAsFactors = FALSE)
          report(sprintf("%s/%s/%s rep %d: ok and valid", rec$family, rec$cell, rec$method,
                         rec$rep), as.numeric(!(rec$ok && valid)), 0.5)
          if (rec$ok && rec$rep == 1L) {
            e <- rec$est
            z <- (e$est - e$truth) / e$se_hessian
            cat(sprintf("      common z: %s\n", paste(sprintf("%s %.2f", e$name[e$block != "asc"],
                                                              z[e$block != "asc"]),
                                                      collapse = ", ")))
            za <- z[e$block == "asc"]
            cat(sprintf("      constants: %d, z mean %.2f sd %.2f; pooled 95%% coverage %.3f\n",
                        length(za), mean(za, na.rm = TRUE), stats::sd(za, na.rm = TRUE),
                        mean(abs(za) < stats::qnorm(0.975), na.rm = TRUE)))
            cat(sprintf("      derived: %s\n", paste(sprintf(
              "%s %.4g (truth %.4g)", rec$derived$name, rec$derived$est, rec$derived$truth),
              collapse = "; ")))
          }
        }
        cat(sprintf("      %s %s rep %d: %.1f s\n", fam$id, cl, r, now() - t0))
        flush(stdout())
      }
    }
  }
  print(do.call(rbind, tab), row.names = FALSE)
}

cat(sprintf("\n%d checks, %d FAIL; %.1f s\n", length(verdicts), sum(!verdicts),
            proc.time()[["elapsed"]] - t_start))
if (any(!verdicts)) quit(status = 1L)
