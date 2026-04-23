# _validation/mxl_monte_carlo.R
#
# Driver for the Mixed Logit Monte Carlo validation suite.
#
# Modes:
#   QUICK=TRUE   ~90 min pre-flight; small R, trimmed N-grid and S-grid
#   (default)    full overnight run; see _validation/README.md
#
# Parallelism:
#   OMP_NUM_THREADS=2 externally; future::plan(multisession) over floor(cores/2)
#   workers to avoid oversubscription. Each run_mxlogit() uses OpenMP internally.
#
# Output artifacts (written under _validation/output/):
#   mc_<scenario>.rds
#   asymptotics_<scenario>_raw.csv
#   asymptotics_<scenario>_natural.csv
#   plot_qq_<scenario>.png
#   plot_consistency.png           (Scenario A only)
#   plot_coverage.png              (Scenario A only)
#   plot_sim_bias.png              (Scenario D only)
#   plot_se_ratios.png
#   plot_lr_cdf.png                (Scenario A only)
#   REPORT.md

t0 <- Sys.time()
suppressPackageStartupMessages({
  library(choicer)
  library(data.table)
})

args_env <- function(name, default) {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v) || v == "") default else v
}
QUICK <- isTRUE(tolower(args_env("QUICK", "false")) %in% c("true", "1", "t"))

# Mirror the parallelism guardrails from the plan.
if (Sys.getenv("OMP_NUM_THREADS") == "") Sys.setenv(OMP_NUM_THREADS = "2")
try(choicer:::set_num_threads(as.integer(Sys.getenv("OMP_NUM_THREADS"))),
    silent = TRUE)
n_workers <- max(1L, floor(parallel::detectCores() / 2))

PARALLEL <- requireNamespace("future.apply", quietly = TRUE) &&
  requireNamespace("future", quietly = TRUE)
if (PARALLEL) {
  future::plan(future::multisession, workers = n_workers)
} else {
  message("future.apply / future not available; running serially.")
}

options(choicer.quick_mode = QUICK)

# Resolve repo root / source helpers.
this_dir <- tryCatch(
  dirname(normalizePath(sys.frames()[[1]]$ofile %||% "_validation/mxl_monte_carlo.R")),
  error = function(e) "_validation"
)
if (!file.exists(file.path(this_dir, "mxl_mc_helpers.R"))) {
  # Fallback when sourced from repo root.
  this_dir <- "_validation"
}
source(file.path(this_dir, "mxl_mc_helpers.R"))

out_dir <- file.path(this_dir, "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("[driver] QUICK=", QUICK, " workers=", n_workers,
        " OMP_NUM_THREADS=", Sys.getenv("OMP_NUM_THREADS"))

# ---------------------------------------------------------------------------
# Helper: run one MC scenario and assemble outputs
# ---------------------------------------------------------------------------

run_scenario <- function(scn, tag = scn$id, base_seed = 20260423L,
                         N_override = NULL, S_override = NULL,
                         with_lr = FALSE, lr_restriction = NULL,
                         hessian_first_rep = TRUE) {
  # Choose sim_fun and fit_fun based on scenario id and any override.
  if (!is.null(N_override)) {
    sim_fun <- scn$sim_fun_for_N(N_override)
  } else {
    sim_fun <- scn$sim_fun
  }
  if (!is.null(S_override)) {
    fit_fun <- scn$fit_fun_for_S(S_override)
  } else {
    fit_fun <- scn$fit_fun
  }

  R <- if (!is.null(scn$R_map) && !is.null(N_override)) {
    scn$R_map[[as.character(N_override)]]
  } else {
    scn$R
  }

  run_one <- function(r) {
    seed_r <- base_seed + r - 1L
    sim <- sim_fun(seed = seed_r)
    tic <- Sys.time()
    fit <- tryCatch(fit_fun(sim), error = function(e) e)
    toc <- Sys.time()
    elapsed <- as.numeric(difftime(toc, tic, units = "secs"))

    if (inherits(fit, "error")) {
      return(list(
        rep_id = r, seed = seed_r, error = conditionMessage(fit),
        time_sec = elapsed, reps = NULL, lr = NULL, hess = NULL
      ))
    }

    reps <- tryCatch(augment_replication(fit, sim), error = function(e) NULL)

    lr_rec <- NULL
    if (with_lr && !is.null(lr_restriction)) {
      restr <- restricted_fit(sim, scn, lr_restriction)
      lr_rec <- data.table(
        rep_id = r, seed = seed_r,
        loglik_full = as.numeric(stats::logLik(fit)),
        loglik_restricted = restr$loglik_restricted,
        lr_stat = 2 * (as.numeric(stats::logLik(fit)) - restr$loglik_restricted),
        restricted_status = restr$convergence
      )
    }

    hess_rec <- NULL
    if (hessian_first_rep && r == 1L) {
      hg <- hessian_agreement(fit, sim)
      hess_rec <- data.table(
        rep_id = r, max_rel_err = hg$max_rel_err, note = hg$note
      )
    }

    list(
      rep_id = r, seed = seed_r, error = NA_character_,
      time_sec = elapsed, loglik = as.numeric(stats::logLik(fit)),
      converged = choicer:::.is_converged(fit),
      reps = reps, lr = lr_rec, hess = hess_rec
    )
  }

  message(sprintf("[%s] R=%d %s%s",
                  tag, R,
                  if (!is.null(N_override)) paste0(" N=", N_override) else "",
                  if (!is.null(S_override)) paste0(" S=", S_override) else ""))

  if (PARALLEL) {
    results <- future.apply::future_lapply(
      seq_len(R), run_one, future.seed = TRUE
    )
  } else {
    results <- lapply(seq_len(R), run_one)
  }

  # Assemble long data.table of replications, augmented with rep metadata.
  reps_list <- lapply(results, function(res) {
    if (is.null(res$reps)) return(NULL)
    r <- res$reps
    r[, `:=`(
      rep_id = res$rep_id, seed = res$seed,
      loglik = res$loglik, converged = res$converged,
      time_sec = res$time_sec, error = res$error
    )]
    r
  })
  replications <- rbindlist(reps_list, use.names = TRUE, fill = TRUE)
  if (!nrow(replications)) {
    warning(sprintf("[%s] no successful replications", tag))
    return(list(mc_raw = NULL, mc_natural = NULL, lr = NULL, hess = NULL,
                conv_rate = NA_real_))
  }

  # Split by scale and box into choicer_mc.
  mk_mc <- function(scale_tag) {
    r <- replications[scale == scale_tag]
    if (!nrow(r)) return(NULL)
    r2 <- copy(r); r2[, scale := NULL]
    structure(
      list(
        replications = r2,
        meta = list(
          R = R, seed = base_seed, seeds = base_seed + seq_len(R) - 1L,
          parallel = PARALLEL, timestamp = Sys.time(),
          scenario = tag, scale = scale_tag,
          N_override = N_override, S_override = S_override
        )
      ),
      class = "choicer_mc"
    )
  }

  mc_raw <- mk_mc("raw")
  mc_natural <- mk_mc("natural")

  lr_tbl <- rbindlist(lapply(results, `[[`, "lr"),
                      use.names = TRUE, fill = TRUE)
  hess_tbl <- rbindlist(lapply(results, `[[`, "hess"),
                        use.names = TRUE, fill = TRUE)

  conv_rate <- mean(vapply(results, function(x) isTRUE(x$converged),
                           logical(1)))

  list(mc_raw = mc_raw, mc_natural = mc_natural,
       lr = lr_tbl, hess = hess_tbl, conv_rate = conv_rate)
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

scenario_ids <- c("A", "B", "B2", "C", "D", "F")
scenario_results <- list()

for (id in scenario_ids) {
  t_scn <- Sys.time()
  scn <- scenario_spec(id, quick = QUICK)
  message("[scenario ", id, "] ", scn$purpose)

  if (id == "A") {
    # N-grid; augment each arm with LR test + Hessian-agreement first rep.
    mc_list_raw <- list(); mc_list_nat <- list()
    lr_list <- list(); hess_list <- list()
    for (N in scn$N_grid) {
      r <- run_scenario(scn, tag = paste0("A_N", N), N_override = N,
                        with_lr = TRUE,
                        lr_restriction = scn$meta$lr_restriction,
                        hessian_first_rep = TRUE)
      if (!is.null(r$mc_raw)) {
        r$mc_raw$meta$N <- N
        mc_list_raw[[as.character(N)]] <- r$mc_raw
      }
      if (!is.null(r$mc_natural)) {
        r$mc_natural$meta$N <- N
        mc_list_nat[[as.character(N)]] <- r$mc_natural
      }
      if (!is.null(r$lr) && nrow(r$lr)) {
        r$lr[, N := N]
        lr_list[[as.character(N)]] <- r$lr
      }
      if (!is.null(r$hess) && nrow(r$hess)) {
        r$hess[, N := N]
        hess_list[[as.character(N)]] <- r$hess
      }
    }

    # For the per-scenario asymptotics table we use the largest-N arm.
    biggest <- as.character(max(scn$N_grid))
    asymp_raw <- if (!is.null(mc_list_raw[[biggest]])) mc_asymptotics(mc_list_raw[[biggest]]) else NULL
    asymp_nat <- if (!is.null(mc_list_nat[[biggest]])) mc_asymptotics(mc_list_nat[[biggest]]) else NULL

    # Persist.
    saveRDS(mc_list_raw, file.path(out_dir, sprintf("mc_%s_raw.rds", id)))
    saveRDS(mc_list_nat, file.path(out_dir, sprintf("mc_%s_natural.rds", id)))
    lr_all <- rbindlist(lr_list, use.names = TRUE, fill = TRUE)
    hess_all <- rbindlist(hess_list, use.names = TRUE, fill = TRUE)
    saveRDS(list(lr = lr_all, hess = hess_all),
            file.path(out_dir, sprintf("aux_%s.rds", id)))
    if (!is.null(asymp_raw))
      data.table::fwrite(asymp_raw, file.path(out_dir, sprintf("asymptotics_%s_raw.csv", id)))
    if (!is.null(asymp_nat))
      data.table::fwrite(asymp_nat, file.path(out_dir, sprintf("asymptotics_%s_natural.csv", id)))

    # QQ plot on the largest-N arm.
    if (!is.null(mc_list_raw[[biggest]])) {
      plot_qq(mc_list_raw[[biggest]]$replications,
              scenario_id = id,
              outfile = file.path(out_dir, sprintf("plot_qq_%s.png", id)))
    }

    # Consistency / coverage plots across N.
    cons_tbl <- rbindlist(lapply(names(mc_list_raw), function(nN) {
      mc <- mc_list_raw[[nN]]
      a <- tryCatch(mc_asymptotics(mc), error = function(e) NULL)
      if (is.null(a)) return(NULL)
      data.table(N = as.integer(nN), parameter = a$parameter,
                 group = a$group, abs_bias = abs(a$bias),
                 rmse = sqrt(a$bias^2 + a$sd_emp^2))
    }), use.names = TRUE, fill = TRUE)
    if (nrow(cons_tbl)) {
      plot_consistency(cons_tbl,
                       outfile = file.path(out_dir, "plot_consistency.png"))
    }
    cov_tbl <- rbindlist(lapply(names(mc_list_raw), function(nN) {
      mc <- mc_list_raw[[nN]]
      a <- tryCatch(mc_asymptotics(mc), error = function(e) NULL)
      if (is.null(a)) return(NULL)
      data.table(N = as.integer(nN), parameter = a$parameter,
                 cov95 = a$cov95, cov95_lower = a$cov95_lower,
                 cov95_upper = a$cov95_upper)
    }), use.names = TRUE, fill = TRUE)
    if (nrow(cov_tbl)) {
      plot_coverage(cov_tbl,
                    outfile = file.path(out_dir, "plot_coverage.png"))
    }

    # LR CDF plot + KS + size test.
    if (nrow(lr_all)) {
      plot_lr_cdf(lr_all$lr_stat,
                  outfile = file.path(out_dir, "plot_lr_cdf.png"),
                  df = 1L)
    }

    lr_clean <- lr_all$lr_stat[is.finite(lr_all$lr_stat) & lr_all$lr_stat >= 0]
    lr_ks_p <- if (length(lr_clean) >= 4) {
      suppressWarnings(ks.test(lr_clean, "pchisq", df = 1)$p.value)
    } else NA_real_
    lr_hits <- sum(lr_clean > qchisq(0.95, 1))
    lr_size <- if (length(lr_clean)) lr_hits / length(lr_clean) else NA_real_
    lr_wilson <- choicer:::.wilson_ci(
      p = lr_size, n = length(lr_clean), level = 0.95
    )

    hess_err <- if (nrow(hess_all)) median(hess_all$max_rel_err, na.rm = TRUE) else NA_real_

    # Convergence rate across all arms (fraction of reps that converged).
    conv_rate <- if (length(mc_list_raw)) {
      reps_all <- rbindlist(lapply(mc_list_raw, function(mc) mc$replications),
                            use.names = TRUE, fill = TRUE)
      mean(reps_all[, any(converged, na.rm = TRUE), by = rep_id]$V1)
    } else NA_real_

    scenario_results[[id]] <- list(
      asymptotics_raw = asymp_raw,
      asymptotics_natural = asymp_nat,
      meta = list(purpose = scn$purpose),
      extras = list(
        lr_ks_p = lr_ks_p, lr_size = lr_size,
        lr_size_lower = lr_wilson$lower, lr_size_upper = lr_wilson$upper,
        hessian_err = hess_err,
        conv_rate = conv_rate
      )
    )

  } else if (id == "D") {
    # S-grid at fixed N for simulation-bias regression.
    mc_list_raw <- list(); mc_list_nat <- list()
    hess_list <- list()
    for (S in scn$S_grid) {
      r <- run_scenario(scn, tag = paste0("D_S", S), S_override = S,
                        with_lr = FALSE, hessian_first_rep = TRUE)
      if (!is.null(r$mc_raw)) {
        r$mc_raw$meta$S <- S
        mc_list_raw[[as.character(S)]] <- r$mc_raw
      }
      if (!is.null(r$mc_natural)) {
        r$mc_natural$meta$S <- S
        mc_list_nat[[as.character(S)]] <- r$mc_natural
      }
      if (!is.null(r$hess) && nrow(r$hess)) {
        r$hess[, S := S]
        hess_list[[as.character(S)]] <- r$hess
      }
    }
    biggest <- as.character(max(scn$S_grid))
    asymp_raw <- if (!is.null(mc_list_raw[[biggest]])) mc_asymptotics(mc_list_raw[[biggest]]) else NULL
    asymp_nat <- if (!is.null(mc_list_nat[[biggest]])) mc_asymptotics(mc_list_nat[[biggest]]) else NULL

    saveRDS(mc_list_raw, file.path(out_dir, sprintf("mc_%s_raw.rds", id)))
    saveRDS(mc_list_nat, file.path(out_dir, sprintf("mc_%s_natural.rds", id)))
    hess_all <- rbindlist(hess_list, use.names = TRUE, fill = TRUE)
    saveRDS(list(hess = hess_all), file.path(out_dir, sprintf("aux_%s.rds", id)))
    if (!is.null(asymp_raw))
      data.table::fwrite(asymp_raw, file.path(out_dir, sprintf("asymptotics_%s_raw.csv", id)))
    if (!is.null(asymp_nat))
      data.table::fwrite(asymp_nat, file.path(out_dir, sprintf("asymptotics_%s_natural.csv", id)))

    if (!is.null(mc_list_raw[[biggest]])) {
      plot_qq(mc_list_raw[[biggest]]$replications,
              scenario_id = id,
              outfile = file.path(out_dir, sprintf("plot_qq_%s.png", id)))
    }

    simbias_tbl <- rbindlist(lapply(names(mc_list_raw), function(sS) {
      mc <- mc_list_raw[[sS]]
      a <- tryCatch(mc_asymptotics(mc), error = function(e) NULL)
      if (is.null(a)) return(NULL)
      data.table(S = as.integer(sS), parameter = a$parameter,
                 abs_bias = abs(a$bias))
    }), use.names = TRUE, fill = TRUE)
    if (nrow(simbias_tbl)) {
      plot_sim_bias(simbias_tbl, outfile = file.path(out_dir, "plot_sim_bias.png"))
    }

    hess_err <- if (nrow(hess_all)) median(hess_all$max_rel_err, na.rm = TRUE) else NA_real_
    scenario_results[[id]] <- list(
      asymptotics_raw = asymp_raw,
      asymptotics_natural = asymp_nat,
      meta = list(purpose = scn$purpose),
      extras = list(hessian_err = hess_err)
    )

  } else {
    # Single-(N, S) scenarios.
    r <- run_scenario(scn, tag = id, with_lr = FALSE,
                      hessian_first_rep = TRUE)
    asymp_raw <- if (!is.null(r$mc_raw)) mc_asymptotics(r$mc_raw) else NULL
    asymp_nat <- if (!is.null(r$mc_natural)) mc_asymptotics(r$mc_natural) else NULL

    if (!is.null(r$mc_raw))
      saveRDS(r$mc_raw, file.path(out_dir, sprintf("mc_%s_raw.rds", id)))
    if (!is.null(r$mc_natural))
      saveRDS(r$mc_natural, file.path(out_dir, sprintf("mc_%s_natural.rds", id)))
    if (!is.null(r$hess) && nrow(r$hess))
      saveRDS(list(hess = r$hess),
              file.path(out_dir, sprintf("aux_%s.rds", id)))
    if (!is.null(asymp_raw))
      data.table::fwrite(asymp_raw,
                         file.path(out_dir, sprintf("asymptotics_%s_raw.csv", id)))
    if (!is.null(asymp_nat))
      data.table::fwrite(asymp_nat,
                         file.path(out_dir, sprintf("asymptotics_%s_natural.csv", id)))

    if (!is.null(r$mc_raw)) {
      plot_qq(r$mc_raw$replications,
              scenario_id = id,
              outfile = file.path(out_dir, sprintf("plot_qq_%s.png", id)))
    }

    hess_err <- if (!is.null(r$hess) && nrow(r$hess)) {
      median(r$hess$max_rel_err, na.rm = TRUE)
    } else NA_real_

    scenario_results[[id]] <- list(
      asymptotics_raw = asymp_raw,
      asymptotics_natural = asymp_nat,
      meta = list(purpose = scn$purpose),
      extras = list(hessian_err = hess_err, conv_rate = r$conv_rate)
    )
  }

  message(sprintf("[scenario %s] done in %.1f min", id,
                  as.numeric(difftime(Sys.time(), t_scn, units = "mins"))))
}

# ---------------------------------------------------------------------------
# Cross-scenario plots and final report
# ---------------------------------------------------------------------------

# SE ratio plot across scenarios on the raw scale.
se_rows <- list()
for (id in names(scenario_results)) {
  a <- scenario_results[[id]]$asymptotics_raw
  if (is.null(a) || !nrow(a)) next
  se_rows[[id]] <- data.table(
    scenario = id, parameter = a$parameter,
    hess_ratio = a$se_ratio,
    bhhh_ratio = NA_real_  # BHHH SE ratio is computed separately; see below
  )
}
se_tbl <- rbindlist(se_rows, use.names = TRUE, fill = TRUE)
if (nrow(se_tbl)) {
  plot_se_ratios(se_tbl, outfile = file.path(out_dir, "plot_se_ratios.png"))
}

build_report(scenario_results, outfile = file.path(out_dir, "REPORT.md"))

elapsed_total <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
message(sprintf("[driver] done in %.1f min total", elapsed_total))
message("Report written to ", file.path(out_dir, "REPORT.md"))
