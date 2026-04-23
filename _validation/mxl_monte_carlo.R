# _validation/mxl_monte_carlo.R
#
# Driver for the Mixed Logit Monte Carlo validation suite.
#
# Modes:
#   QUICK=TRUE   ~90 min pre-flight; small R, trimmed N-grid and S-grid
#   (default)    full overnight run; see _validation/README.md
#
# Run tagging and checkpointing:
#   RUN_TAG=<name>     Resume into or create the folder
#                      _validation/output/<name>/. If the folder exists,
#                      any scenario with an existing mc_<id>_raw.rds is
#                      skipped and rehydrated from disk (useful after a
#                      spot preemption).
#   SCENARIOS=A,B,F    Only run this comma-separated subset. Combined with
#                      RUN_TAG=<existing> to re-run specific scenarios
#                      into an existing folder.
#
#   Default behavior (no RUN_TAG): fresh folder named by timestamp
#   YYYYMMDD-HHMMSS. _validation/output/latest always points at the most
#   recent run.
#
# Parallelism:
#   OMP_NUM_THREADS=2 externally; future::plan(multisession) over floor(cores/2)
#   workers to avoid oversubscription. Each run_mxlogit() uses OpenMP internally.
#
# Output artifacts (written under _validation/output/<RUN_TAG>/):
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
  # Prefer devtools::load_all() when running from a development checkout so
  # the driver sees the in-tree source (including any unreleased edits to
  # mc_asymptotics() etc.). Fall back to library(choicer) when running from
  # an installed package (e.g., on a CI machine without devtools).
  loaded <- FALSE
  if (requireNamespace("devtools", quietly = TRUE)) {
    pkg_root <- tryCatch(rprojroot::find_package_root_file(), error = function(e) NULL)
    if (!is.null(pkg_root) && file.exists(file.path(pkg_root, "DESCRIPTION"))) {
      devtools::load_all(pkg_root, quiet = TRUE)
      loaded <- TRUE
    }
  }
  if (!loaded) library(choicer)
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

# Run tag: controls the output subfolder under _validation/output/.
#   - Unset (default): fresh run, use timestamp YYYYMMDD-HHMMSS.
#   - RUN_TAG=<name>: resume into an existing tagged folder (so partial
#     outputs from a prior spot-preempted run are reused) OR create a
#     named fresh folder.
# The `latest` symlink / LATEST pointer file in _validation/output always
# points at the most recent run so ad-hoc tooling works without a tag.
root_out <- file.path(this_dir, "output")
if (!dir.exists(root_out)) dir.create(root_out, recursive = TRUE)
run_tag <- args_env("RUN_TAG", format(Sys.time(), "%Y%m%d-%H%M%S"))
out_dir <- file.path(root_out, run_tag)
resuming <- dir.exists(out_dir)
if (!resuming) dir.create(out_dir, recursive = TRUE)

# Update the `latest` pointer. Prefer a symlink; fall back to a text file
# on filesystems that reject symlinks (rare on *nix, but e.g. some S3-
# backed mounts). Guard against edge cases where the target exists as a
# regular file (e.g., old flat layout) or as a dangling symlink.
latest_link <- file.path(root_out, "latest")
existing_target <- suppressWarnings(Sys.readlink(latest_link))
if (file.exists(latest_link) || (!is.na(existing_target) && nzchar(existing_target))) {
  try(unlink(latest_link, force = TRUE), silent = TRUE)
}
symlink_ok <- tryCatch(
  isTRUE(file.symlink(run_tag, latest_link)),
  error = function(e) FALSE, warning = function(w) FALSE
)
if (!symlink_ok) writeLines(run_tag, file.path(root_out, "LATEST"))

# Scenario filter: SCENARIOS=A,B,B2 runs the intersection with the six
# defined scenarios. Unset = run all. Combines with RUN_TAG=<existing> to
# re-run individual scenarios into an existing folder after a preemption.
scenario_filter <- args_env("SCENARIOS", "")
scenario_filter <- if (nzchar(scenario_filter)) {
  toupper(trimws(strsplit(scenario_filter, ",")[[1]]))
} else NULL

message("[driver] QUICK=", QUICK, " workers=", n_workers,
        " OMP_NUM_THREADS=", Sys.getenv("OMP_NUM_THREADS"),
        " RUN_TAG=", run_tag, " resuming=", resuming,
        " out_dir=", out_dir)
if (!is.null(scenario_filter)) {
  message("[driver] SCENARIOS filter: ", paste(scenario_filter, collapse = ","))
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Compute mc_asymptotics() against the BHHH per-rep SE column (`se_bhhh`,
# populated by augment_replication() in mxl_mc_helpers.R). Returns NULL
# gracefully if the column is absent or every value is NA -- BHHH failure
# on a subset of reps is expected on weak-identification scenarios.
.safe_mc_asymptotics_bhhh <- function(mc) {
  if (is.null(mc) || !inherits(mc, "choicer_mc")) return(NULL)
  if (!"se_bhhh" %in% names(mc$replications)) return(NULL)
  if (all(is.na(mc$replications$se_bhhh))) return(NULL)
  tryCatch(mc_asymptotics(mc, se_col = "se_bhhh"),
           error = function(e) { message("BHHH asymptotics failed: ",
                                         conditionMessage(e)); NULL })
}

# Checkpointing: a scenario is "done" if its primary RDS output already
# exists in out_dir. We use mc_<id>_raw.rds as the completion sentinel
# because every scenario writes it last (after the per-rep loop completes).
.scenario_is_done <- function(id, out_dir) {
  file.exists(file.path(out_dir, sprintf("mc_%s_raw.rds", id)))
}

# Scenario selection: apply SCENARIOS filter AND the RUN_TAG-based
# completed-scenario skip. Returns TRUE if this scenario should run.
.should_run_scenario <- function(id, out_dir, filter = NULL) {
  if (!is.null(filter) && !(id %in% filter)) {
    message(sprintf("[scenario %s] skipped (not in SCENARIOS filter)", id))
    return(FALSE)
  }
  if (.scenario_is_done(id, out_dir)) {
    message(sprintf("[scenario %s] skipped (already complete in %s)",
                    id, basename(out_dir)))
    return(FALSE)
  }
  TRUE
}

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
  if (!.should_run_scenario(id, out_dir, filter = scenario_filter)) {
    # Resume path: rehydrate outputs from disk so cross-scenario aggregation
    # (SE-ratio plot, REPORT.md) still sees them. If the RDS is missing
    # we silently skip aggregation for this scenario.
    if (.scenario_is_done(id, out_dir)) {
      mc_raw_path <- file.path(out_dir, sprintf("mc_%s_raw.rds", id))
      mc_nat_path <- file.path(out_dir, sprintf("mc_%s_natural.rds", id))
      mc_raw <- tryCatch(readRDS(mc_raw_path), error = function(e) NULL)
      mc_nat <- if (file.exists(mc_nat_path)) {
        tryCatch(readRDS(mc_nat_path), error = function(e) NULL)
      } else NULL
      # For Scenario A / D the on-disk RDS is a named list of per-arm
      # choicer_mc objects; pick the biggest arm for the summary table,
      # matching the fresh-run behavior further down.
      pick_biggest <- function(mc_list) {
        if (is.null(mc_list) || !length(mc_list)) return(NULL)
        if (inherits(mc_list, "choicer_mc")) return(mc_list)
        mc_list[[as.character(max(as.numeric(names(mc_list))))]]
      }
      mc_big_raw <- pick_biggest(mc_raw)
      mc_big_nat <- pick_biggest(mc_nat)
      scenario_results[[id]] <- list(
        asymptotics_raw      = if (!is.null(mc_big_raw)) mc_asymptotics(mc_big_raw) else NULL,
        asymptotics_raw_bhhh = .safe_mc_asymptotics_bhhh(mc_big_raw),
        asymptotics_natural  = if (!is.null(mc_big_nat)) mc_asymptotics(mc_big_nat) else NULL,
        meta   = list(purpose = scenario_spec(id, quick = QUICK)$purpose,
                      resumed = TRUE),
        extras = list()
      )
    }
    next
  }

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
    asymp_raw_bhhh <- .safe_mc_asymptotics_bhhh(mc_list_raw[[biggest]])

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

    # Per-arm convergence rate: group by (N_label, rep_id) so reps that
    # converged in one N-arm but not another are counted per-arm, not
    # collapsed across arms via any().
    conv_rate <- if (length(mc_list_raw)) {
      reps_all <- rbindlist(
        lapply(names(mc_list_raw), function(nN) {
          r <- data.table::copy(mc_list_raw[[nN]]$replications)
          r[, N_label := nN]
          r
        }),
        use.names = TRUE, fill = TRUE
      )
      mean(reps_all[, any(converged, na.rm = TRUE), by = .(N_label, rep_id)]$V1)
    } else NA_real_

    scenario_results[[id]] <- list(
      asymptotics_raw = asymp_raw,
      asymptotics_raw_bhhh = asymp_raw_bhhh,
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
    asymp_raw_bhhh <- .safe_mc_asymptotics_bhhh(mc_list_raw[[biggest]])

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

    # Claim 5: per-parameter OLS of |bias| on 1/S. Expected slope > 0 and
    # t-stat meaningfully away from zero if simulation bias is O(1/S).
    simbias_slopes <- if (nrow(simbias_tbl) && length(unique(simbias_tbl$S)) >= 2) {
      simbias_tbl[, {
        d <- data.frame(abs_bias = abs_bias, inv_S = 1 / S)
        d <- d[is.finite(d$abs_bias) & is.finite(d$inv_S), ]
        if (nrow(d) >= 2) {
          fit <- stats::lm(abs_bias ~ inv_S, data = d)
          cf <- stats::coef(summary(fit))
          if (nrow(cf) >= 2) {
            list(slope = cf[2, 1], slope_se = cf[2, 2],
                 slope_t = cf[2, 3], slope_p = cf[2, 4])
          } else list(slope = NA_real_, slope_se = NA_real_,
                      slope_t = NA_real_, slope_p = NA_real_)
        } else list(slope = NA_real_, slope_se = NA_real_,
                    slope_t = NA_real_, slope_p = NA_real_)
      }, by = parameter]
    } else NULL

    hess_err <- if (nrow(hess_all)) median(hess_all$max_rel_err, na.rm = TRUE) else NA_real_
    scenario_results[[id]] <- list(
      asymptotics_raw = asymp_raw,
      asymptotics_raw_bhhh = asymp_raw_bhhh,
      asymptotics_natural = asymp_nat,
      meta = list(purpose = scn$purpose),
      extras = list(hessian_err = hess_err,
                    simbias_slopes = simbias_slopes)
    )

  } else {
    # Single-(N, S) scenarios.
    r <- run_scenario(scn, tag = id, with_lr = FALSE,
                      hessian_first_rep = TRUE)
    asymp_raw <- if (!is.null(r$mc_raw)) mc_asymptotics(r$mc_raw) else NULL
    asymp_nat <- if (!is.null(r$mc_natural)) mc_asymptotics(r$mc_natural) else NULL
    asymp_raw_bhhh <- .safe_mc_asymptotics_bhhh(r$mc_raw)

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
      asymptotics_raw_bhhh = asymp_raw_bhhh,
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

# SE ratio plot across scenarios on the raw scale. Both Hessian-based and
# BHHH-based ratios are surfaced so Claim 4 (information-matrix equality)
# is fully adjudicated: under correct specification, both should converge
# to 1.
se_rows <- list()
for (id in names(scenario_results)) {
  a <- scenario_results[[id]]$asymptotics_raw
  if (is.null(a) || !nrow(a)) next
  ab <- scenario_results[[id]]$asymptotics_raw_bhhh
  # Align BHHH rows to the Hessian parameter order (same (parameter, group)
  # keys since they come from the same choicer_mc); left-join by parameter.
  bhhh_ratio <- if (!is.null(ab) && nrow(ab)) {
    ab[match(a$parameter, ab$parameter), se_ratio]
  } else rep(NA_real_, nrow(a))
  se_rows[[id]] <- data.table(
    scenario = id, parameter = a$parameter,
    hess_ratio = a$se_ratio,
    bhhh_ratio = bhhh_ratio
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
