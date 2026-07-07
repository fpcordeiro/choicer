#!/usr/bin/env Rscript
# Bitwise sanity gate: run a tiny fixed-seed config through (A) the
# CURRENTLY INSTALLED choicer (non-instrumented, default library path) and
# (B) the HB_PROFILE-instrumented copy (lib.loc = --lib-dir), in TWO
# separate subprocesses so no session ever has two versions of the choicer
# shared object loaded at once. Both processes consume the IDENTICAL
# prepared inputs (X, Z, M, choice_pos, alt_of_row, Ti, priors, inits) --
# generated once, here, with the currently installed package -- so the only
# difference between the two kernel calls is the -DHB_PROFILE compile flag.
#
# Asserts bdraw / wdraw / deltadraw / thetadraw / sigma_d2draw (+ sigma2draw
# for HMNP) are IDENTICAL (bit-for-bit: `identical()`, backed up by a
# max-abs-diff == 0 check) between (A) and (B), for both hmnl_gibbs and
# hmnp_gibbs. Also measures timer overhead: total wall time with vs without
# HB_PROFILE at a moderately larger R (same tiny N/J), both same-process
# comparable since instrumented-vs-not is still two subprocesses.
#
# Usage:
#   Rscript _benchmarks/hb_profiling/01_bitwise_gate.R \
#     --lib-dir <instrumented templib> --ref-lib-dir <non-instrumented templib>
#
# NOTE on Build A: the spec allows "the currently installed (or a
# non-instrumented temp) build" for Build A. The system-installed choicer at
# this repo's default library path predates the HMNL/HMNP feature (stale --
# `simulate_hmnl_data` etc. do not exist there), so this script uses a
# non-instrumented TEMP build compiled from the same source tree instead
# (see 00_build_instrumented.R's sibling: a plain `R CMD INSTALL`, no patch
# applied, into --ref-lib-dir).

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) == 0 || i == length(args)) return(default)
  args[i + 1]
}
this_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(this_file) == 1) dirname(this_file) else "."
lib_dir <- get_arg("--lib-dir",
  "/Users/fernando/.claude/plugins/data/statsclaw-statsclaw/workspace/choicer/tmp/hb_prof_build/templib")
ref_lib_dir <- get_arg("--ref-lib-dir",
  "/Users/fernando/.claude/plugins/data/statsclaw-statsclaw/workspace/choicer/tmp/hb_prof_build/templib_ref")
tmp_dir <- file.path(script_dir, "results", "tmp")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

library(choicer, lib.loc = ref_lib_dir)   # non-instrumented reference build
cat("Build A (non-instrumented): choicer ", as.character(packageVersion("choicer")),
   " @ ", find.package("choicer"), "\n", sep = "")

set_num_threads(5L)

# --- Tiny fixed-seed config (per spec): N=200, J=10, K=4, P=3, R=200,
# burn=100 ---------------------------------------------------------------
N <- 200L; J <- 10L; K <- 4L; P <- 3L; R <- 200L; burn <- 100L; SEED <- 42L

run_gate_for_model <- function(model) {
  cat("\n--- Gate:", model, "---\n")
  beta_true <- rep(0.4, K); theta_true <- rep(0.2, P)
  x_cols <- paste0("x", seq_len(K)); z_cols <- paste0("z", seq_len(P - 1))

  sim <- if (model == "hmnl") {
    simulate_hmnl_data(N = N, T = 1, J = J, beta = beta_true,
                       theta = theta_true, seed = SEED)
  } else {
    simulate_hmnp_data(N = N, T = 1, J = J, beta = beta_true,
                       theta = theta_true, seed = SEED)
  }
  prep_fun <- if (model == "hmnl") prepare_hmnl_data else prepare_hmnp_data
  d <- prep_fun(sim$data, "task", "alt", "choice", x_cols, person_col = NULL,
               alt_covariate_cols = z_cols)

  b_bar <- rep(0, K); A <- 0.01 * diag(K); nu <- K + 3; V <- nu * diag(K)
  theta_bar <- rep(0, d$P); A_theta <- 0.01 * diag(d$P)
  sd_prior <- list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3)

  inputs <- list(model = model, X = d$X, Z = d$Z, M = d$M,
                choice_pos = d$choice_pos, alt_of_row = d$alt_of_row,
                Ti = d$Ti, b_bar = b_bar, A = A, nu = nu, V = V,
                theta_bar = theta_bar, A_theta = A_theta, sd_prior = sd_prior,
                R = R, burn = burn, thin = 1L, seed = SEED + 1000,
                keep_beta_i = 0L)
  if (model == "hmnl") {
    inputs$rc_dist <- as.integer(d$rc_dist)
    inputs$beta_pooled <- rep(0, K)
    inputs$delta_init <- rep(0, d$J)
    inputs$theta_init <- rep(0, d$P)
    inputs$s_init <- 2.38 / sqrt(K)
    inputs$accept_target <- 0.234
  } else {
    inputs$delta_init <- rep(0, d$J)
    inputs$theta_init <- rep(0, d$P)
    inputs$a0 <- 3; inputs$s0 <- 3
  }
  inputs_path <- file.path(tmp_dir, paste0("gate_inputs_", model, ".rds"))
  saveRDS(inputs, inputs_path)

  # --- (A) run here, non-instrumented -----------------------------------
  outA <- if (model == "hmnl") {
    hmnl_gibbs(X = inputs$X, Z = inputs$Z, M = inputs$M,
              choice_pos = inputs$choice_pos, include_outside_option = TRUE,
              alt_of_row = inputs$alt_of_row, Ti = inputs$Ti,
              rc_dist = inputs$rc_dist, beta_pooled = inputs$beta_pooled,
              delta_init = inputs$delta_init, theta_init = inputs$theta_init,
              b_bar = inputs$b_bar, A = inputs$A, nu = inputs$nu,
              V = inputs$V, theta_bar = inputs$theta_bar,
              A_theta = inputs$A_theta, sd_prior = inputs$sd_prior,
              R = inputs$R, burn = inputs$burn, thin = inputs$thin,
              seed = inputs$seed, keep_beta_i = inputs$keep_beta_i,
              s_init = inputs$s_init, accept_target = inputs$accept_target)
  } else {
    hmnp_gibbs(X = inputs$X, Z = inputs$Z, M = inputs$M,
              choice_pos = inputs$choice_pos, include_outside_option = TRUE,
              alt_of_row = inputs$alt_of_row, Ti = inputs$Ti,
              delta_init = inputs$delta_init, theta_init = inputs$theta_init,
              b_bar = inputs$b_bar, A = inputs$A, nu = inputs$nu,
              V = inputs$V, theta_bar = inputs$theta_bar,
              A_theta = inputs$A_theta, sd_prior = inputs$sd_prior,
              a0 = inputs$a0, s0 = inputs$s0, R = inputs$R,
              burn = inputs$burn, thin = inputs$thin, seed = inputs$seed,
              keep_beta_i = inputs$keep_beta_i)
  }
  stopifnot(is.null(outA$profile))   # non-instrumented build: no profile field
  outA_path <- file.path(tmp_dir, paste0("gate_outA_", model, ".rds"))
  saveRDS(outA, outA_path)

  # --- (B) run in a subprocess against the instrumented build -----------
  worker_script <- tempfile(fileext = ".R")
  outB_path <- file.path(tmp_dir, paste0("gate_outB_", model, ".rds"))
  writeLines(sprintf('
    library(choicer, lib.loc = "%s")
    choicer::set_num_threads(5L)
    inputs <- readRDS("%s")
    if (inputs$model == "hmnl") {
      out <- hmnl_gibbs(X = inputs$X, Z = inputs$Z, M = inputs$M,
                choice_pos = inputs$choice_pos, include_outside_option = TRUE,
                alt_of_row = inputs$alt_of_row, Ti = inputs$Ti,
                rc_dist = inputs$rc_dist, beta_pooled = inputs$beta_pooled,
                delta_init = inputs$delta_init, theta_init = inputs$theta_init,
                b_bar = inputs$b_bar, A = inputs$A, nu = inputs$nu,
                V = inputs$V, theta_bar = inputs$theta_bar,
                A_theta = inputs$A_theta, sd_prior = inputs$sd_prior,
                R = inputs$R, burn = inputs$burn, thin = inputs$thin,
                seed = inputs$seed, keep_beta_i = inputs$keep_beta_i,
                s_init = inputs$s_init, accept_target = inputs$accept_target)
    } else {
      out <- hmnp_gibbs(X = inputs$X, Z = inputs$Z, M = inputs$M,
                choice_pos = inputs$choice_pos, include_outside_option = TRUE,
                alt_of_row = inputs$alt_of_row, Ti = inputs$Ti,
                delta_init = inputs$delta_init, theta_init = inputs$theta_init,
                b_bar = inputs$b_bar, A = inputs$A, nu = inputs$nu,
                V = inputs$V, theta_bar = inputs$theta_bar,
                A_theta = inputs$A_theta, sd_prior = inputs$sd_prior,
                a0 = inputs$a0, s0 = inputs$s0, R = inputs$R,
                burn = inputs$burn, thin = inputs$thin, seed = inputs$seed,
                keep_beta_i = inputs$keep_beta_i)
    }
    stopifnot(!is.null(out$profile))   # instrumented build: profile MUST be present
    saveRDS(out, "%s")
  ', lib_dir, inputs_path, outB_path), worker_script)
  status <- system2("Rscript", shQuote(worker_script))
  stopifnot(status == 0)
  outB <- readRDS(outB_path)

  # --- Compare ------------------------------------------------------------
  fields <- c("bdraw", "wdraw", "deltadraw", "thetadraw", "sigma_d2draw")
  if (model == "hmnp") fields <- c(fields, "sigma2draw")
  cmp <- lapply(fields, function(f) {
    a <- outA[[f]]; b <- outB[[f]]
    list(field = f, identical = identical(a, b),
        max_abs_diff = suppressWarnings(max(abs(a - b))))
  })
  cmp_dt <- rbindlist(lapply(cmp, function(x) {
    data.table(field = x$field, identical = x$identical,
              max_abs_diff = x$max_abs_diff)
  }))
  cat("Comparison (A = non-instrumented vs B = instrumented):\n")
  print(cmp_dt)
  all_identical <- all(cmp_dt$identical) && all(cmp_dt$max_abs_diff == 0)
  cat("ALL FIELDS BITWISE IDENTICAL: ", all_identical, "\n")

  list(model = model, all_identical = all_identical, comparison = cmp_dt,
      profile_B = outB$profile)
}

gate_hmnl <- run_gate_for_model("hmnl")
gate_hmnp <- run_gate_for_model("hmnp")

cat("\n=== Timer overhead check (R=2000, same tiny N/J, wall time only) ===\n")
overhead_check <- function(model) {
  R_over <- 2000L
  x_cols <- paste0("x", seq_len(K)); z_cols <- paste0("z", seq_len(P - 1))
  sim <- if (model == "hmnl") {
    simulate_hmnl_data(N = N, T = 1, J = J, beta = rep(0.4, K),
                       theta = rep(0.2, P), seed = SEED)
  } else {
    simulate_hmnp_data(N = N, T = 1, J = J, beta = rep(0.4, K),
                       theta = rep(0.2, P), seed = SEED)
  }
  prep_fun <- if (model == "hmnl") prepare_hmnl_data else prepare_hmnp_data
  d <- prep_fun(sim$data, "task", "alt", "choice", x_cols, person_col = NULL,
               alt_covariate_cols = z_cols)
  b_bar <- rep(0, K); A <- 0.01 * diag(K); nu <- K + 3; V <- nu * diag(K)
  theta_bar <- rep(0, d$P); A_theta <- 0.01 * diag(d$P)
  sd_prior <- list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3)

  t_noninstr <- system.time({
    if (model == "hmnl") {
      hmnl_gibbs(X = d$X, Z = d$Z, M = d$M, choice_pos = d$choice_pos,
                include_outside_option = TRUE, alt_of_row = d$alt_of_row,
                Ti = d$Ti, rc_dist = as.integer(d$rc_dist),
                beta_pooled = rep(0, K), delta_init = rep(0, d$J),
                theta_init = rep(0, d$P), b_bar = b_bar, A = A, nu = nu,
                V = V, theta_bar = theta_bar, A_theta = A_theta,
                sd_prior = sd_prior, R = R_over, burn = R_over %/% 2,
                thin = 1L, seed = SEED + 2000, keep_beta_i = 0L,
                s_init = 2.38 / sqrt(K), accept_target = 0.234)
    } else {
      hmnp_gibbs(X = d$X, Z = d$Z, M = d$M, choice_pos = d$choice_pos,
                include_outside_option = TRUE, alt_of_row = d$alt_of_row,
                Ti = d$Ti, delta_init = rep(0, d$J), theta_init = rep(0, d$P),
                b_bar = b_bar, A = A, nu = nu, V = V, theta_bar = theta_bar,
                A_theta = A_theta, sd_prior = sd_prior, a0 = 3, s0 = 3,
                R = R_over, burn = R_over %/% 2, thin = 1L,
                seed = SEED + 2000, keep_beta_i = 0L)
    }
  })["elapsed"]

  saveRDS(list(X = d$X, Z = d$Z, M = d$M, choice_pos = d$choice_pos,
              alt_of_row = d$alt_of_row, Ti = d$Ti, rc_dist = as.integer(d$rc_dist),
              b_bar = b_bar, A = A, nu = nu, V = V, theta_bar = theta_bar,
              A_theta = A_theta, sd_prior = sd_prior, J = d$J, P = d$P,
              model = model, R = R_over),
         file.path(tmp_dir, paste0("overhead_inputs_", model, ".rds")))

  worker_script <- tempfile(fileext = ".R")
  writeLines(sprintf('
    library(choicer, lib.loc = "%s")
    choicer::set_num_threads(5L)
    inputs <- readRDS("%s")
    K <- ncol(inputs$X)
    t <- system.time({
      if (inputs$model == "hmnl") {
        out <- hmnl_gibbs(X = inputs$X, Z = inputs$Z, M = inputs$M,
                  choice_pos = inputs$choice_pos, include_outside_option = TRUE,
                  alt_of_row = inputs$alt_of_row, Ti = inputs$Ti,
                  rc_dist = inputs$rc_dist, beta_pooled = rep(0, K),
                  delta_init = rep(0, inputs$J), theta_init = rep(0, inputs$P),
                  b_bar = inputs$b_bar, A = inputs$A, nu = inputs$nu,
                  V = inputs$V, theta_bar = inputs$theta_bar,
                  A_theta = inputs$A_theta, sd_prior = inputs$sd_prior,
                  R = inputs$R, burn = floor(inputs$R / 2), thin = 1L,
                  seed = %d, keep_beta_i = 0L, s_init = 2.38 / sqrt(K),
                  accept_target = 0.234)
      } else {
        out <- hmnp_gibbs(X = inputs$X, Z = inputs$Z, M = inputs$M,
                  choice_pos = inputs$choice_pos, include_outside_option = TRUE,
                  alt_of_row = inputs$alt_of_row, Ti = inputs$Ti,
                  delta_init = rep(0, inputs$J), theta_init = rep(0, inputs$P),
                  b_bar = inputs$b_bar, A = inputs$A, nu = inputs$nu,
                  V = inputs$V, theta_bar = inputs$theta_bar,
                  A_theta = inputs$A_theta, sd_prior = inputs$sd_prior,
                  a0 = 3, s0 = 3, R = inputs$R, burn = floor(inputs$R / 2),
                  thin = 1L, seed = %d, keep_beta_i = 0L)
      }
    })["elapsed"]
    cat("INSTRUMENTED_WALL_TIME=", t, "\\n", sep = "")
    cat("PROFILE_T_TOTAL=", out$profile$t_total, "\\n", sep = "")
  ', lib_dir, file.path(tmp_dir, paste0("overhead_inputs_", model, ".rds")),
    SEED + 2000, SEED + 2000), worker_script)
  out_lines <- system2("Rscript", shQuote(worker_script), stdout = TRUE)
  wall_instr <- as.numeric(sub("INSTRUMENTED_WALL_TIME=", "",
                               grep("INSTRUMENTED_WALL_TIME=", out_lines, value = TRUE)))
  prof_total <- as.numeric(sub("PROFILE_T_TOTAL=", "",
                               grep("PROFILE_T_TOTAL=", out_lines, value = TRUE)))

  data.table(model = model, R = R_over,
            wall_noninstrumented_sec = unname(t_noninstr),
            wall_instrumented_sec = wall_instr,
            profile_t_total_sec = prof_total,
            overhead_pct = 100 * (wall_instr - unname(t_noninstr)) / unname(t_noninstr))
}

overhead_dt <- rbindlist(lapply(c("hmnl", "hmnp"), overhead_check))
print(overhead_dt)

cat("\n=== GATE SUMMARY ===\n")
cat("HMNL bitwise identical: ", gate_hmnl$all_identical, "\n")
cat("HMNP bitwise identical: ", gate_hmnp$all_identical, "\n")

saveRDS(list(hmnl = gate_hmnl, hmnp = gate_hmnp, overhead = overhead_dt),
       file.path(script_dir, "results", "bitwise_gate_result.rds"))
fwrite(overhead_dt, file.path(script_dir, "results", "timer_overhead.csv"))
