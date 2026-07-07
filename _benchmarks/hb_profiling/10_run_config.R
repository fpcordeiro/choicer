#!/usr/bin/env Rscript
# Worker: run ONE profiling config against the HB_PROFILE-instrumented build
# and write a one-row result (JSON) to --out. Meant to be launched as a fresh
# subprocess (so `/usr/bin/time -l` on the parent call captures this
# process's peak RSS) -- never call this from inside a long-lived session.
#
# The kernels (hmnl_gibbs / hmnp_gibbs) are called DIRECTLY rather than via
# run_hmnlogit()/run_hmnprobit(): the R wrappers do not surface the `profile`
# list element (adding that would mean touching package R source, which is
# out of the simulator's write surface). This worker reproduces the
# wrappers' prior defaults and initialization logic (pooled MNL MLE for HMNL,
# shrunk log-share delta init, OLS theta init) so the profiled run matches
# real usage; init cost is timed and reported SEPARATELY from the Gibbs
# kernel's own `profile` block totals.
#
# Args (key=value, any order):
#   model=hmnl|hmnp|mnprobit   N= J= K= P= R= burn= thin=1 threads=5
#   keep_beta_i=none|means|draws   seed=1   lib_dir=<instrumented lib path>
#   config_id=<label>   out=<path to write JSON result>
#   trace=0

suppressPackageStartupMessages({
  library(data.table)
})

args_raw <- commandArgs(trailingOnly = TRUE)
kv <- strsplit(args_raw, "=", fixed = TRUE)
opt <- setNames(vapply(kv, `[`, character(1), 2), vapply(kv, `[`, character(1), 1))

geti <- function(name, default = NULL) {
  if (is.na(opt[name]) || is.null(opt[[name]])) return(default)
  as.integer(opt[[name]])
}
getc <- function(name, default = NULL) {
  v <- opt[[name]]
  if (is.null(v) || is.na(v)) default else v
}

model      <- getc("model", "hmnl")
N          <- geti("N", 500L)
J          <- geti("J", 20L)
K          <- geti("K", 10L)
P          <- geti("P", 10L)
R          <- geti("R", 100L)
burn       <- geti("burn", R %/% 3L)
thin       <- geti("thin", 1L)
threads    <- geti("threads", 5L)
keep_beta_i <- getc("keep_beta_i", "means")
seed       <- geti("seed", 1L)
lib_dir    <- getc("lib_dir", NULL)
config_id  <- getc("config_id", "unnamed")
out_path   <- getc("out", NULL)
trace      <- geti("trace", 0L)

stopifnot(!is.null(lib_dir), !is.null(out_path))
stopifnot(model %in% c("hmnl", "hmnp", "mnprobit"))

library(choicer, lib.loc = lib_dir)
choicer::set_num_threads(threads)

t_setup_start <- Sys.time()

result <- list(
  config_id = config_id, model = model, N = N, J = J, K = K, P = P,
  R = R, burn = burn, thin = thin, threads = threads,
  keep_beta_i = keep_beta_i, seed = seed,
  choicer_version = as.character(utils::packageVersion("choicer")),
  choicer_lib = find.package("choicer")
)

keep_code <- switch(keep_beta_i, none = 0L, means = 1L, draws = 2L)

if (model %in% c("hmnl", "hmnp")) {
  beta_true  <- rep(0.4, K)
  theta_true <- rep(0.2, P)   # length P: intercept + (P-1) z-covariates
  x_cols <- paste0("x", seq_len(K))
  z_cols <- if (P > 1) paste0("z", seq_len(P - 1)) else character(0)

  sim <- if (model == "hmnl") {
    simulate_hmnl_data(N = N, T = 1, J = J, beta = beta_true,
                       theta = theta_true, seed = seed)
  } else {
    simulate_hmnp_data(N = N, T = 1, J = J, beta = beta_true,
                       theta = theta_true, seed = seed)
  }

  prep_fun <- if (model == "hmnl") prepare_hmnl_data else prepare_hmnp_data
  d <- prep_fun(sim$data, "task", "alt", "choice", x_cols,
               person_col = NULL,
               alt_covariate_cols = if (length(z_cols)) z_cols else NULL)

  b_bar <- rep(0, K); A <- 0.01 * diag(K)
  nu <- K + 3; V <- nu * diag(K)
  theta_bar <- rep(0, d$P); A_theta <- 0.01 * diag(d$P)
  sd_prior <- list(half_cauchy = TRUE, s_d = 1, c0 = 3, d0 = 3)

  if (model == "hmnl") {
    rc <- as.integer(d$rc_dist)
    weights <- rep(1, d$n_tasks)
    eval_f <- function(theta_par) {
      mnl_loglik_gradient_parallel(
        theta_par, d$X, d$alt_idx, d$choice_pos, d$M,
        weights, use_asc = FALSE, include_outside_option = TRUE
      )
    }
    opt_r <- tryCatch(run_optimizer("nloptr", rep(0, K), eval_f),
                      error = function(e) NULL)
    beta_pooled <- if (!is.null(opt_r) && all(is.finite(opt_r$par))) {
      opt_r$par
    } else rep(0, K)

    am <- d$alt_mapping
    n_in <- am[am$alt_int > 0, ][["N_CHOICES"]]
    n_out <- am[am$alt_int == 0, ][["N_CHOICES"]]
    delta_init <- log((n_in + 0.5) / (n_out + 0.5))
    theta_init <- as.numeric(qr.solve(d$Z, delta_init))

    t_setup_end <- Sys.time()

    out <- hmnl_gibbs(
      X = d$X, Z = d$Z, M = d$M, choice_pos = d$choice_pos,
      include_outside_option = TRUE, alt_of_row = d$alt_of_row, Ti = d$Ti,
      rc_dist = rc, beta_pooled = beta_pooled,
      delta_init = delta_init, theta_init = theta_init,
      b_bar = b_bar, A = A, nu = nu, V = V,
      theta_bar = theta_bar, A_theta = A_theta, sd_prior = sd_prior,
      R = R, burn = burn, thin = thin, seed = seed + 1000,
      keep_beta_i = keep_code, s_init = 2.38 / sqrt(K),
      accept_target = 0.234, trace = trace
    )
  } else {
    delta_init <- rep(0, d$J); theta_init <- rep(0, d$P)
    t_setup_end <- Sys.time()
    out <- hmnp_gibbs(
      X = d$X, Z = d$Z, M = d$M, choice_pos = d$choice_pos,
      include_outside_option = TRUE, alt_of_row = d$alt_of_row, Ti = d$Ti,
      delta_init = delta_init, theta_init = theta_init,
      b_bar = b_bar, A = A, nu = nu, V = V,
      theta_bar = theta_bar, A_theta = A_theta, sd_prior = sd_prior,
      a0 = 3, s0 = 3, R = R, burn = burn, thin = thin, seed = seed + 1000,
      keep_beta_i = keep_code, trace = trace
    )
  }

  result$t_setup_sec <- as.numeric(difftime(t_setup_end, t_setup_start,
                                            units = "secs"))
  result$profile <- out$profile
  result$R_keep <- out$R_keep
} else {
  # run_mnprobit contrast: wall time only, NOT instrumented (the currently
  # installed build is fine here -- see spec: "no instrumentation needed").
  beta_true <- rep(0.4, K)
  sim <- simulate_mnp_data(N = N, J = J, beta = beta_true, seed = seed)
  t0 <- Sys.time()
  fit <- run_mnprobit(sim$data, "task", "alt", "choice",
                      paste0("x", seq_len(K)),
                      mcmc = list(R = R, burn = burn, seed = seed + 1000))
  t1 <- Sys.time()
  result$t_setup_sec <- 0
  result$profile <- list(t_total = as.numeric(difftime(t1, t0, units = "secs")))
  result$R_keep <- fit$mcmc$R - fit$mcmc$burn
}

jsonlite_ok <- requireNamespace("jsonlite", quietly = TRUE)
if (jsonlite_ok) {
  writeLines(jsonlite::toJSON(result, auto_unbox = TRUE, digits = 10), out_path)
} else {
  saveRDS(result, sub("\\.json$", ".rds", out_path))
}
message("worker done: ", config_id)
