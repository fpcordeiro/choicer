# _benchmarks/halton_benchmark.R
#
# Benchmarks for the on-the-fly Halton draw generator (Feature 1).
# Three benchmark groups:
#
#   1. Generation speed: HaltonGen (C++ on-the-fly, scramble="none") vs
#      get_halton_normals() (randtoolbox) for various (S, N, K_w) configs.
#
#   2. Per-gradient overhead: mxl_loglik_gradient_parallel store vs generate,
#      at small m_i (J=3) and large m_i (J=50), showing overhead is diluted
#      by large m_i.
#
#   3. End-to-end fit wall-clock + peak RSS: run_mxlogit store vs generate on
#      a large-N config (N=2000, J=5, K_w=3, S=100).
#
# Not part of the installed package (_benchmarks/ is .Rbuildignore'd).
# Run from package root:
#   Rscript _benchmarks/halton_benchmark.R
#
# Requires: choicer installed or loaded via devtools::load_all()

rm(list = ls(all.names = TRUE))
gc()

options(pkg.build_extra_flags = FALSE)
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})
source("_benchmarks/bench_helpers.R")

library(data.table)

cat("\n====================================================================\n")
cat("  Halton On-the-Fly Generator Benchmark\n")
cat("  choicer version:", as.character(packageVersion("choicer")), "\n")
cat("====================================================================\n\n")

# ---------------------------------------------------------------------------
# Helper: format time nicely
# ---------------------------------------------------------------------------
.fmt_ms <- function(t_sec) sprintf("%.2f ms", t_sec * 1000)
.fmt_us <- function(t_sec) sprintf("%.2f us", t_sec * 1e6)

# ---------------------------------------------------------------------------
# 1. Generation speed: HaltonGen vs get_halton_normals
# ---------------------------------------------------------------------------

cat("----------------------------------------------------------------------\n")
cat("1. GENERATION SPEED: HaltonGen (C++, scramble=none) vs get_halton_normals\n")
cat("----------------------------------------------------------------------\n")

gen_speed_configs <- list(
  list(S = 100L, N = 500L,  K_w = 3L),
  list(S = 100L, N = 2000L, K_w = 3L),
  list(S = 200L, N = 2000L, K_w = 5L),
  list(S = 100L, N = 5000L, K_w = 3L)
)

for (cfg in gen_speed_configs) {
  S <- cfg$S; N <- cfg$N; K_w <- cfg$K_w

  # Warm up
  invisible(get_halton_normals(S, N, K_w))
  invisible(choicer:::halton_generate_normal(S, N, K_w, seed = 0L, scramble = 0L))

  # Repeated timing
  R_reps <- 10L

  t_store <- system.time(
    for (r in seq_len(R_reps)) get_halton_normals(S, N, K_w)
  )["elapsed"] / R_reps

  t_gen <- system.time(
    for (r in seq_len(R_reps))
      choicer:::halton_generate_normal(S, N, K_w, seed = 0L, scramble = 0L)
  )["elapsed"] / R_reps

  total_draws <- S * N * K_w
  cat(sprintf(
    "  S=%d, N=%d, K_w=%d (%d total draws):\n",
    S, N, K_w, total_draws
  ))
  cat(sprintf(
    "    randtoolbox get_halton_normals : %s\n", .fmt_ms(t_store)
  ))
  cat(sprintf(
    "    HaltonGen C++ (scramble=none)  : %s\n", .fmt_ms(t_gen)
  ))
  cat(sprintf(
    "    Speedup (store/gen): %.2fx\n", t_store / t_gen
  ))
  cat("\n")
}

# Owen scramble overhead vs no scramble
cat("  Owen scramble overhead (S=100, N=2000, K_w=5):\n")
S <- 100L; N <- 2000L; K_w <- 5L; R_reps <- 10L
t_none <- system.time(
  for (r in seq_len(R_reps))
    choicer:::halton_generate_normal(S, N, K_w, seed = 42L, scramble = 0L)
)["elapsed"] / R_reps
t_owen <- system.time(
  for (r in seq_len(R_reps))
    choicer:::halton_generate_normal(S, N, K_w, seed = 42L, scramble = 1L)
)["elapsed"] / R_reps
cat(sprintf("    scramble=none : %s\n", .fmt_ms(t_none)))
cat(sprintf("    scramble=owen : %s\n", .fmt_ms(t_owen)))
cat(sprintf("    Owen overhead : %.1fx\n", t_owen / t_none))
cat("\n")

# ---------------------------------------------------------------------------
# 2. Per-gradient overhead: store vs generate, small vs large m_i
# ---------------------------------------------------------------------------

cat("----------------------------------------------------------------------\n")
cat("2. PER-GRADIENT OVERHEAD: store vs generate (scramble=none)\n")
cat("   Small m_i (J=3) vs large m_i (J=50)\n")
cat("----------------------------------------------------------------------\n")

.bench_gradient <- function(J, N, S, K_w, R_reps = 20L) {
  # Simulate data
  set.seed(1)
  sim <- simulate_mxl_data(
    N = N, J = J,
    beta  = rep(0.5, 1),
    delta = rep(0, J),
    Sigma = matrix(0.5, K_w, K_w) + 0.5 * diag(K_w),
    seed  = 1,
    outside_option = FALSE,
    vary_choice_set = FALSE
  )

  K_x <- 1L
  inp <- prepare_mxl_data(
    sim$data, "id", "alt", "choice", "x1",
    paste0("w", seq_len(K_w)),
    rc_correlation = FALSE
  )
  L_size <- K_w
  n_asc  <- J - 1L
  theta  <- c(runif(K_x, -0.3, 0.3), log(runif(K_w, 0.3, 0.8)),
              runif(n_asc, -0.2, 0.2))

  eta       <- get_halton_normals(S, inp$N, K_w)
  eta_empty <- array(0, dim = c(K_w, 0L, 0L))

  rc_dist <- rep(0L, K_w)

  # Warm up
  mxl_loglik_gradient_parallel(theta, inp$X, inp$W, inp$alt_idx, inp$choice_idx,
    inp$M, inp$weights, eta, rc_dist, FALSE, FALSE, TRUE, FALSE)
  mxl_loglik_gradient_parallel(theta, inp$X, inp$W, inp$alt_idx, inp$choice_idx,
    inp$M, inp$weights, eta_empty, rc_dist, FALSE, FALSE, TRUE, FALSE,
    gen_seed=0L, gen_scramble=0L, gen_S=S)

  t_store <- system.time(
    for (r in seq_len(R_reps))
      mxl_loglik_gradient_parallel(theta, inp$X, inp$W, inp$alt_idx, inp$choice_idx,
        inp$M, inp$weights, eta, rc_dist, FALSE, FALSE, TRUE, FALSE)
  )["elapsed"] / R_reps

  t_gen <- system.time(
    for (r in seq_len(R_reps))
      mxl_loglik_gradient_parallel(theta, inp$X, inp$W, inp$alt_idx, inp$choice_idx,
        inp$M, inp$weights, eta_empty, rc_dist, FALSE, FALSE, TRUE, FALSE,
        gen_seed=0L, gen_scramble=0L, gen_S=S)
  )["elapsed"] / R_reps

  list(J=J, N=N, S=S, K_w=K_w, t_store=t_store, t_gen=t_gen)
}

for (J in c(3L, 50L)) {
  cfg <- .bench_gradient(J = J, N = 300L, S = 50L, K_w = 2L)
  cat(sprintf(
    "  J=%d, N=%d, S=%d, K_w=%d:\n", cfg$J, cfg$N, cfg$S, cfg$K_w
  ))
  cat(sprintf("    store mode : %s\n", .fmt_ms(cfg$t_store)))
  cat(sprintf("    gen mode   : %s\n", .fmt_ms(cfg$t_gen)))
  overhead_pct <- 100 * (cfg$t_gen - cfg$t_store) / cfg$t_store
  cat(sprintf("    gen overhead: %+.1f%%\n", overhead_pct))
  cat("\n")
}

# ---------------------------------------------------------------------------
# 3. End-to-end fit: store vs generate, large N
# ---------------------------------------------------------------------------

cat("----------------------------------------------------------------------\n")
cat("3. END-TO-END FIT: store vs generate (scramble=none and owen)\n")
cat("   N=2000, J=5, K_w=3, S=100\n")
cat("----------------------------------------------------------------------\n")

set.seed(42)
sim_large <- simulate_mxl_data(
  N   = 2000L, J = 5L,
  beta  = c(0.8),
  delta = rep(0, 5L),
  Sigma = diag(c(0.4, 0.5, 0.3)),
  seed  = 42,
  outside_option = FALSE,
  vary_choice_set = FALSE
)
dt_large <- sim_large$data
S_fit <- 100L

# Measure peak RSS before fits
gc()
rss_before <- gc()["Vcells", 2] * 8 / 1024^2  # MB

# Store mode
tic <- proc.time()["elapsed"]
fit_store <- run_mxlogit(
  data = dt_large, id_col = "id", alt_col = "alt", choice_col = "choice",
  covariate_cols = "x1", random_var_cols = c("w1", "w2", "w3"),
  S = S_fit, rc_correlation = FALSE, rc_mean = FALSE, rc_dist = rep(0L, 3L),
  use_asc = TRUE, include_outside_option = FALSE, draws = "store",
  control = list(print_level = 0L, maxeval = 200L)
)
t_store <- proc.time()["elapsed"] - tic
gc()
rss_store <- gc()["Vcells", 2] * 8 / 1024^2  # MB after GC

# Generate mode (scramble=none, seed=0 for fair comparison)
tic <- proc.time()["elapsed"]
fit_gen_none <- run_mxlogit(
  data = dt_large, id_col = "id", alt_col = "alt", choice_col = "choice",
  covariate_cols = "x1", random_var_cols = c("w1", "w2", "w3"),
  S = S_fit, rc_correlation = FALSE, rc_mean = FALSE, rc_dist = rep(0L, 3L),
  use_asc = TRUE, include_outside_option = FALSE,
  draws = "generate", seed = 0L, scramble = "none",
  control = list(print_level = 0L, maxeval = 200L)
)
t_gen_none <- proc.time()["elapsed"] - tic
gc()
rss_gen <- gc()["Vcells", 2] * 8 / 1024^2

# Generate mode (Owen scramble)
tic <- proc.time()["elapsed"]
fit_gen_owen <- run_mxlogit(
  data = dt_large, id_col = "id", alt_col = "alt", choice_col = "choice",
  covariate_cols = "x1", random_var_cols = c("w1", "w2", "w3"),
  S = S_fit, rc_correlation = FALSE, rc_mean = FALSE, rc_dist = rep(0L, 3L),
  use_asc = TRUE, include_outside_option = FALSE,
  draws = "generate", seed = 42L, scramble = "owen",
  control = list(print_level = 0L, maxeval = 200L)
)
t_gen_owen <- proc.time()["elapsed"] - tic

cat(sprintf("  Store mode (randtoolbox)      : %.2f sec  (loglik=%.4f)\n",
            t_store,     as.numeric(logLik(fit_store))))
cat(sprintf("  Generate mode (scramble=none) : %.2f sec  (loglik=%.4f)\n",
            t_gen_none,  as.numeric(logLik(fit_gen_none))))
cat(sprintf("  Generate mode (scramble=owen) : %.2f sec  (loglik=%.4f)\n",
            t_gen_owen,  as.numeric(logLik(fit_gen_owen))))
cat(sprintf("\n  Peak RSS after store fit: ~%.1f MB  (eta cube: S*N*K_w*8 = %.1f MB)\n",
            rss_store,
            S_fit * 2000 * 3 * 8 / 1024^2))
cat(sprintf("  Peak RSS after gen fit:   ~%.1f MB  (eta cube eliminated)\n", rss_gen))
cat(sprintf("\n  Wall-clock gen/store ratio: %.2fx\n", t_gen_none / t_store))

cat("\n====================================================================\n")
cat("  Benchmark complete.\n")
cat("====================================================================\n")
