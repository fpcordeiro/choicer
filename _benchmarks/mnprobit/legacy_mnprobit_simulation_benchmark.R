# MNP cross-package benchmark
# Compares choicer::run_mnprobit() against two reference implementations of
# the Bayesian multinomial probit on the simulate_mnp_data() DGP used by
# inst/simulations/mnp_simulation.R:
#   - bayesm::rmnpGibbs(): the same McCulloch-Rossi (1994) Gibbs sampler,
#     run on the identical differenced design matrix and priors as choicer
#     (built by prepare_mnp_data()) -- the closest apples-to-apples reference.
#   - MNP::mnp(): the Imai-van Dyk (2005) marginal data augmentation sampler,
#     with its own design construction and identification scheme; priors are
#     matched where the parameterizations align.
#
# One estimation run per package on a common dataset and common chain length
# (this is a smoke-level comparison, not a Monte Carlo study). All estimates
# are posterior means of the identified draws (normalized per draw by
# sigma_11), so the packages and `sim$true_params` are directly comparable.
# `se` is the posterior SD and `ess` the effective sample size of the kept
# draws -- for MCMC, time per effective draw is the meaningful speed metric,
# not wall time alone.
#
# Reference packages are not in DESCRIPTION; install them manually first:
#   install.packages(c("bayesm", "MNP"))
# Run from package root: Rscript _benchmarks/mnprobit_simulation_benchmark.R

rm(list = ls(all.names = TRUE))
gc()

# devtools::load_all() lets pkgbuild inject debug flags (-O0 -UNDEBUG) by
# overriding variables in the user Makevars, which would benchmark an
# unoptimized choicer against -O2 CRAN builds of the reference packages (and
# can clobber a custom OpenMP setup in ~/.R/Makevars on macOS). Disable the
# injection and force a clean rebuild so the comparison is fair (pkgbuild
# does not invalidate its cache on flag changes, hence the unconditional
# clean).
options(pkg.build_extra_flags = FALSE)
pkgbuild::clean_dll()
devtools::load_all()
source("_benchmarks/bench_helpers.R")

# choicer parallelizes the data-augmentation step with OpenMP (draws are
# bitwise invariant to thread count). bayesm and MNP are single-threaded, so
# the default here is 1 thread for a fair wall-time comparison; raise it to
# see the multithreaded speed.
N_THREADS <- 1L
set_num_threads(N_THREADS)

# Simulated DGP (same preset as inst/simulations/mnp_simulation.R) -------------
sim <- simulate_mnp_data(N = 5000, J = 3, seed = 123)
print(sim)

# Common chain settings ---------------------------------------------------------
mcmc <- list(R = 20000, burn = 5000, thin = 5)
cat(sprintf(
  "\nChain: R = %d, burn = %d, thin = %d (%d kept draws); choicer threads = %d\n",
  mcmc$R, mcmc$burn, mcmc$thin, (mcmc$R - mcmc$burn) %/% mcmc$thin, N_THREADS
))

# Run benchmarks ----------------------------------------------------------------
set.seed(42)
res <- benchmark_fit(sim, packages = c("choicer", "bayesm", "MNP"), mcmc = mcmc)
print(res)

# Side-by-side: truth vs posterior mean / SD ------------------------------------
p <- sim$settings$J - 1L
truth <- data.table(
  parameter = c(paste0("x", seq_len(sim$settings$K_x)),
                paste0("ASC_", setdiff(seq_len(sim$settings$J), sim$settings$base_alt)),
                .sigma_names(p)),
  truth = c(sim$true_params$beta, sim$true_params$delta,
            unlist(lapply(seq_len(p), function(i) sim$true_params$Sigma[i, seq_len(i)])))
)

cmp <- dcast(res, parameter ~ package, value.var = c("estimate", "se"))
cmp <- cmp[match(truth$parameter, parameter)]
cmp <- cbind(truth, cmp[, -"parameter"])
num_cols <- setdiff(names(cmp), "parameter")
cmp[, (num_cols) := lapply(.SD, round, 3), .SDcols = num_cols]
cat("\n--- Posterior means (identified scale) vs truth ---\n")
print(cmp)

# Sampling efficiency -----------------------------------------------------------
# min ESS is the binding constraint (Sigma_11 is the normalization and has no
# ESS). sec_per_1k_ess = wall time to obtain 1000 effective draws of the
# slowest-mixing parameter.
eff <- res[, .(
  time_sec = round(time_sec[1L], 1),
  min_ess  = round(min(ess, na.rm = TRUE)),
  med_ess  = round(stats::median(ess, na.rm = TRUE))
), by = package]
# grouped [.data.table keeps the choicer_benchmark class; drop it before any
# further [.data.table call, or the schema-specific print method is
# dispatched on this summary (:= auto-prints through custom print methods)
setattr(eff, "class", class(eff)[!class(eff) %in% "choicer_benchmark"])
eff[, sec_per_1k_ess := round(1000 * time_sec / min_ess, 1)]
cat("\n--- Sampling efficiency ---\n")
print(eff)
