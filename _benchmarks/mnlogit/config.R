mnlogit_benchmark_config <- function() {
  list(
    n_runs = 3L,
    n_vec = c(1000L, 5000L, 10000L),
    j_vec = c(5L, 10L, 25L),
    fixed_n = 10000L,
    fixed_j = 10L,
    grid_mode = "two_sweeps",
    packages = c("choicer", "mlogit", "logitr", "gmnl", "apollo", "mixl"),
    skip_missing = TRUE,
    coef_abs_tol = 1e-3,
    coef_rel_tol = 1e-3,
    loglik_rel_tol = 1e-6,
    n_threads = 1L,
    seed = 12345L,
    tag = NULL,
    beta = c(0.8, -0.6),
    optimizer = "nloptr",
    optimizer_control = list(
      algorithm = "NLOPT_LD_LBFGS",
      xtol_rel = 1e-8,
      maxeval = 1000L,
      print_level = 0L
    ),
    plot_metric = "mean_fit_time_sec",
    clean_dll = TRUE
  )
}

mnlogit_normalize_config <- function(config) {
  config$n_runs <- as.integer(config$n_runs[[1L]])
  config$n_vec <- as.integer(config$n_vec)
  config$j_vec <- as.integer(config$j_vec)
  config$fixed_n <- as.integer(config$fixed_n[[1L]])
  config$fixed_j <- as.integer(config$fixed_j[[1L]])
  config$n_threads <- as.integer(config$n_threads[[1L]])
  config$seed <- as.integer(config$seed[[1L]])
  config$packages <- unique(tolower(config$packages))
  config$packages[config$packages == "survival"] <- "clogit"
  config$grid_mode <- config$grid_mode[[1L]]
  config$plot_metric <- config$plot_metric[[1L]]
  config
}
