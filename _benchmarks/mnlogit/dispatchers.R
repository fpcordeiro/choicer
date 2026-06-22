mnlogit_dependency <- function(package) {
  switch(package,
    choicer = "choicer",
    mlogit = "mlogit",
    logitr = "logitr",
    gmnl = "gmnl",
    apollo = "apollo",
    mixl = "mixl",
    clogit = "survival",
    package
  )
}

mnlogit_raw_row <- function(package, spec, run, status, converged = NA,
                            setup_time_sec = NA_real_, fit_time_sec = NA_real_,
                            total_time_sec = NA_real_, loglik = NA_real_,
                            n_params = NA_integer_, error = NA_character_,
                            warnings = NA_character_, timeout_sec = NA_real_,
                            process_elapsed_sec = NA_real_,
                            worker_exit_status = NA_integer_) {
  dep <- mnlogit_dependency(package)
  data.table::data.table(
    benchmark = "mnlogit",
    sweep = spec$sweep,
    dimension = spec$dimension,
    dimension_value = spec$dimension_value,
    spec_id = spec$spec_id,
    N = as.integer(spec$N),
    J = as.integer(spec$J),
    run = as.integer(run),
    package = package,
    package_dependency = dep,
    package_version = bench_package_version(dep),
    status = status,
    converged = as.logical(converged),
    setup_time_sec = as.numeric(setup_time_sec),
    fit_time_sec = as.numeric(fit_time_sec),
    total_time_sec = as.numeric(total_time_sec),
    loglik = as.numeric(loglik),
    n_params = as.integer(n_params),
    timeout_sec = as.numeric(timeout_sec),
    process_elapsed_sec = as.numeric(process_elapsed_sec),
    worker_exit_status = as.integer(worker_exit_status),
    error = error,
    warnings = warnings
  )
}

mnlogit_canonicalize <- function(package, estimates, x_cols, J, raw_loglik,
                                 spec, run, status = "ok") {
  if (is.null(estimates) || !length(estimates)) return(data.table::data.table())
  raw_names <- names(estimates)
  estimates <- as.numeric(estimates)
  if (is.null(raw_names)) raw_names <- rep("", length(estimates))

  out <- data.table::data.table(
    benchmark = "mnlogit",
    sweep = spec$sweep,
    dimension = spec$dimension,
    dimension_value = spec$dimension_value,
    spec_id = spec$spec_id,
    N = as.integer(spec$N),
    J = as.integer(spec$J),
    run = as.integer(run),
    package = package,
    raw_parameter = raw_names,
    parameter = raw_names,
    group = "other",
    estimate = estimates,
    loglik = as.numeric(raw_loglik),
    status = status
  )

  out[raw_parameter %in% x_cols, `:=`(parameter = raw_parameter, group = "beta")]
  out[grepl("^beta_x[0-9]+$", raw_parameter), `:=`(
    parameter = sub("^beta_", "", raw_parameter),
    group = "beta"
  )]
  out[grepl("^b_x[0-9]+$", raw_parameter), `:=`(
    parameter = sub("^b_", "", raw_parameter),
    group = "beta"
  )]

  out[grepl("^ASC_[0-9]+$", raw_parameter), `:=`(
    parameter = raw_parameter,
    group = "asc"
  )]
  out[grepl("^asc_[0-9]+$", raw_parameter), `:=`(
    parameter = sub("^asc_", "ASC_", raw_parameter),
    group = "asc"
  )]
  out[grepl("^alt_factor[0-9]+$", raw_parameter), `:=`(
    parameter = sub("^alt_factor", "ASC_", raw_parameter),
    group = "asc"
  )]
  out[grepl("^alt_factor[T.]?[0-9]+$", raw_parameter), `:=`(
    parameter = paste0("ASC_", sub("^alt_factor[T.]?", "", raw_parameter)),
    group = "asc"
  )]
  out[grepl("\\(Intercept\\):[0-9]+", raw_parameter), `:=`(
    parameter = paste0("ASC_", sub(".*\\(Intercept\\):", "", raw_parameter)),
    group = "asc"
  )]
  out[grepl("\\(intercept\\):[0-9]+", raw_parameter, ignore.case = TRUE), `:=`(
    parameter = paste0("ASC_", sub(".*\\(intercept\\):", "", raw_parameter, ignore.case = TRUE)),
    group = "asc"
  )]
  out[grepl("^[0-9]+:\\(Intercept\\)", raw_parameter), `:=`(
    parameter = paste0("ASC_", sub(":.*", "", raw_parameter)),
    group = "asc"
  )]
  out[grepl("^[0-9]+:\\(intercept\\)", raw_parameter, ignore.case = TRUE), `:=`(
    parameter = paste0("ASC_", sub(":.*", "", raw_parameter)),
    group = "asc"
  )]

  expected <- c(x_cols, paste0("ASC_", 2:J))
  out <- out[parameter %in% expected]
  out[, parameter_order := match(parameter, expected)]
  data.table::setorder(out, parameter_order)
  out[, parameter_order := NULL]
  out[]
}

mnlogit_choicer_converged <- function(fit) {
  code <- fit$convergence
  if (is.null(code) || length(code) != 1L || !is.finite(code)) return(FALSE)
  nm <- fit$optimizer$name %||% "nloptr"
  switch(nm,
    nloptr = as.integer(code) %in% c(1L, 2L, 3L, 4L),
    optim = as.integer(code) == 0L,
    as.integer(code) %in% c(0L, 1L, 2L, 3L, 4L)
  )
}

mnlogit_setup_choicer <- function(dt, x_cols, J, config) {
  choicer::prepare_mnl_data(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = x_cols,
    include_outside_option = FALSE
  )
}

mnlogit_fit_choicer <- function(setup, x_cols, J, config) {
  fit <- choicer::run_mnlogit(
    input_data = setup,
    optimizer = config$optimizer,
    control = config$optimizer_control,
    use_asc = TRUE,
    keep_data = FALSE,
    scale_vars = "none"
  )
  list(
    estimates = stats::coef(fit),
    loglik = as.numeric(stats::logLik(fit)),
    converged = mnlogit_choicer_converged(fit),
    n_params = length(stats::coef(fit))
  )
}

mnlogit_setup_mlogit <- function(dt, x_cols, J, config) {
  dt <- data.table::copy(dt)
  dt[, alt_factor := factor(as.character(alt), levels = as.character(seq_len(J)))]
  list(
    data = dfidx::dfidx(as.data.frame(dt), idx = c("id", "alt_factor")),
    formula = stats::as.formula(paste("choice ~", paste(x_cols, collapse = " + "), "| 1"))
  )
}

mnlogit_fit_mlogit <- function(setup, x_cols, J, config) {
  fit <- mlogit::mlogit(
    setup$formula,
    data = setup$data,
    reflevel = "1",
    print.level = 0,
    hessian = FALSE
  )
  est <- stats::coef(fit)
  list(
    estimates = est,
    loglik = as.numeric(stats::logLik(fit)),
    converged = TRUE,
    n_params = length(est)
  )
}

mnlogit_setup_logitr <- function(dt, x_cols, J, config) {
  dt <- data.table::copy(dt)
  dt[, alt_factor := factor(as.character(alt), levels = as.character(seq_len(J)))]
  as.data.frame(dt)
}

mnlogit_fit_logitr <- function(setup, x_cols, J, config) {
  old_core_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  if (is.na(parallel::detectCores())) {
    Sys.setenv("_R_CHECK_LIMIT_CORES_" = "true")
  }
  on.exit({
    if (is.na(old_core_limit)) {
      Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
    } else {
      Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_core_limit)
    }
  }, add = TRUE)

  fit <- logitr::logitr(
    data = setup,
    outcome = "choice",
    obsID = "id",
    pars = c(x_cols, "alt_factor"),
    scaleInputs = FALSE,
    vcov = FALSE,
    predict = FALSE,
    numCores = 1L,
    options = config$optimizer_control
  )
  est <- fit$coefficients
  list(
    estimates = est,
    loglik = as.numeric(fit$logLik),
    converged = isTRUE(fit$status > 0),
    n_params = length(est)
  )
}

mnlogit_setup_gmnl <- function(dt, x_cols, J, config) {
  dt <- data.table::copy(dt)
  dt[, alt_factor := factor(as.character(alt), levels = as.character(seq_len(J)))]
  list(
    data = mlogit::mlogit.data(
      as.data.frame(dt),
      choice = "choice",
      shape = "long",
      alt.var = "alt_factor",
      chid.var = "id"
    ),
    formula = stats::as.formula(paste("choice ~", paste(x_cols, collapse = " + "), "| 1"))
  )
}

mnlogit_fit_gmnl <- function(setup, x_cols, J, config) {
  fit <- gmnl::gmnl(
    setup$formula,
    data = setup$data,
    model = "mnl",
    reflevel = "1",
    print.init = FALSE,
    gradient = TRUE
  )
  est <- stats::coef(fit)
  list(
    estimates = est,
    loglik = as.numeric(stats::logLik(fit)),
    converged = TRUE,
    n_params = length(est)
  )
}

mnlogit_setup_clogit <- function(dt, x_cols, J, config) {
  dt <- data.table::copy(dt)
  dt[, alt_factor := stats::relevel(
    factor(as.character(alt), levels = as.character(seq_len(J))),
    ref = "1"
  )]
  as.data.frame(dt)
}

mnlogit_fit_clogit <- function(setup, x_cols, J, config) {
  form <- stats::as.formula(paste(
    "choice ~", paste(c(x_cols, "alt_factor"), collapse = " + "), "+ survival::strata(id)"
  ))
  fit <- survival::clogit(
    form,
    data = setup,
    method = "exact",
    control = survival::coxph.control(iter.max = as.integer(config$optimizer_control$maxeval %||% 1000L))
  )
  est <- stats::coef(fit)
  list(
    estimates = est,
    loglik = as.numeric(stats::logLik(fit)),
    converged = isTRUE(fit$converged %||% TRUE),
    n_params = length(est)
  )
}

mnlogit_setup_apollo <- function(dt, x_cols, J, config) {
  suppressPackageStartupMessages(library(apollo))
  wide <- mnlogit_wide_data(dt, x_cols, J)
  apollo_mnl <- apollo::apollo_mnl
  apollo_prepareProb <- apollo::apollo_prepareProb
  alternatives <- stats::setNames(seq_len(J), paste0("alt", seq_len(J)))
  avail <- stats::setNames(as.list(rep(1L, J)), names(alternatives))
  apollo_beta <- stats::setNames(
    rep(0, length(x_cols) + J - 1L),
    c(x_cols, paste0("ASC_", 2:J))
  )
  apollo_control <- list(
    modelName = "choicer_mnlogit_benchmark",
    modelDescr = "choicer MNL benchmark",
    indivID = "id",
    nCores = as.integer(config$n_threads),
    noDiagnostics = TRUE,
    outputDirectory = tempdir()
  )

  apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    db <- apollo_inputs$database
    V <- vector("list", J)
    names(V) <- names(alternatives)
    for (j in seq_len(J)) {
      util <- rep(0, nrow(db))
      for (x in x_cols) {
        util <- util + apollo_beta[[x]] * db[[paste0(x, "_", j)]]
      }
      if (j > 1L) util <- util + apollo_beta[[paste0("ASC_", j)]]
      V[[j]] <- util
    }
    mnl_settings <- list(
      alternatives = alternatives,
      avail = avail,
      choiceVar = db$choice_alt,
      utilities = V
    )
    P <- list(model = apollo_mnl(mnl_settings, functionality))
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }

  inputs <- apollo::apollo_validateInputs(
    apollo_beta = apollo_beta,
    apollo_fixed = character(),
    database = as.data.frame(wide),
    apollo_control = apollo_control,
    silent = TRUE
  )

  list(
    beta = apollo_beta,
    fixed = character(),
    probabilities = apollo_probabilities,
    inputs = inputs
  )
}

mnlogit_fit_apollo <- function(setup, x_cols, J, config) {
  fit <- apollo::apollo_estimate(
    setup$beta,
    setup$fixed,
    setup$probabilities,
    setup$inputs,
    estimate_settings = list(
      estimationRoutine = "bfgs",
      maxIterations = as.integer(config$optimizer_control$maxeval %||% 1000L),
      hessianRoutine = "none",
      printLevel = 0L,
      silent = TRUE
    )
  )
  est <- fit$estimate %||% stats::coef(fit)
  list(
    estimates = est,
    loglik = as.numeric(fit$maximum %||% stats::logLik(fit)),
    converged = isTRUE(fit$successfulEstimation %||% TRUE),
    n_params = length(est)
  )
}

mnlogit_setup_mixl <- function(dt, x_cols, J, config) {
  wide <- mnlogit_wide_data(dt, x_cols, J)
  util_lines <- vapply(seq_len(J), function(j) {
    rhs <- paste(sprintf("@beta_%s * $%s_%d", x_cols, x_cols, j), collapse = " + ")
    if (j > 1L) rhs <- paste(sprintf("@ASC_%d", j), rhs, sep = " + ")
    sprintf("U_%d = %s;", j, rhs)
  }, character(1L))
  utility_script <- paste(util_lines, collapse = "\n")
  model_spec <- mixl::specify_model(
    utility_script,
    dataset = as.data.frame(wide),
    model_name = paste0("choicer_mnlogit_benchmark_", Sys.getpid()),
    disable_multicore = TRUE
  )
  list(
    data = as.data.frame(wide),
    model_spec = model_spec,
    start = stats::setNames(
      rep(0, length(x_cols) + J - 1L),
      c(paste0("beta_", x_cols), paste0("ASC_", 2:J))
    ),
    availabilities = mixl::generate_default_availabilities(
      as.data.frame(wide),
      model_spec$num_utility_functions
    )
  )
}

mnlogit_fit_mixl <- function(setup, x_cols, J, config) {
  fit <- mixl::estimate(
    setup$model_spec,
    setup$start,
    setup$data,
    availabilities = setup$availabilities,
    nDraws = 0L,
    num_threads = as.integer(config$n_threads)
  )
  est <- tryCatch(stats::coef(fit), error = function(e) NULL)
  est <- est %||% fit$coefficients %||% fit$estimate
  ll <- tryCatch(as.numeric(stats::logLik(fit)), error = function(e) fit$logLik %||% fit$ll)
  list(
    estimates = est,
    loglik = as.numeric(ll),
    converged = isTRUE(fit$converged %||% TRUE),
    n_params = length(est)
  )
}

mnlogit_setup_dispatch <- function(package) {
  get(paste0("mnlogit_setup_", package), mode = "function")
}

mnlogit_fit_dispatch <- function(package) {
  get(paste0("mnlogit_fit_", package), mode = "function")
}

mnlogit_run_package <- function(package, dt, x_cols, spec, run, config) {
  dep <- mnlogit_dependency(package)
  if (!requireNamespace(dep, quietly = TRUE)) {
    if (!isTRUE(config$skip_missing)) {
      stop("Required benchmark package is not installed: ", dep)
    }
    raw <- mnlogit_raw_row(
      package = package,
      spec = spec,
      run = run,
      status = "not_installed",
      error = sprintf("Package '%s' is not installed.", dep)
    )
    return(list(raw = raw, coef = data.table::data.table()))
  }

  setup_time <- fit_time <- total_time <- NA_real_
  warning_text <- NA_character_
  total_start <- proc.time()[["elapsed"]]

  result <- tryCatch({
    setup_fun <- mnlogit_setup_dispatch(package)
    fit_fun <- mnlogit_fit_dispatch(package)

    setup_q <- bench_quiet(bench_time(setup_fun(dt, x_cols, spec$J, config)))
    setup_time <- setup_q$value$elapsed
    setup <- setup_q$value$value

    fit_q <- bench_quiet(bench_time(fit_fun(setup, x_cols, spec$J, config)))
    fit_time <- fit_q$value$elapsed
    fit <- fit_q$value$value

    all_warnings <- c(setup_q$warnings, fit_q$warnings)
    warning_text <- if (length(all_warnings)) paste(unique(all_warnings), collapse = " | ") else NA_character_
    total_time <- proc.time()[["elapsed"]] - total_start

    raw <- mnlogit_raw_row(
      package = package,
      spec = spec,
      run = run,
      status = "ok",
      converged = fit$converged,
      setup_time_sec = setup_time,
      fit_time_sec = fit_time,
      total_time_sec = total_time,
      loglik = fit$loglik,
      n_params = fit$n_params,
      warnings = warning_text
    )
    coef <- mnlogit_canonicalize(
      package = package,
      estimates = fit$estimates,
      x_cols = x_cols,
      J = spec$J,
      raw_loglik = fit$loglik,
      spec = spec,
      run = run,
      status = "ok"
    )
    list(raw = raw, coef = coef)
  }, error = function(e) {
    total_time <- proc.time()[["elapsed"]] - total_start
    raw <- mnlogit_raw_row(
      package = package,
      spec = spec,
      run = run,
      status = "error",
      setup_time_sec = setup_time,
      fit_time_sec = fit_time,
      total_time_sec = total_time,
      error = conditionMessage(e),
      warnings = warning_text
    )
    list(raw = raw, coef = data.table::data.table())
  })

  result
}
