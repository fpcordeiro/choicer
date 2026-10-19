# The fit-time variance of run_mnlogit(), run_mxlogit() and run_nestlogit().
# Its inputs (the weights against se_method and the sampling provenance, the
# cluster labels) are checked before any evaluation of the likelihood, and
# an error in the variance step itself, after the optimizer (a mixed logit
# kernel's memory check, a failed allocation), returns the estimates without
# a variance (vcov NULL, NA standard errors, the error in variance_error)
# and a warning of class choicer_vcov_failed instead of discarding them. The
# variance kernels are mocked to fail with a memory-check message; one test
# reaches the mixed logit's real memory check.

# Ninety choice situations of three alternatives, three per decision maker
vf_data <- function(seed = 7L) {
  set.seed(seed)
  dt <- data.table::data.table(id = rep(1:90, each = 3L), alt = rep(1:3, 90L))
  dt[, `:=`(x1 = stats::rnorm(.N), x2 = stats::rnorm(.N))]
  dt[, nest := c("A", "A", "B")[alt]]
  dt[, person := (id - 1L) %/% 3L + 1L]
  dt[, V := 0.8 * x1 - 0.5 * x2]
  dt[, choice := as.integer(alt == sample(alt, 1L, prob = exp(V) / sum(exp(V)))),
     by = id]
  dt[, V := NULL]
  dt[]
}

vf_fit <- function(model, dt, ...) {
  suppressMessages(switch(model,
    mnl = run_mnlogit(data = dt, id_col = "id", alt_col = "alt",
                      choice_col = "choice", covariate_cols = c("x1", "x2"), ...),
    mxl = run_mxlogit(data = dt, id_col = "id", alt_col = "alt",
                      choice_col = "choice", covariate_cols = "x1",
                      random_var_cols = "x2", rc_mean = TRUE, S = 20L, ...),
    mxl_noasc = run_mxlogit(data = dt, id_col = "id", alt_col = "alt",
                            choice_col = "choice", covariate_cols = "x1",
                            random_var_cols = "x2", rc_mean = TRUE, S = 20L,
                            use_asc = FALSE, ...),
    nl = run_nestlogit(data = dt, id_col = "id", alt_col = "alt",
                       choice_col = "choice", covariate_cols = c("x1", "x2"),
                       nest_col = "nest", ...)
  ))
}

# Counts the evaluations of the three models' likelihood (the kernel behind
# the optimizer and the "bhhh" pass) in `counter$n` until the calling test
# ends
vf_count_evals <- function(env = parent.frame()) {
  counter <- new.env()
  counter$n <- 0L
  count <- function(name) {
    real <- get(name, envir = asNamespace("choicer"))
    function(...) {
      counter$n <- counter$n + 1L
      real(...)
    }
  }
  local_mocked_bindings(
    mnl_loglik_gradient_parallel = count("mnl_loglik_gradient_parallel"),
    mxl_loglik_gradient_parallel = count("mxl_loglik_gradient_parallel"),
    nl_loglik_gradient_parallel = count("nl_loglik_gradient_parallel"),
    .env = env
  )
  counter
}

# Prepared inputs of each model ("mxl_panel": decision makers by `person`)
vf_prepare <- function(model, dt, ...) {
  suppressWarnings(switch(model,
    mnl = prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"), ...),
    mxl = prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2", ...),
    mxl_panel = prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2",
                                 person_col = "person", ...),
    nl = prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest", ...)
  ))
}

vf_fit_input <- function(model, input, ...) {
  run <- switch(sub("_panel", "", model), mnl = run_mnlogit,
                mxl = run_mxlogit, nl = run_nestlogit)
  args <- list(input_data = input, ...)
  if (startsWith(model, "mxl")) {
    args$eta_draws <- get_halton_normals(5L, length(choicer:::.unit_first(input)), 1L)
  }
  suppressMessages(do.call(run, args))
}

test_that("fit-time cluster labels are checked before any likelihood evaluation", {
  evals <- vf_count_evals()
  dt <- vf_data()
  dt[, team := (person - 1L) %/% 2L + 1L]  # clusters that nest decision makers
  for (model in c("mnl", "mxl", "nl", "mxl_panel")) {
    input <- vf_prepare(model, dt, cluster_col = "team")
    # well-formed labels fit, with their variance
    evals$n <- 0L
    fit <- vf_fit_input(model, input, se_method = "cluster")
    expect_s3_class(fit, "choicer_fit")
    expect_null(fit$variance_error)
    expect_gt(evals$n, 0L)

    bad <- input
    bad$cluster[2L] <- NA
    evals$n <- 0L
    expect_error(vf_fit_input(model, bad, se_method = "cluster"),
                 if (model == "mxl_panel") {
                   "must be constant within each decision maker"
                 } else {
                   "`cluster` contains missing values."
                 }, fixed = TRUE, label = model)
    expect_identical(evals$n, 0L, label = paste(model, "NA"))

    bad <- input
    bad$cluster <- bad$cluster[-1L]
    evals$n <- 0L
    expect_error(vf_fit_input(model, bad, se_method = "cluster"),
                 if (model == "mxl_panel") {
                   "`cluster_col`: got 89 values for 90 choice situations."
                 } else {
                   "`cluster` has length 89 but the model has 90 choice situations."
                 }, fixed = TRUE, label = model)
    expect_identical(evals$n, 0L, label = paste(model, "length"))
  }
  # labels with an NA for a whole decision maker reach the NA check
  bad <- vf_prepare("mxl_panel", dt, cluster_col = "team")
  bad$cluster[1:3] <- NA
  evals$n <- 0L
  expect_error(vf_fit_input("mxl_panel", bad, se_method = "cluster"),
               "`cluster` contains missing values.", fixed = TRUE)
  expect_identical(evals$n, 0L)
  # a panel whose Ti does not split the situations is the kernels' to
  # report, at their first call, whatever the labels
  bad <- vf_prepare("mxl_panel", dt, cluster_col = "team")
  bad$Ti <- bad$Ti[-1L]
  bad$cluster[2L] <- NA
  evals$n <- 0L
  expect_error(vf_fit_input("mxl_panel", bad, se_method = "cluster"),
               "sum(Ti) (87) does not match the number of choice situations (90).",
               fixed = TRUE)
  expect_identical(evals$n, 1L)  # the "bhhh" pass
})

test_that("the weights' guard runs before any likelihood evaluation", {
  evals <- vf_count_evals()
  dt <- vf_data()
  # a prepared input that carries choice-based-sampling provenance with
  # uniform weights stops
  flagged <- data.table::copy(dt)
  flagged[, w := 2]
  data.table::setattr(flagged, "choice_sampling", structure(
    list(scheme = "wesml", weight_name = "w"), class = "choicer_sampling"))
  for (model in c("mnl", "mxl", "nl")) {
    input <- vf_prepare(model, flagged, weights_col = "w")
    evals$n <- 0L
    expect_error(vf_fit_input(model, input, se_method = "sandwich"),
                 "`input_data` is flagged as a WESML choice-based sample",
                 fixed = TRUE, label = model)
    expect_identical(evals$n, 0L, label = model)
    # before the cluster labels, as when both ran after the optimizer
    input <- vf_prepare(model, flagged, weights_col = "w", cluster_col = "person")
    input$cluster[2L] <- NA
    evals$n <- 0L
    expect_error(vf_fit_input(model, input, se_method = "cluster"),
                 "`input_data` is flagged as a WESML choice-based sample",
                 fixed = TRUE, label = model)
    expect_identical(evals$n, 0L, label = model)
  }
  # non-uniform weights without the sandwich or clustered variance warn
  weighted <- data.table::copy(dt)
  weighted[, w := 1 + id %% 2]
  for (model in c("mnl", "mxl", "nl")) {
    for (se_method in c("hessian", "bhhh")) {
      evals$n <- 0L
      at <- integer(0)
      fit <- withCallingHandlers(
        vf_fit(model, weighted, weights_col = "w", se_method = se_method),
        warning = function(w) {
          if (startsWith(conditionMessage(w), "Non-uniform weights detected")) {
            at <<- c(at, evals$n)
            invokeRestart("muffleWarning")
          }
        })
      expect_s3_class(fit, "choicer_fit")
      expect_null(fit$variance_error)
      expect_identical(at, 0L, label = paste(model, se_method))
      expect_gt(evals$n, 0L)
    }
  }
})

# The models and their se_method; mxl_noasc is the mixed logit without
# constants (use_asc = FALSE), whose parameter vector has no asc block
vf_methods <- list(
  mnl = c("hessian", "bhhh", "sandwich", "cluster"),
  mxl = c("hessian", "bhhh", "sandwich", "cluster"),
  mxl_noasc = c("hessian", "bhhh", "sandwich", "cluster"),
  nl = c("hessian", "numeric", "bhhh", "sandwich", "cluster")
)

# The kernel behind each model's fit-time variance (the meat's, for the
# sandwich), and the vcov(type = ) that recomputes it
vf_kernel <- function(model, se_method) {
  switch(paste(sub("_noasc$", "", model), se_method),
    "mnl hessian" = "mnl_loglik_hessian_parallel",
    "mnl bhhh" = , "mnl sandwich" = "mnl_bhhh_parallel",
    "mnl cluster" = "mnl_scores_parallel",
    "nl hessian" = "nl_loglik_hessian_parallel",
    "nl numeric" = "nl_loglik_numeric_hessian",
    "nl bhhh" = , "nl sandwich" = "nl_bhhh_parallel",
    "nl cluster" = "nl_scores_parallel",
    "mxl hessian" = "mxl_hessian_parallel",
    "mxl bhhh" = , "mxl sandwich" = "mxl_bhhh_parallel",
    "mxl cluster" = "mxl_cluster_meat_parallel")
}
vf_type <- function(se_method) {
  switch(se_method, sandwich = "robust", numeric = "hessian", se_method)
}

vf_memory_stop <- paste(
  "The Hessian needs about 9.80 GB with 2 threads (4.90 GB per thread and",
  "70 kB shared), more than the 8.00 GB of physical memory. Use",
  "set_num_threads(1), or raise options(choicer.max_memory =).")

# The value of `expr` and the choicer_vcov_failed warnings it raised (muffled)
vf_capture <- function(expr) {
  warnings <- list()
  value <- withCallingHandlers(expr, choicer_vcov_failed = function(w) {
    warnings[[length(warnings) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

# A fit whose variance kernel stops with `msg` (NULL: not mocked here); the
# mock ends when this returns
vf_fit_failing <- function(model, se_method, dt, ..., msg = vf_memory_stop) {
  if (!is.null(msg)) {
    mock <- stats::setNames(list(function(...) stop(msg)),
                            vf_kernel(model, se_method))
    do.call(local_mocked_bindings, c(mock, list(.env = environment())))
  }
  vf_capture(vf_fit(model, dt, se_method = se_method, ...))
}

# A fit object without its variance, its record of a failed variance step
# and its run time
vf_strip <- function(fit) {
  fit$vcov <- NULL
  fit$se <- NULL
  fit$variance_error <- NULL
  fit$optimizer$elapsed_time <- NULL
  fit
}

# What a fit whose variance step failed with `msg` holds instead
vf_expect_no_vcov <- function(fit, ref, msg = vf_memory_stop, label = NULL) {
  nm <- names(ref$coefficients)
  expect_true("vcov" %in% names(fit), label = label)
  expect_null(fit$vcov, label = label)
  expect_identical(fit$se, stats::setNames(rep(NA_real_, length(nm)), nm),
                   label = label)
  expect_identical(fit$variance_error, msg, label = label)
}

test_that("a failed variance step keeps the estimates, without a variance and with a warning", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)  # so that the two fits agree bit for bit
  dt <- vf_data()
  for (model in names(vf_methods)) {
    for (se_method in vf_methods[[model]]) {
      label <- paste(model, se_method)
      cluster_col <- if (se_method == "cluster") "person"
      ref <- vf_fit(model, dt, se_method = se_method, cluster_col = cluster_col)
      expect_true(all(is.finite(ref$se)), label = label)
      expect_false("variance_error" %in% names(ref))
      res <- vf_fit_failing(model, se_method, dt, cluster_col = cluster_col)
      fit <- res$value

      # One warning, quoting the error (whose own text says how much memory
      # and how many threads) and naming the post-hoc route
      expect_length(res$warnings, 1L)
      w <- res$warnings[[1L]]
      expect_null(conditionCall(w))
      expect_identical(conditionMessage(w), sprintf(paste0(
        "The variance could not be computed (se_method = \"%s\"), so the ",
        "estimates are returned without it, with NA standard errors: %s To ",
        "compute it from the fit, call vcov(fit, type = \"%s\")."),
        se_method, vf_memory_stop, vf_type(se_method)), label = label)

      # The fit is the successful one but for its variance, which it lacks
      expect_identical(vf_strip(fit), vf_strip(ref), label = label)
      vf_expect_no_vcov(fit, ref, label = label)

      # vcov(type = ) computes the variance from the fit, as from the
      # successful one
      expect_identical(vcov(fit, type = vf_type(se_method)),
                       vcov(ref, type = vf_type(se_method)), label = label)
      expect_equal(vcov(fit, type = vf_type(se_method)), ref$vcov,
                   tolerance = 1e-10, label = label)
    }
  }
})

test_that("summary(), print() and post-estimation work on a fit without a variance", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  dt <- vf_data()
  for (model in names(vf_methods)) {
    ref <- vf_fit(model, dt)
    fit <- vf_fit_failing(model, "hessian", dt)$value

    # The methods below never run the failed step again (ensure_vcov()
    # leaves the fit alone): the kernel, no longer mocked, is not called
    hess_calls <- 0L
    kernel <- vf_kernel(model, "hessian")
    real <- get(kernel, envir = asNamespace("choicer"))
    do.call(local_mocked_bindings, stats::setNames(list(function(...) {
      hess_calls <<- hess_calls + 1L
      real(...)
    }), kernel))

    expect_null(vcov(fit), label = model)
    expect_identical(coef(fit), coef(ref))
    expect_identical(logLik(fit), logLik(ref))
    expect_identical(AIC(fit), AIC(ref))
    expect_identical(BIC(fit), BIC(ref))
    expect_identical(nobs(fit), nobs(ref))

    # summary(): the estimates of the successful fit (the mixed logit's
    # covariance block as Sigma), every standard error NA, and a note
    s <- summary(fit)
    s_ref <- summary(ref)
    expect_identical(s$coefficients$Estimate, s_ref$coefficients$Estimate,
                     label = model)
    expect_identical(rownames(s$coefficients), rownames(s_ref$coefficients))
    expect_true(all(is.na(s$coefficients$Std_Error)), label = model)
    expect_true(all(is.na(s$coefficients$Pr_z)))
    expect_identical(s$gof, s_ref$gof)
    expect_identical(s$variance_error, vf_memory_stop)
    expect_false("variance_error" %in% names(s_ref))
    out <- capture.output(print(s))
    expect_true(any(grepl(paste0("not computed at fit time .*",
                                 "vcov\\(fit, type = \"hessian\"\\) computes them"),
                          out)), label = model)
    expect_false(any(grepl("not computed", capture.output(print(s_ref)))))
    expect_output(print(fit), "Log-likelihood")

    # Post-estimation: the point estimates of the successful fit, NA
    # standard errors where they need the variance
    expect_identical(predict(fit), predict(ref))
    expect_identical(predict(fit, type = "shares"), predict(ref, type = "shares"))
    expect_identical(logsum(fit), logsum(ref))
    w_fit <- wtp(fit, "x1")
    w_ref <- wtp(ref, "x1")
    expect_identical(w_fit$Estimate, w_ref$Estimate, label = model)
    expect_true(all(is.na(w_fit$Std_Error)), label = model)
    expect_true(all(is.finite(w_ref$Std_Error)))
    cs <- consumer_surplus(fit, "x1")
    expect_identical(cs$mean_cs, consumer_surplus(ref, "x1")$mean_cs)
    expect_true(is.na(cs$se_mean_cs))
    expect_output(print(cs), "Consumer surplus")
    truth <- list(beta = c(0.8, -0.5)[seq_along(ref$param_map$beta)])
    rt <- recovery_table(fit, truth)
    expect_identical(rt$estimate, recovery_table(ref, truth)$estimate)
    expect_true(all(is.na(rt$se)))
    expect_identical(hess_calls, 0L, label = model)
  }
})

test_that("summary() of a mixed logit without a variance shows the covariance", {
  # A singular information (x2 duplicates x1, so the fixed coefficient and
  # the random coefficient's mean are not separately identified) leaves vcov
  # NULL, as a failed variance step does. summary() reported the Cholesky
  # parameters themselves under the Sigma labels; it now reports Sigma.
  dt <- vf_data()
  dt[, x2 := x1]
  fit <- suppressMessages(vf_fit("mxl", dt))
  expect_null(fit$vcov)
  expect_false("variance_error" %in% names(fit))
  s <- suppressMessages(summary(fit, gof = FALSE))
  est <- stats::setNames(s$coefficients$Estimate, rownames(s$coefficients))
  expect_equal(est[["Sigma_11"]], exp(2 * fit$coefficients[["L_11"]]))
  expect_identical(est[["Mu_x2"]], fit$coefficients[["Mu_x2"]])
  expect_true(all(is.na(s$coefficients$Std_Error)))
})

test_that("a fit without its data says that computing the variance takes a refit", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  dt <- vf_data()
  for (model in names(vf_methods)) {
    ref <- vf_fit(model, dt, keep_data = FALSE)
    res <- vf_fit_failing(model, "hessian", dt, keep_data = FALSE)
    expect_length(res$warnings, 1L)
    expect_match(conditionMessage(res$warnings[[1L]]),
                 paste("The fit keeps no data (keep_data = FALSE): refit with",
                       "keep_data = TRUE to compute it with",
                       "vcov(fit, type = \"hessian\")"), fixed = TRUE)
    fit <- res$value
    expect_output(print(summary(fit)), paste(
      "refit with keep_data = TRUE for vcov\\(fit, type = \"hessian\"\\)",
      "to compute them"))
    # the data slot is kept, NULL, as in the successful fit
    expect_identical(names(fit), c(names(ref), "variance_error"), label = model)
    expect_null(fit$data)
    a <- vf_strip(fit)
    b <- vf_strip(ref)
    a$call <- b$call <- NULL  # called with and without se_method
    expect_identical(a, b, label = model)
    vf_expect_no_vcov(fit, ref, label = model)
  }
})

test_that("the error is quoted whole, with a full stop", {
  dt <- vf_data()
  res <- vf_fit_failing("mnl", "hessian", dt,
                        msg = "cannot allocate vector of size 2.3 Gb")
  expect_match(conditionMessage(res$warnings[[1L]]),
               "standard errors: cannot allocate vector of size 2.3 Gb. To compute",
               fixed = TRUE)
  expect_identical(res$value$variance_error, "cannot allocate vector of size 2.3 Gb.")
  res <- vf_fit_failing("mnl", "hessian", dt, msg = "")
  expect_identical(res$value$variance_error, "an error of class simpleError.")
})

test_that("an inversion that fails for want of memory takes the same route", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  dt <- vf_data()
  ref <- vf_fit("mnl", dt)
  # The fit's messages, recorded (vf_fit() suppresses them)
  fit_msgs <- function() {
    msgs <- character()
    res <- vf_capture(withCallingHandlers(
      run_mnlogit(data = dt, id_col = "id", alt_col = "alt",
                  choice_col = "choice", covariate_cols = c("x1", "x2")),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }))
    c(res, list(messages = msgs))
  }
  # An allocation in the inversion was reported as a singular information
  # (vcov NULL, which ensure_vcov() then recomputed on every summary())
  oom <- "cannot allocate a working copy of 9.9 Gb for the inversion"
  local_mocked_bindings(equilibrated_inverse = function(...) stop(oom))
  res <- fit_msgs()
  expect_length(res$warnings, 1L)
  expect_false(any(grepl("singular", res$messages)))
  expect_identical(vf_strip(res$value), vf_strip(ref))
  vf_expect_no_vcov(res$value, ref, msg = paste0(oom, "."))
  # post hoc, the error is the caller's
  expect_error(vcov(ref, type = "hessian"), oom, fixed = TRUE)
  # a singular information is still a singular information: NA standard
  # errors with a message, and no warning
  local_mocked_bindings(equilibrated_inverse = function(...) {
    stop("Lapack routine dgetrf: system is exactly singular: U[2,2] = 0")
  })
  res <- fit_msgs()
  expect_length(res$warnings, 0L)
  expect_true(any(grepl("likely singular", res$messages)))
  expect_null(res$value$vcov)
  expect_false("variance_error" %in% names(res$value))
})

test_that("a sandwich that fails for want of memory takes the same route", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  dt <- vf_data()
  ref <- vf_fit_failing("mnl", "sandwich", dt, msg = NULL)$value  # the same call
  oom <- "cannot allocate vector of size 9.9 Gb"
  local_mocked_bindings(equilibrated_sandwich = function(...) stop(oom))
  res <- vf_fit_failing("mnl", "sandwich", dt, msg = NULL)
  expect_length(res$warnings, 1L)
  expect_identical(vf_strip(res$value), vf_strip(ref))
  vf_expect_no_vcov(res$value, ref, msg = paste0(oom, "."))
  expect_error(vcov(ref, type = "robust"), oom, fixed = TRUE)
})

test_that("errors outside the variance step still stop the fit", {
  dt <- vf_data()
  failing_optimizer <- function(theta_init, eval_f, ...) stop("the optimizer failed")
  for (model in names(vf_methods)) {
    expect_error(vf_fit(model, dt, optimizer = failing_optimizer),
                 "the optimizer failed", label = model)
  }
  # the "bhhh" pass of scale_vars = "bhhh", before the optimizer, alone
  for (model in c("mnl", "mxl")) {
    kernel <- paste0(model, "_loglik_gradient_parallel")
    real <- get(kernel, envir = asNamespace("choicer"))
    expect_error(with_mocked_bindings(
      vf_fit(model, dt),
      !!kernel := function(..., opg_diag = FALSE) {
        if (opg_diag) stop("the bhhh pass failed")
        real(..., opg_diag = opg_diag)
      }
    ), "the bhhh pass failed", label = model)
  }
})

test_that("the mixed logit's memory check stopping the variance keeps the estimates", {
  skip_on_cran()
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)  # the check's estimate depends on the threads
  # 400 alternatives, four in each choice situation. Both kernels hold a
  # draw batch of up to 4 MiB per thread; the Hessian also holds its
  # 401 x 401 result (1.3 MB), so a limit between the two needs stops the
  # variance and nothing before it
  set.seed(19)
  N <- 400L
  dt <- data.table::data.table(id = rep(seq_len(N), each = 4L))
  dt[, alt := sample.int(400L, 4L), by = id]
  dt[, `:=`(x1 = stats::rnorm(.N), x2 = stats::rnorm(.N))]
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  input <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2")
  eta <- get_halton_normals(10L, N, 1L)
  fit_at <- function(limit) {
    old <- options(choicer.max_memory = limit)
    on.exit(options(old))
    suppressMessages(vf_capture(run_mxlogit(
      input_data = input, eta_draws = eta, control = list(maxeval = 20L)
    )))
  }
  # What the gradient (with the B pass's diagonal) and the Hessian need, from
  # their own messages at a limit of one byte
  need <- function(kernel, ...) {
    old <- options(choicer.max_memory = 1)
    on.exit(options(old))
    msg <- tryCatch(kernel(...), error = conditionMessage)
    m <- regmatches(msg, regexec("needs about ([0-9.]+) (kB|MB|GB) even with one thread",
                                 msg))[[1L]]
    as.numeric(m[2L]) * c(kB = 1e3, MB = 1e6, GB = 1e9)[[m[3L]]]
  }
  ref <- fit_at(NULL)$value
  args <- list(theta = unname(ref$coefficients), X = input$X, W = input$W,
               alt_idx = input$alt_idx, choice_idx = input$choice_idx,
               M = input$M, weights = input$weights, eta_draws = eta,
               rc_dist = 0L, rc_correlation = FALSE, rc_mean = FALSE)
  grad_need <- do.call(need, c(list(choicer:::mxl_loglik_gradient_parallel),
                               args, list(opg_diag = TRUE)))
  hess_need <- do.call(need, c(list(choicer:::mxl_hessian_parallel), args))
  # the messages round to a tenth of a megabyte
  expect_gt(hess_need - grad_need, 1e6)

  res <- fit_at((grad_need + hess_need) / 2)
  expect_length(res$warnings, 1L)
  expect_match(conditionMessage(res$warnings[[1L]]),
               "standard errors: The Hessian needs about [0-9.]+ MB even with one thread")
  fit <- res$value
  expect_identical(vf_strip(fit), vf_strip(ref))
  expect_null(fit$vcov)
  expect_match(fit$variance_error, "^The Hessian needs about [0-9.]+ MB even with one thread")
  # With the limit lifted, vcov(type = "hessian") is the successful fit's
  # variance, bit for bit (an information that 20 evaluations leave
  # indefinite says so in a message)
  expect_identical(suppressMessages(vcov(fit, type = "hessian")), ref$vcov)
})
