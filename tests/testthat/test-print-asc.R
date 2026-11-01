# Deterministic fitted summaries: printing must not depend on estimation noise.
asc_print_fit <- function(model = "mnl", labels = sprintf("alt%02d", 1:25)) {
  pm <- list(beta = 1L)
  extra <- switch(model, mxl = "L_11", nl = "lambda_A", character())
  if (model == "mxl") pm$sigma <- 2L
  if (model == "nl") pm$lambda <- 2L
  pm$asc <- length(c("cost", extra)) + seq_along(labels)
  # Deliberately avoid ASC_ names; use the parameter map and alternative labels.
  nm <- c("ASC_covariate", extra, paste0("constant_", labels, recycle0 = TRUE))
  co <- setNames(rep(0.5, length(nm)), nm)
  fit <- list(model = model, coefficients = co, se = rep(0.1, length(co)),
              vcov = diag(0.01, length(co)), param_map = pm,
              alt_mapping = data.table(alt_int = 0:length(labels),
                                       mode = c("base", labels)),
              n_params = length(co), nobs = 200, loglik = -100,
              convergence = 0L, message = "converged",
              optimizer = list(name = "optim", elapsed_time = 0),
              rc_dist = 0L, rc_correlation = FALSE, rc_mean = FALSE,
              sigma = matrix(1, 1, 1), include_outside_option = TRUE)
  class(fit) <- c(paste0("choicer_", model), "choicer_fit")
  fit
}

asc_print_truth <- function(fit) {
  list(beta = 0.5, L_params = 0.5, lambdas = 0.5,
       delta = rep(0.5, length(fit$param_map$asc)))
}

asc_output <- function(x, ...) paste(capture.output(print(x, ...)), collapse = "\n")
asc_lines <- function(x, ...) {
  out <- capture.output(print(x, ...))
  trimws(out[grepl("constant_", out)])
}

test_that("frequentist summaries limit only ASC printing and retain all estimates", {
  for (model in c("mnl", "mxl", "nl")) {
    fit <- asc_print_fit(model)
    s <- summary(fit, gof = FALSE)
    expect_length(asc_lines(s), 20L)
    expect_match(asc_output(s), "20 shown, 5 hidden")
    expect_match(asc_output(s), "ASC_covariate")
    expect_equal(nrow(s$coefficients), length(coef(fit)))
    expect_length(asc_lines(s, print_asc = TRUE), 25L)
    expect_false(grepl("hidden", asc_output(s, print_asc = TRUE)))
    expect_length(asc_lines(s, print_asc = FALSE), 0L)
    expect_match(asc_output(s, print_asc = FALSE), "0 shown, 25 hidden")
    small <- summary(fit, gof = FALSE, print_asc = 1L)
    expect_length(asc_lines(small), 1L)
    expect_equal(small$coefficients, s$coefficients)
    expect_match(asc_output(small), "constant_alt01")
    expect_equal(coef(fit), fit$coefficients)
    expect_match(asc_output(s), "Log-likelihood:")
    if (model == "mxl") expect_match(asc_output(s), "Sigma_11")
    if (model == "nl") expect_match(asc_output(s), "lambda_A")
  }
})

test_that("labels and full-table positions select ASCs in existing order", {
  fit <- asc_print_fit(labels = c("car", "train", "bus", "ferry"))
  s <- summary(fit, gof = FALSE, print_asc = c("ferry", "car", "car"))
  lines <- asc_lines(s)
  expect_length(lines, 2L)
  expect_match(lines[1], "constant_car")
  expect_match(lines[2], "constant_ferry")
  expect_identical(lines, asc_lines(s, print_asc = c(5L, 2L, 5L)))
  expect_match(asc_output(s), "2 shown, 2 hidden")
  expect_length(asc_lines(s, print_asc = "train"), 1L)
  for (selector in list(0L, FALSE, integer(), character())) {
    expect_length(asc_lines(s, print_asc = selector), 0L)
  }
  for (selector in list(4L, 50L, 1e20, TRUE)) {
    expect_length(asc_lines(s, print_asc = selector), 4L)
    expect_false(grepl("hidden", asc_output(s, print_asc = selector)))
  }
  expect_identical(asc_lines(s, print_asc = 2), asc_lines(s, print_asc = 2L))
})

test_that("real fit labels exclude the reference with and without an outside option", {
  d <- create_small_mnl_data()
  d[, alt := c("car", "train", "bus")[alt]]
  for (outside in c(FALSE, TRUE)) {
    fit <- suppressMessages(run_mnlogit(
      d, "id", "alt", "choice", c("x1", "x2"),
      include_outside_option = outside, outside_opt_label = "car"
    ))
    s <- summary(fit, gof = FALSE, print_asc = "train")
    expect_match(asc_output(s), "ASC_train")
    expect_false(grepl("ASC_bus", asc_output(s)))
    expect_error(summary(fit, gof = FALSE, print_asc = "car"), "without an ASC")
    rt <- recovery_table(fit, list(beta = c(0, 0), delta = c(0, 0)),
                         print_asc = "train")
    expect_match(asc_output(rt), "ASC_train")
    expect_false(grepl("ASC_bus", asc_output(rt)))
  }
})

test_that("invalid ASC selectors fail clearly at construction and printing", {
  fit <- asc_print_fit()
  s <- summary(fit, gof = FALSE)
  rt <- recovery_table(fit, asc_print_truth(fit))
  invalid <- list(-1L, 1.5, NA, NA_integer_, NA_character_, Inf, NaN,
                  c(TRUE, FALSE), logical(), c(2, 3.5), c(0L, 2L),
                  c(1L, 2L), c(2L, 99L), "unknown", "constant_alt01",
                  "base", list(2L), NULL, factor("alt01"), matrix(1L), 1+1i)
  for (selector in invalid) {
    expect_error(summary(fit, gof = FALSE, print_asc = selector), "print_asc")
    expect_error(print(s, print_asc = selector), "print_asc")
    expect_error(recovery_table(fit, asc_print_truth(fit), print_asc = selector),
                 "print_asc")
    expect_error(print(rt, print_asc = selector), "print_asc")
  }
})

test_that("empty ASC blocks and ASC-only coefficient tables print cleanly", {
  fit <- asc_print_fit(labels = character())
  for (selector in list(20L, TRUE, FALSE, character(), integer())) {
    s <- summary(fit, gof = FALSE, print_asc = selector)
    expect_match(asc_output(s), "ASC_covariate")
    expect_false(grepl("hidden", asc_output(s)))
    rt <- recovery_table(fit, list(beta = 0.5), print_asc = selector)
    expect_false(grepl("hidden", asc_output(rt)))
  }
  fit <- asc_print_fit()
  fit$coefficients <- fit$coefficients[-1L]
  fit$se <- fit$se[-1L]
  fit$param_map <- list(asc = seq_along(fit$coefficients))
  s <- summary(fit, gof = FALSE, print_asc = FALSE)
  expect_match(asc_output(s), "0 shown, 25 hidden")
  expect_equal(nrow(s$coefficients), 25L)
})

test_that("recovery selectors respect missing truth blocks, subsetting and print options", {
  fit <- asc_print_fit("nl")
  truth <- asc_print_truth(fit)
  rt <- recovery_table(fit, truth, print_asc = 2L)
  original <- data.table::copy(rt)
  expect_length(asc_lines(rt), 2L)
  expect_match(asc_output(rt), "2 shown, 23 hidden")
  expect_match(asc_output(rt), "lambda_A")
  expect_equal(nrow(rt), 27L)
  expect_equal(asc_lines(rt, print_asc = c(27L, 3L)),
               asc_lines(rt, print_asc = c("alt25", "alt01")))
  # A low data.table print limit must not silently undo print_asc = TRUE.
  old_options <- options(datatable.print.nrows = 2L, datatable.print.topn = 1L)
  on.exit(options(old_options), add = TRUE)
  expect_length(asc_lines(rt, print_asc = TRUE), 25L)
  expect_equal(rt, original)
  # Labels remain attached to parameters even after rows move.
  reordered <- rt[c(27L, 3L, 1L)]
  expect_match(asc_output(reordered, print_asc = "alt25"), "constant_alt25")
  expect_false(grepl("constant_alt01", asc_output(reordered, print_asc = "alt25")))
  partial <- recovery_table(fit, list(delta = truth$delta), print_asc = c(1L, 25L))
  expect_length(asc_lines(partial), 2L)
  expect_match(asc_output(partial), "constant_alt25")
})

test_that("MNP summaries and recovery keep covariance rows visible", {
  fit <- asc_print_fit("mnp", c("car", "train", "bus"))
  class(fit) <- "choicer_mnp"
  fit$base_alt <- "base"
  fit$mcmc <- list(R_keep = 10L, R = 10L, burn = 0L, thin = 1L, seed = 1L)
  fit$draws <- list(
    beta = matrix(seq_len(40) / 40, 10, 4,
                  dimnames = list(NULL, names(fit$coefficients))),
    sigma = matrix(rep(c(1, 0, 1, 0, 0, 1), each = 10), 10, 6,
                   dimnames = list(NULL, c("Sigma_11", "Sigma_21", "Sigma_22",
                                          "Sigma_31", "Sigma_32", "Sigma_33"))))
  fit$sigma <- diag(3)
  s <- summary(fit, print_asc = "train")
  expect_length(asc_lines(s), 1L)
  expect_match(asc_output(s), "constant_train")
  expect_match(asc_output(s), "Sigma_33")
  rt <- recovery_table(fit, c(asc_print_truth(fit), list(Sigma = diag(3))),
                       print_asc = FALSE)
  expect_length(asc_lines(rt), 0L)
  expect_match(asc_output(rt), "Sigma_33")
  expect_equal(nrow(rt), 10L)
})

test_that("hierarchical Bayes selectors cover delta and xi without hiding diagnostics", {
  for (model in c("hmnl", "hmnp")) {
    sim <- if (model == "hmnl") {
      simulate_hmnl_data(N = 20, T = 2, J = 4, theta = 0.5, seed = 42)
    } else {
      simulate_hmnp_data(N = 20, T = 2, J = 4, theta = 0.5, seed = 42)
    }
    run <- if (model == "hmnl") run_hmnlogit else run_hmnprobit
    fit <- suppressMessages(suppressWarnings(run(sim$data, "task", "alt", "choice",
                                 c("x1", "x2"), person_col = "pid",
                                 mcmc = list(R = 40, burn = 20, seed = 7))))
    labels <- as.character(fit$delta$alternative)
    s <- summary(fit, print_asc = labels[c(4, 2)])
    expect_match(asc_output(s), "2 shown, 2 hidden")
    expect_equal(asc_output(s), asc_output(s, print_asc = c(4L, 2L)))
    expect_match(asc_output(s, print_asc = FALSE), "Convergence diagnostics")
    expect_equal(nrow(s$delta), 4L)
    expect_equal(nrow(s$xi), 4L)
    expect_false(grepl("hidden", asc_output(s, print_asc = TRUE)))
    rt <- recovery_table(fit, sim, print_asc = labels[2])
    expect_match(asc_output(rt), "1 shown, 3 hidden")
    positions <- which(rt$group == "delta")
    expect_equal(asc_output(rt, print_asc = positions[c(4, 2)]),
                 asc_output(rt, print_asc = labels[c(4, 2)]))
    expect_equal(sum(rt$group == "delta"), 4L)
  }
})

test_that("Monte Carlo recovery retains labels and all rows across failed replications", {
  fit <- asc_print_fit(labels = c("car", "train", "bus"))
  truth <- asc_print_truth(fit)
  fit_fun <- function(sim) {
    if (sim$seed == 1L) stop("failed replication")
    fit
  }
  # Supply the truth through the usual choicer_sim interface.
  sim_fun <- function(seed) structure(list(true_params = truth, seed = seed),
                                     class = "choicer_sim")
  mc <- monte_carlo(sim_fun, fit_fun, R = 3L, progress = FALSE)
  s <- recovery_table(mc, print_asc = "train")
  expect_equal(nrow(s$summary), 4L)
  expect_length(asc_lines(s), 1L)
  expect_match(asc_output(s), "constant_train")
  expect_match(asc_output(s), "1 shown, 2 hidden")
  expect_length(asc_lines(s, print_asc = TRUE), 3L)
  expect_length(asc_lines(summary(mc, print_asc = FALSE)), 0L)
})
