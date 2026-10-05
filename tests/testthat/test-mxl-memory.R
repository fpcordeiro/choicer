# The mixed logit estimation kernels' memory check. Before allocating, each
# kernel compares an upper estimate of its working memory (scratch for each
# thread plus what the threads share) with options(choicer.max_memory) or,
# when that is unset, the machine's physical memory, and stops with an error
# that gives the number of threads that fits. Kernel calls go through
# mxlp_call() (helper-mxl-panel.R).

# The estimate that a call with a tiny limit reports for one thread, in bytes.
mxm_one_thread <- function(kernel, fx) {
  old <- options(choicer.max_memory = 1)
  on.exit(options(old), add = TRUE)
  msg <- tryCatch(mxlp_call(kernel, fx), error = conditionMessage)
  m <- regmatches(msg, regexec("needs about ([0-9.]+) (kB|MB|GB) even with one thread", msg))[[1]]
  as.numeric(m[2L]) * c(kB = 1e3, MB = 1e6, GB = 1e9)[[m[3L]]]
}

test_that("the physical memory is read on every platform", {
  ram <- test_physical_memory()
  expect_gt(ram, 1e8)
  expect_lt(ram, 1e16)
})

test_that("every kernel stops before allocating more than choicer.max_memory", {
  fx <- mxlp_fixture("memory", 3101)
  old <- options(choicer.max_memory = 1)
  on.exit(options(old), add = TRUE)
  for (k in mxlp_kernels) {
    expect_error(mxlp_call(k, fx),
                 "needs about [0-9.]+ MB even with one thread .* more than the 0 kB of options\\(choicer.max_memory\\)",
                 label = k)
  }
  # Inf lifts the check, and a large limit does not bind: the results are
  # those without the option (one thread, so that they repeat bit for bit)
  on.exit(mxlp_threads(2L), add = TRUE)
  mxlp_threads(1L)
  options(choicer.max_memory = NULL)
  ref <- lapply(setNames(mxlp_kernels, mxlp_kernels), function(k) mxlp_call(k, fx))
  for (lim in list(Inf, 1e15)) {
    options(choicer.max_memory = lim)
    for (k in mxlp_kernels) {
      expect_identical(mxlp_call(k, fx), ref[[k]], label = paste(k, "with", format(lim)))
    }
  }
})

test_that("choicer.max_memory must be a positive number of bytes", {
  fx <- mxlp_fixture("memory", 3101)
  old <- options(choicer.max_memory = NULL)
  on.exit(options(old), add = TRUE)
  bad_values <- list(0, -1, NA_real_, NaN, NA_integer_, "1e9", c(1e9, 2e9), TRUE)
  if (requireNamespace("bit64", quietly = TRUE)) {
    bad_values <- c(bad_values, list(bit64::NA_integer64_, bit64::as.integer64(-5)))
  }
  for (bad in bad_values) {
    options(choicer.max_memory = bad)
    expect_error(mxlp_call("gradient", fx),
                 "options(choicer.max_memory =) must be a positive number of bytes",
                 fixed = TRUE)
  }
  options(choicer.max_memory = 2000000000L)  # an integer number of bytes
  expect_no_error(mxlp_call("gradient", fx))
  if (requireNamespace("bit64", quietly = TRUE)) {
    # an integer64 number of bytes is read as its value, not its bit pattern
    options(choicer.max_memory = bit64::as.integer64(8e9))
    expect_no_error(mxlp_call("gradient", fx))
    options(choicer.max_memory = bit64::as.integer64(1))
    expect_error(mxlp_call("gradient", fx), "more than the 0 kB of", fixed = TRUE)
  }
  # something classed integer64 but not stored as one is read by its storage
  # (an integer number of bytes) or refused, without reading it as doubles
  options(choicer.max_memory = structure(5L, class = "integer64"))
  expect_error(mxlp_call("gradient", fx), "more than the 0 kB of", fixed = TRUE)
  options(choicer.max_memory = structure("5", class = "integer64"))
  expect_error(mxlp_call("gradient", fx),
               "options(choicer.max_memory =) must be a positive number of bytes",
               fixed = TRUE)
})

test_that("the Hessian counts a split unit's per-situation draws", {
  # One situation of 2,000 alternatives and 600 draws, without ASCs: the
  # unit is split into draw batches, and pass 2 forms the situation's
  # 2,000 x 600 W_t Gamma (9.6 MB), outside the batch budget.
  fx <- mxlp_cross_section(mxlp_fixture("split unit", 3103, J = 2000L, S = 600L,
                                        use_asc = FALSE))
  expect_gt(mxm_one_thread("hessian", fx), 8 * 2000 * 600)
})

test_that("the memory error gives the number of threads that fits", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  skip_if(isTRUE(thread_info()$omp_get_thread_limit < 2L), "thread limit below 2")
  on.exit(mxlp_threads(2L), add = TRUE)
  mxlp_threads(2L)
  # 1,000 alternatives: the BHHH's result (8 MB) and each thread's triangle
  # (4 MB) are a large part of the estimate, so two threads need about twice
  # what one does
  fx <- mxlp_fixture("memory, many alternatives", 3102, J = 1000L)
  one <- mxm_one_thread("bhhh", fx)
  expect_gt(one, 5e6)
  old <- options(choicer.max_memory = 1.3 * one)
  on.exit(options(old), add = TRUE)
  expect_error(mxlp_call("bhhh", fx),
               "with 2 threads .* Use set_num_threads\\(1\\), or raise")
  mxlp_threads(1L)
  B <- mxlp_call("bhhh", fx)
  options(choicer.max_memory = NULL)
  expect_identical(B, mxlp_call("bhhh", fx))
})

test_that("a call beyond the physical memory stops before allocating", {
  skip_on_cran()  # builds a 10^6-row design
  skip_if_not(test_physical_memory() > 0, "physical memory unknown")
  skip_if(test_physical_memory() > 4e12, "more than 4 TB of memory")
  old <- options(choicer.max_memory = NULL)
  on.exit(options(old), add = TRUE)
  # One choice situation offering 10^6 alternatives: the BHHH's 10^6 x 10^6
  # result alone would take 8 TB.
  J <- 1e6L
  fx <- list(X = matrix(stats::rnorm(2 * J), J, 2L), W = matrix(0, J, 1L),
             alt_idx = seq_len(J), choice_idx = 1L, M = J, weights = 1,
             Ti = NULL, eta = get_halton_normals(2L, 1L, 1L), K_w = 1L, S = 2L,
             theta = c(0.1, -0.1, 0, rep(0, J - 1L)), rc_dist = 0L,
             rc_correlation = FALSE, rc_mean = FALSE, use_asc = TRUE,
             include_outside_option = FALSE)
  expect_error(mxlp_call("bhhh", fx), "GB of physical memory")
})
