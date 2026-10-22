# The multinomial and nested logit estimation kernels' memory check (the
# mixed logit's, test-mxl-memory.R). Before allocating, each kernel compares
# an upper estimate of its working memory (each thread's buffers plus what the
# threads share) with options(choicer.max_memory) or, when that is unset, the
# machine's physical memory, and stops with an error that gives the number of
# threads that fits. The check follows the kernels' input checks.

# N choice situations, each offering m of J alternatives (distinct codes, the
# first situation offering every one), K covariates, random choices.
mnm_fixture <- function(seed, N = 40L, J = 12L, m = 5L, K = 2L) {
  set.seed(seed)
  alts <- c(list(seq_len(J)),
            lapply(seq_len(N - 1L), function(i) sort(sample.int(J, m))))
  M <- lengths(alts)
  rows <- sum(M)
  list(X = matrix(stats::rnorm(rows * K), rows, K),
       alt_idx = as.integer(unlist(alts)),
       choice_idx = vapply(M, function(mm) sample.int(mm, 1L), integer(1)),
       M = as.integer(M), weights = stats::runif(N, 0.5, 2),
       theta = c(stats::rnorm(K, sd = 0.3), stats::rnorm(J - 1L, sd = 0.3)))
}

mnm_kernels <- c("gradient", "gradient_opg", "hessian", "bhhh", "scores", "meat")

mnm_call <- function(kernel, f) {
  ns <- asNamespace("choicer")
  switch(kernel,
    gradient = ns$mnl_loglik_gradient_parallel(f$theta, f$X, f$alt_idx,
                                               f$choice_idx, f$M, f$weights),
    gradient_opg = ns$mnl_loglik_gradient_parallel(f$theta, f$X, f$alt_idx,
                                                   f$choice_idx, f$M, f$weights,
                                                   opg_diag = TRUE),
    hessian = ns$mnl_loglik_hessian_parallel(f$theta, f$X, f$alt_idx,
                                             f$choice_idx, f$M, f$weights),
    bhhh = ns$mnl_bhhh_parallel(f$theta, f$X, f$alt_idx, f$choice_idx, f$M,
                                f$weights),
    scores = ns$mnl_scores_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
                                    f$M),
    meat = ns$mnl_cluster_meat_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
                                        f$M, f$weights,
                                        rep(1:4, length.out = length(f$M))))
}

# The estimate that a call with a tiny limit reports for one thread, in bytes.
mnm_one_thread <- function(kernel, f) {
  old <- options(choicer.max_memory = 1)
  on.exit(options(old), add = TRUE)
  msg <- tryCatch(mnm_call(kernel, f), error = conditionMessage)
  m <- regmatches(msg, regexec("needs about ([0-9.]+) (kB|MB|GB) even with one thread", msg))[[1]]
  as.numeric(m[2L]) * c(kB = 1e3, MB = 1e6, GB = 1e9)[[m[3L]]]
}

test_that("every MNL kernel stops before allocating more than choicer.max_memory", {
  f <- mnm_fixture(4101)
  old <- options(choicer.max_memory = 1)
  on.exit(options(old), add = TRUE)
  for (k in mnm_kernels) {
    expect_error(mnm_call(k, f),
                 "needs about [0-9.]+ kB even with one thread .* more than the 0 kB of options\\(choicer.max_memory\\)",
                 label = k)
  }
  # Inf lifts the check, and a large limit does not bind: the results are
  # those without the option (one thread, so that they repeat bit for bit)
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  options(choicer.max_memory = NULL)
  ref <- lapply(setNames(mnm_kernels, mnm_kernels), function(k) mnm_call(k, f))
  for (lim in list(Inf, 1e15)) {
    options(choicer.max_memory = lim)
    for (k in mnm_kernels) {
      expect_identical(mnm_call(k, f), ref[[k]], label = paste(k, "with", format(lim)))
    }
  }
})

test_that("the MNL kernels check their inputs before their memory", {
  f <- mnm_fixture(4102)
  old <- options(choicer.max_memory = 1)
  on.exit(options(old), add = TRUE)
  bad_choice <- f
  bad_choice$choice_idx[2L] <- 99L
  for (k in setdiff(mnm_kernels, "hessian")) {
    expect_error(mnm_call(k, bad_choice), "Invalid chosen alternative index",
                 label = k)
  }
  bad_m <- f
  bad_m$M[1L] <- 0L
  for (k in mnm_kernels) {
    expect_error(mnm_call(k, bad_m), "M must be positive", label = k)
  }
  # The Hessian does not read the choices: its memory check comes next
  expect_error(mnm_call("hessian", bad_choice), "needs about")
})

test_that("the MNL memory estimates grow with the parameters and the threads", {
  f_small <- mnm_fixture(4103, J = 12L)
  f_large <- mnm_fixture(4103, J = 400L)
  for (k in c("hessian", "bhhh")) {
    small <- mnm_one_thread(k, f_small)
    large <- mnm_one_thread(k, f_large)
    # at least the n x n result (8 n^2 bytes; n = 401)
    expect_gt(large, 8 * 401^2, label = k)
    expect_lt(small, large, label = k)
  }
  expect_gt(mnm_one_thread("scores", f_large), 16 * 40 * 401)
})

test_that("the MNL memory error gives the number of threads that fits", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  skip_if(isTRUE(thread_info()$omp_get_thread_limit < 2L), "thread limit below 2")
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(2L)
  # 1,000 alternatives: one thread adds into the 8 MB result itself, while
  # two threads also keep a 4 MB packed triangle each, so two threads need
  # about twice what one does
  f <- mnm_fixture(4104, J = 1000L, m = 4L)
  one <- mnm_one_thread("bhhh", f)
  expect_gt(one, 7.5e6)
  old <- options(choicer.max_memory = 1.15 * one)
  on.exit(options(old), add = TRUE)
  expect_error(mnm_call("bhhh", f),
               "with 2 threads .* Use set_num_threads\\(1\\), or raise")
  set_omp_threads(1L)
  B <- mnm_call("bhhh", f)
  options(choicer.max_memory = NULL)
  expect_identical(B, mnm_call("bhhh", f))
})

test_that("an MNL call beyond the physical memory stops before allocating", {
  skip_on_cran()  # builds a 10^6-row design
  skip_if_not(test_physical_memory() > 0, "physical memory unknown")
  skip_if(test_physical_memory() > 4e12, "more than 4 TB of memory")
  old <- options(choicer.max_memory = NULL)
  on.exit(options(old), add = TRUE)
  # One choice situation offering 10^6 alternatives: the BHHH's 10^6 x 10^6
  # result alone would take 8 TB.
  J <- 1e6L
  f <- list(X = matrix(stats::rnorm(2 * J), J, 2L), alt_idx = seq_len(J),
            choice_idx = 1L, M = J, weights = 1,
            theta = c(0.1, -0.1, rep(0, J - 1L)))
  for (k in c("hessian", "bhhh")) {
    expect_error(mnm_call(k, f), "GB of physical memory", label = k)
  }
})
