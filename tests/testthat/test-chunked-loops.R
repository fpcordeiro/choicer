# The multinomial and nested logits' loops over choice situations hand each
# thread chunks of about CHOICER_CHUNK_ROWS stacked rows (256 by default),
# at least one situation: with situations of four rows a chunk holds dozens
# of them, with wide ones a single one, and a design of a few situations is
# one chunk, which one thread of the team takes. One thread's results do
# not depend on the chunks; with more, each thread sums its own chunks, so
# the results agree with one thread's up to the order of the sums.

chl_design <- function(seed, N, m, J = 8L, K = 3L, nests = FALSE) {
  set.seed(seed)
  alts <- lapply(seq_len(N), function(i) sort(sample.int(J, m)))
  M <- lengths(alts)
  rows <- sum(M)
  list(X = matrix(stats::rnorm(rows * K), rows, K),
       alt_idx = as.integer(unlist(alts)), M = as.integer(M),
       choice_idx = vapply(M, function(k) sample.int(k, 1L), integer(1)),
       weights = stats::runif(N, 0.5, 2),
       nest_idx = as.integer((seq_len(J) - 1L) %/% 2L + 1L),
       theta_mnl = c(stats::rnorm(K, sd = 0.3), stats::rnorm(J - 1L, sd = 0.3)),
       theta_nl = c(stats::rnorm(K, sd = 0.3), stats::runif(J %/% 2L, 0.5, 0.9),
                    stats::rnorm(J - 1L, sd = 0.3)))
}

chl_kernels <- function(d) {
  ns <- asNamespace("choicer")
  list(
    mnl_grad = ns$mnl_loglik_gradient_parallel(d$theta_mnl, d$X, d$alt_idx, d$choice_idx,
                                               d$M, d$weights, opg_diag = TRUE),
    mnl_hess = ns$mnl_loglik_hessian_parallel(d$theta_mnl, d$X, d$alt_idx, d$choice_idx,
                                              d$M, d$weights),
    mnl_bhhh = ns$mnl_bhhh_parallel(d$theta_mnl, d$X, d$alt_idx, d$choice_idx, d$M, d$weights),
    mnl_scores = ns$mnl_scores_parallel(d$theta_mnl, d$X, d$alt_idx, d$choice_idx, d$M),
    mnl_pred = ns$mnl_predict(d$theta_mnl, d$X, d$alt_idx, d$M),
    mnl_shares = ns$mnl_predict_shares(d$theta_mnl, d$X, d$alt_idx, d$M, d$weights),
    mnl_elas = ns$mnl_elasticities_parallel(d$theta_mnl, d$X, d$alt_idx, d$choice_idx, d$M,
                                            d$weights, 1L),
    mnl_dr = ns$mnl_diversion_ratios_parallel(d$theta_mnl, d$X, d$alt_idx, d$M, d$weights),
    nl_grad = ns$nl_loglik_gradient_parallel(d$theta_nl, d$X, d$alt_idx, d$choice_idx,
                                             d$nest_idx, d$M, d$weights, opg_diag = TRUE),
    nl_hess = ns$nl_loglik_hessian_parallel(d$theta_nl, d$X, d$alt_idx, d$choice_idx,
                                            d$nest_idx, d$M, d$weights),
    nl_bhhh = ns$nl_bhhh_parallel(d$theta_nl, d$X, d$alt_idx, d$choice_idx, d$nest_idx,
                                  d$M, d$weights),
    nl_scores = ns$nl_scores_parallel(d$theta_nl, d$X, d$alt_idx, d$choice_idx,
                                      d$nest_idx, d$M),
    nl_pred = ns$nl_predict(d$theta_nl, d$X, d$alt_idx, d$M, d$nest_idx),
    nl_shares = ns$nl_predict_shares(d$theta_nl, d$X, d$alt_idx, d$M, d$weights,
                                     d$nest_idx),
    nl_elas = ns$nl_elasticities_parallel(d$theta_nl, d$X, d$alt_idx, d$choice_idx,
                                          d$nest_idx, d$M, d$weights, 1L),
    nl_dr = ns$nl_diversion_ratios_parallel(d$theta_nl, d$X, d$alt_idx, d$nest_idx, d$M,
                                            d$weights))
}

test_that("chunked loops agree with one thread on small, wide and few situations", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  on.exit(set_omp_threads(2L), add = TRUE)
  designs <- list(
    small = chl_design(7101, N = 3000L, m = 4L),           # chunks of 64 situations
    wide = chl_design(7102, N = 200L, m = 300L, J = 400L),  # chunks of one
    few = chl_design(7103, N = 5L, m = 4L))                # one chunk
  for (nm in names(designs)) {
    set_omp_threads(1L)
    one <- chl_kernels(designs[[nm]])
    set_omp_threads(2L)
    two <- chl_kernels(designs[[nm]])
    for (k in names(one)) {
      expect_equal(two[[k]], one[[k]], tolerance = 1e-12, label = paste(nm, k))
    }
  }
})
