# The multinomial and nested logits' loops over choice situations hand each
# thread chunks of about CHOICER_CHUNK_ROWS stacked rows (256 by default),
# at least one situation: with situations of four rows a chunk holds dozens
# of them, with wide ones a single one, and a design of a few situations is
# one chunk, which one thread of the team takes. One thread's results do
# not depend on the chunks; with more, each thread sums its own chunks, so
# the results agree with one thread's up to the order of the sums.

chl_design <- function(seed, N, m, J = 8L, K = 3L, ioo = FALSE, use_asc = TRUE,
                       blp = FALSE) {
  set.seed(seed)
  alts <- lapply(seq_len(N), function(i) sort(sample.int(J, m)))
  M <- lengths(alts)
  rows <- sum(M)
  n_delta <- if (!use_asc) 0L else if (ioo) J else J - 1L
  list(X = matrix(stats::rnorm(rows * K), rows, K),
       alt_idx = as.integer(unlist(alts)), M = as.integer(M),
       choice_idx = vapply(M, function(k) if (ioo) sample(0:k, 1L) else sample.int(k, 1L),
                           integer(1)),
       weights = stats::runif(N, 0.5, 2),
       nest_idx = as.integer((seq_len(J) - 1L) %/% 2L + 1L),
       ioo = ioo, use_asc = use_asc, blp = blp, J = J, K = K,
       theta_mnl = c(stats::rnorm(K, sd = 0.3), stats::rnorm(n_delta, sd = 0.3)),
       theta_nl = c(stats::rnorm(K, sd = 0.3), stats::runif(J %/% 2L, 0.5, 0.9),
                    stats::rnorm(n_delta, sd = 0.3)))
}

chl_kernels <- function(d) {
  ns <- asNamespace("choicer")
  o <- list(use_asc = d$use_asc, include_outside_option = d$ioo)
  mnl <- function(f, ...) do.call(ns[[f]], c(list(d$theta_mnl, d$X, d$alt_idx), list(...), o))
  nl <- function(f, ...) do.call(ns[[f]], c(list(d$theta_nl, d$X, d$alt_idx), list(...), o))
  out <- list(
    mnl_grad = mnl("mnl_loglik_gradient_parallel", d$choice_idx, d$M, d$weights,
                   opg_diag = TRUE),
    mnl_hess = mnl("mnl_loglik_hessian_parallel", d$choice_idx, d$M, d$weights),
    mnl_bhhh = mnl("mnl_bhhh_parallel", d$choice_idx, d$M, d$weights),
    mnl_scores = mnl("mnl_scores_parallel", d$choice_idx, d$M),
    mnl_pred = mnl("mnl_predict", d$M),
    mnl_shares = mnl("mnl_predict_shares", d$M, d$weights),
    mnl_elas = mnl("mnl_elasticities_parallel", d$choice_idx, d$M, d$weights, 1L),
    mnl_dr = mnl("mnl_diversion_ratios_parallel", d$M, d$weights),
    nl_grad = nl("nl_loglik_gradient_parallel", d$choice_idx, d$nest_idx, d$M, d$weights,
                 opg_diag = TRUE),
    nl_hess = nl("nl_loglik_hessian_parallel", d$choice_idx, d$nest_idx, d$M, d$weights),
    nl_bhhh = nl("nl_bhhh_parallel", d$choice_idx, d$nest_idx, d$M, d$weights),
    nl_scores = nl("nl_scores_parallel", d$choice_idx, d$nest_idx, d$M),
    nl_pred = nl("nl_predict", d$M, d$nest_idx),
    nl_shares = nl("nl_predict_shares", d$M, d$weights, d$nest_idx),
    nl_elas = nl("nl_elasticities_parallel", d$choice_idx, d$nest_idx, d$M, d$weights, 1L),
    nl_dr = nl("nl_diversion_ratios_parallel", d$nest_idx, d$M, d$weights))
  if (d$blp) {
    # the contractions reuse their shares' per-thread buffers at every iteration
    # (the step alone, whose path the sums' last bits do not move, as they can
    # move Anderson acceleration's)
    target <- rep(1 / d$J, d$J)
    out$mnl_blp <- blp_contraction(rep(0, d$J), target, d$X, d$theta_mnl[seq_len(d$K)],
                                   d$alt_idx, d$M, d$weights, tol = 1e-10,
                                   method = "simple")
    out$nl_blp <- nl_blp_contraction(rep(0, d$J), target, d$X, d$theta_nl[seq_len(d$K)],
                                     d$theta_nl[d$K + seq_len(d$J %/% 2L)], d$alt_idx,
                                     d$nest_idx, d$M, d$weights, tol = 1e-10,
                                     method = "simple")
  }
  out
}

test_that("chunked loops agree with one thread on small, wide and few situations", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  on.exit(set_omp_threads(2L), add = TRUE)
  designs <- list(
    small = chl_design(7101, N = 3000L, m = 4L, blp = TRUE),  # chunks of 64 situations
    wide = chl_design(7102, N = 200L, m = 300L, J = 400L),     # chunks of one
    few = chl_design(7103, N = 5L, m = 4L),                    # one chunk
    outside = chl_design(7104, N = 2000L, m = 3L, ioo = TRUE),
    no_constants = chl_design(7105, N = 2000L, m = 4L, use_asc = FALSE))
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
