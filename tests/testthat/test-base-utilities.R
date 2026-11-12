# The multinomial and nested logit kernels form each choice situation's base
# utilities, X_i beta plus its rows' constants, inside their loops over the
# situations, on X in place (BaseUtil, src/choicer_internal.h): sixteen rows
# at a time, then eight, four, two and one, the rows that whole blocks leave
# over as one more block ending at the last row. mnl_predict() and
# nl_predict() return those utilities, so they show every path: situations of
# every size from 1 to 40 rows and a few larger ones, against X %*% beta plus
# the constants, and the same rows cut into other situations or taken alone
# in designs of one row, which must not move a bit (each row's sum is its
# own, whatever block forms it).

# Rows stacked into situations of sizes `M`, K covariates on mixed scales, J
# alternatives in nests of four (codes drawn with replacement where a
# situation has more rows than J; the utilities do not care).
bu_design <- function(seed, M, K, J = 12L) {
  stopifnot(J %% 4L == 0L)
  set.seed(seed)
  rows <- sum(M)
  X <- matrix(stats::rnorm(rows * K), rows, K) *
    rep(10^seq(-2, 2, length.out = K), each = rows)
  list(X = X, M = as.integer(M), J = J, K = K,
       alt_idx = as.integer(unlist(lapply(M, function(m) {
         sample.int(J, m, replace = m > J)
       }))),
       beta = stats::rnorm(K), delta = stats::rnorm(J),
       nest_idx = as.integer((seq_len(J) - 1L) %/% 4L + 1L),
       lambda = stats::runif(J %/% 4L, 0.4, 0.9))
}

# The utilities the two kernels return, for constants or none and an outside
# option or none, with the reference: X beta, plus each row's constant (the
# first alternative's is 0 without an outside option).
bu_utilities <- function(d, use_asc, ioo) {
  ns <- asNamespace("choicer")
  delta_free <- if (!use_asc) numeric(0) else if (ioo) d$delta else d$delta[-1L]
  delta_all <- if (!use_asc) rep(0, d$J) else if (ioo) d$delta else c(0, d$delta[-1L])
  mnl <- ns$mnl_predict(c(d$beta, delta_free), d$X, d$alt_idx, d$M,
                        use_asc = use_asc, include_outside_option = ioo)
  nl <- ns$nl_predict(c(d$beta, d$lambda, delta_free), d$X, d$alt_idx, d$M,
                      d$nest_idx, use_asc = use_asc, include_outside_option = ioo)
  list(mnl = drop(mnl$utility), nl = drop(nl$utility),
       ref = drop(d$X %*% d$beta) + delta_all[d$alt_idx],
       p_mnl = drop(mnl$choice_prob), p_nl = drop(nl$choice_prob))
}

test_that("the kernels' utilities are X beta plus the constants for every situation size", {
  sizes <- c(1:40, 47L, 48L, 49L, 63L, 64L, 65L, 100L, 3L, 1L, 17L)
  for (K in c(1L, 2L, 3L, 7L)) {
    d <- bu_design(8100 + K, sizes, K)
    for (use_asc in c(TRUE, FALSE)) {
      for (ioo in c(FALSE, TRUE)) {
        u <- bu_utilities(d, use_asc, ioo)
        lab <- sprintf("K = %d, use_asc = %s, outside = %s", K, use_asc, ioo)
        expect_equal(u$mnl, u$ref, tolerance = 1e-13, label = paste("MNL,", lab))
        # one helper serves both models
        expect_identical(u$nl, u$mnl, label = paste("NL,", lab))
        # and the probabilities are the softmax of those utilities
        if (!ioo) {
          sit <- rep(seq_along(d$M), d$M)
          p_ref <- unlist(lapply(split(u$mnl, sit), function(v) {
            e <- exp(v - max(v))
            e / sum(e)
          }), use.names = FALSE)
          expect_equal(u$p_mnl, p_ref, tolerance = 1e-12, label = paste("MNL P,", lab))
        }
      }
    }
  }
})

test_that("a row's utility does not depend on the situation it is stacked in", {
  rows <- 211L
  cuts <- list(
    ones = rep(1L, rows),
    whole = rows,
    twos = c(rep(2L, 105L), 1L),
    threes = c(1L, rep(3L, 70L)),
    fives = c(rep(5L, 42L), 1L),
    sevens = c(4L, rep(7L, 29L), 4L),
    nines = c(rep(9L, 23L), 4L),
    mixed = c(16L, 8L, 4L, 2L, 1L, 31L, 33L, 15L, 17L, 23L, 40L, 21L))
  expect_true(all(vapply(cuts, sum, integer(1)) == rows))
  for (K in c(1L, 3L, 10L)) {
    d <- bu_design(8200 + K, rows, K, J = 40L)
    ref <- list()  # the first cut's utilities, with and without an outside option
    for (nm in names(cuts)) {
      d$M <- as.integer(cuts[[nm]])
      for (ioo in c(FALSE, TRUE)) {
        u <- bu_utilities(d, use_asc = TRUE, ioo = ioo)
        key <- paste(K, ioo)
        if (is.null(ref[[key]])) ref[[key]] <- u$mnl
        expect_identical(u$mnl, ref[[key]],
                         label = sprintf("MNL, K = %d, outside = %s, cut %s", K, ioo, nm))
        expect_identical(u$nl, ref[[key]],
                         label = sprintf("NL, K = %d, outside = %s, cut %s", K, ioo, nm))
      }
    }
    # ... and each row alone, a design of one row, whose columns lie one value
    # apart (a single BLAS call over such a design could sum it otherwise)
    for (ioo in c(FALSE, TRUE)) {
      alone <- vapply(seq_len(rows), function(i) {
        di <- d
        di$X <- d$X[i, , drop = FALSE]
        di$alt_idx <- d$alt_idx[i]
        di$M <- 1L
        u <- bu_utilities(di, use_asc = TRUE, ioo = ioo)
        c(u$mnl, u$nl)
      }, numeric(2))
      lab <- sprintf("K = %d, outside = %s, designs of one row", K, ioo)
      expect_identical(alone[1L, ], ref[[paste(K, ioo)]], label = paste("MNL,", lab))
      expect_identical(alone[2L, ], ref[[paste(K, ioo)]], label = paste("NL,", lab))
    }
  }
})

test_that("the utilities do not depend on the thread count", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  on.exit(set_omp_threads(2L), add = TRUE)
  d <- bu_design(8301, rep(c(3L, 9L, 20L, 37L), 150L), 4L)
  set_omp_threads(1L)
  one <- bu_utilities(d, TRUE, FALSE)
  set_omp_threads(2L)
  two <- bu_utilities(d, TRUE, FALSE)
  expect_identical(two$mnl, one$mnl)
  expect_identical(two$nl, one$nl)
  # each situation's probabilities come from its own utilities alone
  expect_identical(two$p_mnl, one$p_mnl)
  expect_identical(two$p_nl, one$p_nl)
})

test_that("a zero coefficient and non-finite coefficients pass through the sums", {
  d <- bu_design(8401, c(1L, 2L, 5L, 8L, 13L, 16L, 21L), 3L)
  ns <- asNamespace("choicer")
  theta <- c(d$beta, d$delta[-1L])
  theta[2L] <- 0
  u0 <- drop(ns$mnl_predict(theta, d$X, d$alt_idx, d$M)$utility)
  expect_equal(u0, drop(d$X[, -2L] %*% d$beta[-2L]) + c(0, d$delta[-1L])[d$alt_idx],
               tolerance = 1e-13)
  # all coefficients zero: each row's utility is its constant, exactly
  theta[seq_len(3L)] <- 0
  expect_identical(drop(ns$mnl_predict(theta, d$X, d$alt_idx, d$M)$utility),
                   c(0, d$delta[-1L])[d$alt_idx])
  theta[1L] <- NaN
  expect_true(all(is.nan(ns$mnl_predict(theta, d$X, d$alt_idx, d$M)$utility)))
  theta[1L] <- Inf
  u_inf <- drop(ns$mnl_predict(theta, d$X, d$alt_idx, d$M)$utility)
  expect_identical(u_inf, ifelse(d$X[, 1L] > 0, Inf, -Inf))
})
