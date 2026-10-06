# The gradient kernels' OPG diagonal: with opg_diag = TRUE they also return
# sum_u w_u s_u o s_u, the diagonal of the BHHH matrix with the same
# weights, and leave the objective and gradient as they are. Each term is
# formed and rounded as the BHHH kernel forms it, so at one thread the two
# agree to the last bit where the compiler fuses neither (clang, GCC at -O2);
# the tests allow 1e-13 relative, and 1e-12 at two threads, whose partial
# sums are added in the order the threads finish.

# Largest deviation of the diagonal d from the BHHH matrix's (or a vector's
# entries), relative above one and absolute below.
opg_dev <- function(d, B) {
  b <- if (is.matrix(B)) diag(B) else B
  max(abs(as.numeric(d) - b) / pmax(1, abs(b)))
}

# The flag-on call's objective and gradient, and its diagonal.
opg_split <- function(g) {
  d <- g$opg_diag
  g$opg_diag <- NULL
  list(g = g, d = d)
}

test_that("the mixed logit's OPG diagonal is the BHHH diagonal", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  cfgs <- list(
    list(rc_dist = c(0L, 0L), corr = TRUE, mean = FALSE, ioo = FALSE, W = "row"),
    list(rc_dist = c(0L, 1L), corr = FALSE, mean = TRUE, ioo = TRUE, W = "row"),
    list(rc_dist = c(0L, 0L), corr = TRUE, mean = TRUE, ioo = FALSE, W = "alt"))
  for (i in seq_along(cfgs)) {
    cf <- cfgs[[i]]
    fx <- mxlp_fixture(paste0("opg", i), 900L + i, rc_dist = cf$rc_dist,
                       rc_correlation = cf$corr, rc_mean = cf$mean,
                       include_outside_option = cf$ioo, W_layout = cf$W,
                       weight_type = "person")
    for (panel in c(TRUE, FALSE)) {
      f <- if (panel) fx else mxlp_cross_section(fx, 910L + i)
      for (generate in c(FALSE, TRUE)) {
        what <- sprintf("[%d, %s, %s]", i,
                        if (panel) "panel" else "cross-section",
                        if (generate) "generate" else "store")
        g0 <- mxlp_call("gradient", f, generate = generate)
        g1 <- opg_split(mxlp_call("gradient", f, generate = generate,
                                  opg_diag = TRUE))
        expect_identical(names(g0), c("objective", "gradient", "overflow"))
        expect_identical(g1$g, g0, label = paste("objective and gradient", what))
        expect_length(g1$d, length(f$theta))
        B <- mxlp_call("bhhh", f, generate = generate)
        expect_lt(opg_dev(g1$d, B), 1e-13, label = paste("diagonal", what))
      }
      # Draws in batches of two.
      g1 <- opg_split(mxlp_call("gradient", f, draw_batch = 2L, opg_diag = TRUE))
      expect_lt(opg_dev(g1$d, mxlp_call("bhhh", f, draw_batch = 2L)), 1e-13)
    }
  }
})

# The fixture with the last row of situation `t` repeated (row-level W), so
# that the alternative appears twice in the situation.
opg_dup_row <- function(fx, t) {
  r <- cumsum(fx$M)[t]
  ord <- append(seq_len(nrow(fx$X)), r, after = r)
  fx$X <- fx$X[ord, , drop = FALSE]
  fx$W <- fx$W[ord, , drop = FALSE]
  fx$alt_idx <- fx$alt_idx[ord]
  fx$M[t] <- fx$M[t] + 1L
  fx
}

test_that("the mixed logit's OPG diagonal holds without constants and with a repeated row", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  fx <- mxlp_fixture("opg_noasc", 921L, rc_correlation = TRUE,
                     use_asc = FALSE, weight_type = "person")
  fd <- mxlp_fixture("opg_dup", 922L, rc_mean = TRUE,
                     include_outside_option = TRUE, weight_type = "person")
  t_dup <- cumsum(fd$Ti)[fd$probe_units[2L]]  # the interior unit's last
  fd <- opg_dup_row(fd, t_dup)
  for (f in list(fx, mxlp_cross_section(fx, 923L), fd,
                 mxlp_cross_section(fd, 924L))) {
    g1 <- opg_split(mxlp_call("gradient", f, opg_diag = TRUE))
    expect_identical(g1$g, mxlp_call("gradient", f), label = f$name)
    expect_lt(opg_dev(g1$d, mxlp_call("bhhh", f)), 1e-13, label = f$name)
  }
})

test_that("the mixed logit's OPG diagonal is returned as computed", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  fx <- mxlp_fixture("opg_edge", 931L, rc_correlation = TRUE)
  th <- fx$theta
  th[1L] <- 4
  fx$X[5L, 1L] <- 1e308  # utilities overflow: the objective's sentinel
  g <- mxlp_call("gradient", fx, theta = th, opg_diag = TRUE)
  expect_true(g$overflow)
  expect_identical(g$objective, 1e10)
  expect_true(all(g$gradient == 0))
  expect_identical(opg_split(g)$g, mxlp_call("gradient", fx, theta = th))
  expect_true(anyNA(g$opg_diag))
  B <- mxlp_call("bhhh", fx, theta = th)
  expect_identical(is.na(as.numeric(g$opg_diag)), is.na(diag(B)))
})

# Multinomial and nested logit designs: choice sets of 2..J alternatives
# (the first situation sees all), non-uniform finite weights, and, with
# `dup`, one row repeated per situation.
opg_choice_data <- function(seed, N, J, K, ioo, dup = FALSE) {
  set.seed(seed)
  alts <- lapply(seq_len(N), function(i) {
    a <- if (i == 1L) seq_len(J) else sort(sample.int(J, sample(2:J, 1L)))
    if (dup) a <- sort(c(a, a[sample.int(length(a), 1L)]))
    a
  })
  M <- lengths(alts)
  lo <- if (ioo) 0L else 1L
  list(X = matrix(stats::rnorm(sum(M) * K), sum(M), K),
       alt_idx = as.integer(unlist(alts)), M = as.integer(M),
       choice_idx = vapply(M, function(m) sample(lo:m, 1L), integer(1)),
       weights = stats::runif(N, 0.5, 2))
}

test_that("the multinomial logit's OPG diagonal is the BHHH diagonal", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  for (ioo in c(FALSE, TRUE)) {
    for (use_asc in c(TRUE, FALSE)) {
      for (dup in c(FALSE, TRUE)) {
        d <- opg_choice_data(940L + ioo + 2L * dup, 40L, 6L, 3L, ioo, dup)
        n_delta <- if (!use_asc) 0L else if (ioo) 6L else 5L
        set.seed(950L + ioo)
        a <- list(stats::rnorm(3L + n_delta, sd = 0.4), d$X, d$alt_idx,
                  d$choice_idx, d$M, d$weights, use_asc, ioo)
        what <- sprintf("[oo = %s, asc = %s, dup = %s]", ioo, use_asc, dup)
        g0 <- do.call(mnl_loglik_gradient_parallel, a)
        g1 <- opg_split(do.call(mnl_loglik_gradient_parallel,
                                c(a, list(opg_diag = TRUE))))
        expect_identical(names(g0), c("objective", "gradient"))
        expect_identical(g1$g, g0, label = paste("objective and gradient", what))
        B <- do.call(choicer:::mnl_bhhh_parallel, a)
        expect_lt(opg_dev(g1$d, B), 1e-13, label = paste("diagonal", what))
      }
    }
  }
  # One parameter: the BHHH kernel's 1 x 1 product.
  d <- opg_choice_data(946L, 40L, 4L, 1L, FALSE)
  a <- list(0.3, d$X, d$alt_idx, d$choice_idx, d$M, d$weights, FALSE, FALSE)
  g1 <- do.call(mnl_loglik_gradient_parallel, c(a, list(opg_diag = TRUE)))
  expect_lt(opg_dev(g1$opg_diag, do.call(choicer:::mnl_bhhh_parallel, a)), 1e-13)
})

test_that("the nested logit's OPG diagonal is the BHHH diagonal", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  nest <- c(1L, 1L, 2L, 2L, 2L, 3L, 4L, 4L)  # nest 3 is a singleton
  for (ioo in c(FALSE, TRUE)) {
    for (use_asc in c(TRUE, FALSE)) {
      for (dup in c(FALSE, TRUE)) {
        d <- opg_choice_data(960L + ioo + 2L * dup, 50L, 8L, 2L, ioo, dup)
        n_delta <- if (!use_asc) 0L else if (ioo) 8L else 7L
        set.seed(970L + ioo)
        theta <- c(stats::rnorm(2L, sd = 0.4), stats::runif(3L, 0.55, 1),
                   stats::rnorm(n_delta, sd = 0.4))
        a <- list(theta, d$X, d$alt_idx, d$choice_idx, nest, d$M, d$weights,
                  use_asc, ioo)
        what <- sprintf("[oo = %s, asc = %s, dup = %s]", ioo, use_asc, dup)
        g0 <- do.call(nl_loglik_gradient_parallel, a)
        g1 <- opg_split(do.call(nl_loglik_gradient_parallel,
                                c(a, list(opg_diag = TRUE))))
        expect_identical(g1$g, g0, label = paste("objective and gradient", what))
        B <- do.call(choicer:::nl_bhhh_parallel, a)
        expect_lt(opg_dev(g1$d, B), 1e-13, label = paste("diagonal", what))
      }
    }
  }
})

test_that("the multinomial and nested logits' OPG diagonals are returned as computed", {
  mxlp_threads(1L)
  on.exit(mxlp_threads(2L), add = TRUE)
  d <- opg_choice_data(948L, 30L, 4L, 2L, FALSE)
  d$X[5L, 1L] <- 1e308  # utilities overflow
  th <- c(4, 0.2, 0.1, -0.1, 0.3)
  a <- list(th, d$X, d$alt_idx, d$choice_idx, d$M, d$weights, TRUE, FALSE)
  msg <- "Non-finite log-probability"
  expect_warning(g1 <- opg_split(do.call(mnl_loglik_gradient_parallel,
                                         c(a, list(opg_diag = TRUE)))), msg)
  expect_warning(g0 <- do.call(mnl_loglik_gradient_parallel, a), msg)
  expect_identical(g1$g, g0)
  B <- do.call(choicer:::mnl_bhhh_parallel, a)
  expect_true(anyNA(g1$d))
  expect_identical(is.na(as.numeric(g1$d)), is.na(diag(B)))
  ok <- !is.na(diag(B))
  expect_lt(opg_dev(g1$d[ok], diag(B)[ok]), 1e-13)
  nest <- c(1L, 1L, 2L, 2L)
  a <- list(c(4, 0.2, 0.7, 0.8, 0.1, -0.1, 0.3), d$X, d$alt_idx, d$choice_idx,
            nest, d$M, d$weights, TRUE, FALSE)
  g1 <- opg_split(do.call(nl_loglik_gradient_parallel,
                          c(a, list(opg_diag = TRUE))))
  expect_identical(g1$g, do.call(nl_loglik_gradient_parallel, a))
  B <- do.call(choicer:::nl_bhhh_parallel, a)
  expect_true(anyNA(g1$d))
  expect_identical(is.na(as.numeric(g1$d)), is.na(diag(B)))
  ok <- !is.na(diag(B))
  expect_lt(opg_dev(g1$d[ok], diag(B)[ok]), 1e-13)
})

test_that("the OPG diagonal is the BHHH diagonal at two threads", {
  mxlp_threads(2L)
  # About 1,500 situations, so that both threads take part.
  fx <- mxlp_fixture("opg_threads", 981L, rc_correlation = TRUE,
                     include_outside_option = TRUE, weight_type = "person",
                     U = 600L)
  for (f in list(fx, mxlp_cross_section(fx, 982L))) {
    g <- mxlp_call("gradient", f, opg_diag = TRUE)
    expect_lt(opg_dev(g$opg_diag, mxlp_call("bhhh", f)), 1e-12)
  }
  for (ioo in c(FALSE, TRUE)) {
    d <- opg_choice_data(983L + ioo, 1500L, 6L, 3L, ioo, dup = TRUE)
    set.seed(985L + ioo)
    n_delta <- if (ioo) 6L else 5L
    a <- list(stats::rnorm(3L + n_delta, sd = 0.4), d$X, d$alt_idx,
              d$choice_idx, d$M, d$weights, TRUE, ioo)
    g <- do.call(mnl_loglik_gradient_parallel, c(a, list(opg_diag = TRUE)))
    expect_lt(opg_dev(g$opg_diag, do.call(choicer:::mnl_bhhh_parallel, a)),
              1e-12)
    theta <- c(a[[1L]][1:3], 0.6, 0.8, a[[1L]][-(1:3)])
    a <- list(theta, d$X, d$alt_idx, d$choice_idx, c(1L, 1L, 1L, 2L, 2L, 2L),
              d$M, d$weights, TRUE, ioo)
    g <- do.call(nl_loglik_gradient_parallel, c(a, list(opg_diag = TRUE)))
    expect_lt(opg_dev(g$opg_diag, do.call(choicer:::nl_bhhh_parallel, a)),
              1e-12)
  }
})
