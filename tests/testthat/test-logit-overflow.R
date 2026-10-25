# The multinomial and nested logit gradient kernels' overflow flag: TRUE where
# the objective or the gradient is not the negated log-likelihood's at theta
# (a situation's log-probability clamped to -1e10, the total replaced by the
# 1e10 sentinel, a gradient entry set to 0), FALSE elsewhere. run_polish()
# reads it to fail a trial point; the optimizers see the list without it.

# N situations over J alternatives in sets of 3 to 5 (`lo`..`hi`), nests
# {1, 2, 3}, {4, 5} and {6} for J = 6, else J / 2 nests of two; with `ioo`,
# the outside option can be chosen.
lo_data <- function(seed, ioo = FALSE, N = 40L, J = 6L, lo = 3L, hi = 5L) {
  set.seed(seed)
  M <- sample(lo:hi, N, replace = TRUE)
  alt <- unlist(lapply(M, function(m) sort(sample.int(J, m))))
  list(X = matrix(stats::rnorm(sum(M) * 2L), ncol = 2L), alt_idx = as.integer(alt),
       choice_idx = vapply(M, function(m) sample.int(m + ioo, 1L) - ioo, 1L),
       M = as.integer(M), weights = stats::runif(N, 0.5, 2),
       nest_idx = if (J == 6L) c(1L, 1L, 1L, 2L, 2L, 3L) else rep(seq_len(J / 2L), each = 2L),
       ioo = ioo, J = J)
}

lo_theta <- function(model, d, seed = 1L) {
  set.seed(seed)
  K_l <- sum(tabulate(d$nest_idx) > 1L)
  c(stats::rnorm(2L, sd = 0.4), if (model == "nl") seq(0.5, 0.8, length.out = K_l),
    stats::rnorm(if (d$ioo) d$J else d$J - 1L, sd = 0.4))
}

lo_grad <- function(model, theta, d, ...) {
  if (model == "mnl") {
    mnl_loglik_gradient_parallel(theta, d$X, d$alt_idx, d$choice_idx, d$M,
                                 d$weights, TRUE, d$ioo, ...)
  } else {
    nl_loglik_gradient_parallel(theta, d$X, d$alt_idx, d$choice_idx,
                                d$nest_idx, d$M, d$weights, TRUE, d$ioo, ...)
  }
}

# A row of situation s that is not chosen
lo_other_row <- function(d, s) {
  c(0L, cumsum(d$M))[s] + (if (d$choice_idx[s] == 1L) 2L else 1L)
}

# Situations `ss` made to clamp alone: two rows, alternative 4 (nest 2) and
# the chosen alternative 6 (row 2, alone in nest 3), utilities +0.9e308 and
# -0.9e308 through x1 (1e300 times +-0.9e8; x1 is 0 elsewhere). The chosen
# log-probability is -0.9e308 - 0.9e308 = -Inf, while every probability,
# score and the rest of the objective stay finite (in NL while nest 2's
# dissimilarity exceeds 0.9e308 / DBL_MAX, about 0.5: lo_theta() gives 0.8).
lo_clamp <- function(d, ss) {
  rows <- split(seq_along(d$alt_idx), rep(seq_along(d$M), d$M))
  keep <- unlist(lapply(seq_along(d$M), function(s) {
    if (s %in% ss) rows[[s]][1:2] else rows[[s]]
  }))
  X <- d$X[keep, , drop = FALSE]; alt <- d$alt_idx[keep]
  M <- d$M; M[ss] <- 2L
  start <- c(0L, cumsum(M))
  for (s in ss) {
    r <- start[s] + 1:2
    alt[r] <- c(4L, 6L)
    X[r, 1L] <- c(0.9e8, -0.9e8)
  }
  X[-unlist(lapply(ss, function(s) start[s] + 1:2)), 1L] <- 0
  d$X <- X; d$alt_idx <- as.integer(alt); d$M <- as.integer(M)
  d$choice_idx[ss] <- 2L
  d
}

test_that("the flag is FALSE where the objective is the log-likelihood", {
  for (model in c("mnl", "nl")) for (ioo in c(FALSE, TRUE)) {
    d <- lo_data(101L + ioo, ioo)
    for (od in c(FALSE, TRUE)) {
      r <- lo_grad(model, lo_theta(model, d), d, opg_diag = od)
      lab <- sprintf("%s oo = %s opg = %s", model, ioo, od)
      expect_named(r, c("objective", "gradient", "overflow",
                        if (od) "opg_diag"), label = lab)
      expect_identical(r$overflow, FALSE, label = lab)
    }
  }
})

test_that("a clamped log-probability alone sets the flag", {
  # finite utilities, probabilities and scores: only the chosen
  # log-probability is not finite, and only the clamp marks it
  for (model in c("mnl", "nl")) for (ioo in c(FALSE, TRUE)) for (od in c(FALSE, TRUE)) {
    d <- lo_clamp(lo_data(111L + ioo, ioo), 2L)
    th <- lo_theta(model, d)
    th[1L] <- 1e300
    lab <- sprintf("%s oo = %s opg = %s", model, ioo, od)
    r <- NULL
    if (model == "mnl") {
      expect_warning(r <- lo_grad(model, th, d, opg_diag = od), "clamped to -1e10")
    } else {
      expect_no_warning(r <- lo_grad(model, th, d, opg_diag = od))
    }
    expect_true(r$overflow, label = lab)
    expect_true(all(is.finite(r$gradient)) && all(r$gradient != 0), label = lab)
    # the clamped situation contributes its weight times 1e10, the others
    # their (finite) terms: not the sentinel
    expect_gt(r$objective, 1e10 * d$weights[2L] * (1 - 1e-12))
    expect_lt(r$objective, 1e10 * d$weights[2L] + 1e6)
  }
})

test_that("the 1e10 sentinel sets the flag", {
  # an infinite weight makes the total infinite: the objective is the
  # sentinel and the gradient zero
  for (model in c("mnl", "nl")) for (ioo in c(FALSE, TRUE)) {
    d <- lo_data(116L + ioo, ioo)
    d$weights[1L] <- Inf
    r <- lo_grad(model, lo_theta(model, d), d)
    lab <- sprintf("%s oo = %s", model, ioo)
    expect_identical(r$objective, 1e10, label = lab)
    expect_true(all(r$gradient == 0), label = lab)
    expect_true(r$overflow, label = lab)
  }
})

test_that("scores whose sum overflows set the flag, the objective finite", {
  # an unchosen row of each of ten situations at 1e308 times 1e-308: finite
  # utilities, but the coefficient's score sums past the largest double
  for (model in c("mnl", "nl")) for (ioo in c(FALSE, TRUE)) {
    d <- lo_data(121L + ioo, ioo)
    th <- lo_theta(model, d)
    r0 <- lo_grad(model, th, d)
    for (s in 1:10) d$X[lo_other_row(d, s), 2L] <- 1e308
    th[2L] <- 1e-308
    r <- lo_grad(model, th, d)
    lab <- sprintf("%s oo = %s", model, ioo)
    expect_true(is.finite(r$objective) && r$objective < 1e10, label = lab)
    expect_identical(r$gradient[2L], 0, label = lab)
    expect_true(r$overflow, label = lab)
    expect_false(r0$overflow, label = lab)
  }
})

test_that("a nest whose inclusive value overflows sets the flag", {
  # An unchosen row of nest 1 (dissimilarity 0.5) at a utility of 1e308:
  # finite, but its utility over lambda is Inf, so the nest's log I_k is
  # -Inf and the nest leaves the denominator. The chosen log-probability
  # stays finite and too high (above the one at 0.4e308, where the nest
  # dominates), while the nest's probabilities, NaN, reach the gradient,
  # which the kernel sets to 0: only the flag tells. With 3 nests over sets
  # of 3 to 5 the kernel takes every nest; with 10 nests of two over sets of
  # 2 to 3, the nests each situation offers.
  for (J in c(6L, 20L)) for (ioo in c(FALSE, TRUE)) {
    d <- if (J == 6L) lo_data(131L + ioo, ioo) else lo_data(133L + ioo, ioo, J = 20L, lo = 2L, hi = 3L)
    th <- lo_theta("nl", d)
    th[3L] <- 0.5
    th[1L] <- 1
    d$X[, 1L] <- 0
    lab <- sprintf("J = %d, oo = %s", J, ioo)
    expect_false(lo_grad("nl", th, d)$overflow, label = lab)
    # a situation that offers nest 1 and chooses outside it
    start <- c(0L, cumsum(d$M))
    nest1 <- d$nest_idx[d$alt_idx] == 1L
    chosen_row <- start[seq_along(d$M)] + d$choice_idx
    s <- which(vapply(seq_along(d$M), function(i) {
      rows <- (start[i] + 1L):start[i + 1L]
      any(nest1[rows]) && (d$choice_idx[i] == 0L || !nest1[chosen_row[i]])
    }, NA))[1L]
    expect_false(is.na(s), label = lab)
    rows <- (start[s] + 1L):start[s + 1L]
    r_up <- rows[nest1[rows]][1L]
    d_inf <- d; d_inf$X[r_up, 1L] <- 1e308
    d_big <- d; d_big$X[r_up, 1L] <- 0.4e308
    r <- lo_grad("nl", th, d_inf)
    expect_true(r$overflow, label = paste("dropped nest,", lab))
    expect_lt(r$objective, 1e6)
    expect_identical(r$gradient[1:2], c(0, 0), label = lab)
    expect_lt(r$objective, lo_grad("nl", th, d_big)$objective)
    # a NaN dissimilarity takes the kernel to the path over every nest,
    # where it clamps every situation
    th_nan <- th
    th_nan[3L] <- NaN
    expect_true(lo_grad("nl", th_nan, d)$overflow, label = paste("NaN,", lab))
  }
})

test_that("the flag does not depend on the number of threads", {
  # clamped situations spread over the data (the dynamic schedule decides
  # which thread meets them): a flag the reduction loses would show
  on.exit(set_omp_threads(2L), add = TRUE)
  for (model in c("mnl", "nl")) {
    d <- lo_data(141L, TRUE, N = 400L)
    th <- lo_theta(model, d)
    d2 <- lo_clamp(d, c(150L, 390L))
    th2 <- th
    th2[1L] <- 1e300
    flags <- vapply(1:2, function(t) {
      set_omp_threads(t)
      c(lo_grad(model, th, d)$overflow,
        suppressWarnings(lo_grad(model, th2, d2))$overflow)
    }, logical(2))
    expect_identical(flags[, 1L], c(FALSE, TRUE), label = model)
    expect_identical(flags[, 2L], flags[, 1L], label = model)
  }
})

test_that("the optimizers see the objective and gradient alone", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  set.seed(151)
  N <- 60L
  dt <- data.table::data.table(id = rep(seq_len(N), each = 4L), alt = rep(1:4, N))
  dt[, `:=`(x1 = stats::rnorm(.N), x2 = stats::rnorm(.N))]
  dt[, choice := as.integer(alt == sample.int(4L, 1L)), by = id]
  dt[, nest := c(1L, 1L, 2L, 2L)[alt]]
  seen <- NULL
  local_mocked_bindings(run_optimizer = function(optimizer, theta_init, eval_f,
                                                 lower = NULL, upper = NULL,
                                                 control = list()) {
    seen <<- eval_f(theta_init)
    stop("captured")
  })
  for (sv in c("none", "bhhh")) {
    seen <- NULL
    expect_error(run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"),
                             scale_vars = sv), "captured")
    expect_named(seen, c("objective", "gradient"))
    seen <- NULL
    expect_error(run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest",
                               scale_vars = sv), "captured")
    expect_named(seen, c("objective", "gradient"))
  }
})

test_that("a custom optimizer gets the objective and gradient alone", {
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  set.seed(152)
  N <- 60L
  dt <- data.table::data.table(id = rep(seq_len(N), each = 4L), alt = rep(1:4, N))
  dt[, `:=`(x1 = stats::rnorm(.N), x2 = stats::rnorm(.N))]
  dt[, choice := as.integer(alt == sample.int(4L, 1L)), by = id]
  dt[, nest := c(1L, 1L, 2L, 2L)[alt]]
  nms <- NULL
  opt <- function(theta_init, eval_f, lower, upper) {
    r <- eval_f(theta_init)
    nms <<- names(r)
    list(par = theta_init, value = r$objective)
  }
  suppressMessages(run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"),
                               optimizer = opt))
  expect_identical(nms, c("objective", "gradient"))
  nms <- NULL
  suppressMessages(run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest",
                                 optimizer = opt))
  expect_identical(nms, c("objective", "gradient"))
})
