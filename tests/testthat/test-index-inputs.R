# Index inputs of the multinomial and nested logit kernels.
#
# The kernels read each row's alternative code, each situation's choice, the
# alternatives' nests and the situation sizes in place from integer vectors
# (a double vector is coerced once per call, by Rcpp), validate them before
# any parallel region, and never modify them. The checks run in a fixed
# order, so a call with two broken inputs reports the first of them.

# A small design with unbalanced choice sets over J alternatives; situation 1
# offers every alternative.
idx_fixture <- function(seed, N = 12L, J = 4L, K = 2L, ioo = FALSE,
                        use_asc = TRUE) {
  set.seed(seed)
  sets <- c(list(seq_len(J)),
            lapply(seq_len(N - 1L),
                   function(i) sort(sample.int(J, sample(2:J, 1L)))))
  M <- lengths(sets)
  n <- sum(M)
  lo <- if (ioo) 0L else 1L
  n_delta <- if (!use_asc) 0L else if (ioo) J else J - 1L
  f <- list(theta = rnorm(K + n_delta, sd = 0.3),
            X = matrix(rnorm(n * K), n, K), alt_idx = unlist(sets),
            choice_idx = vapply(M, function(m) sample(lo:m, 1L), integer(1)),
            M = M, weights = runif(N, 0.5, 2), J = J, K = K, N = N,
            ioo = ioo, use_asc = use_asc)
  # BLP target: the shares the model predicts with ASCs of 0.1, which the
  # contraction recovers in a few iterations
  f$target <- as.numeric(mnl_predict_shares(
    c(f$theta[seq_len(K)], rep(0.1, if (ioo) J else J - 1L)), f$X, f$alt_idx,
    f$M, f$weights, TRUE, ioo))
  f
}

mnl_kernels <- c("gradient", "bhhh", "scores", "predict", "shares", "blp",
                 "hessian", "elasticities", "diversion")

mnl_call <- function(k, f) {
  switch(k,
    gradient = mnl_loglik_gradient_parallel(f$theta, f$X, f$alt_idx,
      f$choice_idx, f$M, f$weights, f$use_asc, f$ioo),
    bhhh = mnl_bhhh_parallel(f$theta, f$X, f$alt_idx, f$choice_idx, f$M,
      f$weights, f$use_asc, f$ioo),
    scores = mnl_scores_parallel(f$theta, f$X, f$alt_idx, f$choice_idx, f$M,
      f$use_asc, f$ioo),
    predict = mnl_predict(f$theta, f$X, f$alt_idx, f$M, f$use_asc, f$ioo),
    shares = mnl_predict_shares(f$theta, f$X, f$alt_idx, f$M, f$weights,
      f$use_asc, f$ioo),
    blp = blp_contraction(rep(0, f$J), f$target, f$X, f$theta[seq_len(f$K)],
      f$alt_idx, f$M, f$weights, f$ioo),
    hessian = mnl_loglik_hessian_parallel(f$theta, f$X, f$alt_idx,
      f$choice_idx, f$M, f$weights, f$use_asc, f$ioo),
    elasticities = mnl_elasticities_parallel(f$theta, f$X, f$alt_idx,
      f$choice_idx, f$M, f$weights, 1L, f$use_asc, f$ioo),
    diversion = mnl_diversion_ratios_parallel(f$theta, f$X, f$alt_idx, f$M,
      f$weights, f$use_asc, f$ioo),
    stop("unknown kernel: ", k))
}

idx_err <- function(k, f) {
  warns <- character()
  msg <- tryCatch(
    withCallingHandlers({
      mnl_call(k, f)
      "<no error>"
    }, warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = conditionMessage)
  c(msg, warns)
}

with_input <- function(f, ...) {
  a <- list(...)
  for (nm in names(a)) f[nm] <- list(a[[nm]])
  f
}

test_that("the MNL kernels give identical results for double-typed indices", {
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L) # one thread: two calls sum in the same order
  for (ioo in c(FALSE, TRUE)) {
    for (use_asc in c(TRUE, FALSE)) {
      f <- idx_fixture(11L + ioo, ioo = ioo, use_asc = use_asc)
      fd <- with_input(f, alt_idx = as.double(f$alt_idx),
                       choice_idx = as.double(f$choice_idx),
                       M = as.double(f$M))
      for (k in mnl_kernels) {
        expect_identical(mnl_call(k, fd), mnl_call(k, f), label = sprintf(
          "%s (outside option %s, ASCs %s) with double-typed indices",
          k, ioo, use_asc))
      }
    }
  }
})

test_that("MNL input errors give the same message in every kernel", {
  # One broken input per case, in the checks' order; the last cases break two
  # checks at once and pin which one reports first.
  f <- idx_fixture(21L)                 # no outside option, ASCs
  fo <- idx_fixture(22L, ioo = TRUE)    # outside option
  n <- sum(f$M)
  N <- f$N
  J <- f$J
  choice_kernels <- c("gradient", "bhhh", "scores")
  weight_kernels <- c("gradient", "bhhh", "shares", "blp", "hessian",
                      "elasticities", "diversion")
  suffix <- c(bhhh = " (mnl_bhhh_parallel)", scores = " (mnl_scores_parallel)")
  bad_choice <- function(k) {
    if (!k %in% choice_kernels) return("<no error>")
    paste0("Invalid chosen alternative index for individual 3",
           if (k %in% names(suffix)) suffix[[k]] else "")
  }
  code_msg <- function(found) {
    sprintf("alt_idx must use 1-based alternative indices (found %s).", found)
  }
  cases <- list(
    list("M = 0", with_input(f, M = replace(f$M, 3L, 0L)),
         "M must be positive for every individual (M[3] = 0)."),
    list("M = NA", with_input(f, M = replace(f$M, 3L, NA_integer_)),
         "M must be positive for every individual (M[3] = -2147483648)."),
    list("X one row short", with_input(f, X = f$X[-1L, , drop = FALSE]),
         sprintf("X has %d rows but sum(M) is %d.", n - 1L, n)),
    list("alt_idx one short", with_input(f, alt_idx = f$alt_idx[-1L]),
         sprintf(paste("alt_idx length (%d) does not match the number of",
                       "rows of X (%d)."), n - 1L, n)),
    list("weights one short", with_input(f, weights = f$weights[-1L]),
         function(k) if (k %in% weight_kernels) {
           sprintf("weights length (%d) does not match N (%d)", N - 1L, N)
         } else "<no error>"),
    list("choice_idx one short", with_input(f, choice_idx = f$choice_idx[-1L]),
         function(k) if (k %in% choice_kernels) {
           sprintf("choice_idx length (%d) does not match N (%d)", N - 1L, N)
         } else "<no error>"),
    list("alternative code 0", with_input(f, alt_idx = replace(f$alt_idx, 5L, 0L)),
         code_msg(0)),
    list("alternative code NA",
         with_input(f, alt_idx = replace(f$alt_idx, 5L, NA_integer_)),
         code_msg("NA")),
    list("alternative code -1",
         with_input(f, alt_idx = replace(f$alt_idx, 5L, -1L)), code_msg(-1)),
    list("alternative code 3e9 (a double beyond the integer range)",
         with_input(f, alt_idx = replace(as.double(f$alt_idx), 5L, 3e9)),
         c(code_msg("NA"), "NAs introduced by coercion to integer range")),
    list("alternative code past the delta block",
         with_input(f, alt_idx = replace(f$alt_idx, 5L, J + 1L)),
         sprintf(paste("Theta's delta (ASC) block implies %d alternatives but",
                       "alt_idx references alternative %d."), J, J + 1L)),
    list("choice past the set",
         with_input(f, choice_idx = replace(f$choice_idx, 4L, f$M[4L] + 1L)),
         bad_choice),
    list("choice 0 without an outside option",
         with_input(f, choice_idx = replace(f$choice_idx, 4L, 0L)), bad_choice),
    list("choice NA", with_input(f, choice_idx = replace(f$choice_idx, 4L, NA_integer_)),
         bad_choice),
    list("choice NA, outside option",
         with_input(fo, choice_idx = replace(fo$choice_idx, 4L, NA_integer_)),
         bad_choice),
    list("choice -1, outside option",
         with_input(fo, choice_idx = replace(fo$choice_idx, 4L, -1L)),
         bad_choice),
    list("choice past the set, outside option",
         with_input(fo, choice_idx = replace(fo$choice_idx, 4L, fo$M[4L] + 1L)),
         bad_choice),
    list("design one row short and a code 0: the design's height first",
         with_input(f, X = f$X[-1L, , drop = FALSE],
                    alt_idx = replace(f$alt_idx, 5L, 0L)),
         sprintf("X has %d rows but sum(M) is %d.", n - 1L, n)),
    list("a code 0 and an NA choice: the codes first",
         with_input(f, alt_idx = replace(f$alt_idx, 5L, 0L),
                    choice_idx = replace(f$choice_idx, 4L, NA_integer_)),
         code_msg(0))
  )
  for (cs in cases) {
    for (k in mnl_kernels) {
      expected <- if (is.function(cs[[3L]])) cs[[3L]](k) else cs[[3L]]
      expect_identical(idx_err(k, cs[[2L]]), expected,
                       label = sprintf("%s: %s", cs[[1L]], k))
    }
  }
})

test_that("MNL kernels check their own arguments before the layout", {
  f <- idx_fixture(23L)
  m0 <- replace(f$M, 3L, 0L)
  # BLP: the target's length, then its signs
  expect_identical(
    idx_err("blp", with_input(f, M = m0, target = f$target[-1L])),
    paste("Error: target_shares must have the same length as the total",
          "number of alternatives."))
  expect_identical(
    idx_err("blp", with_input(f, M = m0,
                              target = replace(f$target[-1L], 1L, 0))),
    paste("Error: target_shares must have the same length as the total",
          "number of alternatives."))
  expect_identical(
    idx_err("blp", with_input(f, M = m0, target = replace(f$target, 2L, 0))),
    paste("Error: all target_shares must be strictly positive (log(share)",
          "is undefined otherwise)."))
  # Elasticities: the variable's column, then theta, then the layout
  elas <- function(theta, var, M) {
    tryCatch({
      mnl_elasticities_parallel(theta, f$X, f$alt_idx, f$choice_idx, M,
                                f$weights, var, TRUE, FALSE)
      "<no error>"
    }, error = conditionMessage)
  }
  expect_identical(elas(f$theta[1L], 0L, m0),
                   "elast_var_idx is out of bounds for design matrix X.")
  expect_identical(elas(f$theta[1L], 1L, m0), paste(
    "Theta vector too short: missing beta parameters (expected at least 2,",
    "got 1)."))
  expect_identical(elas(f$theta, 1L, m0),
                   "M must be positive for every individual (M[3] = 0).")
})

test_that("MNL kernels validate and add ASCs in parallel above 10^6 rows", {
  # The code scan and the ASC pass split their rows across threads above
  # 10^6 rows: a bad code in the last row is still reported, and the
  # utilities do not depend on the thread count.
  n <- 1000001L
  set.seed(51)
  X <- matrix(rnorm(n), n, 1L)
  M <- rep(1L, n)
  alt <- rep(1L, n)
  theta <- c(0.3, -0.2) # beta, and the ASC of the one inside alternative
  for (bad in list(NA_integer_, 0L)) {
    expect_error(mnl_predict(theta, X, replace(alt, n, bad), M, TRUE, TRUE),
                 sprintf("alt_idx must use 1-based alternative indices (found %s).",
                         if (is.na(bad)) "NA" else bad), fixed = TRUE)
  }
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  u1 <- mnl_predict(theta, X, alt, M, TRUE, TRUE)$utility
  set_num_threads(2L)
  u2 <- mnl_predict(theta, X, alt, M, TRUE, TRUE)$utility
  expect_identical(u2, u1)
  expect_equal(as.numeric(u1), as.numeric(X) * 0.3 - 0.2)
})

test_that("MNL kernels on empty input keep Armadillo's errors", {
  f <- idx_fixture(31L)
  e <- with_input(f, X = f$X[0L, , drop = FALSE], alt_idx = integer(0),
                  choice_idx = integer(0), M = integer(0), weights = numeric(0))
  expect_error(mnl_call("shares", e), "Mat::max(): object has no elements",
               fixed = TRUE, class = "std::logic_error")
  en <- with_input(e, use_asc = FALSE, theta = f$theta[seq_len(f$K)])
  expect_error(mnl_call("elasticities", en), "max(): object has no elements",
               fixed = TRUE, class = "std::logic_error")
  expect_error(mnl_call("diversion", en), "max(): object has no elements",
               fixed = TRUE, class = "std::logic_error")
  expect_identical(mnl_call("gradient", e),
                   list(objective = -0, gradient = matrix(0, length(f$theta), 1L)))
  # BLP stops at the weights before it reads delta, so an empty delta block
  # is never indexed (main read past it, or failed an Armadillo bounds check
  # with an outside option)
  for (ioo in c(FALSE, TRUE)) {
    expect_error(blp_contraction(numeric(0), if (ioo) 1 else numeric(0), e$X,
                                 f$theta[seq_len(f$K)], e$alt_idx, e$M,
                                 e$weights, ioo),
                 "Error: Sum of weights must be positive.", fixed = TRUE)
  }
})

test_that("MNL kernels never modify their index vectors", {
  for (ioo in c(FALSE, TRUE)) {
    f <- idx_fixture(41L + ioo, ioo = ioo)
    before <- unserialize(serialize(f, NULL))
    for (k in mnl_kernels) mnl_call(k, f)
    expect_identical(f, before)
  }
  # A compact sequence (ALTREP) is read like any integer vector: one
  # alternative per situation, no ASCs.
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  n <- 25L
  set.seed(42)
  X <- matrix(rnorm(2L * n), n, 2L)
  alt <- seq_len(n)
  plain <- seq_len(n) + 0L
  args <- function(a) list(c(0.2, -0.1), X, a, rep(1L, n), rep(1L, n),
                           rep(1, n), FALSE, FALSE)
  expect_identical(do.call(mnl_loglik_gradient_parallel, args(alt)),
                   do.call(mnl_loglik_gradient_parallel, args(plain)))
  expect_identical(alt, seq_len(n))
})

# Nested logit ----------------------------------------------------------------

# Five alternatives in three nests, the third a singleton (lambda fixed at
# 1); situation 1 offers every alternative.
nl_idx_fixture <- function(seed, N = 12L, ioo = FALSE, use_asc = TRUE) {
  set.seed(seed)
  J <- 5L
  K <- 2L
  nest <- c(1L, 1L, 2L, 2L, 3L)
  lambda <- c(0.7, 0.8) # the two nests of more than one alternative
  sets <- c(list(seq_len(J)),
            lapply(seq_len(N - 1L),
                   function(i) sort(sample.int(J, sample(2:J, 1L)))))
  M <- lengths(sets)
  n <- sum(M)
  lo <- if (ioo) 0L else 1L
  n_delta <- if (!use_asc) 0L else if (ioo) J else J - 1L
  f <- list(theta = c(rnorm(K, sd = 0.3), lambda, rnorm(n_delta, sd = 0.3)),
            X = matrix(rnorm(n * K), n, K), alt_idx = unlist(sets),
            choice_idx = vapply(M, function(m) sample(lo:m, 1L), integer(1)),
            nest_idx = nest, M = M, weights = runif(N, 0.5, 2),
            lambda_full = c(lambda, 1), J = J, K = K, N = N, ioo = ioo,
            use_asc = use_asc)
  # BLP target: the shares the model predicts with ASCs of 0.1
  f$target <- as.numeric(nl_predict_shares(
    c(f$theta[seq_len(K)], lambda, rep(0.1, if (ioo) J else J - 1L)), f$X,
    f$alt_idx, f$M, f$weights, f$nest_idx, TRUE, ioo))
  f
}

nl_kernels <- c("gradient", "bhhh", "scores", "hessian", "num_hessian",
                "predict", "shares", "elasticities", "diversion", "blp")

nl_call <- function(k, f) {
  switch(k,
    gradient = nl_loglik_gradient_parallel(f$theta, f$X, f$alt_idx,
      f$choice_idx, f$nest_idx, f$M, f$weights, f$use_asc, f$ioo),
    bhhh = nl_bhhh_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
      f$nest_idx, f$M, f$weights, f$use_asc, f$ioo),
    scores = nl_scores_parallel(f$theta, f$X, f$alt_idx, f$choice_idx,
      f$nest_idx, f$M, f$use_asc, f$ioo),
    hessian = nl_loglik_hessian_parallel(f$theta, f$X, f$alt_idx,
      f$choice_idx, f$nest_idx, f$M, f$weights, f$use_asc, f$ioo),
    num_hessian = nl_loglik_numeric_hessian(f$theta, f$X, f$alt_idx,
      f$choice_idx, f$nest_idx, f$M, f$weights, f$use_asc, f$ioo),
    predict = nl_predict(f$theta, f$X, f$alt_idx, f$M, f$nest_idx, f$use_asc,
      f$ioo),
    shares = nl_predict_shares(f$theta, f$X, f$alt_idx, f$M, f$weights,
      f$nest_idx, f$use_asc, f$ioo),
    elasticities = nl_elasticities_parallel(f$theta, f$X, f$alt_idx,
      f$choice_idx, f$nest_idx, f$M, f$weights, 1L, f$use_asc, f$ioo),
    diversion = nl_diversion_ratios_parallel(f$theta, f$X, f$alt_idx,
      f$nest_idx, f$M, f$weights, f$use_asc, f$ioo),
    blp = nl_blp_contraction(rep(0, f$J), f$target, f$X, f$theta[seq_len(f$K)],
      f$lambda_full, f$alt_idx, f$nest_idx, f$M, f$weights, f$ioo),
    stop("unknown kernel: ", k))
}

nl_err <- function(k, f) {
  warns <- character()
  msg <- tryCatch(
    withCallingHandlers({
      nl_call(k, f)
      "<no error>"
    }, warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = conditionMessage)
  c(msg, warns)
}

test_that("the NL kernels give identical results for double-typed indices", {
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L) # one thread: two calls sum in the same order
  for (ioo in c(FALSE, TRUE)) {
    for (use_asc in c(TRUE, FALSE)) {
      f <- nl_idx_fixture(61L + ioo, ioo = ioo, use_asc = use_asc)
      fd <- with_input(f, alt_idx = as.double(f$alt_idx),
                       choice_idx = as.double(f$choice_idx),
                       nest_idx = as.double(f$nest_idx), M = as.double(f$M))
      for (k in nl_kernels) {
        expect_identical(nl_call(k, fd), nl_call(k, f), label = sprintf(
          "%s (outside option %s, ASCs %s) with double-typed indices",
          k, ioo, use_asc))
      }
    }
  }
})

test_that("NL input errors give the same message in every kernel", {
  f <- nl_idx_fixture(71L)                       # no outside option, ASCs
  fo <- nl_idx_fixture(72L, ioo = TRUE)          # outside option
  fn <- nl_idx_fixture(73L, use_asc = FALSE)     # no ASCs
  n <- sum(f$M)
  N <- f$N
  J <- f$J
  choice_kernels <- c("gradient", "bhhh", "scores", "hessian", "num_hessian")
  weight_kernels <- c("gradient", "bhhh", "hessian", "num_hessian", "shares",
                      "elasticities", "diversion", "blp")
  suffix <- c(bhhh = " (nl_bhhh_parallel)", scores = " (nl_scores_parallel)")
  bad_choice <- function(k) {
    if (!k %in% choice_kernels) return("<no error>")
    paste0("Invalid chosen alternative index for individual 3",
           if (k %in% names(suffix)) suffix[[k]] else "")
  }
  code_msg <- function(found) {
    sprintf("alt_idx must use 1-based alternative indices (found %s).", found)
  }
  delta_msg <- sprintf(paste("Theta's delta (ASC) block implies %d",
                             "alternatives but alt_idx references alternative",
                             "%d."), J, J + 1L)
  bad_nest <- function(found) function(k) {
    if (k == "blp") {
      sprintf("nest_idx must use 1-based nest indices (found %s).", found)
    } else "Invalid nest index found in nest_idx."
  }
  cases <- list(
    list("M = 0", with_input(f, M = replace(f$M, 3L, 0L)),
         "M must be positive for every individual (M[3] = 0)."),
    list("X one row short", with_input(f, X = f$X[-1L, , drop = FALSE]),
         sprintf("X has %d rows but sum(M) is %d.", n - 1L, n)),
    list("alt_idx one short", with_input(f, alt_idx = f$alt_idx[-1L]),
         sprintf(paste("alt_idx length (%d) does not match the number of",
                       "rows of X (%d)."), n - 1L, n)),
    list("weights one short", with_input(f, weights = f$weights[-1L]),
         function(k) if (k %in% weight_kernels) {
           sprintf("weights length (%d) does not match N (%d)", N - 1L, N)
         } else "<no error>"),
    list("choice_idx one short", with_input(f, choice_idx = f$choice_idx[-1L]),
         function(k) if (k %in% choice_kernels) {
           sprintf("choice_idx length (%d) does not match N (%d)", N - 1L, N)
         } else "<no error>"),
    list("alternative code 0",
         with_input(f, alt_idx = replace(f$alt_idx, 5L, 0L)), code_msg(0)),
    list("alternative code NA",
         with_input(f, alt_idx = replace(f$alt_idx, 5L, NA_integer_)),
         code_msg("NA")),
    list("alternative code -1",
         with_input(f, alt_idx = replace(f$alt_idx, 5L, -1L)), code_msg(-1)),
    list("alternative code 3e9 (a double beyond the integer range)",
         with_input(f, alt_idx = replace(as.double(f$alt_idx), 5L, 3e9)),
         c(code_msg("NA"), "NAs introduced by coercion to integer range")),
    list("alternative code past the delta block",
         with_input(f, alt_idx = replace(f$alt_idx, 5L, J + 1L)), delta_msg),
    list("alternative code past nest_idx, no ASCs",
         with_input(fn, alt_idx = replace(fn$alt_idx, 5L, J + 1L)),
         function(k) if (k == "blp") delta_msg else sprintf(paste(
           "nest_idx has %d entries but alt_idx references alternative %d",
           "(one nest index per global alternative is required)."), J, J + 1L)),
    list("nest code 0", with_input(f, nest_idx = replace(f$nest_idx, 2L, 0L)),
         bad_nest(0)),
    list("nest code NA",
         with_input(f, nest_idx = replace(f$nest_idx, 2L, NA_integer_)),
         bad_nest("NA")),
    list("nest code -1", with_input(f, nest_idx = replace(f$nest_idx, 2L, -1L)),
         bad_nest(-1)),
    list("every nest code NA",
         with_input(f, nest_idx = rep(NA_integer_, length(f$nest_idx))),
         bad_nest("NA")),
    list("nest code 3e9 (a double beyond the integer range)",
         with_input(f, nest_idx = replace(as.double(f$nest_idx), 2L, 3e9)),
         function(k) c(bad_nest("NA")(k),
                       "NAs introduced by coercion to integer range")),
    list("choice past the set",
         with_input(f, choice_idx = replace(f$choice_idx, 4L, f$M[4L] + 1L)),
         bad_choice),
    list("choice 0 without an outside option",
         with_input(f, choice_idx = replace(f$choice_idx, 4L, 0L)), bad_choice),
    list("choice NA",
         with_input(f, choice_idx = replace(f$choice_idx, 4L, NA_integer_)),
         bad_choice),
    list("choice NA, outside option",
         with_input(fo, choice_idx = replace(fo$choice_idx, 4L, NA_integer_)),
         bad_choice),
    list("choice -1, outside option",
         with_input(fo, choice_idx = replace(fo$choice_idx, 4L, -1L)),
         bad_choice),
    list("design one row short and a nest code 0: the nest codes first",
         with_input(f, X = f$X[-1L, , drop = FALSE],
                    nest_idx = replace(f$nest_idx, 2L, 0L)),
         function(k) if (k == "blp") {
           sprintf("X has %d rows but sum(M) is %d.", n - 1L, n)
         } else "Invalid nest index found in nest_idx."),
    list("a code 0 and an NA choice: the codes first",
         with_input(f, alt_idx = replace(f$alt_idx, 5L, 0L),
                    choice_idx = replace(f$choice_idx, 4L, NA_integer_)),
         code_msg(0))
  )
  for (cs in cases) {
    for (k in nl_kernels) {
      expected <- if (is.function(cs[[3L]])) cs[[3L]](k) else cs[[3L]]
      expect_identical(nl_err(k, cs[[2L]]), expected,
                       label = sprintf("%s: %s", cs[[1L]], k))
    }
  }
  # The messages are Rcpp's exceptions, as on main
  expect_error(nl_call("gradient",
                       with_input(f, nest_idx = replace(f$nest_idx, 2L, 0L))),
               "Invalid nest index found in nest_idx.", fixed = TRUE,
               class = "Rcpp::exception")
  # An empty nest_idx keeps arma::max()'s error, in every kernel
  e <- with_input(f, nest_idx = integer(0))
  for (k in nl_kernels) {
    expect_error(nl_call(k, e), "max(): object has no elements", fixed = TRUE,
                 class = "std::logic_error", label = k)
  }
})

test_that("NL kernels on empty input keep Armadillo's errors", {
  f <- nl_idx_fixture(81L)
  e <- with_input(f, X = f$X[0L, , drop = FALSE], alt_idx = integer(0),
                  choice_idx = integer(0), M = integer(0), weights = numeric(0))
  expect_error(nl_call("shares", e), "Mat::max(): object has no elements",
               fixed = TRUE, class = "std::logic_error")
  en <- with_input(e, use_asc = FALSE, theta = f$theta[seq_len(f$K + 2L)])
  expect_error(nl_call("elasticities", en), "max(): object has no elements",
               fixed = TRUE, class = "std::logic_error")
  expect_error(nl_call("diversion", en), "max(): object has no elements",
               fixed = TRUE, class = "std::logic_error")
  # BLP stops at the weights before it reads delta, so an empty delta block
  # is never indexed (main read past it, or failed an Armadillo bounds check
  # with an outside option)
  for (ioo in c(FALSE, TRUE)) {
    expect_error(nl_blp_contraction(numeric(0), if (ioo) 1 else numeric(0),
                                    e$X, f$theta[seq_len(f$K)], f$lambda_full,
                                    e$alt_idx, f$nest_idx, e$M, e$weights, ioo),
                 "Error: Sum of weights must be positive.", fixed = TRUE)
  }
})

test_that("NL kernels never modify their index vectors", {
  for (ioo in c(FALSE, TRUE)) {
    f <- nl_idx_fixture(91L + ioo, ioo = ioo)
    before <- unserialize(serialize(f, NULL))
    for (k in nl_kernels) nl_call(k, f)
    expect_identical(f, before)
  }
})
