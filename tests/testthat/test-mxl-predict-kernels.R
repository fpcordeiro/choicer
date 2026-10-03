# The mixed logit prediction kernels read their alternative codes in place and
# work one choice situation at a time in per-thread buffers. These tests pin
# their inputs and errors, check them against a brute-force R oracle written
# from the estimand, and check that a situation's outputs do not depend on the
# thread count.

# --- Fixtures ----------------------------------------------------------------

# A stacked design of N situations over J alternatives. Every situation sees M[t]
# distinct alternatives (codes from `codes`, all of 1..J by default); W is
# row-aligned or alternative-level (J x K_w).
mxp_data <- function(seed, M, J, K_x = 2L, K_w = 2L, S = 5L, ioo = FALSE,
                     use_asc = TRUE, w_type = "row", rc_dist = NULL,
                     rc_corr = FALSE, rc_mean = TRUE, codes = seq_len(J)) {
  set.seed(seed)
  N <- length(M)
  alt_idx <- as.integer(unlist(lapply(M, function(m) sort(sample(codes, m)))))
  n <- sum(M)
  X <- matrix(rnorm(n * K_x), n, K_x)
  W <- if (w_type == "row") matrix(rnorm(n * K_w), n, K_w) else
    matrix(rnorm(J * K_w), J, K_w)
  rc_dist <- if (is.null(rc_dist)) rep(0L, K_w) else as.integer(rc_dist)
  L_size <- if (rc_corr) K_w * (K_w + 1L) / 2L else K_w
  n_delta <- if (!use_asc) 0L else if (ioo) J else J - 1L
  theta <- rnorm(K_x + (if (rc_mean) K_w else 0L) + L_size + n_delta, sd = 0.4)
  eta <- get_halton_normals(S, N, K_w)
  list(theta = theta, X = X, W = W, alt_idx = alt_idx, M = as.integer(M),
       eta = eta, rc_dist = rc_dist, rc_corr = rc_corr, rc_mean = rc_mean,
       use_asc = use_asc, ioo = ioo, weights = runif(N, 0.5, 2))
}

mxp_predict <- function(d, ...) {
  choicer:::mxl_predict(d$theta, d$X, d$W, d$alt_idx, d$M, d$eta, d$rc_dist,
                        d$rc_corr, d$rc_mean, d$use_asc, d$ioo, ...)
}
mxp_logsum <- function(d, ...) {
  choicer:::mxl_logsum(d$theta, d$X, d$W, d$alt_idx, d$M, d$eta, d$rc_dist,
                       d$rc_corr, d$rc_mean, d$use_asc, d$ioo, ...)
}

# The generate-mode draws of situation t (1-based): Halton block t, indices
# (t - 1) S + 1, ..., t S.
mxp_gen_cube <- function(N, S, K_w, seed, scramble) {
  eta <- array(0, c(K_w, S, N))
  for (t in seq_len(N)) {
    eta[, , t] <- choicer:::halton_fill_block((t - 1) * S + 1, S, K_w, seed,
                                              scramble)
  }
  eta
}

# Brute-force predictions from the estimand: for situation t and draw s, the
# random coefficients mu~ + gamma_ts (gamma_ts = L eta_ts, exponentiated in
# log-normal rows; mu~ = exp(mu) there, mu otherwise, 0 without rc_mean), the
# utilities V = X beta + W (mu~ + gamma) + delta, and the logit probabilities
# with the outside option's utility 0 when present.
mxp_oracle <- function(d, eta = d$eta) {
  X <- d$X; W <- d$W; K_x <- ncol(X); K_w <- ncol(W); th <- d$theta
  beta <- th[seq_len(K_x)]; pos <- K_x
  mu <- rep(0, K_w)
  if (d$rc_mean) {
    mu <- th[pos + seq_len(K_w)]
    pos <- pos + K_w
  }
  mu_f <- ifelse(d$rc_dist == 1L & d$rc_mean, exp(mu), mu)
  L <- matrix(0, K_w, K_w)
  if (d$rc_corr) {
    for (i in seq_len(K_w)) for (j in seq_len(i)) {
      pos <- pos + 1L
      L[i, j] <- if (i == j) exp(th[pos]) else th[pos]
    }
  } else if (K_w > 0L) {
    diag(L) <- exp(th[pos + seq_len(K_w)])
    pos <- pos + K_w
  }
  delta <- if (!d$use_asc) NULL else if (d$ioo) th[-seq_len(pos)] else
    c(0, th[-seq_len(pos)])
  alt_level <- nrow(W) != nrow(X)
  ends <- cumsum(d$M); starts <- ends - d$M + 1L
  S <- dim(eta)[2]
  prob <- util <- numeric(sum(d$M))
  prob_out <- ls <- numeric(length(d$M))
  for (t in seq_along(d$M)) {
    r <- starts[t]:ends[t]
    W_t <- if (alt_level) W[d$alt_idx[r], , drop = FALSE] else W[r, , drop = FALSE]
    base <- drop(X[r, , drop = FALSE] %*% beta)
    if (d$use_asc) base <- base + delta[d$alt_idx[r]]
    P <- 0; U <- 0; Po <- 0; Ls <- 0
    for (s in seq_len(S)) {
      g <- if (K_w > 0L) drop(L %*% eta[, s, t]) else numeric(0)
      g[d$rc_dist == 1L] <- exp(g[d$rc_dist == 1L])
      V <- base + drop(W_t %*% (mu_f + g))
      V_all <- if (d$ioo) c(0, V) else V
      e <- exp(V_all - max(V_all))
      p <- e / sum(e)
      if (d$ioo) {
        Po <- Po + p[1]
        p <- p[-1]
      }
      P <- P + p; U <- U + V
      Ls <- Ls + max(V_all) + log(sum(e))
    }
    prob[r] <- P / S; util[r] <- U / S
    prob_out[t] <- Po / S; ls[t] <- Ls / S
  }
  list(prob = prob, util = util, prob_out = prob_out, logsum = ls)
}

mxp_expect_oracle <- function(d, gen = NULL) {
  if (is.null(gen)) {
    p <- mxp_predict(d)
    ls <- mxp_logsum(d)
    o <- mxp_oracle(d)
  } else {
    p <- mxp_predict(d, gen_seed = gen$seed, gen_scramble = gen$scramble,
                     gen_S = gen$S)
    ls <- mxp_logsum(d, gen_seed = gen$seed, gen_scramble = gen$scramble,
                     gen_S = gen$S)
    o <- mxp_oracle(d, mxp_gen_cube(length(d$M), gen$S, ncol(d$W), gen$seed,
                                    gen$scramble))
  }
  expect_equal(as.numeric(p$choice_prob), o$prob, tolerance = 1e-12)
  expect_equal(as.numeric(p$utility), o$util, tolerance = 1e-12)
  if (d$ioo) {
    expect_equal(as.numeric(p$choice_prob_outside), o$prob_out, tolerance = 1e-12)
  }
  expect_equal(as.numeric(ls), o$logsum, tolerance = 1e-12)
}

# --- Values ------------------------------------------------------------------

test_that("mxl_predict and mxl_logsum match a brute-force oracle", {
  set.seed(1)
  M <- sample(1:5, 12, replace = TRUE)
  for (ioo in c(FALSE, TRUE)) {
    for (w_type in c("row", "alt")) {
      d <- mxp_data(10 + ioo, M, J = 5L, K_w = 2L, ioo = ioo, w_type = w_type,
                    rc_dist = c(0L, 1L), rc_corr = TRUE)
      mxp_expect_oracle(d)
      mxp_expect_oracle(d, gen = list(seed = 17L, scramble = 1L, S = 6L))
    }
  }
  # Without ASCs or random-coefficient means, and with codes starting at 3
  d <- mxp_data(20, rep(3L, 9), J = 6L, K_w = 3L, use_asc = FALSE,
                rc_mean = FALSE, codes = 3:6)
  mxp_expect_oracle(d)
})

test_that("tiny designs and a single draw match the oracle", {
  # Designs of 1 to 5 rows: Armadillo forms X beta by dgemv('T') for one row
  # and by its own code for square designs of up to 4 x 4.
  for (n in 1:5) {
    for (K_x in unique(c(n, 2L))) {
      d <- mxp_data(30 + n, n, J = 5L, K_x = K_x, K_w = 1L)
      mxp_expect_oracle(d)
    }
  }
  d <- mxp_data(40, c(1L, 4L, 2L), J = 4L, S = 1L, ioo = TRUE)
  mxp_expect_oracle(d)
  # m = K_w = S = 3 in every situation: every product is a tiny square one
  d <- mxp_data(42, rep(3L, 6), J = 3L, K_w = 3L, S = 3L, rc_corr = TRUE)
  mxp_expect_oracle(d)
})

test_that("a situation's predictions do not depend on the thread count", {
  set.seed(2)
  d <- mxp_data(50, sample(2:20, 300, replace = TRUE), J = 25L, K_w = 3L,
                S = 7L, ioo = TRUE, rc_dist = c(0L, 1L, 0L), rc_corr = TRUE)
  on.exit(set_num_threads(2L), add = TRUE)
  gen <- function(k) k(d, gen_seed = 3L, gen_scramble = 1L, gen_S = 7L)
  set_num_threads(1L)
  p1 <- mxp_predict(d); l1 <- mxp_logsum(d)
  gp1 <- gen(mxp_predict); gl1 <- gen(mxp_logsum)
  set_num_threads(2L)
  expect_identical(mxp_predict(d), p1)
  expect_identical(mxp_logsum(d), l1)
  expect_identical(gen(mxp_predict), gp1)
  expect_identical(gen(mxp_logsum), gl1)
})

test_that("double-typed alternative codes give the integer result", {
  d <- mxp_data(60, c(2L, 3L, 4L, 2L, 3L), J = 4L, ioo = TRUE)
  dd <- d
  dd$alt_idx <- as.double(d$alt_idx)
  expect_identical(mxp_predict(dd), mxp_predict(d))
  expect_identical(mxp_logsum(dd), mxp_logsum(d))
})

test_that("predictions return n x 1 matrices and leave their inputs alone", {
  d <- mxp_data(70, c(2L, 3L, 4L), J = 4L, ioo = TRUE)
  before <- serialize(d, NULL)
  p <- mxp_predict(d)
  expect_named(p, c("choice_prob", "utility", "choice_prob_outside"))
  expect_identical(dim(p$choice_prob), c(9L, 1L))
  expect_identical(dim(p$utility), c(9L, 1L))
  expect_identical(dim(p$choice_prob_outside), c(3L, 1L))
  expect_type(p$choice_prob, "double")
  expect_identical(dim(mxp_logsum(d)), c(3L, 1L))
  expect_identical(serialize(d, NULL), before)
  expect_named(mxp_predict(mxp_data(71, c(2L, 3L), J = 3L)),
               c("choice_prob", "utility"))
})

test_that("predictions with no choice situations are empty", {
  d <- mxp_data(80, c(2L, 3L), J = 3L, ioo = TRUE)
  d$X <- d$X[0, , drop = FALSE]; d$W <- d$W[0, , drop = FALSE]
  d$alt_idx <- integer(0); d$M <- integer(0); d$eta <- d$eta[, , 0, drop = FALSE]
  p <- mxp_predict(d)
  expect_identical(dim(p$choice_prob), c(0L, 1L))
  expect_identical(dim(p$choice_prob_outside), c(0L, 1L))
  expect_identical(dim(mxp_logsum(d)), c(0L, 1L))
})

# --- Errors ------------------------------------------------------------------

test_that("prediction kernels report malformed alternative codes", {
  d <- mxp_data(90, c(2L, 3L, 4L, 2L), J = 4L)
  for (k in list(mxp_predict, mxp_logsum)) {
    for (gen in c(FALSE, TRUE)) {
      call_k <- function(alt) {
        d$alt_idx <- alt
        if (gen) k(d, gen_seed = 1L, gen_scramble = 1L, gen_S = 5L) else k(d)
      }
      expect_error(call_k(replace(d$alt_idx, 3L, NA_integer_)),
                   "alt_idx must use 1-based alternative indices (found NA).",
                   fixed = TRUE)
      expect_error(call_k(replace(d$alt_idx, 3L, -1L)),
                   "alt_idx must use 1-based alternative indices (found -1).",
                   fixed = TRUE)
      expect_error(call_k(replace(d$alt_idx, 3L, 0L)),
                   "alt_idx must use 1-based alternative indices (found 0).",
                   fixed = TRUE)
      expect_warning(
        expect_error(call_k(replace(as.double(d$alt_idx), 3L, 3e9)),
                     "(found NA)", fixed = TRUE),
        "NAs introduced by coercion to integer range")
      expect_error(call_k(replace(d$alt_idx, 3L, 5L)),
                   "Theta's delta (ASC) block implies 4 alternatives but alt_idx references alternative 5.",
                   fixed = TRUE)
    }
  }
})

test_that("prediction kernels check the design, the draws and W", {
  d <- mxp_data(91, c(2L, 3L, 4L, 2L), J = 4L)
  for (k in list(mxp_predict, mxp_logsum)) {
    e <- d; e$M[2] <- 0L
    expect_error(k(e), "M must be positive for every individual (M[2] = 0).",
                 fixed = TRUE)
    e <- d; e$X <- e$X[-1, , drop = FALSE]
    expect_error(k(e), "X has 10 rows but sum(M) is 11.", fixed = TRUE)
    e <- d; e$alt_idx <- e$alt_idx[-1]
    expect_error(k(e), "alt_idx length (10) does not match the number of rows of X (11).",
                 fixed = TRUE)
    e <- d; e$eta <- get_halton_normals(5L, 5L, 2L)
    expect_error(k(e), "eta_draws 3rd dimension (5) does not match N (4)",
                 fixed = TRUE)
    e <- d; e$eta <- get_halton_normals(5L, 4L, 3L)
    expect_error(k(e), "eta_draws 1st dimension (3) does not match K_w (2)",
                 fixed = TRUE)
    e <- d; e$theta <- e$theta[1:3]
    expect_error(k(e), "Theta vector too short")
    # An alternative-level W must cover every code, in both draw modes
    e <- d; e$W <- matrix(1, 3L, 2L)
    msg <- "W must be row-aligned with X (11 rows) or contain one row per global alternative (at least 4 rows); got 3 rows."
    expect_error(k(e), msg, fixed = TRUE)
    expect_error(k(e, gen_seed = 1L, gen_scramble = 1L, gen_S = 5L), msg,
                 fixed = TRUE)
    expect_error(k(d, gen_seed = 1L, gen_scramble = 1L, gen_S = 0L),
                 "gen_S must be positive when gen_seed >= 0", fixed = TRUE)
  }
  # Generated draws have one Halton base per random coefficient, up to 128
  wide <- mxp_data(92, c(2L, 3L), J = 3L, K_w = 129L, S = 1L)
  for (k in list(mxp_predict, mxp_logsum)) {
    expect_error(k(wide, gen_seed = 1L, gen_scramble = 1L, gen_S = 2L),
                 "K_w exceeds the primes table size (128)", fixed = TRUE)
  }
})
