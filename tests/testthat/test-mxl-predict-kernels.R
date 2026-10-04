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
       use_asc = use_asc, ioo = ioo, weights = runif(N, 0.5, 2), J = J)
}

mxp_predict <- function(d, ...) {
  choicer:::mxl_predict(d$theta, d$X, d$W, d$alt_idx, d$M, d$eta, d$rc_dist,
                        d$rc_corr, d$rc_mean, d$use_asc, d$ioo, ...)
}
mxp_logsum <- function(d, ...) {
  choicer:::mxl_logsum(d$theta, d$X, d$W, d$alt_idx, d$M, d$eta, d$rc_dist,
                       d$rc_corr, d$rc_mean, d$use_asc, d$ioo, ...)
}
mxp_shares <- function(d, ...) {
  choicer:::mxl_predict_shares(d$theta, d$X, d$W, d$alt_idx, d$M, d$weights,
                               d$eta, d$rc_dist, d$rc_corr, d$rc_mean,
                               d$use_asc, d$ioo, ...)
}
mxp_elas <- function(d, var, random, ...) {
  choicer:::mxl_elasticities_parallel(d$theta, d$X, d$W, d$alt_idx,
                                      rep(1L, length(d$M)), d$M, d$weights,
                                      d$eta, d$rc_dist, var, random, d$rc_corr,
                                      d$rc_mean, d$use_asc, d$ioo, ...)
}
mxp_dr <- function(d, var, random, ...) {
  choicer:::mxl_diversion_ratios_parallel(d$theta, d$X, d$W, d$alt_idx, d$M,
                                          d$weights, d$eta, d$rc_dist, var,
                                          random, d$rc_corr, d$rc_mean,
                                          d$use_asc, d$ioo, ...)
}
# BLP inversion from `delta0` towards `target`, at the design's beta, mu and L
mxp_blp <- function(d, target, delta0, beta = NULL, ...) {
  K_x <- ncol(d$X); K_w <- ncol(d$W)
  pos <- K_x + (if (d$rc_mean) K_w else 0L)
  L_size <- if (d$rc_corr) K_w * (K_w + 1L) / 2L else K_w
  if (is.null(beta)) beta <- d$theta[seq_len(K_x)]
  mu <- if (d$rc_mean) d$theta[K_x + seq_len(K_w)] else rep(0, K_w)
  mxl_blp_contraction(delta0, target, d$X, d$W, beta, mu,
                      d$theta[pos + seq_len(L_size)], d$alt_idx, d$M,
                      d$weights, d$eta, d$rc_dist, d$rc_corr, d$rc_mean,
                      d$ioo, ...)
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

# The oracle's weighted shares: one per alternative code up to the largest
# (all J with ASCs), the outside option first when present
mxp_shares_oracle <- function(d, eta = d$eta) {
  o <- mxp_oracle(d, eta)
  J_in <- if (d$use_asc) d$J else max(d$alt_idx)
  sh <- numeric(J_in + d$ioo)
  ends <- cumsum(d$M); starts <- ends - d$M + 1L
  for (t in seq_along(d$M)) {
    if (d$ioo) sh[1] <- sh[1] + d$weights[t] * o$prob_out[t]
    for (r in starts[t]:ends[t]) {
      j <- d$alt_idx[r] + d$ioo
      sh[j] <- sh[j] + d$weights[t] * o$prob[r]
    }
  }
  sh / sum(d$weights)
}

# Brute-force aggregate elasticities and diversion ratios with respect to
# variable `var` of X (fixed) or W (random): per situation and draw, the
# realized coefficient b (beta_var, or mu~_var + gamma_var), the
# probabilities P (outside option first when present) and the variable's
# values x (0 for the outside option); then
#   E_t(j, m) = mean_s b x_m P_j (1{j = m} - P_m) / mean_s P_j
#   DR(k, j)  = sum_t w_t mean_s b P_j P_k / sum_t w_t mean_s b P_j (1 - P_j)
# on the alternatives' global slots, E averaged with the weights.
mxp_elas_dr_oracle <- function(d, var, random, eta = d$eta) {
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
  } else {
    diag(L) <- exp(th[pos + seq_len(K_w)])
    pos <- pos + K_w
  }
  delta <- if (!d$use_asc) NULL else if (d$ioo) th[-seq_len(pos)] else
    c(0, th[-seq_len(pos)])
  alt_level <- nrow(W) != nrow(X)
  J_in <- if (d$use_asc) d$J else max(d$alt_idx)
  J_tot <- J_in + d$ioo
  E <- num <- matrix(0, J_tot, J_tot)
  den <- numeric(J_tot)
  ends <- cumsum(d$M); starts <- ends - d$M + 1L
  S <- dim(eta)[2]
  for (t in seq_along(d$M)) {
    r <- starts[t]:ends[t]
    W_t <- if (alt_level) W[d$alt_idx[r], , drop = FALSE] else W[r, , drop = FALSE]
    base <- drop(X[r, , drop = FALSE] %*% beta)
    if (d$use_asc) base <- base + delta[d$alt_idx[r]]
    g_slot <- c(if (d$ioo) 1L, d$alt_idx[r] + d$ioo)
    x <- c(if (d$ioo) 0, if (random) W_t[, var] else X[r, var])
    n <- length(g_slot)
    acc_E <- matrix(0, n, n); acc_P <- numeric(n)
    acc_num <- matrix(0, n, n); acc_den <- numeric(n)
    for (s in seq_len(S)) {
      g <- drop(L %*% eta[, s, t])
      g[d$rc_dist == 1L] <- exp(g[d$rc_dist == 1L])
      V <- base + drop(W_t %*% (mu_f + g))
      V_all <- if (d$ioo) c(0, V) else V
      P <- exp(V_all - max(V_all)); P <- P / sum(P)
      b <- if (random) mu_f[var] + g[var] else beta[var]
      acc_P <- acc_P + P
      acc_E <- acc_E + b * outer(P, x) * (diag(n) - matrix(P, n, n, byrow = TRUE))
      acc_den <- acc_den + b * P * (1 - P)
      acc_num <- acc_num + b * outer(P, P)   # [j, k] = b P_j P_k
    }
    P_bar <- acc_P / S
    E_t <- (acc_E / S) / P_bar
    E_t[P_bar <= 1e-12, ] <- 0
    E[g_slot, g_slot] <- E[g_slot, g_slot] + d$weights[t] * E_t
    nk <- t(acc_num / S); diag(nk) <- 0   # [k, j]
    num[g_slot, g_slot] <- num[g_slot, g_slot] + d$weights[t] * nk
    den[g_slot] <- den[g_slot] + d$weights[t] * acc_den / S
  }
  DR <- sweep(num, 2, den, "/")
  DR[, abs(den) <= 1e-15] <- 0
  diag(DR) <- 0
  list(elas = E / sum(d$weights), dr = DR)
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
  if (is.null(gen)) {
    expect_equal(as.numeric(mxp_shares(d)), mxp_shares_oracle(d),
                 tolerance = 1e-12)
  } else {
    expect_equal(
      as.numeric(mxp_shares(d, gen_seed = gen$seed, gen_scramble = gen$scramble,
                            gen_S = gen$S)),
      mxp_shares_oracle(d, mxp_gen_cube(length(d$M), gen$S, ncol(d$W),
                                        gen$seed, gen$scramble)),
      tolerance = 1e-12)
  }
}

# Elasticities and diversion ratios with respect to the first fixed and the
# first random coefficient, against the oracle
mxp_expect_elas_dr_oracle <- function(d) {
  for (random in c(FALSE, TRUE)) {
    o <- mxp_elas_dr_oracle(d, 1L, random)
    expect_equal(unname(mxp_elas(d, 1L, random)), o$elas, tolerance = 1e-12)
    expect_equal(unname(mxp_dr(d, 1L, random)), o$dr, tolerance = 1e-12)
  }
}

# --- Values ------------------------------------------------------------------

test_that("predictions, log-sums and shares match a brute-force oracle", {
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
  d$ioo <- TRUE
  mxp_expect_oracle(d)
})

test_that("tiny designs and a single draw match the oracle", {
  # Designs of 1 to 5 rows: Armadillo forms X beta by dgemv('T') for one row
  # and by its own code for square designs of up to 4 x 4.
  for (n in 1:5) {
    for (K_x in unique(c(n, 2L))) {
      d <- mxp_data(30 + n, n, J = 5L, K_x = K_x, K_w = 1L)
      mxp_expect_oracle(d)
      mxp_expect_elas_dr_oracle(d)
    }
  }
  d <- mxp_data(40, c(1L, 4L, 2L), J = 4L, S = 1L, ioo = TRUE)
  mxp_expect_oracle(d)
  mxp_expect_elas_dr_oracle(d)
  # m = K_w = S = 3 in every situation: every product is a tiny square one
  d <- mxp_data(42, rep(3L, 6), J = 3L, K_w = 3L, S = 3L, rc_corr = TRUE)
  mxp_expect_oracle(d)
  mxp_expect_elas_dr_oracle(d)
})

test_that("the BLP contraction inverts the simulated shares", {
  set.seed(4)
  M <- sample(2:5, 40, replace = TRUE)
  for (ioo in c(FALSE, TRUE)) {
    for (gen in c(FALSE, TRUE)) {
      w_type <- if (gen) "alt" else "row"
      d <- mxp_data(15 + ioo, M, J = 5L, K_w = 2L, ioo = ioo, w_type = w_type,
                    rc_dist = c(1L, 0L), rc_corr = TRUE)
      g <- if (gen) list(gen_seed = 9L, gen_scramble = 1L, gen_S = 5L) else list()
      target <- as.numeric(do.call(mxp_shares, c(list(d), g)))
      n_free <- if (ioo) 5L else 4L
      delta <- if (ioo) d$theta[length(d$theta) - n_free + seq_len(n_free)] else
        c(0, d$theta[length(d$theta) - n_free + seq_len(n_free)])
      est <- do.call(mxp_blp, c(list(d, target, rep(0, 5L), tol = 1e-13), g))
      expect_equal(as.numeric(est), delta, tolerance = 1e-8)
    }
  }
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
  # Shares add situations up in each thread, so only the rounding of the sum
  # depends on the thread count
  set_num_threads(1L)
  s1 <- mxp_shares(d); e1 <- mxp_elas(d, 2L, TRUE); r1 <- mxp_dr(d, 1L, FALSE)
  set_num_threads(2L)
  expect_equal(mxp_shares(d), s1, tolerance = 1e-13)
  expect_equal(mxp_elas(d, 2L, TRUE), e1, tolerance = 1e-13)
  expect_equal(mxp_dr(d, 1L, FALSE), r1, tolerance = 1e-13)
})

test_that("double-typed alternative codes give the integer result", {
  d <- mxp_data(60, c(2L, 3L, 4L, 2L, 3L), J = 4L, ioo = TRUE)
  dd <- d
  dd$alt_idx <- as.double(d$alt_idx)
  expect_identical(mxp_predict(dd), mxp_predict(d))
  expect_identical(mxp_logsum(dd), mxp_logsum(d))
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)  # one thread: sums over situations in a fixed order
  expect_identical(mxp_shares(dd), mxp_shares(d))
  target <- as.numeric(mxp_shares(d))
  expect_identical(mxp_blp(dd, target, rep(0, 4L)), mxp_blp(d, target, rep(0, 4L)))
  expect_identical(mxp_elas(dd, 1L, TRUE), mxp_elas(d, 1L, TRUE))
  expect_identical(mxp_dr(dd, 2L, FALSE), mxp_dr(d, 2L, FALSE))
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
  blp_k <- function(d, ...) mxp_blp(d, rep(0.25, 4L), rep(0, 4L), ...)
  elas_k <- function(d, ...) mxp_elas(d, 1L, TRUE, ...)
  dr_k <- function(d, ...) mxp_dr(d, 2L, FALSE, ...)
  for (k in list(mxp_predict, mxp_logsum, mxp_shares, blp_k, elas_k, dr_k)) {
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
      if (!identical(k, blp_k)) {  # BLP sizes its ASCs from the codes
        expect_error(call_k(replace(d$alt_idx, 3L, 5L)),
                     "Theta's delta (ASC) block implies 4 alternatives but alt_idx references alternative 5.",
                     fixed = TRUE)
      }
    }
  }
})

test_that("prediction kernels check the design, the draws and W", {
  d <- mxp_data(91, c(2L, 3L, 4L, 2L), J = 4L)
  blp_k <- function(d, ...) mxp_blp(d, rep(0.25, 4L), rep(0, 4L), ...)
  elas_k <- function(d, ...) mxp_elas(d, 1L, TRUE, ...)
  dr_k <- function(d, ...) mxp_dr(d, 2L, FALSE, ...)
  for (k in list(mxp_predict, mxp_logsum, mxp_shares, blp_k, elas_k, dr_k)) {
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
    if (!identical(k, blp_k)) {
      e <- d; e$theta <- e$theta[1:3]
      expect_error(k(e), "Theta vector too short")
    }
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
  blp_w <- function(d, ...) mxp_blp(d, rep(1 / 3, 3L), rep(0, 3L), ...)
  elas_w <- function(d, ...) mxp_elas(d, 1L, FALSE, ...)
  dr_w <- function(d, ...) mxp_dr(d, 1L, TRUE, ...)
  for (k in list(mxp_predict, mxp_logsum, mxp_shares, blp_w, elas_w, dr_w)) {
    expect_error(k(wide, gen_seed = 1L, gen_scramble = 1L, gen_S = 2L),
                 "K_w exceeds the primes table size (128)", fixed = TRUE)
  }
})

test_that("shares and the BLP contraction keep their order of checks", {
  d <- mxp_data(93, c(2L, 3L, 4L, 2L), J = 4L)
  target <- rep(0.25, 4L)
  expect_error(mxp_blp(d, rep(0.2, 5L), rep(0, 4L)),
               "target_shares must have length 4", fixed = TRUE)
  expect_error(mxp_blp(d, target, rep(0, 6L)),
               "delta must have length 4 (full) or 3 (free, with baseline omitted).",
               fixed = TRUE)
  # The weights are checked before X beta is formed: a beta of the wrong
  # length stops there only when the weights are fine
  e <- d; e$weights[] <- 0
  expect_error(mxp_blp(e, target, rep(0, 4L), beta = c(0.1, 0.2, 0.3)),
               "Sum of weights must be positive", fixed = TRUE)
  expect_error(mxp_blp(d, target, rep(0, 4L), beta = c(0.1, 0.2, 0.3)),
               "incompatible matrix dimensions", fixed = TRUE)
  expect_error(mxp_shares(e), "Sum of weights must be positive", fixed = TRUE)
  # Without ASCs the largest code numbers the shares, and an empty design
  # keeps Armadillo's error
  n <- mxp_data(94, rep(5L, 6), J = 7L, use_asc = FALSE, codes = 3:7)
  expect_length(mxp_shares(n), 7L)
  expect_identical(as.numeric(mxp_shares(n))[1:2], c(0, 0))
  n$X <- n$X[0, , drop = FALSE]; n$W <- n$W[0, , drop = FALSE]
  n$alt_idx <- integer(0); n$M <- integer(0); n$weights <- numeric(0)
  n$eta <- n$eta[, , 0, drop = FALSE]
  expect_error(mxp_shares(n), "max(): object has no elements", fixed = TRUE)
  expect_error(mxp_blp(n, numeric(0), numeric(0)),
               "max(): object has no elements", fixed = TRUE)
  # The largest code, 2^31 - 1, plus the outside option overflows an int
  big <- mxp_data(95, rep(3L, 4), J = 5L, use_asc = FALSE, ioo = TRUE)
  big$alt_idx[2] <- .Machine$integer.max
  msg <- "alt_idx references alternative 2147483647, which with the outside option is more alternatives than an int can count."
  expect_error(mxp_shares(big), msg, fixed = TRUE)
  expect_error(mxp_elas(big, 1L, FALSE), msg, fixed = TRUE)
  expect_error(mxp_dr(big, 1L, TRUE), msg, fixed = TRUE)
})

test_that("elasticities and diversion ratios match a brute-force oracle", {
  set.seed(6)
  M <- sample(2:5, 15, replace = TRUE)
  for (ioo in c(FALSE, TRUE)) {
    for (w_type in c("row", "alt")) {
      d <- mxp_data(25 + ioo, M, J = 5L, K_w = 2L, ioo = ioo, w_type = w_type,
                    rc_dist = c(0L, 1L), rc_corr = TRUE)
      for (random in c(FALSE, TRUE)) {
        var <- if (random) 2L else 1L
        o <- mxp_elas_dr_oracle(d, var, random)
        expect_equal(unname(mxp_elas(d, var, random)), o$elas, tolerance = 1e-12)
        expect_equal(unname(mxp_dr(d, var, random)), o$dr, tolerance = 1e-12)
        g <- mxp_gen_cube(length(d$M), 6L, 2L, 21L, 1L)
        og <- mxp_elas_dr_oracle(d, var, random, eta = g)
        expect_equal(unname(mxp_elas(d, var, random, gen_seed = 21L,
                                     gen_scramble = 1L, gen_S = 6L)),
                     og$elas, tolerance = 1e-12)
        expect_equal(unname(mxp_dr(d, var, random, gen_seed = 21L,
                                   gen_scramble = 1L, gen_S = 6L)),
                     og$dr, tolerance = 1e-12)
      }
    }
  }
  # Without ASCs or random-coefficient means, and with codes starting at 3
  d <- mxp_data(20, rep(3L, 9), J = 6L, K_w = 3L, use_asc = FALSE,
                rc_mean = FALSE, codes = 3:6)
  mxp_expect_elas_dr_oracle(d)
  d$ioo <- TRUE
  mxp_expect_elas_dr_oracle(d)
})

test_that("elasticities and diversion ratios of no situations or of one", {
  d <- mxp_data(97, c(2L, 3L), J = 3L, ioo = TRUE)
  e <- d
  e$X <- e$X[0, , drop = FALSE]; e$W <- e$W[0, , drop = FALSE]
  e$alt_idx <- integer(0); e$M <- integer(0); e$weights <- numeric(0)
  e$eta <- e$eta[, , 0, drop = FALSE]
  # With ASCs the alternatives are known: zero matrices over all of them
  expect_identical(mxp_elas(e, 1L, FALSE), matrix(0, 4L, 4L))
  expect_identical(mxp_dr(e, 1L, TRUE), matrix(0, 4L, 4L))
  # Without, the largest code numbers them, and there is none
  e$use_asc <- FALSE
  e$theta <- e$theta[seq_len(ncol(e$X) + 2L * ncol(e$W))]
  expect_error(mxp_elas(e, 1L, FALSE), "max(): object has no elements",
               fixed = TRUE)
  expect_error(mxp_dr(e, 1L, TRUE), "max(): object has no elements",
               fixed = TRUE)
  # One situation runs on two threads, one of them idle; its accumulators add
  # zeros
  o <- d
  o$X <- o$X[1:2, , drop = FALSE]; o$W <- o$W[1:2, , drop = FALSE]
  o$alt_idx <- o$alt_idx[1:2]; o$M <- o$M[1]; o$weights <- o$weights[1]
  o$eta <- o$eta[, , 1, drop = FALSE]
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  e1 <- mxp_elas(o, 1L, FALSE); r1 <- mxp_dr(o, 2L, TRUE); s1 <- mxp_shares(o)
  set_num_threads(2L)
  expect_identical(mxp_elas(o, 1L, FALSE), e1)
  expect_identical(mxp_dr(o, 2L, TRUE), r1)
  expect_identical(mxp_shares(o), s1)
})

test_that("elasticities and diversion ratios check their variable first", {
  d <- mxp_data(96, c(2L, 3L, 4L), J = 4L)
  expect_error(mxp_elas(d, 3L, FALSE),
               "elast_var_idx (3) is out of bounds for X matrix (K_x=2).", fixed = TRUE)
  expect_error(mxp_dr(d, 0L, TRUE),
               "elast_var_idx (0) is out of bounds for W matrix (K_w=2).", fixed = TRUE)
  e <- d; e$theta <- e$theta[1:2]  # a short theta is reported only afterwards
  expect_error(mxp_elas(e, 3L, FALSE), "out of bounds for X matrix", fixed = TRUE)
  e <- d; e$X <- e$X[, 0, drop = FALSE]
  expect_error(mxp_dr(e, 1L, FALSE),
               "the model has no fixed coefficients (K_x = 0)", fixed = TRUE)
  # The unused choice_idx is never converted (one thread: sums over
  # situations in a fixed order)
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  expect_identical(
    choicer:::mxl_elasticities_parallel(d$theta, d$X, d$W, d$alt_idx, NULL,
                                        d$M, d$weights, d$eta, d$rc_dist, 1L,
                                        FALSE, d$rc_corr, d$rc_mean, d$use_asc,
                                        d$ioo),
    mxp_elas(d, 1L, FALSE))
})
