# Brute-force oracle and fixtures for the panel mixed logit kernels.
#
# The oracle is written from the estimand alone, never from the C++ kernels, so
# test-mxl-panel-kernels.R compares two independent computations. It is
# deliberately slow: plain nested loops over decision makers, draws, choice
# situations and alternatives.
#
# Prepared (stacked) layout, as produced by prepare_mxl_data():
#   * Situation t owns M[t] consecutive inside rows; alt_idx[r] is the global
#     1-based alternative code of row r. choice_idx[t] is the 1-based position
#     of the chosen row within t, or 0 when the outside option (utility 0) is
#     chosen, which requires include_outside_option = TRUE.
#   * W is row-aligned (sum(M) x K_w) or alternative-level (J x K_w, row j
#     holds alternative j's attributes).
#   * Ti[u] is the number of consecutive situations of decision maker u;
#     Ti = NULL makes every situation its own decision maker.
#
# Parameters, theta = [beta (K_x), mu (K_w, if rc_mean), L block, delta]:
#   * L block: with rc_correlation the lower triangle of L packed row-major
#     (l11, l21, l22, l31, l32, l33, ...), otherwise the K_w diagonal entries;
#     diagonal entries on the log scale, off-diagonal entries raw.
#   * delta (if use_asc): one ASC per inside alternative with the outside
#     option; otherwise alternatives 2..J, with delta_1 = 0.
#
# Estimand. Decision maker u uses draw block eta[, , u] (K_w x S). For draw s,
# gamma_us = L eta[, s, u] and the random coefficient is
#   c_us,k = mu~_k + gamma_us,k        (normal, rc_dist = 0)
#   c_us,k = mu~_k + exp(gamma_us,k)   (shifted log-normal, rc_dist = 1)
# with mu~_k = mu_k (normal) or exp(mu_k) (log-normal), and mu~ = 0 when
# rc_mean = FALSE. The utility of inside row r is
#   V_r = X_r beta + W_r c_us + delta_alt(r).
# Then
#   lambda_us = sum_{t in u} log P_ts(chosen_t)
#   l_u       = log((1 / S) sum_s exp(lambda_us))
#   loglik    = sum_u w_u l_u,  w_u = weight of u's first situation
#   omega_us  = exp(lambda_us) / sum_s exp(lambda_us)
#   E[c_u | y_u]  ~ sum_s omega_us c_us
#   SD[c_u | y_u] ~ sqrt(sum_s omega_us (c_us - E[c_u | y_u])^2)

# --- Parameter unpacking -----------------------------------------------------

# Lower-triangular Cholesky factor L from its packed parameter block.
mxlp_chol <- function(L_params, K_w, rc_correlation) {
  L <- matrix(0, K_w, K_w)
  if (rc_correlation) {
    stopifnot(length(L_params) == K_w * (K_w + 1L) / 2L)
    k <- 0L
    for (i in seq_len(K_w)) {
      for (j in seq_len(i)) {
        k <- k + 1L
        L[i, j] <- if (i == j) exp(L_params[k]) else L_params[k]
      }
    }
  } else {
    stopifnot(length(L_params) == K_w)
    for (i in seq_len(K_w)) L[i, i] <- exp(L_params[i])
  }
  L
}

# Split theta into beta, the location shift mu~, L and the full-length delta.
mxlp_unpack <- function(theta, fx) {
  K_x <- ncol(fx$X)
  K_w <- fx$K_w
  J <- fx$J
  pos <- 0L
  take <- function(n) {
    idx <- pos + seq_len(n)
    pos <<- pos + n
    theta[idx]
  }
  beta <- take(K_x)
  mu <- if (fx$rc_mean) take(K_w) else rep(0, K_w)
  L_size <- if (fx$rc_correlation) (K_w * (K_w + 1L)) %/% 2L else K_w
  L <- mxlp_chol(take(L_size), K_w, fx$rc_correlation)
  delta <- rep(0, J)
  if (fx$use_asc) {
    if (fx$include_outside_option) {
      delta <- take(J)
    } else {
      delta[2:J] <- take(J - 1L)
    }
  }
  stopifnot(pos == length(theta))

  mu_tilde <- rep(0, K_w)
  if (fx$rc_mean) {
    for (k in seq_len(K_w)) {
      mu_tilde[k] <- if (fx$rc_dist[k] == 1L) exp(mu[k]) else mu[k]
    }
  }
  list(beta = beta, mu_tilde = mu_tilde, L = L, delta = delta)
}

# --- Oracle ------------------------------------------------------------------

# Situation offsets of the decision makers: unit u owns situations
# (off[u] + 1):off[u + 1].
mxlp_unit_offsets <- function(Ti, N) {
  if (is.null(Ti)) Ti <- rep(1L, N)
  stopifnot(all(Ti >= 1L), sum(Ti) == N)
  c(0L, cumsum(as.integer(Ti)))
}

# Weight of each decision maker: the weight of its first situation.
mxlp_unit_weights <- function(fx, weights = fx$weights, Ti = fx[["Ti"]]) {
  off <- mxlp_unit_offsets(Ti, length(fx$M))
  weights[off[-length(off)] + 1L]
}

# lambda_us (length S) and c_us (K_w x S) for one decision maker u.
mxlp_unit_draws <- function(par, fx, u, off, eta) {
  X <- fx$X
  W <- fx$W
  alt_idx <- fx$alt_idx
  choice_idx <- fx$choice_idx
  rc_dist <- fx$rc_dist
  outside <- fx$include_outside_option
  K_w <- fx$K_w
  S <- dim(eta)[2L]
  row_off <- c(0L, cumsum(fx$M))
  alt_level_W <- nrow(W) != length(alt_idx)
  lambda <- numeric(S)
  coef <- matrix(0, K_w, S)
  for (s in seq_len(S)) {
    gamma <- as.numeric(par$L %*% eta[, s, u])
    c_s <- numeric(K_w)
    for (k in seq_len(K_w)) {
      c_s[k] <- par$mu_tilde[k] +
        (if (rc_dist[k] == 1L) exp(gamma[k]) else gamma[k])
    }
    coef[, s] <- c_s
    for (t in (off[u] + 1L):off[u + 1L]) {
      rows <- (row_off[t] + 1L):row_off[t + 1L]
      V <- numeric(length(rows))
      for (m in seq_along(rows)) {
        r <- rows[m]
        a <- alt_idx[r]
        w_r <- if (alt_level_W) W[a, ] else W[r, ]
        V[m] <- sum(X[r, ] * par$beta) + sum(w_r * c_s) + par$delta[a]
      }
      # Log denominator by log-sum-exp; the outside option has utility 0.
      V_all <- if (outside) c(V, 0) else V
      v_max <- max(V_all)
      log_denom <- v_max + log(sum(exp(V_all - v_max)))
      j <- choice_idx[t]
      if (j == 0L) stopifnot(outside)
      v_chosen <- if (j == 0L) 0 else V[j]
      lambda[s] <- lambda[s] + v_chosen - log_denom
    }
  }
  list(lambda = lambda, coef = coef)
}

# Brute-force panel simulated log-likelihood, posterior draw weights and
# conditional tastes. Returns a list with `loglik` (sum_u w_u l_u),
# `loglik_u` and `w_u` (length U), `lambda` and `omega` (U x S), and the
# conditional-taste `mean` and `sd` (K_w x U).
mxlp_oracle <- function(theta, fx, Ti = fx[["Ti"]], eta = fx$eta,
                        weights = fx$weights) {
  par <- mxlp_unpack(theta, fx)
  N <- length(fx$M)
  off <- mxlp_unit_offsets(Ti, N)
  U <- length(off) - 1L
  S <- dim(eta)[2L]
  K_w <- fx$K_w
  stopifnot(dim(eta)[1L] == K_w, dim(eta)[3L] == U, length(weights) == N)

  w_u <- weights[off[-(U + 1L)] + 1L]
  stopifnot(all(weights == rep(w_u, diff(off))))  # decision-maker weights

  lambda <- matrix(0, U, S)
  omega <- matrix(0, U, S)
  loglik_u <- numeric(U)
  t_mean <- matrix(0, K_w, U)
  t_sd <- matrix(0, K_w, U)
  for (u in seq_len(U)) {
    d <- mxlp_unit_draws(par, fx, u, off, eta)
    lam_max <- max(d$lambda)
    e <- exp(d$lambda - lam_max)
    loglik_u[u] <- lam_max + log(sum(e) / S)
    lambda[u, ] <- d$lambda
    omega[u, ] <- e / sum(e)
    for (k in seq_len(K_w)) {
      t_mean[k, u] <- sum(omega[u, ] * d$coef[k, ])
      t_sd[k, u] <- sqrt(sum(omega[u, ] * (d$coef[k, ] - t_mean[k, u])^2))
    }
  }
  list(loglik = sum(w_u * loglik_u), loglik_u = loglik_u, w_u = w_u,
       lambda = lambda, omega = omega, mean = t_mean, sd = t_sd)
}

# Weight-free log-likelihood l_u of one decision maker (numDeriv target for
# the score rows).
mxlp_oracle_unit_loglik <- function(theta, fx, u, Ti = fx[["Ti"]],
                                    eta = fx$eta) {
  par <- mxlp_unpack(theta, fx)
  off <- mxlp_unit_offsets(Ti, length(fx$M))
  lambda <- mxlp_unit_draws(par, fx, u, off, eta)$lambda
  lam_max <- max(lambda)
  lam_max + log(mean(exp(lambda - lam_max)))
}

# --- Fixtures ----------------------------------------------------------------

# Small prepared panel in the stacked layout, plus a moderate theta.
#
# U decision makers with unbalanced T_u in 1..4 (every value occurs), choice
# sets of varying size over J global alternatives (situation 1 sees all of
# them), and, with the outside option, real outside choosers: every situation
# of one multi-situation decision maker plus about a quarter of the rest.
# Weights are 1 (`weight_type = "unit"`) or non-uniform but constant within
# each decision maker (`"person"`). Draws are get_halton_normals(S, U, K_w),
# one K_w x S block per decision maker.
mxlp_fixture <- function(name, seed, rc_dist = c(0L, 0L),
                         rc_correlation = FALSE, rc_mean = FALSE,
                         use_asc = TRUE, include_outside_option = FALSE,
                         W_layout = c("row", "alt"),
                         weight_type = c("unit", "person"),
                         U = 14L, J = 4L, K_x = 2L, S = 12L) {
  W_layout <- match.arg(W_layout)
  weight_type <- match.arg(weight_type)
  rc_dist <- as.integer(rc_dist)
  K_w <- length(rc_dist)
  stopifnot(U >= 5L, J >= 3L)
  set.seed(seed)

  Ti <- sample(c(1:4, sample.int(4L, U - 4L, replace = TRUE)))
  N <- sum(Ti)
  unit_of <- rep(seq_len(U), Ti)

  min_M <- if (include_outside_option) 1L else 2L
  alts <- vector("list", N)
  for (t in seq_len(N)) {
    size <- if (t == 1L) J else min_M - 1L + sample.int(J - min_M + 1L, 1L)
    alts[[t]] <- sort(sample.int(J, size))
  }
  M <- lengths(alts)
  alt_idx <- unlist(alts)
  n_rows <- sum(M)

  choice_idx <- vapply(M, function(m) sample.int(m, 1L), integer(1L))
  if (include_outside_option) {
    u_out <- which(Ti >= 2L & seq_len(U) > 1L)[1L]
    choice_idx[unit_of == u_out] <- 0L
    choice_idx[stats::runif(N) < 0.25] <- 0L
  }

  X <- matrix(stats::rnorm(n_rows * K_x), n_rows, K_x)
  W <- if (W_layout == "row") {
    matrix(stats::rnorm(n_rows * K_w), n_rows, K_w)
  } else {
    matrix(stats::rnorm(J * K_w), J, K_w)
  }

  w_u <- if (weight_type == "person") stats::runif(U, 0.5, 2) else rep(1, U)
  weights <- rep(w_u, Ti)

  L_params <- numeric(0)
  if (rc_correlation) {
    for (i in seq_len(K_w)) {
      for (j in seq_len(i)) {
        L_params <- c(L_params, if (i == j) log(stats::runif(1L, 0.3, 0.8))
                                else stats::runif(1L, -0.3, 0.3))
      }
    }
  } else {
    L_params <- log(stats::runif(K_w, 0.3, 0.8))
  }
  n_delta <- if (!use_asc) 0L else if (include_outside_option) J else J - 1L
  theta <- c(stats::rnorm(K_x, sd = 0.3),
             if (rc_mean) stats::rnorm(K_w, sd = 0.3),
             L_params,
             stats::rnorm(n_delta, sd = 0.3))

  # Units probed by the score-row checks: the first, the interior unit with
  # the most situations, and the last. At least three units have T_u >= 2,
  # so one of them is interior.
  interior <- 2:(U - 1L)
  u_multi <- interior[which.max(Ti[interior])]

  list(
    name = name, X = X, W = W, alt_idx = as.integer(alt_idx),
    choice_idx = as.integer(choice_idx), M = as.integer(M),
    weights = weights, Ti = as.integer(Ti), eta = get_halton_normals(S, U, K_w),
    theta = theta, rc_dist = rc_dist, rc_correlation = rc_correlation,
    rc_mean = rc_mean, use_asc = use_asc,
    include_outside_option = include_outside_option,
    K_w = K_w, J = J, S = S, U = U, N = N,
    probe_units = c(1L, u_multi, U)
  )
}

# A long panel whose likelihood underflows double precision: decision maker 2
# makes T_long choices, each of the alternative with the lowest fixed utility,
# so lambda_us is about -1000, past the exp() underflow threshold (about -745),
# and exp(lambda_us) is exactly 0 for every draw. Random-coefficient covariates
# are small, so the posterior draw weights stay spread out and the likelihood
# stays smooth at numDeriv's step sizes.
mxlp_long_panel_fixture <- function(seed = 211, T_long = 150L, S = 12L) {
  set.seed(seed)
  J <- 3L
  K_w <- 2L
  Ti <- c(2L, as.integer(T_long), 1L)
  U <- length(Ti)
  N <- sum(Ti)
  M <- rep(J, N)
  alt_idx <- rep(seq_len(J), N)
  n_rows <- N * J
  X <- matrix(stats::rnorm(n_rows, sd = 2), n_rows, 1L)
  W <- matrix(stats::rnorm(n_rows * K_w, sd = 0.25), n_rows, K_w)
  # beta, L (l11, l21, l22), delta_2, delta_3
  theta <- c(2, log(0.5), 0.2, log(0.4), 0.3, -0.3)
  choice_idx <- apply(matrix(X[, 1L], nrow = J), 2L, which.min)
  list(
    name = "long panel", X = X, W = W, alt_idx = as.integer(alt_idx),
    choice_idx = as.integer(choice_idx), M = as.integer(M),
    weights = rep(1, N), Ti = Ti, eta = get_halton_normals(S, U, K_w),
    theta = theta, rc_dist = c(0L, 1L), rc_correlation = TRUE,
    rc_mean = FALSE, use_asc = TRUE, include_outside_option = FALSE,
    K_w = K_w, J = J, S = S, U = U, N = N, probe_units = seq_len(U)
  )
}

# A small panel (or, with panel = FALSE, its cross-section) in which the chosen
# alternative of one situation sits `gap` utils below its competitors, through
# the first fixed covariate. The situation is the first one of an interior
# multi-situation decision maker; `edge_situation` and `edge_unit` locate it.
mxlp_edge_fixture <- function(gap, panel = TRUE) {
  fx <- mxlp_fixture("edge", seed = 301, rc_dist = c(0L, 1L),
                     rc_correlation = TRUE, rc_mean = TRUE, use_asc = TRUE,
                     include_outside_option = TRUE, W_layout = "row",
                     weight_type = "person", U = 6L, J = 3L, S = 10L)
  fx$theta[1L] <- 0.5
  u <- fx$probe_units[2L]
  t <- sum(fx$Ti[seq_len(u - 1L)]) + 1L
  if (!panel) fx <- mxlp_cross_section(fx)
  row_off <- c(0L, cumsum(fx$M))
  rows <- (row_off[t] + 1L):row_off[t + 1L]
  if (fx$choice_idx[t] == 0L) fx$choice_idx[t] <- 1L
  fx$X[rows, 1L] <- 0
  fx$X[rows[fx$choice_idx[t]], 1L] <- gap / fx$theta[1L]
  fx$name <- sprintf("gap %g, %s", gap, if (panel) "panel" else "cross-section")
  fx$edge_situation <- t
  fx$edge_unit <- if (panel) u else t
  fx
}

# log P(chosen) of situation t at every draw of its decision maker's block.
mxlp_situation_logp <- function(theta, fx, t) {
  off <- mxlp_unit_offsets(fx[["Ti"]], length(fx$M))
  u <- max(which(off < t))  # the decision maker that owns situation t
  mxlp_unit_draws(mxlp_unpack(theta, fx), fx, 1L, c(t - 1L, t),
                  fx$eta[, , u, drop = FALSE])$lambda
}

# The fixture with decision maker u removed: its situations, their rows and
# its draw block.
mxlp_drop_unit <- function(fx, u) {
  off <- mxlp_unit_offsets(fx[["Ti"]], length(fx$M))
  keep_t <- setdiff(seq_along(fx$M), (off[u] + 1L):off[u + 1L])
  row_off <- c(0L, cumsum(fx$M))
  keep_r <- unlist(lapply(keep_t, function(t) {
    (row_off[t] + 1L):row_off[t + 1L]
  }))
  if (nrow(fx$W) == length(fx$alt_idx)) fx$W <- fx$W[keep_r, , drop = FALSE]
  fx$X <- fx$X[keep_r, , drop = FALSE]
  fx$alt_idx <- fx$alt_idx[keep_r]
  fx$choice_idx <- fx$choice_idx[keep_t]
  fx$M <- fx$M[keep_t]
  fx$weights <- fx$weights[keep_t]
  if (!is.null(fx[["Ti"]])) fx$Ti <- fx$Ti[-u]
  fx$eta <- fx$eta[, , -u, drop = FALSE]
  fx$N <- length(keep_t)
  fx$U <- dim(fx$eta)[3L]
  fx
}

# The same data as a cross-section: every situation its own decision maker,
# one draw block per situation, and situation weights that vary freely.
mxlp_cross_section <- function(fx, seed = 1L) {
  set.seed(seed)
  fx[["Ti"]] <- NULL
  fx$U <- fx$N
  fx$eta <- get_halton_normals(fx$S, fx$N, fx$K_w)
  fx$weights <- stats::runif(fx$N, 0.5, 2)
  fx$probe_units <- c(1L, 2L, fx$N)
  fx
}

# Cell grid for the kernel tests. Every level of every factor appears in at
# least two cells: rc_dist normal / log-normal / mixed, rc_correlation,
# rc_mean, the outside option (with real outside choosers), use_asc,
# row-aligned vs alternative-level W, and unit vs person-constant weights.
# The two K_w = 3 correlated cells pin the row-major Cholesky packing, which
# coincides with column-major packing when K_w = 2.
mxlp_cell_configs <- function() {
  cell <- function(seed, rc_dist, rc_correlation, rc_mean, use_asc,
                   include_outside_option, W_layout, weight_type) {
    list(seed = seed, rc_dist = rc_dist, rc_correlation = rc_correlation,
         rc_mean = rc_mean, use_asc = use_asc,
         include_outside_option = include_outside_option,
         W_layout = W_layout, weight_type = weight_type)
  }
  list(
    cell(101, c(0L, 0L),     FALSE, FALSE, TRUE,  FALSE, "row", "unit"),
    cell(102, c(0L, 0L),     TRUE,  FALSE, TRUE,  TRUE,  "row", "person"),
    cell(103, c(0L, 0L),     FALSE, TRUE,  TRUE,  FALSE, "alt", "person"),
    cell(104, c(1L, 1L),     FALSE, FALSE, FALSE, TRUE,  "row", "unit"),
    cell(105, c(0L, 1L),     TRUE,  TRUE,  FALSE, FALSE, "alt", "unit"),
    cell(106, c(1L, 0L),     FALSE, TRUE,  TRUE,  TRUE,  "alt", "person"),
    cell(107, c(0L, 0L, 1L), TRUE,  TRUE,  FALSE, TRUE,  "row", "person"),
    cell(108, c(0L, 0L, 0L), TRUE,  FALSE, FALSE, FALSE, "row", "unit"),
    cell(109, c(1L, 1L),     TRUE,  TRUE,  TRUE,  TRUE,  "alt", "unit"),
    cell(110, c(0L, 0L),     FALSE, TRUE,  FALSE, FALSE, "row", "person")
  )
}

# Readable cell label, e.g. "n+ln|corr|mean|outside|asc|W=alt|w=person".
mxlp_cell_name <- function(cfg) {
  paste0(
    paste(ifelse(cfg$rc_dist == 1L, "ln", "n"), collapse = "+"),
    if (cfg$rc_correlation) "|corr" else "|diag",
    if (cfg$rc_mean) "|mean" else "",
    if (cfg$include_outside_option) "|outside" else "",
    if (cfg$use_asc) "|asc" else "|noasc",
    "|W=", cfg$W_layout,
    "|w=", cfg$weight_type
  )
}

mxlp_build_cells <- function() {
  lapply(mxlp_cell_configs(), function(cfg) {
    do.call(mxlp_fixture, c(list(name = mxlp_cell_name(cfg)), cfg))
  })
}

# --- Kernel calls and comparisons -------------------------------------------

mxlp_kernels <- c("gradient", "hessian", "bhhh", "scores", "tastes")

# Call one of the five MXL panel kernels on a fixture. `Ti = NULL` omits the
# argument (the kernels' cross-sectional default); `generate = TRUE` switches
# to on-the-fly Halton draws with identity scrambling, which reproduce
# get_halton_normals() exactly.
mxlp_call <- function(kernel, fx, theta = fx$theta, weights = fx$weights,
                      Ti = fx[["Ti"]], eta = fx$eta, generate = FALSE) {
  args <- list(theta = theta, X = fx$X, W = fx$W, alt_idx = fx$alt_idx,
               choice_idx = fx$choice_idx, M = fx$M)
  if (kernel %in% c("gradient", "hessian", "bhhh")) args$weights <- weights
  args$eta_draws <- if (generate) array(0, dim = c(fx$K_w, 0L, 0L)) else eta
  args <- c(args, list(
    rc_dist = fx$rc_dist, rc_correlation = fx$rc_correlation,
    rc_mean = fx$rc_mean, use_asc = fx$use_asc,
    include_outside_option = fx$include_outside_option
  ))
  if (generate) {
    args <- c(args, list(gen_seed = 0L, gen_scramble = 0L, gen_S = fx$S))
  }
  if (!is.null(Ti)) args$Ti <- as.integer(Ti)
  fn <- switch(kernel,
    gradient = mxl_loglik_gradient_parallel,
    hessian  = mxl_hessian_parallel,
    bhhh     = mxl_bhhh_parallel,
    scores   = mxl_scores_parallel,
    tastes   = mxl_conditional_tastes_parallel,
    stop("unknown kernel: ", kernel)
  )
  do.call(fn, args)
}

# Flatten a kernel result (matrix or list of arrays) to one numeric vector.
mxlp_flat <- function(x) {
  if (!is.list(x)) return(as.numeric(x))
  unlist(lapply(x, as.numeric), use.names = FALSE)
}

# Largest elementwise error |a - e| / max(1, |e|): relative for entries above
# one, absolute below. Inf on a length mismatch or any non-finite value, so a
# comparison can never pass by accident.
mxlp_rel_err <- function(actual, expected) {
  actual <- mxlp_flat(actual)
  expected <- mxlp_flat(expected)
  if (length(expected) == 0L || length(actual) != length(expected) ||
      !all(is.finite(actual)) || !all(is.finite(expected))) {
    return(Inf)
  }
  max(abs(actual - expected) / pmax(1, abs(expected)))
}

# expect_lt() on mxlp_rel_err(); `what` names the cell and the check, and the
# failure message reports the error attained.
mxlp_expect_close <- function(actual, expected, tol, what) {
  err <- mxlp_rel_err(actual, expected)
  expect_lt(err, tol, label = sprintf("%s: max rel err %.2e", what, err))
}
