# Brute-force reference for run_mxlogit(concentrate_out_asc = "blp"): the
# share equations on the estimation draws, written from the estimand alone
# with the panel oracle's parameter unpacking (helper-mxl-panel.R).
#
# Decision maker u uses draw block eta[, , u]; at draw s its coefficients are
# mxlp_unit_draws()'s c_us, and each of its situations t has logit
# probabilities P_ts (outside option, utility 0, when present). The share
# equations weight the draws equally:
#   counts_j = sum_t w_t (1/S) sum_s P_tsj       (inside alternatives)
#   Q_usj    = sum_{t in u} P_tsj                (a unit's counts at draw s)
# and their Jacobian on the free constants is what
# mxl_blp_jacobian_parallel() returns (here by numerical differentiation of
# counts_j), K = sum_u w_u Cov_s(Q_us) its variance mode's.

# counts (J) and, with `per_draw`, Q (U x J x S)
blpr_counts <- function(theta, fx, Ti = fx[["Ti"]], eta = fx$eta,
                        weights = fx$weights, per_draw = FALSE) {
  par <- mxlp_unpack(theta, fx)
  N <- length(fx$M)
  off <- mxlp_unit_offsets(Ti, N)
  U <- length(off) - 1L
  S <- dim(eta)[2L]
  J <- fx$J
  row_off <- c(0L, cumsum(fx$M))
  alt_level_W <- nrow(fx$W) != length(fx$alt_idx)
  counts <- numeric(J)
  Q <- if (per_draw) array(0, c(U, J, S))
  for (u in seq_len(U)) {
    for (s in seq_len(S)) {
      gamma <- as.numeric(par$L %*% eta[, s, u])
      c_s <- par$mu_tilde + ifelse(fx$rc_dist == 1L, exp(gamma), gamma)
      for (t in (off[u] + 1L):off[u + 1L]) {
        rows <- (row_off[t] + 1L):row_off[t + 1L]
        a <- fx$alt_idx[rows]
        Wr <- if (alt_level_W) fx$W[a, , drop = FALSE] else fx$W[rows, , drop = FALSE]
        V <- as.numeric(fx$X[rows, , drop = FALSE] %*% par$beta + Wr %*% c_s) +
          par$delta[a]
        v_max <- max(c(V, if (fx$include_outside_option) 0))
        e <- exp(V - v_max)
        P <- e / (sum(e) + if (fx$include_outside_option) exp(-v_max) else 0)
        for (k in seq_along(a)) {
          counts[a[k]] <- counts[a[k]] + weights[t] * P[k] / S
          if (per_draw) Q[u, a[k], s] <- Q[u, a[k], s] + P[k]
        }
      }
    }
  }
  if (per_draw) list(counts = counts, Q = Q) else counts
}

# The free constants' alternatives
blpr_free <- function(fx) if (fx$include_outside_option) seq_len(fx$J) else 2:fx$J

# K = sum_u w_u Cov_s(Q_us) on the free alternatives
blpr_K <- function(theta, fx, Ti = fx[["Ti"]], eta = fx$eta, weights = fx$weights) {
  Q <- blpr_counts(theta, fx, Ti, eta, weights, per_draw = TRUE)$Q
  fr <- blpr_free(fx)
  off <- mxlp_unit_offsets(Ti, length(fx$M))
  w_u <- weights[off[-length(off)] + 1L]
  K <- matrix(0, length(fr), length(fr))
  for (u in seq_len(dim(Q)[1L])) {
    q <- Q[u, fr, , drop = TRUE]
    q <- matrix(q, length(fr))
    qc <- q - rowMeans(q)
    K <- K + w_u[u] * tcrossprod(qc) / ncol(q)
  }
  K
}

# Call mxl_blp_jacobian_parallel() on a fixture
blpr_call <- function(fx, theta = fx$theta, weights = fx$weights, Ti = fx[["Ti"]],
                      eta = fx$eta, ...) {
  choicer:::mxl_blp_jacobian_parallel(
    theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, weights, eta,
    fx$rc_dist, fx$rc_correlation, fx$rc_mean, TRUE, fx$include_outside_option,
    Ti = Ti, ...)
}

# --- The variance (run_mxlogit(concentrate_out_asc = "blp")) ----------------
#
# The estimate solves the stacked equations sum_u w_u (s_c,u + D' s_d,u) = 0
# (the projected score; D = -C_d^{-1} C_c held fixed) and sum_u w_u c_u = 0
# (the share residuals), so its variance is G^{-1} Omega G^{-T} with G their
# Jacobian (here by numDeriv) and Omega their sums' variance.

# Each situation's probabilities averaged over its unit's draws (N x J, zero
# where an alternative is not offered) and its choice (N x J)
blpr_pbar <- function(theta, fx, Ti = fx[["Ti"]], eta = fx$eta) {
  par <- mxlp_unpack(theta, fx)
  N <- length(fx$M)
  off <- mxlp_unit_offsets(Ti, N)
  S <- dim(eta)[2L]
  row_off <- c(0L, cumsum(fx$M))
  alt_level_W <- nrow(fx$W) != length(fx$alt_idx)
  Pbar <- matrix(0, N, fx$J)
  Y <- matrix(0, N, fx$J)
  for (u in seq_len(length(off) - 1L)) {
    for (s in seq_len(S)) {
      gamma <- as.numeric(par$L %*% eta[, s, u])
      c_s <- par$mu_tilde + ifelse(fx$rc_dist == 1L, exp(gamma), gamma)
      for (t in (off[u] + 1L):off[u + 1L]) {
        rows <- (row_off[t] + 1L):row_off[t + 1L]
        a <- fx$alt_idx[rows]
        Wr <- if (alt_level_W) fx$W[a, , drop = FALSE] else fx$W[rows, , drop = FALSE]
        V <- as.numeric(fx$X[rows, , drop = FALSE] %*% par$beta + Wr %*% c_s) +
          par$delta[a]
        v_max <- max(c(V, if (fx$include_outside_option) 0))
        e <- exp(V - v_max)
        P <- e / (sum(e) + if (fx$include_outside_option) exp(-v_max) else 0)
        # row by row: a situation may list an alternative twice
        for (k in seq_along(a)) Pbar[t, a[k]] <- Pbar[t, a[k]] + P[k] / S
      }
    }
  }
  ch <- fx$choice_idx > 0L
  Y[cbind(which(ch), fx$alt_idx[row_off[which(ch)] + fx$choice_idx[ch]])] <- 1
  list(Pbar = Pbar, Y = Y)
}

# Each situation's market (the market of its first alternative)
blpr_market_of_situations <- function(fx, market_of_alt) {
  market_of_alt[fx$alt_idx[c(0L, cumsum(fx$M))[seq_along(fx$M)] + 1L]]
}

# The units' share residuals on the free constants (U x J_d): for the
# sample's shares (tg NULL), c_u = sum_{t in u} (Pbar_t - y_t); for supplied
# targets (tg$s, each free constant's share within its market, and
# tg$market_of_alt), c_u = sum_{t in u} (Pbar_t - s_m(t)) over the free
# constants of t's market
blpr_residuals <- function(theta, fx, tg = NULL, Ti = fx[["Ti"]], eta = fx$eta) {
  pb <- blpr_pbar(theta, fx, Ti, eta)
  fr <- blpr_free(fx)
  off <- mxlp_unit_offsets(Ti, length(fx$M))
  unit <- rep(seq_len(length(off) - 1L), diff(off))
  if (is.null(tg)) {
    r <- pb$Pbar[, fr, drop = FALSE] - pb$Y[, fr, drop = FALSE]
  } else {
    mk <- blpr_market_of_situations(fx, tg$market_of_alt)
    tgt <- t(vapply(mk, function(m) ifelse(tg$market_of_alt[fr] == m, tg$s, 0),
                    numeric(length(fr))))
    r <- pb$Pbar[, fr, drop = FALSE] - matrix(tgt, ncol = length(fr))
  }
  unname(rowsum(r, unit, reorder = FALSE))
}

# The stacked equations at theta, D held fixed (theta_c first)
blpr_stacked <- function(theta, fx, D, tg = NULL) {
  K_c <- ncol(D)
  g <- -choicer:::mxl_loglik_gradient_parallel(
    theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$weights, fx$eta,
    fx$rc_dist, fx$rc_correlation, fx$rc_mean, TRUE, fx$include_outside_option,
    Ti = fx[["Ti"]])$gradient
  c(g[seq_len(K_c)] + as.numeric(crossprod(D, g[-seq_len(K_c)])),
    colSums(mxlp_unit_weights(fx) * blpr_residuals(theta, fx, tg)))
}

# G^{-1} Omega G^{-T} at fx$theta. Omega: "hessian", the model's (It and
# C_d + K, or sum_u w_u c_u c_u' for supplied targets); "sandwich", sum_u
# w_u^2 v_u v_u' from the units' scores and residuals; "cluster", by
# cluster (one code per unit). Vs and W_free: the targets' covariance and
# their markets' weights.
blpr_vcov_bf <- function(fx, type, tg = NULL, cluster = NULL, Vs = NULL,
                         W_free = NULL) {
  theta <- fx$theta
  n <- length(theta)
  J_d <- length(blpr_free(fx))
  ci <- seq_len(n - J_d)
  di <- n - J_d + seq_len(J_d)
  jac <- blpr_call(fx, variance = TRUE)
  Cd <- jac$jacobian[di, di]
  D <- -solve(Cd, jac$jacobian[di, ci])
  G <- numDeriv::jacobian(function(th) blpr_stacked(th, fx, D, tg), theta)
  w_u <- mxlp_unit_weights(fx)
  cu <- blpr_residuals(theta, fx, tg)
  if (type == "hessian") {
    I <- choicer:::mxl_hessian_parallel(
      theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$weights, fx$eta,
      fx$rc_dist, fx$rc_correlation, fx$rc_mean, TRUE, fx$include_outside_option,
      Ti = fx[["Ti"]])
    It <- I[ci, ci] + I[ci, di] %*% D + t(D) %*% I[di, ci] + t(D) %*% I[di, di] %*% D
    Omega <- matrix(0, n, n)
    Omega[ci, ci] <- It
    Omega[di, di] <- if (is.null(tg)) Cd + jac$K else crossprod(sqrt(w_u) * cu)
  } else {
    sc <- choicer:::mxl_scores_parallel(
      theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M, fx$eta, fx$rc_dist,
      fx$rc_correlation, fx$rc_mean, TRUE, fx$include_outside_option,
      Ti = fx[["Ti"]])
    v <- cbind(sc[, ci, drop = FALSE] + sc[, di, drop = FALSE] %*% D, cu)
    Omega <- if (type == "sandwich") {
      crossprod(w_u * v)
    } else {
      crossprod(rowsum(w_u * v, cluster, reorder = FALSE))
    }
  }
  if (!is.null(Vs)) Omega[di, di] <- Omega[di, di] + (W_free %o% W_free) * Vs
  Gi <- solve(G)
  V <- Gi %*% Omega %*% t(Gi)
  (V + t(V)) / 2
}

# The variance's block formulas with solve() throughout, from the
# information I (Hessian or BHHH), the Jacobian blocks C_c (J_d x K_c) and
# C_d, and Omega's blocks
blpr_assemble <- function(I, Cc, Cd, Om11, Om12, Omcc) {
  p <- ncol(Cc)
  ci <- seq_len(p)
  di <- p + seq_len(nrow(Cd))
  Cdi <- solve(Cd)
  D <- -Cdi %*% Cc
  It <- I[ci, ci] + I[ci, di] %*% D + t(D) %*% I[di, ci] + t(D) %*% I[di, di] %*% D
  R <- (I[ci, di] + t(D) %*% I[di, di]) %*% Cdi
  Iti <- solve(It)
  Vc <- Iti %*% (Om11 + R %*% t(Om12) + Om12 %*% t(R) + R %*% Omcc %*% t(R)) %*% Iti
  X <- Iti %*% (Om12 + R %*% Omcc) %*% Cdi
  Vcd <- Vc %*% t(D) - X
  Vd <- D %*% Vc %*% t(D) - D %*% X - t(X) %*% t(D) + Cdi %*% Omcc %*% Cdi
  V <- rbind(cbind(Vc, Vcd), cbind(t(Vcd), Vd))
  (V + t(V)) / 2
}

# Two markets: two panel fixtures with an outside option and three
# alternatives each, stacked, the second's alternatives coded 4 to 6, so
# that no choice situation offers both markets' alternatives and every
# decision maker lies in one market. `...` goes to mxlp_fixture().
blpr_two_markets <- function(seed, ...) {
  a <- mxlp_fixture("market 1", seed, J = 3L, include_outside_option = TRUE, ...)
  b <- mxlp_fixture("market 2", seed + 1, J = 3L, include_outside_option = TRUE, ...)
  K_c <- length(a$theta) - a$J
  fx <- a
  fx$name <- "two markets"
  fx$X <- rbind(a$X, b$X)
  fx$W <- rbind(a$W, b$W)
  fx$alt_idx <- c(a$alt_idx, b$alt_idx + a$J)
  fx$choice_idx <- c(a$choice_idx, b$choice_idx)
  fx$M <- c(a$M, b$M)
  fx$weights <- c(a$weights, b$weights)
  fx$Ti <- c(a$Ti, b$Ti)
  fx$U <- a$U + b$U
  fx$N <- a$N + b$N
  fx$J <- a$J + b$J
  fx$eta <- get_halton_normals(a$S, fx$U, a$K_w)
  fx$theta <- c(a$theta, b$theta[-seq_len(K_c)])
  fx$market_of_alt <- rep(1:2, each = 3L)
  fx$unit_market <- rep(1:2, c(a$U, b$U))
  fx
}
