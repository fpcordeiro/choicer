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
