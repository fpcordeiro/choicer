# A reference for the mixed logit Hessian written from the estimand alone
# (the per-draw formulas of the math article's sections 4.1-4.3 and the Louis
# identity of section 4.5), never from the C++ kernel: for each decision
# maker u and draw s, every alternative's utility gradient z and second
# derivative H_V as dense n_params vectors and matrices, the per-draw score
# g_us and Hessian H_us summed over u's situations, then
#   H_u = sum_s omega_us (H_us + (g_us - g_bar)(g_us - g_bar)'),
# and the negated sum -sum_u w_u H_u, which is what mxl_hessian_parallel()
# returns. It forms what the kernel's Identity D avoids, an outer product per
# alternative and draw, so it is for small fixtures only. The layouts,
# parameters and unpacking are helper-mxl-panel.R's (mxlp_unpack()).
mxlh_reference <- function(theta, fx, eta = fx$eta, weights = fx$weights,
                           Ti = fx[["Ti"]]) {
  par <- mxlp_unpack(theta, fx)
  X <- fx$X
  W <- fx$W
  K_x <- ncol(X)
  K_w <- fx$K_w
  n <- length(theta)
  outside <- fx$include_outside_option
  alt_level_W <- nrow(W) != length(fx$alt_idx)
  lognormal <- fx$rc_dist == 1L
  # The Cholesky parameters r = (p, q), in theta's order, with dL/dl and
  # d2L/dl2 (exp(l) on the diagonal, 1 and 0 off it)
  pq <- if (fx$rc_correlation) {
    do.call(rbind, lapply(seq_len(K_w), function(p) cbind(p, seq_len(p))))
  } else {
    cbind(seq_len(K_w), seq_len(K_w))
  }
  dL <- ifelse(pq[, 1] == pq[, 2], diag(par$L)[pq[, 1]], 1)
  d2L <- ifelse(pq[, 1] == pq[, 2], diag(par$L)[pq[, 1]], 0)
  i_mu <- K_x
  i_L <- K_x + if (fx$rc_mean) K_w else 0L
  i_d <- i_L + nrow(pq)
  mu_raw <- if (fx$rc_mean) theta[i_mu + seq_len(K_w)] else rep(0, K_w)
  dmu <- ifelse(lognormal, exp(mu_raw), 1)    # d mu~ / d mu
  d2mu <- ifelse(lognormal, exp(mu_raw), 0)   # d2 mu~ / d mu2
  free <- function(code) {                    # theta index of a row's ASC, or NA
    if (!fx$use_asc) return(NA_integer_)
    k <- if (outside) code else code - 1L
    if (k >= 1L) i_d + k else NA_integer_
  }
  row_off <- c(0L, cumsum(fx$M))
  off <- mxlp_unit_offsets(Ti, length(fx$M))
  U <- length(off) - 1L
  S <- dim(eta)[2L]
  w_u <- weights[off[-(U + 1L)] + 1L]
  H <- matrix(0, n, n)
  for (u in seq_len(U)) {
    lam <- numeric(S)
    g_s <- matrix(0, n, S)
    H_s <- array(0, c(n, n, S))
    for (s in seq_len(S)) {
      e_s <- eta[, s, u]
      gam <- as.numeric(par$L %*% e_s)
      f1 <- ifelse(lognormal, exp(gam), 1)
      f2 <- ifelse(lognormal, exp(gam), 0)
      coef <- par$mu_tilde + ifelse(lognormal, exp(gam), gam)
      for (t in (off[u] + 1L):off[u + 1L]) {
        rows <- (row_off[t] + 1L):row_off[t + 1L]
        m <- length(rows)
        Wt <- if (alt_level_W) W[fx$alt_idx[rows], , drop = FALSE] else W[rows, , drop = FALSE]
        V <- as.numeric(X[rows, , drop = FALSE] %*% par$beta + Wt %*% coef) +
          par$delta[fx$alt_idx[rows]]
        V_all <- if (outside) c(V, 0) else V
        P_all <- exp(V_all - max(V_all))
        P_all <- P_all / sum(P_all)
        P <- P_all[seq_len(m)]
        j <- fx$choice_idx[t]                  # 0: the outside option
        lam[s] <- lam[s] + log(if (j == 0L) P_all[m + 1L] else P[j])
        Z <- matrix(0, m, n)                   # utility gradients, one row per alternative
        if (K_x > 0L) Z[, seq_len(K_x)] <- X[rows, , drop = FALSE]
        if (fx$rc_mean) Z[, i_mu + seq_len(K_w)] <- sweep(Wt, 2L, dmu, "*")
        for (r in seq_len(nrow(pq))) {
          p <- pq[r, 1]; q <- pq[r, 2]
          Z[, i_L + r] <- Wt[, p] * f1[p] * dL[r] * e_s[q]
        }
        for (a in seq_len(m)) {
          k <- free(fx$alt_idx[rows[a]])
          if (!is.na(k)) Z[a, k] <- Z[a, k] + 1
        }
        d <- -P
        if (j > 0L) d[j] <- d[j] + 1
        zbar <- colSums(P * Z)
        Hts <- -crossprod(Z, P * Z) + tcrossprod(zbar)
        for (a in seq_len(m)) {               # second derivatives of utility
          if (fx$rc_mean) {
            for (p in seq_len(K_w)) {
              Hts[i_mu + p, i_mu + p] <- Hts[i_mu + p, i_mu + p] + d[a] * Wt[a, p] * d2mu[p]
            }
          }
          for (r in seq_len(nrow(pq))) {
            for (r2 in seq_len(nrow(pq))) {
              p <- pq[r, 1]
              if (pq[r2, 1] != p) next
              a1 <- dL[r] * e_s[pq[r, 2]]
              a2 <- dL[r2] * e_s[pq[r2, 2]]
              h <- f2[p] * a1 * a2 + if (r == r2) f1[p] * d2L[r] * e_s[pq[r, 2]] else 0
              Hts[i_L + r, i_L + r2] <- Hts[i_L + r, i_L + r2] + d[a] * Wt[a, p] * h
            }
          }
        }
        g_s[, s] <- g_s[, s] + colSums(d * Z)
        H_s[, , s] <- H_s[, , s] + Hts
      }
    }
    omega <- exp(lam - max(lam))
    omega <- omega / sum(omega)
    # A draw of weight zero contributes nothing, though its terms may not be
    # finite (a choice probability of zero): leave it out.
    keep <- which(omega > 0)
    g_bar <- as.numeric(g_s[, keep, drop = FALSE] %*% omega[keep])
    H_u <- matrix(0, n, n)
    for (s in keep) {
      dev <- g_s[, s] - g_bar
      H_u <- H_u + omega[s] * (H_s[, , s] + tcrossprod(dev))
    }
    H <- H - w_u[u] * H_u
  }
  H
}
