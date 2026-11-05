# ============================================================================
# Q1.3 (plan §5, claim C4): choicer's mixed logit by maximum simulated
# likelihood, S = 4000 generate-mode draws, against the suite's exact
# maximum likelihood on the same data (D4-40, N = 2000). With K_w = 2 random coefficients the
# integral over eta ~ N(0, I_2) is a product Gauss-Hermite rule, exact for
# any transformation of eta up to the rule's error, which the node count
# makes negligible (chosen by the log-likelihood's convergence over a list
# of grids, the result certified on the finest). The reference is maximized
# in R with its analytic gradient from choicer's estimate, by steps in
# choicer's inverse Hessian, which sets the path only (the fixed point is the
# reference's g = 0), and it is certified in its own metric: the decrement
# g' B^-1 g, B the reference's outer product of scores, below 1e-12 (BHHH
# steps alone reach it too, more slowly). Pass: max |theta_ML - theta_MSL| /
# SE < 0.05 over every parameter, SE choicer's Hessian SE. Definitions only.
# Usage: ref_q13() returns check rows (test, check, value, threshold, pass;
# pass NA for information rows).
# ============================================================================

# Gauss-Hermite rule for E f(z), z ~ N(0, 1) (Golub-Welsch): nodes z, weights w
# summing to 1.
gh_normal <- function(Q) {
  Q <- as.integer(Q)
  jac <- matrix(0, Q, Q)
  if (Q > 1L) jac[cbind(1:(Q - 1L), 2:Q)] <- jac[cbind(2:Q, 1:(Q - 1L))] <- sqrt(1:(Q - 1L) / 2)
  e <- eigen(jac, symmetric = TRUE)
  list(z = sqrt(2) * e$values, w = e$vectors[1L, ]^2)
}

# The product rule on eta in R^2, Q[1] nodes on the first coordinate and
# Q[2] (default Q[1]) on the second: nodes (2 x Q1 Q2) and weights.
gh_product2 <- function(Q) {
  Q <- rep_len(as.integer(Q), 2L)
  g1 <- gh_normal(Q[1L])
  g2 <- gh_normal(Q[2L])
  list(z = rbind(rep(g1$z, Q[2L]), rep(g2$z, each = Q[1L])),
       w = rep(g1$w, Q[2L]) * rep(g2$w, each = Q[1L]))
}

# The exact (quadrature) mixed logit at theta for design d (lib/simulate.R's
# design(), cross-section, with an outside option or not), chosen rows
# `chosen` (chosen_rows()), the rule `rule` (gh_product2()) and situation
# weights `w`. what = "ll" gives the log-likelihood; "grad" adds its
# gradient; "scores" adds the N x p matrix of per-situation scores.
#   ll = sum_n w_n log sum_q omega_q P_{n c(n) q},
# P_{nrq} the logit at node q, coefficients rc_coefs() at eta = z_q. With
# R_nq = omega_q P_{ncq} / L_n and Pi_r = sum_q R_{n(r)q} P_rq, the score of a
# coefficient on a column x (X, the constants' indicators) is
# sum_r (y_r - Pi_r) x_r over the situation's rows; of a parameter t moving
# random coefficient a by G_tq at node q, sum_q R_nq G_tq (W_{c,a} -
# sum_r P_rq W_ra).
quad_mxl <- function(theta, d, spec, chosen, rule, w = NULL, what = c("ll", "grad", "scores"),
                     chunk_values = 4e6) {
  what <- match.arg(what)
  spec <- .spec("mxl", spec)
  par <- unpack_theta("mxl", spec, theta)
  N <- .n_sit(d$sit)
  n <- length(d$sit)
  if (is.null(w)) w <- rep(1, N)
  K_w <- ncol(d$W)
  nq <- length(rule$w)
  g <- par$L %*% rule$z                                  # K_w x Q^2
  B <- rc_coefs(par$mu, par$L, rule$z, par$rc_dist, par$rc_mean)
  base <- .utility(d, par)
  inside <- chosen > 0L
  # The nodes in chunks of at most chunk_values row x node values. Each
  # chunk's logit: P (rows x nodes) and Pc, the chosen alternative's (the
  # outside option's for its choosers).
  chunks <- split(seq_len(nq), ceiling(seq_len(nq) / max(1, floor(chunk_values / n))))
  logit <- function(q) {
    V <- base + d$W %*% B[, q, drop = FALSE]
    if (max(V) > 600) stop("quad_mxl(): utilities above 600; exp() would overflow")
    E <- exp(V)
    D <- rowsum(E, d$sit, reorder = TRUE)
    if (par$outside) D <- D + 1
    P <- E / D[d$sit, , drop = FALSE]
    Pc <- if (par$outside) 1 / D else matrix(0, N, length(q))
    Pc[inside, ] <- P[chosen[inside], , drop = FALSE]
    list(P = P, Pc = Pc)
  }
  Pc <- matrix(0, N, nq)
  for (q in chunks) Pc[, q] <- logit(q)$Pc
  Lnq <- sweep(Pc, 2L, rule$w, `*`)
  Ln <- rowSums(Lnq)
  ll <- sum(w * log(Ln))
  if (what == "ll") return(list(ll = ll))

  # Derivative pieces, chunk by chunk: R = Lnq / Ln (N x nodes), RP = R P
  # (rows x nodes) summed into Pi; per random coefficient a, the chosen
  # row's W_a (0 for the outside option) against sum_r P W_ra.
  R <- Lnq / Ln
  pos <- theta_blocks("mxl", spec)
  nm <- theta_names("mxl", spec)
  Gs <- .quad_G(par, g, B, rule$z, pos)                  # list: pos, a (column), G (Q^2)
  Wc <- matrix(0, N, K_w)
  Wc[inside, ] <- d$W[chosen[inside], , drop = FALSE]
  wr <- w[d$sit]
  Pi <- numeric(n)
  S <- matrix(0, N, length(nm), dimnames = list(NULL, nm))
  gs <- numeric(length(nm))                              # random-coefficient parameters
  for (q in chunks) {
    P <- logit(q)$P
    Rq <- R[, q, drop = FALSE]
    RP <- Rq[d$sit, , drop = FALSE] * P
    Pi <- Pi + rowSums(RP)
    for (a in seq_len(K_w)) {
      ts <- Gs[vapply(Gs, `[[`, 1L, "a") == a]
      if (!length(ts)) next
      if (what == "grad") {
        dif <- colSums(w * Rq * Wc[, a]) - colSums(wr * RP * d$W[, a])
        for (t in ts) gs[t$pos] <- gs[t$pos] + sum(t$G[q] * dif)
      } else {
        M <- Rq * (Wc[, a] - rowsum(P * d$W[, a], d$sit, reorder = TRUE))
        for (t in ts) S[, t$pos] <- S[, t$pos] + drop(M %*% t$G[q])
      }
    }
  }
  # Coefficients on X and the constants: (y - Pi) over the situation's rows.
  y <- numeric(n)
  y[chosen[inside]] <- 1
  res <- y - Pi
  free_asc <- if (length(pos$asc)) match(.asc_labels(spec), .alt_levels(spec)) else integer(0)
  if (what == "grad") {
    gr <- gs
    gr[pos$beta] <- colSums(wr * res * d$X)
    if (length(pos$asc)) {
      ga <- numeric(length(.alt_levels(spec)))
      s <- rowsum(wr * res, d$alt, reorder = TRUE)
      ga[as.integer(rownames(s))] <- s[, 1L]
      gr[pos$asc] <- ga[free_asc]
    }
    return(list(ll = ll, grad = stats::setNames(gr, nm)))
  }
  S[, pos$beta] <- rowsum(res * d$X, d$sit, reorder = TRUE)
  if (length(pos$asc)) {
    Sa <- matrix(0, N, length(.alt_levels(spec)))
    Sa[cbind(d$sit, d$alt)] <- res
    S[, pos$asc] <- Sa[, free_asc]
  }
  S <- w * S
  list(ll = ll, grad = colSums(S), scores = S)
}

# For each parameter moving the random coefficients (mu, sigma): its
# position, the coefficient a it moves and G_q = d b_a(z_q) / d t, where b_a
# = f_a(mu_a, g_a), g = L z, f = mu + g (normal) or exp(mu) + exp(g)
# (log-normal; exp(g) without rc_mean), sigma's entry (i, j) of L moving g_i
# by z_j (times L_ii on the log diagonal).
.quad_G <- function(par, g, B, z, pos) {
  K_w <- nrow(par$L)
  ln <- par$rc_dist == 1L
  fprime <- g
  fprime[] <- 1
  fprime[ln, ] <- exp(g[ln, , drop = FALSE])
  out <- list()
  if (length(pos$mu)) {
    for (a in seq_len(K_w)) {
      out[[length(out) + 1L]] <- list(pos = pos$mu[a], a = a,
                                      G = rep(if (ln[a]) exp(par$mu[[a]]) else 1, ncol(z)))
    }
  }
  if (par$correlated) {
    ij <- cbind(rep(seq_len(K_w), seq_len(K_w)), sequence(seq_len(K_w)))
  } else {
    ij <- cbind(seq_len(K_w), seq_len(K_w))
  }
  for (k in seq_len(nrow(ij))) {
    i <- ij[k, 1L]
    j <- ij[k, 2L]
    G <- fprime[i, ] * z[j, ] * (if (i == j) par$L[i, i] else 1)
    out[[length(out) + 1L]] <- list(pos = pos$sigma[k], a = i, G = G)
  }
  out
}

# Maximizes the quadrature log-likelihood from theta by steps M g with step
# halving, M the fixed `metric` (an inverse curvature: choicer's vcov at its
# estimate, which sets only the path, since the fixed point g = 0 is the
# reference's) or, without one, BHHH's (S'S)^-1 at each point; stops when
# the decrement g' M g is below `tol`. A step is accepted when the
# log-likelihood does not fall or, once the decrement is below `tol_ll` (the
# predicted gain, about half of it, is then within the rounding of a sum
# over thousands of situations), when the decrement falls. Returns theta,
# ll, grad, iterations, and the certificate: `decrement`, the reference's
# own BHHH decrement g' (S'S)^-1 g at theta, with its BHHH matrix `B`.
quad_mxl_fit <- function(theta, d, spec, chosen, rule, metric = NULL, tol = 1e-12,
                         tol_ll = 1e-8, max_iter = 100L) {
  bhhh <- is.null(metric)
  eval <- function(th) quad_mxl(th, d, spec, chosen, rule, what = if (bhhh) "scores" else "grad")
  direction <- function(s) {
    if (bhhh) solve(crossprod(s$scores), s$grad) else drop(metric %*% s$grad)
  }
  cur <- eval(theta)
  step <- direction(cur)
  dec <- sum(cur$grad * step)
  it <- 0L
  while (dec >= tol && it < max_iter) {
    if (!is.finite(dec)) stop("quad_mxl_fit(): the step is not finite")
    t <- 1
    repeat {
      nxt <- eval(theta + t * step)
      step_n <- direction(nxt)
      dec_n <- sum(nxt$grad * step_n)
      if (is.finite(nxt$ll) && if (dec < tol_ll) dec_n < dec else nxt$ll >= cur$ll) break
      t <- t / 2
      if (t < 1e-10) stop("quad_mxl_fit(): no ascent along the step")
    }
    theta <- theta + t * step
    cur <- nxt
    step <- step_n
    dec <- dec_n
    it <- it + 1L
  }
  sc <- if (bhhh) cur else quad_mxl(theta, d, spec, chosen, rule, what = "scores")
  B <- crossprod(sc$scores)
  list(theta = theta, ll = cur$ll, grad = cur$grad, iterations = it,
       decrement = sum(sc$grad * solve(B, sc$grad)), B = B)
}

# Q1.3: D4's MXL variant (log-normal on -distance, normal on quality,
# K_w = 2) on `design` (D4-40, F3's, by default) at N students, the
# never-chosen alternatives dropped (drop_unbounded_alts()); choicer at S
# draws (generate, polish = "newton") against the quadrature ML. The rule is
# the first of `grids` (nodes on -distance's and quality's eta) whose
# log-likelihood at choicer's estimate is within `ll_tol` of the last's
# (the finest); the reference is then certified on the finest rule too: its
# BHHH decrement there bounds the remaining error, |dtheta|/SE of about
# sqrt(decrement). With s_q = 1.75 the quality integrand is a steep logistic
# in eta, so Gauss-Hermite needs about a hundred nodes on it (measured
# 2026-11-02 on D4-small, N = 2000: ll within 3e-6 at (96, 96), 1e-11
# between 96 and 128 nodes on distance). Returns check rows; `verbose`
# prints progress.
ref_q13 <- function(N = 2000L, S = 4000L, design = GEO_D4_40,
                    grids = list(c(32L, 64L), c(48L, 96L), c(64L, 128L), c(96L, 192L)),
                    ll_tol = 1e-6, threshold = 0.05,
                    seed = vs_seed("Q1.3", paste0("N", N), 1), verbose = FALSE) {
  t_start <- proc.time()[["elapsed"]]
  say <- function(...) if (verbose) cat(sprintf(...), "\n", sep = "")
  sim <- d4_simulate(N, seed, "default", design = design)
  dr <- drop_unbounded_alts(sim$data)
  dt <- dr$data
  alts <- sort(unique(dt$alt))
  spec <- d4_spec("default", sim$meta$schools, alts)

  # choicer's MSL fit
  t0 <- proc.time()[["elapsed"]]
  fit <- suppressMessages(choicer::run_mxlogit(
    dt, "id", "alt", "choice", D4_X, D4_W, rc_dist = c(1L, 0L), rc_mean = FALSE,
    include_outside_option = TRUE, S = S, draws = "generate", seed = seed + 1L,
    polish = "newton"))
  t_fit <- proc.time()[["elapsed"]] - t0
  th_msl <- coef(fit)
  if (!identical(names(th_msl), theta_names("mxl", spec))) {
    stop("ref_q13(): choicer's parameters differ from the suite's names")
  }
  se <- fit$se
  say("choicer: N = %d, %d rows, %d parameters, S = %d: %.1f s, polish converged %s",
      max(dt$id) - dr$counts[["situations"]], nrow(dt), length(th_msl), S, t_fit,
      isTRUE(fit$optimizer$polish$converged))

  # The rule: the first grid whose log-likelihood at choicer's estimate is
  # within ll_tol of the finest grid's.
  d <- design(dt, spec)
  chosen <- chosen_rows(dt$choice, d$sit)
  t0 <- proc.time()[["elapsed"]]
  ll_g <- vapply(grids, function(Q) quad_mxl(th_msl, d, spec, chosen, gh_product2(Q))$ll, 1)
  G <- length(grids)
  gap <- abs(ll_g - ll_g[G])
  k <- which(gap[-G] < ll_tol)[1L]
  if (is.na(k)) k <- G
  lab <- function(Q) paste(Q, collapse = "x")
  say("quadrature ll at the MSL estimate by grid: %s (%.0f s)",
      paste(sprintf("%s: %.10f", vapply(grids, lab, ""), ll_g), collapse = "; "),
      proc.time()[["elapsed"]] - t0)

  # The path in choicer's metric (its inverse Hessian; BHHH alone took 23
  # steps at the plan's design), the certificate in the reference's own BHHH
  # metric, on the chosen grid and on the finest.
  t0 <- proc.time()[["elapsed"]]
  V <- fit$vcov
  ref <- quad_mxl_fit(th_msl, d, spec, chosen, gh_product2(grids[[k]]),
                      metric = if (!is.null(V) && all(is.finite(V))) unname(V), tol = 1e-14)
  fine <- quad_mxl(ref$theta, d, spec, chosen, gh_product2(grids[[G]]), what = "grad")
  dec_fine <- sum(fine$grad * solve(ref$B, fine$grad))
  t_ref <- proc.time()[["elapsed"]] - t0
  z <- abs(ref$theta - th_msl) / se
  worst <- names(which.max(z))
  say("reference: grid %s, %d steps, BHHH decrement %.2e (finest grid %.2e), %.1f s; max |dtheta|/SE %.4f (%s)",
      lab(grids[[k]]), ref$iterations, ref$decrement, dec_fine, t_ref, max(z), worst)

  # Comparisons by .ref_row(): a NaN fails. q_row(): a verdict, or information.
  rbind(
    q_row("Q1.3", sprintf("choicer polish converged (decrement over the mean weight; S = %d, %.0f s)",
                          S, t_fit),
          fit$optimizer$polish$decrement %||% NA_real_,
          pass = isTRUE(fit$optimizer$polish$converged)),
    .ref_row("Q1.3", sprintf("quadrature: |ll(%s) - ll(%s)| at the MSL estimate", lab(grids[[k]]),
                             lab(grids[[G]])), gap[k], ll_tol),
    .ref_row("Q1.3", sprintf("reference converged: its BHHH decrement after %d steps (%.0f s)",
                             ref$iterations, t_ref), ref$decrement, 1e-12),
    .ref_row("Q1.3", sprintf("reference on the finest rule (%s): BHHH decrement",
                             lab(grids[[G]])), dec_fine, 1e-6),
    q_row("Q1.3", sprintf("ll(ML) - ll(MSL estimate), quadrature, %d parameters", length(th_msl)),
          ref$ll - ll_g[k]),
    .ref_row("Q1.3", sprintf("max |theta_ML - theta_MSL| / SE (worst %s; N = %d, S = %d)", worst,
                             length(unique(dt$id)), S), max(z), threshold),
    q_row("Q1.3", sprintf("median |theta_ML - theta_MSL| / SE; %.0f s in all",
                          proc.time()[["elapsed"]] - t_start), stats::median(z)))
}
