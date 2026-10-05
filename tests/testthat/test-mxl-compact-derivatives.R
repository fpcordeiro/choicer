# The mixed logit derivative kernels on designs with many alternatives and
# sparse choice sets. A decision maker's score and Hessian are nonzero only on
# the continuous parameters (beta, mu, L) and the free ASCs of the
# alternatives they saw, so the kernels work in that compact block and add it
# into the n x n result; these designs (30 alternatives, choice sets of 5-10)
# make the blocks small next to n and exercise the maps: the normalized
# reference alternative, the outside option, panels whose situations offer
# different alternatives, a situation listing an alternative twice, the
# reference never offered, extra trailing ASCs, draw batches and threads.
# Kernel calls go through mxlp_call() (helper-mxl-panel.R).

# A prepared many-alternative design in the stacked layout, with the fields
# mxlp_call() reads. U decision makers with 1-3 situations each (panel) or one
# (cross-section); each situation offers `m_range` distinct alternatives drawn
# from 1..J_codes (from 2..J_codes when `absent_ref`), plus one repeated row
# when `dup`; weights vary across decision makers.
cdx_fixture <- function(name, seed, panel = TRUE,
                        include_outside_option = FALSE, J = 30L,
                        m_range = c(5L, 10L), U = if (panel) 10L else 24L,
                        rc_dist = c(0L, 1L), rc_correlation = TRUE,
                        rc_mean = TRUE, use_asc = TRUE, K_x = 2L, S = 12L,
                        dup = FALSE, absent_ref = FALSE, J_codes = J,
                        ref_only = FALSE, shuffle = FALSE) {
  set.seed(seed)
  rc_dist <- as.integer(rc_dist)
  K_w <- length(rc_dist)
  Ti <- if (panel) sample(1:3, U, replace = TRUE) else rep(1L, U)
  N <- sum(Ti)
  unit <- rep(seq_len(U), Ti)
  pool <- if (absent_ref) 2:J_codes else seq_len(J_codes)
  alts <- lapply(seq_len(N), function(t) {
    # `ref_only`: every situation of unit 1, and about a third of the others,
    # offer only the reference alternative (no free ASC)
    if (ref_only && (unit[t] == 1L || stats::runif(1L) < 1 / 3)) return(1L)
    m <- sample(m_range[1L]:m_range[2L], 1L)
    a <- sort(pool[sample.int(length(pool), min(m, length(pool)))])
    if (dup) a <- sort(c(a, a[sample.int(length(a), 1L)]))
    # `shuffle`: rows in no particular order (repeats not adjacent)
    if (shuffle) a <- a[sample.int(length(a))]
    a
  })
  M <- lengths(alts)
  n_rows <- sum(M)
  lo <- if (include_outside_option) 0L else 1L
  choice_idx <- vapply(M, function(m) sample(lo:m, 1L), integer(1L))
  w_u <- stats::runif(U, 0.5, 2)
  L_size <- if (rc_correlation) (K_w * (K_w + 1L)) %/% 2L else K_w
  n_delta <- if (!use_asc) 0L else if (include_outside_option) J else J - 1L
  theta <- c(stats::rnorm(K_x, sd = 0.4), if (rc_mean) stats::rnorm(K_w, sd = 0.3),
             stats::rnorm(L_size, sd = 0.3) - 0.5, stats::rnorm(n_delta, sd = 0.4))
  list(
    name = name, X = matrix(stats::rnorm(n_rows * K_x), n_rows, K_x),
    W = matrix(stats::rnorm(n_rows * K_w), n_rows, K_w),
    alt_idx = as.integer(unlist(alts)), choice_idx = as.integer(choice_idx),
    M = as.integer(M), weights = rep(w_u, Ti),
    Ti = if (panel) as.integer(Ti) else NULL,
    eta = get_halton_normals(S, U, K_w), theta = theta, rc_dist = rc_dist,
    rc_correlation = rc_correlation, rc_mean = rc_mean, use_asc = use_asc,
    include_outside_option = include_outside_option, K_w = K_w, K_x = K_x,
    J = J, S = S, U = U, N = N
  )
}

# Index of the first ASC in theta, and the global parameter indices of the
# compact block of decision maker u: the continuous parameters and the free
# ASCs of the alternatives in u's rows.
cdx_delta_start <- function(fx) {
  L_size <- if (fx$rc_correlation) (fx$K_w * (fx$K_w + 1L)) %/% 2L else fx$K_w
  fx$K_x + (if (fx$rc_mean) fx$K_w else 0L) + L_size + 1L
}
cdx_block <- function(fx, u) {
  Ti <- if (is.null(fx$Ti)) rep(1L, fx$N) else fx$Ti
  sit <- which(rep(seq_along(Ti), Ti) == u)
  row_off <- c(0L, cumsum(fx$M))
  codes <- fx$alt_idx[unlist(lapply(sit, function(t) (row_off[t] + 1L):row_off[t + 1L]))]
  free <- if (fx$include_outside_option) codes else codes[codes > 1L] - 1L
  d0 <- cdx_delta_start(fx)
  c(seq_len(d0 - 1L), sort(unique(free)) + d0 - 1L)
}

cdx_cells <- function() {
  list(
    cdx_fixture("xsec", 1001, panel = FALSE),
    cdx_fixture("xsec, outside option", 1002, panel = FALSE,
                include_outside_option = TRUE),
    cdx_fixture("panel", 1003),
    cdx_fixture("panel, outside option", 1004, include_outside_option = TRUE),
    cdx_fixture("panel, repeated alternative", 1005, dup = TRUE),
    cdx_fixture("panel, reference never offered", 1006, absent_ref = TRUE),
    cdx_fixture("panel, trailing ASCs", 1007, J_codes = 24L),
    cdx_fixture("panel, normal, uncorrelated, no means", 1008,
                rc_dist = c(0L, 0L), rc_correlation = FALSE, rc_mean = FALSE),
    cdx_fixture("panel, no ASCs", 1009, use_asc = FALSE),
    cdx_fixture("one decision maker", 1010, U = 1L),
    cdx_fixture("panel, situations with the reference only", 1020,
                ref_only = TRUE),
    cdx_fixture("xsec, shuffled rows, repeated alternatives", 1021,
                panel = FALSE, dup = TRUE, shuffle = TRUE)
  )
}

test_that("the compact Hessian matches numDeriv of the gradient on sparse designs", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  for (fx in cdx_cells()) {
    num <- numDeriv::jacobian(
      function(th) mxlp_call("gradient", fx, theta = th)$gradient, fx$theta)
    for (gen in c(FALSE, TRUE)) {
      H <- mxlp_call("hessian", fx, generate = gen)
      ref <- if (gen) numDeriv::jacobian(
        function(th) mxlp_call("gradient", fx, theta = th, generate = TRUE)$gradient,
        fx$theta) else num
      mxlp_expect_close(H, ref, 1e-6, sprintf("[%s%s] Hessian vs numDeriv",
                                             fx$name, if (gen) ", generate" else ""))
    }
    # Pass 2 forming each situation's rows of W Gamma in turn
    if (!is.null(fx$Ti)) {
      mxlp_expect_close(mxlp_call("hessian", fx, draw_batch = 5L), num, 1e-6,
                        sprintf("[%s] Hessian with draw batches vs numDeriv",
                                fx$name))
    }
  }
})

test_that("relabeling the alternatives permutes the ASC block of the Hessian", {
  for (ioo in c(FALSE, TRUE)) {
    fx <- cdx_fixture("relabel", 1011, include_outside_option = ioo)
    # Without an outside option alternative 1 is the normalized reference, so
    # it keeps its code; the others are permuted.
    J <- fx$J
    perm <- if (ioo) sample.int(J) else c(1L, 1L + sample.int(J - 1L))
    fy <- fx
    fy$alt_idx <- perm[fx$alt_idx]
    d0 <- cdx_delta_start(fx)
    n_delta <- length(fx$theta) - d0 + 1L
    # new free ASC of old free ASC k
    new_of_old <- if (ioo) perm else perm[-1L] - 1L
    fy$theta[d0 - 1L + new_of_old] <- fx$theta[d0 - 1L + seq_len(n_delta)]
    p <- c(seq_len(d0 - 1L), d0 - 1L + new_of_old)
    H <- mxlp_call("hessian", fx)
    G <- mxlp_call("hessian", fy)
    mxlp_expect_close(G[p, p], H, 1e-12,
                      sprintf("[outside option = %s] relabeled Hessian", ioo))
  }
})

test_that("the Hessian is exactly symmetric, linear in the weights and stable across threads", {
  on.exit(set_num_threads(2L), add = TRUE)
  cells <- cdx_cells()
  for (fx in cells[vapply(cells, function(f) f$name %in% c(
    "xsec", "panel, outside option", "panel, repeated alternative",
    "panel, situations with the reference only",
    "xsec, shuffled rows, repeated alternatives"), TRUE)]) {
    for (cfg in list(list(), list(generate = TRUE), list(draw_batch = 3L))) {
      what <- sprintf("[%s%s]", fx$name, if (length(cfg)) paste0(", ", names(cfg)) else "")
      set_num_threads(2L)
      H2 <- do.call(mxlp_call, c(list("hessian", fx), cfg))
      expect_identical(H2, t(H2), label = paste(what, "symmetric"))
      mxlp_expect_close(do.call(mxlp_call, c(list("hessian", fx,
                                                 weights = 2 * fx$weights), cfg)),
                        2 * H2, 1e-12, paste(what, "H(2w) vs 2 H(w)"))
      set_num_threads(1L)
      H1 <- do.call(mxlp_call, c(list("hessian", fx), cfg))
      expect_identical(H1, t(H1), label = paste(what, "symmetric, 1 thread"))
      mxlp_expect_close(H1, H2, 1e-12, paste(what, "1 vs 2 threads"))
    }
  }
})

test_that("non-finite blocks spread over the result as the dense products did", {
  # Draws 1-5 of decision maker u give its choice probability exactly zero
  # (its one chosen alternative's random utility overflows to -Inf), so its
  # likelihood is finite but its Hessian pieces at those zero-weight draws
  # are NaN. The dense kernel spread a NaN row of the score stashes over the
  # whole row and column of the result (0 * NaN); the result keeps that:
  # every NaN lies in a row and column that are NaN throughout.
  fx <- cdx_fixture("zero-weight draws", 1012, rc_dist = 0L)
  u <- 3L
  Ti <- fx$Ti
  first <- cumsum(c(1L, Ti[-length(Ti)]))
  row_off <- c(0L, cumsum(fx$M))
  rows_u <- (row_off[first[u]] + 1L):row_off[first[u] + Ti[u]]
  fx$W[rows_u, 1L] <- 0
  fx$W[row_off[first[u]] + fx$choice_idx[first[u]], 1L] <- 10
  fx$eta[1L, 1:5, u] <- -1e308
  H <- mxlp_call("hessian", fx)
  bad <- which(apply(is.nan(H), 1L, all))
  # Only the Cholesky parameter's row was NaN in the dense kernel: the
  # overflowing draw enters the unit's scores through it alone.
  expect_identical(unname(bad), cdx_delta_start(fx) - 1L)
  expect_true(all(is.nan(H[bad, ])) && all(is.nan(H[, bad])))
  expect_true(all(is.finite(H[-bad, -bad])))

  # A non-finite weight multiplied the zeros outside the unit's block: those
  # entries are NaN. (R rejects such weights; direct calls can pass them.)
  fx <- cdx_fixture("infinite weight", 1013)
  u <- 2L
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  fx$weights[first[u] - 1L + seq_len(fx$Ti[u])] <- Inf
  H <- mxlp_call("hessian", fx)
  blk <- cdx_block(fx, u)
  out <- matrix(TRUE, nrow(H), ncol(H))
  out[blk, blk] <- FALSE
  expect_true(all(is.nan(H[out])))
  expect_false(any(is.finite(H[blk, blk])))

  # A missing weight: w * 0 is missing everywhere (NA or NaN, by platform).
  fx$weights[first[u] - 1L + seq_len(fx$Ti[u])] <- NA_real_
  expect_true(all(is.na(mxlp_call("hessian", fx))))
  fx$weights[first[u] - 1L + seq_len(fx$Ti[u])] <- NaN
  expect_true(all(is.na(mxlp_call("hessian", fx))))

  # A decision maker whose utilities overflow is left out, so its weight,
  # finite or not, spreads nothing (bit for bit at one thread).
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  fx <- cdx_fixture("skipped unit", 1022)
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  row_off <- c(0L, cumsum(fx$M))
  u <- 3L
  fx$theta[1L] <- 2
  fx$X[row_off[first[u]] + fx$choice_idx[first[u]], 1L] <- -1e308
  H1 <- mxlp_call("hessian", fx)
  expect_true(all(is.finite(H1)))
  fx$weights[first[u] - 1L + seq_len(fx$Ti[u])] <- Inf
  expect_identical(mxlp_call("hessian", fx), H1,
                   label = "skipped unit with an infinite weight")
})
