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

cdx_unit_weights <- function(fx) {
  fx$weights[if (is.null(fx$Ti)) seq_len(fx$N) else cumsum(c(1L, fx$Ti[-length(fx$Ti)]))]
}

# --- BHHH --------------------------------------------------------------------

# Main's dense BHHH update, sum_u w_u s_u s_u', formed entry by entry as the
# kernel forms it: the product of the two scores, times the weight, added in
# unit order (0 * NaN and 0 * Inf included).
cdx_bhhh_ref <- function(S_u, w_u) {
  Reduce(`+`, lapply(seq_along(w_u), function(v) w_u[v] * outer(S_u[v, ], S_u[v, ])),
         matrix(0, ncol(S_u), ncol(S_u)))
}

test_that("the compact BHHH is the weighted cross-product of the unit scores", {
  on.exit(set_num_threads(2L), add = TRUE)
  for (fx in cdx_cells()) {
    w_u <- cdx_unit_weights(fx)
    for (nt in 1:2) {
      set_num_threads(nt)
      S_u <- mxlp_call("scores", fx)
      B <- mxlp_call("bhhh", fx)
      expect_identical(B, t(B), label = sprintf("[%s] BHHH symmetric", fx$name))
      mxlp_expect_close(B, crossprod(sqrt(w_u) * S_u), 1e-12,
                        sprintf("[%s, %d threads] BHHH vs crossprod(sqrt(w) S)",
                                fx$name, nt))
    }
    mxlp_expect_close(mxlp_call("bhhh", fx, generate = TRUE),
                      crossprod(sqrt(w_u) * mxlp_call("scores", fx, generate = TRUE)),
                      1e-12, sprintf("[%s] BHHH vs scores (generate)", fx$name))
    mxlp_expect_close(mxlp_call("bhhh", fx, draw_batch = 3L),
                      crossprod(sqrt(w_u) * mxlp_call("scores", fx, draw_batch = 3L)),
                      1e-12, sprintf("[%s] BHHH vs scores (draw batches)", fx$name))
  }
})

test_that("at one thread the BHHH is the dense update bit for bit", {
  # Each entry is the product of the two scores, times the weight, added in
  # unit order, rounded at each step as R rounds the reference on every
  # toolchain (main's update rounded the same way under clang and GCC -O2).
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  for (fx in cdx_cells()) {
    expect_identical(mxlp_call("bhhh", fx),
                     cdx_bhhh_ref(mxlp_call("scores", fx), cdx_unit_weights(fx)),
                     label = sprintf("[%s] BHHH vs the dense update", fx$name))
  }
})

test_that("non-finite scores and weights spread over the BHHH as before", {
  # Patterns of NaN and infinities are those of the dense update.
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  same_pattern <- function(B, ref, what) {
    expect_identical(is.nan(B), is.nan(ref), label = paste(what, "NaN"))
    expect_identical(is.infinite(B), is.infinite(ref), label = paste(what, "Inf"))
    expect_identical(B[is.infinite(B)], ref[is.infinite(ref)],
                     label = paste(what, "infinities and their signs"))
    expect_identical(B[is.finite(B)], ref[is.finite(ref)],
                     label = paste(what, "finite entries"))
  }
  # A decision maker whose utilities overflow has a NaN score.
  fx <- cdx_fixture("overflow", 1014)
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  row_off <- c(0L, cumsum(fx$M))
  u <- 2L
  fx$theta[1L] <- 2
  fx$X[row_off[first[u]] + fx$choice_idx[first[u]], 1L] <- -1e308
  S_u <- mxlp_call("scores", fx)
  expect_true(any(is.nan(S_u[u, ])))
  same_pattern(mxlp_call("bhhh", fx), cdx_bhhh_ref(S_u, cdx_unit_weights(fx)),
               "overflowed unit")

  # A score entry past the largest double with finite utilities: a covariate
  # of 1e308 on a non-chosen alternative with a coefficient of 1e-306 in
  # every situation of a two-situation decision maker, so beta_1's score
  # overflows to -Inf while the decision maker's other entries stay finite.
  fx <- cdx_fixture("infinite score", 1023)
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  row_off <- c(0L, cumsum(fx$M))
  u <- which(fx$Ti >= 2L)[1L]
  fx$theta[1L] <- 1e-306
  for (t in first[u] - 1L + seq_len(fx$Ti[u])) {
    other <- setdiff(seq_len(fx$M[t]), fx$choice_idx[t])[1L]
    fx$X[row_off[t] + other, 1L] <- 1e308
  }
  S_u <- mxlp_call("scores", fx)
  expect_true(is.infinite(S_u[u, 1L]) && all(is.finite(S_u[u, -1L])))
  same_pattern(mxlp_call("bhhh", fx), cdx_bhhh_ref(S_u, cdx_unit_weights(fx)),
               "infinite score")

  # An infinite weight: w * 0 is NaN outside the unit's block, w * s
  # infinite inside.
  fx <- cdx_fixture("infinite weight", 1015)
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  u <- 4L
  fx$weights[first[u] - 1L + seq_len(fx$Ti[u])] <- Inf
  B <- mxlp_call("bhhh", fx)
  same_pattern(B, cdx_bhhh_ref(mxlp_call("scores", fx), cdx_unit_weights(fx)),
               "infinite weight")
  blk <- cdx_block(fx, u)
  out <- matrix(TRUE, nrow(B), ncol(B))
  out[blk, blk] <- FALSE
  expect_true(all(is.nan(B[out])))
  expect_false(any(is.finite(B[blk, blk])))

  # Missing weights: every entry is missing (NA or NaN, by platform).
  for (w_bad in c(NA_real_, NaN)) {
    fx$weights[first[u] - 1L + seq_len(fx$Ti[u])] <- w_bad
    expect_true(all(is.na(mxlp_call("bhhh", fx))))
  }
})

# --- Cluster meat -------------------------------------------------------------

# The cluster meat kernel on a fixture, with labels per likelihood unit coded
# as the R route codes them.
cdx_meat <- function(fx, cl, generate = FALSE, weights = fx$weights) {
  lab <- as.character(cl)
  args <- list(theta = fx$theta, X = fx$X, W = fx$W, alt_idx = fx$alt_idx,
               choice_idx = fx$choice_idx, M = fx$M, weights = weights,
               cluster = match(lab, unique(lab)),
               eta_draws = if (generate) array(0, dim = c(fx$K_w, 0L, 0L)) else fx$eta,
               rc_dist = fx$rc_dist, rc_correlation = fx$rc_correlation,
               rc_mean = fx$rc_mean, use_asc = fx$use_asc,
               include_outside_option = fx$include_outside_option)
  if (generate) args <- c(args, list(gen_seed = 0L, gen_scramble = 0L, gen_S = fx$S))
  if (!is.null(fx$Ti)) args$Ti <- fx$Ti
  do.call(mxl_cluster_meat_parallel, args)
}
cdx_units <- function(fx) if (is.null(fx$Ti)) fx$N else length(fx$Ti)

test_that("the cluster meat equals crossprod(rowsum(w * S, cluster))", {
  on.exit(set_num_threads(2L), add = TRUE)
  for (fx in cdx_cells()) {
    U <- cdx_units(fx)
    set.seed(1100)
    labels <- list(singletons = seq_len(U), one = rep("a", U),
                   coarse = sample(1:4, U, replace = TRUE),
                   two = rep(c(2.5, 1.5), length.out = U),
                   shuffled = sample(seq_len(U)) * 7L)
    S_u <- mxlp_call("scores", fx)
    w_u <- cdx_unit_weights(fx)
    for (nt in 1:2) {
      set_num_threads(nt)
      for (lb in names(labels)) {
        cl <- labels[[lb]]
        mxlp_expect_close(cdx_meat(fx, cl),
                          crossprod(rowsum(w_u * S_u, group = as.character(cl))),
                          1e-12, sprintf("[%s, %s, %d threads] cluster meat",
                                         fx$name, lb, nt))
      }
      # Singletons: the robust meat, the BHHH kernel with squared weights
      mxlp_expect_close(cdx_meat(fx, seq_len(U)),
                        mxlp_call("bhhh", fx, weights = fx$weights^2), 1e-12,
                        sprintf("[%s, %d threads] singletons vs BHHH(w^2)",
                                fx$name, nt))
    }
    mxlp_expect_close(cdx_meat(fx, labels$coarse, generate = TRUE),
                      crossprod(rowsum(w_u * mxlp_call("scores", fx, generate = TRUE),
                                       group = as.character(labels$coarse))),
                      1e-12, sprintf("[%s] cluster meat (generate)", fx$name))
  }
})

test_that("the cluster meat does not depend on how labels are coded or ordered", {
  fx <- cdx_fixture("codes", 1016, U = 40L)
  set.seed(1101)
  cl <- sample(1:6, 40L, replace = TRUE)
  M1 <- cdx_meat(fx, cl)
  mxlp_expect_close(cdx_meat(fx, letters[7L - cl]), M1, 1e-12, "relabeled")
  mxlp_expect_close(cdx_meat(fx, factor(cl, levels = 6:1)), M1, 1e-12, "factor")
  expect_identical(M1, t(M1))
  # A cluster cut by every chunk boundary, at many threads' worth of chunks
  mxlp_expect_close(cdx_meat(fx, rep(1L, 40L)),
                    tcrossprod(colSums(cdx_unit_weights(fx) * mxlp_call("scores", fx))),
                    1e-12, "one cluster")
})

test_that("a non-finite score spreads over the cluster meat as crossprod() did", {
  fx <- cdx_fixture("overflow", 1017, U = 12L)
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  row_off <- c(0L, cumsum(fx$M))
  u <- 5L
  fx$theta[1L] <- 2
  fx$X[row_off[first[u]] + fx$choice_idx[first[u]], 1L] <- -1e308
  cl <- rep(1:4, length.out = 12L)
  S_u <- mxlp_call("scores", fx)
  ref <- crossprod(rowsum(cdx_unit_weights(fx) * S_u, group = as.character(cl)))
  M <- cdx_meat(fx, cl)
  expect_identical(is.nan(M), is.nan(ref))
  expect_identical(is.finite(M), is.finite(ref))
})

test_that("non-finite weights and scores, and overflowing sums, spread as crossprod() did", {
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  same_pattern <- function(M, ref, what) {
    expect_identical(is.nan(M), is.nan(ref), label = paste(what, "NaN"))
    expect_identical(is.infinite(M), is.infinite(ref), label = paste(what, "Inf"))
    expect_identical(M[is.infinite(M)], ref[is.infinite(ref)],
                     label = paste(what, "infinities and their signs"))
    expect_equal(M[is.finite(M)], ref[is.finite(ref)], tolerance = 1e-12,
                 label = paste(what, "finite entries"))
  }
  ref_meat <- function(fx, cl, w = fx$weights) {
    first <- if (is.null(fx$Ti)) seq_len(fx$N) else cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
    crossprod(rowsum(w[first] * mxlp_call("scores", fx), group = as.character(cl)))
  }
  fx <- cdx_fixture("weights", 1024, U = 12L)
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  cl <- rep(1:4, length.out = 12L)
  for (w_bad in c(Inf, -Inf)) {
    w <- fx$weights
    w[first[5L] - 1L + seq_len(fx$Ti[5L])] <- w_bad
    same_pattern(cdx_meat(fx, cl, weights = w), ref_meat(fx, cl, w),
                 sprintf("weight %s", w_bad))
  }
  for (w_bad in c(NA_real_, NaN)) {
    w <- fx$weights
    w[first[5L] - 1L + seq_len(fx$Ti[5L])] <- w_bad
    expect_true(all(is.na(cdx_meat(fx, cl, weights = w))))
  }
  # An infinite score entry with finite utilities (as in the BHHH test)
  fx <- cdx_fixture("infinite score", 1023)
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  row_off <- c(0L, cumsum(fx$M))
  u <- which(fx$Ti >= 2L)[1L]
  fx$theta[1L] <- 1e-306
  for (t in first[u] - 1L + seq_len(fx$Ti[u])) {
    other <- setdiff(seq_len(fx$M[t]), fx$choice_idx[t])[1L]
    fx$X[row_off[t] + other, 1L] <- 1e308
  }
  cl <- rep(1:3, length.out = length(fx$Ti))
  same_pattern(cdx_meat(fx, cl), ref_meat(fx, cl), "infinite score")
  # Finite scores whose cluster sum overflows: two cross-sectional units of
  # one cluster, each with a beta_1 score near -1e308
  fx <- cdx_fixture("overflowing sum", 1025, panel = FALSE, U = 10L)
  row_off <- c(0L, cumsum(fx$M))
  fx$theta[1L] <- 1e-306
  for (t in 1:2) {
    other <- setdiff(seq_len(fx$M[t]), fx$choice_idx[t])[1L]
    fx$X[row_off[t] + other, 1L] <- 1e308
  }
  S_u <- mxlp_call("scores", fx)
  expect_true(all(is.finite(S_u)))
  cl <- c(1L, 1L, 2:9)
  ref <- ref_meat(fx, cl)
  expect_true(any(is.nan(ref)))
  same_pattern(cdx_meat(fx, cl), ref, "overflowing cluster sum")
})

test_that("unused codes and skewed units do not change the cluster meat", {
  on.exit(set_num_threads(2L), add = TRUE)
  fx <- cdx_fixture("codes", 1026, U = 30L)
  set.seed(1103)
  cl <- sample(1:5, 30L, replace = TRUE)
  args <- function(codes) list(fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx,
                               fx$M, fx$weights, codes, fx$eta, fx$rc_dist,
                               fx$rc_correlation, fx$rc_mean, Ti = fx$Ti)
  set_num_threads(1L)
  # codes in the same order, with gaps: the same sort, the same result
  expect_identical(do.call(mxl_cluster_meat_parallel, args(3L * cl + 2L)),
                   do.call(mxl_cluster_meat_parallel, args(cl)))
  set_num_threads(2L)
  # one decision maker holding most rows. Sorted first, it has the first
  # chunk to itself and the next few chunks hold one unit each; sorted last,
  # it leaves the trailing chunks empty
  fy <- cdx_fixture("skewed", 1027, U = 12L, m_range = c(25L, 30L))
  fy$eta <- get_halton_normals(fy$S, 12L, fy$K_w)
  n_sit <- sum(fy$Ti)
  for (Ti in list(c(n_sit - 11L, rep(1L, 11L)), c(rep(1L, 11L), n_sit - 11L))) {
    fy$Ti <- Ti
    fy$weights <- rep(stats::runif(12L, 0.5, 2), Ti)
    w_u <- fy$weights[cumsum(c(1L, Ti[-12L]))]
    for (cl in list(seq_len(12L), rep(1:2, 6L), rep(1L, 12L))) {
      mxlp_expect_close(cdx_meat(fy, cl),
                        crossprod(rowsum(w_u * mxlp_call("scores", fy),
                                         group = as.character(cl))),
                        1e-12, "skewed units")
    }
  }
})

test_that("the cluster meat checks its labels after the other inputs", {
  fx <- cdx_fixture("labels", 1018, U = 6L)
  expect_error(cdx_meat(fx, 1:5), "cluster length (5) does not match the number of likelihood units (6)",
               fixed = TRUE)
  bad <- list(theta = fx$theta, X = fx$X, W = fx$W, alt_idx = fx$alt_idx,
              choice_idx = fx$choice_idx, M = fx$M, weights = fx$weights,
              cluster = c(1L, NA, 2L, 2L, 3L, 3L), eta_draws = fx$eta,
              rc_dist = fx$rc_dist, rc_correlation = fx$rc_correlation,
              rc_mean = fx$rc_mean, Ti = fx$Ti)
  expect_error(do.call(mxl_cluster_meat_parallel, bad),
               "cluster codes must be positive integers (unit 2 has NA).", fixed = TRUE)
  bad$cluster <- c(1L, 0L, 2L, 2L, 3L, 3L)
  expect_error(do.call(mxl_cluster_meat_parallel, bad),
               "cluster codes must be positive integers (unit 2 has 0).", fixed = TRUE)
  bad$weights <- bad$weights[-1L]
  expect_error(do.call(mxl_cluster_meat_parallel, bad),
               "weights length", fixed = TRUE)
})

# --- Shared accumulator -------------------------------------------------------

# With more than one thread, the kernels add into per-thread packed triangles
# while T n (n + 1) / 2 doubles fit in acc_bytes (2 GiB by default) and share
# the result above it: the continuous rows through per-thread buffers, the
# ASC block with atomic additions. acc_bytes = 0 forces the shared result.

test_that("the kernels share the result exactly when the triangles exceed acc_bytes", {
  acc <- function(n, K_c, T, bytes) test_mxl_acc(n, K_c, T, bytes)
  tri <- function(n, T) 8 * T * n * (n + 1) / 2
  # the T triangles fit at exactly acc_bytes
  expect_false(acc(100, 5, 4, tri(100, 4))$shared)
  expect_true(acc(100, 5, 4, tri(100, 4) - 1)$shared)
  # the default 2 GiB at eleven threads: per-thread up to 6,985 parameters
  expect_false(acc(6985, 19, 11, 2^31)$shared)
  expect_true(acc(6986, 19, 11, 2^31)$shared)
  # NaN, negative and zero budgets share; Inf never does; one thread never
  for (bytes in c(NaN, -1, 0)) {
    expect_true(acc(100, 5, 2, bytes)$shared)
    expect_false(acc(100, 5, 1, bytes)$shared)
  }
  expect_false(acc(100, 5, 64, Inf)$shared)
  # A thread's buffer packs each column's first min(s, j + 1) rows: the
  # triangle (s = n) or the continuous rows (s = K_c), never longer than the
  # triangle, and the same when every parameter is continuous.
  expect_equal(diff(acc(100, 5, 2, Inf)$col_off), 1:100)
  for (K_c in c(1, 7, 99, 100)) {
    off <- acc(100, K_c, 2, 0)$col_off
    expect_equal(diff(off), pmin(K_c, 1:100))
    expect_lte(off[101L], 100 * 101 / 2)
  }
  expect_identical(acc(100, 100, 2, 0)$col_off, acc(100, 100, 2, Inf)$col_off)
  # offsets past 2^32 doubles, and the continuous rows' length
  expect_equal(acc(1e5, 19, 2, Inf)$col_off[1e5 + 1], 1e5 * (1e5 + 1) / 2)
  expect_equal(acc(1e5, 19, 2, 0)$col_off[1e5 + 1], 1e5 * 19 - 19 * 18 / 2)
})

test_that("a shared result equals the per-thread accumulators' sum", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  on.exit(set_num_threads(2L), add = TRUE)
  cells <- cdx_cells()
  for (fx in cells[vapply(cells, function(f) f$name %in% c(
    "xsec, outside option", "panel", "panel, repeated alternative",
    "panel, trailing ASCs", "panel, no ASCs", "one decision maker",
    "xsec, shuffled rows, repeated alternatives"), TRUE)]) {
    cl <- (seq_len(cdx_units(fx)) - 1L) %/% 3L + 1L
    for (k in c("hessian", "bhhh", "meat")) {
      for (gen in c(FALSE, TRUE)) {
        what <- sprintf("[%s%s] %s", fx$name, if (gen) ", generate" else "", k)
        kcall <- function(...) {
          mxlp_call(k, fx, generate = gen, cluster = if (k == "meat") cl, ...)
        }
        set_num_threads(2L)
        shared <- kcall(acc_bytes = 0)
        expect_identical(shared, t(shared), label = paste(what, "symmetric"))
        mxlp_expect_close(shared, kcall(), 1e-12, paste(what, "shared vs per-thread"))
        # one thread adds into the result itself, whatever the budget
        set_num_threads(1L)
        expect_identical(kcall(acc_bytes = 0), kcall(), label = paste(what, "one thread"))
      }
    }
  }
})

test_that("non-finite entries spread over a shared result as over per-thread ones", {
  skip_if_not(isTRUE(thread_info()$openmp_enabled), "needs OpenMP")
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(2L)
  same_pattern <- function(fx, k, what, cl = NULL) {
    sh <- mxlp_call(k, fx, cluster = cl, acc_bytes = 0)
    ref <- mxlp_call(k, fx, cluster = cl)
    what <- paste(what, k)
    expect_true(any(!is.finite(ref)), label = paste(what, "has non-finite entries"))
    expect_identical(is.na(sh), is.na(ref), label = paste(what, "NA"))
    expect_identical(is.nan(sh), is.nan(ref), label = paste(what, "NaN"))
    expect_identical(is.infinite(sh), is.infinite(ref), label = paste(what, "Inf"))
    expect_identical(sh[is.infinite(sh)], ref[is.infinite(ref)],
                     label = paste(what, "infinities and their signs"))
    # (a non-finite weight leaves no finite entry)
    if (any(is.finite(ref))) {
      mxlp_expect_close(sh[is.finite(ref)], ref[is.finite(ref)], 1e-12,
                        paste(what, "finite entries"))
    }
  }
  # draws at which a decision maker's choice probability is zero: NaN rows of
  # the Hessian (as in the non-finite Hessian test above)
  fx <- cdx_fixture("zero-weight draws", 1012, rc_dist = 0L)
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  row_off <- c(0L, cumsum(fx$M))
  rows_u <- (row_off[first[3L]] + 1L):row_off[first[3L] + fx$Ti[3L]]
  fx$W[rows_u, 1L] <- 0
  fx$W[row_off[first[3L]] + fx$choice_idx[first[3L]], 1L] <- 10
  fx$eta[1L, 1:5, 3L] <- -1e308
  same_pattern(fx, "hessian", "zero-weight draws")
  # utilities that overflow: a NaN score
  fx <- cdx_fixture("overflow", 1014)
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  row_off <- c(0L, cumsum(fx$M))
  fx$theta[1L] <- 2
  fx$X[row_off[first[2L]] + fx$choice_idx[first[2L]], 1L] <- -1e308
  for (k in c("bhhh", "meat")) same_pattern(fx, k, "overflowed unit", cl = c(1L, 1L, 2:9))
  # an infinite score entry with finite utilities
  fx <- cdx_fixture("infinite score", 1023)
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  row_off <- c(0L, cumsum(fx$M))
  u <- which(fx$Ti >= 2L)[1L]
  fx$theta[1L] <- 1e-306
  for (t in first[u] - 1L + seq_len(fx$Ti[u])) {
    other <- setdiff(seq_len(fx$M[t]), fx$choice_idx[t])[1L]
    fx$X[row_off[t] + other, 1L] <- 1e308
  }
  for (k in c("bhhh", "meat")) {
    same_pattern(fx, k, "infinite score", cl = rep(1:3, length.out = length(fx$Ti)))
  }
  # finite scores whose cluster sum overflows
  fx <- cdx_fixture("overflowing sum", 1025, panel = FALSE, U = 10L)
  row_off <- c(0L, cumsum(fx$M))
  fx$theta[1L] <- 1e-306
  for (t in 1:2) {
    other <- setdiff(seq_len(fx$M[t]), fx$choice_idx[t])[1L]
    fx$X[row_off[t] + other, 1L] <- 1e308
  }
  same_pattern(fx, "meat", "overflowing cluster sum", cl = c(1L, 1L, 2:9))
  # infinite, then missing, weights on one decision maker
  fx <- cdx_fixture("weights", 1015)
  first <- cumsum(c(1L, fx$Ti[-length(fx$Ti)]))
  for (w_bad in c(Inf, NA_real_)) {
    fx$weights[first[4L] - 1L + seq_len(fx$Ti[4L])] <- w_bad
    for (k in c("hessian", "bhhh", "meat")) {
      same_pattern(fx, k, sprintf("weight %s", w_bad), cl = rep(1:5, 2L))
    }
  }
})

test_that("the score matrix refuses to exceed max_bytes", {
  fx <- cdx_fixture("guard", 1019, U = 5L)
  expect_error(
    mxl_scores_parallel(fx$theta, fx$X, fx$W, fx$alt_idx, fx$choice_idx, fx$M,
                        fx$eta, fx$rc_dist, fx$rc_correlation, fx$rc_mean,
                        Ti = fx$Ti, max_bytes = 100),
    "The score matrix would take 0.00 GiB (5 likelihood units x", fixed = TRUE)
  expect_equal(dim(mxl_scores_parallel(fx$theta, fx$X, fx$W, fx$alt_idx,
                                       fx$choice_idx, fx$M, fx$eta, fx$rc_dist,
                                       fx$rc_correlation, fx$rc_mean, Ti = fx$Ti)),
               c(5L, length(fx$theta)))
})

# --- The variance routes --------------------------------------------------------

cdx_route_fit <- function(panel, draws, cluster = TRUE, seed = 1102) {
  sim <- simulate_mxl_data(N = if (panel) 60L else 120L, J = 6L,
                           T = if (panel) 2L else 1L, beta = c(0.8, -0.5),
                           Sigma = diag(c(0.6, 0.4)), outside_option = FALSE,
                           vary_choice_set = TRUE, seed = seed)
  dt <- data.table::as.data.table(sim$data)
  set.seed(seed)
  if (panel) {
    dt[, grp := pid %% 7L]
    w_p <- stats::runif(max(dt$pid), 0.5, 2)
    dt[, w := w_p[pid]]
  } else {
    dt[, grp := id %% 9L]
    w_i <- stats::runif(max(dt$id), 0.5, 2)
    dt[, w := w_i[id]]
  }
  args <- list(data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
               covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
               person_col = if (panel) "pid" else NULL, weights_col = "w",
               S = 20L, draws = draws, seed = 3L)
  if (cluster) args <- c(args, list(se_method = "cluster", cluster_col = "grp"))
  list(fit = suppressMessages(suppressWarnings(do.call(run_mxlogit, args))), dt = dt)
}

test_that("mixed logit variances never form the score matrix and build their draws once", {
  skip_on_cran()
  for (cfg in list(list(panel = TRUE, draws = "store"),
                   list(panel = FALSE, draws = "generate"))) {
    what <- sprintf("[%s, %s]", if (cfg$panel) "panel" else "cross-section", cfg$draws)
    fit <- cdx_route_fit(cfg$panel, cfg$draws)$fit
    expect_true(is.call(fit$call), label = paste(what, "the fit keeps its call"))
    # The old route, from the score matrix (non-uniform weights, so that w,
    # w^2 and sqrt(w) give different meats)
    d <- fit$data
    S_u <- choicer:::compute_scores(fit)
    w_u <- d$weights[choicer:::.unit_first(d)]
    expect_gt(length(unique(w_u)), 1L)
    A <- choicer:::.compute_bread(fit)
    old <- function(B) choicer:::.sandwich_combine(A, B)$vcov
    cl <- choicer:::.to_units(d$cluster, d, "cluster")
    ref <- list(
      robust = old(crossprod(w_u * S_u)),
      cluster = old(crossprod(rowsum(w_u * S_u, group = as.character(cl)))),
      bhhh = choicer:::invert_hessian(crossprod(sqrt(w_u) * S_u))$vcov)
    n_cubes <- 0L
    with_mocked_bindings({
      for (type in names(ref)) {
        n_cubes <- 0L
        V <- unname(vcov(fit, type = type))
        expect_equal(n_cubes, if (cfg$draws == "store") 1L else 0L,
                     label = sprintf("%s draws built for type = %s", what, type))
        expect_equal(V, unname(ref[[type]]), tolerance = 1e-10,
                     label = sprintf("%s vcov(type = %s) vs the score-matrix route",
                                     what, type))
      }
      expect_equal(unname(wesml_vcov(fit)), unname(ref$robust), tolerance = 1e-10,
                   label = paste(what, "wesml_vcov()"))
    },
    mxl_scores_parallel = function(...) stop("the score matrix was formed"),
    get_halton_normals = function(...) {
      n_cubes <<- n_cubes + 1L
      choicer:::.halton_cube(...)
    })
    expect_equal(unname(fit$vcov), unname(ref$cluster), tolerance = 1e-10,
                 label = paste(what, "fit-time cluster vcov"))
  }
})

test_that("mixed logit cluster labels are checked before any draws are built", {
  skip_on_cran()
  rf <- cdx_route_fit(TRUE, "store", cluster = FALSE)
  fit <- rf$fit
  dt <- rf$dt
  expect_error(vcov(fit, type = "cluster"),
               "Cluster-robust standard errors need cluster labels", fixed = TRUE)
  ids <- fit$data$situation_ids
  per_sit <- dt[!duplicated(id), setNames(pid %% 7L, id)]
  n_cubes <- 0L
  local_mocked_bindings(get_halton_normals = function(...) {
    n_cubes <<- n_cubes + 1L
    choicer:::.halton_cube(...)
  })
  bad_na <- per_sit
  # every situation of the first decision maker (prepared order: by person)
  bad_na[as.character(ids[seq_len(fit$data$Ti[1L])])] <- NA
  expect_error(vcov(fit, type = "cluster", cluster = bad_na),
               "`cluster` contains missing values.", fixed = TRUE)
  expect_equal(n_cubes, 0L, label = "no draws for labels with an NA")
  expect_error(suppressWarnings(vcov(fit, type = "cluster",
                                     cluster = unname(per_sit)[-1L])),
               "`cluster` has length", fixed = TRUE)
  varying <- per_sit
  varying[as.character(ids[1L])] <- 99L  # a decision maker's first situation only
  expect_error(vcov(fit, type = "cluster", cluster = varying),
               "must be constant within each decision maker", fixed = TRUE)
  expect_equal(n_cubes, 0L)
})
