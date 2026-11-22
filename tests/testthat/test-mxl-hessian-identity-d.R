# The mixed logit Hessian against a reference written from the estimand
# (mxlh_reference(), helper-mxl-hessian.R): the kernel forms every sum over a
# situation's alternatives and draws as matrix products of its probabilities
# (Identity D, §4.4.4 of the math article), on covariates centered on the
# alternative of largest posterior probability, and a decision maker outside
# Identity D's range by the direct formula, the per-alternative, per-draw
# outer products of the utility's derivatives; the reference forms those
# outer products itself. The fixtures are helper-mxl-compact.R's and
# helper-mxl-panel.R's.

mxlh_cells <- function() {
  c(cdx_cells(), list(
    cdx_fixture("xsec, log-normal, uncorrelated, with means", 1101, panel = FALSE,
                rc_dist = c(1L, 1L), rc_correlation = FALSE, rc_mean = TRUE),
    cdx_fixture("xsec, outside option, log-normal, correlated", 1102, panel = FALSE,
                include_outside_option = TRUE, rc_dist = c(1L, 1L)),
    cdx_fixture("panel, three correlated coefficients", 1103,
                rc_dist = c(0L, 1L, 0L)),
    cdx_fixture("panel, no fixed coefficients", 1104, K_x = 0L),
    cdx_fixture("panel, no fixed coefficients, no means", 1105, K_x = 0L,
                rc_mean = FALSE),
    cdx_fixture("panel, repeated alternatives, outside option", 1106, dup = TRUE,
                include_outside_option = TRUE),
    mxlp_fixture("panel, alternative-level W", 1107, rc_dist = c(0L, 1L),
                 rc_correlation = TRUE, rc_mean = TRUE, W_layout = "alt",
                 weight_type = "person")
  ))
}

# The draws generate mode forms for each likelihood unit (gen_seed 0,
# identity permutations, as mxlp_call(generate = TRUE) asks for them).
mxlh_generated_eta <- function(fx) {
  U <- if (is.null(fx$Ti)) length(fx$M) else length(fx$Ti)
  eta <- array(0, c(fx$K_w, fx$S, U))
  for (u in seq_len(U)) {
    eta[, , u] <- halton_fill_block((u - 1) * fx$S + 1, fx$S, fx$K_w, 0, 0L)
  }
  eta
}

# Whether mxl_hessian_parallel() forms each likelihood unit of fx by Identity
# D (TRUE) or, outside its range, by the direct formula (FALSE): the kernel's
# own check, per unit.
mxlh_identity_d_units <- function(fx, theta = fx$theta, eta = fx$eta) {
  args <- list(theta = theta, X = fx$X, W = fx$W, alt_idx = fx$alt_idx,
               choice_idx = fx$choice_idx, M = fx$M, weights = fx$weights,
               eta_draws = eta, rc_dist = fx$rc_dist,
               rc_correlation = fx$rc_correlation, rc_mean = fx$rc_mean,
               use_asc = fx$use_asc,
               include_outside_option = fx$include_outside_option)
  if (!is.null(fx[["Ti"]])) args$Ti <- as.integer(fx$Ti)
  do.call(mxl_hessian_identity_d_units, args)
}

# The likelihood unit of each situation of fx
mxlh_unit_of <- function(fx) {
  if (is.null(fx[["Ti"]])) seq_along(fx$M) else rep(seq_along(fx$Ti), fx$Ti)
}

# fx with one more fixed covariate, of coefficient zero, whose spread over the
# alternatives of one situation (the first with two) is `spread`: no utility
# moves, and that situation's decision maker (`unit`) lies inside Identity
# D's range (spreads up to 2^64) or outside it.
mxlh_with_spread <- function(fx, spread) {
  t <- which(fx$M >= 2L)[1L]
  x <- numeric(nrow(fx$X))
  x[c(0L, cumsum(fx$M))[t] + 2L] <- spread
  K_x <- ncol(fx$X)
  fx$X <- unname(cbind(fx$X, x))
  fx$theta <- append(fx$theta, 0, after = K_x)
  fx$K_x <- K_x + 1L
  fx$unit <- mxlh_unit_of(fx)[t]
  fx
}

# Two alternatives per situation, one random coefficient, no fixed
# covariates or constants: `w` gives W's rows, `eta` the K_w x S x U draws.
mxlh_binary <- function(w, chosen, weights, eta, theta, rc_dist = 0L,
                        rc_mean = FALSE, Ti = NULL) {
  n <- length(chosen)
  list(name = "binary", X = matrix(numeric(), 2L * n, 0L), W = matrix(w),
       alt_idx = rep(1:2, n), choice_idx = as.integer(chosen), M = rep(2L, n),
       weights = weights, Ti = if (is.null(Ti)) NULL else as.integer(Ti),
       eta = eta, theta = theta, rc_dist = as.integer(rc_dist),
       rc_correlation = FALSE, rc_mean = rc_mean, use_asc = FALSE,
       include_outside_option = FALSE, K_w = 1L, K_x = 0L, J = 2L,
       S = dim(eta)[2L])
}

test_that("the Hessian is the per-draw Louis identity on every cell and path", {
  on.exit(set_omp_threads(2L), add = TRUE)
  for (fx in mxlh_cells()) {
    expect_true(all(mxlh_identity_d_units(fx)),
                label = sprintf("[%s] every unit in Identity D's range", fx$name))
    ref <- mxlh_reference(fx$theta, fx)
    set_omp_threads(1L)
    H1 <- mxlp_call("hessian", fx)
    mxlp_expect_close(H1, ref, 1e-12, sprintf("[%s] Hessian vs the reference", fx$name))
    # Draw batches of two draws: every situation's probabilities formed again
    mxlp_expect_close(mxlp_call("hessian", fx, draw_batch = 2L), H1, 1e-13,
                      sprintf("[%s] draw batches", fx$name))
    ref_gen <- mxlh_reference(fx$theta, fx, eta = mxlh_generated_eta(fx))
    mxlp_expect_close(mxlp_call("hessian", fx, generate = TRUE), ref_gen, 1e-12,
                      sprintf("[%s] generate mode vs the reference", fx$name))
    set_omp_threads(2L)
    mxlp_expect_close(mxlp_call("hessian", fx), H1, 1e-13,
                      sprintf("[%s] 2 threads", fx$name))
  }
})

test_that("a common level in a covariate leaves the Hessian as it is", {
  # Without an outside option, adding a constant to a column of X or W for
  # every row shifts every alternative's utility in a situation alike, so the
  # probabilities, the scores and the Hessian do not change. Formed from the
  # covariates' levels, the terms Identities B and C subtract grow with the
  # square of the constant while their difference does not; centered, the
  # kernel loses only the rounding of the shifted covariates themselves.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (fx in list(cdx_fixture("xsec", 1111, panel = FALSE),
                  cdx_fixture("panel", 1112, rc_dist = c(0L, 1L)))) {
    H <- mxlp_call("hessian", fx)
    fy <- fx
    fy$X[, 1L] <- fy$X[, 1L] + 1e5
    fy$W[, 1L] <- fy$W[, 1L] + 1e3
    mxlp_expect_close(mxlp_call("hessian", fy), H, 1e-9,
                      sprintf("[%s] Hessian with X[, 1] + 1e5 and W[, 1] + 1e3", fx$name))
  }
})

test_that("draws of weight zero leave the other entries exact", {
  # Draws 1-5 of decision maker u give its chosen alternative a utility of
  # -Inf (W[, 1] is 10 there and 0 elsewhere in u's rows, eta_1 = -1e308), so
  # those draws have weight zero, while the Cholesky parameters' draw factors
  # overflow at them, past Identity D's range: u takes the direct formula,
  # which leaves such draws out of Identity C and H_V as Identity D does. The
  # rows and columns of the Cholesky parameters whose factors are not finite
  # there are NaN throughout, as the dense products spread them; every other
  # entry is the Hessian without those draws, which the reference gives (it
  # leaves draws of weight zero out).
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (ioo in c(FALSE, TRUE)) {
    fx <- cdx_fixture("zero-weight draws", 1121 + ioo, panel = FALSE,
                      include_outside_option = ioo, rc_dist = c(0L, 0L))
    u <- 3L
    row_off <- c(0L, cumsum(fx$M))
    if (fx$choice_idx[u] == 0L) fx$choice_idx[u] <- 1L
    rows_u <- (row_off[u] + 1L):row_off[u + 1L]
    fx$W[rows_u, 1L] <- 0
    fx$W[row_off[u] + fx$choice_idx[u], 1L] <- 10
    fx$eta[1L, 1:5, u] <- -1e308
    expect_identical(mxlh_identity_d_units(fx), seq_along(fx$M) != u)
    H <- mxlp_call("hessian", fx)
    bad <- which(apply(is.nan(H), 1L, all))
    expect_gt(length(bad), 0L)
    # Only Cholesky parameters (after beta and mu; K_w = 2, correlated)
    expect_true(all(bad %in% (fx$K_x + fx$K_w + seq_len(3L))))
    ref <- mxlh_reference(fx$theta, fx)
    expect_true(all(is.finite(ref)))
    mxlp_expect_close(H[-bad, -bad], ref[-bad, -bad], 1e-12,
                      sprintf("[outside option = %s] the entries outside the NaN rows", ioo))
  }
})

test_that("a draw factor whose square overflows leaves the Hessian finite", {
  # Decision maker u's covariates of the log-normal coefficient are zero, and
  # at its first draw gamma_2 = 400, so that the coefficient's draw factors
  # e_r = exp(gamma_2) dL_2q eta_q (about 1e176) are finite but square past
  # the double range. Identity D's (zeta, zeta) terms omega e_r e_r2 y, with
  # y = sum_a P w_2 w_p2 = 0, would be Inf x 0; such factors lie outside its
  # range, so u takes the direct formula, whose products (w e_r)(w e_r2) are
  # zero, as the exact terms are (and the reference's products).
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (ioo in c(FALSE, TRUE)) {
    fx <- cdx_fixture("an overflowing draw factor", 1141 + ioo, panel = FALSE,
                      include_outside_option = ioo, rc_dist = c(0L, 1L))
    u <- 3L
    row_off <- c(0L, cumsum(fx$M))
    fx$W[(row_off[u] + 1L):row_off[u + 1L], 2L] <- 0
    L <- mxlp_unpack(fx$theta, fx)$L
    fx$eta[2L, 1L, u] <- (400 - L[2L, 1L] * fx$eta[1L, 1L, u]) / L[2L, 2L]
    expect_identical(mxlh_identity_d_units(fx), seq_along(fx$M) != u)
    H <- mxlp_call("hessian", fx)
    expect_true(all(is.finite(H)))
    mxlp_expect_close(H, mxlh_reference(fx$theta, fx), 1e-12,
                      sprintf("[outside option = %s] Hessian vs the reference", ioo))
  }
})

test_that("an overflowing draw factor against small covariates stays exact", {
  # Two normal coefficients. At decision maker u's first draw eta_1 = 1e160
  # and L_21 = 1e-160, so the off-diagonal entry's draw factor e_21 = eta_1 is
  # 1e160 while gamma_2 = L_21 eta_1 + L_22 eta_2 stays of order one. u's
  # covariates of the first coefficient are zero and those of the second of
  # order 1e-10, so omega e_21 e_21 passes the double range against
  # y = sum_a P w_2 w_2 of order 1e-20 while the exact term, near 1e298, is
  # finite. The factor lies outside Identity D's range: u takes the direct
  # formula, which squares w_2 e_21 itself, and whose H_V has no f'' term for
  # a normal coefficient (0 x e_21^2 would be NaN).
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  fx <- cdx_fixture("an overflowing draw factor, small covariates", 1171,
                    panel = FALSE, rc_dist = c(0L, 0L))
  u <- 3L
  rows_u <- (c(0L, cumsum(fx$M))[u] + 1L):cumsum(fx$M)[u]
  fx$W[rows_u, 1L] <- 0
  fx$W[rows_u, 2L] <- 1e-10 * fx$W[rows_u, 2L]
  fx$eta[1L, 1L, u] <- 1e160
  i_21 <- fx$K_x + fx$K_w + 2L      # theta's L_21 (after L_11)
  fx$theta[i_21] <- 1e-160
  expect_identical(mxlh_identity_d_units(fx), seq_along(fx$M) != u)
  H <- mxlp_call("hessian", fx)
  ref <- mxlh_reference(fx$theta, fx)
  expect_true(all(is.finite(H)))
  expect_gt(abs(ref[i_21, i_21]), 1e250)
  expect_lt(abs(H[i_21, i_21] / ref[i_21, i_21] - 1), 1e-12)
  mxlp_expect_close(H, ref, 1e-12, "the kernel vs the reference")
})

test_that("covariates of extreme magnitude leave the Hessian as it is", {
  # W[, 1] multiplied by 2^530 (or 2^-530) and its coefficient's Cholesky
  # entry L_11 divided by as much leave every utility, and the Hessian in
  # theta's parameters (log L_11 on the diagonal), as they are: no means, so
  # that no parameter multiplies w_1 alone. Such covariates' products w_1 w_p
  # pass the double range. The covariates' spreads (2^530) or the draw
  # factors L_11 eta (2^-530) lie outside Identity D's range, so every
  # decision maker takes the direct formula, which multiplies w_1 by its draw
  # factor first, and whose H_V has no f'' term (0 x Inf where L_11 is of
  # order 2^530).
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (fx in list(cdx_fixture("xsec", 1151, panel = FALSE, rc_mean = FALSE,
                              rc_dist = c(0L, 0L)),
                  cdx_fixture("panel, outside option", 1152, rc_mean = FALSE,
                              rc_dist = c(0L, 0L), include_outside_option = TRUE))) {
    H <- mxlp_call("hessian", fx)
    i_11 <- fx$K_x + 1L
    for (sh in c(530, -530)) {
      fy <- fx
      fy$W[, 1L] <- fy$W[, 1L] * 2^sh
      fy$theta[i_11] <- fy$theta[i_11] - sh * log(2)
      expect_false(any(mxlh_identity_d_units(fy)))
      mxlp_expect_close(mxlp_call("hessian", fy), H, 1e-12,
                        sprintf("[%s] W[, 1] x 2^%d", fx$name, sh))
      mxlp_expect_close(mxlp_call("hessian", fy, draw_batch = 2L), H, 1e-12,
                        sprintf("[%s] W[, 1] x 2^%d, draw batches", fx$name, sh))
    }
  }
})

test_that("covariates scaled in opposite directions leave the Hessian as it is", {
  # W[, 1] x 2^200 and W[, 2] x 2^-200 (every |w| in [1/2, 1) before), with
  # L_11 / 2^200, L_21 x 2^200 and L_22 x 2^200, leave every utility as it
  # is; the Hessian is the same in log L_11 and log L_22, and its L_21 rows
  # and columns are divided by 2^200. The spreads of W[, 1] and the draw
  # factors of L_21 and L_22 lie outside Identity D's range, so every
  # decision maker takes the direct formula, whose cross products of the two
  # coefficients' derivatives are of order one.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  fx <- cdx_fixture("xsec, opposite scales", 1181, panel = FALSE,
                    rc_mean = FALSE, rc_dist = c(0L, 0L))
  set.seed(1182)
  n_w <- length(fx$W)
  fx$W[] <- sample(c(-1, 1), n_w, replace = TRUE) * stats::runif(n_w, 0.5, 0.99)
  H <- mxlp_call("hessian", fx)
  i <- fx$K_x + 1:3                  # log L_11, L_21, log L_22
  fy <- fx
  fy$W[, 1L] <- fy$W[, 1L] * 2^200
  fy$W[, 2L] <- fy$W[, 2L] * 2^-200
  fy$theta[i] <- fy$theta[i] + c(-200 * log(2), 0, 200 * log(2))
  fy$theta[i[2L]] <- fy$theta[i[2L]] * 2^200
  d <- rep(1, length(fx$theta))
  d[i[2L]] <- 2^200                  # the L_21 rows back to fx's scale
  expect_true(all(mxlh_identity_d_units(fx)))
  expect_false(any(mxlh_identity_d_units(fy)))
  for (batch in c(0L, 2L)) {
    Hy <- mxlp_call("hessian", fy, draw_batch = batch)
    mxlp_expect_close(Hy * outer(d, d), H, 1e-12,
                      sprintf("opposite scales, draw batch %d", batch))
  }
})

test_that("a decision maker outside Identity D's range takes the direct formula, on every cell", {
  # A fixed covariate of coefficient zero spread by 2^70 over one situation's
  # alternatives puts that situation's decision maker outside the range
  # (2^64) without moving a utility: its Hessian is the direct formula's,
  # the others' Identity D's, and the sum the reference's. At a spread of
  # 2^64 it stays inside.
  on.exit(set_omp_threads(2L), add = TRUE)
  for (fx0 in mxlh_cells()) {
    for (k in c(64L, 70L)) {
      fx <- mxlh_with_spread(fx0, 2^k)
      what <- sprintf("[%s] a spread of 2^%d", fx$name, k)
      inside <- rep(TRUE, max(mxlh_unit_of(fx)))
      if (k > 64L) inside[fx$unit] <- FALSE
      expect_identical(mxlh_identity_d_units(fx), inside,
                       label = paste(what, "units in Identity D's range"))
      set_omp_threads(1L)
      H1 <- mxlp_call("hessian", fx)
      mxlp_expect_close(H1, mxlh_reference(fx$theta, fx), 1e-12,
                        paste(what, "vs the reference"))
      mxlp_expect_close(mxlp_call("hessian", fx, draw_batch = 2L), H1, 1e-13,
                        paste(what, "in draw batches"))
      set_omp_threads(2L)
      mxlp_expect_close(mxlp_call("hessian", fx), H1, 1e-13, paste(what, "on 2 threads"))
    }
  }
})

test_that("a mean derivative past the range keeps the mean's information", {
  # A log-normal coefficient's mean mu = -log(1e-160) gives mu*' = e^mu =
  # 1e160 against covariates of 1e-160: the mean's utility derivatives w mu*'
  # are of order one while mu*' squared overflows. Every decision maker takes
  # the direct formula; four of them, at draws of 0 and 1.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  w <- 1e-160
  p <- stats::plogis(c(1, 2))
  fx <- mxlh_binary(rep(c(0, w), 4L), rep(1:2, 2L),
                    c(1 - p[1L], p[1L], 1 - p[2L], p[2L]),
                    array(rep(c(0, 0, 1, 1), each = 2L), c(1L, 2L, 4L)),
                    c(-log(w), log(-log(w))), rc_dist = 1L, rc_mean = TRUE)
  expect_false(any(mxlh_identity_d_units(fx)))
  H <- mxlp_call("hessian", fx)
  expect_true(all(is.finite(H)))
  mxlp_expect_close(H, mxlh_reference(fx$theta, fx), 1e-12, "the kernel vs the reference")
})

test_that("a large covariate in one task leaves a panel's other tasks", {
  # Two decision makers with two tasks each; the second tasks carry the
  # information. The first tasks' covariates are 0 and -2^510, so that their
  # first alternative's probability is one: a spread of 2^510, past Identity
  # D's range. The second tasks' covariates are 0 and 2^-30, and the weights
  # scale their information to P (1 - P) in all.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  w <- 2^-30
  p <- stats::plogis(w)
  fx <- mxlh_binary(rep(c(0, -2^510, 0, w), 2L), c(1L, 1L, 1L, 2L),
                    rep(c(1 - p, p) / w^2, each = 2L), array(1, c(1L, 2L, 2L)), 0,
                    Ti = c(2L, 2L))
  expect_false(any(mxlh_identity_d_units(fx)))
  H <- mxlp_call("hessian", fx)
  mxlp_expect_close(H, mxlh_reference(fx$theta, fx), 1e-12, "a spread of 2^510")
  expect_equal(H[1L, 1L], p * (1 - p), tolerance = 1e-12)
  # Second tasks' covariates 0 and 1, weights below one; the first tasks'
  # covariates (c, c) carry no information at any level c. A common level
  # spreads nothing, so the decision makers stay in the range, and centered,
  # it leaves the Hessian as it is.
  p <- stats::plogis(1)
  fy <- mxlh_binary(rep(c(0, 0, 0, 1), 2L), c(1L, 1L, 1L, 2L),
                    rep(c(1 - p, p), each = 2L), array(1, c(1L, 2L, 2L)), 0,
                    Ti = c(2L, 2L))
  H0 <- mxlp_call("hessian", fy)
  mxlp_expect_close(H0, mxlh_reference(fy$theta, fy), 1e-12, "no common level")
  for (c0 in c(2^510, 2^600, 3e150)) {
    fy$W[c(1L, 2L, 5L, 6L), 1L] <- c0
    expect_true(all(mxlh_identity_d_units(fy)))
    mxlp_expect_close(mxlp_call("hessian", fy), H0, 1e-12,
                      sprintf("a common level of %.3g", c0))
  }
})

test_that("a covariate spread that overflows takes the direct formula", {
  # Covariates of -1e308 and 1e308 against L = 1e-308: the utilities are
  # near -1 and 1, while the covariates' difference overflows.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  p <- stats::plogis(2)
  fx <- mxlh_binary(rep(c(-1e308, 1e308), 2L), 1:2, c(1 - p, p),
                    array(1, c(1L, 1L, 2L)), log(1e-308))
  expect_false(any(mxlh_identity_d_units(fx)))
  H <- mxlp_call("hessian", fx)
  expect_true(all(is.finite(H)))
  mxlp_expect_close(H, mxlh_reference(fx$theta, fx), 1e-12, "the kernel vs the reference")
})

test_that("a covariate below the range's floor keeps a small parameter's information", {
  # One decision maker with two mirrored binary tasks: a fixed covariate (0, 1)
  # and a log-normal coefficient's covariates (0, c), L = 37, draws 1 and 0.97,
  # so draw factors near 4e17 and 1e17. At c = 2^-540 the coefficient's
  # utility derivatives are near 2^-482 and its information near 4e-291,
  # while c^2 underflows: below the floor of Identity D's range (2^-128), the
  # decision maker takes the direct formula. Small as it is, that information
  # is correlated with beta's: with it beta's standard error is 3.09, without
  # it 1.06, and with its sign reversed the information is not positive
  # definite. At 2^-129 the decision maker is outside the range too, at 2^-128
  # inside. The information in closed form: at every draw each task is a
  # binary logit at P = 1/2 (beta = 0, and c e_s is far below a rounding of
  # 1/2), so it adds P (1 - P) z z' with z = (1, c e_s), the derivatives of
  # the second alternative's utility; the draws weigh 1/2 each, the scores
  # of the two tasks cancel, and so do their H_V terms. (The reference
  # helper adds those H_V terms, near 1e-144, to each task's 1e-291 before
  # they cancel, and loses it.)
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  se_beta <- function(I) sqrt(1 / (I[1L, 1L] - I[1L, 2L]^2 / I[2L, 2L]))
  for (k in c(540L, 129L, 128L)) {
    c0 <- 2^-k
    fx <- list(name = "two mirrored tasks", X = matrix(c(0, 1, 0, 1)),
               W = matrix(c(0, c0, 0, c0)), alt_idx = rep(1:2, 2L),
               choice_idx = c(1L, 2L), M = c(2L, 2L), weights = c(1, 1), Ti = 2L,
               eta = array(c(1, 0.97), c(1L, 2L, 1L)), theta = c(0, log(37)),
               rc_dist = 1L, rc_correlation = FALSE, rc_mean = FALSE,
               use_asc = FALSE, include_outside_option = FALSE, K_w = 1L,
               K_x = 1L, J = 2L, S = 2L)
    expect_identical(mxlh_identity_d_units(fx), k <= 128L)
    H <- mxlp_call("hessian", fx)
    eta <- fx$eta[1L, , 1L]
    z <- c0 * exp(37 * eta) * 37 * eta      # d utility / d log L, per draw
    I <- 0.5 * matrix(c(1, mean(z), mean(z), mean(z^2)), 2L)
    expect_lt(max(abs(H / I - 1)), 1e-12, label = sprintf("c = 2^-%d: entrywise error", k))
    expect_equal(se_beta(H), se_beta(I), tolerance = 1e-12)
    expect_gt(se_beta(I), 3)
  }
})

test_that("a choice whose probability underflows at every draw counts", {
  # Utilities (0, 0, -800) with the third alternative chosen: its probability
  # is zero in floating point at every draw, while its log, -800 - log 2, and
  # the draw weights are finite. The information on the coefficient of the
  # second alternative's indicator is P_2 (1 - P_2) = 1/4 (W is zero, so the
  # draws play no part).
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  fx <- cdx_fixture("an underflowing choice", 1161, panel = FALSE, U = 1L,
                    J = 3L, m_range = c(3L, 3L), rc_dist = 0L, use_asc = FALSE)
  fx$X[] <- c(0, 1, 0, 0, 0, -800)  # alternative 2's indicator; -800 on 3
  fx$W[] <- 0
  fx$theta[1:2] <- c(0, 1)
  fx$choice_idx[] <- 3L
  fx$weights[] <- 1
  ref <- mxlh_reference(fx$theta, fx)
  expect_equal(ref[1L, 1L], 0.25)
  mxlp_expect_close(mxlp_call("hessian", fx), ref, 1e-12, "the kernel vs the reference")
  # The same among other decision makers: u's chosen utility lowered by 800
  fx <- cdx_fixture("xsec", 1162, panel = FALSE)
  u <- 3L
  row_c <- c(0L, cumsum(fx$M))[u] + fx$choice_idx[u]
  fx$theta[1L] <- 1
  fx$X[row_c, 1L] <- fx$X[row_c, 1L] - 800
  mxlp_expect_close(mxlp_call("hessian", fx), mxlh_reference(fx$theta, fx), 1e-12,
                    "a decision maker whose choice underflows, among others")
})

test_that("a covariate that varies within no situation has exactly zero rows", {
  # Without an outside option a covariate with one value per situation (a
  # decision maker's characteristic entered alone) cancels from every utility
  # difference: its parameters are not identified and their rows of the exact
  # Hessian are zero. Centered on one of the situation's own alternatives, the
  # kernel's covariates are exact zeros there, so those rows are exact zeros
  # too and the information is singular. With an outside option the covariate
  # moves the inside alternatives against it and is identified.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (ioo in c(FALSE, TRUE)) {
    fx <- cdx_fixture("situation-level covariates", 1131 + ioo,
                      include_outside_option = ioo)
    sit <- rep(seq_along(fx$M), fx$M)
    fx$X[, 1L] <- stats::rnorm(length(fx$M))[sit]
    fx$W[, 1L] <- stats::rnorm(length(fx$M))[sit]
    unid <- c(1L, fx$K_x + 1L, fx$K_x + fx$K_w + 1L)  # beta_1, mu_1, l_11
    H <- mxlp_call("hessian", fx)
    mxlp_expect_close(H, mxlh_reference(fx$theta, fx), 1e-12,
                      sprintf("[outside option = %s] Hessian vs the reference", ioo))
    if (ioo) {
      expect_true(all(diag(H)[unid] != 0))
    } else {
      expect_true(all(H[unid, ] == 0) && all(H[, unid] == 0))
    }
  }
})
