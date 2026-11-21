# The mixed logit Hessian against a reference written from the estimand
# (mxlh_reference(), helper-mxl-hessian.R): the kernel forms every sum over a
# situation's alternatives and draws as matrix products of its probabilities
# (Identity D, §4.4.4 of the math article), on covariates centered on the
# alternative of largest posterior probability; the reference forms the
# per-draw outer products themselves. The fixtures are helper-mxl-compact.R's
# and helper-mxl-panel.R's.

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

test_that("the Hessian is the per-draw Louis identity on every cell and path", {
  on.exit(set_omp_threads(2L), add = TRUE)
  for (fx in mxlh_cells()) {
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
  # overflow at them. The rows and columns of the Cholesky parameters whose
  # factors are not finite there are NaN throughout, as the dense products
  # spread them; every other entry is the Hessian without those draws, which
  # the reference gives (it leaves draws of weight zero out).
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
  # the double range. The (zeta, zeta) terms omega e_r e_r2 y with y = sum_a P
  # w_2 w_p2 = 0 are then Inf x 0 in plain products, while the exact terms
  # are zero (as are the products (w e_r)(w e_r2) the reference forms).
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
  # y = sum_a P w_2 w_2 of order 1e-20: the exact term, near 1e298, is
  # finite, and the kernel sums such entries again from the factors'
  # significands and exponents.
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
  # pass the double range; the kernel forms them shifted by a power of two.
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
      mxlp_expect_close(mxlp_call("hessian", fy), H, 1e-12,
                        sprintf("[%s] W[, 1] x 2^%d", fx$name, sh))
      mxlp_expect_close(mxlp_call("hessian", fy, draw_batch = 2L), H, 1e-12,
                        sprintf("[%s] W[, 1] x 2^%d, draw batches", fx$name, sh))
    }
  }
})

test_that("shifts that cancel in a pair leave the Hessian as it is", {
  # Every |w| lies in [1/2, 1), so with W[, 1] x 2^200 and W[, 2] x 2^-200
  # each unit shifts the two columns by 2^200 and 2^-200 exactly: their
  # cross products enter the (zeta, zeta) sums unscaled (shifts of sum zero),
  # their own products scaled. L_11 / 2^200, L_21 x 2^200 and L_22 x 2^200
  # leave every utility as it is; the Hessian is the same in log L_11 and
  # log L_22, and its L_21 rows and columns are divided by 2^200.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  fx <- cdx_fixture("xsec, opposite shifts", 1181, panel = FALSE,
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
  for (batch in c(0L, 2L)) {
    Hy <- mxlp_call("hessian", fy, draw_batch = batch)
    mxlp_expect_close(Hy * outer(d, d), H, 1e-12,
                      sprintf("opposite shifts, draw batch %d", batch))
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
