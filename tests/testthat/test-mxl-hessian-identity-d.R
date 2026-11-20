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
