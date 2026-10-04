# tests/testthat/test-halton-generator.R
#
# Phase-A deterministic correctness gates for src/halton.h.
#
# Tested via the thin Rcpp wrappers in src/halton_test_exports.cpp:
#   choicer:::halton_radical_inverse(n, base)
#   choicer:::halton_inv_normal_cdf(p)
#   choicer:::halton_generate_uniform(n, dim, seed, scramble)
#   choicer:::halton_generate_normal(S, N, K_w, seed, scramble)
#   choicer:::halton_fill_block(n0, S, K_w, seed, scramble)
#   choicer:::halton_reference_block(n0, S, K_w, seed, scramble)
#
# These wrappers are @noRd (internal); they live in the package namespace
# and are accessible via ::: in testthat / devtools::test.
#
# Do NOT add hash/bit-equality assertions on randomized draws across
# platforms.  Estimation-level tests are Phase B.  (Section 6 compares two
# code paths of one build bit for bit, not draws across platforms.)

# ---------------------------------------------------------------------------
# 1. radical_inverse: known values
# ---------------------------------------------------------------------------

test_that("radical_inverse returns correct known values", {
  # phi_2(1) = 0.5   (binary: 0.1  -> 1/2)
  expect_equal(choicer:::halton_radical_inverse(1, 2), 0.5)

  # phi_2(2) = 0.25  (binary: 10   -> 0.01 = 1/4)
  expect_equal(choicer:::halton_radical_inverse(2, 2), 0.25)

  # phi_2(3) = 0.75  (binary: 11   -> 0.11 = 3/4)
  expect_equal(choicer:::halton_radical_inverse(3, 2), 0.75)

  # phi_3(1) = 1/3
  expect_equal(choicer:::halton_radical_inverse(1, 3), 1 / 3)

  # phi_3(2) = 2/3
  expect_equal(choicer:::halton_radical_inverse(2, 3), 2 / 3)

  # phi_2(0) = 0 (loop body never executes)
  expect_equal(choicer:::halton_radical_inverse(0, 2), 0.0)

  # phi_5(4) = 4/5
  expect_equal(choicer:::halton_radical_inverse(4, 5), 4 / 5)
})

# ---------------------------------------------------------------------------
# 2. inv_normal_cdf: max absolute error < 1e-12 vs stats::qnorm
# ---------------------------------------------------------------------------

test_that("inv_normal_cdf achieves max abs error < 1e-12 vs stats::qnorm", {
  # 1600-point dense grid spanning the central, intermediate-tail, and
  # far-tail regimes.
  p_central  <- seq(1e-10, 1 - 1e-10, length.out = 1600)
  our_vals   <- vapply(p_central, choicer:::halton_inv_normal_cdf, numeric(1))
  r_vals     <- stats::qnorm(p_central)
  max_err    <- max(abs(our_vals - r_vals))

  expect_true(
    max_err < 1e-12,
    label = paste0("max abs error = ", max_err, " (gate: < 1e-12)")
  )

  # NIT-1 (Phase-A review): explicitly probe the E/F far-tail branch.
  # The far-tail branch activates when sqrt(-log(min(p,1-p))) > 5,
  # which requires p < exp(-25) ~ 1.39e-11.  The grid above starts at
  # 1e-10, so the E/F coefficients are not exercised by the grid test.
  # These two points exercise the far-tail path directly.
  far_tail_p  <- c(1e-12, 1e-13)
  far_tail_r  <- stats::qnorm(far_tail_p)
  far_tail_us <- vapply(far_tail_p, choicer:::halton_inv_normal_cdf, numeric(1))
  far_tail_err <- max(abs(far_tail_us - far_tail_r))
  expect_true(
    far_tail_err < 1e-12,
    label = paste0("far-tail max abs error = ", far_tail_err, " (gate: < 1e-12)")
  )
})

test_that("inv_normal_cdf handles edge-case sentinels without NaN", {
  # p <= 0 returns -8.29 (sentinel, not -Inf)
  expect_equal(choicer:::halton_inv_normal_cdf(0.0),  -8.29)
  expect_equal(choicer:::halton_inv_normal_cdf(-0.5), -8.29)

  # p >= 1 returns +8.29 (sentinel, not +Inf)
  expect_equal(choicer:::halton_inv_normal_cdf(1.0),   8.29)
  expect_equal(choicer:::halton_inv_normal_cdf(1.5),   8.29)
})

# ---------------------------------------------------------------------------
# 3. compat mode (scramble=0): randtoolbox::halton's points (bit for bit
#    where both are compiled with the same floating-point contraction; the
#    tolerance covers builds that differ, test-mxl-predict-store-points.R)
# ---------------------------------------------------------------------------

test_that("halton_generate_uniform (scramble=0) matches randtoolbox", {
  skip_if_not_installed("randtoolbox")

  tol <- 1e-12  # the same algorithm: equal up to the last bit

  # Dims 1..10, n = 100
  for (d in 1:10) {
    rt   <- randtoolbox::halton(100, d, normal = FALSE)
    ours <- choicer:::halton_generate_uniform(100, d, seed = 0, scramble = 0)
    expect_true(
      max(abs(rt - ours)) <= tol,
      label = paste0("compat match, dim=", d, ", n=100")
    )
  }

  # Dims 1, 5, 10 with n = 2000
  for (d in c(1, 5, 10)) {
    rt   <- randtoolbox::halton(2000, d, normal = FALSE)
    ours <- choicer:::halton_generate_uniform(2000, d, seed = 0, scramble = 0)
    expect_true(
      max(abs(rt - ours)) <= tol,
      label = paste0("compat match, dim=", d, ", n=2000")
    )
  }
})

# ---------------------------------------------------------------------------
# 4. Position-wise permutation mode (scramble=1): reproducibility and sanity
# ---------------------------------------------------------------------------

test_that("permuted mode is bitwise reproducible with the same seed", {
  a1 <- choicer:::halton_generate_normal(100, 10, 3, seed = 42, scramble = 1)
  a2 <- choicer:::halton_generate_normal(100, 10, 3, seed = 42, scramble = 1)
  expect_identical(a1, a2)
})

test_that("permuted mode produces different draws with a different seed", {
  a <- choicer:::halton_generate_normal(100, 10, 3, seed = 42,  scramble = 1)
  b <- choicer:::halton_generate_normal(100, 10, 3, seed = 999, scramble = 1)
  expect_false(isTRUE(all.equal(a, b)))
})

test_that("permuted mode leaves implicit trailing zero digits unchanged", {
  # For n=1 in base 2 the implementation visits one explicit digit and then
  # stops. Across seeds the first coordinate can therefore only be 0 or 1/2;
  # it does not acquire a randomized binary tail. This pins the precise
  # position-wise construction and guards against calling it nested Owen.
  first_points <- vapply(0:63, function(seed) {
    choicer:::halton_generate_uniform(
      n = 1L, dim = 1L, seed = seed, scramble = 1L
    )[1L, 1L]
  }, numeric(1L))
  expect_true(all(first_points %in% c(0, 0.5)))
  expect_setequal(unique(first_points), c(0, 0.5))
})

test_that("permuted mode per-dim mean is approximately 0 and var approximately 1", {
  # Use N=2000 individuals with S=1 draw each per dimension.
  # Low-discrepancy sequences have good Monte Carlo statistics at N=2000.
  # NIT-2 (Phase-A review): the gate |mean| < 0.1 and |var-1| < 0.3 is
  # intentionally loose (sanity-only, not precision). The real correctness
  # gate is the quantile-match test (slice-indexing) below, which uses a
  # 1e-10 tolerance. The loose mean/var gate will not pass a broken
  # implementation that the quantile test would catch.
  S <- 1L; N <- 2000L; K_w <- 5L
  m <- choicer:::halton_generate_normal(S, N, K_w, seed = 42, scramble = 1)
  # m is K_w x (S*N) = 5 x 2000
  for (k in seq_len(K_w)) {
    row_vals <- m[k, ]
    expect_true(
      abs(mean(row_vals)) < 0.1,
      label = paste0("permuted dim ", k, " mean approx 0: mean=", round(mean(row_vals), 4))
    )
    expect_true(
      abs(var(row_vals) - 1.0) < 0.3,
      label = paste0("permuted dim ", k, " var approx 1: var=", round(var(row_vals), 4))
    )
  }
})

test_that("permuted-mode indexing uses n = (i-1)*S + s + 1", {
  # In compat mode (scramble=0), verify that the K_w x (S*N) output of
  # halton_generate_normal is consistent with the per-point formula:
  #   column (i-1)*S + s (0-based) = inv_normal_cdf(phi_b((i-1)*S + s + 1))
  # We check a representative selection of (i, s, k) triples.
  skip_if_not_installed("randtoolbox")

  S <- 10L; N <- 3L; K_w <- 3L
  full_normal <- choicer:::halton_generate_normal(S, N, K_w, seed = 0, scramble = 0)
  # full_normal is K_w x (S*N); randtoolbox provides reference uniform values.
  rt_uniform <- randtoolbox::halton(S * N, K_w, normal = FALSE)

  check_pairs <- list(
    list(i = 1L, s = 0L), list(i = 1L, s = 9L),
    list(i = 2L, s = 0L), list(i = 2L, s = 5L),
    list(i = 3L, s = 0L), list(i = 3L, s = 9L)
  )
  primes_r <- c(2, 3, 5, 7, 11, 13)  # first 6 primes (K_w <= 3 here)

  for (pair in check_pairs) {
    i <- pair$i; s <- pair$s
    n_global  <- (i - 1L) * S + s + 1L   # 1-based global Halton index
    col_idx_R <- n_global                 # same value: 1-based column in S*N output
    for (k in seq_len(K_w)) {
      u_rt <- rt_uniform[n_global, k]
      z_expected <- stats::qnorm(u_rt)
      z_actual   <- full_normal[k, col_idx_R]
      expect_true(
        abs(z_expected - z_actual) < 1e-10,
        label = paste0("slice indexing: i=", i, " s=", s, " k=", k)
      )
    }
  }
})

# ---------------------------------------------------------------------------
# 5. halton_generate_normal output layout documentation test
# ---------------------------------------------------------------------------

test_that("halton_generate_normal has documented output layout K_w x (S*N)", {
  S <- 7L; N <- 5L; K_w <- 3L
  m <- choicer:::halton_generate_normal(S, N, K_w, seed = 1, scramble = 0)

  # Dimensions: K_w rows, S*N columns
  expect_equal(nrow(m), K_w)
  expect_equal(ncol(m), S * N)

  # All values are finite real numbers
  expect_true(all(is.finite(m)))

  # Spot-check column layout: column 1 is individual 1 draw 0 (n=1),
  # column S+1 is individual 2 draw 0 (n=S+1), dimension 1 (base 2).
  col_i1_s0 <- m[1, 1]
  col_i2_s0 <- m[1, S + 1L]

  z_i1 <- stats::qnorm(choicer:::halton_radical_inverse(1,      2))
  z_i2 <- stats::qnorm(choicer:::halton_radical_inverse(S + 1L, 2))

  expect_equal(col_i1_s0, z_i1, tolerance = 1e-10)
  expect_equal(col_i2_s0, z_i2, tolerance = 1e-10)
})

test_that("the full 128-prime Halton table is addressable", {
  # Dimension indices are 0..127 in C++, so K_w = 128 is valid and uses the
  # final prime (719). This guards the former >= off-by-one rejection.
  u <- choicer:::halton_generate_uniform(
    n = 2L, dim = 128L, seed = 0, scramble = 0L
  )
  expect_equal(dim(u), c(2L, 128L))
  expect_equal(u[1L, 128L], 1 / 719, tolerance = 1e-15)
  expect_equal(u[2L, 128L], 2 / 719, tolerance = 1e-15)

  expect_error(
    choicer:::halton_generate_uniform(1L, 129L, seed = 0, scramble = 0L),
    "between 1 and 128"
  )
  expect_error(
    choicer:::halton_generate_normal(1L, 1L, 129L, seed = 0, scramble = 0L),
    "between 1 and 128"
  )

  # Exercise a production generate-mode kernel at K_w = 128, not only the
  # test wrapper. W is zero, so the expected shares are exactly symmetric.
  K_w <- 128L
  shares <- mxl_predict_shares(
    theta = c(0, rep(log(0.1), K_w)),
    X = matrix(0, 2L, 1L), W = matrix(0, 2L, K_w),
    alt_idx = 1:2, M = 2L, weights = 1,
    eta_draws = array(0, dim = c(K_w, 0L, 0L)),
    rc_dist = rep(0L, K_w), rc_correlation = FALSE,
    rc_mean = FALSE, use_asc = FALSE, include_outside_option = FALSE,
    gen_seed = 0L, gen_scramble = 0L, gen_S = 1L
  )
  expect_equal(as.numeric(shares), c(0.5, 0.5), tolerance = 1e-15)
})

# ---------------------------------------------------------------------------
# 6. Block generator: HaltonGen::fill_block() against the per-index loop
# ---------------------------------------------------------------------------

test_that("fill_block reproduces the per-index draws bit for bit", {
  # fill_block() forms a block's uniforms first -- digits by division by a
  # compile-time constant for bases up to 31, base 2 by an exact integer
  # odometer while the block's last index is below 2^53 -- and then applies
  # inv_normal_cdf. Generate-mode draws, hence fits, must not change, so the
  # block must equal inv_normal_cdf(scrambled_halton_uniform(n, k)) exactly.
  # The starts put blocks just below, at and across 2^53 (the odometer's
  # limit), across 2^32 and 2^52 (where the odometer grows a top digit, at
  # 2^52 its last), and up to 64 binary digits; K_w = 12 and 20 reach bases
  # 37 and 71, past the compile-time bases. 32-bit x87 builds keep doubles in
  # 80-bit registers, where the two code paths may round differently.
  skip_if(R.version$arch %in% c("i386", "i686"), "x87 excess precision")
  starts <- function(S) {
    c(1, 1e4, 2.5e8, 2^32 - S %/% 2, 1e10, 2^52 - S %/% 2, 2^53 - S,
      2^53 - S + 1, 2^53 - S %/% 2, 1e18, 2^64 - 2048)
  }
  differs <- function(n0, S, K_w, scramble, seed = 2718) {
    blk <- choicer:::halton_fill_block(n0, S, K_w, seed, scramble)
    ref <- choicer:::halton_reference_block(n0, S, K_w, seed, scramble)
    if (identical(blk, ref, num.eq = FALSE)) character(0) else
      sprintf("n0 = %.17g, S = %d, K_w = %d, scramble = %d", n0, S, K_w,
              scramble)
  }
  for (S in c(1L, 2L, 7L, 100L, 1000L)) {
    bad <- character(0)
    for (K_w in c(1L, 3L, 12L, 20L)) {
      for (scramble in 0:1) {
        for (n0 in starts(S)) bad <- c(bad, differs(n0, S, K_w, scramble))
      }
    }
    expect_identical(bad, character(0), label = sprintf("blocks at S = %d", S))
  }

  # The whole primes table (bases up to 719; near 2^64 every base reaches its
  # top digit), and a block that ends at the last index, 2^64 - 1.
  bad <- character(0)
  for (scramble in 0:1) {
    for (n0 in c(1, 2^53 - 7, 1e18, 2^64 - 2048)) {
      bad <- c(bad, differs(n0, 7L, 128L, scramble))
    }
    bad <- c(bad, differs(2^64 - 2048, 2048L, 3L, scramble))
  }
  expect_identical(bad, character(0), label = "K_w = 128 and the last index")
})

test_that("the permutation table keeps every digit position of a 64-bit index", {
  # Both sides of the block test read the same table, so its extent is checked
  # here: dimension k keeps one permutation per base-b digit of 2^64 - 1, the
  # largest index. That count is 64 in base 2 and ceiling(64 / log2(b)) for
  # an odd prime b, where 64 / log2(b) is never an integer (the margin below
  # keeps the floating-point ceiling unambiguous).
  primes <- function(n) {
    p <- integer(0); k <- 2L
    while (length(p) < n) {
      if (all(k %% p[p * p <= k] != 0L)) p <- c(p, k)
      k <- k + 1L
    }
    p
  }
  b <- primes(128L)
  ratio <- 64 / log2(b[-1L])
  expect_true(all(abs(ratio - round(ratio)) > 1e-6))
  digits <- c(64, ceiling(ratio))
  lay <- choicer:::halton_table_layout(128L)
  expect_identical(lay$digits, as.integer(digits))
  expect_identical(lay$offsets, cumsum(c(0, utils::head(digits * b, -1L))))
  expect_identical(lay$size, sum(digits * b))
})

test_that("the block exports reject invalid blocks", {
  msg <- "n0 must be an integer in \\[1, 2\\^64\\)"
  expect_error(choicer:::halton_fill_block(0, 1L, 1L, 0, 0L), msg)
  expect_error(choicer:::halton_fill_block(1.5, 1L, 1L, 0, 0L), msg)
  expect_error(choicer:::halton_fill_block(NaN, 1L, 1L, 0, 0L), msg)
  expect_error(choicer:::halton_fill_block(2^64, 1L, 1L, 0, 0L), msg)
  expect_error(choicer:::halton_reference_block(2^64 - 2048, 2049L, 1L, 0, 0L),
               "passes 2\\^64 - 1")
  expect_error(choicer:::halton_fill_block(1, 0L, 1L, 0, 0L),
               "S must be positive")
  expect_error(choicer:::halton_reference_block(1, 1L, 129L, 0, 0L),
               "between 1 and 128")
  expect_error(choicer:::halton_fill_block(1, 1L, 1L, -1, 1L),
               "seed must be an integer in \\[0, 2\\^64\\)")
})
