# Tests for utility functions

# --- convertTime tests ---

test_that("convertTime formats seconds correctly", {
  time <- c(user = 0, system = 0, elapsed = 45)
  result <- convertTime(time)
  expect_equal(result, "0h:0m:45s")
})
test_that("convertTime formats minutes correctly", {
  time <- c(user = 0, system = 0, elapsed = 125)  # 2min 5sec
  result <- convertTime(time)
  expect_equal(result, "0h:2m:5s")
})

test_that("convertTime formats hours correctly", {
  time <- c(user = 0, system = 0, elapsed = 3661)  # 1h 1m 1s
  result <- convertTime(time)
  expect_equal(result, "1h:1m:1s")
})

test_that("convertTime handles zero", {
  time <- c(user = 0, system = 0, elapsed = 0)
  result <- convertTime(time)
  expect_equal(result, "0h:0m:0s")
})

test_that("convertTime handles sub-second times", {
  time <- c(user = 0, system = 0, elapsed = 0.5)
  result <- convertTime(time)
  expect_equal(result, "0h:0m:0.5s")
})

# --- vech tests ---

test_that("vech extracts lower triangle correctly for 2x2 matrix", {
  M <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  # M = [1 2]
  #     [3 4]
  # Lower triangle (column-major): 1, 3, 4
  result <- vech_col(M)
  expect_equal(result, c(1, 3, 4))
})

test_that("vech extracts lower triangle correctly for 3x3 matrix", {
  M <- matrix(1:9, nrow = 3, byrow = TRUE)
  # M = [1 2 3]
  #     [4 5 6]
  #     [7 8 9]
  # Lower triangle (column-major): 1, 4, 7, 5, 8, 9
  result <- vech_col(M)
  expect_equal(result, c(1, 4, 7, 5, 8, 9))
})

test_that("vech handles 1x1 matrix", {
  M <- matrix(5, nrow = 1, ncol = 1)
  result <- vech_col(M)
  expect_equal(result, 5)
})

test_that("vech returns correct length", {
  for (n in 1:5) {
    M <- matrix(runif(n * n), nrow = n)
    result <- vech_col(M)
    expected_length <- n * (n + 1) / 2
    expect_length(result, expected_length)
  }
})

# --- get_halton_normals tests ---

test_that("get_halton_normals returns correct dimensions", {
  S <- 50
  N <- 30
  K_w <- 3

  eta <- get_halton_normals(S, N, K_w)

  expect_equal(dim(eta), c(K_w, S, N))
})

test_that("get_halton_normals produces finite values", {
  S <- 100
  N <- 50
  K_w <- 2

  eta <- get_halton_normals(S, N, K_w)

  expect_true(all(is.finite(eta)))
})

test_that("get_halton_normals produces approximately standard normal draws", {
  S <- 500
  N <- 100
  K_w <- 2

  eta <- get_halton_normals(S, N, K_w)

  # Flatten and check moments
  all_draws <- as.vector(eta)

  # Mean should be close to 0
  expect_equal(mean(all_draws), 0, tolerance = 0.1)

  # SD should be close to 1
  expect_equal(sd(all_draws), 1, tolerance = 0.1)
})

test_that("get_halton_normals handles K_w = 1", {
  S <- 50
  N <- 20
  K_w <- 1

  eta <- get_halton_normals(S, N, K_w)

  expect_equal(dim(eta), c(1, S, N))
})

test_that("get_halton_normals is deterministic", {
  S <- 30
  N <- 20
  K_w <- 2

  eta1 <- get_halton_normals(S, N, K_w)
  eta2 <- get_halton_normals(S, N, K_w)

  expect_equal(eta1, eta2)
})

test_that("get_halton_normals gives unit i the points (i - 1) S + 1, ..., i S", {
  # Reference: one halton() call over the whole sequence, cut into consecutive
  # blocks of S points, one per unit.
  full_cube <- function(S, N, K_w) {
    h <- matrix(randtoolbox::halton(S * N, K_w, normal = TRUE), ncol = K_w)
    array(t(h), dim = c(K_w, S, N))
  }
  for (p in list(c(1, 1, 1), c(7, 5, 1), c(50, 9, 2), c(3, 11, 5), c(1, 4, 3))) {
    S <- p[1]; N <- p[2]; K_w <- p[3]
    ref <- full_cube(S, N, K_w)
    expect_identical(get_halton_normals(S, N, K_w), ref)
    # Blocks of 1, 2, 3, N - 1, N and N + 1 units, including a partial last block.
    for (units in unique(pmax(1, c(1, 2, 3, N - 1, N, N + 1)))) {
      expect_identical(.halton_cube(S, N, K_w, block = units * S * K_w), ref)
    }
  }
})

test_that("get_halton_normals validates sizes and the sequence index limit", {
  expect_error(get_halton_normals(0, 10, 2), "`S` must be a single positive whole number")
  expect_error(get_halton_normals(10, 2.5, 2), "`N` must be a single positive whole number")
  expect_error(get_halton_normals(10, 10, NA), "`K_w` must be a single positive whole number")
  expect_error(get_halton_normals(c(10, 20), 10, 2), "`S` must be")
  # None of these allocates: the guards run first.
  expect_error(get_halton_normals(100, 3e7, 1), "more than 2\\^31 - 1")
  expect_error(get_halton_normals(2, 2^30, 1), "more than 2\\^31 - 1")       # S * N = 2^31
  expect_error(get_halton_normals(50000L, 50000L, 1L), "more than 2\\^31 - 1") # no integer NA
  expect_error(get_halton_normals(100, 1.5e7, 3), "more than 2\\^32 - 1")    # K_w * S * N
  expect_error(get_halton_normals(2^16, 2^14, 4), "more than 2\\^32 - 1")    # exactly 2^32
})

# --- check_collinearity / remove_nullspace_cols tests ---

test_that("check_collinearity returns list with correct elements", {
  X <- matrix(rnorm(20), nrow = 5, ncol = 4)
  colnames(X) <- c("a", "b", "c", "d")

  result <- check_collinearity(X)

  expect_true("mat" %in% names(result))
  expect_true("dropped" %in% names(result))
})

test_that("check_collinearity identifies linearly dependent columns", {
  # Create columns where c = a + b (perfectly collinear)
  set.seed(123)
  a <- rnorm(10)
  b <- rnorm(10)
  X <- cbind(
    a = a,
    b = b,
    c = a + b  # Linear combination
  )

  result <- check_collinearity(X)

  # Should drop one column
  expect_true(ncol(result$mat) < 3 || length(result$dropped) > 0)
})

test_that("check_collinearity handles all independent columns", {
  X <- diag(5)
  colnames(X) <- letters[1:5]

  result <- check_collinearity(X)

  expect_equal(ncol(result$mat), 5)
  expect_equal(length(result$dropped), 0)
})

test_that("check_collinearity handles single-column matrix", {
  X <- matrix(1:5, ncol = 1)
  colnames(X) <- "a"

  result <- check_collinearity(X)

  expect_equal(ncol(result$mat), 1)
})

test_that(".qr_rank_pivot makes qr()'s rank and drop decisions by row chunks", {
  set.seed(7)
  n <- 3000
  x1 <- rnorm(n); x2 <- rnorm(n) * 1e4; z <- rnorm(n)
  alt <- sample.int(4, n, replace = TRUE)
  D <- sapply(1:4, function(j) as.numeric(alt == j))
  rare <- numeric(n); rare[c(3, 77)] <- 1
  ints <- cbind(a = sample(0:5, n, replace = TRUE),
                b = sample(0:3, n, replace = TRUE))
  designs <- list(
    lincomb = cbind(x1, x2, x3 = 0.7 * x1 - 3e-5 * x2),
    near_dependent = cbind(x1, x2, x3 = x1 + 1e-10 * z),
    near_independent = cbind(x1, x2, x3 = x1 + 1e-4 * z),
    dummy_trap = cbind(int = 1, D, x1),
    duplicated = cbind(a = x1, b = z, c = x1),
    zero = cbind(a = x1, z = 0, b = z),
    two_constants = cbind(a = 1, b = 2, c = z),
    rare_dummy = cbind(d = rare, x = x1, dx = 2 * rare),
    integer = cbind(ints, c = ints[, "a"] + 2L * ints[, "b"]),
    independent = cbind(x1, x2, z)
  )
  expect_type(designs$integer, "integer")
  dropped <- function(q, p) sort(q$pivot[seq_len(p)[-seq_len(q$rank)]])
  for (nm in names(designs)) {
    m <- designs[[nm]]
    ref <- qr(m, tol = 1e-7)
    for (rows in c(1, 2, ncol(m), 7, 64, n)) {
      q <- .qr_rank_pivot(m, rows = rows)
      expect_identical(q$rank, ref$rank, info = paste(nm, rows))
      expect_identical(dropped(q, ncol(m)), dropped(ref, ncol(m)),
                       info = paste(nm, rows))
    }
  }
})

test_that("remove_nullspace_cols matches qr() on both sides of the switch", {
  set.seed(8)
  n <- 300000                                      # a default chunk + remainder
  m <- cbind(a = rnorm(n), b = rnorm(n), c = rnorm(n))
  m <- cbind(m, d = m[, "a"] - 2 * m[, "c"])      # n * p >= 1e6: by chunks
  expect_identical(remove_nullspace_cols(m), m[, c("a", "b", "c")])
  small <- m[seq_len(240000), ]                    # n * p < 1e6: qr()
  expect_identical(remove_nullspace_cols(small), small[, c("a", "b", "c")])
  m[n - 5, "b"] <- NA                              # in the second chunk
  expect_error(remove_nullspace_cols(m), "NA/NaN/Inf in foreign function call")
  m[n - 5, "b"] <- Inf
  expect_error(remove_nullspace_cols(m), "NA/NaN/Inf in foreign function call")
})

test_that(".n_distinct_by and .first_by match the per-group expressions they replace", {
  dt <- data.table::data.table(id = c(3L, 1L, 3L, 1L, 2L, 1L, 2L),
                               p = c(5L, 6L, 5L, 7L, NA, 6L, 8L),
                               w = c(1.5, 2, 1.5, 2, 3, 2, NA),
                               s = c("a", NA, "a", "b", "c", NA, "c"))
  ids <- unique(dt$id)
  for (col in c("p", "w", "s", "id")) {
    expect_identical(.n_distinct_by(dt, col, "id"),
                     dt[, data.table::uniqueN(get(col)), by = "id"][["V1"]])
    expect_identical(.first_by(dt, col, "id", ids),
                     dt[, get(col)[1L], by = "id"][["V1"]][match(ids, unique(dt$id))])
  }
})

test_that("per-situation columns work whatever the id and value columns are called", {
  # With an id column named V1 the old per-group results were read back as
  # the ids; a comma in a name broke `by = c(id, col)`.
  dt <- data.table::data.table(V1 = c(1L, 1L, 2L, 2L, 3L, 3L),
                               w = c(0.5, 0.5, 2, 2, 4, 4))
  data.table::setnames(dt, "w", "w,adj")
  expect_identical(.n_distinct_by(dt, "w,adj", "V1"), c(1L, 1L, 1L))
  expect_identical(.first_by(dt, "w,adj", "V1", 3:1), c(4, 2, 0.5))
  expect_identical(.collapse_situation_col(dt, "w,adj", "V1", 1:3), c(0.5, 2, 4))
})

# --- OpenMP thread control tests ---

test_that("get_num_threads returns valid output", {
  # This should not error and return something
  result <- get_num_threads()

  # Result might be printed; just check it doesn't error
  expect_true(TRUE)
})

test_that("thread_info returns structured OpenMP state", {
  info <- thread_info()

  expect_type(info, "list")
  expect_named(info, c(
    "openmp_enabled", "_OPENMP", "omp_get_num_threads",
    "omp_get_max_threads", "omp_get_num_procs", "omp_get_thread_limit",
    "OMP_THREAD_LIMIT", "OMP_NUM_THREADS"
  ))
  expect_type(info$openmp_enabled, "logical")

  if (isTRUE(info$openmp_enabled)) {
    expect_type(info$`_OPENMP`, "integer")
    expect_type(info$omp_get_num_threads, "integer")
    expect_type(info$omp_get_max_threads, "integer")
    expect_type(info$omp_get_num_procs, "integer")
    expect_type(info$omp_get_thread_limit, "integer")
    expect_gte(info$omp_get_num_threads, 1L)
    expect_gte(info$omp_get_max_threads, 1L)
    expect_gte(info$omp_get_num_procs, 1L)
  } else {
    expect_true(is.na(info$`_OPENMP`))
    expect_true(is.na(info$omp_get_num_threads))
  }
})

test_that("set_num_threads accepts valid input", {
  # Should not error
  expect_no_error(set_num_threads(1))
  expect_no_error(set_num_threads(2))
})

test_that("set_num_threads updates reported OpenMP max threads", {
  info_before <- thread_info()
  on.exit({
    if (isTRUE(info_before$openmp_enabled)) {
      set_num_threads(info_before$omp_get_max_threads)
    }
  }, add = TRUE)

  expect_error(set_num_threads(0), "`n_threads` must be a positive integer", fixed = TRUE)
  expect_no_error(set_num_threads(1))

  info_after <- thread_info()
  if (isTRUE(info_after$openmp_enabled)) {
    expect_equal(info_after$omp_get_max_threads, 1L)
    expect_equal(info_after$omp_get_num_threads, 1L)
  }
})

test_that("choicer is built with Armadillo's 64-bit word", {
  # src/Makevars defines ARMA_64BIT_WORD: the kernels view R's matrices and
  # cubes without copying them, and under a 32-bit word the element count of
  # one of more than 2^32 - 1 values wraps. Armadillo keeps the 32-bit word
  # where pointers are 32-bit, so the expectation follows the platform.
  expect_identical(choicer:::test_arma_word_bytes(), .Machine$sizeof.pointer)
  if (.Machine$sizeof.pointer == 8L) {
    expect_identical(choicer:::test_arma_view_n_elem(), 6442450941)
  }
})
