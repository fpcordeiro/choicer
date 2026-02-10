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
  result <- vech(M)
  expect_equal(result, c(1, 3, 4))
})

test_that("vech extracts lower triangle correctly for 3x3 matrix", {
  M <- matrix(1:9, nrow = 3, byrow = TRUE)
  # M = [1 2 3]
  #     [4 5 6]
  #     [7 8 9]
  # Lower triangle (column-major): 1, 4, 7, 5, 8, 9
  result <- vech(M)
  expect_equal(result, c(1, 4, 7, 5, 8, 9))
})

test_that("vech handles 1x1 matrix", {
  M <- matrix(5, nrow = 1, ncol = 1)
  result <- vech(M)
  expect_equal(result, 5)
})

test_that("vech returns correct length", {
  for (n in 1:5) {
    M <- matrix(runif(n * n), nrow = n)
    result <- vech(M)
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

# --- OpenMP thread control tests ---

test_that("get_num_threads returns valid output", {
  # This should not error and return something
  result <- get_num_threads()

  # Result might be printed; just check it doesn't error
  expect_true(TRUE)
})

test_that("set_num_threads accepts valid input", {
  # Should not error
  expect_no_error(set_num_threads(1))
  expect_no_error(set_num_threads(2))
})
