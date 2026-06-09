# Tests for counterfactual prediction: predict(..., newdata = )
# - data.frame path (prepare_newdata)
# - list path ("modified design matrix" for policy simulation)
# - NL nest_idx resolution (top-level field with $data fallback)

# =============================================================================
# Helpers: fit small models and keep the raw data alongside
# =============================================================================

make_mnl <- function(keep_data = TRUE, scale_vars = "none") {
  dt <- create_small_mnl_data()
  fit <- run_mnlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    control = list(maxeval = 50L),
    keep_data = keep_data,
    scale_vars = scale_vars
  )
  list(fit = fit, dt = dt)
}

make_mxl <- function(keep_data = TRUE) {
  dt <- create_small_mxl_data()
  fit <- run_mxlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = "x1",
    random_var_cols = c("w1", "w2"),
    S = 10L,
    control = list(maxeval = 50L),
    keep_data = keep_data
  )
  list(fit = fit, dt = dt)
}

make_nl <- function(keep_data = TRUE) {
  dt <- create_small_nl_data()
  fit <- run_nestlogit(
    data = dt,
    id_col = "id",
    alt_col = "alt",
    choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    nest_col = "nest",
    control = list(maxeval = 50L),
    keep_data = keep_data
  )
  list(fit = fit, dt = dt)
}

# MNL dataset with an explicit outside option (alt = 0, zero covariates);
# some ids choose the outside option (choice = 0 on all inside rows).
create_mnl_outside_data <- function(seed = 7, N = 40, J_inside = 3) {
  set.seed(seed)
  dt <- data.table(
    id = rep(1:N, each = J_inside + 1),
    alt = rep(0:J_inside, N),
    x1 = rnorm(N * (J_inside + 1)),
    x2 = runif(N * (J_inside + 1), -1, 1)
  )
  dt[alt == 0, c("x1", "x2") := 0]
  dt[, choice := 0L]
  dt[, choice := {
    pick <- sample.int(J_inside + 1, 1) - 1L
    as.integer(alt == pick)
  }, by = id]
  dt[]
}

# Per-id probability sums for a stacked choice_prob vector
per_id_sums <- function(choice_prob, M) {
  ends <- cumsum(M)
  starts <- c(1L, head(ends, -1) + 1L)
  vapply(seq_along(M),
         function(i) sum(choice_prob[starts[i]:ends[i]]),
         numeric(1))
}

# =============================================================================
# 1. Round-trip identity: newdata = original data reproduces stored-data path
# =============================================================================

test_that("MNL round-trip: newdata = original data matches predict(fit)", {
  m <- make_mnl()

  expect_equal(
    predict(m$fit, type = "probabilities", newdata = m$dt),
    predict(m$fit, type = "probabilities"),
    tolerance = 1e-12
  )
  expect_equal(
    predict(m$fit, type = "shares", newdata = m$dt),
    predict(m$fit, type = "shares"),
    tolerance = 1e-12
  )
})

test_that("MXL round-trip: newdata = original data matches predict(fit)", {
  m <- make_mxl()

  expect_equal(
    predict(m$fit, type = "probabilities", newdata = m$dt),
    predict(m$fit, type = "probabilities"),
    tolerance = 1e-12
  )
  expect_equal(
    predict(m$fit, type = "shares", newdata = m$dt),
    predict(m$fit, type = "shares"),
    tolerance = 1e-12
  )
})

test_that("NL round-trip: newdata = original data matches predict(fit)", {
  m <- make_nl()

  expect_equal(
    predict(m$fit, type = "probabilities", newdata = m$dt),
    predict(m$fit, type = "probabilities"),
    tolerance = 1e-12
  )
  expect_equal(
    predict(m$fit, type = "shares", newdata = m$dt),
    predict(m$fit, type = "shares"),
    tolerance = 1e-12
  )
})

# =============================================================================
# 2. Row order and irrelevant columns are immaterial
# =============================================================================

test_that("shuffled rows and extra columns in newdata do not change predictions", {
  m <- make_mnl()

  set.seed(99)
  dt_shuffled <- m$dt[sample(.N)]
  dt_shuffled[, junk := rnorm(.N)]
  dt_shuffled[, junk_chr := "a"]

  expect_equal(
    predict(m$fit, type = "probabilities", newdata = dt_shuffled),
    predict(m$fit, type = "probabilities"),
    tolerance = 1e-12
  )
  expect_equal(
    predict(m$fit, type = "shares", newdata = dt_shuffled),
    predict(m$fit, type = "shares"),
    tolerance = 1e-12
  )
})

# =============================================================================
# 3. Natural-scale invariance: scale_vars fit predicts identically on newdata
# =============================================================================

test_that("MNL fit with scale_vars = 'sd' round-trips on newdata", {
  m <- make_mnl(scale_vars = "sd")

  expect_equal(
    predict(m$fit, type = "probabilities", newdata = m$dt),
    predict(m$fit, type = "probabilities"),
    tolerance = 1e-12
  )
  expect_equal(
    predict(m$fit, type = "shares", newdata = m$dt),
    predict(m$fit, type = "shares"),
    tolerance = 1e-12
  )
})

# =============================================================================
# 4. keep_data = FALSE: newdata works, stored-data path keeps erroring
# =============================================================================

test_that("keep_data = FALSE fits predict with newdata but error without", {
  m_full <- make_mnl(keep_data = TRUE)
  m_slim <- make_mnl(keep_data = FALSE)

  expect_error(predict(m_slim$fit), "keep_data")

  expect_equal(
    predict(m_slim$fit, type = "probabilities", newdata = m_slim$dt),
    predict(m_full$fit, type = "probabilities"),
    tolerance = 1e-12
  )
})

test_that("keep_data = FALSE MXL and NL fits predict with newdata", {
  m_mxl <- make_mxl(keep_data = FALSE)
  expect_error(predict(m_mxl$fit), "keep_data")
  p_mxl <- predict(m_mxl$fit, type = "probabilities", newdata = m_mxl$dt)
  expect_length(p_mxl$choice_prob, nrow(m_mxl$dt))

  m_nl <- make_nl(keep_data = FALSE)
  expect_error(predict(m_nl$fit), "keep_data")
  p_nl <- predict(m_nl$fit, type = "probabilities", newdata = m_nl$dt)
  expect_length(p_nl$choice_prob, nrow(m_nl$dt))
})

# =============================================================================
# 5. Counterfactual direction and probability normalization
# =============================================================================

test_that("increasing x1 for alternative 2 moves its share with sign(coef)", {
  m <- make_mnl()
  base_shares <- predict(m$fit, type = "shares")

  dt_cf <- copy(m$dt)
  dt_cf[alt == 2, x1 := x1 + 1]

  cf_shares <- predict(m$fit, type = "shares", newdata = dt_cf)
  b_x1 <- coef(m$fit)[["x1"]]

  expect_equal(
    sign(cf_shares[[2]] - base_shares[[2]]),
    sign(b_x1)
  )

  # Probabilities still sum to 1 per individual
  preds <- predict(m$fit, type = "probabilities", newdata = dt_cf)
  M <- m$fit$data$M
  expect_equal(per_id_sums(preds$choice_prob, M), rep(1, length(M)),
               tolerance = TOL_PROB)
})

# =============================================================================
# 6. Input validation errors (data.frame path)
# =============================================================================

test_that("newdata with unseen alternative label errors", {
  m <- make_mnl()
  dt_bad <- copy(m$dt)
  dt_bad[1, alt := 99L]

  expect_error(
    predict(m$fit, newdata = dt_bad),
    "not seen at fit time"
  )
})

test_that("newdata missing a covariate column errors", {
  m <- make_mnl()
  dt_bad <- copy(m$dt)
  dt_bad[, x2 := NULL]

  expect_error(
    predict(m$fit, newdata = dt_bad),
    "Missing columns in newdata: x2"
  )
})

test_that("newdata with NA in a covariate errors", {
  m <- make_mnl()
  dt_bad <- copy(m$dt)
  dt_bad[3, x1 := NA_real_]

  expect_error(
    predict(m$fit, newdata = dt_bad),
    "missing values"
  )
})

test_that("newdata with duplicated (id, alt) pairs errors", {
  m <- make_mnl()
  dt_bad <- rbind(m$dt, m$dt[1])

  expect_error(
    predict(m$fit, newdata = dt_bad),
    "duplicated"
  )
})

test_that("non-numeric covariate in newdata errors", {
  m <- make_mnl()
  dt_bad <- copy(m$dt)
  dt_bad[, x1 := as.character(x1)]

  expect_error(
    predict(m$fit, newdata = dt_bad),
    "numeric"
  )
})

test_that("weights of wrong length error; weights ignored without newdata", {
  m <- make_mnl()

  expect_error(
    predict(m$fit, type = "shares", newdata = m$dt, weights = c(1, 2)),
    "'weights' must be a numeric vector"
  )

  # Ignored (no error) when newdata is NULL
  expect_equal(
    predict(m$fit, type = "shares", weights = c(1, 2)),
    predict(m$fit, type = "shares")
  )
})

test_that("prediction weights are used for share aggregation", {
  m <- make_mnl()
  N <- length(m$fit$data$M)

  # All weight on the first individual -> shares equal its probabilities
  w <- c(1, rep(0, N - 1))
  shares_w <- predict(m$fit, type = "shares", newdata = m$dt, weights = w)
  preds <- predict(m$fit, type = "probabilities", newdata = m$dt)
  M1 <- m$fit$data$M[1]

  expect_equal(as.numeric(shares_w)[seq_len(M1)],
               preds$choice_prob[seq_len(M1)],
               tolerance = 1e-10)
})

# =============================================================================
# 7. MXL: newdata with a subset of ids (different N)
# =============================================================================

test_that("MXL predicts on a subset of ids with probabilities summing to 1", {
  m <- make_mxl()
  nd <- m$dt[id <= 10]

  preds <- predict(m$fit, type = "probabilities", newdata = nd)
  expect_length(preds$choice_prob, nrow(nd))

  M <- nd[, .N, by = id][["N"]]
  expect_equal(per_id_sums(preds$choice_prob, M), rep(1, length(M)),
               tolerance = 1e-6)

  shares <- predict(m$fit, type = "shares", newdata = nd)
  expect_equal(sum(shares), 1, tolerance = 1e-6)
})

# =============================================================================
# 8. NL nest_idx resolution: top-level field, $data fallback, refit error
# =============================================================================

test_that("NL newdata prediction uses top-level nest_idx with $data fallback", {
  m <- make_nl()

  # Fresh fits carry nest_idx top-level
  expect_false(is.null(m$fit$nest_idx))
  expect_identical(m$fit$nest_idx, m$fit$data$nest_idx)

  p_new <- predict(m$fit, type = "probabilities", newdata = m$dt)

  # Simulate an old-style fit: no top-level nest_idx, fallback to $data
  fit_old <- m$fit
  fit_old$nest_idx <- NULL
  p_old <- predict(fit_old, type = "probabilities", newdata = m$dt)
  expect_equal(p_old, p_new, tolerance = 1e-12)

  # Both missing -> informative error
  fit_none <- fit_old
  fit_none$data <- NULL
  expect_error(
    predict(fit_none, newdata = m$dt),
    "Refit to enable newdata prediction"
  )
})

# =============================================================================
# 9. Outside option: inside probabilities + outside probability sum to 1
# =============================================================================

test_that("MNL with outside option: inside + outside probabilities sum to 1", {
  dt <- create_mnl_outside_data()
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"),
    outside_opt_label = 0L,
    include_outside_option = TRUE,
    control = list(maxeval = 50L)
  )

  # Round trip: outside rows in newdata are dropped internally
  preds <- predict(fit, type = "probabilities", newdata = dt)
  expect_equal(preds, predict(fit, type = "probabilities"),
               tolerance = 1e-12)

  # Per id: P_outside = 1 / (1 + sum(exp(V_inside))) and the full
  # choice set (inside + outside) must sum to 1.
  M <- fit$data$M
  ends <- cumsum(M)
  starts <- c(1L, head(ends, -1) + 1L)
  full_sums <- vapply(seq_along(M), function(i) {
    rows <- starts[i]:ends[i]
    p_out <- 1 / (1 + sum(exp(preds$utility[rows])))
    sum(preds$choice_prob[rows]) + p_out
  }, numeric(1))
  expect_equal(full_sums, rep(1, length(M)), tolerance = TOL_PROB)

  # Shares include the outside option (first element) and sum to 1
  shares <- predict(fit, type = "shares", newdata = dt)
  J_inside <- max(fit$alt_mapping$alt_int)
  expect_length(shares, J_inside + 1)
  expect_equal(sum(shares), 1, tolerance = 1e-8)
})

# =============================================================================
# 10. List path ("modified design matrix")
# =============================================================================

test_that("list path matches the data.frame path on the same inputs", {
  m <- make_mnl()
  d <- m$fit$data
  nd_list <- list(X = d$X, alt_idx = d$alt_idx, M = d$M)

  expect_equal(
    predict(m$fit, type = "probabilities", newdata = nd_list),
    predict(m$fit, type = "probabilities", newdata = m$dt),
    tolerance = 1e-12
  )
  expect_equal(
    predict(m$fit, type = "shares", newdata = nd_list),
    predict(m$fit, type = "shares", newdata = m$dt),
    tolerance = 1e-12
  )
})

test_that("MXL list path works and requires W", {
  m <- make_mxl()
  d <- m$fit$data

  expect_equal(
    predict(m$fit, type = "probabilities",
            newdata = list(X = d$X, W = d$W, alt_idx = d$alt_idx, M = d$M)),
    predict(m$fit, type = "probabilities"),
    tolerance = 1e-12
  )

  expect_error(
    predict(m$fit, newdata = list(X = d$X, alt_idx = d$alt_idx, M = d$M)),
    "must include 'W'"
  )
})

test_that("list path validates dimensions and alternative codes", {
  m <- make_mnl()
  d <- m$fit$data

  # Wrong number of X columns
  expect_error(
    predict(m$fit, newdata = list(X = d$X[, 1, drop = FALSE],
                                  alt_idx = d$alt_idx, M = d$M)),
    "expects 2"
  )

  # Mismatched column names
  X_bad <- d$X
  colnames(X_bad) <- c("x2", "x1")
  expect_error(
    predict(m$fit, newdata = list(X = X_bad, alt_idx = d$alt_idx, M = d$M)),
    "column names must match"
  )

  # sum(M) != nrow(X)
  expect_error(
    predict(m$fit, newdata = list(X = d$X, alt_idx = d$alt_idx,
                                  M = d$M + 1L)),
    "sum\\(newdata\\$M\\)"
  )

  # alt_idx out of range
  alt_bad <- d$alt_idx
  alt_bad[1] <- 99L
  expect_error(
    predict(m$fit, newdata = list(X = d$X, alt_idx = alt_bad, M = d$M)),
    "alt_idx"
  )

  # Missing required element
  expect_error(
    predict(m$fit, newdata = list(X = d$X, alt_idx = d$alt_idx)),
    "missing element"
  )
})

# =============================================================================
# 10. Review fixes: weight realignment, finiteness, NA alt labels, list checks
# =============================================================================

test_that("prediction weights follow first-appearance id order in newdata", {
  m <- make_mnl()
  N <- length(m$fit$data$M)

  # Reverse the id blocks: id N appears first in newdata
  nd_rev <- m$dt[order(-id, alt)]
  # All weight on the first id to appear in newdata (id N)
  w <- c(1, rep(0, N - 1))
  shares_rev <- predict(m$fit, type = "shares", newdata = nd_rev, weights = w)

  # Same weighting expressed against sorted newdata: all weight on id N
  w_sorted <- c(rep(0, N - 1), 1)
  shares_sorted <- predict(m$fit, type = "shares", newdata = m$dt,
                           weights = w_sorted)
  expect_equal(shares_rev, shares_sorted, tolerance = 1e-12)
})

test_that("newdata with non-finite covariates errors (data.frame path)", {
  m <- make_mnl()
  nd <- data.table::copy(m$dt)
  nd[1, x1 := Inf]
  expect_error(predict(m$fit, newdata = nd), "finite")
})

test_that("NA alternative labels error even with an outside option", {
  dt <- create_mnl_outside_data()
  fit <- run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), include_outside_option = TRUE,
    outside_opt_label = 0L, control = list(maxeval = 50L)
  )
  nd <- data.table::copy(dt)
  nd[2, alt := NA]
  expect_error(predict(fit, newdata = nd), "missing values")
})

test_that("list path rejects non-integer M and repeated alt codes in a block", {
  m <- make_mnl()
  d <- m$fit$data

  M_frac <- as.numeric(d$M)
  M_frac[1] <- M_frac[1] + 0.7
  M_frac[2] <- M_frac[2] - 0.7
  expect_error(
    predict(m$fit, newdata = list(X = d$X, alt_idx = d$alt_idx, M = M_frac)),
    "positive integers"
  )

  alt_dup <- d$alt_idx
  alt_dup[2] <- alt_dup[1]
  expect_error(
    predict(m$fit, newdata = list(X = d$X, alt_idx = alt_dup, M = d$M)),
    "repeat an alternative"
  )
})
