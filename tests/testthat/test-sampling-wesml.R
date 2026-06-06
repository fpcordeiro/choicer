# Tests for choice-based sampling + WESML weighting + robust sandwich SEs.
# Fast tests (weight math, sampling quotas, validation) run everywhere; tests
# that fit a model are gated with skip_on_cran().

# --- helpers -----------------------------------------------------------------

# Deterministic 6-id, 2-alt data set: ids 1,2,3,6 choose alt 1; ids 4,5 choose 2.
make_tiny <- function() {
  data.table(
    id     = rep(1:6, each = 2),
    alt    = rep(c(1L, 2L), 6),
    x1     = rnorm(12),
    w1     = rnorm(12),
    choice = c(1, 0,  1, 0,  1, 0,  0, 1,  0, 1,  1, 0)
  )
}

# Population with unequal choice strata (alt 3 rare), then a choice-based sample.
make_cb_sample <- function(seed = 1, N = 900L, J = 3L, n_per = 70L) {
  set.seed(seed)
  pop <- data.table(id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N))
  pop[, x1 := rnorm(.N)]
  pop[, w1 := rnorm(.N)]
  pr <- c(0.5, 0.35, 0.15)[seq_len(J)]
  pop[, choice := {
    k <- sample(J, 1L, prob = pr)
    as.integer(seq_len(.N) == k)
  }, by = id]
  sample_by_choice(pop, "id", "alt", "choice", n_per_alt = n_per, seed = seed + 1L)
}

fit_sandwich <- function(dat, S = 30L, se_method = "sandwich", weights_col = ".wesml_weight") {
  suppressWarnings(suppressMessages(
    run_mxlogit(dat, "id", "alt", "choice", "x1", "w1",
                S = S, weights_col = weights_col, se_method = se_method)
  ))
}

# --- 1. weight correctness ---------------------------------------------------

test_that("wesml_weights returns Q(j)/H(j) and supports id-keyed + attach", {
  set.seed(1)
  dt <- make_tiny()
  Q <- c("1" = 0.3, "2" = 0.7)             # already sums to 1
  w <- wesml_weights(dt, "id", "alt", "choice", Q = Q, normalize = FALSE)

  expect_s3_class(w, "data.table")
  expect_equal(nrow(w), 6L)
  expect_setequal(names(w), c("id", ".wesml_weight"))

  # H(1) = 4/6, H(2) = 2/6  ->  w1 = 0.3/(4/6) = 0.45 ; w2 = 0.7/(2/6) = 2.1
  expect_equal(w[order(id)][[".wesml_weight"]],
               c(0.45, 0.45, 0.45, 2.1, 2.1, 0.45), tolerance = 1e-10)
  expect_equal(unname(attr(w, "H")[c("1", "2")]), c(4/6, 2/6), tolerance = 1e-10)

  # attach = TRUE -> row-level column, constant within id
  wa <- wesml_weights(dt, "id", "alt", "choice", Q = Q, attach = TRUE)
  expect_equal(nrow(wa), nrow(dt))
  per_id_u <- wa[, data.table::uniqueN(.wesml_weight), by = id][["V1"]]
  expect_true(all(per_id_u == 1L))
})

# --- 2. weights_col alignment survives reorder / id drop ---------------------

test_that("weights_col collapses by id and survives row shuffle + id drop", {
  set.seed(2)
  N <- 40L; J <- 3L
  dt <- data.table(id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N))
  dt[, x1 := rnorm(.N)]
  dt[, w1 := rnorm(.N)]
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  dt[, wt := as.numeric(id)]               # ground-truth weight == id

  shuffled <- dt[sample(.N)][id != 7L]      # reorder rows, drop id 7 entirely
  prepared <- suppressWarnings(
    prepare_mxl_data(shuffled, "id", "alt", "choice", "x1", "w1",
                     weights_col = "wt")
  )
  expect_equal(prepared$weights, as.numeric(setdiff(seq_len(N), 7L)))
  expect_equal(prepared$N, N - 1L)
})

# --- 3. constant-within-id validation ----------------------------------------

test_that("weights_col that varies within an id errors", {
  set.seed(3)
  dt <- make_tiny()
  dt[, wt := as.numeric(id)]
  dt[1L, wt := 999]                         # break constancy for id 1
  expect_error(
    prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1", weights_col = "wt"),
    "constant within"
  )
  # weights + weights_col are mutually exclusive
  expect_error(
    prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1",
                     weights = rep(1, 6), weights_col = "wt"),
    "only one of"
  )
})

# --- 4. stratum-key coercion (numeric alt vs character Q names) --------------

test_that("character Q names match numeric alternatives; mismatches error", {
  set.seed(4)
  dt <- make_tiny()
  expect_silent(
    wesml_weights(dt, "id", "alt", "choice", Q = c("1" = 0.4, "2" = 0.6))
  )
  expect_error(                              # missing a realized stratum
    wesml_weights(dt, "id", "alt", "choice", Q = c("1" = 1.0)),
    "match the chosen strata"
  )
  expect_error(                              # unexpected extra stratum
    wesml_weights(dt, "id", "alt", "choice",
                  Q = c("1" = 0.3, "2" = 0.5, "3" = 0.2)),
    "match the chosen strata"
  )
})

# --- 7. sample_by_choice quotas / integrity / errors -------------------------

test_that("sample_by_choice honors quotas, keeps situations intact, errors on over-request", {
  set.seed(7)
  N <- 500L; J <- 3L
  pop <- data.table(id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N))
  pop[, x1 := rnorm(.N)]
  pop[, w1 := rnorm(.N)]
  pop[, choice := {
    k <- sample(J, 1L, prob = c(0.6, 0.3, 0.1))
    as.integer(seq_len(.N) == k)
  }, by = id]

  s <- sample_by_choice(pop, "id", "alt", "choice", n_per_alt = 40L, seed = 11L)
  # exactly 40 choosers per stratum
  cnt <- s[choice == 1, .N, by = alt][order(alt)][["N"]]
  expect_equal(cnt, c(40L, 40L, 40L))
  # group integrity: every sampled id retains all J rows
  expect_true(all(s[, .N, by = id][["N"]] == J))
  # provenance attribute
  cs <- attr(s, "choice_sampling")
  expect_equal(cs$scheme, "wesml")
  expect_equal(cs$meat, "robust")

  # attached weights reproduce wesml_weights() on the subsample with attached Q
  w2 <- wesml_weights(s, "id", "alt", "choice", Q = attr(s, "Q"))
  setkey(w2, id)
  per_id <- s[s[, .I[1L], by = id][["V1"]]][order(id)]
  expect_equal(per_id[[".wesml_weight"]], w2[order(id)][[".wesml_weight"]],
               tolerance = 1e-10)

  # over-request without replacement -> error
  expect_error(
    sample_by_choice(pop, "id", "alt", "choice", n_per_alt = 99999L, seed = 1L),
    "without replacement"
  )
  # exactly one of n_per_alt / frac_per_alt
  expect_error(
    sample_by_choice(pop, "id", "alt", "choice", n_per_alt = 10L, frac_per_alt = 0.1),
    "exactly one"
  )
})

# --- 5/6/8. sandwich estimation: uniform-equivalence, invariance, parity -----

test_that("sandwich SEs: uniform-weight equivalence, scale-invariance, wesml_vcov parity", {
  skip_on_cran()
  s <- make_cb_sample(seed = 21)

  fit <- fit_sandwich(s, S = 30L, se_method = "sandwich")
  expect_false(is.null(fit$vcov))
  expect_true(all(is.finite(fit$se)))

  # (8) standalone wesml_vcov() matches the eagerly stored sandwich vcov
  expect_equal(unname(wesml_vcov(fit, "vcov")), unname(fit$vcov), tolerance = 1e-7)

  # (6) scale-invariance: rescale stored weights by 10 -> sandwich vcov unchanged;
  #     inverse-Hessian variance scales by 1/10 and SE by 1/sqrt(10).
  fit10 <- fit
  fit10$data$weights <- 10 * fit$data$weights
  expect_equal(wesml_vcov(fit, "vcov"), wesml_vcov(fit10, "vcov"), tolerance = 1e-6)

  ed <- get_halton_normals(fit$draws_info$S, fit$draws_info$N, fit$draws_info$K_w)
  Hof <- function(wts) mxl_hessian_parallel(
    theta = fit$coefficients, X = fit$data$X, W = fit$data$W,
    alt_idx = fit$data$alt_idx, choice_idx = fit$data$choice_idx,
    M = fit$data$M, weights = wts, eta_draws = ed, rc_dist = fit$rc_dist,
    rc_correlation = fit$rc_correlation, rc_mean = fit$rc_mean,
    use_asc = fit$use_asc, include_outside_option = fit$include_outside_option)
  se1  <- sqrt(diag(solve(Hof(fit$data$weights))))
  se10 <- sqrt(diag(solve(Hof(10 * fit$data$weights))))
  expect_equal(mean(se10 / se1), 1 / sqrt(10), tolerance = 1e-3)

  # (5) under uniform weights the sandwich reduces to the usual robust variance,
  # which agrees with the inverse-Hessian for well-identified parameters. (The
  # RC-variance coordinate can be near-flat, inflating the Hessian SE, so compare
  # only the well-conditioned coordinates.)
  uni <- create_identified_mxl_data(seed = 31, N = 400L, J = 3L)
  fu <- suppressMessages(run_mxlogit(uni, "id", "alt", "choice", "x1", "w1",
                                     S = 40L, se_method = "sandwich"))
  fh <- suppressMessages(run_mxlogit(uni, "id", "alt", "choice", "x1", "w1",
                                     S = 40L, se_method = "hessian"))
  expect_true(all(is.finite(fu$se)))
  ok <- is.finite(fh$se) & fh$se < 5
  expect_gte(sum(ok), 2L)
  expect_equal(unname(fu$se[ok]), unname(fh$se[ok]), tolerance = 0.3)
})

# --- 9. provenance metadata + print label ------------------------------------

test_that("provenance: wesml from sample_by_choice, user from raw weights; labels print", {
  skip_on_cran()
  s <- make_cb_sample(seed = 41)

  fit_w <- fit_sandwich(s, S = 30L)
  expect_equal(fit_w$choice_sampling$scheme, "wesml")

  out <- capture.output(print(summary(fit_w)))
  expect_true(any(grepl("Sandwich \\(robust\\)", out)))
  expect_true(any(grepl("WESML choice-based", out)))

  # Raw user weights (no provenance attribute) -> weighting = "user"
  raw <- data.table::copy(s)
  data.table::setattr(raw, "choice_sampling", NULL)
  fit_u <- fit_sandwich(raw, S = 30L)
  expect_equal(fit_u$choice_sampling$scheme, "user")
})

# --- 10. narrow non-uniform-weight warning -----------------------------------

test_that("non-uniform weights warn under hessian but not under sandwich/uniform", {
  skip_on_cran()
  s <- make_cb_sample(seed = 51)

  expect_warning(
    run_mxlogit(s, "id", "alt", "choice", "x1", "w1", S = 20L,
                weights_col = ".wesml_weight", se_method = "hessian"),
    "sandwich"
  )
  expect_no_warning(
    suppressMessages(run_mxlogit(s, "id", "alt", "choice", "x1", "w1", S = 20L,
                                 weights_col = ".wesml_weight", se_method = "sandwich"))
  )
})

# --- input validation: zero targets, collisions, duplicate/bad names ---------

test_that("sample_by_choice rejects zero-target strata (preserves WESML target)", {
  set.seed(61)
  N <- 300L; J <- 3L
  pop <- data.table(id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N))
  pop[, x1 := rnorm(.N)]
  pop[, w1 := rnorm(.N)]
  pop[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]

  # explicit zero quota for a populated stratum -> error
  expect_error(
    sample_by_choice(pop, "id", "alt", "choice",
                     n_per_alt = c("1" = 20L, "2" = 20L, "3" = 0L)),
    "Zero sampling target"
  )
  # a fraction that rounds to zero for every stratum -> error
  expect_error(
    sample_by_choice(pop, "id", "alt", "choice", frac_per_alt = 0.001),
    "Zero sampling target"
  )
})

test_that("weight_name collisions are rejected", {
  set.seed(62)
  dt <- make_tiny()
  data.table::set(dt, j = ".wesml_weight", value = 1)   # pre-existing column

  expect_error(
    wesml_weights(dt, "id", "alt", "choice", Q = c("1" = .4, "2" = .6),
                  attach = TRUE),
    "already exists"
  )
  expect_error(
    sample_by_choice(dt, "id", "alt", "choice", n_per_alt = 2L),
    "already exists"
  )
})

test_that("malformed shares / quotas are rejected", {
  set.seed(63)
  dt <- make_tiny()
  # duplicate names in Q
  expect_error(
    wesml_weights(dt, "id", "alt", "choice", Q = c("1" = .2, "1" = .3, "2" = .5)),
    "duplicate names"
  )

  pop <- data.table(id = rep(seq_len(30), each = 3L), alt = rep(1:3, 30))
  pop[, x1 := rnorm(.N)]
  pop[, w1 := rnorm(.N)]
  pop[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  # duplicate names in n_per_alt
  expect_error(
    sample_by_choice(pop, "id", "alt", "choice",
                     n_per_alt = c("1" = 5L, "1" = 6L, "2" = 5L, "3" = 5L)),
    "duplicate names"
  )
  # non-finite quota
  expect_error(
    sample_by_choice(pop, "id", "alt", "choice", n_per_alt = NA_integer_),
    "finite"
  )
})

test_that("shares must be strictly positive and fixed quotas must be whole numbers", {
  set.seed(64)
  dt <- make_tiny()
  # A realized stratum with zero population share is incoherent for WESML.
  expect_error(
    wesml_weights(dt, "id", "alt", "choice", Q = c("1" = 0, "2" = 1)),
    "positive"
  )

  pop <- data.table(id = rep(seq_len(30), each = 3L), alt = rep(1:3, 30))
  pop[, x1 := rnorm(.N)]
  pop[, w1 := rnorm(.N)]
  pop[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  # Fractional fixed quotas must error (use frac_per_alt instead of truncating).
  expect_error(
    sample_by_choice(pop, "id", "alt", "choice", n_per_alt = 1.9),
    "whole number"
  )
  expect_error(
    sample_by_choice(pop, "id", "alt", "choice",
                     n_per_alt = c("1" = 5, "2" = 5.5, "3" = 5)),
    "whole number"
  )
  # Whole-valued doubles are still fine.
  expect_silent(
    sample_by_choice(pop, "id", "alt", "choice", n_per_alt = 5, seed = 1L)
  )
})
