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

# --- 11. provenance auto-use, error fallback, uniform downgrade --------------

test_that("WESML provenance auto-uses recorded weight column without weights_col", {
  skip_on_cran()
  skip_on_ci()
  s <- make_cb_sample(seed = 71)

  expect_message(
    fit <- suppressWarnings(
      run_mxlogit(s, "id", "alt", "choice", "x1", "w1", S = 25L,
                  se_method = "sandwich")
    ),
    "WESML"
  )
  expect_true(length(unique(fit$data$weights)) > 1L)
  expect_equal(fit$choice_sampling$scheme, "wesml")
  expect_true(isTRUE(fit$choice_sampling$weights_applied))

  # Identical to an explicit weights_col fit.
  fit_explicit <- fit_sandwich(s, S = 25L)
  expect_equal(unname(fit$coefficients), unname(fit_explicit$coefficients),
               tolerance = 1e-10)
})

test_that("WESML provenance without the weight column errors", {
  skip_on_cran()
  skip_on_ci()
  s <- make_cb_sample(seed = 72)

  dt2 <- data.table::copy(s)
  cs  <- attr(s, "choice_sampling")
  dt2[, (".wesml_weight") := NULL]
  # [, := NULL] keeps attributes, but re-set defensively.
  data.table::setattr(dt2, "choice_sampling", cs)

  expect_error(
    run_mxlogit(dt2, "id", "alt", "choice", "x1", "w1", S = 20L,
                se_method = "sandwich"),
    "provenance"
  )
})

test_that("uniform weights under WESML provenance downgrade the fit", {
  skip_on_cran()
  skip_on_ci()
  s <- make_cb_sample(seed = 73)
  N <- length(unique(s$id))

  expect_warning(
    fit <- suppressMessages(
      run_mxlogit(s, "id", "alt", "choice", "x1", "w1", S = 20L,
                  weights = rep(1, N), se_method = "sandwich")
    ),
    "NOT a WESML"
  )
  expect_equal(fit$choice_sampling$scheme, "wesml")
  expect_false(isTRUE(fit$choice_sampling$weights_applied))

  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("NOT applied", out)))
})

test_that("wesml_vcov warns under uniform weights", {
  skip_on_cran()
  skip_on_ci()
  uni <- create_identified_mxl_data(seed = 74, N = 200L, J = 3L)
  fit <- suppressMessages(
    run_mxlogit(uni, "id", "alt", "choice", "x1", "w1", S = 25L,
                se_method = "hessian")
  )
  expect_warning(wesml_vcov(fit), "uniform")
})

test_that("provenance survives the advanced input_data pathway", {
  skip_on_cran()
  skip_on_ci()
  s <- make_cb_sample(seed = 75)

  prep <- suppressWarnings(
    prepare_mxl_data(s, "id", "alt", "choice", "x1", "w1",
                     weights_col = ".wesml_weight")
  )
  expect_false(is.null(attr(prep, "choice_sampling")))
  expect_equal(attr(prep, "choice_sampling")$scheme, "wesml")

  ed  <- get_halton_normals(25L, prep$N, ncol(prep$W))
  fit <- suppressWarnings(suppressMessages(
    run_mxlogit(input_data = prep, eta_draws = ed, se_method = "sandwich")
  ))
  expect_equal(fit$choice_sampling$scheme, "wesml")
  expect_true(isTRUE(fit$choice_sampling$weights_applied))
})

test_that("scale_vars='sd' sandwich vcov matches recompute and is scale-invariant", {
  skip_on_cran()
  skip_on_ci()
  s <- make_cb_sample(seed = 76)
  ctrl <- list(xtol_rel = 1e-8, maxeval = 500L)

  fit_sd <- suppressWarnings(suppressMessages(
    run_mxlogit(s, "id", "alt", "choice", "x1", "w1", S = 30L,
                weights_col = ".wesml_weight", se_method = "sandwich",
                scale_vars = "sd", keep_data = TRUE, control = ctrl)
  ))
  # Eager (back-transformed) vcov matches the natural-space recompute.
  recomputed <- choicer:::compute_sandwich_vcov(fit_sd)
  expect_equal(unname(fit_sd$vcov), unname(recomputed$vcov), tolerance = 1e-6)

  # Invariance vs scale_vars = "none".
  fit_none <- suppressWarnings(suppressMessages(
    run_mxlogit(s, "id", "alt", "choice", "x1", "w1", S = 30L,
                weights_col = ".wesml_weight", se_method = "sandwich",
                scale_vars = "none", keep_data = TRUE, control = ctrl)
  ))
  expect_equal(unname(fit_sd$coefficients), unname(fit_none$coefficients),
               tolerance = 1e-5)
  expect_equal(unname(fit_sd$se), unname(fit_none$se), tolerance = 1e-5)
  expect_equal(unname(fit_sd$vcov), unname(fit_none$vcov), tolerance = 1e-5)
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

# =============================================================================
# MNL + NL WESML sandwich inference (parity with the MXL implementation).
# =============================================================================

# Choice-based MNL sample: alt 3 oversampled relative to its population share.
make_cb_mnl <- function(seed = 81, N = 1200L, J = 3L, n_per = 90L) {
  set.seed(seed)
  pop <- data.table(id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N))
  pop[, x1 := rnorm(.N)]
  pop[, x2 := rnorm(.N)]
  beta <- c(1.0, -0.6)
  pop[, V := drop(as.matrix(.SD) %*% beta), .SDcols = c("x1", "x2")]
  pop[, prob := exp(V) / sum(exp(V)), by = id]
  pop[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L, prob = prob)), by = id]
  pop[, c("V", "prob") := NULL]
  sample_by_choice(pop, "id", "alt", "choice", n_per_alt = n_per, seed = seed + 1L)
}

# Choice-based NL sample (4 alts, 2 nests), alt 4 oversampled.
make_cb_nl <- function(seed = 91, N = 1400L, J = 4L, n_per = 120L) {
  set.seed(seed)
  pop <- data.table(id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N))
  pop[, x1 := rnorm(.N)]
  pop[, x2 := rnorm(.N)]
  pop[, nest := ifelse(alt <= 2, "A", "B")]
  beta <- c(0.8, -0.5)
  pop[, V := drop(as.matrix(.SD) %*% beta), .SDcols = c("x1", "x2")]
  pop[, prob := exp(V) / sum(exp(V)), by = id]
  pop[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L, prob = prob)), by = id]
  pop[, c("V", "prob") := NULL]
  sample_by_choice(pop, "id", "alt", "choice", n_per_alt = n_per, seed = seed + 1L)
}

# --- (e) numeric: BHHH(weights=1) == R-side sum of outer products of scores,
#         and summed per-individual scores == negated total gradient -----------

test_that("mnl_bhhh_parallel matches outer products and scores sum to -gradient", {
  set.seed(101)
  N <- 40L; J <- 3L
  dt <- data.table(id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N),
                   x1 = rnorm(N * J), x2 = rnorm(N * J))
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  d <- prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"))
  th <- c(0.3, -0.2, 0.1, -0.1)        # 2 beta + 2 ASC

  B <- mnl_bhhh_parallel(th, d$X, d$alt_idx, d$choice_idx, d$M, d$weights)
  g <- mnl_loglik_gradient_parallel(th, d$X, d$alt_idx, d$choice_idx, d$M, d$weights)

  np <- length(th); Bman <- matrix(0, np, np); gsum <- rep(0, np)
  for (i in seq_len(N)) {
    w <- numeric(N); w[i] <- 1
    gi <- -mnl_loglik_gradient_parallel(th, d$X, d$alt_idx, d$choice_idx,
                                        d$M, w)$gradient
    Bman <- Bman + gi %*% t(gi)
    gsum <- gsum + gi
  }
  expect_equal(B, Bman, tolerance = 1e-10)
  expect_equal(gsum, -g$gradient, tolerance = 1e-10)
})

test_that("nl_bhhh_parallel matches outer products and scores sum to -gradient", {
  set.seed(107)
  N <- 50L; J <- 4L
  dt <- data.table(id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N),
                   x1 = rnorm(N * J), x2 = rnorm(N * J))
  dt[, nest := ifelse(alt <= 2, "A", "B")]
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  d <- prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
  Kl <- sum(table(d$nest_idx) > 1)
  th <- c(0.4, -0.3, rep(0.7, Kl), 0.1, -0.2, 0.05)   # 2 beta, lambda, 3 ASC

  B <- nl_bhhh_parallel(th, d$X, d$alt_idx, d$choice_idx, d$nest_idx, d$M, d$weights)
  g <- nl_loglik_gradient_parallel(th, d$X, d$alt_idx, d$choice_idx, d$nest_idx,
                                   d$M, d$weights)

  np <- length(th); Bman <- matrix(0, np, np); gsum <- rep(0, np)
  for (i in seq_len(N)) {
    w <- numeric(N); w[i] <- 1
    gi <- -nl_loglik_gradient_parallel(th, d$X, d$alt_idx, d$choice_idx,
                                       d$nest_idx, d$M, w)$gradient
    Bman <- Bman + gi %*% t(gi)
    gsum <- gsum + gi
  }
  expect_equal(B, Bman, tolerance = 1e-9)
  expect_equal(gsum, -g$gradient, tolerance = 1e-9)
})

# --- (a) sandwich runs and differs from hessian under non-uniform weights -----

test_that("MNL sandwich SEs run and differ from hessian under WESML weights", {
  skip_on_cran()
  s <- make_cb_mnl(seed = 81)

  fit_s <- suppressWarnings(suppressMessages(
    run_mnlogit(s, "id", "alt", "choice", c("x1", "x2"),
                weights_col = ".wesml_weight", se_method = "sandwich")
  ))
  fit_h <- suppressWarnings(suppressMessages(
    run_mnlogit(s, "id", "alt", "choice", c("x1", "x2"),
                weights_col = ".wesml_weight", se_method = "hessian")
  ))
  expect_false(is.null(fit_s$vcov))
  expect_true(all(is.finite(fit_s$se)))
  # Same point estimates, different SEs.
  expect_equal(unname(fit_s$coefficients), unname(fit_h$coefficients),
               tolerance = 1e-8)
  expect_false(isTRUE(all.equal(unname(fit_s$se), unname(fit_h$se),
                                tolerance = 1e-3)))
})

test_that("NL sandwich SEs run and differ from hessian under WESML weights", {
  skip_on_cran()
  s <- make_cb_nl(seed = 91)

  fit_s <- suppressWarnings(suppressMessages(
    run_nestlogit(s, "id", "alt", "choice", c("x1", "x2"), nest_col = "nest",
                  weights_col = ".wesml_weight", se_method = "sandwich")
  ))
  fit_h <- suppressWarnings(suppressMessages(
    run_nestlogit(s, "id", "alt", "choice", c("x1", "x2"), nest_col = "nest",
                  weights_col = ".wesml_weight", se_method = "hessian")
  ))
  expect_false(is.null(fit_s$vcov))
  expect_true(all(is.finite(fit_s$se)))
  expect_equal(unname(fit_s$coefficients), unname(fit_h$coefficients),
               tolerance = 1e-8)
  expect_false(isTRUE(all.equal(unname(fit_s$se), unname(fit_h$se),
                                tolerance = 1e-3)))
})

# --- (b) provenance guard errors when WESML data passed without weights -------

test_that("MNL/NL provenance guard errors when the weight column is missing", {
  skip_on_cran()
  s <- make_cb_mnl(seed = 82)
  cs <- attr(s, "choice_sampling")
  dt2 <- data.table::copy(s)
  dt2[, (".wesml_weight") := NULL]
  data.table::setattr(dt2, "choice_sampling", cs)

  expect_error(
    run_mnlogit(dt2, "id", "alt", "choice", c("x1", "x2"), se_method = "sandwich"),
    "provenance"
  )

  s_nl <- make_cb_nl(seed = 92)
  cs_nl <- attr(s_nl, "choice_sampling")
  dt3 <- data.table::copy(s_nl)
  dt3[, (".wesml_weight") := NULL]
  data.table::setattr(dt3, "choice_sampling", cs_nl)

  expect_error(
    run_nestlogit(dt3, "id", "alt", "choice", c("x1", "x2"), nest_col = "nest",
                  se_method = "sandwich"),
    "provenance"
  )
})

# --- (c) uniform-weight warning fires (non-sandwich) --------------------------

test_that("MNL/NL non-uniform weights warn under hessian but not sandwich", {
  skip_on_cran()
  s <- make_cb_mnl(seed = 83)
  expect_warning(
    suppressMessages(run_mnlogit(s, "id", "alt", "choice", c("x1", "x2"),
                                 weights_col = ".wesml_weight", se_method = "hessian")),
    "sandwich"
  )
  expect_no_warning(
    suppressMessages(run_mnlogit(s, "id", "alt", "choice", c("x1", "x2"),
                                 weights_col = ".wesml_weight", se_method = "sandwich"))
  )

  s_nl <- make_cb_nl(seed = 93)
  expect_warning(
    suppressMessages(run_nestlogit(s_nl, "id", "alt", "choice", c("x1", "x2"),
                                   nest_col = "nest", weights_col = ".wesml_weight",
                                   se_method = "hessian")),
    "sandwich"
  )
})

test_that("MNL/NL bhhh + non-uniform weights warn that BHHH is not a WESML correction", {
  skip_on_cran()
  s <- make_cb_mnl(seed = 87)
  mnl_warn <- NULL
  withCallingHandlers(
    suppressMessages(run_mnlogit(s, "id", "alt", "choice", c("x1", "x2"),
                                 weights_col = ".wesml_weight", se_method = "bhhh")),
    warning = function(w) {
      mnl_warn <<- c(mnl_warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # The bhhh-specific message fires (mentions BHHH/OPG and the w^1 vs w^2 meat).
  expect_true(any(grepl("BHHH/OPG", mnl_warn)))
  expect_true(any(grepl("w\\^2", mnl_warn)))
  expect_true(any(grepl("sandwich", mnl_warn)))
  # The generic steering message does NOT fire redundantly alongside it.
  expect_false(any(grepl("If these are sampling/WESML", mnl_warn)))

  s_nl <- make_cb_nl(seed = 97)
  nl_warn <- NULL
  withCallingHandlers(
    suppressMessages(run_nestlogit(s_nl, "id", "alt", "choice", c("x1", "x2"),
                                   nest_col = "nest", weights_col = ".wesml_weight",
                                   se_method = "bhhh")),
    warning = function(w) {
      nl_warn <<- c(nl_warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("BHHH/OPG", nl_warn)))
  expect_true(any(grepl("w\\^2", nl_warn)))
  expect_true(any(grepl("sandwich", nl_warn)))
  expect_false(any(grepl("If these are sampling/WESML", nl_warn)))
})

test_that("MNL/NL wesml_vcov warns under uniform weights", {
  skip_on_cran()
  dt <- create_small_mnl_data(seed = 84)
  fit <- suppressMessages(run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2")))
  expect_warning(wesml_vcov(fit), "uniform")

  dnl <- create_small_nl_data(seed = 94)
  fit_nl <- suppressMessages(suppressWarnings(
    run_nestlogit(dnl, "id", "alt", "choice", c("x1", "x2"), nest_col = "nest")
  ))
  expect_warning(wesml_vcov(fit_nl), "uniform")
})

# --- (d) wesml_vcov() round-trips and equals the eager sandwich vcov ----------

test_that("MNL wesml_vcov() equals the eager fit-path sandwich vcov", {
  skip_on_cran()
  s <- make_cb_mnl(seed = 85)
  fit <- suppressWarnings(suppressMessages(
    run_mnlogit(s, "id", "alt", "choice", c("x1", "x2"),
                weights_col = ".wesml_weight", se_method = "sandwich",
                keep_data = TRUE)
  ))
  expect_equal(unname(wesml_vcov(fit, "vcov")), unname(fit$vcov), tolerance = 1e-8)
  # Scale-invariance: rescale stored weights -> sandwich vcov unchanged.
  fit10 <- fit
  fit10$data$weights <- 10 * fit$data$weights
  expect_equal(wesml_vcov(fit, "vcov"), wesml_vcov(fit10, "vcov"), tolerance = 1e-7)
})

test_that("NL wesml_vcov() equals the eager fit-path sandwich vcov", {
  skip_on_cran()
  s <- make_cb_nl(seed = 95)
  fit <- suppressWarnings(suppressMessages(
    run_nestlogit(s, "id", "alt", "choice", c("x1", "x2"), nest_col = "nest",
                  weights_col = ".wesml_weight", se_method = "sandwich",
                  keep_data = TRUE)
  ))
  expect_equal(unname(wesml_vcov(fit, "vcov")), unname(fit$vcov), tolerance = 1e-8)
  fit10 <- fit
  fit10$data$weights <- 10 * fit$data$weights
  expect_equal(wesml_vcov(fit, "vcov"), wesml_vcov(fit10, "vcov"), tolerance = 1e-7)
})

# --- provenance / print labels -----------------------------------------------

test_that("MNL/NL provenance recorded and summary prints sandwich + WESML labels", {
  skip_on_cran()
  s <- make_cb_mnl(seed = 86)
  fit <- suppressWarnings(suppressMessages(
    run_mnlogit(s, "id", "alt", "choice", c("x1", "x2"),
                weights_col = ".wesml_weight", se_method = "sandwich")
  ))
  expect_equal(fit$choice_sampling$scheme, "wesml")
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Sandwich \\(robust\\)", out)))
  expect_true(any(grepl("WESML choice-based", out)))

  s_nl <- make_cb_nl(seed = 96)
  fit_nl <- suppressWarnings(suppressMessages(
    run_nestlogit(s_nl, "id", "alt", "choice", c("x1", "x2"), nest_col = "nest",
                  weights_col = ".wesml_weight", se_method = "sandwich")
  ))
  expect_equal(fit_nl$choice_sampling$scheme, "wesml")
  out_nl <- capture.output(print(summary(fit_nl)))
  expect_true(any(grepl("Sandwich \\(robust\\)", out_nl)))
  expect_true(any(grepl("WESML choice-based", out_nl)))
})

# --- (g) Fix A: weights must be finite and strictly positive -----------------

test_that("prepare_*_data reject non-positive and non-finite weights", {
  mnl <- data.table(
    id     = rep(1:4, each = 2),
    alt    = rep(c(1L, 2L), 4),
    x1     = rnorm(8),
    choice = c(1, 0, 1, 0, 0, 1, 0, 1)
  )
  nl <- data.table::copy(mnl)
  nl[, nest := ifelse(alt <= 1, "A", "B")]
  mxl <- data.table::copy(mnl)
  mxl[, w1 := rnorm(8)]

  # Zero weight.
  expect_error(
    prepare_mnl_data(mnl, "id", "alt", "choice", "x1", weights = c(1, 0, 1, 1)),
    "strictly positive"
  )
  expect_error(
    prepare_nl_data(nl, "id", "alt", "choice", "x1", "nest",
                    weights = c(1, 0, 1, 1)),
    "strictly positive"
  )
  expect_error(
    prepare_mxl_data(mxl, "id", "alt", "choice", "x1", "w1",
                     weights = c(1, 0, 1, 1)),
    "strictly positive"
  )

  # Negative weight.
  expect_error(
    prepare_mnl_data(mnl, "id", "alt", "choice", "x1", weights = c(1, -2, 1, 1)),
    "strictly positive"
  )
  expect_error(
    prepare_nl_data(nl, "id", "alt", "choice", "x1", "nest",
                    weights = c(1, -2, 1, 1)),
    "strictly positive"
  )
  expect_error(
    prepare_mxl_data(mxl, "id", "alt", "choice", "x1", "w1",
                     weights = c(1, -2, 1, 1)),
    "strictly positive"
  )

  # Non-finite weight (NA and Inf).
  expect_error(
    prepare_mnl_data(mnl, "id", "alt", "choice", "x1", weights = c(1, NA, 1, 1)),
    "finite"
  )
  expect_error(
    prepare_nl_data(nl, "id", "alt", "choice", "x1", "nest",
                    weights = c(1, Inf, 1, 1)),
    "finite"
  )
  expect_error(
    prepare_mxl_data(mxl, "id", "alt", "choice", "x1", "w1",
                     weights = c(1, NA, 1, 1)),
    "finite"
  )

  # Uniform default weights are accepted.
  expect_silent(prepare_mnl_data(mnl, "id", "alt", "choice", "x1"))
})

# --- (h) Fix B: advanced input_data + WESML provenance + uniform weights errors

test_that("advanced-mode WESML provenance with uniform weights errors (MNL/NL/MXL)", {
  skip_on_cran()
  skip_on_ci()

  # MNL: prepare WITHOUT weights (uniform), then attach WESML provenance.
  s_mnl <- make_cb_mnl(seed = 84)
  cs_mnl <- attr(s_mnl, "choice_sampling")
  prep_mnl <- prepare_mnl_data(s_mnl, "id", "alt", "choice", c("x1", "x2"))
  attr(prep_mnl, "choice_sampling") <- cs_mnl
  expect_error(
    run_mnlogit(input_data = prep_mnl, se_method = "sandwich"),
    "choice_sampling"
  )
  # Stripping the provenance lets the unweighted fit proceed.
  attr(prep_mnl, "choice_sampling") <- NULL
  expect_silent(suppressMessages(
    run_mnlogit(input_data = prep_mnl, se_method = "hessian")
  ))

  # NL.
  s_nl <- make_cb_nl(seed = 94)
  cs_nl <- attr(s_nl, "choice_sampling")
  prep_nl <- prepare_nl_data(s_nl, "id", "alt", "choice", c("x1", "x2"), "nest")
  attr(prep_nl, "choice_sampling") <- cs_nl
  expect_error(
    run_nestlogit(input_data = prep_nl, se_method = "sandwich"),
    "choice_sampling"
  )
  attr(prep_nl, "choice_sampling") <- NULL
  expect_silent(suppressMessages(
    run_nestlogit(input_data = prep_nl, se_method = "hessian")
  ))

  # MXL.
  s_mxl <- make_cb_sample(seed = 76)
  cs_mxl <- attr(s_mxl, "choice_sampling")
  prep_mxl <- prepare_mxl_data(s_mxl, "id", "alt", "choice", "x1", "w1")
  attr(prep_mxl, "choice_sampling") <- cs_mxl
  ed <- get_halton_normals(20L, prep_mxl$N, ncol(prep_mxl$W))
  expect_error(
    run_mxlogit(input_data = prep_mxl, eta_draws = ed, se_method = "sandwich"),
    "choice_sampling"
  )
  attr(prep_mxl, "choice_sampling") <- NULL
  expect_silent(suppressMessages(
    run_mxlogit(input_data = prep_mxl, eta_draws = ed, se_method = "hessian")
  ))
})

test_that("advanced-mode WESML provenance with non-uniform weights still fits", {
  skip_on_cran()
  skip_on_ci()

  # Weights baked into input_data (non-uniform) + provenance => valid WESML fit.
  s_mnl <- make_cb_mnl(seed = 85)
  prep_mnl <- prepare_mnl_data(s_mnl, "id", "alt", "choice", c("x1", "x2"),
                               weights_col = ".wesml_weight")
  expect_false(is.null(attr(prep_mnl, "choice_sampling")))
  fit_mnl <- suppressWarnings(suppressMessages(
    run_mnlogit(input_data = prep_mnl, se_method = "sandwich")
  ))
  expect_equal(fit_mnl$choice_sampling$scheme, "wesml")
  expect_true(isTRUE(fit_mnl$choice_sampling$weights_applied))

  s_nl <- make_cb_nl(seed = 95)
  prep_nl <- prepare_nl_data(s_nl, "id", "alt", "choice", c("x1", "x2"), "nest",
                             weights_col = ".wesml_weight")
  fit_nl <- suppressWarnings(suppressMessages(
    run_nestlogit(input_data = prep_nl, se_method = "sandwich")
  ))
  expect_equal(fit_nl$choice_sampling$scheme, "wesml")
  expect_true(isTRUE(fit_nl$choice_sampling$weights_applied))
})
