# Score kernels + clustered / robust / BHHH variance assembly (v0.2.0, Day 1).
#
# The per-situation score kernels (mnl/nl/mxl_scores_parallel) are internal;
# tests reach them via choicer:::. The R-side assembly is exercised through
# vcov(fit, type = , cluster = ) and the cluster_col fit-time convenience.

# --- Shared fixtures ---------------------------------------------------------

make_panel_mnl <- function(N_persons = 30L, T_per = 4L, J = 3L, seed = 42) {
  set.seed(seed)
  N <- N_persons * T_per
  dt <- data.table::data.table(
    id = rep(seq_len(N), each = J),
    alt = rep(seq_len(J), N)
  )
  dt[, `:=`(x1 = stats::rnorm(.N), x2 = stats::rnorm(.N))]
  dt[, person := rep(rep(seq_len(N_persons), each = T_per), each = J)[seq_len(.N)]]
  beta_true <- c(0.8, -0.5)
  dt[, V := 0.8 * x1 - 0.5 * x2]
  dt[, prob := exp(V) / sum(exp(V)), by = id]
  dt[, choice := as.integer(alt == sample(alt, 1, prob = prob)), by = id]
  dt[, c("V", "prob") := NULL]
  dt
}

fit_mnl <- function(dt, ...) {
  suppressMessages(run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), ...
  ))
}

# Per-situation cluster labels named by choice-situation id. Because the
# labels carry their ids as names, vcov() realigns them to the prepared order
# safely, regardless of the vector's own order.
situation_clusters <- function(dt) {
  m <- unique(dt[, .(id, person)])
  data.table::setorder(m, id)
  stats::setNames(m$person, m$id)
}

# --- Scores vs numDeriv per-situation gradients ------------------------------

test_that("MNL score rows match numDeriv per-situation gradients", {
  skip_if_not_installed("numDeriv")
  dt <- make_panel_mnl(N_persons = 8L, T_per = 2L)
  fit <- fit_mnl(dt)
  d <- fit$data
  S <- choicer:::mnl_scores_parallel(
    theta = fit$coefficients, X = d$X, alt_idx = d$alt_idx,
    choice_idx = d$choice_idx, M = d$M, use_asc = fit$use_asc,
    include_outside_option = fit$include_outside_option
  )
  expect_equal(dim(S), c(fit$nobs, fit$n_params))

  # Situation i's positive loglik via an indicator weight vector.
  loglik_i <- function(theta, i) {
    w <- rep(0, fit$nobs); w[i] <- 1
    -mnl_loglik_gradient_parallel(
      theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
      M = d$M, weights = w, use_asc = fit$use_asc,
      include_outside_option = fit$include_outside_option
    )$objective
  }
  for (i in c(1L, 5L, fit$nobs)) {
    g_num <- numDeriv::grad(loglik_i, fit$coefficients, i = i)
    expect_equal(unname(S[i, ]), g_num, tolerance = 1e-6)
  }
})

test_that("MXL score rows match numDeriv per-situation gradients", {
  skip_if_not_installed("numDeriv")
  set.seed(7)
  N <- 40L; J <- 3L
  dt <- data.table::data.table(
    id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N)
  )
  dt[, `:=`(x1 = stats::rnorm(.N), w1 = stats::rnorm(.N))]
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  fit <- suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = "w1", S = 30L, rc_mean = TRUE
  ))
  d <- fit$data
  gp <- choicer:::.mxl_gen_params(fit$draws_info)
  S_mat <- choicer:::mxl_scores_parallel(
    theta = fit$coefficients, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx, M = d$M,
    eta_draws = gp$eta_draws, rc_dist = fit$rc_dist,
    rc_correlation = fit$rc_correlation, rc_mean = fit$rc_mean,
    use_asc = fit$use_asc,
    include_outside_option = fit$include_outside_option,
    gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S
  )
  loglik_i <- function(theta, i) {
    w <- rep(0, fit$nobs); w[i] <- 1
    -mxl_loglik_gradient_parallel(
      theta = theta, X = d$X, W = d$W, alt_idx = d$alt_idx,
      choice_idx = d$choice_idx, M = d$M, weights = w,
      eta_draws = gp$eta_draws, rc_dist = fit$rc_dist,
      rc_correlation = fit$rc_correlation, rc_mean = fit$rc_mean,
      use_asc = fit$use_asc,
      include_outside_option = fit$include_outside_option,
      gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S
    )$objective
  }
  for (i in c(1L, 17L, N)) {
    g_num <- numDeriv::grad(loglik_i, fit$coefficients, i = i)
    expect_equal(unname(S_mat[i, ]), g_num, tolerance = 1e-6)
  }
})

test_that("NL score rows match numDeriv per-situation gradients", {
  skip_if_not_installed("numDeriv")
  set.seed(11)
  N <- 40L; J <- 4L
  dt <- data.table::data.table(
    id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N)
  )
  dt[, `:=`(x1 = stats::rnorm(.N), x2 = stats::rnorm(.N))]
  dt[, nest := ifelse(alt <= 2, "A", "B")]
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  fit <- suppressMessages(run_nestlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), nest_col = "nest"
  ))
  d <- fit$data
  # Uniformly random choices leave lambda weakly identified, so its MLE can
  # approach the zero boundary differently across platforms. The numerical
  # oracle uses central differences and must stay in the lambda > 0 domain;
  # test the score identity at a fixed, feasible interior point instead.
  theta_eval <- fit$coefficients
  theta_eval[fit$param_map$lambda] <- 0.7
  S_mat <- choicer:::nl_scores_parallel(
    theta = theta_eval, X = d$X, alt_idx = d$alt_idx,
    choice_idx = d$choice_idx, nest_idx = d$nest_idx, M = d$M,
    use_asc = fit$use_asc,
    include_outside_option = fit$include_outside_option
  )
  loglik_i <- function(theta, i) {
    w <- rep(0, fit$nobs); w[i] <- 1
    -nl_loglik_gradient_parallel(
      theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
      nest_idx = d$nest_idx, M = d$M, weights = w, use_asc = fit$use_asc,
      include_outside_option = fit$include_outside_option
    )$objective
  }
  for (i in c(1L, 23L, N)) {
    g_num <- numDeriv::grad(loglik_i, theta_eval, i = i)
    expect_equal(unname(S_mat[i, ]), g_num, tolerance = 1e-6)
  }
})

# --- Assembly identities -----------------------------------------------------

test_that("crossprod(sqrt(w) S) reproduces the BHHH kernel", {
  dt <- make_panel_mnl()
  fit <- fit_mnl(dt)
  d <- fit$data
  S <- choicer:::compute_scores(fit)
  B_kernel <- mnl_bhhh_parallel(
    theta = fit$coefficients, X = d$X, alt_idx = d$alt_idx,
    choice_idx = d$choice_idx, M = d$M, weights = d$weights,
    use_asc = fit$use_asc,
    include_outside_option = fit$include_outside_option
  )
  expect_equal(crossprod(sqrt(d$weights) * S), B_kernel,
               ignore_attr = TRUE, tolerance = 1e-10)
})

test_that("cluster of singletons equals the robust variance", {
  dt <- make_panel_mnl()
  fit <- fit_mnl(dt)
  V_rob <- vcov(fit, type = "robust")
  # Each situation its own cluster, named by situation id.
  singletons <- stats::setNames(seq_len(fit$nobs), fit$data$situation_ids)
  V_singleton <- vcov(fit, type = "cluster", cluster = singletons)
  expect_equal(V_singleton, V_rob, tolerance = 1e-10)
})

test_that("a single cluster gives the outer product of the summed score", {
  dt <- make_panel_mnl()
  fit <- fit_mnl(dt)
  S <- choicer:::compute_scores(fit)
  w <- fit$data$weights
  g <- colSums(w * S)
  B_single <- choicer:::.score_meat(S, w, "cluster",
                                    cluster = rep(1L, fit$nobs))
  expect_equal(B_single, tcrossprod(g), ignore_attr = TRUE,
               tolerance = 1e-10)
  # At the MLE the summed weighted score is ~0, so the meat nearly vanishes.
  expect_lt(max(abs(g)), 1e-3)
})

test_that("cluster meat is invariant to label ordering and coding", {
  dt <- make_panel_mnl()
  fit <- fit_mnl(dt)
  cl_int <- situation_clusters(dt)         # named by situation id
  V1 <- vcov(fit, cluster = cl_int)  # cluster= alone implies type = "cluster"
  # Relabel: reverse the integer coding (names ride along)
  V2 <- vcov(fit, type = "cluster", cluster = max(cl_int) + 1L - cl_int)
  # Recode as unordered character labels, keeping the id names
  pool <- sample(paste0("grp_", seq_len(max(cl_int))))
  V3 <- vcov(fit, type = "cluster",
             cluster = stats::setNames(pool[cl_int], names(cl_int)))
  expect_equal(V1, V2, tolerance = 1e-12)
  expect_equal(V1, V3, tolerance = 1e-12)
})

test_that("post-hoc robust vcov matches wesml_vcov and post-hoc bhhh matches the kernel", {
  dt <- make_panel_mnl()
  fit <- fit_mnl(dt)
  expect_equal(vcov(fit, type = "robust"),
               suppressWarnings(wesml_vcov(fit)), tolerance = 1e-8)
  V_bhhh <- vcov(fit, type = "bhhh")
  d <- fit$data
  B_kernel <- mnl_bhhh_parallel(
    theta = fit$coefficients, X = d$X, alt_idx = d$alt_idx,
    choice_idx = d$choice_idx, M = d$M, weights = d$weights,
    use_asc = fit$use_asc,
    include_outside_option = fit$include_outside_option
  )
  expect_equal(V_bhhh, solve(B_kernel), ignore_attr = TRUE,
               tolerance = 1e-8)
  # type = "hessian" reproduces the as-fitted default vcov
  expect_equal(vcov(fit, type = "hessian"), vcov(fit), tolerance = 1e-8)
})

test_that("clustered SEs exceed robust SEs under strong within-cluster dependence", {
  # Duplicate every situation 4x within a person-cluster: scores are perfectly
  # dependent within clusters, so clustered variances must be larger.
  set.seed(3)
  N_base <- 60L; J <- 3L
  base <- data.table::data.table(
    id = rep(seq_len(N_base), each = J), alt = rep(seq_len(J), N_base)
  )
  base[, `:=`(x1 = stats::rnorm(.N), x2 = stats::rnorm(.N))]
  base[, V := 0.8 * x1 - 0.5 * x2]
  base[, prob := exp(V) / sum(exp(V)), by = id]
  base[, choice := as.integer(alt == sample(alt, 1, prob = prob)), by = id]
  base[, c("V", "prob") := NULL]
  reps <- data.table::rbindlist(lapply(0:3, function(r) {
    b <- data.table::copy(base)
    b[, id := id + r * N_base]
    b
  }))
  reps[, person := ((id - 1L) %% N_base) + 1L]
  fit <- fit_mnl(reps)
  cl <- situation_clusters(reps)
  se_rob <- sqrt(diag(vcov(fit, type = "robust")))
  se_cl <- sqrt(diag(vcov(fit, type = "cluster", cluster = cl)))
  expect_true(all(se_cl > se_rob))
  # With j perfectly-dependent copies, the cluster meat is j^2 x the
  # per-copy meat while robust is j x, so V_cl ~ 4 * V_rob here.
  expect_equal(unname(se_cl / se_rob), rep(2, fit$n_params),
               tolerance = 0.05)
})

# --- Fit-time convenience (cluster_col / se_method = "cluster") --------------

test_that("cluster_col at fit time matches post-hoc clustering", {
  dt <- make_panel_mnl()
  fit_plain <- fit_mnl(dt)
  fit_cl <- fit_mnl(dt, cluster_col = "person")  # implies se_method = "cluster"
  expect_identical(fit_cl$se_method, "cluster")
  cl <- situation_clusters(dt)
  V_posthoc <- vcov(fit_plain, type = "cluster", cluster = cl)
  expect_equal(vcov(fit_cl), V_posthoc, tolerance = 1e-6)
  # Stored labels allow type = "cluster" without re-supplying cluster
  expect_equal(vcov(fit_cl, type = "cluster"), V_posthoc, tolerance = 1e-6)
  # Point estimates are unaffected by the SE method
  expect_equal(coef(fit_cl), coef(fit_plain))
})

test_that("cluster_col works for NL and MXL fits", {
  # Well-conditioned fixture: 50 clusters identify the 7 NL parameters with
  # room to spare, so vcov entries are small and well-scaled. (A small,
  # ill-conditioned fixture makes the fit-time-vs-post-hoc comparison fragile:
  # the parallel Hessian bread differs call-to-call by ~1e-8, which a poorly
  # identified vcov amplifies into a large relative wobble.)
  set.seed(21)
  N <- 200L; J <- 4L
  dt <- data.table::data.table(
    id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N)
  )
  dt[, `:=`(x1 = stats::rnorm(.N), x2 = stats::rnorm(.N))]
  dt[, nest := ifelse(alt <= 2, "A", "B")]
  dt[, person := rep(seq_len(50L), each = 4L)[id]]
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]

  fit_nl <- suppressMessages(run_nestlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), nest_col = "nest",
    cluster_col = "person"
  ))
  expect_identical(fit_nl$se_method, "cluster")
  expect_true(all(is.finite(fit_nl$se)))
  cl <- situation_clusters(dt)
  # Fit-time and post-hoc are the same estimator; they differ only by the
  # parallel-reduction noise in the bread, so allow a margin above it.
  expect_equal(vcov(fit_nl), vcov(fit_nl, type = "cluster", cluster = cl),
               tolerance = 1e-5)

  fit_mxl <- suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = "x2", S = 30L, rc_mean = TRUE,
    cluster_col = "person"
  ))
  expect_identical(fit_mxl$se_method, "cluster")
  expect_equal(vcov(fit_mxl), vcov(fit_mxl, type = "cluster", cluster = cl),
               tolerance = 1e-5)
})

# --- Error handling ----------------------------------------------------------

test_that("cluster errors are informative", {
  dt <- make_panel_mnl()
  fit <- fit_mnl(dt)
  expect_error(vcov(fit, type = "cluster"),
               "need cluster labels")
  expect_error(vcov(fit, type = "cluster", cluster = 1:5),
               "choice situations")
  cl <- situation_clusters(dt)
  cl[1] <- NA
  expect_error(vcov(fit, type = "cluster", cluster = cl),
               "missing values")
  expect_error(fit_mnl(dt, se_method = "cluster"),
               "needs cluster labels")
  expect_error(
    suppressMessages(run_mnlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), cluster_col = "nope"
    )),
    "Missing columns"
  )
  # cluster_col must be constant within a choice situation
  dt_bad <- data.table::copy(dt)
  dt_bad[, person := seq_len(.N)]
  expect_error(
    suppressMessages(run_mnlogit(
      data = dt_bad, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), cluster_col = "person"
    )),
    "constant within"
  )
  # keep_data = FALSE blocks post-hoc recomputation
  fit_nodata <- fit_mnl(dt, keep_data = FALSE)
  expect_error(vcov(fit_nodata, type = "robust"), "keep_data")
})

# --- Named-cluster realignment (the alignment footgun) -----------------------

test_that("named cluster realigns by id, so its order does not matter", {
  dt <- make_panel_mnl()
  fit <- fit_mnl(dt)
  cl <- situation_clusters(dt)                 # named by id, prepared order
  V_ref <- vcov(fit, type = "cluster", cluster = cl)
  # Shuffle the vector: names still identify each situation.
  V_shuf <- vcov(fit, type = "cluster", cluster = cl[sample.int(length(cl))])
  expect_equal(V_ref, V_shuf, tolerance = 1e-12)
})

test_that("a named cluster recovers the fit-time cluster_col result even when the data is not id-ordered", {
  # Scramble row order so first-appearance (data) order != prepared (id) order.
  dt <- make_panel_mnl()
  set.seed(99)
  dt_scrambled <- dt[sample.int(nrow(dt))]
  fit_col <- fit_mnl(dt_scrambled, cluster_col = "person")  # ground truth
  fit_plain <- fit_mnl(dt_scrambled)

  # A per-situation vector built in DATA order (first appearance of each id).
  first_seen <- dt_scrambled[!duplicated(id), .(id, person)]
  data_order <- stats::setNames(first_seen$person, first_seen$id)
  # Named -> realigned by id -> matches the fit-time truth despite its order.
  expect_equal(vcov(fit_plain, type = "cluster", cluster = data_order),
               vcov(fit_col), tolerance = 1e-6)

  # The same labels stripped of names and taken positionally are NOT the
  # fit-time truth (this is exactly the silent-misalignment hazard that names
  # remove); the call warns about the prepared-order assumption.
  V_positional <- suppressWarnings(
    vcov(fit_plain, type = "cluster", cluster = unname(data_order))
  )
  expect_false(isTRUE(all.equal(V_positional, vcov(fit_col), tolerance = 1e-6)))
})

test_that("unnamed cluster warns about the prepared-order assumption", {
  dt <- make_panel_mnl()
  fit <- fit_mnl(dt)
  cl <- unname(situation_clusters(dt))         # prepared order, but unnamed
  expect_warning(vcov(fit, type = "cluster", cluster = cl),
                 "prepared \\(id-sorted\\) order")
  # In prepared order it still gives the right answer (equal to the named form).
  V_unnamed <- suppressWarnings(vcov(fit, type = "cluster", cluster = cl))
  V_named <- vcov(fit, type = "cluster", cluster = situation_clusters(dt))
  expect_equal(V_unnamed, V_named, tolerance = 1e-12)
})

test_that("cluster resolution rejects row-level and incomplete vectors", {
  dt <- make_panel_mnl()
  fit <- fit_mnl(dt)
  # Row-level (length sum(M)) labels: the most common mistake.
  row_level <- dt$person
  expect_error(vcov(fit, type = "cluster", cluster = row_level),
               "stacked")
  # Named but not covering every situation.
  cl <- situation_clusters(dt)
  bad <- cl[-1]
  expect_error(vcov(fit, type = "cluster", cluster = bad),
               "do not cover")
  # Duplicate names.
  dup <- cl
  names(dup)[2] <- names(dup)[1]
  expect_error(vcov(fit, type = "cluster", cluster = dup),
               "duplicate names")
})

# --- Coverage: outside option, weights, and scale_vars with clustering -------

test_that("clustering works with an outside option", {
  set.seed(11)
  N <- 80L; J <- 3L
  dt <- data.table::data.table(
    id = rep(seq_len(N), each = J), alt = rep(seq_len(J), N)
  )
  dt[, `:=`(x1 = stats::rnorm(.N), x2 = stats::rnorm(.N))]
  dt[, person := rep(seq_len(20L), each = 4L)[id]]
  # About a third of situations pick the (implicit) outside option.
  dt[, V := 0.6 * x1 - 0.4 * x2]
  dt[, prob := exp(V) / (1 + sum(exp(V))), by = id]
  dt[, choice := 0L]
  dt[, choice := {
    p <- c(1 - sum(prob), prob)
    pick <- sample.int(J + 1L, 1L, prob = p)
    out <- rep(0L, .N); if (pick > 1L) out[pick - 1L] <- 1L; out
  }, by = id]
  fit <- suppressMessages(run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), include_outside_option = TRUE,
    cluster_col = "person"
  ))
  expect_true(all(is.finite(fit$se)))
  cl <- situation_clusters(dt)
  expect_equal(vcov(fit), vcov(fit, type = "cluster", cluster = cl),
               tolerance = 1e-6)
})

test_that("clustering composes with weights and with scale_vars", {
  dt <- make_panel_mnl()
  dt[, wt := 0.5 + (person %% 3)]              # non-uniform, constant within id

  # Weighted + clustered: fit-time equals post-hoc.
  fit_w <- suppressWarnings(suppressMessages(run_mnlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), weights_col = "wt", cluster_col = "person"
  )))
  cl <- situation_clusters(dt)
  expect_equal(vcov(fit_w), vcov(fit_w, type = "cluster", cluster = cl),
               tolerance = 1e-6)

  # scale_vars is unwound to natural units, so a scaled clustered fit matches
  # the natural-scale post-hoc cluster variance of an unscaled fit.
  fit_scaled <- fit_mnl(dt, scale_vars = "sd", cluster_col = "person")
  fit_plain <- fit_mnl(dt)
  expect_equal(vcov(fit_scaled),
               vcov(fit_plain, type = "cluster", cluster = cl),
               tolerance = 1e-6)
})
