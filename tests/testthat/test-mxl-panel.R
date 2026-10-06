# Integration tests for the panel mixed logit: run_mxlogit(person_col = ).
#
# With `person_col`, each decision maker draws one taste vector that all of
# their choice situations share, and the simulated likelihood integrates the
# product of their logit probabilities (Revelt and Train 1998). The kernels are
# checked against a brute-force oracle in test-mxl-panel-kernels.R. This file
# checks the R layer around them:
#   1. prepare_mxl_data(person_col = ): row order, Ti / person_ids, errors;
#   2. every se_method on a panel fit;
#   3. post-hoc vcov() and ensure_vcov() against the fit-time variance in
#      both draw modes. In generate mode, a recomputation that drops `Ti`
#      would silently return cross-sectional derivatives;
#   4. equivalences: single-situation panels, clustering on the decision
#      maker, scale_vars;
#   5. the panel likelihood is the one maximized;
#   6. post-estimation on a panel fit;
#   7. conditional_tastes(), against the oracle in helper-mxl-panel.R;
#   8. simulate_mxl_data(T = );
#   9. parameter and taste recovery.
# Shared fits are memoized and small; only recovery and the scale_vars
# invariance check are skipped on CRAN.

# --- Fixtures ----------------------------------------------------------------

# Regexp fragments of the panel validation errors (not full messages).
pnl_rx_within <- "constant within|vary within|varies within|differ"
pnl_rx_cluster <- paste0(pnl_rx_within, "|nest")

# Fit-time se_method -> the vcov(type = ) that recomputes it post hoc.
pnl_se_types <- c(hessian = "hessian", bhhh = "bhhh", sandwich = "robust",
                  cluster = "cluster")

# Memoize a fixture builder: build() runs on first use only.
pnl_memo <- function(build) {
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- build()
    cache
  }
}

# Long-format panel with person-level tastes. Decision maker p makes T_p
# choices among J alternatives (T_p drawn from T_range unless given). The
# utility is -x1 + b_p1 w1 + b_p2 w2 + ASC, with b_p1 ~ N(0.5, 1) and
# b_p2 ~ N(-0.5, 0.8^2) fixed across p's situations. Situation ids interleave
# decision makers (the t-th situation of every person is numbered before any
# (t + 1)-th), so the prepared (person, id) order differs from the id order,
# and rows are returned shuffled.
pnl_data <- function(n_persons = 90L, T_range = 3:5, T_p = NULL, J = 3L,
                     seed = 11L) {
  set.seed(seed)
  if (is.null(T_p)) {
    T_p <- T_range[sample.int(length(T_range), n_persons, replace = TRUE)]
  }
  T_p <- as.integer(T_p)
  n_persons <- length(T_p)
  person <- rep(seq_len(n_persons), T_p)      # situations, person-major
  N <- length(person)
  id <- integer(N)
  id[order(sequence(T_p), person)] <- seq_len(N)
  dt <- data.table::data.table(
    id = rep(id, each = J),
    person = rep(person, each = J),
    alt = rep(seq_len(J), N)
  )
  dt[, `:=`(x1 = stats::rnorm(.N), w1 = stats::rnorm(.N),
            w2 = stats::rnorm(.N))]
  taste <- cbind(0.5 + stats::rnorm(n_persons),
                 -0.5 + 0.8 * stats::rnorm(n_persons))
  asc <- c(0, 0.4, -0.3, 0.2)[seq_len(J)]
  dt[, v := -x1 + taste[person, 1L] * w1 + taste[person, 2L] * w2 + asc[alt]]
  dt[, choice := {
    p <- exp(v - max(v))
    as.integer(seq_len(.N) == sample.int(.N, 1L, prob = p))
  }, by = id]
  dt[, v := NULL]
  dt[sample.int(.N)]
}

pnl_prepare <- function(dt, ...) {
  prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("w1", "w2"), ...)
}

pnl_fit <- function(dt, ..., S = 30L) {
  suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    rc_mean = TRUE, S = S, ...
  ))
}

# Custom optimizer that evaluates the objective once and returns theta_init,
# so a fit's eager variance is computed at a chosen parameter vector without
# re-optimizing.
pnl_at <- function(theta_init, eval_f, ...) {
  list(par = theta_init, value = eval_f(theta_init)$objective,
       convergence = 0L)
}

# One label per choice situation, named by situation id.
pnl_labels <- function(dt, col) {
  m <- unique(dt[, c("id", col), with = FALSE])
  stats::setNames(m[[col]], m$id)
}

# Call an MXL kernel on a fit's stored data at its coefficients, with the
# fit's own draws. `panel = TRUE` passes the stored Ti (one draw block per
# decision maker); `panel = FALSE` treats every situation as its own unit.
pnl_kernel <- function(kernel, fit, panel = TRUE) {
  d <- fit$data
  di <- fit$draws_info
  generate <- identical(di$mode, "generate")
  n_blocks <- if (panel) length(d$Ti) else length(d$M)
  args <- list(theta = fit$coefficients, X = d$X, W = d$W,
               alt_idx = d$alt_idx, choice_idx = d$choice_idx, M = d$M)
  if (kernel %in% c("gradient", "hessian", "bhhh")) args$weights <- d$weights
  args$eta_draws <- if (generate) {
    array(0, dim = c(di$K_w, 0L, 0L))
  } else {
    get_halton_normals(di$S, n_blocks, di$K_w)
  }
  args <- c(args, list(
    rc_dist = fit$rc_dist, rc_correlation = fit$rc_correlation,
    rc_mean = fit$rc_mean, use_asc = fit$use_asc,
    include_outside_option = fit$include_outside_option
  ))
  if (generate) {
    args <- c(args, list(
      gen_seed = as.integer(di$seed),
      gen_scramble = if (identical(di$scramble, "none")) 0L else 1L,
      gen_S = as.integer(di$S)
    ))
  }
  if (panel) args$Ti <- as.integer(d$Ti)
  fn <- switch(kernel,
    gradient = mxl_loglik_gradient_parallel,
    hessian  = mxl_hessian_parallel,
    bhhh     = mxl_bhhh_parallel,
    tastes   = mxl_conditional_tastes_parallel,
    stop("unknown kernel: ", kernel)
  )
  do.call(fn, args)
}

# A store-mode fit's stored data as an mxlp_oracle() fixture
# (helper-mxl-panel.R), with the fit's estimation draws.
pnl_oracle_fixture <- function(fit) {
  d <- fit$data
  n_units <- if (is.null(d$Ti)) length(d$M) else length(d$Ti)
  list(
    X = d$X, W = d$W, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
    M = d$M, weights = d$weights, Ti = d$Ti,
    eta = get_halton_normals(fit$draws_info$S, n_units, fit$draws_info$K_w),
    rc_dist = fit$rc_dist, rc_correlation = fit$rc_correlation,
    rc_mean = fit$rc_mean, use_asc = fit$use_asc,
    include_outside_option = fit$include_outside_option,
    K_w = ncol(d$W), J = max(fit$alt_mapping$alt_int)
  )
}

pnl_n_persons <- 90L

# Base panel: 90 decision makers, 3-5 situations each, J = 3, S = 30, with
# `grp` clustering decision makers in threes.
.pnl_base <- pnl_memo(function() {
  dt <- pnl_data(n_persons = pnl_n_persons, T_range = 3:5, seed = 11L)
  dt[, grp := (person - 1L) %/% 3L + 1L]
  list(dt = dt, fit = pnl_fit(dt, person_col = "person"))
})

# The same data fitted as a cross-section.
.pnl_cs <- pnl_memo(function() pnl_fit(.pnl_base()$dt))

# Every se_method in both draw modes, at the base panel estimate.
.pnl_se_fits <- pnl_memo(function() {
  b <- .pnl_base()
  theta <- unname(coef(b$fit))
  fits <- list(store = list(), generate = list())
  for (draws in names(fits)) {
    for (m in names(pnl_se_types)) {
      fits[[draws]][[m]] <- pnl_fit(
        b$dt, person_col = "person", draws = draws,
        seed = if (draws == "generate") 20260922L,
        se_method = m, cluster_col = if (m == "cluster") "grp",
        theta_init = theta, optimizer = pnl_at
      )
    }
  }
  fits
})

# Brute-force panel likelihood and conditional tastes at the base estimate.
.pnl_oracle <- pnl_memo(function() {
  f <- .pnl_base()$fit
  mxlp_oracle(unname(coef(f)), pnl_oracle_fixture(f))
})

# Recovery design: 400 decision makers x 8 situations, J = 4 inside
# alternatives plus physical outside-option rows, correlated 2 x 2 Sigma.
.pnl_recovery <- pnl_memo(function() {
  sim <- simulate_mxl_data(N = 400, T = 8, J = 4, seed = 42)
  fit <- suppressMessages(run_mxlogit(
    data = sim$data, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    outside_opt_label = 0L, include_outside_option = FALSE,
    rc_correlation = TRUE, person_col = "pid", S = 100L
  ))
  list(sim = sim, fit = fit, tastes = conditional_tastes(fit))
})

# --- 1. Data preparation -----------------------------------------------------

test_that("prepare_mxl_data(person_col) orders situations by decision maker", {
  T_p <- c(3L, 1L, 4L, 2L, 3L, 1L, 4L, 2L)
  dt <- pnl_data(T_p = T_p, seed = 3L)
  # Varying choice sets: drop an unchosen alternative from some situations.
  dt <- dt[!(alt == 3L & choice == 0L & id %% 3L == 0L)]
  prep <- pnl_prepare(dt, person_col = "person")

  ref <- dt[order(person, id, alt)]
  sit <- unique(ref[, .(person, id)])
  # The fixture interleaves decision makers in the id sequence, so the
  # prepared order is not the id order.
  expect_false(identical(sit$id, sort(sit$id)))

  expect_equal(prep$situation_ids, sit$id)
  expect_equal(prep$Ti, T_p)
  expect_equal(prep$person_ids, seq_along(T_p))
  expect_identical(prep$data_spec$person_col, "person")
  expect_equal(prep$N, nrow(sit))
  expect_equal(sum(prep$Ti), prep$N)
  expect_equal(unname(prep$X[, "x1"]), ref$x1)
  expect_equal(unname(prep$W), unname(as.matrix(ref[, .(w1, w2)])))
  expect_equal(prep$alt_idx, ref$alt)
  expect_equal(prep$M, ref[, .N, by = .(person, id)]$N)
  expect_equal(prep$choice_idx,
               ref[, which(choice == 1L), by = .(person, id)]$V1)
})

test_that("prepare_mxl_data() without person_col is the cross-sectional preparation", {
  dt <- pnl_data(T_p = c(3L, 1L, 4L, 2L, 3L, 1L, 4L, 2L), seed = 3L)
  p_default <- pnl_prepare(dt)
  p_null <- pnl_prepare(dt, person_col = NULL)
  expect_identical(p_null, p_default)
  expect_null(p_default$Ti)
  expect_null(p_default$person_ids)
  expect_null(p_default$data_spec$person_col)

  # Situations stay in id order.
  ref <- dt[order(id, alt)]
  expect_equal(p_default$situation_ids, sort(unique(dt$id)))
  expect_equal(unname(p_default$X[, "x1"]), ref$x1)
  expect_true(all(c("X", "W", "alt_idx", "choice_idx", "M", "N", "weights",
                    "cluster", "situation_ids", "include_outside_option",
                    "rc_correlation", "alt_mapping", "dropped_cols",
                    "data_spec") %in% names(p_default)))

  # The panel preparation holds the same situations in another order.
  p_panel <- pnl_prepare(dt, person_col = "person")
  expect_equal(p_panel$alt_mapping, p_default$alt_mapping)
  expect_setequal(p_panel$situation_ids, p_default$situation_ids)
})

test_that("prepare_mxl_data(person_col) rejects situation ids shared by decision makers", {
  dt <- pnl_data(T_p = rep(c(3L, 2L, 4L), 4L), seed = 5L)

  # The common mistake: situation ids that restart within each decision
  # maker (1, 2, ...), so every id belongs to several decision makers.
  restarted <- data.table::copy(dt)
  restarted[, id := match(id, sort(unique(id))), by = person]
  expect_error(pnl_prepare(restarted, person_col = "person"), "unique")

  # One situation whose rows are split between two decision makers (it
  # still has exactly one chosen row).
  dt_split <- data.table::copy(dt)
  target <- dt_split[person == 2L, min(id)]
  row <- dt_split[, which(id == target & choice == 0L)[1L]]
  dt_split[row, person := 3L]
  expect_error(pnl_prepare(dt_split, person_col = "person"), "unique")

  expect_error(pnl_prepare(dt, person_col = "nope"), "Missing columns")
})

test_that("weights and cluster labels must be constant within a decision maker", {
  dt <- pnl_data(T_p = rep(c(3L, 2L, 4L), 4L), seed = 5L)
  prep <- pnl_prepare(dt, person_col = "person")
  person_of <- rep(prep$person_ids, prep$Ti)          # prepared order

  # Situation-level weights: constant within an id, not within a person.
  dt[, w_sit := 1 + id / max(id)]
  expect_error(pnl_prepare(dt, weights_col = "w_sit", person_col = "person"),
               pnl_rx_within)
  expect_error(pnl_fit(dt, person_col = "person", weights_col = "w_sit"),
               pnl_rx_within)
  # The cross-section still accepts them.
  expect_no_error(pnl_prepare(dt, weights_col = "w_sit"))

  # Decision-maker weights by column, aligned by id.
  dt[, w_dm := 0.5 + person %% 3L]
  w_expected <- 0.5 + person_of %% 3L
  expect_equal(
    pnl_prepare(dt, weights_col = "w_dm", person_col = "person")$weights,
    w_expected
  )
  # By position they are ambiguous here: pnl_data() interleaves ids across
  # decision makers, so the prepared (person-first) order is not id order and
  # a vector built in id order would be silently misaligned.
  expect_error(pnl_prepare(dt, weights = w_expected, person_col = "person"),
               "weights_col")
  # With ids that are contiguous within decision makers the two orders
  # coincide, and positional weights are accepted and checked.
  dt_c <- data.table::copy(dt)
  key <- unique(dt_c[, .(person, id)])[order(person, id)]
  key[, id_c := .I]
  dt_c[key, id := i.id_c, on = .(person, id)]
  prep_c <- pnl_prepare(dt_c, person_col = "person")
  w_c <- 0.5 + rep(prep_c$person_ids, prep_c$Ti) %% 3L
  expect_equal(
    pnl_prepare(dt_c, weights = w_c, person_col = "person")$weights, w_c
  )
  expect_error(pnl_prepare(dt_c, weights = seq_len(prep_c$N) / prep_c$N,
                           person_col = "person"),
               pnl_rx_within)

  # Clusters must nest decision makers.
  dt[, cl_sit := id]
  expect_error(pnl_prepare(dt, cluster_col = "cl_sit", person_col = "person"),
               pnl_rx_cluster)
  expect_error(pnl_fit(dt, person_col = "person", cluster_col = "cl_sit"),
               pnl_rx_cluster)
  dt[, cl_grp := (person - 1L) %/% 3L + 1L]
  expect_equal(
    pnl_prepare(dt, cluster_col = "cl_grp", person_col = "person")$cluster,
    (person_of - 1L) %/% 3L + 1L
  )
})

test_that("WESML provenance cannot be combined with person_col", {
  dt <- pnl_data(T_p = rep(c(3L, 2L, 4L), 4L), seed = 5L)
  strata <- sort(unique(dt[choice == 1L, alt]))
  Q <- stats::setNames(rep(1 / length(strata), length(strata)), strata)
  dt_w <- wesml_weights(dt, "id", "alt", "choice", Q = Q, attach = TRUE)
  expect_false(is.null(attr(dt_w, "choice_sampling")))

  rx <- "WESML|choice-based"
  expect_error(pnl_prepare(dt_w, weights_col = ".wesml_weight",
                           person_col = "person"), rx)
  expect_error(pnl_prepare(dt_w, person_col = "person"), rx)
  # run_mxlogit() adopts the recorded weight column, then must refuse.
  expect_error(pnl_fit(dt_w, person_col = "person", S = 5L), rx)
  # The cross-sectional WESML preparation is unaffected.
  expect_no_error(pnl_prepare(dt_w, weights_col = ".wesml_weight"))
})

test_that("the advanced workflow takes one draw slice per decision maker", {
  b <- .pnl_base()
  f <- b$fit
  prep <- pnl_prepare(b$dt, person_col = "person")
  K_w <- ncol(prep$W)
  eta <- get_halton_normals(30L, length(prep$Ti), K_w)

  # person_col belongs to the convenience workflow.
  expect_error(
    run_mxlogit(input_data = prep, eta_draws = eta, S = 30L,
                person_col = "person"),
    "convenience"
  )

  # input_data from prepare_mxl_data(person_col = ) gives the panel fit.
  fit_adv <- suppressMessages(run_mxlogit(
    input_data = prep, eta_draws = eta, S = 30L, rc_mean = TRUE,
    theta_init = unname(coef(f)), optimizer = pnl_at
  ))
  expect_equal(fit_adv$n_persons, length(prep$Ti))
  expect_equal(fit_adv$data$Ti, prep$Ti)
  expect_equal(fit_adv$loglik, f$loglik, tolerance = 1e-10)
  expect_equal(vcov(fit_adv), vcov(f), tolerance = 1e-8)

  # A cube with one slice per choice situation is the cross-sectional cube.
  eta_situations <- get_halton_normals(30L, prep$N, K_w)
  expect_error(
    suppressMessages(run_mxlogit(input_data = prep, eta_draws = eta_situations,
                                 S = 30L, rc_mean = TRUE)),
    "eta_draws"
  )
})

# --- 2. Standard errors ------------------------------------------------------

test_that("every se_method gives finite standard errors on a panel fit", {
  fits <- .pnl_se_fits()
  n_sit <- nobs(.pnl_base()$fit)
  for (draws in names(fits)) {
    for (m in names(pnl_se_types)) {
      f <- fits[[draws]][[m]]
      what <- sprintf("[draws = %s, se_method = %s]", draws, m)
      expect_identical(f$se_method, m, label = paste(what, "se_method"))
      expect_identical(f$draws_info$mode, draws)
      expect_equal(dim(f$vcov), c(f$n_params, f$n_params))
      expect_true(all(is.finite(f$se)) && all(f$se > 0),
                  label = paste(what, "finite positive SEs"))
      expect_equal(f$n_persons, pnl_n_persons)
      expect_equal(nobs(f), n_sit)
    }
  }
})

# --- 3. Post-hoc variance ----------------------------------------------------

for (draws in c("store", "generate")) {
  test_that(sprintf("post-hoc vcov() reproduces the fit-time variance (draws = \"%s\")", draws), {
    fits <- .pnl_se_fits()[[draws]]
    for (m in names(pnl_se_types)) {
      f <- fits[[m]]
      what <- sprintf("[draws = %s, se_method = %s]", draws, m)
      expect_equal(vcov(f, type = pnl_se_types[[m]]), f$vcov,
                   tolerance = 1e-8, label = paste(what, "vcov(type = )"))
      lazy <- f
      lazy$vcov <- NULL
      lazy$se <- NULL
      lazy <- choicer:::ensure_vcov(lazy)
      expect_equal(lazy$vcov, f$vcov, tolerance = 1e-8,
                   label = paste(what, "ensure_vcov()"))
      expect_equal(lazy$se, f$se, tolerance = 1e-8,
                   label = paste(what, "ensure_vcov() SEs"))
    }

    # Every estimator can be recomputed from a fit made with another one.
    f <- fits$hessian
    expect_equal(vcov(f, type = "bhhh"), fits$bhhh$vcov, tolerance = 1e-8)
    expect_equal(vcov(f, type = "robust"), fits$sandwich$vcov,
                 tolerance = 1e-8)
    expect_equal(vcov(f, type = "cluster",
                      cluster = pnl_labels(.pnl_base()$dt, "grp")),
                 fits$cluster$vcov, tolerance = 1e-8)

    # The fit-time information matrices are the panel ones: they match the
    # kernels called with Ti, and not the cross-sectional kernels.
    expect_equal(unname(f$vcov), unname(solve(pnl_kernel("hessian", f))),
                 tolerance = 1e-8)
    expect_false(isTRUE(all.equal(
      unname(f$vcov), unname(solve(pnl_kernel("hessian", f, panel = FALSE))),
      tolerance = 1e-3
    )))
    expect_equal(unname(fits$bhhh$vcov),
                 unname(solve(pnl_kernel("bhhh", fits$bhhh))),
                 tolerance = 1e-8)
  })
}

# --- 4. Equivalences ---------------------------------------------------------

test_that("a panel of single-situation decision makers reproduces the cross-sectional fit", {
  # One thread makes the reductions deterministic, so both optimizer runs
  # see the same objective.
  set_num_threads(1L)
  on.exit(set_num_threads(2L), add = TRUE)
  dt <- data.table::copy(.pnl_base()$dt)
  dt[, pid1 := id]
  f_cs <- pnl_fit(dt)
  f_p1 <- pnl_fit(dt, person_col = "pid1")

  expect_equal(f_p1$n_persons, nobs(f_cs))
  expect_equal(f_p1$data$Ti, rep(1L, nobs(f_cs)))
  expect_equal(f_p1$data$person_ids, f_cs$data$situation_ids)
  expect_equal(coef(f_p1), coef(f_cs), tolerance = 1e-8)
  expect_equal(f_p1$loglik, f_cs$loglik, tolerance = 1e-8)
  expect_equal(vcov(f_p1), vcov(f_cs), tolerance = 1e-8)
  expect_equal(vcov(f_p1, type = "robust"), vcov(f_cs, type = "robust"),
               tolerance = 1e-8)
  ct_p1 <- conditional_tastes(f_p1)
  ct_cs <- conditional_tastes(f_cs)
  expect_equal(ct_p1$mean, ct_cs$mean, tolerance = 1e-10)
  expect_equal(ct_p1$sd, ct_cs$sd, tolerance = 1e-10)
})

test_that("clustering a panel fit on its decision makers gives the robust variance", {
  # One thread: the two calls below then share a bit-identical bread.
  set_num_threads(1L)
  on.exit(set_num_threads(2L), add = TRUE)
  b <- .pnl_base()
  f <- b$fit
  # Each decision maker is already one score row.
  expect_equal(vcov(f, type = "cluster", cluster = pnl_labels(b$dt, "person")),
               vcov(f, type = "robust"), tolerance = 1e-10)
  # Fit-time clustering on the decision maker is the same estimator.
  f_cl <- pnl_fit(b$dt, person_col = "person", cluster_col = "person",
                  theta_init = unname(coef(f)), optimizer = pnl_at)
  expect_identical(f_cl$se_method, "cluster")
  expect_equal(vcov(f_cl), vcov(f, type = "robust"), tolerance = 1e-8)
})

test_that("post-hoc cluster labels that vary within a decision maker are rejected", {
  f <- .pnl_base()$fit
  ids <- f$data$situation_ids
  # Every situation its own cluster: valid for a cross-section, not here.
  by_situation <- stats::setNames(seq_along(ids), ids)
  expect_error(vcov(f, type = "cluster", cluster = by_situation),
               pnl_rx_cluster)
})

test_that("scale_vars is invariant on a panel fit", {
  skip_on_cran()
  skip_on_ci()
  dt <- pnl_data(n_persons = 250L, T_range = 4L, seed = 23L)
  # Put the covariates on scales far from 1 so the scaling does something.
  dt[, `:=`(x1 = 10 * x1, w2 = 5 * w2)]
  ctrl <- list(xtol_rel = 1e-12, maxeval = 5000L)
  f_none <- pnl_fit(dt, person_col = "person", rc_correlation = TRUE,
                    S = 40L, control = ctrl)
  f_sd <- pnl_fit(dt, person_col = "person", rc_correlation = TRUE,
                  S = 40L, control = ctrl, scale_vars = "sd")

  expect_true(f_none$convergence > 0)
  expect_true(f_sd$convergence > 0)
  expect_equal(f_sd$n_persons, 250L)
  expect_lt(max(abs(coef(f_sd) - coef(f_none))), 1e-5)
  expect_lt(max(abs(f_sd$se - f_none$se)), 1e-5)
  expect_lt(abs(f_sd$loglik - f_none$loglik), 1e-8)
  expect_lt(max(abs(vcov(f_sd, type = "robust") -
                    vcov(f_none, type = "robust"))), 1e-5)
})

# --- 5. The panel likelihood is the one maximized ----------------------------

test_that("run_mxlogit(person_col) maximizes the panel, not the cross-sectional, likelihood", {
  b <- .pnl_base()
  f <- b$fit
  f_cs <- .pnl_cs()

  # The stored log-likelihood is the panel objective at the estimate, which
  # is a stationary point of it...
  panel <- pnl_kernel("gradient", f)
  expect_equal(f$loglik, -panel$objective, tolerance = 1e-10)
  expect_lt(max(abs(panel$gradient)), 1e-3)
  # ... and equals the brute-force Revelt-Train likelihood, one draw block
  # per decision maker.
  expect_equal(f$loglik, .pnl_oracle()$loglik, tolerance = 1e-8)

  # The cross-sectional likelihood is a different function of the same data.
  cross <- -pnl_kernel("gradient", f, panel = FALSE)$objective
  expect_gt(abs(f$loglik - cross), 1)
  expect_gt(abs(f$loglik - f_cs$loglik), 1e-3)
  expect_false(isTRUE(all.equal(coef(f), coef(f_cs), tolerance = 1e-4)))
  expect_null(f_cs$n_persons)
  expect_null(f_cs$data$Ti)
  expect_null(f_cs$data$person_ids)
})

# --- 6. Post-estimation ------------------------------------------------------

test_that("a panel fit counts choice situations and decision makers", {
  b <- .pnl_base()
  f <- b$fit
  n_sit <- data.table::uniqueN(b$dt$id)
  expect_equal(nobs(f), n_sit)
  expect_equal(attr(logLik(f), "nobs"), n_sit)
  expect_equal(f$draws_info$N, n_sit)
  expect_equal(f$n_persons, pnl_n_persons)
  expect_length(f$data$Ti, pnl_n_persons)
  expect_equal(sum(f$data$Ti), n_sit)
  expect_equal(f$data$person_ids, seq_len(pnl_n_persons))
  expect_identical(f$data_spec$person_col, "person")
  expect_true(is.finite(AIC(f)) && is.finite(BIC(f)))
})

test_that("print() and summary() report the decision makers of a panel fit", {
  b <- .pnl_base()
  line <- sprintf("Respondents: %d", pnl_n_persons)
  out <- capture.output(print(b$fit))
  expect_true(any(grepl(line, out, fixed = TRUE)))
  out_summary <- capture.output(print(summary(b$fit)))
  expect_true(any(grepl(line, out_summary, fixed = TRUE)))

  # Cross-sectional output has no respondent line.
  f_cs <- pnl_fit(b$dt, theta_init = unname(coef(b$fit)), optimizer = pnl_at)
  expect_false(any(grepl("Respondents", capture.output(print(f_cs)))))
  expect_false(any(grepl("Respondents",
                         capture.output(print(summary(f_cs))))))
})

test_that("post-estimation methods work on a panel fit", {
  b <- .pnl_base()
  f <- b$fit
  n_sit <- nobs(f)
  J <- nrow(f$alt_mapping)

  p <- predict(f)
  expect_length(p$choice_prob, nrow(b$dt))
  expect_true(all(is.finite(p$choice_prob)))
  situation <- rep(seq_len(n_sit), f$data$M)
  expect_equal(as.numeric(rowsum(p$choice_prob, situation)), rep(1, n_sit),
               tolerance = 1e-10)
  # Round trip: the (shuffled) estimation data as newdata reproduce the
  # in-sample predictions value for value, in the same order.
  expect_equal(predict(f, newdata = b$dt), p, tolerance = 1e-12)

  shares <- predict(f, type = "shares")
  expect_length(shares, J)
  expect_equal(sum(shares), 1, tolerance = 1e-10)
  expect_equal(predict(f, type = "shares", newdata = b$dt), shares,
               tolerance = 1e-12)

  el <- elasticities(f, elast_var = "x1")
  expect_equal(dim(el), c(J, J))
  expect_true(all(is.finite(el)))
  expect_true(all(is.finite(elasticities(f, elast_var = "w1",
                                         is_random_coef = TRUE))))

  dr <- diversion_ratios(f, wrt_var = "x1")
  expect_equal(dim(dr), c(J, J))
  expect_true(all(is.finite(dr)))

  # Share inversion at the fitted shares returns the fitted ASCs.
  delta <- blp(f, target_shares = shares)
  expect_equal(as.numeric(delta) - as.numeric(delta)[1L],
               unname(c(0, coef(f)[f$param_map$asc])), tolerance = 1e-6)

  g <- gof(f)
  expect_s3_class(g, "choicer_gof")
  expect_equal(g$nobs, n_sit)
  expect_true(is.finite(g$mcfadden_r2) && g$mcfadden_r2 > 0)
  expect_true(g$hit_rate >= 0 && g$hit_rate <= 1)

  w <- wtp(f, price_var = "x1")
  expect_s3_class(w, "choicer_wtp")
  expect_equal(nrow(w), 2L)
  expect_true(all(is.finite(w$Estimate)) && all(is.finite(w$Std_Error)))

  ls <- logsum(f)
  expect_length(ls, n_sit)
  expect_true(all(is.finite(ls)))

  cs <- consumer_surplus(f, price_var = "x1")
  expect_s3_class(cs, "choicer_cs")
  expect_length(cs$cs, n_sit)
  expect_true(is.finite(cs$mean_cs))
})

# --- 7. Conditional tastes ---------------------------------------------------

test_that("conditional_tastes() returns one taste summary per decision maker", {
  f <- .pnl_base()$fit
  ct <- conditional_tastes(f)
  expect_s3_class(ct, "choicer_tastes")

  dims <- list(c("w1", "w2"), as.character(f$data$person_ids))
  expect_equal(dim(ct$mean), c(2L, pnl_n_persons))
  expect_equal(dim(ct$sd), c(2L, pnl_n_persons))
  expect_identical(dimnames(ct$mean), dims)
  expect_identical(dimnames(ct$sd), dims)
  expect_true(all(is.finite(ct$mean)))
  expect_true(all(is.finite(ct$sd)) && all(ct$sd >= 0))
  expect_equal(ct$n_units, pnl_n_persons)
  expect_equal(ct$S, 30L)
  expect_equal(ct$weights, rep(1, pnl_n_persons))
  expect_type(ct$unit, "character")
  expect_length(ct$unit, 1L)

  pop <- ct$population
  expect_s3_class(pop, "data.frame")
  expect_equal(nrow(pop), 2L)
  expect_true(all(c("model_mean", "model_sd", "mean_cond", "sd_cond",
                    "revealed", "residual") %in% names(pop)))

  out <- capture.output(res <- print(ct))
  expect_identical(res, ct)
  expect_true(any(grepl("w1", out, fixed = TRUE)))
  expect_true(any(grepl("w2", out, fixed = TRUE)))

  # Generate mode draws decision maker u's tastes from Halton block u too.
  f_gen <- .pnl_se_fits()$generate$hessian
  ct_gen <- conditional_tastes(f_gen)
  expect_equal(dim(ct_gen$mean), c(2L, pnl_n_persons))
  ref <- pnl_kernel("tastes", f_gen)
  mxlp_expect_close(unname(ct_gen$mean), ref$mean, 1e-12,
                    "generate-mode conditional mean vs kernel with Ti")
  mxlp_expect_close(unname(ct_gen$sd), ref$sd, 1e-12,
                    "generate-mode conditional SD vs kernel with Ti")
})

test_that("conditional tastes equal the brute-force oracle at the estimate", {
  f <- .pnl_base()$fit
  ct <- conditional_tastes(f)
  orc <- .pnl_oracle()
  mxlp_expect_close(unname(ct$mean), orc$mean, 1e-10,
                    "panel conditional mean vs oracle")
  mxlp_expect_close(unname(ct$sd), orc$sd, 1e-10,
                    "panel conditional SD vs oracle")

  # Population table (normal random coefficients, uniform weights).
  pop <- ct$population
  var_model <- diag(f$sigma)
  expect_equal(pop$model_mean, unname(coef(f)[f$param_map$mu]),
               tolerance = 1e-10)
  expect_equal(pop$model_sd, unname(sqrt(var_model)), tolerance = 1e-10)
  expect_equal(pop$mean_cond, rowMeans(orc$mean), tolerance = 1e-10)
  expect_equal(pop$residual, unname(rowMeans(orc$sd^2) / var_model),
               tolerance = 1e-8)
  # Var(E[beta | y]) with either the 1/U or the 1/(U - 1) convention.
  U <- ncol(orc$mean)
  v_pop <- rowMeans((orc$mean - rowMeans(orc$mean))^2)
  v_smp <- v_pop * U / (U - 1)
  v_rev <- pop$revealed * unname(var_model)
  expect_true(all(v_rev >= v_pop * (1 - 1e-8) & v_rev <= v_smp * (1 + 1e-8)))
  expect_true(all(pop$sd_cond >= sqrt(v_pop) * (1 - 1e-8) &
                  pop$sd_cond <= sqrt(v_smp) * (1 + 1e-8)))
})

test_that("conditional_tastes() on a cross-sectional fit has one unit per situation", {
  f_cs <- .pnl_cs()
  ct <- conditional_tastes(f_cs)
  n_sit <- nobs(f_cs)
  expect_s3_class(ct, "choicer_tastes")
  expect_equal(dim(ct$mean), c(2L, n_sit))
  expect_identical(colnames(ct$mean), as.character(f_cs$data$situation_ids))
  expect_equal(ct$n_units, n_sit)
  expect_false(identical(ct$unit, conditional_tastes(.pnl_base()$fit)$unit))

  orc <- mxlp_oracle(unname(coef(f_cs)), pnl_oracle_fixture(f_cs))
  mxlp_expect_close(unname(ct$mean), orc$mean, 1e-10,
                    "cross-section conditional mean vs oracle")
  mxlp_expect_close(unname(ct$sd), orc$sd, 1e-10,
                    "cross-section conditional SD vs oracle")
})

test_that("conditional_tastes() needs the stored data", {
  b <- .pnl_base()
  f_slim <- pnl_fit(b$dt, person_col = "person", keep_data = FALSE,
                    theta_init = unname(coef(b$fit)), optimizer = pnl_at)
  expect_error(conditional_tastes(f_slim), "keep_data")
})

# --- 8. Simulator ------------------------------------------------------------

test_that("simulate_mxl_data(T = 1) keeps the cross-sectional layout", {
  sim <- simulate_mxl_data(N = 40, J = 3, seed = 1)
  expect_identical(names(sim$data),
                   c("id", "alt", "w1", "w2", "x1", "x2", "choice"))
  expect_false("pid" %in% names(sim$data))
  expect_identical(data.table::key(sim$data), c("id", "alt"))
  expect_equal(sort(unique(sim$data$id)), 1:40)
  expect_identical(simulate_mxl_data(N = 40, J = 3, seed = 1, T = 1L), sim)
  expect_equal(sim$settings$T, 1)
  expect_equal(dim(sim$true_params$gamma_i), c(2L, 40L))
  expect_true(all(c("beta", "delta", "Sigma", "L_params", "mu", "rc_dist",
                    "rc_correlation") %in% names(sim$true_params)))

  sim_inside <- simulate_mxl_data(N = 40, J = 3, seed = 1,
                                  outside_option = FALSE)
  expect_identical(names(sim_inside$data),
                   c("id", "w1", "w2", "x1", "x2", "alt", "choice"))
})

test_that("simulate_mxl_data(T > 1) simulates T situations per decision maker", {
  N <- 30L
  T_ <- 4L
  sim <- simulate_mxl_data(N = N, J = 3, seed = 2, T = T_)
  dt <- sim$data
  expect_true("pid" %in% names(dt))
  expect_false(anyNA(dt$pid))

  # Ids 1..N*T; pid constant within an id, in runs of length T.
  sit <- unique(dt[, .(id, pid)])[order(id)]
  expect_equal(sit$id, seq_len(N * T_))
  expect_equal(sit$pid, rep(seq_len(N), each = T_))
  expect_equal(dt[, sum(choice), by = id]$V1, rep(1L, N * T_))

  # Outside-option rows carry their situation's pid and zero covariates.
  out <- dt[alt == 0L]
  expect_equal(nrow(out), N * T_)
  expect_equal(out$pid, (out$id - 1L) %/% T_ + 1L)
  expect_true(all(as.matrix(out[, .(w1, w2, x1, x2)]) == 0))

  expect_equal(dim(sim$true_params$gamma_i), c(2L, N))
  expect_equal(sim$settings$T, T_)
  expect_equal(sim$settings$N, N)

  prep <- prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "x2"),
                           c("w1", "w2"), outside_opt_label = 0L,
                           person_col = "pid")
  expect_equal(prep$Ti, rep(T_, N))
  expect_equal(prep$person_ids, seq_len(N))
})

test_that("simulate_mxl_data() draws one taste vector per decision maker", {
  # The tastes are the first draws after set.seed(): L %*% z, one column of
  # z per decision maker, whatever T is. This pins both the meaning of
  # true_params$gamma_i and the unchanged random-number stream at T = 1.
  N <- 25L
  Sigma <- matrix(c(1.0, 0.5, 0.5, 1.5), 2L)
  set.seed(9)
  gamma_expected <- t(chol(Sigma)) %*% matrix(stats::rnorm(N * 2L), 2L, N)
  for (T_ in c(1L, 4L)) {
    sim <- simulate_mxl_data(N = N, J = 3, seed = 9, T = T_)
    expect_equal(unname(sim$true_params$gamma_i), gamma_expected,
                 tolerance = 1e-12, label = sprintf("gamma_i (T = %d)", T_))
  }
})

test_that("simulate_mxl_data() rejects invalid T", {
  expect_error(simulate_mxl_data(N = 10, J = 3, T = 0L))
  expect_error(simulate_mxl_data(N = 10, J = 3, T = c(2L, 3L)))
})

# --- 9. Recovery -------------------------------------------------------------

test_that("a panel fit recovers the simulated parameters and tastes", {
  skip_on_cran()
  r <- .pnl_recovery()
  fit <- r$fit
  expect_true(fit$convergence %in% 1:4)
  expect_equal(fit$n_persons, 400L)
  expect_equal(nobs(fit), 3200L)

  rt <- recovery_table(fit, r$sim)
  expect_setequal(unique(rt$group), c("beta", "sigma", "asc"))
  expect_true(all(is.finite(rt$se)))
  missed <- rt$parameter[!(abs(rt$estimate - rt$true) < pmax(3 * rt$se, 0.15))]
  expect_identical(missed, character(0))

  # The conditional means track the realized tastes (aligned by pid).
  ct <- r$tastes
  gamma <- r$sim$true_params$gamma_i[, as.integer(colnames(ct$mean)),
                                     drop = FALSE]
  for (k in seq_len(nrow(ct$mean))) {
    expect_gt(stats::cor(ct$mean[k, ], gamma[k, ]), 0.6,
              label = sprintf("cor(E[beta | y], gamma) for %s",
                              rownames(ct$mean)[k]))
  }
})

test_that("revealed and residual taste variance add up to the model variance", {
  skip_on_cran()
  # Law of total variance: Var(E[beta | y]) + E[Var(beta | y)] = Var(beta).
  pop <- .pnl_recovery()$tastes$population
  expect_true(all(abs(pop$revealed + pop$residual - 1) < 0.15))
})

# --- 10. Review follow-ups ---------------------------------------------------

# Person-constant, non-uniform decision-maker weights at the base estimate.
.pnl_wt_fits <- pnl_memo(function() {
  b <- .pnl_base()
  dt <- data.table::copy(b$dt)
  dt[, w_dm := 0.5 + (person %% 4L) / 2]
  theta <- unname(coef(b$fit))
  fits <- list()
  for (m in names(pnl_se_types)) {
    # Non-uniform weights make hessian/bhhh fits warn that they are not WESML
    # corrections; that advice is beside the point here.
    fits[[m]] <- suppressWarnings(pnl_fit(
      dt, person_col = "person", weights_col = "w_dm", se_method = m,
      cluster_col = if (m == "cluster") "grp",
      theta_init = theta, optimizer = pnl_at
    ))
  }
  list(dt = dt, fits = fits)
})

test_that("decision-maker weights enter scores, variances and tastes per unit", {
  wf <- .pnl_wt_fits()
  fits <- wf$fits
  for (m in names(pnl_se_types)) {
    f <- fits[[m]]
    what <- sprintf("[weighted, se_method = %s]", m)
    expect_gt(length(unique(f$data$weights)), 1L)
    expect_equal(vcov(f, type = pnl_se_types[[m]]), f$vcov, tolerance = 1e-8,
                 label = paste(what, "vcov(type = )"))
    lazy <- f
    lazy$vcov <- NULL
    lazy$se <- NULL
    expect_equal(choicer:::ensure_vcov(lazy)$vcov, f$vcov, tolerance = 1e-8,
                 label = paste(what, "ensure_vcov()"))
  }

  # The C++ unit collapse (sandwich meat from the BHHH kernel with squared
  # weights) equals the R collapse (clustering the unit scores on persons).
  set_num_threads(1L)
  on.exit(set_num_threads(2L), add = TRUE)
  f <- fits$sandwich
  expect_equal(vcov(f, type = "cluster", cluster = pnl_labels(wf$dt, "person")),
               f$vcov, tolerance = 1e-10)

  # The population table weights units by their decision-maker weight.
  ct <- conditional_tastes(f)
  w_u <- f$data$weights[cumsum(c(1L, head(f$data$Ti, -1L)))]
  expect_equal(unname(ct$weights), w_u)
  expect_equal(ct$population$mean_cond,
               unname(drop(ct$mean %*% (w_u / sum(w_u)))), tolerance = 1e-12)
})

test_that("population moments of shifted log-normal coefficients match the closed forms", {
  dt <- .pnl_base()$dt
  for (rc_mean in c(TRUE, FALSE)) {
    theta <- c(-0.5, if (rc_mean) c(0.2, -0.3), log(c(0.6, 0.8)), 0.1, -0.2)
    f <- suppressMessages(run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = "x1", random_var_cols = c("w1", "w2"),
      person_col = "person", rc_dist = c(1L, 0L), rc_mean = rc_mean,
      S = 30L, theta_init = theta, optimizer = pnl_at
    ))
    pop <- conditional_tastes(f)$population
    s2 <- unname(diag(f$sigma))
    mu <- if (rc_mean) unname(coef(f)[f$param_map$mu]) else c(0, 0)
    # w1 is shifted log-normal: exp(mu) 1{rc_mean} + exp(L eta).
    shift1 <- if (rc_mean) exp(mu[1]) else 0
    expect_equal(pop$model_mean[1], shift1 + exp(s2[1] / 2), tolerance = 1e-12)
    expect_equal(pop$model_sd[1], sqrt((exp(s2[1]) - 1) * exp(s2[1])),
                 tolerance = 1e-12)
    # w2 is normal: mu 1{rc_mean} + L eta.
    expect_equal(pop$model_mean[2], mu[2], tolerance = 1e-12)
    expect_equal(pop$model_sd[2], sqrt(s2[2]), tolerance = 1e-12)
  }
})

test_that("scale_vars leaves a panel fit's variances unchanged at a fixed theta", {
  b <- .pnl_base()
  theta <- unname(coef(b$fit))
  unscaled <- .pnl_se_fits()$store
  for (sv in c("sd", "mad")) {
    for (m in names(pnl_se_types)) {
      what <- sprintf("[scale_vars = %s, se_method = %s]", sv, m)
      f <- pnl_fit(b$dt, person_col = "person", scale_vars = sv,
                   se_method = m, cluster_col = if (m == "cluster") "grp",
                   theta_init = theta, optimizer = pnl_at)
      # The fit-time vcov, like the post-hoc one, is computed at the natural
      # estimates on the stored natural-scale data.
      expect_equal(vcov(f, type = pnl_se_types[[m]]), f$vcov,
                   tolerance = 1e-8, label = paste(what, "post hoc"))
      expect_equal(f$vcov, unscaled[[m]]$vcov, tolerance = 1e-8,
                   label = paste(what, "vs unscaled"))
    }
  }
})

test_that("unnamed post-hoc cluster labels are rejected when a panel reorders situations", {
  b <- .pnl_base()
  f <- b$fit
  ids <- f$data$situation_ids
  # pnl_data() interleaves ids across decision makers: prepared order is not
  # id order, so unnamed labels cannot be aligned safely.
  expect_false(identical(order(ids, method = "radix"), seq_along(ids)))
  named <- pnl_labels(b$dt, "grp")
  expect_error(vcov(f, type = "cluster", cluster = unname(named)), "[Nn]ame")
  expect_equal(vcov(f, type = "cluster", cluster = named),
               .pnl_se_fits()$store$cluster$vcov, tolerance = 1e-8)

  # A cross-sectional fit keeps the unnamed path (with its warning).
  f_cs <- .pnl_cs()
  cl_cs <- named[as.character(f_cs$data$situation_ids)]
  expect_warning(vcov(f_cs, type = "cluster", cluster = unname(cl_cs)),
                 "Unnamed")
})

test_that("newdata without the person column is flagged for panel fits", {
  b <- .pnl_base()
  f <- b$fit
  expect_message(predict(f, newdata = b$dt[, !"person"]), "person")
  expect_no_message(predict(f, newdata = b$dt))
})

test_that("the advanced workflow records the draws actually used", {
  b <- .pnl_base()
  prep <- pnl_prepare(b$dt, person_col = "person")
  eta40 <- get_halton_normals(40L, length(prep$Ti), ncol(prep$W))
  # S is left at its default (100): the cube's 40 draws must be recorded.
  f <- suppressMessages(run_mxlogit(
    input_data = prep, eta_draws = eta40, rc_mean = TRUE,
    theta_init = unname(coef(b$fit)), optimizer = pnl_at
  ))
  expect_equal(f$draws_info$S, 40L)
  lazy <- f
  lazy$vcov <- NULL
  lazy$se <- NULL
  expect_equal(choicer:::ensure_vcov(lazy)$vcov, f$vcov, tolerance = 1e-8)
})

test_that("conditional_tastes() explains which fits it supports", {
  fit_mnl <- suppressMessages(run_mnlogit(
    .pnl_base()$dt, "id", "alt", "choice", "x1",
    control = list(maxeval = 20L)
  ))
  expect_error(conditional_tastes(fit_mnl), "run_mxlogit")
})
