# Store-mode post-estimation forms the draws on the fly (gen_scramble = 2):
# the generator's identity-permutation uniforms are randtoolbox::halton()'s
# points, and R's qnorm() maps them as get_halton_normals() does, so the
# kernels' results equal the K_w x S x N cube's without the cube: per
# situation at any thread count, sums over situations at one thread. Bit for
# bit where randtoolbox and choicer are compiled with the same floating-point
# contraction (as on CRAN's platforms); otherwise some uniforms differ in
# their last bit, so the kernel tests compare against a cube of choicer's own
# uniforms. blp() keeps the draws it forms on the fly across its iterations,
# within its keep_draws_bytes. run_mxlogit() warns before it builds an
# estimation cube above .mxl_cube_budget() (1 GiB).

# --- Generator -------------------------------------------------------------------

test_that("the store-mode points are randtoolbox's", {
  # Uniforms at the start of the sequence, in a middle block and in the last
  # block a store-mode cube can reach (it ends at 2^31 - 1); randtoolbox's
  # normals are R's qnorm() of its uniforms
  same <- close <- TRUE
  for (K_w in c(1L, 2L, 13L, 128L)) {
    for (n0 in c(1, 123457, 2^31 - 8)) {
      # A double, as get_halton_normals() passes it: randtoolbox adds it to
      # its integer offset, which an integer would overflow at 2^31
      S <- if (n0 > 2^30) 8 else 9
      u <- choicer:::halton_fill_uniforms(n0, S, K_w, 0, 0L)
      ref <- t(matrix(randtoolbox::halton(S, K_w, normal = FALSE, start = n0),
                      S, K_w))
      z <- t(matrix(randtoolbox::halton(S, K_w, normal = TRUE, start = n0),
                    S, K_w))
      lab <- sprintf("K_w %d, n0 %g", K_w, n0)
      expect_equal(u, ref, tolerance = 1e-15, label = lab)
      expect_identical(qnorm(ref), z, label = lab)
      close <- close && max(abs(u - ref)) <= .Machine$double.eps
      same <- same && identical(u, ref)
    }
  }
  # Bit for bit when both packages fuse the radical inverse's multiply-adds
  # alike; otherwise some uniforms differ in their last bit
  skip_if(!same && close, paste("randtoolbox and choicer were compiled with",
                                "different floating-point contraction"))
  expect_true(same)
})

test_that("fill_uniforms() past the cube's reach: its uniforms and normals", {
  # Store-mode prediction reads the generator's uniforms past 2^31 - 1: per
  # index against the radical inverse (below 2^53, where R's doubles hold
  # every index), and through fill_block(), which is checked against a
  # per-index reference up to 2^64 - 1 (test-halton-generator.R)
  skip_if(R.version$arch %in% c("i386", "i686"), "x87 excess precision")
  primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
              61, 67, 71)
  for (n0 in c(2^31 - 3, 2^32 - 3, 2^53 - 9)) {
    for (K_w in c(1L, 3L, 20L)) {
      ref <- outer(seq_len(K_w), 0:7, function(k, s) {
        mapply(choicer:::halton_radical_inverse, n0 + s, primes[k])
      })
      expect_identical(choicer:::halton_fill_uniforms(n0, 8L, K_w, 0, 0L), ref,
                       label = sprintf("uniforms, K_w %d, n0 %g", K_w, n0))
    }
  }
  inv <- function(u) {
    array(vapply(u, choicer:::halton_inv_normal_cdf, 0), dim(u))
  }
  for (n0 in c(2^31 - 3, 2^32 - 3, 2^53 - 4, 1e18)) {
    for (K_w in c(1L, 3L, 20L)) {
      expect_identical(inv(choicer:::halton_fill_uniforms(n0, 8L, K_w, 0, 0L)),
                       choicer:::halton_fill_block(n0, 8L, K_w, 0, 0L),
                       label = sprintf("normals, K_w %d, n0 %g", K_w, n0))
    }
  }
})

# The store-mode cube of N situations from choicer's own uniforms and R's
# qnorm(): what the kernels form on the fly, and what get_halton_normals()
# builds where randtoolbox and choicer are compiled alike
mxc_points_cube <- function(S, N, K_w) {
  u <- choicer:::halton_fill_uniforms(1, S * N, K_w, 0, 0L)
  array(qnorm(u), dim = c(K_w, S, N))
}

# --- Kernels -------------------------------------------------------------------

mxc_data <- function(seed, N, J = 5L, K_w = 2L, S = 6L, ioo = FALSE,
                     use_asc = TRUE, w_type = "row") {
  set.seed(seed)
  M <- c(J, sample(2:J, N - 1L, replace = TRUE))  # every alternative appears
  alt_idx <- as.integer(unlist(lapply(M, function(m) sort(sample(J, m)))))
  n <- sum(M)
  K_x <- 2L
  W <- if (w_type == "row") matrix(rnorm(n * K_w), n, K_w) else
    matrix(rnorm(J * K_w), J, K_w)
  n_delta <- if (!use_asc) 0L else if (ioo) J else J - 1L
  L_size <- K_w * (K_w + 1L) / 2L
  list(theta = rnorm(K_x + K_w + L_size + n_delta, sd = 0.4),
       X = matrix(rnorm(n * K_x), n, K_x), W = W, alt_idx = alt_idx, M = M,
       weights = runif(N, 0.5, 2), rc_dist = c(0L, rep(1L, K_w - 1L)),
       eta = mxc_points_cube(S, N, K_w), K_w = K_w, S = S, N = N, J = J,
       ioo = ioo, use_asc = use_asc)
}

# Call kernel `k` on `d` with the whole cube, or with the store-mode points
# formed on the fly (points = TRUE)
mxc_call <- function(k, d, points = FALSE, var = 1L, random = FALSE) {
  eta <- d$eta
  gen <- list()
  if (points) {
    eta <- array(0, dim = c(d$K_w, 0L, 0L))
    gen <- list(gen_seed = 0L, gen_scramble = 2L, gen_S = d$S)
  }
  flags <- list(rc_correlation = TRUE, rc_mean = TRUE, use_asc = d$use_asc,
                include_outside_option = d$ioo)
  switch(k,
    pred = do.call(choicer:::mxl_predict,
                   c(list(d$theta, d$X, d$W, d$alt_idx, d$M, eta, d$rc_dist),
                     flags, gen)),
    logsum = do.call(choicer:::mxl_logsum,
                     c(list(d$theta, d$X, d$W, d$alt_idx, d$M, eta, d$rc_dist),
                       flags, gen)),
    shares = do.call(choicer:::mxl_predict_shares,
                     c(list(d$theta, d$X, d$W, d$alt_idx, d$M, d$weights, eta,
                            d$rc_dist), flags, gen)),
    elas = do.call(choicer:::mxl_elasticities_parallel,
                   c(list(d$theta, d$X, d$W, d$alt_idx, NULL, d$M, d$weights,
                          eta, d$rc_dist, var, random), flags, gen)),
    dr = do.call(choicer:::mxl_diversion_ratios_parallel,
                 c(list(d$theta, d$X, d$W, d$alt_idx, d$M, d$weights, eta,
                        d$rc_dist, var, random), flags, gen)),
    blp = {
      K_w <- d$K_w
      target <- as.numeric(mxc_call("shares", d))
      do.call(mxl_blp_contraction,
              c(list(rep(0, d$J), target, d$X, d$W, d$theta[1:2],
                     d$theta[2L + seq_len(K_w)],
                     d$theta[2L + K_w + seq_len(K_w * (K_w + 1L) / 2L)],
                     d$alt_idx, d$M, d$weights, eta, d$rc_dist,
                     rc_correlation = TRUE, rc_mean = TRUE,
                     include_outside_option = d$ioo), gen))
    })
}

test_that("store-mode points on the fly give the whole cube's results", {
  on.exit(set_num_threads(2L), add = TRUE)
  for (ioo in c(FALSE, TRUE)) {
    for (w_type in c("row", "alt")) {
      d <- mxc_data(100 + ioo, 23L, ioo = ioo, w_type = w_type)
      set_num_threads(1L)
      for (k in c("pred", "logsum", "shares", "elas", "dr", "blp")) {
        expect_identical(mxc_call(k, d, points = TRUE), mxc_call(k, d),
                         label = sprintf("%s, ioo %d, W %s", k, ioo, w_type))
      }
      expect_identical(mxc_call("elas", d, TRUE, 2L, TRUE),
                       mxc_call("elas", d, FALSE, 2L, TRUE))
      expect_identical(mxc_call("dr", d, TRUE, 2L, TRUE),
                       mxc_call("dr", d, FALSE, 2L, TRUE))
      # Per situation at two threads (R's qnorm() runs in the threads); sums
      # over situations to rounding
      whole_pred <- mxc_call("pred", d)
      whole_ls <- mxc_call("logsum", d)
      set_num_threads(2L)
      expect_identical(mxc_call("pred", d, points = TRUE), whole_pred)
      expect_identical(mxc_call("logsum", d, points = TRUE), whole_ls)
      expect_equal(mxc_call("shares", d, points = TRUE), mxc_call("shares", d),
                   tolerance = 1e-13)
    }
  }
})

test_that("store-mode points without ASCs, with one draw, and one coefficient", {
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  d <- mxc_data(110, 9L, use_asc = FALSE)
  for (k in c("pred", "logsum", "shares", "elas", "dr")) {
    expect_identical(mxc_call(k, d, TRUE), mxc_call(k, d))
  }
  d <- mxc_data(111, 9L, S = 1L)
  for (k in c("pred", "logsum", "shares")) {
    expect_identical(mxc_call(k, d, TRUE), mxc_call(k, d))
  }
  d <- mxc_data(112, 9L, K_w = 1L)
  for (k in c("pred", "logsum", "shares", "elas")) {
    expect_identical(mxc_call(k, d, TRUE), mxc_call(k, d))
  }
})

# --- R methods -------------------------------------------------------------

mxc_fit <- function(person = FALSE, draws = "store") {
  set.seed(140)
  N <- 60L; J <- 4L; T <- if (person) 3L else 1L
  dt <- data.table::data.table(
    id = rep(seq_len(N * T), each = J), alt = rep(seq_len(J), N * T),
    person = rep(seq_len(N), each = J * T))
  dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
  dt[, choice := 0L]
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  fit <- suppressMessages(suppressWarnings(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = c("w1", "w2"),
    person_col = if (person) "person" else NULL, S = 20L, draws = draws,
    seed = if (draws == "generate") 3L else NULL,
    control = list(maxeval = 30L))))
  # A counterfactual: the second alternative 0.5 dearer in x1
  cf <- data.table::copy(dt)[alt == 2L, x1 := x1 + 0.5]
  list(fit = fit, cf = cf)
}

mxc_methods <- function(fit, cf) {
  list(
    probs = predict(fit),
    shares = predict(fit, type = "shares"),
    cf_probs = predict(fit, newdata = cf),
    logsum = logsum(fit),
    cf_logsum = logsum(fit, newdata = cf),
    cs = consumer_surplus(fit, price_var = "x1"),
    cf_cs = consumer_surplus(fit, price_var = "x1", newdata = cf),
    elas = elasticities(fit, "w1", is_random_coef = TRUE),
    dr = diversion_ratios(fit, "x1"),
    gof = gof(fit))
}

test_that("store-mode methods form their draws on the fly, as the cube's", {
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  cube_draws <- function(draws_info, N, what) {
    list(eta_draws = mxc_points_cube(draws_info$S, N, draws_info$K_w),
         gen_seed = -1L, gen_scramble = 1L, gen_S = 0L)
  }
  for (person in c(FALSE, TRUE)) {
    fx <- mxc_fit(person)
    # A cube of the same draws. blp() inverts the counterfactual's shares, so
    # that it iterates (the fit's own shares stop it at its first check)
    with_mocked_bindings(
      {
        whole <- mxc_methods(fx$fit, fx$cf)
        target <- as.numeric(predict(fx$fit, newdata = fx$cf, type = "shares"))
        delta <- blp(fx$fit, target_shares = target)
      },
      .mxl_pred_draws = cube_draws)
    expect_output(start <- blp(fx$fit, target_shares = target, max_iter = 0L),
                  "Maximum iterations reached", fixed = TRUE)
    expect_false(identical(delta, start))
    # On the fly: no site may build the cube
    with_mocked_bindings(
      {
        expect_identical(mxc_methods(fx$fit, fx$cf), whole)
        expect_identical(blp(fx$fit, target_shares = target), delta)
      },
      get_halton_normals = function(...) stop("the whole cube was built"),
      .halton_cube = function(...) stop("the whole cube was built"))
  }
})

test_that("the draw arguments of store-mode post-estimation", {
  info <- list(mode = "store", S = 10L, K_w = 2L, N = 50L)
  gp <- choicer:::.mxl_pred_draws(info, 50L, "predict()")
  expect_identical(gp$gen_seed, 0L)
  expect_identical(gp$gen_scramble, 2L)
  expect_identical(gp$gen_S, 10L)
  expect_identical(dim(gp$eta_draws), c(2L, 0L, 0L))
  # No cube, so no limit of 2^31 - 1 points on the situations predicted
  expect_identical(choicer:::.mxl_pred_draws(info, 3e8, "predict()")$gen_S, 10L)
  # Past the generator's 128 primes, the cube
  wide <- choicer:::.mxl_pred_draws(list(mode = "store", S = 2L, K_w = 129L),
                                    3L, "predict()")
  expect_identical(wide$gen_seed, -1L)
  expect_identical(wide$eta_draws, get_halton_normals(2L, 3L, 129L))
  expect_error(choicer:::halton_fill_uniforms(1, 1L, 129L, 0, 0L),
               "between 1 and 128", fixed = TRUE)
  # Generate mode keeps its generator
  gen <- choicer:::.mxl_pred_draws(list(mode = "generate", S = 10L, K_w = 2L,
                                        seed = 4L, scramble = "permuted"),
                                   10^9, "predict()")
  expect_identical(gen$gen_seed, 4L)
  expect_identical(gen$gen_scramble, 1L)
  expect_error(choicer:::.mxl_pred_draws(NULL, 10L, "elasticities()"),
               "elasticities() requires draws_info from a fitted MXL model.",
               fixed = TRUE)
})

test_that("draw codes and draw metadata are checked", {
  d <- mxc_data(180, 5L)
  none <- array(0, dim = c(d$K_w, 0L, 0L))
  msg <- function(v) {
    sprintf("gen_scramble must be 0, 1 or 2 when gen_seed >= 0; got %d.", v)
  }
  for (v in c(-1L, 3L)) {
    expect_error(choicer:::mxl_predict(d$theta, d$X, d$W, d$alt_idx, d$M, none,
                                        d$rc_dist, gen_seed = 0L,
                                        gen_scramble = v, gen_S = d$S),
                 msg(v), fixed = TRUE)
  }
  K_w <- d$K_w
  expect_error(
    mxl_blp_contraction(rep(0, d$J), rep(1 / d$J, d$J), d$X, d$W, d$theta[1:2],
                        d$theta[2L + seq_len(K_w)],
                        d$theta[2L + K_w + seq_len(K_w * (K_w + 1L) / 2L)],
                        d$alt_idx, d$M, d$weights, none, d$rc_dist,
                        gen_seed = 0L, gen_scramble = 3L, gen_S = d$S),
    msg(3L), fixed = TRUE)
  # The estimation kernels have no store-mode points: there 2 would be the
  # identity points through HaltonGen's own inverse normal CDF
  fit <- mxc_fit()$fit
  e <- fit$data
  expect_error(
    choicer:::mxl_loglik_gradient_parallel(
      fit$coefficients, e$X, e$W, e$alt_idx, e$choice_idx, e$M, e$weights,
      array(0, dim = c(2L, 0L, 0L)), fit$rc_dist, fit$rc_correlation,
      fit$rc_mean, fit$use_asc, fit$include_outside_option, gen_seed = 0L,
      gen_scramble = 2L, gen_S = 20L),
    "gen_scramble must be 0 or 1 when gen_seed >= 0; got 2.", fixed = TRUE)
  # A fit's draw metadata, as get_halton_normals() checked it
  info <- list(mode = "store", S = 10L, K_w = 2L)
  for (bad in list(list(S = NA), list(S = 2.5), list(S = 0L))) {
    expect_error(choicer:::.mxl_pred_draws(utils::modifyList(info, bad), 5L,
                                           "predict()"),
                 "`S` must be a single positive whole number.", fixed = TRUE)
  }
  expect_error(choicer:::.mxl_pred_draws(utils::modifyList(info, list(K_w = NULL)),
                                         5L, "predict()"),
               "`K_w` must be a single positive whole number.", fixed = TRUE)
})

# --- Draws kept across blp()'s iterations ----------------------------------------

# mxl_blp_contraction_cached() on `d`, the draws formed on the fly by `gen`,
# kept within `budget` bytes, and mxl_blp_contraction(), which keeps nothing
mxc_blp_kept <- function(d, gen, budget, max_iter = 1000L) {
  K_w <- d$K_w
  target <- as.numeric(mxc_call("shares", d))
  args <- list(rep(0, d$J), target, d$X, d$W, d$theta[1:2],
               d$theta[2L + seq_len(K_w)],
               d$theta[2L + K_w + seq_len(K_w * (K_w + 1L) / 2L)],
               d$alt_idx, d$M, d$weights, array(0, dim = c(K_w, 0L, 0L)),
               d$rc_dist, rc_correlation = TRUE, rc_mean = TRUE,
               include_outside_option = d$ioo, tol = 1e-8, max_iter = max_iter)
  list(formed = do.call(mxl_blp_contraction, c(args, gen)),
       kept = do.call(choicer:::mxl_blp_contraction_cached,
                      c(args, gen, cache_bytes = budget)))
}

test_that("draws kept across blp()'s iterations give the same deltas", {
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  gens <- list(points = list(gen_seed = 0L, gen_scramble = 2L),
               permuted = list(gen_seed = 9L, gen_scramble = 1L),
               plain = list(gen_seed = 9L, gen_scramble = 0L))
  # Every evaluation of the shares after the first reads all N situations'
  # draws from the cache
  read_all <- function(r, d) {
    expect_true(r$kept)
    expect_gt(r$passes, 2)
    expect_identical(r$reads, d$N * (r$passes - 1))
  }
  for (ioo in c(FALSE, TRUE)) {
    for (w_type in c("row", "alt")) {
      d <- mxc_data(160 + ioo, 23L, ioo = ioo, w_type = w_type)
      for (g in names(gens)) {
        r <- mxc_blp_kept(d, c(gens[[g]], gen_S = d$S), Inf)
        lab <- sprintf("%s, ioo %d, W %s", g, ioo, w_type)
        read_all(r$kept, d)
        expect_identical(r$kept$delta, r$formed, label = lab)
      }
    }
  }
  # One draw, and one coefficient
  for (d in list(mxc_data(162, 9L, S = 1L), mxc_data(163, 9L, K_w = 1L))) {
    r <- mxc_blp_kept(d, list(gen_seed = 0L, gen_scramble = 2L, gen_S = d$S), Inf)
    read_all(r$kept, d)
    expect_identical(r$kept$delta, r$formed)
  }
  # max_iter = 0 leaves no later evaluation, so nothing is kept
  d <- mxc_data(164, 9L)
  expect_output(
    r <- mxc_blp_kept(d, list(gen_seed = 0L, gen_scramble = 2L, gen_S = d$S),
                      Inf, max_iter = 0L),
    "Maximum iterations reached without convergence", fixed = TRUE)
  expect_false(r$kept$kept)
  expect_identical(c(r$kept$passes, r$kept$reads), c(1, 0))
  expect_identical(r$kept$delta, r$formed)
  # At two threads, the sums over situations to rounding
  set_num_threads(2L)
  d <- mxc_data(165, 23L)
  r <- mxc_blp_kept(d, list(gen_seed = 0L, gen_scramble = 2L, gen_S = d$S), Inf)
  read_all(r$kept, d)
  expect_equal(r$kept$delta, r$formed, tolerance = 1e-12)
})

test_that("blp() keeps the draws within keep_draws_bytes, and only on the fly", {
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  d <- mxc_data(170, 23L)
  gen <- list(gen_seed = 0L, gen_scramble = 2L, gen_S = d$S)
  bytes <- 8 * d$K_w * d$S * d$N
  for (b in list(list(bytes, TRUE), list(bytes - 1, FALSE), list(0, FALSE),
                 list(NaN, FALSE), list(-Inf, FALSE))) {
    r <- mxc_blp_kept(d, gen, b[[1]])
    lab <- format(b[[1]])
    expect_identical(r$kept$kept, b[[2]], label = lab)
    expect_identical(r$kept$reads,
                     if (b[[2]]) d$N * (r$kept$passes - 1) else 0, label = lab)
    expect_identical(r$kept$delta, r$formed)
  }
  # A cube is read at every pass, not copied
  K_w <- d$K_w
  target <- as.numeric(mxc_call("shares", d))
  r <- choicer:::mxl_blp_contraction_cached(
    rep(0, d$J), target, d$X, d$W, d$theta[1:2], d$theta[2L + seq_len(K_w)],
    d$theta[2L + K_w + seq_len(K_w * (K_w + 1L) / 2L)], d$alt_idx, d$M,
    d$weights, d$eta, d$rc_dist, TRUE, TRUE, d$ioo, 1e-8, 1000L, -1L, 1L, 0L,
    Inf)
  expect_false(r$kept)
  expect_identical(r$reads, 0)
  expect_identical(r$delta, mxc_call("blp", d))

  # blp() passes keep_draws_bytes, by default 4 GiB; a counterfactual target
  # makes it iterate
  real <- choicer:::mxl_blp_contraction_cached
  seen <- NULL
  spy <- function(...) {
    r <- real(...)
    seen <<- rbind(seen, data.frame(budget = list(...)$cache_bytes,
                                    kept = r$kept, passes = r$passes,
                                    reads = r$reads))
    r
  }
  for (draws in c("store", "generate")) {
    fx <- mxc_fit(draws = draws)
    fit <- fx$fit
    N <- length(fit$data$M)
    target <- as.numeric(predict(fit, newdata = fx$cf, type = "shares"))
    seen <- NULL
    with_mocked_bindings(
      {
        delta <- blp(fit, target_shares = target)
        expect_identical(blp(fit, target_shares = target, keep_draws_bytes = 0),
                         delta)
      },
      mxl_blp_contraction_cached = spy)
    expect_identical(seen$budget, c(2^32, 0), label = draws)
    expect_identical(seen$kept, c(TRUE, FALSE), label = draws)
    expect_identical(seen$passes[1], seen$passes[2])
    expect_gt(seen$passes[1], 2)
    expect_identical(seen$reads, c(N * (seen$passes[1] - 1), 0), label = draws)
  }
  # An integer64 budget is read as its value, not its bits
  skip_if_not_installed("bit64")
  seen <- NULL
  with_mocked_bindings(
    blp(fit, target_shares = target, keep_draws_bytes = bit64::as.integer64(2^33)),
    mxl_blp_contraction_cached = spy)
  expect_identical(seen$budget, 2^33)
  expect_true(seen$kept)
  for (bad in list(-1, NA_real_, c(1, 2), "1", NULL)) {
    expect_error(blp(fit, target_shares = target, keep_draws_bytes = bad),
                 "`keep_draws_bytes` must be a single non-negative number of bytes",
                 fixed = TRUE)
  }
})

# --- Large store-mode fits -----------------------------------------------------

# A fitted object without its wall-clock fields
mxc_strip <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (!is.null(names(x))) {
      x <- x[!names(x) %in% c("elapsed_time", "elapsed", "time_elapsed")]
    }
    x[] <- lapply(x, mxc_strip)
  }
  x
}

# Long-format data: N choice situations of J alternatives, T per decision maker
mxc_long <- function(seed, N, J = 4L, T = 1L) {
  set.seed(seed)
  dt <- data.table::data.table(id = rep(seq_len(N), each = J),
                               alt = rep(seq_len(J), N),
                               person = rep(seq_len(N / T), each = J * T))
  dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
  dt[, choice := 0L]
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  dt
}

test_that("a store-mode fit above the cube budget warns, then fits as before", {
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)  # fits at one thread are bit-reproducible
  dt <- mxc_long(150, 60L)
  args <- list(data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
               covariate_cols = "x1", random_var_cols = c("w1", "w2"), S = 20L,
               control = list(maxeval = 20L))
  fit <- function(...) {
    a <- args
    v <- list(...)
    a[names(v)] <- v
    suppressMessages(do.call(run_mxlogit, a))
  }
  expect_no_warning(plain <- fit(), message = "Halton draws")
  local_mocked_bindings(.mxl_cube_budget = function() 1000)
  # 8 x 2 x 20 x 60 = 19,200 bytes
  expect_warning(
    warned <- fit(),
    paste0("draws = \"store\" will hold 18.75 KiB of Halton draws (60 choice ",
           "situations x 20 draws x 2 random coefficients) while the model is ",
           "fitted, and vcov(type = ), wesml_vcov() and conditional_tastes() ",
           "rebuild them whole. draws = \"generate\" forms the draws on the fly ",
           "and stores none (with scramble = \"none\", the same Halton points)."),
    fixed = TRUE)
  expect_identical(mxc_strip(warned), mxc_strip(plain))
  # A panel fit's draws are per decision maker
  expect_warning(
    fit(data = mxc_long(151, 60L, T = 3L), person_col = "person"),
    "6.25 KiB of Halton draws (20 decision makers x 20 draws",
    fixed = TRUE)
  # An invalid S is still reported by get_halton_normals()
  expect_error(fit(S = NA), "`S` must be a single positive whole number.",
               fixed = TRUE)
  expect_error(fit(S = 2.5), "`S` must be a single positive whole number.",
               fixed = TRUE)
  # Generate mode and the advanced workflow build no cube of their own
  expect_no_warning(fit(draws = "generate", seed = 1L), message = "Halton draws")
  d <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("w1", "w2"))
  expect_no_warning(
    suppressMessages(run_mxlogit(input_data = d,
                                 eta_draws = get_halton_normals(20L, d$N, 2L),
                                 control = list(maxeval = 20L))),
    message = "Halton draws")
})

test_that("the store-mode warning starts past the budget and leaves the stop", {
  expect_warning(
    choicer:::.warn_store_cube(100, 2e6, 3L, panel = TRUE),
    "will hold 4.47 GiB of Halton draws (2,000,000 decision makers x 100 draws x 3",
    fixed = TRUE)
  expect_silent(choicer:::.warn_store_cube(100, 1e6, 1L, panel = FALSE))
  # Exactly 1 GiB is within the budget; one more situation is not
  expect_silent(choicer:::.warn_store_cube(1, 2^27, 1L, panel = FALSE))
  expect_warning(choicer:::.warn_store_cube(1, 2^27 + 1, 1L, panel = FALSE),
                 "(134,217,729 choice situations x 1 draw x 1 random coefficient)",
                 fixed = TRUE)
  # The most points get_halton_normals() builds warn; past them it stops
  expect_warning(choicer:::.warn_store_cube(1, 2^31 - 1, 1L, panel = FALSE),
                 "16 GiB of Halton draws", fixed = TRUE)
  expect_silent(choicer:::.warn_store_cube(100, 3e7, 2L, panel = FALSE))
})
