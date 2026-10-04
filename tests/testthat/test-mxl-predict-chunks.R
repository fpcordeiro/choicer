# Store-mode post-estimation regenerates the Halton draws a chunk of choice
# situations at a time when the K_w x S x N cube would exceed its budget
# (.mxl_cube_budget(), 1 GiB). The chunks' draws are the cube's slices, bit for
# bit, so chunked predictions equal the whole-cube ones: per situation at any
# thread count, sums over situations at one thread.

# --- Kernels -------------------------------------------------------------------

# The draw source of a chunked kernel call: points start, ..., start + n - 1 of
# the K_w-dimensional sequence, as get_halton_normals() reads them
mxc_block <- function(K_w) {
  force(K_w)
  function(start, n) randtoolbox::halton(n = n, dim = K_w, normal = TRUE,
                                         start = start)
}

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
       eta = get_halton_normals(S, N, K_w), K_w = K_w, S = S, N = N, J = J,
       ioo = ioo, use_asc = use_asc)
}

# Call kernel `k` on `d` with the whole cube (chunk = NULL) or in chunks
mxc_call <- function(k, d, chunk = NULL, var = 1L, random = FALSE) {
  eta <- d$eta
  extra <- list()
  if (!is.null(chunk)) {
    eta <- array(0, dim = c(d$K_w, d$S, 0L))
    extra <- list(draw_block = mxc_block(d$K_w), chunk_size = chunk)
  }
  flags <- list(rc_correlation = TRUE, rc_mean = TRUE, use_asc = d$use_asc,
                include_outside_option = d$ioo)
  switch(k,
    pred = do.call(choicer:::mxl_predict,
                   c(list(d$theta, d$X, d$W, d$alt_idx, d$M, eta, d$rc_dist),
                     flags, extra)),
    logsum = do.call(choicer:::mxl_logsum,
                     c(list(d$theta, d$X, d$W, d$alt_idx, d$M, eta, d$rc_dist),
                       flags, extra)),
    shares = do.call(choicer:::mxl_predict_shares,
                     c(list(d$theta, d$X, d$W, d$alt_idx, d$M, d$weights, eta,
                            d$rc_dist), flags, extra)),
    elas = do.call(choicer:::mxl_elasticities_parallel,
                   c(list(d$theta, d$X, d$W, d$alt_idx, NULL, d$M, d$weights,
                          eta, d$rc_dist, var, random), flags, extra)),
    dr = do.call(choicer:::mxl_diversion_ratios_parallel,
                 c(list(d$theta, d$X, d$W, d$alt_idx, d$M, d$weights, eta,
                        d$rc_dist, var, random), flags, extra)),
    blp = {
      K_w <- d$K_w
      target <- as.numeric(mxc_call("shares", d))
      args <- list(rep(0, d$J), target, d$X, d$W, d$theta[1:2],
                   d$theta[2L + seq_len(K_w)],
                   d$theta[2L + K_w + seq_len(K_w * (K_w + 1L) / 2L)],
                   d$alt_idx, d$M, d$weights, eta, d$rc_dist,
                   rc_correlation = TRUE, rc_mean = TRUE,
                   include_outside_option = d$ioo)
      if (is.null(chunk)) do.call(mxl_blp_contraction, args) else
        do.call(choicer:::mxl_blp_contraction_chunked, c(args, extra))
    })
}

test_that("chunked store-mode draws give the whole cube's results", {
  on.exit(set_num_threads(2L), add = TRUE)
  N <- 23L
  for (ioo in c(FALSE, TRUE)) {
    for (w_type in c("row", "alt")) {
      d <- mxc_data(100 + ioo, N, ioo = ioo, w_type = w_type)
      set_num_threads(1L)
      whole <- lapply(c(pred = "pred", logsum = "logsum", shares = "shares",
                        elas = "elas", dr = "dr", blp = "blp"),
                      function(k) mxc_call(k, d))
      for (chunk in c(1L, 2L, 3L, 7L, N - 1L, N, N + 1L)) {
        set_num_threads(1L)
        for (k in names(whole)) {
          expect_identical(mxc_call(k, d, chunk), whole[[k]],
                           label = sprintf("%s, chunks of %d", k, chunk))
        }
        # Per situation at two threads too; sums over situations to rounding
        set_num_threads(2L)
        expect_identical(mxc_call("pred", d, chunk), whole$pred)
        expect_identical(mxc_call("logsum", d, chunk), whole$logsum)
        expect_equal(mxc_call("shares", d, chunk), whole$shares,
                     tolerance = 1e-13)
        expect_equal(mxc_call("dr", d, chunk, 2L, TRUE),
                     mxc_call("dr", d, NULL, 2L, TRUE), tolerance = 1e-13)
      }
    }
  }
})

test_that("chunks without ASCs, with one draw, and with one coefficient", {
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  d <- mxc_data(110, 9L, use_asc = FALSE)
  for (k in c("pred", "logsum", "shares", "elas", "dr")) {
    expect_identical(mxc_call(k, d, 4L), mxc_call(k, d))
  }
  d <- mxc_data(111, 9L, S = 1L)
  for (k in c("pred", "logsum", "shares")) {
    expect_identical(mxc_call(k, d, 2L), mxc_call(k, d))
  }
  d <- mxc_data(112, 9L, K_w = 1L)
  for (k in c("pred", "logsum", "shares", "elas")) {
    expect_identical(mxc_call(k, d, 3L), mxc_call(k, d))
  }
})

test_that("chunked kernels check their draw source", {
  d <- mxc_data(120, 6L)
  e <- array(0, dim = c(2L, 6L, 0L))
  call_pred <- function(draw_block, chunk_size, eta = e) {
    choicer:::mxl_predict(d$theta, d$X, d$W, d$alt_idx, d$M, eta, d$rc_dist,
                          TRUE, TRUE, TRUE, FALSE, draw_block = draw_block,
                          chunk_size = chunk_size)
  }
  expect_error(call_pred("halton", 2), "draw_block must be NULL or a function.",
               fixed = TRUE)
  for (cs in list(0, -1, 2.5, NA_real_)) {
    expect_error(call_pred(mxc_block(2L), cs),
                 "chunk_size must be a positive whole number when draw_block is given.",
                 fixed = TRUE)
  }
  # The cube still gives K_w
  expect_error(call_pred(mxc_block(2L), 2, eta = array(0, dim = c(3L, 6L, 0L))),
               "eta_draws 1st dimension (3) does not match K_w (2)", fixed = TRUE)
  # Each block is checked
  expect_error(call_pred(function(start, n) matrix(0, n, 3L), 2),
               "draw_block(start, n) must return n points of 2 coordinates as doubles: 24 values for start = 1 and n = 12; got double of length 36.",
               fixed = TRUE)
  expect_error(call_pred(function(start, n) matrix(0L, n, 2L), 2),
               "got integer of length 24.", fixed = TRUE)
  # An error in draw_block is R's
  expect_error(call_pred(function(start, n) stop("no draws here"), 2),
               "no draws here", fixed = TRUE)
  # Generate mode ignores draw_block
  g <- choicer:::mxl_predict(d$theta, d$X, d$W, d$alt_idx, d$M,
                             array(0, c(2L, 0L, 0L)), d$rc_dist, TRUE, TRUE,
                             TRUE, FALSE, 5L, 1L, 6L)
  expect_identical(
    choicer:::mxl_predict(d$theta, d$X, d$W, d$alt_idx, d$M,
                          array(0, c(2L, 0L, 0L)), d$rc_dist, TRUE, TRUE, TRUE,
                          FALSE, 5L, 1L, 6L, draw_block = "ignored",
                          chunk_size = -1),
    g)
})

test_that("chunks ask draw_block for whole situations in order", {
  d <- mxc_data(130, 10L)
  asked <- NULL
  block <- function(start, n) {
    asked <<- rbind(asked, c(start, n))
    mxc_block(2L)(start, n)
  }
  e <- array(0, dim = c(2L, 6L, 0L))
  choicer:::mxl_logsum(d$theta, d$X, d$W, d$alt_idx, d$M, e, d$rc_dist, TRUE,
                       TRUE, TRUE, FALSE, draw_block = block, chunk_size = 4)
  # Chunks of 4, 4 and 2 situations of 6 draws: points 1-24, 25-48, 49-60
  expect_identical(asked, rbind(c(1, 24), c(25, 24), c(49, 12)))
})

test_that("a chunk spanning several draw_block() calls places each block", {
  skip_on_cran()  # about 0.5 GB at its peak
  # 2^22 + 1 draws per situation: draw_block() is asked for about 2^22 values
  # at a time, so each call returns one situation and a chunk of three takes
  # three calls, as production chunks (2^30 bytes) take about 32; the four
  # situations' draws pass 2^24 values, after which R collects its garbage
  S <- 2^22 + 1
  set.seed(160)
  M <- c(3L, 2L, 3L, 2L)
  alt_idx <- as.integer(unlist(lapply(M, function(m) sort(sample(3L, m)))))
  n <- sum(M)
  X <- matrix(rnorm(n * 2L), n, 2L)
  W <- matrix(rnorm(n), n, 1L)
  theta <- c(rnorm(2L), 0.3, log(0.5), rnorm(2L, sd = 0.3))
  whole <- choicer:::mxl_predict(theta, X, W, alt_idx, M,
                                 get_halton_normals(S, 4L, 1L), 0L, TRUE, TRUE,
                                 TRUE, FALSE)
  asked <- NULL
  block <- function(start, n) {
    asked <<- rbind(asked, c(start, n))
    randtoolbox::halton(n = n, dim = 1L, normal = TRUE, start = start)
  }
  chunked <- choicer:::mxl_predict(theta, X, W, alt_idx, M,
                                   array(0, dim = c(1L, S, 0L)), 0L, TRUE,
                                   TRUE, TRUE, FALSE, draw_block = block,
                                   chunk_size = 3)
  expect_identical(chunked, whole)
  expect_identical(asked, rbind(c(1, S), c(S + 1, S), c(2 * S + 1, S),
                                c(3 * S + 1, S)))
})

test_that("chunk loading keeps every R object it holds protected", {
  skip_on_cran()  # gctorture(): a collection at every R allocation
  # 2^14 draws per situation: each block (256 KB) is large enough that a
  # block released too early is unmapped, not left readable
  S <- 2^14
  d <- mxc_data(170, 4L, S = S)
  # The points of d$eta, precomputed: the source barely allocates, so the
  # collections test the kernel's own objects
  H <- randtoolbox::halton(n = 4 * S, dim = 2L, normal = TRUE)
  src <- function(start, n) H[start:(start + n - 1), , drop = FALSE]
  whole <- mxc_call("pred", d)
  gctorture(TRUE)
  chunked <- tryCatch(
    choicer:::mxl_predict(d$theta, d$X, d$W, d$alt_idx, d$M,
                          array(0, dim = c(2L, S, 0L)), d$rc_dist, TRUE, TRUE,
                          TRUE, FALSE, draw_block = src, chunk_size = 1),
    finally = gctorture(FALSE))
  expect_identical(chunked, whole)
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

test_that("store-mode methods give the same results in chunks", {
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  for (person in c(FALSE, TRUE)) {
    fx <- mxc_fit(person)
    whole <- mxc_methods(fx$fit, fx$cf)
    target <- as.numeric(whole$shares)
    delta <- blp(fx$fit, target_shares = target)
    # 8 * K_w * S = 320 bytes per situation: chunks of 3 situations. No site
    # may build the whole cube.
    with_mocked_bindings(
      {
        expect_identical(mxc_methods(fx$fit, fx$cf), whole)
        expect_message(d2 <- blp(fx$fit, target_shares = target),
                       "regenerates them on one thread, a chunk of situations at a time",
                       fixed = TRUE)
        expect_identical(d2, delta)
      },
      .mxl_cube_budget = function() 1000,
      get_halton_normals = function(...) stop("the whole cube was built"),
      .halton_cube = function(...) stop("the whole cube was built"))
  }
})

test_that("the budget decides between the cube and chunks", {
  info <- list(mode = "store", S = 10L, K_w = 2L, N = 50L)  # 160 B per situation
  local_mocked_bindings(.mxl_cube_budget = function() 8000)
  at <- choicer:::.mxl_pred_draws(info, 50L, "predict()")
  expect_null(at$draw_block)
  expect_identical(at$eta_draws, get_halton_normals(10L, 50L, 2L))
  over <- choicer:::.mxl_pred_draws(info, 51L, "predict()")
  expect_true(is.function(over$draw_block))
  expect_identical(over$chunk_size, 50)
  expect_identical(dim(over$eta_draws), c(2L, 10L, 0L))
  # A situation larger than the budget makes chunks of one
  local_mocked_bindings(.mxl_cube_budget = function() 100)
  expect_identical(choicer:::.mxl_pred_draws(info, 2L, "predict()")$chunk_size, 1)
  # Generate mode never chunks
  gen <- choicer:::.mxl_pred_draws(list(mode = "generate", S = 10L, K_w = 2L,
                                        seed = 4L, scramble = "permuted"),
                                   10^9, "predict()")
  expect_null(gen$draw_block)
  expect_identical(gen$gen_seed, 4L)
  # Chunks reach at most 2^31 - 1 points of the sequence, as the cube does
  expect_error(choicer:::.mxl_pred_draws(info, 3e8, "predict()"),
               "more than 2^31 - 1", fixed = TRUE)
  expect_error(choicer:::.mxl_pred_draws(NULL, 10L, "elasticities()"),
               "elasticities() requires draws_info from a fitted MXL model.",
               fixed = TRUE)
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
