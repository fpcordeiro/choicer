# ============================================================================
# Phase 3 batteries (scale audit: the MXL prediction family on a situation
# layout, chunked store-mode draws), sourced by refactor_identity_check.R
# after p2_batteries.R (whose p2_val / p2_cond / p2_leaf / p2_strip helpers
# they use). Every call exists with the same R signature on main and on the
# Phase 3 builds; the chunk arguments are passed only where the kernel has
# them, so main computes the unchunked result under the same key.
#   * pred_x    : the six kernels on shapes the mxl battery lacks: designs of
#                 1-5 rows (square K_x, the 4/5 boundary), single-row
#                 situations, sets straddling 16, duplicate codes, codes
#                 starting above 1 without ASCs, alternative-level W (and
#                 J == rows), K_w in {0, 1, 3}, S in {0, 1, 3, 4, 7},
#                 m = K_w = S <= 4, odd K_w * S, non-finite utilities, and a
#                 > 10^6-row design; store and generate mode.
#   * pred_panel: person_col fits (store and generate) and every MXL
#                 post-estimation method on them, plus newdata round trips.
#   * pred_dbl  : double-typed alt_idx; equality with the integer call.
#   * err_p3    : every validation path of the six kernels (condition class,
#                 message, call, warnings).
#   * pred_chunk: store-mode draws in chunks of 1, 2, 3, 7, N - 1, N and
#                 N + 1 situations (main: the same calls unchunked).
# Generate mode with a short alternative-level W is in err_p3: on main it
# stops on the primary thread in the W * mu gather (an Armadillo bounds error),
# before the parallel region (p3_abort.R confirms it in subprocesses).
# ============================================================================

P3_K <- c("pred", "logsum", "shares", "elas_fixed", "elas_random",
          "dr_fixed", "dr_random", "blp")

p3_has_chunks <- function() {
  "draw_block" %in% names(formals(choicer:::mxl_predict))
}

# Builds whose prediction kernels form store-mode draws on the fly
# (gen_scramble = 2): the pred_chunk cases run that way there, against main's
# cube under the same keys.
p3_has_store_points <- function() {
  exists(".MXL_STORE_POINTS", envir = asNamespace("choicer"), inherits = FALSE)
}

# The store-mode draw source of a chunked call: n points of the K_w-dimensional
# Halton sequence from index `start`, mapped to standard normals (what
# get_halton_normals() reads a block of units at a time).
p3_block_fun <- function(K_w) {
  force(K_w)
  function(start, n) randtoolbox::halton(n = n, dim = K_w, normal = TRUE,
                                         start = start)
}

# One design. Returns everything the kernels take; `eta` is the full store
# cube (or the generate-mode placeholder).
p3_case <- function(seed, M, J, K_x, K_w, S, codes = NULL, dup = FALSE,
                    w_type = "row", ioo = FALSE, use_asc = TRUE,
                    rc_corr = FALSE, rc_mean = TRUE, rc_dist = NULL,
                    gen = FALSE, scramble = 1L, X_fun = NULL) {
  set.seed(seed)
  N <- length(M)
  alt_list <- lapply(seq_len(N), function(i) {
    pool <- if (is.null(codes)) seq_len(J) else codes
    sort(pool[sample.int(length(pool), M[i], replace = dup)])
  })
  alt_idx <- as.integer(unlist(alt_list))
  n <- sum(M)
  X <- matrix(rnorm(n * K_x), n, K_x)
  if (!is.null(X_fun)) X <- X_fun(X, M)
  W <- if (w_type == "row") {
    matrix(rnorm(n * K_w), n, K_w)
  } else {
    matrix(rnorm(J * K_w), J, K_w)
  }
  rc_dist <- if (is.null(rc_dist)) rep(0L, K_w) else as.integer(rc_dist)
  L_size <- if (rc_corr) K_w * (K_w + 1) / 2 else K_w
  n_par <- K_x + (if (rc_mean) K_w else 0L) + L_size +
    n_delta_free(J, use_asc, ioo)
  theta <- rnorm(n_par, sd = 0.3)
  weights <- runif(N, 0.5, 2)
  eta <- if (gen) {
    array(0, dim = c(K_w, 0L, 0L))
  } else if (K_w == 0L || S == 0L) {
    array(0, dim = c(K_w, S, N))
  } else {
    get_halton_normals(S, N, K_w)
  }
  list(theta = theta, X = X, W = W, alt_idx = alt_idx, M = as.integer(M),
       choice_idx = rep(1L, N), weights = weights, eta = eta,
       rc_dist = rc_dist, rc_corr = rc_corr, rc_mean = rc_mean,
       use_asc = use_asc, ioo = ioo, J = J, K_x = K_x, K_w = K_w, S = S,
       L_size = L_size, gen_seed = if (gen) 11L else -1L,
       gen_scramble = if (gen) as.integer(scramble) else 1L,
       gen_S = if (gen) as.integer(S) else 0L)
}

# Chunk arguments (and the placeholder cube) when the build has them and the
# case asks for chunks.
p3_draws <- function(f) {
  if (!is.null(f$chunk) && f$gen_seed < 0L && p3_has_store_points()) {
    return(list(eta = array(0, dim = c(f$K_w, 0L, 0L)), extra = list(),
                gen = list(gen_seed = 0L, gen_scramble = 2L,
                           gen_S = as.integer(f$S))))
  }
  if (is.null(f$chunk) || !p3_has_chunks()) {
    return(list(eta = f$eta, extra = list()))
  }
  list(eta = array(0, dim = c(f$K_w, f$S, 0L)),
       extra = list(draw_block = p3_block_fun(f$K_w), chunk_size = f$chunk))
}

p3_call <- function(k, f) {
  d <- p3_draws(f)
  head <- list(f$theta, f$X, f$W, f$alt_idx, f$M)
  flags <- list(rc_correlation = f$rc_corr, rc_mean = f$rc_mean,
                use_asc = f$use_asc, include_outside_option = f$ioo,
                gen_seed = f$gen_seed, gen_scramble = f$gen_scramble,
                gen_S = f$gen_S)
  if (!is.null(d$gen)) flags[names(d$gen)] <- d$gen
  var_r <- f$var_r %||% 1L
  switch(k,
    pred = do.call(mxl_predict,
                   c(head, list(eta_draws = d$eta, rc_dist = f$rc_dist), flags,
                     d$extra)),
    logsum = do.call(mxl_logsum,
                     c(head, list(eta_draws = d$eta, rc_dist = f$rc_dist), flags,
                       d$extra)),
    shares = do.call(mxl_predict_shares,
                     c(head, list(weights = f$weights, eta_draws = d$eta,
                                  rc_dist = f$rc_dist), flags, d$extra)),
    elas_fixed = do.call(mxl_elasticities_parallel,
                         c(list(f$theta, f$X, f$W, f$alt_idx, f$choice_idx, f$M,
                                weights = f$weights, eta_draws = d$eta,
                                rc_dist = f$rc_dist,
                                elast_var_idx = f$var_x %||% 1L,
                                is_random_coef = FALSE), flags, d$extra)),
    elas_random = do.call(mxl_elasticities_parallel,
                          c(list(f$theta, f$X, f$W, f$alt_idx, f$choice_idx, f$M,
                                 weights = f$weights, eta_draws = d$eta,
                                 rc_dist = f$rc_dist, elast_var_idx = var_r,
                                 is_random_coef = TRUE), flags, d$extra)),
    dr_fixed = do.call(mxl_diversion_ratios_parallel,
                       c(head, list(weights = f$weights, eta_draws = d$eta,
                                    rc_dist = f$rc_dist,
                                    elast_var_idx = f$var_x %||% 1L,
                                    is_random_coef = FALSE), flags, d$extra)),
    dr_random = do.call(mxl_diversion_ratios_parallel,
                        c(head, list(weights = f$weights, eta_draws = d$eta,
                                     rc_dist = f$rc_dist, elast_var_idx = var_r,
                                     is_random_coef = TRUE), flags, d$extra)),
    # the solver tier's leaf (blp_batteries.R), its arguments formed inside
    # the recorder, which records their errors and warnings as the call's
    blp = golden_blp("mxl_blp_contraction", p3_blp_args(f, d, flags),
                     ref = if (is.null(f$target) && is.null(f$delta) &&
                               is.null(f$max_iter)) p3_blp_ref(f) else FALSE,
                     expect = if (is.null(f$max_iter)) "any" else "capped"))
}

# The BLP call of a case: BLP's alternatives are 1..max(alt_idx), from zeros,
# with the case's own shares over them as the target (or the case's own
# delta and target).
p3_blp_args <- function(f, d, flags) {
  K_x <- ncol(f$X); K_w <- ncol(f$W)
  th <- f$theta
  beta <- f$beta %||% th[seq_len(K_x)]
  mu <- f$mu %||% (if (f$rc_mean) th[K_x + seq_len(K_w)] else rep(0, K_w))
  L_params <- f$L_params %||%
    th[K_x + (if (f$rc_mean) K_w else 0L) + seq_len(f$L_size)]
  J_in <- if (length(f$alt_idx)) max(as.integer(f$alt_idx)) else f$J
  target <- f$target %||%
    as.numeric(p3_call("shares", f))[seq_len(J_in + f$ioo)]
  c(list(delta = f$delta %||% rep(0, J_in), target_shares = target,
         X = f$X, W = f$W, beta = beta, mu = mu,
         L_params = L_params, alt_idx = f$alt_idx, M = f$M,
         weights = f$weights, eta_draws = d$eta,
         rc_dist = f$rc_dist, rc_correlation = f$rc_corr,
         rc_mean = f$rc_mean, include_outside_option = f$ioo,
         tol = 1e-10, max_iter = f$max_iter %||% 200L,
         gen_seed = flags$gen_seed, gen_scramble = flags$gen_scramble,
         gen_S = flags$gen_S), d$extra)
}

# The exact reference of a case's BLP call on its own shares: theta's
# constants (zeros without them), where theta covers BLP's alternatives;
# otherwise NULL (a Newton reference).
p3_blp_ref <- function(f) {
  tryCatch({
    J_in <- max(as.integer(f$alt_idx))
    n_head <- ncol(f$X) + (if (f$rc_mean) ncol(f$W) else 0L) + f$L_size
    n_delta <- length(f$theta) - n_head
    if (f$use_asc && n_delta != J_in - !f$ioo) NULL else
      gb_delta_star(f$theta, n_head, J_in, f$use_asc, f$ioo)
  }, warning = function(w) NULL, error = function(e) NULL)
}

# Each kernel's value, warnings and error, plus anything it printed (the BLP
# contraction prints a warning when it stops at max_iter).
p3_outputs <- function(f, kernels = P3_K) {
  res <- list()
  for (k in kernels) {
    if (k == "blp") {
      res$blp <- p3_call(k, f)
      next
    }
    printed <- character()
    res[[k]] <- p2_val({
      printed <- utils::capture.output(v <- p3_call(k, f))
      v
    })
    if (length(printed)) res[[paste0(k, "_stdout")]] <- printed
  }
  res
}

# ---------------------------------------------------------------------------
# pred_x
# ---------------------------------------------------------------------------

run_pred_x <- function() {
  out <- list()
  add <- function(key, f, kernels = P3_K) out[[key]] <<- p3_outputs(f, kernels)

  # Designs of 1-5 rows: X * beta takes Armadillo's dgemv('T') path at one
  # row, its tiny-square emulation when rows == K_x <= 4, and dgemv('N')
  # otherwise; W * mu likewise with K_w.
  parts <- list(n1 = list(1L), n2 = list(2L, c(1L, 1L)), n3 = list(3L, c(1L, 2L)),
                n4 = list(4L, c(2L, 2L)), n5 = list(5L, c(2L, 3L)))
  seed <- 9100L
  for (nn in names(parts)) {
    n <- as.integer(substring(nn, 2))
    J <- if (n <= 3L) 3L else 5L
    for (ip in seq_along(parts[[nn]])) {
      M <- parts[[nn]][[ip]]
      for (K_x in unique(c(n, 2L))) {
        for (K_w in unique(c(1L, 2L, min(n, 4L)))) {
          for (cfg in c("asc", "asc_oo", "noasc")) {
            seed <- seed + 1L
            f <- p3_case(seed, M, J, K_x, K_w, S = 5L,
                         ioo = cfg == "asc_oo", use_asc = cfg != "noasc",
                         rc_corr = K_w > 1L, rc_mean = TRUE)
            add(sprintf("tiny|%s|M=%s|Kx=%d|Kw=%d|%s", nn,
                        paste(M, collapse = "+"), K_x, K_w, cfg), f)
          }
        }
      }
    }
  }

  # Single-row situations inside a larger design, odd K_w * S (21 values per
  # draw block); store and generate mode.
  set.seed(9200)
  M1 <- sample(c(1L, 1L, 2L, 3L, 6L), 25, replace = TRUE)
  for (ioo in c(FALSE, TRUE)) {
    for (mode in c("store", "gen0", "gen1")) {
      f <- p3_case(9201L + ioo, M1, J = 6L, K_x = 2L, K_w = 3L, S = 7L,
                   ioo = ioo, rc_corr = TRUE, rc_dist = c(0L, 1L, 0L),
                   gen = mode != "store", scramble = if (mode == "gen0") 0L else 1L)
      f$var_r <- 2L
      add(sprintf("m1|oo=%d|%s", ioo, mode), f)
    }
  }

  # Sets straddling 16 rows (14-20 of 30 alternatives).
  set.seed(9300)
  M16 <- sample(14:20, 12, replace = TRUE)
  for (ioo in c(FALSE, TRUE)) {
    for (w_type in c("row", "alt")) {
      f <- p3_case(9301L + ioo, M16, J = 30L, K_x = 3L, K_w = 2L, S = 4L,
                   w_type = w_type, ioo = ioo, rc_corr = TRUE)
      add(sprintf("wide|oo=%d|W=%s", ioo, w_type), f)
    }
  }

  # Duplicate codes within a situation.
  for (ioo in c(FALSE, TRUE)) {
    f <- p3_case(9401L + ioo, c(4L, 5L, 6L, 4L, 5L, 6L, 4L, 5L), J = 4L,
                 K_x = 2L, K_w = 2L, S = 6L, dup = TRUE, ioo = ioo)
    add(sprintf("dup|oo=%d", ioo), f)
  }

  # Codes starting above 1 without ASCs (alternatives 1-2 never appear).
  for (ioo in c(FALSE, TRUE)) {
    f <- p3_case(9501L + ioo, rep(3L, 15), J = 7L, K_x = 2L, K_w = 2L, S = 5L,
                 codes = 3:7, ioo = ioo, use_asc = FALSE)
    add(sprintf("codes_above1|oo=%d", ioo), f)
  }

  # Alternative-level W, and J equal to the number of rows (which the kernels
  # read as row-aligned).
  for (ioo in c(FALSE, TRUE)) {
    f <- p3_case(9601L + ioo, rep(3L, 10), J = 5L, K_x = 2L, K_w = 2L, S = 5L,
                 w_type = "alt", ioo = ioo, rc_dist = c(1L, 0L))
    add(sprintf("altW|oo=%d", ioo), f)
    g <- p3_case(9611L + ioo, c(2L, 3L), J = 5L, K_x = 2L, K_w = 2L, S = 5L,
                 w_type = "alt", ioo = ioo)
    add(sprintf("altW_J_eq_rows|oo=%d", ioo), g)
  }

  # K_w in {0, 1, 3}.
  for (K_w in c(0L, 1L, 3L)) {
    for (ioo in c(FALSE, TRUE)) {
      set.seed(9700L + 10L * K_w + ioo)
      f <- p3_case(9701L + 10L * K_w + ioo, sample(2:4, 12, replace = TRUE),
                   J = 4L, K_x = 2L, K_w = K_w, S = 5L, ioo = ioo,
                   rc_corr = K_w > 1L, rc_mean = K_w > 0L)
      if (K_w == 0L) f$mu <- numeric(0)
      add(sprintf("Kw=%d|oo=%d", K_w, ioo), f)
    }
  }

  # S in {0, 1, 3, 4, 7} (S = 0: store mode only).
  for (S in c(0L, 1L, 3L, 4L, 7L)) {
    for (mode in c("store", "gen")) {
      if (S == 0L && mode == "gen") next
      set.seed(9800L + S)
      f <- p3_case(9801L + S, sample(2:4, 12, replace = TRUE), J = 4L,
                   K_x = 2L, K_w = 2L, S = S, ioo = TRUE, rc_corr = TRUE,
                   gen = mode == "gen")
      add(sprintf("S=%d|%s", S, mode), f)
    }
  }

  # m = K_w = S <= 4: every product in a situation is Armadillo's tiny-square
  # emulation.
  for (m in 2:4) {
    for (ioo in c(FALSE, TRUE)) {
      f <- p3_case(9901L + 10L * m + ioo, rep(m, 9), J = m, K_x = 2L, K_w = m,
                   S = m, ioo = ioo, rc_corr = TRUE)
      add(sprintf("mKS=%d|oo=%d", m, ioo), f)
    }
  }

  # Non-finite and extreme utilities: a utility of ~3000 in situation 1, +Inf
  # in situation 2, NaN in situation 3, then finite situations (a thread's
  # buffers must not carry them over).
  nonfinite <- function(X, M) {
    r0 <- c(0L, cumsum(M))
    X[r0[1] + 1L, 1L] <- 1e4
    X[r0[2] + 2L, 1L] <- Inf
    X[r0[3] + 1L, 2L] <- NaN
    X
  }
  for (ioo in c(FALSE, TRUE)) {
    f <- p3_case(9951L + ioo, rep(3L, 8), J = 4L, K_x = 2L, K_w = 2L, S = 5L,
                 ioo = ioo, X_fun = nonfinite)
    f$theta[1:2] <- c(0.3, 0.2)
    add(sprintf("nonfinite|oo=%d", ioo), f,
        kernels = setdiff(P3_K, "blp"))
  }

  # More than 10^6 rows: the layout's validation pass runs in parallel there.
  f <- p3_case(9990L, rep(5L, 210000L), J = 5L, K_x = 2L, K_w = 1L, S = 2L,
               rc_mean = TRUE)
  add("big|rows=1.05e6", f, kernels = c("pred", "logsum", "shares"))
  out
}

# ---------------------------------------------------------------------------
# pred_panel
# ---------------------------------------------------------------------------

run_pred_panel <- function() {
  out <- list()
  sim <- simulate_mxl_data(N = 60, J = 4, T = 3L, seed = 11)
  dt <- sim$data
  for (mode in c("store", "generate")) {
    fit <- suppressMessages(run_mxlogit(
      data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
      outside_opt_label = 0L, include_outside_option = TRUE,
      person_col = "pid", S = 20L, draws = mode,
      seed = if (mode == "generate") 5L else NULL))
    res <- list(coef = p2_leaf(coef(fit)))
    res$pred <- p2_val(predict(fit))
    res$pred_new <- p2_val(predict(fit, newdata = dt))
    res$shares <- p2_val(predict(fit, type = "shares"))
    res$logsum <- p2_val(logsum(fit))
    res$logsum_new <- p2_val(logsum(fit, newdata = dt))
    res$cs <- p2_val(consumer_surplus(fit, price_var = "x1"))
    res$elas_fixed <- p2_val(elasticities(fit, "x1"))
    res$elas_random <- p2_val(elasticities(fit, "w1", is_random_coef = TRUE))
    res$dr_fixed <- p2_val(diversion_ratios(fit, "x1"))
    res$dr_random <- p2_val(diversion_ratios(fit, "w2", is_random_coef = TRUE))
    shares <- predict(fit, type = "shares")
    res$blp <- golden_blp("blp", list(fit, target_shares = as.numeric(shares)),
                          ref = gb_fit_delta(fit))
    res$gof <- p2_val(gof(fit))
    out[[sprintf("panel|%s", mode)]] <- res
  }
  out
}

# ---------------------------------------------------------------------------
# pred_dbl
# ---------------------------------------------------------------------------

run_pred_dbl <- function() {
  out <- list()
  for (ioo in c(FALSE, TRUE)) {
    for (mode in c("store", "gen")) {
      set.seed(9650L + ioo)
      f <- p3_case(9601L + ioo, sample(2:5, 20, replace = TRUE), J = 5L,
                   K_x = 2L, K_w = 2L, S = 6L, ioo = ioo, rc_corr = TRUE,
                   gen = mode == "gen")
      a <- p3_outputs(f)
      fd <- f
      fd$alt_idx <- as.double(f$alt_idx)
      b <- p3_outputs(fd)
      out[[sprintf("dbl|oo=%d|%s", ioo, mode)]] <-
        list(dbl = b, same_as_int = identical(a, b))
    }
  }
  out
}

# ---------------------------------------------------------------------------
# err_p3
# ---------------------------------------------------------------------------

run_err_p3 <- function() {
  out <- list()
  J <- 4L; K_x <- 2L; K_w <- 2L; S <- 4L
  base <- function(ioo, use_asc = TRUE, gen = FALSE) {
    f <- p3_case(9701L + ioo + 2L * use_asc, c(2L, 3L, 4L, 2L, 3L, 4L, 4L, 3L,
                                               2L, 4L),
                 J = J, K_x = K_x, K_w = K_w, S = S, ioo = ioo,
                 use_asc = use_asc, rc_corr = TRUE, gen = gen)
    f$max_iter <- 5L
    f$target <- rep(1 / (J + ioo), J + ioo)
    f
  }
  f0 <- base(FALSE); fo <- base(TRUE); fn <- base(FALSE, FALSE)
  fg <- base(FALSE, gen = TRUE)
  n <- sum(f0$M)
  w <- function(f, ...) p2_with(f, ...)
  alt_rep <- function(f, i, v) { a <- f$alt_idx; a[i] <- v; a }
  empty <- function(f) {
    w(f, X = f$X[0, , drop = FALSE],
      W = if (nrow(f$W) == sum(f$M)) f$W[0, , drop = FALSE] else f$W,
      alt_idx = integer(0), choice_idx = integer(0), M = integer(0),
      weights = numeric(0), eta = array(0, c(K_w, S, 0L)))
  }
  wide_kw <- function() {
    f <- p3_case(9790L, rep(2L, 6), J = 3L, K_x = 2L, K_w = 129L, S = 3L,
                 gen = TRUE)
    f
  }
  x0 <- function(f) {
    w(f, X = f$X[, 0, drop = FALSE], theta = f$theta[-seq_len(K_x)],
      beta = numeric(0))
  }
  cases <- list(
    ok = f0, ok_oo = fo, ok_noasc = fn, ok_gen = fg,
    M_zero = w(f0, M = replace(f0$M, 3L, 0L)),
    M_na = w(f0, M = replace(f0$M, 3L, NA_integer_)),
    M_neg = w(f0, M = replace(f0$M, 3L, -2L)),
    X_short = w(f0, X = f0$X[-1L, , drop = FALSE]),
    M_long = w(f0, M = c(f0$M, 2L)),
    alt_short = w(f0, alt_idx = f0$alt_idx[-1L]),
    w_short = w(f0, weights = f0$weights[-1L]),
    alt_zero = w(f0, alt_idx = alt_rep(f0, 5L, 0L)),
    alt_na = w(f0, alt_idx = alt_rep(f0, 5L, NA_integer_)),
    alt_neg = w(f0, alt_idx = alt_rep(f0, 5L, -1L)),
    alt_3e9 = w(f0, alt_idx = replace(as.double(f0$alt_idx), 5L, 3e9)),
    alt_3e9_noasc = w(fn, alt_idx = replace(as.double(fn$alt_idx), 5L, 3e9)),
    alt_big = w(f0, alt_idx = alt_rep(f0, 5L, J + 1L)),
    alt_big_noasc = w(fn, alt_idx = alt_rep(fn, 5L, J + 1L)),
    alt_frac = w(f0, alt_idx = as.double(f0$alt_idx) + 0.25),
    alt_lgl = w(f0, alt_idx = rep(TRUE, n)),
    alt_factor = w(f0, alt_idx = factor(f0$alt_idx)),
    alt_null = w(f0, alt_idx = NULL),
    alt_chr = w(f0, alt_idx = as.character(f0$alt_idx)),
    alt_zero_gen = w(fg, alt_idx = alt_rep(fg, 5L, 0L)),
    alt_neg_gen = w(fg, alt_idx = alt_rep(fg, 5L, -1L)),
    alt_big_gen = w(fg, alt_idx = alt_rep(fg, 5L, J + 1L)),
    theta_short = w(f0, theta = f0$theta[1:2]),
    theta_long_noasc = w(fn, theta = c(fn$theta, 0.1)),
    delta_short = w(f0, theta = f0$theta[-length(f0$theta)]),
    rc_dist_short = w(f0, rc_dist = 0L),
    rc_dist_short_gen = w(fg, rc_dist = 0L),
    empty_asc = empty(f0), empty_noasc = empty(fn), empty_oo = empty(fo),
    w_zero = w(f0, weights = rep(0, length(f0$weights))),
    eta_slices = w(f0, eta = get_halton_normals(S, length(f0$M) + 1L, K_w)),
    eta_rows = w(f0, eta = get_halton_normals(S, length(f0$M), K_w + 1L)),
    W_alt_short = w(f0, W = f0$W[seq_len(J - 1L), , drop = FALSE]),
    W_alt_ok = w(f0, W = f0$W[seq_len(J), , drop = FALSE]),
    W_alt_short_gen = w(fg, W = fg$W[seq_len(J - 1L), , drop = FALSE]),
    W_alt_ok_gen = w(fg, W = fg$W[seq_len(J), , drop = FALSE]),
    gen_S_zero = w(fg, gen_S = 0L),
    gen_Kw129 = wide_kw(),
    var_oob = w(f0, var_x = 3L, var_r = 3L),
    var_zero = w(f0, var_x = 0L, var_r = 0L),
    X0 = x0(f0),
    Kw0_random = local({
      g <- p3_case(9795L, rep(3L, 6), J = 3L, K_x = 2L, K_w = 0L, S = 3L,
                   rc_mean = FALSE)
      g$mu <- numeric(0)
      g
    }),
    blp_target_len = w(f0, target = rep(0.25, 3)),
    blp_target_zero = w(f0, target = c(0.5, 0.5, 0, 0)),
    blp_delta_len = w(f0, delta = rep(0, J + 1L)),
    blp_delta_free = w(f0, delta = rep(0, J - 1L)),
    blp_delta_len_oo = w(fo, delta = rep(0, J - 1L)),
    blp_beta_len = w(f0, beta = c(f0$theta[seq_len(K_x)], 0.1)),
    blp_mu_len = w(f0, mu = c(0.1, 0.2, 0.3)),
    blp_L_short = w(f0, L_params = 0.1),
    blp_maxiter = w(f0, max_iter = 1L))
  for (cn in names(cases)) {
    for (k in P3_K) {
      out[[sprintf("err|%s|%s", cn, k)]] <- if (k == "blp") {
        p3_call(k, cases[[cn]])
      } else {
        p2_cond(p3_call(k, cases[[cn]]))
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# pred_chunk
# ---------------------------------------------------------------------------

run_pred_chunk <- function() {
  out <- list()
  N <- 30L
  for (ioo in c(FALSE, TRUE)) {
    for (use_asc in c(TRUE, FALSE)) {
      for (w_type in c("row", "alt")) {
        set.seed(9850L + ioo + 2L * use_asc)
        f <- p3_case(9851L + ioo + 2L * use_asc, sample(2:5, N, replace = TRUE),
                     J = 5L, K_x = 2L, K_w = 2L, S = 6L, w_type = w_type,
                     ioo = ioo, use_asc = use_asc, rc_corr = TRUE,
                     rc_dist = c(0L, 1L))
        f$var_r <- 2L
        kernels <- if (use_asc) P3_K else setdiff(P3_K, "blp")
        for (cs in c(1L, 2L, 3L, 7L, N - 1L, N, N + 1L)) {
          f$chunk <- cs
          out[[sprintf("chunk|oo=%d|asc=%d|W=%s|size=%d", ioo, use_asc,
                       w_type, cs)]] <- p3_outputs(f, kernels)
        }
      }
    }
  }
  out
}
