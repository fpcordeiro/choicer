# One variance path: a fit's variance is the one vcov(fit, type = se_method)
# returns post hoc, computed at the natural estimates in natural units, and
# every fit's information is inverted equilibrated (Jacobi, by powers of two).

# A symmetric positive definite A with a known inverse, and the badly scaled
# D A D with D spanning 18 orders of magnitude: its reciprocal condition
# number is far below the machine epsilon.
vop_bad_pair <- function(seed) {
  set.seed(seed)
  p <- 6L
  Q <- qr.Q(qr(matrix(rnorm(p * p), p)))
  A <- Q %*% diag(seq(1, 10, length.out = p)) %*% t(Q)
  A <- (A + t(A)) / 2
  s <- 10^c(-9, -4, 0, 3, 7, 9)
  list(A = A, s = s, H = A * tcrossprod(s))
}

# max |V - V_true| / sqrt(V_true[i,i] V_true[j,j]): scale-free.
vop_err <- function(V, V_true) {
  max(abs(V - V_true) / sqrt(tcrossprod(diag(V_true))))
}

# The value of `expr` and the messages it emitted.
vop_msgs <- function(expr) {
  msgs <- character()
  value <- withCallingHandlers(expr, message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  list(value = value, messages = msgs)
}

test_that("the Jacobi scale is a power of two that brings the diagonal near one", {
  bp <- vop_bad_pair(1L)
  D <- choicer:::.jacobi_scale(bp$H)
  expect_identical(log2(D), round(log2(D)))
  dd <- diag(bp$H * tcrossprod(D))
  expect_true(all(dd >= 0.5 & dd <= 2))
  # by magnitude: a negative diagonal entry is scaled like a positive one
  expect_identical(choicer:::.jacobi_scale(diag(c(4, 0, -4, NaN, Inf, NA))),
                   c(0.5, 1, 0.5, 1, 1, 1))
  # Clamped to [2^-500, 2^500] (2^532 and 2^-512 before the clamp).
  expect_identical(choicer:::.jacobi_scale(diag(c(1e-320, 1e308))),
                   c(2^500, 2^-500))
})

test_that("the equilibrated inversion recovers a badly scaled inverse solve() refuses", {
  bp <- vop_bad_pair(2L)
  V_true <- solve(bp$A) / tcrossprod(bp$s)
  expect_error(solve(bp$H), "computationally singular")
  eq <- vop_msgs(choicer:::invert_hessian(bp$H))
  expect_length(eq$messages, 0L)
  expect_lt(vop_err(eq$value$vcov, V_true), 1e-12)
  expect_lt(max(abs(eq$value$se / sqrt(diag(V_true)) - 1)), 1e-12)
})

test_that("an equilibrated inversion reports negative variances and singularity", {
  # Away from the optimum: a negative diagonal entry, badly scaled. Scaled by
  # its magnitude it inverts; the negative variance is reported and its
  # standard error set to NA.
  A <- matrix(c(-1, 0.2, 0.1, 0.2, 2, 0.3, 0.1, 0.3, 3), 3)
  s <- c(1e9, 1, 1)
  H <- A * tcrossprod(s)
  r <- vop_msgs(choicer:::invert_hessian(H))
  expect_match(r$messages[1], "not positive definite; 1 variance")
  V_true <- solve(A) / tcrossprod(s)
  expect_lt(max(abs(r$value$vcov - V_true) / sqrt(abs(tcrossprod(diag(V_true))))),
            1e-12)
  expect_true(is.na(r$value$se[1]))
  expect_equal(r$value$se[2:3], sqrt(diag(V_true)[2:3]), tolerance = 1e-12)
  # A singular matrix stays singular.
  S <- matrix(c(1, 2, 2, 4), 2) * tcrossprod(c(1e6, 1))
  r <- vop_msgs(choicer:::invert_hessian(S))
  expect_null(r$value$vcov)
  expect_match(r$messages[1], "likely singular")
  expect_true(all(is.na(r$value$se)))
  r <- vop_msgs(choicer:::.sandwich_combine(S, diag(2)))
  expect_null(r$value$vcov)
  expect_match(r$messages[1], "likely singular")
})

test_that("the equilibrated inverse and sandwich are the two-step formulas", {
  # The scaled matrix factored in one working copy (the sandwich solved
  # twice against the factors, without an explicit inverse), against R's
  # solve() of the scaled matrices and their products.
  bp <- vop_bad_pair(5L)
  set.seed(6)
  B <- crossprod(matrix(rnorm(36), 6)) * tcrossprod(bp$s)
  An <- bp$H                       # and a bread that is not symmetric
  An[2, 5] <- An[2, 5] * (1 + 1e-3)
  for (A in list(bp$H, An)) {
    D <- choicer:::.jacobi_scale(A)
    S <- tcrossprod(D)
    V <- choicer:::equilibrated_inverse(A, D)
    expect_lt(vop_err(V, solve(A * S) * S), 1e-13)
    Ai <- solve(A * S)
    W <- Ai %*% (B * S) %*% Ai
    W <- (W + t(W)) / 2 * S
    V <- choicer:::equilibrated_sandwich(A, B, D)
    expect_identical(V, t(V))
    expect_lt(vop_err(V, W), 1e-13)
  }
  # Unit scales: the plain inverse, and the input is left as it was.
  H0 <- vop_bad_pair(2L)$A + diag(6)
  H0_copy <- H0 + 0
  expect_lt(vop_err(choicer:::equilibrated_inverse(H0, rep(1, 6)), solve(H0)),
            1e-14)
  expect_identical(H0, H0_copy)
})

test_that("the equilibrated helpers hold past one 64 x 64 tile", {
  # The sandwich transposes and symmetrizes by tiles of 64 x 64 (columns of
  # tiles in parallel from order 256): orders at, just past and well past a
  # tile boundary, and past the parallel threshold; a badly scaled bread,
  # symmetric and not, and a positive semi-definite meat of lower rank.
  vop_big <- function(p, seed) {
    set.seed(seed)
    Q <- qr.Q(qr(matrix(rnorm(p * p), p)))
    A <- Q %*% diag(seq(1, 20, length.out = p)) %*% t(Q)
    A <- (A + t(A)) / 2
    s <- 10^runif(p, -4, 4)
    G <- matrix(rnorm(p * (p - 3L)), p)
    list(A = A * tcrossprod(s), B = tcrossprod(G) * tcrossprod(s))
  }
  for (p in c(64L, 65L, 130L, 300L)) {
    bp <- vop_big(p, p)
    An <- bp$A
    An[2, p] <- An[2, p] * (1 + 1e-3)
    # Observed errors are at most about 3e-15 (OpenBLAS, reference LAPACK).
    for (k in 1:2) {
      A <- list(bp$A, An)[[k]]
      what <- sprintf("p = %d, %s bread", p, c("symmetric", "non-symmetric")[k])
      A0 <- A + 0
      B0 <- bp$B + 0
      D <- choicer:::.jacobi_scale(A)
      S <- tcrossprod(D)
      Ai <- solve(A * S)
      W <- Ai %*% (bp$B * S) %*% Ai
      W <- (W + t(W)) / 2 * S
      V <- choicer:::equilibrated_sandwich(A, bp$B, D)
      expect_identical(V, t(V), label = paste("sandwich symmetric,", what))
      expect_lt(vop_err(V, W), 1e-12, label = paste("sandwich error,", what))
      expect_lt(vop_err(choicer:::equilibrated_inverse(A, D), Ai * S), 1e-12,
                label = paste("inverse error,", what))
      expect_identical(A, A0, label = paste("bread unchanged,", what))
      expect_identical(bp$B, B0, label = paste("meat unchanged,", what))
    }
  }
})

test_that("the equilibrated inverse is solve() on the scaled matrix, bit for bit", {
  # The same LAPACK calls (dgetrf, then dgetrs on the identity) on the same
  # values: bit for bit on reference BLAS, OpenBLAS and Accelerate, not
  # promised everywhere.
  skip_on_cran()
  for (seed in 7:9) {
    bp <- vop_bad_pair(seed)
    D <- choicer:::.jacobi_scale(bp$H)
    S <- tcrossprod(D)
    expect_identical(choicer:::equilibrated_inverse(bp$H, D),
                     solve(bp$H * S) * S)
  }
})

# The outcome of an inversion: its value, or its error message with the
# LAPACK routine's name dropped (solve() calls dgesv, the helpers dgetrf).
vop_outcome <- function(expr) {
  tryCatch(list(value = expr), error = function(e) {
    list(error = sub("dgesv|dgetrf", "LU", conditionMessage(e)))
  })
}

test_that("the equilibrated inverse keeps solve()'s checks", {
  # Errors where solve() errors (an exact zero pivot; a reciprocal condition
  # number below the machine epsilon), and the same value where it does not.
  cases <- list(singular = matrix(c(1, 2, 2, 4), 2),
                near_singular = diag(c(1, 1e-17)), one = matrix(4))
  for (nm in names(cases)) {
    M <- cases[[nm]]
    expect_identical(
      vop_outcome(choicer:::equilibrated_inverse(M, rep(1, nrow(M)))),
      vop_outcome(solve(M)), label = nm)
  }
  expect_error(choicer:::equilibrated_inverse(matrix(0, 0, 0), numeric(0)),
               "0-diml")
  expect_error(choicer:::equilibrated_inverse(matrix(1, 2, 3), c(1, 1)),
               "must be square")
  expect_error(choicer:::equilibrated_inverse(diag(2), 1), "length 2")
  expect_error(choicer:::equilibrated_sandwich(diag(2), diag(3), c(1, 1)),
               "order of `A`")
  expect_error(choicer:::equilibrated_sandwich(cases$singular, diag(2), c(1, 1)),
               "exactly singular")
})

test_that("the equilibrated inverse treats non-finite matrices as solve() does", {
  # solve() estimates the condition only for finite factors; a non-finite
  # pivot's fate is the BLAS's (both call the same LU factorization).
  skip_on_cran()
  cases <- list(nan_pivot = diag(c(NaN, 1)),
                nan_offdiag = matrix(c(1, NaN, NaN, 1), 2),
                inf = diag(c(Inf, 1)), na = matrix(NA_real_, 1, 1))
  for (nm in names(cases)) {
    M <- cases[[nm]]
    expect_identical(
      vop_outcome(choicer:::equilibrated_inverse(M, rep(1, nrow(M)))),
      vop_outcome(solve(M)), label = nm)
  }
})

test_that("an equilibrated sandwich recovers a badly scaled A^-1 B A^-1", {
  bp <- vop_bad_pair(3L)
  set.seed(4)
  G <- matrix(rnorm(36), 6)
  B0 <- crossprod(G) / 6
  B <- B0 * tcrossprod(bp$s)
  Ainv <- solve(bp$A)
  V_true <- (Ainv %*% B0 %*% Ainv) / tcrossprod(bp$s)
  expect_error(solve(bp$H), "computationally singular")
  eq <- vop_msgs(choicer:::.sandwich_combine(bp$H, B))
  expect_length(eq$messages, 0L)
  expect_lt(vop_err(eq$value$vcov, V_true), 1e-12)
  expect_identical(eq$value$vcov, t(eq$value$vcov))
})

# Simulated mixed logit data with badly scaled columns, weights and cluster
# labels constant within decision maker.
vop_mxl_data <- function(N, T, seed, bad = FALSE) {
  sim <- simulate_mxl_data(N = N, J = 4L, T = T, beta = c(0.8, -0.5),
                           Sigma = diag(c(1, 0.6)), vary_choice_set = TRUE,
                           seed = seed)
  dt <- data.table::as.data.table(sim$data)
  if (bad) {
    dt[, `:=`(x1 = x1 * 1e7, x2 = x2 * 1e-7)]
  } else {
    dt[, `:=`(x1 = x1 * 10, x2 = x2 / 40, w1 = w1 * 5)]
  }
  unit <- if (T > 1L) dt$pid else dt$id
  u <- match(unit, unique(unit))
  dt[, `:=`(grp = (u - 1L) %% 7L, w = 0.5 + (u %% 5L) / 4)]
  dt
}

vop_mxl_fit <- function(dt, panel = FALSE, control = list(maxeval = 400L),
                        ...) {
  suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = c("x1", "x2"), random_var_cols = c("w1", "w2"),
    person_col = if (panel) "pid", S = 20L, seed = 3L, control = control,
    ...))
}

vop_type <- function(se) switch(se, sandwich = "robust", numeric = "hessian", se)

# The variance of the lazy route: vcov() after the fit's vcov was dropped.
vop_lazy <- function(fit) {
  fit["vcov"] <- list(NULL)
  fit["se"] <- list(NULL)
  suppressMessages(vcov(fit))
}

test_that("a mixed logit's variance is the one vcov() returns post hoc", {
  skip_on_cran()
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  # A design whose unscaled fits' matrices invert (with seed 11 the unscaled
  # BHHH is singular: a variance near zero). The
  # unscaled fits start where 0.2.x did: from the default, the unscaled
  # optimizer drives w2's variance onto its zero plateau, where the
  # information is singular. The layout: x1, x2; L_11, L_21, L_22; constants.
  dt <- vop_mxl_data(200L, 1L, 13L)
  start_0_2 <- c(0, 0, log(0.5), 0, log(0.5),
                 rep(0, data.table::uniqueN(dt$alt) - 1L))
  for (sv in c("none", "sd")) {
    for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
      fit <- vop_mxl_fit(dt, scale_vars = sv, se_method = se,
                         weights_col = if (se == "sandwich") "w",
                         rc_correlation = TRUE,
                         cluster_col = if (se == "cluster") "grp",
                         theta_init = if (sv == "none") start_0_2)
      what <- paste(sv, se)
      expect_false(is.null(fit$vcov), label = what)
      expect_identical(suppressMessages(vcov(fit, type = vop_type(se))),
                       fit$vcov, label = paste("post hoc", what))
      expect_identical(vop_lazy(fit), fit$vcov, label = paste("lazy", what))
    }
  }
  # A panel with generated draws (the fit's own draws' arguments) and a
  # log-normal coefficient.
  dtp <- vop_mxl_data(60L, 3L, 12L)
  for (se in c("hessian", "sandwich", "cluster")) {
    fit <- vop_mxl_fit(dtp, panel = TRUE, scale_vars = "mad",
                       draws = "generate", rc_dist = c(0L, 1L), rc_mean = TRUE,
                       se_method = se, cluster_col = if (se == "cluster") "grp")
    expect_identical(suppressMessages(vcov(fit, type = vop_type(se))),
                     fit$vcov, label = paste("panel", se))
  }
  # wesml_vcov() is the robust route.
  fit <- vop_mxl_fit(dt, scale_vars = "sd", se_method = "sandwich",
                     weights_col = "w")
  expect_identical(unname(wesml_vcov(fit)), unname(fit$vcov))
})

test_that("multinomial and nested logit variances follow their post-hoc routes", {
  skip_on_cran()
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  sim <- simulate_mnl_data(N = 300, J = 4, seed = 21, outside_option = FALSE)
  dm <- data.table::as.data.table(sim$data)
  dm[, x1 := x1 * 100]
  dm[, clus := (id - 1L) %/% 10L]
  for (sv in c("none", "sd")) {
    for (se in c("hessian", "bhhh", "sandwich", "cluster")) {
      fit <- suppressMessages(run_mnlogit(
        dm, "id", "alt", "choice", c("x1", "x2"), scale_vars = sv,
        se_method = se, cluster_col = if (se == "cluster") "clus"))
      post <- suppressMessages(vcov(fit, type = vop_type(se)))
      lazy <- vop_lazy(fit)
      what <- paste("MNL", sv, se)
      if (se %in% c("hessian", "cluster")) {
        expect_identical(post, fit$vcov, label = what)
      } else {
        # The fit's BHHH kernel against the post-hoc score matrix: the same
        # sums in another order.
        expect_equal(post, fit$vcov, tolerance = 1e-10, label = what)
      }
      # The lazy route (compute_hessian()) is the fit's kernel route, except
      # for the sandwich (the post-hoc route).
      if (se == "sandwich") {
        expect_equal(lazy, fit$vcov, tolerance = 1e-10, label = paste("lazy", what))
      } else {
        expect_identical(lazy, fit$vcov, label = paste("lazy", what))
      }
    }
  }
  simn <- simulate_nl_data(N = 400, seed = 22)
  dn <- data.table::as.data.table(simn$data)
  dn[, clus := (id - 1L) %/% 10L]
  for (sv in c("none", "bhhh")) {
    for (se in c("hessian", "numeric", "bhhh", "sandwich", "cluster")) {
      what <- paste("NL", sv, se)
      fit <- suppressMessages(run_nestlogit(
        dn, "id", "j", "choice", c("X", "W"), "nest", se_method = se,
        cluster_col = if (se == "cluster") "clus", scale_vars = sv))
      post <- suppressMessages(vcov(fit, type = vop_type(se)))
      if (se %in% c("hessian", "numeric", "cluster")) {
        expect_identical(post, fit$vcov, label = what)
      } else {
        expect_equal(post, fit$vcov, tolerance = 1e-10, label = what)
      }
    }
  }
})

test_that("scaled fits of badly scaled designs invert where solve() fails", {
  skip_on_cran()
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- vop_mxl_data(300L, 1L, 31L, bad = TRUE)
  fit <- vop_mxl_fit(dt, scale_vars = "sd", control = list(maxeval = 300L))
  expect_true(all(is.finite(fit$se)))
  post <- vop_msgs(vcov(fit, type = "hessian"))
  expect_length(post$messages, 0L)
  expect_identical(post$value, fit$vcov)
  # The natural information itself is singular to a raw solve().
  H <- choicer:::.compute_bread(fit)
  expect_error(solve(H), "singular")

  sim <- simulate_mnl_data(N = 300, J = 4, seed = 32, outside_option = FALSE)
  dm <- data.table::as.data.table(sim$data)
  dm[, `:=`(x1 = x1 * 1e7, x2 = x2 * 1e-7)]
  fm <- suppressMessages(run_mnlogit(dm, "id", "alt", "choice", c("x1", "x2"),
                                     scale_vars = "sd"))
  expect_true(all(is.finite(fm$se)))
  post <- vop_msgs(vcov(fm, type = "hessian"))
  expect_length(post$messages, 0L)
  expect_identical(post$value, fm$vcov)

  # The inversion does not depend on scale_vars: the same fit, marked
  # unscaled, inverts its information as the scaled one did, where solve()
  # refuses it.
  for (f in list(fit, fm)) {
    f0 <- f
    f0$scale_vars <- "none"
    f0["vcov"] <- list(NULL)
    f0["se"] <- list(NULL)
    lazy <- vop_msgs(vcov(f0))
    expect_length(lazy$messages, 0L)
    expect_identical(lazy$value, f$vcov)
    expect_error(solve(choicer:::.compute_bread(f0)), "singular")
  }
})

test_that("fit objects keep their shape and their own draws", {
  skip_on_cran()
  set_omp_threads(1L)
  on.exit(set_omp_threads(2L), add = TRUE)
  dt <- vop_mxl_data(200L, 1L, 41L)
  with_data <- vop_mxl_fit(dt, scale_vars = "sd")
  no_data <- vop_mxl_fit(dt, scale_vars = "sd", keep_data = FALSE)
  expect_identical(names(no_data), names(with_data))
  expect_null(no_data$data)
  expect_identical(no_data$vcov, with_data$vcov)
  expect_identical(unname(with_data$param_scale[with_data$param_map$beta]),
                   unname(1 / with_data$sX))
  expect_identical(names(with_data$param_shift), names(coef(with_data)))
  none <- vop_mxl_fit(dt, scale_vars = "none")
  expect_identical(unname(none$param_scale), rep(1, length(coef(none))))
  expect_identical(unname(none$param_shift), rep(0, length(coef(none))))

  # Objects without the new fields (fitted before them) still work: the
  # variance reads neither the map nor scale_vars.
  old <- with_data
  old$param_scale <- NULL
  old$param_shift <- NULL
  expect_identical(suppressMessages(vcov(old, type = "hessian")),
                   with_data$vcov)

  # A stored-draws fit builds its cube once: the variance reuses it.
  real <- get("get_halton_normals", envir = asNamespace("choicer"))
  n_cubes <- 0L
  local_mocked_bindings(get_halton_normals = function(...) {
    n_cubes <<- n_cubes + 1L
    real(...)
  })
  fs <- vop_mxl_fit(dt, se_method = "sandwich", scale_vars = "sd")
  expect_identical(n_cubes, 1L)

  # The advanced workflow's own draws give the fit's variance (post hoc
  # regenerates Halton draws instead).
  d_adv <- prepare_mxl_data(dt, "id", "alt", "choice", c("x1", "x2"),
                            c("w1", "w2"))
  set.seed(42)
  eta <- array(rnorm(2L * 15L * d_adv$N), dim = c(2L, 15L, d_adv$N))
  fa <- suppressMessages(run_mxlogit(input_data = d_adv, eta_draws = eta,
                                     scale_vars = "sd",
                                     control = list(maxeval = 200L)))
  H <- choicer:::mxl_hessian_parallel(
    theta = coef(fa), X = d_adv$X, W = d_adv$W, alt_idx = d_adv$alt_idx,
    choice_idx = d_adv$choice_idx, M = d_adv$M, weights = d_adv$weights,
    eta_draws = eta, rc_dist = c(0L, 0L), rc_correlation = FALSE,
    rc_mean = FALSE, use_asc = TRUE, include_outside_option = FALSE)
  expect_identical(unname(fa$vcov),
                   unname(choicer:::invert_hessian(H)$vcov))
})
