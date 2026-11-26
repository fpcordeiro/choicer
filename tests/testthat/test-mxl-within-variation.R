# Covariates that vary within no choice situation (prepare_mxl_data()).
#
# A covariate with one value on all of a situation's alternatives, in every
# situation (a decision maker's age entered alone), adds the same amount to
# every alternative's utility in each situation. Without an outside option it
# cancels from every utility difference and is not identified: the preparation
# drops it with a warning (class choicer_unidentified_covariates). With an
# outside option it moves the inside alternatives against it and is kept, with
# a message. Situations are compared one by one, so in a panel a covariate
# that changes between a decision maker's situations but within none of them
# is dropped too.

# N decision makers with T situations each (T = 1: a cross-section) and J
# alternatives, in prepared order. `age` is the decision maker's, `inc`
# changes between a decision maker's situations, `x` and `w` vary over
# alternatives. Choices are drawn from a mixed logit in which age, inc and x
# shift the inside alternatives against the outside option and w's
# coefficient varies over decision makers; without an outside option, age and
# inc cancel.
wv_data <- function(seed, N = 60L, T = 1L, J = 3L, outside = FALSE) {
  set.seed(seed)
  n_sit <- N * T
  d <- data.frame(person = rep(seq_len(N), each = T * J),
                  id = rep(seq_len(n_sit), each = J),
                  alt = rep(seq_len(J), n_sit))
  d$age <- rep(stats::runif(N, 20, 80), each = T * J)
  d$inc <- rep(stats::rnorm(n_sit), each = J)
  d$x <- stats::rnorm(nrow(d))
  d$w <- stats::rnorm(nrow(d))
  b_w <- rep(stats::rnorm(N, 1, 0.5), each = T * J)
  v <- 1 + d$x - 0.02 * d$age + 0.5 * d$inc + b_w * d$w
  u <- matrix(v - log(-log(stats::runif(nrow(d)))), J)
  pick <- apply(u, 2L, which.max)
  if (outside) pick[apply(u, 2L, max) < -log(-log(stats::runif(n_sit)))] <- 0L
  d$choice <- as.integer(d$alt == rep(pick, each = J))
  d
}

# The value of `expr`, with the messages and warnings it signalled, muffled
wv_collect <- function(expr) {
  msgs <- character()
  warns <- list()
  value <- withCallingHandlers(expr,
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns[[length(warns) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  list(value = value, messages = msgs, warnings = warns)
}

# The warnings of class choicer_unidentified_covariates among those collected
wv_unidentified <- function(res) {
  Filter(function(w) inherits(w, "choicer_unidentified_covariates"),
         res$warnings)
}

test_that("design_varies_within() finds the columns that vary within a situation", {
  # The reference compares every row but a situation's first with its first
  # row, a NaN differing from every value
  ref <- function(X, M) {
    first <- rep(c(0L, cumsum(M)[-length(M)]) + 1L, M)
    later <- seq_len(nrow(X)) != first
    vapply(seq_len(ncol(X)), function(k) {
      d <- X[later, k] != X[first[later], k]
      any(is.na(d) | d)
    }, logical(1L))
  }
  set.seed(1)
  M <- c(3L, 1L, 4L, 2L)
  sit <- rep(seq_along(M), M)
  n <- sum(M)
  last <- stats::rnorm(length(M))[sit]
  last[n] <- last[n] + 1  # varies within the last situation only
  X <- cbind(everywhere = stats::rnorm(n),
             situation = stats::rnorm(length(M))[sit],
             last = last,
             constant = rep(5, n),
             signed_zero = c(0, -0, 0, 1, -0, 0, 0, 0, 2, 2),
             nan = c(1, NaN, 1, 2, 3, 3, 3, 3, 4, 4))
  expect_identical(design_varies_within(X, M),
                   c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))
  expect_identical(design_varies_within(X, M), ref(X, M))
  # A situation of one alternative never varies, and an empty one is skipped
  expect_identical(design_varies_within(X[1:3, , drop = FALSE], c(1L, 2L, 0L)),
                   c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(design_varies_within(X[, 0L, drop = FALSE], M), logical())
  expect_identical(design_varies_within(X[0L, , drop = FALSE], integer()),
                   rep(FALSE, ncol(X)))
  # Random designs: situations of one to five alternatives, columns that vary
  # within each situation with some probability
  for (r in 1:20) {
    M <- sample.int(5L, sample.int(30L, 1L), replace = TRUE)
    sit <- rep(seq_along(M), M)
    X <- vapply(c(0, 0.02, 0.2, 1), function(p) {
      x <- stats::rnorm(length(M))[sit]
      bump <- stats::runif(length(sit)) < p
      x[bump] <- x[bump] + 1
      x
    }, numeric(length(sit)))
    expect_identical(design_varies_within(X, M), ref(X, M))
  }
  expect_error(design_varies_within(X, as.numeric(M)), "integer vector")
  expect_error(design_varies_within(X, M[-1L]), "sum to the rows")
  expect_error(design_varies_within(X, c(M[-1L], NA)), "situations' sizes")
  expect_error(design_varies_within(X, c(M, -1L)), "situations' sizes")
  expect_error(design_varies_within(1:3, 3L), "double matrix")
})

test_that("such a covariate is dropped without an outside option, situation by situation", {
  # In covariate_cols and in random_var_cols, in a cross-section (age) and in
  # a panel (inc, which changes between a decision maker's situations and so
  # would pass a check by decision maker). What is left is the preparation
  # without the covariate, bit for bit.
  for (T in c(1L, 3L)) for (role in c("covariate_cols", "random_var_cols")) {
    d <- wv_data(40L + T, T = T)
    person <- if (T > 1L) "person"
    col <- if (T > 1L) "inc" else "age"
    lab <- sprintf("[T = %d, %s]", T, role)
    base <- list(covariate_cols = "x", random_var_cols = "w")
    cols <- base
    cols[[role]] <- c(cols[[role]], col)
    res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice",
                                       cols$covariate_cols, cols$random_var_cols,
                                       person_col = person))
    w <- wv_unidentified(res)
    expect_length(w, 1L)
    expect_match(conditionMessage(w[[1L]]),
                 sprintf("dropped: %s \\(%s\\)\\.", col, role), label = lab)
    expect_length(res$warnings, 1L)
    expect_length(res$messages, 0L)
    a <- res$value
    expect_identical(a$dropped_cols, col, label = lab)
    # The positions in random_var_cols of the columns dropped from it
    expect_identical(a$random_var_dropped,
                     if (role == "random_var_cols") {
                       list(within = 2L, collinear = integer())
                     }, label = lab)
    b <- prepare_mxl_data(d, "id", "alt", "choice", base$covariate_cols,
                          base$random_var_cols, person_col = person)
    expect_null(b$random_var_dropped)
    same <- setdiff(names(b), c("dropped_cols", "data_spec"))
    expect_identical(a[same], b[same], label = lab)
    expect_identical(a$data_spec[[role]], cols[[role]], label = lab)
  }
})

test_that("a covariate that varies within one situation is kept", {
  for (T in c(1L, 3L)) {
    d <- wv_data(50L + T, T = T)
    d$age[2L] <- d$age[2L] + 1  # the first situation's second alternative
    person <- if (T > 1L) "person"
    res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice", c("x", "age"),
                                       c("w", "age"), person_col = person))
    expect_length(res$warnings, 0L)
    expect_length(res$messages, 0L)
    expect_identical(colnames(res$value$X), c("x", "age"))
    expect_identical(colnames(res$value$W), c("w", "age"))
    expect_identical(res$value$dropped_cols, character())
  }
})

test_that("with an outside option such a covariate is kept, with a message", {
  for (T in c(1L, 3L)) for (role in c("covariate_cols", "random_var_cols")) {
    d <- wv_data(60L + T, T = T, outside = TRUE)
    person <- if (T > 1L) "person"
    col <- if (T > 1L) "inc" else "age"
    lab <- sprintf("[T = %d, %s]", T, role)
    cols <- list(covariate_cols = "x", random_var_cols = "w")
    cols[[role]] <- c(cols[[role]], col)
    res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice",
                                       cols$covariate_cols, cols$random_var_cols,
                                       include_outside_option = TRUE,
                                       person_col = person))
    expect_length(res$warnings, 0L)
    expect_length(res$messages, 1L)
    expect_match(res$messages, sprintf("kept: %s \\(%s\\)\\.", col, role),
                 label = lab)
    a <- res$value
    expect_identical(colnames(a$X), cols$covariate_cols, label = lab)
    expect_identical(colnames(a$W), cols$random_var_cols, label = lab)
    expect_identical(a$dropped_cols, character(), label = lab)
  }
})

test_that("the preparation stops when no random coefficient remains", {
  d <- wv_data(71L)
  w <- NULL
  expect_error(withCallingHandlers(
    prepare_mxl_data(d, "id", "alt", "choice", "x", c("age", "inc")),
    choicer_unidentified_covariates = function(cnd) {
      w <<- cnd
      invokeRestart("muffleWarning")
    }), "No random coefficient remains")
  expect_match(conditionMessage(w),
               "dropped: age \\(random_var_cols\\), inc \\(random_var_cols\\)\\.")
  # Two columns of zeros vary within no situation either. With an outside
  # option that check keeps them, and collinearity leaves none, which reached
  # the fit and stopped it on K_w
  d$z1 <- 0
  d$z2 <- 0
  expect_error(suppressWarnings(prepare_mxl_data(d, "id", "alt", "choice", "x",
                                                 c("z1", "z2"))),
               "No random coefficient remains")
  expect_error(suppressMessages(prepare_mxl_data(d, "id", "alt", "choice", "x",
                                                 c("z1", "z2"),
                                                 include_outside_option = TRUE)),
               "No random coefficient remains")
  # Every fixed coefficient may go: a mixed logit needs none, and the
  # preparation is the one without fixed coefficients
  res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice", "age", "w"))
  expect_identical(res$value$dropped_cols, "age")
  expect_length(wv_unidentified(res), 1L)
  b <- prepare_mxl_data(d, "id", "alt", "choice", NULL, "w")
  same <- setdiff(names(b), c("dropped_cols", "data_spec"))
  expect_identical(res$value[same], b[same])
})

test_that("a shuffled panel is checked in prepared order", {
  # The rows come in any order; the check reads the design the preparation
  # sorted, so the result is the sorted data's
  d <- wv_data(73L, T = 3L)
  s <- d[sample.int(nrow(d)), ]
  a <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice", c("x", "inc"),
                                   c("w", "age"), person_col = "person"))
  b <- wv_collect(prepare_mxl_data(s, "id", "alt", "choice", c("x", "inc"),
                                   c("w", "age"), person_col = "person"))
  expect_identical(b$value$dropped_cols, c("inc", "age"))
  expect_identical(b$value[c("X", "W", "M", "Ti", "choice_idx", "alt_idx")],
                   a$value[c("X", "W", "M", "Ti", "choice_idx", "alt_idx")])
  # One situation's covariate varies: kept, wherever its rows fall
  s$inc[which(s$id == s$id[1L])[1L]] <- 99
  b <- wv_collect(prepare_mxl_data(s, "id", "alt", "choice", c("x", "inc"),
                                   c("w", "age"), person_col = "person"))
  expect_identical(b$value$dropped_cols, "age")
})

test_that("bounds named by the Cholesky factor stop after a random coefficient's column is dropped", {
  # L_<i><j> number the kept columns: dropping age would move L_22 from w1
  # to w2. A bound named by a column (Mu_w1) keeps its meaning.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- wv_data(74L, N = 80L)
  d$w1 <- d$w
  d$w2 <- stats::rnorm(nrow(d))
  fit <- function(rv, ...) {
    wv_collect(run_mxlogit(data = d, id_col = "id", alt_col = "alt",
                           choice_col = "choice", covariate_cols = "x",
                           random_var_cols = rv, S = 20L, ...))$value
  }
  expect_error(fit(c("age", "w1", "w2"), upper = c(L_22 = log(0.05))),
               "Bounds named by the Cholesky factor \\(L_22\\) are ambiguous")
  expect_error(fit(c("w1", "age", "w2"), lower = c(L_11 = -5)),
               "dropped age from")
  f1 <- fit(c("age", "w1", "w2"), rc_mean = TRUE, upper = c(Mu_w1 = 0))
  f0 <- fit(c("w1", "w2"), rc_mean = TRUE, upper = c(Mu_w1 = 0))
  expect_identical(coef(f1), coef(f0))
  # Without a drop, as before: L_22 is w2's, which moves no choice, so its
  # bound binds
  f2 <- fit(c("w1", "w2"), lower = c(L_22 = 0))
  expect_gte(coef(f2)[["L_22"]], 0)
  # The preparation's record of its drops decides, by position, not W's
  # names: without them, a W with nothing dropped keeps its named bounds, one
  # with a drop does not, and rc_dist still follows the kept columns
  a <- prepare_mxl_data(d, "id", "alt", "choice", "x", "w1")
  colnames(a$W) <- NULL
  f3 <- wv_collect(run_mxlogit(input_data = a,
                               eta_draws = get_halton_normals(20L, a$N, 1L),
                               lower = c(L_11 = 0)))$value
  expect_gte(coef(f3)[["L_11"]], 0)
  a <- suppressWarnings(prepare_mxl_data(d, "id", "alt", "choice", "x",
                                         c("age", "w1", "w2")))
  colnames(a$W) <- NULL
  draws <- get_halton_normals(20L, a$N, 2L)
  expect_error(run_mxlogit(input_data = a, eta_draws = draws,
                           upper = c(L_22 = 0)),
               "dropped age from `random_var_cols`, so L_<i><j> number the columns kept \\(w1, w2\\)")
  f4 <- wv_collect(run_mxlogit(input_data = a, eta_draws = draws,
                               rc_dist = c(0L, 1L, 0L)))$value
  expect_identical(f4$rc_dist, c(1L, 0L))
})

test_that("with an outside option and one alternative per situation there is no message", {
  # A binary choice against the outside option: every covariate has one value
  # in each situation, and each moves its alternative against the outside
  d <- wv_data(75L, J = 1L, outside = TRUE)
  res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice", c("x", "age"),
                                     "w", include_outside_option = TRUE))
  expect_length(res$messages, 0L)
  expect_length(res$warnings, 0L)
  expect_identical(colnames(res$value$X), c("x", "age"))
})

test_that("the warning names at most ten columns", {
  d <- wv_data(72L)
  extra <- sprintf("c%02d", 1:12)
  for (k in extra) d[[k]] <- rep(stats::rnorm(60L), each = 3L)
  res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice", c("x", extra), "w"))
  msg <- conditionMessage(wv_unidentified(res)[[1L]])
  expect_match(msg, "c10 \\(covariate_cols\\) and 2 more\\.")
  expect_false(grepl("c11", msg, fixed = TRUE))
  expect_identical(res$value$dropped_cols, extra)
})

test_that("run_mxlogit() fits the identified model after the drop, rc_dist following the kept columns", {
  # A fixed coefficient's covariate and a random coefficient's column dropped
  # together, with rc_dist given for random_var_cols: the fit is the one
  # without them, bit for bit at one thread, in a cross-section and a panel.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (T in c(1L, 3L)) {
    d <- wv_data(80L + T, N = 100L, T = T)
    person <- if (T > 1L) "person"
    lab <- sprintf("[T = %d]", T)
    fit <- function(cov, rv, rc_dist) {
      wv_collect(run_mxlogit(data = d, id_col = "id", alt_col = "alt",
                             choice_col = "choice", covariate_cols = cov,
                             random_var_cols = rv, rc_dist = rc_dist, S = 20L,
                             person_col = person))
    }
    res <- fit(c("x", "age"), c("inc", "w"), c(0L, 1L))
    expect_length(wv_unidentified(res), 1L)
    f1 <- res$value
    f0 <- fit("x", "w", 1L)$value
    expect_identical(names(coef(f1)), names(coef(f0)), label = lab)
    expect_identical(coef(f1), coef(f0), label = lab)
    expect_identical(f1$se, f0$se, label = lab)
    expect_identical(f1$loglik, f0$loglik, label = lab)
    expect_identical(f1$rc_dist, 1L, label = lab)
    expect_true(all(is.finite(f1$se)), label = lab)
  }
})

test_that("with an outside option such a covariate is estimated", {
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  for (T in c(1L, 3L)) {
    d <- wv_data(90L + T, N = 300L, T = T, outside = TRUE)
    person <- if (T > 1L) "person"
    col <- if (T > 1L) "inc" else "age"
    res <- wv_collect(run_mxlogit(
      data = d, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = c("x", col), random_var_cols = "w",
      include_outside_option = TRUE, S = 20L, person_col = person))
    expect_match(res$messages, sprintf("kept: %s \\(covariate_cols\\)\\.", col),
                 all = FALSE)
    expect_length(wv_unidentified(res), 0L)
    se <- res$value$se
    expect_true(all(is.finite(se)), label = sprintf("[T = %d]", T))
  }
})

test_that("after a collinearity drop rc_dist takes one entry per kept column", {
  # w2 = 2 * w: the coefficient on w absorbs w2's, so an rc_dist for both
  # columns (a log-normal draw plus twice a normal one) asks for a
  # coefficient that neither distribution gives. It stops, as it did at the
  # first kernel call, now naming the kept columns, in both workflows. One
  # entry per kept column fits the model without w2, as does the default.
  on.exit(set_omp_threads(2L), add = TRUE)
  set_omp_threads(1L)
  d <- wv_data(95L, N = 100L)
  d$w2 <- 2 * d$w
  fit <- function(rv, rc_dist) {
    wv_collect(run_mxlogit(
      data = d, id_col = "id", alt_col = "alt", choice_col = "choice",
      covariate_cols = "x", random_var_cols = rv, rc_dist = rc_dist,
      S = 20L))$value
  }
  stops <- "dropped w2 as collinear with the others.*one entry per column kept \\(w\\)"
  expect_error(suppressMessages(fit(c("w", "w2"), c(1L, 0L))), stops)
  f1 <- fit(c("w", "w2"), 1L)
  f0 <- fit("w", 1L)
  expect_identical(colnames(f1$data$W), "w")
  expect_identical(f1$rc_dist, 1L)
  expect_identical(coef(f1), coef(f0))
  expect_identical(f1$se, f0$se)
  expect_identical(coef(fit(c("w", "w2"), NULL)), coef(fit("w", NULL)))
  a <- suppressMessages(prepare_mxl_data(d, "id", "alt", "choice", "x",
                                         c("w", "w2")))
  expect_identical(a$random_var_dropped,
                   list(within = integer(), collinear = 2L))
  draws <- get_halton_normals(20L, 100L, 1L)
  expect_error(run_mxlogit(input_data = a, eta_draws = draws,
                           rc_dist = c(1L, 0L)), stops)
  f2 <- wv_collect(run_mxlogit(input_data = a, eta_draws = draws,
                               rc_dist = 1L))$value
  expect_identical(coef(f2), coef(f0))
  expect_identical(f2$se, f0$se)
})

test_that("such columns go first, and collinearity is judged within situations", {
  # A task's start time (one value per task, a timestamp) and its end, or two
  # end times: on raw levels, at their scale, the QR finds them collinear and
  # drops the later column. The start goes first, and the collinearity check
  # takes each row less its situation's first row, so every end time is kept,
  # in either order, next to another random coefficient too.
  d <- wv_data(96L)
  d$start <- 1.6e9 + d$id
  d$end <- d$start + stats::runif(nrow(d), 0, 10)
  d$end_b <- d$start + stats::runif(nrow(d), 0, 10)
  qr_kept <- function(cols) {
    colnames(suppressMessages(check_collinearity(as.matrix(d[cols])))$mat)
  }
  expect_identical(qr_kept(c("start", "end")), "start")
  expect_identical(qr_kept(c("end", "start")), "end")
  expect_identical(qr_kept(c("end", "end_b")), "end")
  for (rv in list(c("start", "end"), c("end", "start"), c("w", "start", "end"),
                  c("start", "end", "end_b"))) {
    lab <- paste(rv, collapse = ", ")
    res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice", "x", rv))
    expect_identical(colnames(res$value$W), setdiff(rv, "start"), label = lab)
    expect_identical(res$value$dropped_cols, "start", label = lab)
    expect_length(wv_unidentified(res), 1L)
    expect_length(res$messages, 0L)  # nothing collinear is left
    expect_identical(res$value$random_var_dropped,
                     list(within = match("start", rv), collinear = integer()),
                     label = lab)
  }
  # x and x + age vary alike within every situation: one goes, in either
  # order and in either role, where raw levels kept both. So does one of a
  # set of dummies for every alternative, which sums to one in each row.
  d$xa <- d$x + d$age
  d$wa <- d$w + d$age
  for (cov in list(c("age", "x", "xa"), c("x", "xa", "age"), c("xa", "age", "x"))) {
    lab <- paste(cov, collapse = ", ")
    res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice", cov, "w"))
    kept <- setdiff(cov, "age")
    expect_identical(colnames(res$value$X), kept[1L], label = lab)
    expect_identical(res$value$dropped_cols, c("age", kept[2L]), label = lab)
    expect_length(res$messages, 1L)
    expect_match(res$messages, paste0("dropped due to collinearity of their ",
                                      "variation within choice situations: ",
                                      kept[2L], "[[:space:]]*$"), label = lab)
  }
  res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice", "x",
                                     c("w", "age", "wa")))
  expect_identical(colnames(res$value$W), "w")
  expect_identical(res$value$random_var_dropped,
                   list(within = 2L, collinear = 3L))
  for (j in 1:3) d[[paste0("d", j)]] <- as.numeric(d$alt == j)
  res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice",
                                     c("x", "d1", "d2", "d3"), "w"))
  expect_identical(colnames(res$value$X), c("x", "d1", "d2"))
  expect_identical(res$value$dropped_cols, "d3")
  # With an outside option the check comes after the QR, so that its message
  # names only columns kept
  d$zero <- 0
  res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice",
                                     c("x", "zero", "age"), "w",
                                     include_outside_option = TRUE))
  expect_identical(colnames(res$value$X), c("x", "age"))
  expect_identical(res$value$dropped_cols, "zero")
  expect_match(res$messages, "dropped due to collinearity: zero", all = FALSE)
  kept <- grep("kept:", res$messages, value = TRUE)
  expect_length(kept, 1L)
  expect_match(kept, "kept: age \\(covariate_cols\\)\\. ")
})

test_that("the QR by row chunks judges variation within groups as in memory", {
  # remove_nullspace_cols() factors a million elements or more by row chunks;
  # with `first` each chunk is less its rows' groups' first rows, read from
  # the matrix wherever the chunk boundaries fall. Its rank and pivot are
  # those of qr() on the whole matrix less its first rows.
  set.seed(99)
  M <- sample(2:5, 40L, replace = TRUE)
  g <- rep(seq_along(M), M)
  first <- rep.int(cumsum(c(1, M[-length(M)])), M)
  n <- sum(M)
  m <- cbind(a = stats::rnorm(n), lvl = stats::rnorm(length(M))[g],
             b = stats::rnorm(n))
  m <- cbind(m, ab = m[, "a"] + m[, "lvl"], big = 1e6 + m[, "lvl"] + m[, "b"])
  # design_rows_less_first() writes a block of rows less their groups' first
  # rows, which may lie before the block
  expect_identical(design_rows_less_first(m, 5, 30, first),
                   unname(m[5:30, ] - m[first[5:30], ]))
  expect_identical(design_rows_less_first(m, 1, n, first),
                   unname(m - m[first, ]))
  expect_error(design_rows_less_first(m, 0, 3, first), "rows of `X`")
  expect_error(design_rows_less_first(m, 3, n + 1, first), "rows of `X`")
  expect_error(design_rows_less_first(m, 4, 3, first), "rows of `X`")
  expect_error(design_rows_less_first(m, 1, 3, as.integer(first)),
               "double vector")
  bad <- first
  bad[2L] <- n + 1
  expect_error(design_rows_less_first(m, 1, 3, bad), "hold rows of `X`")
  ref <- qr(m - m[first, , drop = FALSE], tol = 1e-7)
  expect_identical(ref$rank, 2L)
  for (rows in c(7L, 23L, n)) {
    q <- .qr_rank_pivot(m, 1e-7, rows = rows, first = first)
    expect_identical(q$rank, ref$rank, label = rows)
    expect_identical(q$pivot, ref$pivot, label = rows)
  }
  # Within groups only a and b vary on their own; on raw levels lvl and big
  # are columns of their own, and ab is a + lvl
  expect_identical(colnames(remove_nullspace_cols(m, first = first)),
                   c("a", "b"))
  expect_identical(colnames(remove_nullspace_cols(m)), c("a", "lvl", "b", "big"))
})

test_that("the record of the drops is read only while W has the columns kept", {
  # Reordered or replaced columns of a prepared W: the record no longer says
  # which coefficient is which, so an rc_dist for every column of
  # random_var_cols is refused, as without a record
  d <- wv_data(99L)
  d$w2 <- stats::rnorm(nrow(d))
  d$z <- stats::rnorm(nrow(d))
  a <- suppressWarnings(prepare_mxl_data(d, "id", "alt", "choice", "x",
                                         c("age", "w", "w2")))
  draws <- get_halton_normals(5L, a$N, 2L)
  b <- a
  b$W <- b$W[, c("w2", "w")]
  expect_error(run_mxlogit(input_data = b, eta_draws = draws,
                           rc_dist = c(0L, 1L, 0L)), "length K_w")
  b <- a
  b$W[, "w2"] <- d$z  # wv_data() rows are in prepared order
  colnames(b$W)[2L] <- "z"
  expect_error(run_mxlogit(input_data = b, eta_draws = draws,
                           rc_dist = c(0L, 1L, 0L)), "length K_w")
})

test_that("a column named twice in one list is an error", {
  # Collinearity used to drop the repeat without a message, and rc_dist and
  # the bounds named by the Cholesky factor could not say which was meant. A
  # column may be in both lists: a fixed mean and a random deviation from it.
  d <- wv_data(97L)
  d$z <- stats::rnorm(nrow(d))
  expect_error(prepare_mxl_data(d, "id", "alt", "choice", "x", c("w", "w", "z")),
               "`random_var_cols` names a column more than once: w\\.")
  expect_error(prepare_mxl_data(d, "id", "alt", "choice", c("x", "z", "x"), "w"),
               "`covariate_cols` names a column more than once: x\\.")
  a <- prepare_mxl_data(d, "id", "alt", "choice", c("x", "w"), "w")
  expect_identical(colnames(a$X), c("x", "w"))
  expect_identical(colnames(a$W), "w")
})

test_that("the outside-option message qualifies what identifies a random coefficient", {
  # With constants a covariate's variation between situations identifies a
  # fixed coefficient or a random coefficient's mean only, while a random
  # coefficient's spread also shows in the correlation of the inside
  # alternatives' utilities
  d <- wv_data(98L, outside = TRUE)
  res <- wv_collect(prepare_mxl_data(d, "id", "alt", "choice", "x",
                                     c("w", "age"),
                                     include_outside_option = TRUE))
  expect_match(res$messages,
               "identifies a fixed coefficient or a random coefficient's mean; a random coefficient's spread is also identified",
               fixed = TRUE)
})
