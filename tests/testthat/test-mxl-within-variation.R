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
    b <- prepare_mxl_data(d, "id", "alt", "choice", base$covariate_cols,
                          base$random_var_cols, person_col = person)
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
  # So does one that collinearity leaves with none (two columns of zeros),
  # which reached the fit and stopped it on K_w
  d$z1 <- 0
  d$z2 <- 0
  expect_error(suppressMessages(prepare_mxl_data(d, "id", "alt", "choice", "x",
                                                 c("z1", "z2"))),
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

test_that("rc_dist follows the kept columns after a collinearity drop", {
  # A random coefficient's column dropped for collinearity used to leave an
  # rc_dist given for random_var_cols one entry too long, which the first
  # kernel call refused; in both workflows it now follows the kept columns
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
  f1 <- fit(c("w", "w2"), c(1L, 0L))
  f0 <- fit("w", 1L)
  expect_identical(colnames(f1$data$W), "w")
  expect_identical(f1$rc_dist, 1L)
  expect_identical(coef(f1), coef(f0))
  expect_identical(f1$se, f0$se)
  a <- suppressMessages(prepare_mxl_data(d, "id", "alt", "choice", "x",
                                         c("w", "w2")))
  f2 <- wv_collect(run_mxlogit(input_data = a,
                               eta_draws = get_halton_normals(20L, 100L, 1L),
                               rc_dist = c(1L, 0L)))$value
  expect_identical(f2$rc_dist, 1L)
  expect_identical(coef(f2), coef(f0))
  expect_identical(f2$se, f0$se)
})
