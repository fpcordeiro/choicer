# Post-estimation and the WESML helpers must not depend on the names of the
# user's columns. They used to look up function locals inside `[.data.table`,
# where a column of the same name takes their place (`am` and `pos` in
# prepare_newdata(), `spec` in its outside-option filter, `keep_ids` and
# `choice_col` in the sampling helpers, `d` and `am` in the hierarchical Bayes
# label lookups), and to read a table grouped by the id column back under
# data.table's fixed name `V1`. Each test renames one input column at a time
# to such a name and expects the same results back.

# `data` with column `from` renamed to `to`, and `args` with every reference
# to `from` renamed too.
rename_input <- function(data, args, from, to) {
  d <- copy(data)
  setnames(d, from, to)
  args <- lapply(args, function(a) {
    if (is.character(a)) replace(a, a == from, to) else a
  })
  list(data = d, args = args)
}

# Rename each role column to each name in turn and collect, as
# "role -> name: what", every way the result of `check(data, args, from, to)`
# (a named logical vector) fails; errors are reported too.
rename_failures <- function(data, args, renames, check) {
  failures <- character(0)
  for (role in names(renames)) {
    for (to in renames[[role]]) {
      r <- rename_input(data, args, role, to)
      bad <- tryCatch({
        ok <- check(r$data, r$args, role, to)
        names(ok)[!ok]
      }, error = function(e) paste("error:", conditionMessage(e)))
      if (length(bad) > 0) {
        failures <- c(failures, sprintf("%s -> %s: %s", role, to,
                                        paste(bad, collapse = ", ")))
      }
    }
  }
  failures
}

# Ids are read back from tables grouped by the id column (and, in a panel, by
# the decision maker); covariates met the locals and arguments of
# prepare_newdata(), and `inside` is the name its outside-option filter now
# uses.
id_clash_names <- c("N", "V1")
covariate_clash_names <- c("pos", "am", "spec", "inside", "id_col", "alt_col")

# Fit `fit_fun` on the renamed data and check that predict(), logsum() and
# consumer_surplus() on that data as `newdata` reproduce the stored-data path,
# and that the predictions reproduce those of the fit `ref` on the original
# names.
newdata_check <- function(fit_fun, ref, ref_data, price_var) {
  ref_probs <- predict(ref, newdata = ref_data)
  function(data, args, from, to) {
    fit <- suppressMessages(do.call(fit_fun, c(list(data), args)))
    price <- if (price_var == from) to else price_var
    probs <- predict(fit, newdata = data)
    c(
      probabilities = identical(probs, predict(fit)),
      shares = identical(predict(fit, type = "shares", newdata = data),
                         predict(fit, type = "shares")),
      logsum = identical(logsum(fit, newdata = data), logsum(fit)),
      surplus = identical(consumer_surplus(fit, price, newdata = data),
                          consumer_surplus(fit, price)),
      reference = identical(probs, ref_probs)
    )
  }
}

# On one thread every fit and prediction is reproducible bit for bit, so a
# refit on renamed columns must match the reference exactly. On two, the
# kernels' reductions can differ in the last bits from run to run, which
# moves MXL shares and is enough to stop an NL refit elsewhere.
expect_newdata_immune <- function(fit_fun, data, args, renames, price_var) {
  set_num_threads(1L)
  on.exit(set_num_threads(2L), add = TRUE)
  ref <- suppressMessages(do.call(fit_fun, c(list(data), args)))
  check <- newdata_check(fit_fun, ref, data, price_var)
  expect_identical(rename_failures(data, args, renames, check), character(0))
}

test_that("MNL predictions on newdata are unaffected by column names", {
  expect_newdata_immune(
    run_mnlogit, create_small_mnl_data(),
    list(id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = c("x1", "x2")),
    list(id = id_clash_names, x2 = covariate_clash_names), price_var = "x1"
  )
})

test_that("MNL predictions with an outside option are unaffected by column names", {
  # Alternative 0 becomes an explicit outside row, which the newdata path
  # drops with the filter that used to evaluate `spec` inside the data.
  dt <- create_small_mnl_data()
  dt[, alt := alt - 1L]
  dt[alt == 0L, c("x1", "x2") := 0]
  expect_newdata_immune(
    run_mnlogit, dt,
    list(id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = c("x1", "x2"), outside_opt_label = 0L,
         include_outside_option = TRUE),
    list(id = id_clash_names, x2 = covariate_clash_names), price_var = "x1"
  )
})

test_that("MXL predictions on newdata are unaffected by column names", {
  dt <- create_small_mxl_data()
  args <- list(id_col = "id", alt_col = "alt", choice_col = "choice",
               covariate_cols = "x1", random_var_cols = c("w1", "w2"),
               S = 10L)
  expect_newdata_immune(
    run_mxlogit, dt, args,
    list(id = id_clash_names, w2 = covariate_clash_names), price_var = "x1"
  )
  # Panel: ten decision makers with three situations each.
  dt[, pid := (id - 1L) %/% 3L + 1L]
  expect_newdata_immune(
    run_mxlogit, dt, c(args, list(person_col = "pid")),
    list(id = id_clash_names, pid = id_clash_names), price_var = "x1"
  )
})

test_that("NL predictions on newdata are unaffected by column names", {
  expect_newdata_immune(
    run_nestlogit, create_small_nl_data(),
    list(id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = c("x1", "x2"), nest_col = "nest"),
    list(id = id_clash_names, x2 = covariate_clash_names), price_var = "x1"
  )
})

test_that("a covariate named pos that permutes 1..J keeps the alternative codes", {
  # The reported silent case: prepare_newdata() coded each row by the
  # covariate `pos` in place of its position in alt_mapping, so a covariate
  # that permutes 1..J within every situation (a rank, say) passed every
  # check and mislabelled the rows.
  dt <- create_small_mnl_data()
  dt[, pos := as.numeric(sample(.N)), by = id]
  fit <- suppressMessages(
    run_mnlogit(dt, "id", "alt", "choice", c("x1", "pos"))
  )
  expect_identical(predict(fit, newdata = dt), predict(fit))
  expect_identical(predict(fit, type = "shares", newdata = dt),
                   predict(fit, type = "shares"))
})

# A population of 60 choice situations over three named alternatives, ids
# not 1..N, and an unused column `u` (unaligned with the ids) to rename to the
# helpers' local names.
wesml_fixture <- function() {
  set.seed(31)
  dt <- data.table(id = rep(seq_len(60L) * 10L, each = 3L),
                   alt = rep(c("bus", "car", "rail"), 60L))
  dt[, x1 := round(rnorm(.N), 3)]
  dt[, u := round(runif(.N), 3)]
  dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
  dt[]
}

# The three WESML results on renamed data, with the renamed column's name
# restored so they can be compared with the reference.
wesml_results <- function(data, args, from, to, Q, n_per_alt) {
  outputs <- list(
    weights = do.call(wesml_weights, c(list(data), args, list(Q = Q))),
    attached = do.call(wesml_weights,
                       c(list(data), args, list(Q = Q, attach = TRUE))),
    sample = do.call(sample_by_choice,
                     c(list(data), args, list(n_per_alt = n_per_alt,
                                              seed = 1L)))
  )
  lapply(outputs, function(o) {
    if (to %in% names(o)) setnames(o, to, from)
    o
  })
}

expect_wesml_immune <- function(data, args, renames, Q, n_per_alt) {
  ref <- wesml_results(data, args, "id", "id", Q, n_per_alt)
  check <- function(data, args, from, to) {
    got <- wesml_results(data, args, from, to, Q, n_per_alt)
    mapply(identical, got, ref)
  }
  expect_identical(rename_failures(data, args, renames, check), character(0))
}

test_that("WESML weights and samples are unaffected by column names", {
  # Names a used column met: data.table's V1/N, the old working column
  # `.strat`, and the locals the helpers looked up inside the data. The
  # unused column `u` is renamed to those locals as well, since the helpers
  # search every column of `data`.
  locals <- c("choice_col", "id_col", "alt_col", "keep_ids", "sampled",
              "is_chosen", "chosen")
  renames <- list(
    id = c("V1", "N", ".strat", "chosen", "id_col", "keep_ids"),
    alt = c("V1", "N", ".strat"),
    choice = c("V1", "N", "chosen", "choice_col"),
    u = locals
  )
  args <- list(id_col = "id", alt_col = "alt", choice_col = "choice")
  Q <- c(bus = 0.2, car = 0.5, rail = 0.3)
  expect_wesml_immune(wesml_fixture(), args, renames, Q = Q, n_per_alt = 5L)

  # The id-keyed weights keep the id column's attributes (a Stata label from
  # haven, say), as the row subset of the data did.
  dt <- wesml_fixture()
  setattr(dt$id, "label", "Choice situation")
  w <- wesml_weights(dt, "id", "alt", "choice", Q = Q)
  expect_identical(attr(w$id, "label"), "Choice situation")

  # Outside option: situations with no chosen alternative form stratum "o".
  dt <- wesml_fixture()
  dt[id %in% c(40L, 80L), choice := 0L]
  expect_wesml_immune(
    dt, c(args, list(include_outside_option = TRUE, outside_opt_label = "o")),
    renames, Q = c(bus = 0.2, car = 0.4, rail = 0.3, o = 0.1), n_per_alt = 2L
  )
})

test_that("hierarchical Bayes fits label alternatives under any alternative column name", {
  # run_hmnlogit()/run_hmnprobit() and predict() with newdata evaluated `am`,
  # and the stored-data path of predict() and the methods built on it
  # evaluated `d`, inside alt_mapping, where an alternative column of the
  # same name took their place.
  labels <- c("bus", "car", "rail")
  mcmc <- list(R = 200L, burn = 100L, seed = 7L)
  hmnl <- simulate_hmnl_data(N = 30, T = 2, J = 3, seed = 5)$data
  hmnp <- simulate_hmnp_data(N = 30, T = 2, J = 3, seed = 11)$data
  set(hmnl, j = "alt", value = labels[hmnl$alt])
  set(hmnp, j = "alt", value = labels[hmnp$alt])

  for (to in c("d", "am")) {
    for (model in c("hmnl", "hmnp")) {
      data <- copy(if (model == "hmnl") hmnl else hmnp)
      setnames(data, "alt", to)
      run <- if (model == "hmnl") run_hmnlogit else run_hmnprobit
      fit <- suppressWarnings(suppressMessages(
        run(data, "task", to, "choice", c("x1", "x2"), person_col = "pid",
            mcmc = mcmc)
      ))
      expect_identical(colnames(fit$draws$delta), labels)
      set.seed(1)
      in_sample <- predict(fit)
      set.seed(1)
      expect_identical(predict(fit, newdata = data), in_sample)
      expect_identical(in_sample$alternative, c(labels, "(outside)"))
      set.seed(1)
      expect_identical(colnames(elasticities(fit, "x1", n_draws = 10L)),
                       labels)
    }
  }
})
