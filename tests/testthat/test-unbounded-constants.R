# Constants without a finite maximum-likelihood estimate. With constants
# (use_asc = TRUE), run_mnlogit(), run_mxlogit() and run_nestlogit() warn,
# with class choicer_unbounded_constants and before the optimizer, when an
# alternative is never chosen, or chosen in every choice situation that
# offers it: the likelihood keeps rising as its constant falls, or grows. For
# the reference whose constant is normalized (the outside option, or the base
# alternative) the other constants move off together
# (.unbounded_constants(), .warn_unbounded_constants()).

# N choice situations over J alternatives (labelled `labels`); the chosen one
# cycles through `chosen`, so exactly those alternatives are chosen. x1 is a
# fixed and w1 a random covariate; `nest` puts the first half of the
# alternatives in nest "A", the rest in "B".
ucon_data <- function(N = 40L, J = 4L, chosen = seq_len(J),
                      labels = seq_len(J), seed = 1L) {
  set.seed(seed)
  pick <- chosen[(seq_len(N) - 1L) %% length(chosen) + 1L]
  dt <- data.table(id = rep(seq_len(N), each = J),
                   code = rep(seq_len(J), N))
  dt[, choice := as.integer(code == pick[id])]
  dt[, `:=`(alt = labels[code], x1 = rnorm(.N), w1 = rnorm(.N),
            nest = ifelse(code <= J %/% 2L, "A", "B"))]
  dt[, code := NULL]
  dt[]
}

ucon_mnl <- function(dt, ...) {
  suppressMessages(run_mnlogit(dt, "id", "alt", "choice", "x1", ...))
}
ucon_mxl <- function(dt, ...) {
  suppressMessages(run_mxlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", random_var_cols = "w1", S = 10L, ...))
}
ucon_nl <- function(dt, ...) {
  suppressMessages(run_nestlogit(
    data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
    covariate_cols = "x1", nest_col = "nest", ...))
}

UCON <- "choicer_unbounded_constants"

test_that(".unbounded_constants() reads alt_mapping's counts by row", {
  am <- data.table(alt_int = 1:4, alt = c("a", "b", "c", "d"),
                   N_OBS = c(10L, 10L, 10L, 3L),
                   N_CHOICES = c(4L, 6L, 0L, 3L))
  u <- .unbounded_constants(am)
  expect_identical(u$never, 3L)
  expect_identical(u$always, 4L)
  expect_identical(u$ref, 1L)
  expect_false(u$ref_never)
  expect_false(u$ref_always)
  expect_false(u$outside)

  # every constant finite: nothing to report
  am_ok <- data.table(alt_int = 1:3, alt = 1:3, N_OBS = c(5L, 5L, 5L),
                      N_CHOICES = c(1L, 2L, 2L))
  expect_null(.unbounded_constants(am_ok))

  # the outside option's row (alt_int 0) is the reference
  am_o <- data.table(alt_int = 0:2, alt = c(NA, "a", "b"),
                     N_OBS = c(10L, 10L, 10L), N_CHOICES = c(0L, 4L, 6L))
  u_o <- .unbounded_constants(am_o)
  expect_true(u_o$outside)
  expect_identical(u_o$ref, 1L)
  expect_true(u_o$ref_never)
  expect_length(u_o$never, 0L)
  expect_length(u_o$always, 0L)

  # without it, the base alternative (alt_int 1) is
  am_b <- data.table(alt_int = 1:3, alt = 1:3, N_OBS = c(2L, 10L, 10L),
                     N_CHOICES = c(2L, 4L, 4L))
  u_b <- .unbounded_constants(am_b)
  expect_true(u_b$ref_always)
  expect_false(u_b$ref_never)
  expect_length(u_b$always, 0L)

  # a hand-built alt_mapping without the counts is not checked
  expect_null(.unbounded_constants(data.table(alt_int = 1:2, alt = 1:2)))
})

test_that("a never-chosen alternative warns in all three models", {
  dt <- ucon_data(chosen = 1:3)
  msg <- "Never chosen in the data (their constants fall without bound): 4."
  expect_warning(fit <- ucon_mnl(dt), msg, fixed = TRUE, class = UCON)
  expect_s3_class(fit, "choicer_mnl")
  expect_warning(fit <- ucon_mxl(dt), msg, fixed = TRUE, class = UCON)
  expect_s3_class(fit, "choicer_mxl")
  expect_warning(fit <- ucon_nl(dt), msg, fixed = TRUE, class = UCON)
  expect_s3_class(fit, "choicer_nl")

  # without constants there is nothing to warn about
  expect_no_warning(ucon_mnl(dt, use_asc = FALSE), class = UCON)
  expect_no_warning(ucon_mxl(dt, use_asc = FALSE), class = UCON)
  expect_no_warning(ucon_nl(dt, use_asc = FALSE), class = UCON)

  # data where every alternative is chosen, and none wherever offered
  expect_no_warning(ucon_mnl(ucon_data()), class = UCON)
})

test_that("the warning names five alternatives and counts the rest", {
  dt <- ucon_data(N = 60L, J = 9L, chosen = 1:2)
  expect_warning(ucon_mnl(dt),
                 "fall without bound): 3, 4, 5, 6, 7 and 2 more. Drop these",
                 fixed = TRUE, class = UCON)
  # text labels are quoted
  dt <- ucon_data(chosen = 1:3, labels = c("bus", "car", "rail", "walk"))
  expect_warning(ucon_mnl(dt), 'without bound): "walk".', fixed = TRUE,
                 class = UCON)
})

test_that("an alternative chosen wherever it is offered warns", {
  dt <- ucon_data()
  # alternative 5 is offered in situations 1 to 3 only, and chosen there
  dt[id <= 3L, choice := 0L]
  extra <- data.table(id = 1:3, choice = 1L, alt = 5L, x1 = c(0.1, -0.2, 0.3),
                      w1 = c(0.4, 0.5, -0.6), nest = "B")
  dt <- rbind(dt, extra)
  data.table::setorderv(dt, c("id", "alt"))
  msg <- paste0("Chosen in every choice situation that offers them (their ",
                "constants rise without bound): 5.")
  expect_warning(ucon_mnl(dt), msg, fixed = TRUE, class = UCON)
  expect_warning(ucon_mxl(dt), msg, fixed = TRUE, class = UCON)
  expect_warning(ucon_nl(dt), msg, fixed = TRUE, class = UCON)
})

test_that("an outside option never chosen, or always, warns by name", {
  # every situation chooses an alternative in the data
  dt <- ucon_data()
  msg <- "The outside option is never chosen, so every constant rises"
  w <- expect_warning(ucon_mnl(dt, include_outside_option = TRUE),
                      class = UCON)
  expect_match(conditionMessage(w), msg, fixed = TRUE)
  # the remedy: the likelihood's supremum is the model without it
  expect_match(conditionMessage(w), "fit without it (include_outside_option = FALSE)",
               fixed = TRUE)
  # the outside option's label (NA without outside_opt_label) is not shown
  expect_false(grepl("\\bNA\\b", conditionMessage(w)))
  expect_warning(ucon_mxl(dt, include_outside_option = TRUE,
                          outside_opt_label = 0L), msg, fixed = TRUE,
                 class = UCON)

  # no situation does: the outside option is always chosen, and every
  # alternative in the data never
  dt[, choice := 0L]
  expect_warning(ucon_mnl(dt, include_outside_option = TRUE),
                 "The outside option is chosen in every choice situation",
                 fixed = TRUE, class = UCON)
})

test_that("a base alternative never chosen, or always, warns by name", {
  # without the outside option the first alternative is the base
  dt <- ucon_data(chosen = 2:4)
  expect_warning(ucon_mnl(dt),
                 paste0("The base alternative 1 (its constant normalized ",
                        "to 0) is never chosen, so the other constants rise ",
                        "together without bound: drop it, or pool it with a ",
                        "chosen alternative"),
                 fixed = TRUE, class = UCON)
  expect_warning(ucon_nl(dt), "The base alternative 1", fixed = TRUE,
                 class = UCON)

  # the base, offered in situations 1 to 5 only, is chosen in all of them
  dt <- ucon_data(chosen = 2:4)
  dt <- dt[!(alt == 1L & id > 5L)]
  dt[id <= 5L, choice := as.integer(alt == 1L)]
  expect_warning(ucon_mxl(dt),
                 paste0("is chosen in every choice situation that offers it, ",
                        "so the other constants fall together without bound: ",
                        "drop the situations that offer it, or pool it with ",
                        "an alternative that is not chosen wherever offered."),
                 fixed = TRUE, class = UCON)
})

test_that("prepared input_data warns too", {
  dt <- ucon_data(chosen = 1:3)
  inp <- prepare_mnl_data(dt, "id", "alt", "choice", "x1")
  expect_warning(suppressMessages(run_mnlogit(input_data = inp)),
                 class = UCON)
  inp <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", "w1")
  expect_warning(suppressMessages(run_mxlogit(
    input_data = inp, eta_draws = get_halton_normals(10L, inp$N, 1L))),
    class = UCON)
  inp <- prepare_nl_data(dt, "id", "alt", "choice", "x1", "nest")
  expect_warning(suppressMessages(run_nestlogit(input_data = inp)),
                 class = UCON)
})

test_that("factor labels are quoted and dates printed as dates", {
  am <- data.table(alt_int = 1:3, alt = factor(c("x", "y", "z")),
                   N_OBS = c(5L, 5L, 5L), N_CHOICES = c(2L, 3L, 0L))
  w <- expect_warning(.warn_unbounded_constants(am), class = UCON)
  expect_match(conditionMessage(w), 'without bound): "z".', fixed = TRUE)
  am <- data.table(alt_int = 1:3, alt = as.Date("2020-01-01") + 0:2,
                   N_OBS = c(5L, 5L, 5L), N_CHOICES = c(2L, 3L, 0L))
  w <- expect_warning(.warn_unbounded_constants(am), class = UCON)
  expect_match(conditionMessage(w), "without bound): 2020-01-03.", fixed = TRUE)
})

test_that("a finite bound in the direction a constant diverges holds it", {
  # alternative 4 is never chosen: its constant falls without bound unless a
  # lower bound holds it
  dt <- ucon_data(chosen = 1:3)
  expect_no_warning(ucon_mxl(dt, lower = c(ASC_4 = -20)), class = UCON)
  expect_no_warning(ucon_mxl(dt, lower = c(ASC_4 = -20), upper = c(ASC_4 = -20)),
                    class = UCON)
  expect_warning(ucon_mxl(dt, upper = c(ASC_4 = 5)), class = UCON)
  # with two never chosen, the one left unbounded is still named
  dt2 <- ucon_data(chosen = 1:2)
  expect_warning(ucon_mxl(dt2, lower = c(ASC_3 = -20)),
                 "fall without bound): 4. Drop", fixed = TRUE, class = UCON)
  # an always-chosen alternative's constant needs an upper bound
  dt3 <- ucon_data()
  dt3[id <= 3L, choice := 0L]
  dt3 <- rbind(dt3, data.table(id = 1:3, choice = 1L, alt = 5L,
                               x1 = c(0.1, -0.2, 0.3), w1 = c(0.4, 0.5, -0.6),
                               nest = "B"))
  data.table::setorderv(dt3, c("id", "alt"))
  expect_no_warning(ucon_mxl(dt3, upper = c(ASC_5 = 20)), class = UCON)
  expect_warning(ucon_mxl(dt3, lower = c(ASC_5 = -20)), class = UCON)
  # with the outside option, whose row comes first, alternative 4 never
  # chosen and the outside option chosen in every fifth situation
  dt4 <- ucon_data(chosen = 1:3)
  dt4[id %% 5L == 0L, choice := 0L]
  expect_no_warning(ucon_mxl(dt4, include_outside_option = TRUE,
                             outside_opt_label = 0L, lower = c(ASC_4 = -20)),
                    class = UCON)
  expect_warning(ucon_mxl(dt4, include_outside_option = TRUE,
                          outside_opt_label = 0L, lower = c(ASC_3 = -20)),
                 "fall without bound): 4. Drop", fixed = TRUE, class = UCON)
})

test_that("the warning comes after the variance's input checks", {
  dt <- ucon_data(chosen = 1:3)
  dt[, clu := id %% 5L]
  inp <- prepare_mnl_data(dt, "id", "alt", "choice", "x1", cluster_col = "clu")
  inp$cluster[1L] <- NA
  seen <- character()
  expect_error(withCallingHandlers(
    run_mnlogit(input_data = inp, se_method = "cluster"),
    choicer_unbounded_constants = function(w) {
      seen <<- c(seen, "warning")
      invokeRestart("muffleWarning")
    }), "`cluster` contains missing values")
  expect_length(seen, 0L)
})

test_that("the warning comes before the draws and the optimizer", {
  dt <- ucon_data(chosen = 1:3)
  seen <- character()
  record <- function(w) {
    seen <<- c(seen, "warning")
    invokeRestart("muffleWarning")
  }
  stopper <- function(theta_init, eval_f, lower, upper) {
    seen <<- c(seen, "optimizer")
    stop("optimizer reached")
  }
  expect_error(withCallingHandlers(ucon_mnl(dt, optimizer = stopper),
                                   choicer_unbounded_constants = record),
               "optimizer reached")
  expect_identical(seen, c("warning", "optimizer"))

  seen <- character()
  expect_error(withCallingHandlers(ucon_nl(dt, optimizer = stopper),
                                   choicer_unbounded_constants = record),
               "optimizer reached")
  expect_identical(seen, c("warning", "optimizer"))

  # the mixed logit's store-mode draws come after it
  seen <- character()
  local_mocked_bindings(get_halton_normals = function(...) {
    seen <<- c(seen, "draws")
    stop("draws reached")
  })
  expect_error(withCallingHandlers(ucon_mxl(dt),
                                   choicer_unbounded_constants = record),
               "draws reached")
  expect_identical(seen, c("warning", "draws"))
})

test_that(".warn_now() signals the condition and restores warn", {
  old <- options(warn = 0)
  on.exit(options(old))
  cond <- warningCondition("now", class = "choicer_test_now")
  expect_warning(.warn_now(cond), "now", class = "choicer_test_now")
  expect_equal(getOption("warn"), 0)
  # while the condition is signalled warn is 1, which prints a warning that
  # no handler muffles as it occurs
  during <- NULL
  withCallingHandlers(.warn_now(cond), choicer_test_now = function(w) {
    during <<- getOption("warn")
    invokeRestart("muffleWarning")
  })
  expect_equal(during, 1)
  expect_equal(getOption("warn"), 0)
  options(warn = 1)
  expect_warning(.warn_now(cond), class = "choicer_test_now")
  expect_equal(getOption("warn"), 1)
})
