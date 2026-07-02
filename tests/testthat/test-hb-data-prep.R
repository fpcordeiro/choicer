# Tests for the shared hierarchical-Bayes data prep (R/hb_data.R) through its
# exported wrappers:
# - prepare_hmnl_data()
# - prepare_hmnp_data()
# Covers the two-level (person, task) indexing, the implicit outside option,
# the alternative-level design Z, cf_residual handling, rc_dist alignment,
# and the conditional task-constant covariate rule.

# Hand-built 3-person panel with known indexing:
#   P1: task 1 (a, b, c; chose b)  task 2 (a, c; chose the outside row "out")
#   P2: task 3 (b, c; chose c)
#   P3: task 4 (a, b, c; no choice = outside)  task 5 (a, b; chose a)
# Physical outside rows exist for tasks 1 and 2 only; rows are stored in
# reversed order so the prep's (person, task, alt) sort is exercised.
make_hb_panel_data <- function() {
  inside <- data.table(
    pid  = c("P1", "P1", "P1", "P1", "P1", "P2", "P2",
             "P3", "P3", "P3", "P3", "P3"),
    task = c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 4L, 5L, 5L),
    alt  = c("a", "b", "c", "a", "c", "b", "c", "a", "b", "c", "a", "b"),
    x1   = as.numeric(1:12),
    x2   = c(0.3, -0.1, 0.4, 0.2, -0.5, 0.7, 0.1, -0.2, 0.6, -0.3, 0.05, 0.25),
    choice = c(0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L)
  )
  qual_map <- c(a = 0.1, b = 0.2, c = 0.3)
  inside[, qual := qual_map[alt]]
  outside <- data.table(
    pid = c("P1", "P1"), task = c(1L, 2L), alt = c("out", "out"),
    x1 = 0, x2 = 0, choice = c(0L, 1L), qual = 0
  )
  dt <- rbindlist(list(inside, outside), use.names = TRUE)
  dt[rev(seq_len(nrow(dt)))]
}

test_that("prepare_hmnl_data round-trips the hand-built 3-person panel", {
  dt <- make_hb_panel_data()
  d <- prepare_hmnl_data(
    dt, "task", "alt", "choice", c("x1", "x2"),
    person_col = "pid", alt_covariate_cols = "qual",
    outside_opt_label = "out", include_outside_option = TRUE
  )

  expect_s3_class(d, "choicer_data_hmnl")
  expect_true(inherits(d, "list"))

  # Two-level indexing
  expect_equal(d$N_persons, 3L)
  expect_equal(d$person_ids, c("P1", "P2", "P3"))
  expect_equal(d$Ti, c(2L, 1L, 2L))
  expect_equal(d$n_tasks, 5L)
  expect_equal(sum(d$Ti), d$n_tasks)

  # Task-level round trip: M, choice_pos (0 = outside chosen)
  expect_equal(d$M, c(3L, 2L, 2L, 3L, 2L))
  expect_equal(d$choice_pos, c(2L, 0L, 2L, 0L, 1L))
  expect_equal(cumsum(d$M), c(3L, 5L, 7L, 10L, 12L))

  # Row-level round trip: sort restored the (person, task, alt) order
  expect_equal(d$J, 3L)
  expect_equal(d$alt_of_row,
               c(1L, 2L, 3L, 1L, 3L, 2L, 3L, 1L, 2L, 3L, 1L, 2L))
  expect_equal(d$alt_idx, d$alt_of_row)
  expect_equal(nrow(d$X), 12L)          # outside rows removed
  expect_equal(unname(d$X[, "x1"]), as.numeric(1:12))
  expect_equal(d$K_struct, 2L)
  expect_equal(colnames(d$X), c("x1", "x2"))

  # Z build: intercept first, one deduplicated row per alternative
  expect_equal(d$P, 2L)
  expect_equal(colnames(d$Z), c("(Intercept)", "qual"))
  expect_equal(unname(d$Z[, 1]), rep(1, 3))
  expect_equal(unname(d$Z[, 2]), c(0.1, 0.2, 0.3))

  # param_map robust naming
  expect_equal(unname(d$param_map$beta), 1:2)
  expect_equal(names(d$param_map$theta), c("(Intercept)", "qual"))

  # alt_mapping: outside option is alt_int = 0 with tasks 2 and 4 outside
  expect_equal(d$alt_mapping$alt_int[1], 0L)
  expect_equal(d$alt_mapping$alt[1], "out")
  expect_equal(d$alt_mapping$N_OBS[1], 5L)
  expect_equal(d$alt_mapping$N_CHOICES[1], 2L)

  # rc_dist defaults to all-normal aligned with X
  expect_equal(unname(d$rc_dist), c(0L, 0L))
  expect_equal(names(d$rc_dist), c("x1", "x2"))

  # data_spec provenance
  expect_equal(d$data_spec$person_col, "pid")
  expect_equal(d$data_spec$outside_opt_label, "out")
  expect_null(d$data_spec$cf_residual_col)
  expect_true(d$include_outside_option)
})

test_that("the all-zeros choice encodes an outside choice without physical rows", {
  dt <- make_hb_panel_data()[alt != "out"]   # no physical outside rows at all
  d <- prepare_hmnp_data(dt, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", include_outside_option = TRUE)
  expect_s3_class(d, "choicer_data_hmnp")
  expect_equal(d$choice_pos, c(2L, 0L, 2L, 0L, 1L))
  expect_equal(nrow(d$X), 12L)
  # outside label unknown -> NA in the mapping, counts still correct
  expect_true(is.na(d$alt_mapping$alt[1]))
  expect_equal(d$alt_mapping$N_CHOICES[1], 2L)
})

test_that("cf_residual_col is appended to X as an ordinary covariate", {
  dt <- make_hb_panel_data()
  dt[, cfres := x1^2 / 10]
  d <- prepare_hmnp_data(dt, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", outside_opt_label = "out",
                         cf_residual_col = "cfres")
  expect_equal(colnames(d$X), c("x1", "x2", "cfres"))
  expect_equal(d$K_struct, 3L)
  expect_equal(unname(d$X[, "cfres"]), (1:12)^2 / 10)
  expect_equal(d$data_spec$cf_residual_col, "cfres")

  # cf residual must not double as a structural covariate
  expect_error(
    prepare_hmnp_data(dt, "task", "alt", "choice", c("x1", "cfres"),
                      person_col = "pid", outside_opt_label = "out",
                      cf_residual_col = "cfres"),
    "must not also appear"
  )
})

test_that("person_col = NULL makes each task its own respondent", {
  dt <- make_hb_panel_data()
  d <- prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "x2"),
                         outside_opt_label = "out")
  expect_equal(d$N_persons, 5L)
  expect_true(all(d$Ti == 1L))
  expect_equal(d$person_ids, 1:5)      # person ids inherit the task ids
  expect_equal(d$M, c(3L, 2L, 2L, 3L, 2L))
  expect_equal(d$choice_pos, c(2L, 0L, 2L, 0L, 1L))
  expect_null(d$data_spec$person_col)
})

test_that("rc_dist validates against the input and aligns through drops", {
  dt <- make_hb_panel_data()
  dt[, cfres := x1^2 / 10]
  # x3 is an exact (tiny-norm) linear combination -> QR drops it
  dt[, x3 := (x1 + x2) * 1e-9]

  expect_error(
    prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "x2"),
                      person_col = "pid", outside_opt_label = "out",
                      rc_dist = c(0L, 1L, 0L)),
    "must have length 2"
  )
  expect_error(
    prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "x2"),
                      person_col = "pid", outside_opt_label = "out",
                      rc_dist = c(0L, 2L)),
    "0 \\(normal\\) or 1 \\(log-normal\\)"
  )

  expect_message(
    d <- prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "x2", "x3"),
                           person_col = "pid", outside_opt_label = "out",
                           cf_residual_col = "cfres",
                           rc_dist = c(0L, 1L, 0L)),
    "collinearity"
  )
  expect_equal(colnames(d$X), c("x1", "x2", "cfres"))
  expect_equal(d$dropped_cols, "x3")
  # x2 keeps its log-normal flag; the cf residual coordinate is always normal
  expect_equal(unname(d$rc_dist), c(0L, 1L, 0L))
  expect_equal(names(d$rc_dist), c("x1", "x2", "cfres"))
})

test_that("task-constant covariates are kept with an outside option, dropped without", {
  dt <- make_hb_panel_data()
  dt[, xc := as.numeric(task)]         # constant within every task

  # WITH the outside option: identified vs the outside good -> kept + message
  expect_message(
    d <- prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "x2", "xc"),
                           person_col = "pid", outside_opt_label = "out",
                           include_outside_option = TRUE),
    "constant within every choice situation kept"
  )
  expect_true("xc" %in% colnames(d$X))
  expect_null(d$dropped_cols)

  # WITHOUT: unidentified -> dropped with a warning
  dt2 <- make_hb_panel_data()[alt != "out"]
  dt2[task == 2 & alt == "a", choice := 1L]   # give every task an inside choice
  dt2[task == 4 & alt == "c", choice := 1L]
  dt2[, xc := as.numeric(task)]
  expect_warning(
    d2 <- prepare_hmnl_data(dt2, "task", "alt", "choice", c("x1", "x2", "xc"),
                            person_col = "pid",
                            include_outside_option = FALSE),
    "not identified without an outside option"
  )
  expect_false("xc" %in% colnames(d2$X))
  expect_true("xc" %in% d2$dropped_cols)
  expect_equal(unname(d2$rc_dist), c(0L, 0L))
})

test_that("one-choice-per-task enforcement follows the outside-option convention", {
  # Two chosen alternatives in one task always errors
  dt <- make_hb_panel_data()
  dt[task == 1 & alt == "a", choice := 1L]    # task 1 now has b AND a chosen
  expect_error(
    prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "x2"),
                      person_col = "pid", outside_opt_label = "out"),
    "at most one chosen"
  )

  # Without an outside option a zero-choice task errors
  dt2 <- make_hb_panel_data()[alt != "out"]   # tasks 2 and 4 have no choice
  expect_error(
    prepare_hmnp_data(dt2, "task", "alt", "choice", c("x1", "x2"),
                      person_col = "pid", include_outside_option = FALSE),
    "exactly one chosen"
  )

  # Non-0/1 choice column errors (on an inside row: outside rows are removed
  # before the 0/1 check)
  dt3 <- make_hb_panel_data()
  dt3[task == 1 & alt == "b", choice := 2L]
  expect_error(
    prepare_hmnl_data(dt3, "task", "alt", "choice", c("x1", "x2"),
                      person_col = "pid", outside_opt_label = "out"),
    "only 0 and 1"
  )
})

test_that("Z build validates and prunes alternative-level covariates", {
  dt <- make_hb_panel_data()

  # Intercept-only Z when alt_covariate_cols = NULL
  d0 <- prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "x2"),
                          person_col = "pid", outside_opt_label = "out")
  expect_equal(d0$P, 1L)
  expect_equal(colnames(d0$Z), "(Intercept)")
  expect_equal(unname(d0$Z[, 1]), rep(1, 3))

  # Not constant within an alternative -> error
  dtv <- make_hb_panel_data()
  dtv[, qual2 := qual]
  dtv[task == 3 & alt == "b", qual2 := 9]
  expect_error(
    prepare_hmnl_data(dtv, "task", "alt", "choice", c("x1", "x2"),
                      person_col = "pid", alt_covariate_cols = "qual2",
                      outside_opt_label = "out"),
    "constant within each alternative"
  )

  # Constant ACROSS alternatives -> dropped from Z with a message
  dtc <- make_hb_panel_data()
  dtc[, zflat := 7]
  expect_message(
    dz <- prepare_hmnl_data(dtc, "task", "alt", "choice", c("x1", "x2"),
                            person_col = "pid",
                            alt_covariate_cols = c("qual", "zflat"),
                            outside_opt_label = "out"),
    "constant across alternatives"
  )
  expect_equal(colnames(dz$Z), c("(Intercept)", "qual"))
  expect_equal(dz$dropped_z_cols, "zflat")
})

test_that("prep validates choice-set sizes and missing values", {
  # M >= 2 required without an outside option
  dt <- make_hb_panel_data()[alt != "out"]
  dt[task == 2 & alt == "a", choice := 1L]
  dt[task == 4 & alt == "c", choice := 1L]
  dt_small <- dt[!(task == 3 & alt == "b")]   # task 3 down to one alternative
  expect_error(
    prepare_hmnl_data(dt_small, "task", "alt", "choice", c("x1", "x2"),
                      person_col = "pid", include_outside_option = FALSE),
    "at least 2 alternatives"
  )
  # ...but a single inside alternative is fine WITH the outside option
  dt_small2 <- make_hb_panel_data()[!(task == 3 & alt == "b")]
  d <- prepare_hmnl_data(dt_small2, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", outside_opt_label = "out",
                         include_outside_option = TRUE)
  expect_equal(d$M, c(3L, 2L, 1L, 3L, 2L))

  # NA task dropping removes the whole task (and its person if emptied)
  dtn <- make_hb_panel_data()
  dtn[task == 3 & alt == "b", x1 := NA]
  expect_warning(
    dn <- prepare_hmnl_data(dtn, "task", "alt", "choice", c("x1", "x2"),
                            person_col = "pid", outside_opt_label = "out"),
    "1 choice situations containing missing values"
  )
  expect_equal(dn$n_tasks, 4L)
  expect_equal(dn$N_persons, 2L)       # P2 had only task 3
  expect_equal(dn$Ti, c(2L, 2L))

  # Missing columns error
  expect_error(
    prepare_hmnl_data(make_hb_panel_data(), "task", "alt", "choice",
                      c("x1", "nope"), person_col = "pid"),
    "Missing columns"
  )
})

test_that("price-like covariates without a cf residual trigger a reminder message", {
  dt <- make_hb_panel_data()
  dt[, price := x1^2 / 10]
  expect_message(
    prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "price"),
                      person_col = "pid", outside_opt_label = "out"),
    "cf_residual_col"
  )
  # supplying the residual silences the reminder
  dt[, cfres := x2^2]
  expect_no_message(
    prepare_hmnl_data(dt, "task", "alt", "choice", c("x1", "price"),
                      person_col = "pid", outside_opt_label = "out",
                      cf_residual_col = "cfres")
  )
})

test_that("prepare_hmnl_data round-trips simulate_hmnl_data output", {
  sim <- simulate_hmnl_data(N = 12, T = 3, J = 3, seed = 5,
                            include_outside = TRUE)
  d <- prepare_hmnl_data(sim$data, "task", "alt", "choice", c("x1", "x2"),
                         person_col = "pid", alt_covariate_cols = "z1",
                         outside_opt_label = 0L)
  expect_equal(d$n_tasks, 36L)
  expect_equal(d$N_persons, 12L)
  expect_true(all(d$Ti == 3L))
  expect_equal(d$J, 3L)
  expect_equal(d$M, rep(3L, 36L))
  # the prep's Z reproduces the DGP's mean-function design
  expect_equal(unname(d$Z), unname(sim$true_params$Z))
  # outside-chosen tasks are the all-zeros ones
  n_outside <- sim$data[alt == 0L & choice == 1L, .N]
  expect_equal(sum(d$choice_pos == 0L), n_outside)
})
