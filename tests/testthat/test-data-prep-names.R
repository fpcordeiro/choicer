# User column names must not collide with the prepare_*_data() functions' own
# names. The preparations copy the user's index columns (ids, alternatives,
# choices, weights, clusters, decision makers) into a working table, add
# working columns to it and group it by user columns. Before those names
# carried the reserved ".choicer_" prefix, an index column called alt_int,
# idx_in_group, HB_PERSON or task_idx was overwritten, an id column called
# chosen or pos was misread, and a column called like a local variable the
# preparation looked up inside the table (levels, J, outside_opt_label,
# ids_to_drop) took that variable's place. Until the column names were read
# outside the table too, an index column called id_col, choice_col,
# person_col or task_by stopped the preparation with an unrelated error.
# The grid tests rename one input column at a time to such a name and expect
# the same preparation back.

# Names the preparations wrote or read by fixed name, names of the working
# tables' summaries and of former working columns, local variables that a
# same-named column masked, and locals the preparations use as the sole `i`
# of dt[...], which data.table looks up in the calling scope (keep, has_na,
# has_bad, bad_tasks).
clash_names <- c(
  "alt_int", "idx_in_group", "HAS_NA", "HB_PERSON", "task_idx",
  "TASK_HAS_NA", "HAS_BAD", "TASK_HAS_BAD", "N", "V1", "chosen", "pos",
  "levels", "J", "outside_opt_label", "ids_to_drop", "id_col", "choice_col",
  "person_col", "task_by", "keep", "has_na", "has_bad", "bad_tasks", "N_OBS",
  "N_CHOICES", "TAKE_RATE", "MKT_SHARE"
)
# The alternative column reappears in alt_mapping, which keeps these names.
alt_mapping_names <- c("alt_int", "N_OBS", "N_CHOICES", "TAKE_RATE", "MKT_SHARE")

# The full grids take about 8 seconds, so they run only where testthat runs
# CRAN-skipped tests (NOT_CRAN=true: local checks, and the GitHub matrix jobs,
# whose R setup sets it). On CRAN and in the gcc16 job they try one name for
# each way a clash went wrong (an overwritten working column, a result read
# under a shadowed name, and a local masked in j and i, in get() and in by),
# which keeps data.table's handling of the names checked there in about two
# seconds.
not_cran <- Sys.getenv("NOT_CRAN")  # read as testthat::skip_on_cran() does
full_grid <- if (identical(not_cran, "")) interactive() else
  isTRUE(as.logical(not_cran))
grid_names <- if (full_grid) {
  clash_names
} else {
  c("alt_int", "N", "levels", "id_col", "task_by")
}

# Seven choice situations over alternatives a, b, c, stored in reverse so the
# preparation's sort is exercised; ids are unsorted and not 1..N, so an id
# read as a count or a position shows up. Situation 70 has a missing x1 and is
# dropped. Weights and clusters are constant within person, as the panel
# mixed logit requires.
names_fixture <- function() {
  ids <- c(30, 10, 60, 20, 50, 40, 70)
  chosen <- c("a", "b", "a", "c", "b", "c", "a")
  pid <- c("p2", "p1", "p3", "p1", "p3", "p2", "p3")
  dt <- data.table(
    id = rep(ids, each = 3),
    alt = rep(c("a", "b", "c"), 7),
    pid = rep(pid, each = 3),
    x1 = round(seq(-1.2, 2.3, length.out = 21), 3),
    x2 = round(cos(1:21), 3)
  )
  dt[, choice := as.integer(alt == chosen[match(id, ids)])]
  dt[, w := unname(c(p1 = 1.5, p2 = 0.5, p3 = 2)[pid])]
  dt[, cl := unname(c(p1 = "k1", p2 = "k2", p3 = "k2")[pid])]
  dt[, nest := ifelse(alt == "a", "A", "B")]
  dt[id == 70 & alt == "b", x1 := NA]
  dt[rev(seq_len(nrow(dt)))]
}

# The same situations with the outside good chosen in situations 20 and 50.
# `physical = TRUE` adds explicit outside rows labelled "o" (chosen in 20 and
# 50), which prepare_mnl_data() removes before modelling the outside good.
names_outside_fixture <- function(physical = FALSE) {
  dt <- names_fixture()
  dt[id %in% c(20, 50), choice := 0L]
  if (physical) {
    out <- unique(dt[, .(id, pid, w, cl)])
    out[, `:=`(alt = "o", x1 = 0, x2 = 0, nest = "O",
               choice = as.integer(id %in% c(20, 50)))]
    dt <- rbindlist(list(dt, out), use.names = TRUE)
  }
  dt[]
}

# Three respondents; task ids restart within respondent, so tasks are keyed by
# (pid, task). Respondent p1's task 2 chooses the outside good, recorded on a
# physical "o" row; p3's task 1 chooses it with no row. p2's task 2 has a
# missing x1 and p3's task 2 an infinite x2: both tasks are dropped. q is an
# alternative-level covariate and cf a control-function residual.
# `outside = FALSE` gives every task one chosen inside alternative and no "o"
# row; `cross = TRUE` numbers the tasks uniquely, for person_col = NULL.
hb_names_fixture <- function(outside = TRUE, cross = FALSE) {
  pid <- rep(c("p1", "p2", "p3"), each = 6)
  task <- rep(rep(1:2, each = 3), 3)
  dt <- data.table(
    pid = pid, task = task,
    alt = rep(c("a", "b", "c"), 6),
    x1 = round(seq(-1.2, 2.29, length.out = 18), 3),
    x2 = round(sin(1:18), 3),
    cf = round(cos(1:18) / 4, 3)
  )
  dt[, q := unname(c(a = 0.1, b = 0.4, c = 0.2)[alt])]
  # per (pid, task), "" = outside
  chosen <- if (outside) c("b", "", "c", "a", "", "b")
            else c("b", "c", "c", "a", "a", "b")
  dt[, choice := as.integer(alt == chosen[.GRP]), by = .(pid, task)]
  dt[pid == "p2" & task == 2 & alt == "c", x1 := NA]
  dt[pid == "p3" & task == 2 & alt == "a", x2 := Inf]
  if (outside) {
    out <- data.table(pid = "p1", task = 2L, alt = "o", x1 = 0, x2 = 0,
                      cf = 0, q = 0, choice = 1L)
    dt <- rbindlist(list(dt, out), use.names = TRUE)
  }
  if (cross) dt[, task := 10L * match(pid, c("p1", "p2", "p3")) + task]
  dt[rev(seq_len(nrow(dt)))]
}

# A prepared object's values with the user's column names stripped, so a
# preparation of renamed columns can be compared with the reference.
prep_values <- function(p) {
  strip <- function(v) {
    if (is.data.frame(v)) return(unname(as.list(v)))
    if (is.list(v)) return(lapply(unname(v), strip))
    if (is.matrix(v)) dimnames(v) <- NULL
    unname(v)
  }
  p <- unclass(p)
  p$data_spec <- NULL
  lapply(p, strip)
}

# Run a preparation, keeping the text of its warnings and messages.
run_prep <- function(prep, data, args) {
  conditions <- character(0)
  keep <- function(cnd) conditions <<- c(conditions, conditionMessage(cnd))
  value <- withCallingHandlers(
    do.call(prep, c(list(data), args)),
    warning = function(w) { keep(w); invokeRestart("muffleWarning") },
    message = function(m) { keep(m); invokeRestart("muffleMessage") }
  )
  list(value = value, conditions = conditions)
}

# Column names and keys a prepared object shows the user.
visible_names <- function(p) {
  c(names(p$alt_mapping), key(p$alt_mapping),
    colnames(p$X), colnames(p$W), colnames(p$Z))
}

# Ways in which the preparation `got`, made with column `from` renamed to
# `to`, differs from the reference `ref`; empty when it is the reference with
# the renamed column under its new name.
rename_mismatches <- function(got, ref, from, to) {
  rename <- function(x) {
    if (is.list(x)) return(lapply(x, rename))
    if (is.character(x)) x <- replace(x, x == from, to)
    if (!is.null(names(x))) names(x) <- rename(names(x))
    x
  }
  g <- got$value
  r <- ref$value
  checks <- c(
    values = identical(prep_values(g), prep_values(r)),
    X = identical(colnames(g$X), rename(colnames(r$X))),
    W = identical(colnames(g$W), rename(colnames(r$W))),
    Z = identical(colnames(g$Z), rename(colnames(r$Z))),
    alt_mapping = identical(names(g$alt_mapping), rename(names(r$alt_mapping))),
    key = identical(key(g$alt_mapping), rename(key(r$alt_mapping))),
    data_spec = identical(g$data_spec, rename(r$data_spec)),
    conditions = identical(got$conditions, ref$conditions),
    no_leak = !any(startsWith(visible_names(g), ".choicer_"))
  )
  names(checks)[!checks]
}

# Rename each role column to each clash name in turn and expect the reference
# preparation back; failures are listed as "role -> name: checks".
expect_prep_immune <- function(prep, data, args, roles, alt_role) {
  ref <- run_prep(prep, data, args)
  expect_false(any(startsWith(visible_names(ref$value), ".choicer_")))
  failures <- character(0)
  for (role in roles) {
    to_names <- if (role == alt_role) setdiff(grid_names, alt_mapping_names)
                else grid_names
    for (to in to_names) {
      d <- copy(data)
      setnames(d, role, to)
      args_to <- lapply(args, function(a) {
        if (is.character(a)) replace(a, a == role, to) else a
      })
      bad <- tryCatch(
        rename_mismatches(run_prep(prep, d, args_to), ref, role, to),
        error = function(e) paste("error:", conditionMessage(e))
      )
      if (length(bad) > 0) {
        failures <- c(failures, sprintf("%s -> %s: %s", role, to,
                                        paste(bad, collapse = ", ")))
      }
    }
  }
  expect_identical(failures, character(0))
}

test_that("prepare_mnl_data is unaffected by the names of its input columns", {
  expect_prep_immune(
    prepare_mnl_data, names_fixture(),
    list(id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = c("x1", "x2"), weights_col = "w", cluster_col = "cl"),
    roles = c("id", "alt", "choice", "x1", "w", "cl"), alt_role = "alt"
  )
  expect_prep_immune(
    prepare_mnl_data, names_outside_fixture(physical = TRUE),
    list(id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = c("x1", "x2"), outside_opt_label = "o",
         include_outside_option = TRUE),
    roles = c("id", "alt", "choice", "x1"), alt_role = "alt"
  )
})

test_that("prepare_mxl_data is unaffected by the names of its input columns", {
  expect_prep_immune(
    prepare_mxl_data, names_fixture(),
    list(id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = "x1", random_var_cols = "x2", weights_col = "w",
         cluster_col = "cl", person_col = "pid"),
    roles = c("id", "alt", "choice", "x1", "x2", "w", "cl", "pid"),
    alt_role = "alt"
  )
  expect_prep_immune(
    prepare_mxl_data, names_outside_fixture(),
    list(id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = "x1", random_var_cols = "x2",
         outside_opt_label = "o", include_outside_option = TRUE),
    roles = c("id", "alt", "choice", "x2"), alt_role = "alt"
  )
})

test_that("prepare_nl_data is unaffected by the names of its input columns", {
  expect_prep_immune(
    prepare_nl_data, names_fixture(),
    list(id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = c("x1", "x2"), nest_col = "nest",
         weights_col = "w", cluster_col = "cl"),
    roles = c("id", "alt", "choice", "x1", "nest", "w", "cl"),
    alt_role = "alt"
  )
})

test_that("prepare_mnp_data is unaffected by the names of its input columns", {
  expect_prep_immune(
    prepare_mnp_data, names_fixture(),
    list(id_col = "id", alt_col = "alt", choice_col = "choice",
         covariate_cols = c("x1", "x2")),
    roles = c("id", "alt", "choice", "x1"), alt_role = "alt"
  )
})

test_that("prepare_hmnl_data and prepare_hmnp_data are unaffected by the names of their input columns", {
  args <- list(id_col = "task", alt_col = "alt", choice_col = "choice",
               covariate_cols = c("x1", "x2"), person_col = "pid",
               alt_covariate_cols = "q", cf_residual_col = "cf",
               outside_opt_label = "o")
  roles <- c("task", "alt", "choice", "x1", "pid", "q", "cf")
  expect_prep_immune(prepare_hmnl_data, hb_names_fixture(),
                     c(args, list(rc_dist = c(0L, 1L))), roles, "alt")
  expect_prep_immune(prepare_hmnp_data, hb_names_fixture(), args, roles, "alt")
  # Each task its own respondent (person_col = NULL), and no outside good.
  cross_args <- args[setdiff(names(args), "person_col")]
  expect_prep_immune(prepare_hmnl_data, hb_names_fixture(cross = TRUE),
                     cross_args, setdiff(roles, "pid"), "alt")
  inside_args <- c(args[setdiff(names(args), "outside_opt_label")],
                   list(include_outside_option = FALSE))
  expect_prep_immune(prepare_hmnp_data, hb_names_fixture(outside = FALSE),
                     inside_args, roles, "alt")
})

test_that("an alternative column named like an alt_mapping column is an error", {
  fx <- names_fixture()
  hb <- hb_names_fixture()
  for (nm in alt_mapping_names) {
    msg <- paste0("alternative column cannot be named '", nm, "'")
    d <- copy(fx)
    setnames(d, "alt", nm)
    expect_error(prepare_mnl_data(d, "id", nm, "choice", "x1"), msg,
                 fixed = TRUE)
    expect_error(prepare_mxl_data(d, "id", nm, "choice", "x1", "x2"), msg,
                 fixed = TRUE)
    expect_error(prepare_nl_data(d, "id", nm, "choice", "x1", "nest"), msg,
                 fixed = TRUE)
    expect_error(prepare_mnp_data(d, "id", nm, "choice", "x1"), msg,
                 fixed = TRUE)
    h <- copy(hb)
    setnames(h, "alt", nm)
    expect_error(prepare_hmnl_data(h, "task", nm, "choice", "x1",
                                   person_col = "pid"), msg, fixed = TRUE)
    expect_error(prepare_hmnp_data(h, "task", nm, "choice", "x1",
                                   person_col = "pid"), msg, fixed = TRUE)
  }
})

test_that("column names with the reserved .choicer_ prefix are an error", {
  msg <- "Column names beginning with '.choicer_' are reserved"
  d <- names_fixture()
  setnames(d, c("x1", "pid"), c(".choicer_alt_int", ".choicer_person"))
  expect_error(prepare_mnl_data(d, "id", "alt", "choice", ".choicer_alt_int"),
               msg, fixed = TRUE)
  expect_error(prepare_mxl_data(d, "id", "alt", "choice", "x2", "w",
                                person_col = ".choicer_person"),
               msg, fixed = TRUE)
  expect_error(prepare_nl_data(d, "id", "alt", "choice", ".choicer_alt_int",
                               "nest"), msg, fixed = TRUE)
  n <- copy(d)
  setnames(n, "nest", ".choicer_nest")
  expect_error(prepare_nl_data(n, "id", "alt", "choice", "x2", ".choicer_nest"),
               msg, fixed = TRUE)
  expect_error(prepare_mnp_data(d, "id", "alt", "choice", ".choicer_alt_int"),
               msg, fixed = TRUE)
  h <- hb_names_fixture()
  setnames(h, "q", ".choicer_q")
  expect_error(prepare_hmnl_data(h, "task", "alt", "choice", "x1",
                                 alt_covariate_cols = ".choicer_q"),
               msg, fixed = TRUE)
  expect_error(prepare_hmnp_data(h, "task", "alt", "choice", "x1",
                                 alt_covariate_cols = ".choicer_q"),
               msg, fixed = TRUE)
})

test_that("unused columns with the reserved prefix are ignored", {
  fx <- names_fixture()
  extra <- copy(fx)[, .choicer_alt_int := -1]
  expect_identical(
    suppressWarnings(prepare_mnl_data(extra, "id", "alt", "choice", "x1")),
    suppressWarnings(prepare_mnl_data(fx, "id", "alt", "choice", "x1"))
  )
  hb <- hb_names_fixture()
  extra_hb <- copy(hb)[, .choicer_person := "z"]
  expect_identical(
    suppressWarnings(prepare_hmnp_data(extra_hb, "task", "alt", "choice", "x1",
                                       person_col = "pid",
                                       outside_opt_label = "o")),
    suppressWarnings(prepare_hmnp_data(hb, "task", "alt", "choice", "x1",
                                       person_col = "pid",
                                       outside_opt_label = "o"))
  )
})

test_that("an alternative column named J labels the ASCs of a fit", {
  # run_*logit() read the ASC labels as alt_mapping[2:J], where a column
  # named J masked the local J: with labels 1 to 3 the third alternative's
  # ASC was labelled as the base alternative's, and other labels errored.
  # (The estimates are covered by the preparation grids above; separate fits
  # of these flat likelihoods differ in the last digits with the kernels'
  # OpenMP reduction order.)
  expect_same_names <- function(fit_fun, data, ...) {
    d_J <- copy(data)
    setnames(d_J, "alt", "J")
    fit <- suppressMessages(fit_fun(data, "id", "alt", "choice", ...))
    fit_J <- suppressMessages(fit_fun(d_J, "id", "J", "choice", ...))
    expect_identical(names(coef(fit_J)), names(coef(fit)))
  }
  expect_same_names(run_mnlogit, create_small_mnl_data(), c("x1", "x2"))
  expect_same_names(run_mxlogit, create_small_mxl_data(), "x1", "w1", S = 20L)
  expect_same_names(run_nestlogit, create_small_nl_data(), c("x1", "x2"),
                    nest_col = "nest")
})

test_that("task columns named with join operators prepare as under any other name", {
  # The hierarchical preparations drop the situations of a flagged row by an
  # anti-join, and data.table parses the strings of `on =` as conditions.
  hb <- hb_names_fixture()
  args <- list(id_col = "task", alt_col = "alt", choice_col = "choice",
               covariate_cols = c("x1", "x2"), person_col = "pid",
               outside_opt_label = "o")
  ref <- run_prep(prepare_hmnl_data, hb, args)
  for (nm in c("t>0", "t==1", "t<=2", " task")) {
    d <- copy(hb)
    setnames(d, "task", nm)
    got <- run_prep(prepare_hmnl_data, d, modifyList(args, list(id_col = nm)))
    expect_identical(rename_mismatches(got, ref, "task", nm), character(0),
                     label = paste0("mismatches for a task column named '",
                                    nm, "'"))
  }
})

test_that("the preparations work without spare column slots", {
  # options(datatable.alloccol = 0) leaves new tables without room to add a
  # column in place; the working columns must still be added.
  fx <- names_fixture()
  hb <- hb_names_fixture()
  prep_all <- function() {
    suppressWarnings(suppressMessages(list(
      mnl = prepare_mnl_data(fx, "id", "alt", "choice", c("x1", "x2")),
      mxl = prepare_mxl_data(fx, "id", "alt", "choice", "x1", "x2",
                             person_col = "pid"),
      nl = prepare_nl_data(fx, "id", "alt", "choice", c("x1", "x2"), "nest"),
      mnp = prepare_mnp_data(fx, "id", "alt", "choice", c("x1", "x2")),
      hmnl = prepare_hmnl_data(hb, "task", "alt", "choice", c("x1", "x2"),
                               person_col = "pid", outside_opt_label = "o")
    )))
  }
  ref <- prep_all()
  op <- options(datatable.alloccol = 0L)
  got <- tryCatch(prep_all(), finally = options(op))
  expect_identical(lapply(got, prep_values), lapply(ref, prep_values))
})

