# The fit-time variance of run_mnlogit(), run_mxlogit() and run_nestlogit().
# Its inputs (the weights against se_method and the sampling provenance, the
# cluster labels) are checked before any evaluation of the likelihood.

# Ninety choice situations of three alternatives, three per decision maker
vf_data <- function(seed = 7L) {
  set.seed(seed)
  dt <- data.table::data.table(id = rep(1:90, each = 3L), alt = rep(1:3, 90L))
  dt[, `:=`(x1 = stats::rnorm(.N), x2 = stats::rnorm(.N))]
  dt[, nest := c("A", "A", "B")[alt]]
  dt[, person := (id - 1L) %/% 3L + 1L]
  dt[, V := 0.8 * x1 - 0.5 * x2]
  dt[, choice := as.integer(alt == sample(alt, 1L, prob = exp(V) / sum(exp(V)))),
     by = id]
  dt[, V := NULL]
  dt[]
}

vf_fit <- function(model, dt, ...) {
  suppressMessages(switch(model,
    mnl = run_mnlogit(data = dt, id_col = "id", alt_col = "alt",
                      choice_col = "choice", covariate_cols = c("x1", "x2"), ...),
    mxl = run_mxlogit(data = dt, id_col = "id", alt_col = "alt",
                      choice_col = "choice", covariate_cols = "x1",
                      random_var_cols = "x2", rc_mean = TRUE, S = 20L, ...),
    nl = run_nestlogit(data = dt, id_col = "id", alt_col = "alt",
                       choice_col = "choice", covariate_cols = c("x1", "x2"),
                       nest_col = "nest", ...)
  ))
}

# Counts the evaluations of the three models' likelihood (the kernel behind
# the optimizer and the "bhhh" pass) in `counter$n` until the calling test
# ends
vf_count_evals <- function(env = parent.frame()) {
  counter <- new.env()
  counter$n <- 0L
  count <- function(name) {
    real <- get(name, envir = asNamespace("choicer"))
    function(...) {
      counter$n <- counter$n + 1L
      real(...)
    }
  }
  local_mocked_bindings(
    mnl_loglik_gradient_parallel = count("mnl_loglik_gradient_parallel"),
    mxl_loglik_gradient_parallel = count("mxl_loglik_gradient_parallel"),
    nl_loglik_gradient_parallel = count("nl_loglik_gradient_parallel"),
    .env = env
  )
  counter
}

# Prepared inputs of each model ("mxl_panel": decision makers by `person`)
vf_prepare <- function(model, dt, ...) {
  suppressWarnings(switch(model,
    mnl = prepare_mnl_data(dt, "id", "alt", "choice", c("x1", "x2"), ...),
    mxl = prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2", ...),
    mxl_panel = prepare_mxl_data(dt, "id", "alt", "choice", "x1", "x2",
                                 person_col = "person", ...),
    nl = prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest", ...)
  ))
}

vf_fit_input <- function(model, input, ...) {
  run <- switch(sub("_panel", "", model), mnl = run_mnlogit,
                mxl = run_mxlogit, nl = run_nestlogit)
  args <- list(input_data = input, ...)
  if (startsWith(model, "mxl")) {
    args$eta_draws <- get_halton_normals(5L, length(choicer:::.unit_first(input)), 1L)
  }
  suppressMessages(do.call(run, args))
}

test_that("fit-time cluster labels are checked before any likelihood evaluation", {
  evals <- vf_count_evals()
  dt <- vf_data()
  dt[, team := (person - 1L) %/% 2L + 1L]  # clusters that nest decision makers
  for (model in c("mnl", "mxl", "nl", "mxl_panel")) {
    input <- vf_prepare(model, dt, cluster_col = "team")
    # well-formed labels fit, with their variance
    evals$n <- 0L
    fit <- vf_fit_input(model, input, se_method = "cluster")
    expect_s3_class(fit, "choicer_fit")
    expect_null(fit$variance_error)
    expect_gt(evals$n, 0L)

    bad <- input
    bad$cluster[2L] <- NA
    evals$n <- 0L
    expect_error(vf_fit_input(model, bad, se_method = "cluster"),
                 if (model == "mxl_panel") {
                   "must be constant within each decision maker"
                 } else {
                   "`cluster` contains missing values."
                 }, fixed = TRUE, label = model)
    expect_identical(evals$n, 0L, label = paste(model, "NA"))

    bad <- input
    bad$cluster <- bad$cluster[-1L]
    evals$n <- 0L
    expect_error(vf_fit_input(model, bad, se_method = "cluster"),
                 if (model == "mxl_panel") {
                   "`cluster_col`: got 89 values for 90 choice situations."
                 } else {
                   "`cluster` has length 89 but the model has 90 choice situations."
                 }, fixed = TRUE, label = model)
    expect_identical(evals$n, 0L, label = paste(model, "length"))
  }
  # labels with an NA for a whole decision maker reach the NA check
  bad <- vf_prepare("mxl_panel", dt, cluster_col = "team")
  bad$cluster[1:3] <- NA
  evals$n <- 0L
  expect_error(vf_fit_input("mxl_panel", bad, se_method = "cluster"),
               "`cluster` contains missing values.", fixed = TRUE)
  expect_identical(evals$n, 0L)
  # a panel whose Ti does not split the situations is the kernels' to
  # report, at their first call, whatever the labels
  bad <- vf_prepare("mxl_panel", dt, cluster_col = "team")
  bad$Ti <- bad$Ti[-1L]
  bad$cluster[2L] <- NA
  evals$n <- 0L
  expect_error(vf_fit_input("mxl_panel", bad, se_method = "cluster"),
               "sum(Ti) (87) does not match the number of choice situations (90).",
               fixed = TRUE)
  expect_identical(evals$n, 1L)  # the "bhhh" pass
})

test_that("the weights' guard runs before any likelihood evaluation", {
  evals <- vf_count_evals()
  dt <- vf_data()
  # a prepared input that carries choice-based-sampling provenance with
  # uniform weights stops
  flagged <- data.table::copy(dt)
  flagged[, w := 2]
  data.table::setattr(flagged, "choice_sampling", structure(
    list(scheme = "wesml", weight_name = "w"), class = "choicer_sampling"))
  for (model in c("mnl", "mxl", "nl")) {
    input <- vf_prepare(model, flagged, weights_col = "w")
    evals$n <- 0L
    expect_error(vf_fit_input(model, input, se_method = "sandwich"),
                 "`input_data` is flagged as a WESML choice-based sample",
                 fixed = TRUE, label = model)
    expect_identical(evals$n, 0L, label = model)
    # before the cluster labels, as when both ran after the optimizer
    input <- vf_prepare(model, flagged, weights_col = "w", cluster_col = "person")
    input$cluster[2L] <- NA
    evals$n <- 0L
    expect_error(vf_fit_input(model, input, se_method = "cluster"),
                 "`input_data` is flagged as a WESML choice-based sample",
                 fixed = TRUE, label = model)
    expect_identical(evals$n, 0L, label = model)
  }
  # non-uniform weights without the sandwich or clustered variance warn
  weighted <- data.table::copy(dt)
  weighted[, w := 1 + id %% 2]
  for (model in c("mnl", "mxl", "nl")) {
    for (se_method in c("hessian", "bhhh")) {
      evals$n <- 0L
      at <- integer(0)
      fit <- withCallingHandlers(
        vf_fit(model, weighted, weights_col = "w", se_method = se_method),
        warning = function(w) {
          if (startsWith(conditionMessage(w), "Non-uniform weights detected")) {
            at <<- c(at, evals$n)
            invokeRestart("muffleWarning")
          }
        })
      expect_s3_class(fit, "choicer_fit")
      expect_null(fit$variance_error)
      expect_identical(at, 0L, label = paste(model, se_method))
      expect_gt(evals$n, 0L)
    }
  }
})
