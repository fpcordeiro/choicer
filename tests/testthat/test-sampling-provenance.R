# Characterize the shared post-fit WESML contract through the public wrappers.
# A fixed custom optimizer keeps these tests cheap and platform-independent;
# data preparation, covariance calculations, and fitted objects remain real.
provenance_fixture <- function(weights = c(0.5, 1, 2, 0.5, 1, 2)) {
  data.frame(
    id = rep(seq_len(6), each = 3),
    alt = rep(seq_len(3), 6),
    choice = as.integer(rep(seq_len(3), 6) == rep(c(1, 2, 3, 2, 1, 3), each = 3)),
    x = sin(seq_len(18)), w = cos(seq_len(18)),
    nest = rep(c("A", "A", "B"), 6),
    cluster = rep(rep(seq_len(3), each = 2), each = 3),
    weight = rep(weights, each = 3)
  )
}

provenance_metadata <- function() {
  structure(
    list(scheme = "wesml", weight_name = "weight", se_method = "old",
         weights_applied = NA, custom = list(label = "retained", values = 1:3)),
    class = "choicer_sampling", source = "characterization"
  )
}

fit_provenance_case <- function(model, pathway, weights, metadata, se_method) {
  data <- provenance_fixture(weights)
  attr(data, "choice_sampling") <- metadata
  columns <- list(id_col = "id", alt_col = "alt", choice_col = "choice",
                  covariate_cols = "x", weights_col = "weight", cluster_col = "cluster")
  if (model == "mxl") columns$random_var_cols <- "w"
  if (model == "nl") columns$nest_col <- "nest"
  prepare <- switch(model, mnl = prepare_mnl_data, mxl = prepare_mxl_data,
                    nl = prepare_nl_data)
  run <- switch(model, mnl = run_mnlogit, mxl = run_mxlogit, nl = run_nestlogit)
  args <- if (pathway == "prepared") {
    list(input_data = do.call(prepare, c(list(data = data), columns)))
  } else {
    c(list(data = data), columns)
  }
  if (model == "mxl" && pathway == "prepared") {
    args$eta_draws <- get_halton_normals(5L, 6L, 1L)
  }
  original <- serialize(args, NULL)
  events <- character()
  warnings <- list()
  args$optimizer <- function(theta_init, eval_f, ...) {
    events <<- c(events, "optimizer")
    list(par = theta_init, value = eval_f(theta_init)$objective, convergence = 0L)
  }
  args$se_method <- se_method
  if (model == "mxl") args$S <- 5L
  result <- suppressMessages(withCallingHandlers(
    tryCatch(do.call(run, args), error = function(e) {
      events <<- c(events, "error")
      e
    }),
    warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      events <<- c(events, "warning")
      invokeRestart("muffleWarning")
    }
  ))
  args$optimizer <- args$se_method <- args$S <- NULL
  expect_identical(serialize(args, NULL), original)
  list(result = result, warnings = warnings, events = events)
}

provenance_weight_warning <- paste0(
  "Non-uniform weights detected. If these are sampling/WESML ",
  "weights, use se_method = 'sandwich' for valid inference."
)
provenance_bhhh_warning <- paste0(
  "Non-uniform weights detected with se_method = 'bhhh': BHHH/OPG ",
  "standard errors use the w^1 meat (sum w_i s_i s_i')^{-1}, which is ",
  "NOT a valid choice-based-sampling (WESML) correction; the correct ",
  "sandwich meat is w^2. Use se_method = 'sandwich' for valid WESML inference."
)
provenance_uniform_warning <- paste0(
  "WESML provenance is present but the applied weights are uniform; the fit ",
  "is effectively unweighted and is NOT a WESML-corrected estimator."
)

test_that("weighted provenance and warnings agree across models and input pathways", {
  for (model in c("mnl", "mxl", "nl")) {
    methods <- c("hessian", "bhhh", "sandwich", "cluster")
    if (model == "nl") methods <- c(methods, "numeric")
    for (pathway in c("convenience", "prepared")) {
      for (se_method in methods) {
        for (metadata in list(NULL, provenance_metadata())) {
          case <- fit_provenance_case(model, pathway, c(0.5, 1, 2, 0.5, 1, 2),
                                      metadata, se_method)
          expect_s3_class(case$result, paste0("choicer_", model))
          expected <- metadata
          if (is.null(expected)) expected <- list(scheme = "user")
          expected$se_method <- se_method
          expected$weights_applied <- TRUE
          expect_identical(case$result$choice_sampling, expected)
          message <- switch(se_method, hessian = provenance_weight_warning,
                            numeric = provenance_weight_warning,
                            bhhh = provenance_bhhh_warning, character())
          expect_identical(vapply(case$warnings, conditionMessage, character(1)), message)
          expect_true(all(vapply(case$warnings, function(w) is.null(conditionCall(w)),
                                 logical(1))))
          expect_identical(case$events, c("optimizer", rep("warning", length(message))))
        }
      }
    }
  }
})

test_that("uniform weights preserve the distinct convenience and prepared guards", {
  for (model in c("mnl", "mxl", "nl")) {
    for (pathway in c("convenience", "prepared")) {
      plain <- fit_provenance_case(model, pathway, rep(2, 6), NULL, "bhhh")
      expect_s3_class(plain$result, paste0("choicer_", model))
      expect_null(plain$result$choice_sampling)
      expect_length(plain$warnings, 0L)
      expect_identical(plain$events, "optimizer")

      metadata <- provenance_metadata()
      case <- fit_provenance_case(model, pathway, rep(2, 6), metadata, "sandwich")
      if (pathway == "convenience") {
        expect_s3_class(case$result, paste0("choicer_", model))
        metadata$se_method <- "sandwich"
        metadata$weights_applied <- FALSE
        expect_identical(case$result$choice_sampling, metadata)
        expect_length(case$warnings, 1L)
        expect_identical(conditionMessage(case$warnings[[1]]), provenance_uniform_warning)
        expect_null(conditionCall(case$warnings[[1]]))
        expect_identical(case$events, c("optimizer", "warning"))
      } else {
        expect_s3_class(case$result, "error")
        prepare_name <- paste0("prepare_", model, "_data")
        expected <- paste0(
          "`input_data` is flagged as a WESML choice-based sample (it carries ",
          "`choice_sampling` provenance), but the resolved weights are uniform. ",
          "Fitting would produce an invalid unweighted estimator mislabeled as ",
          "WESML. To proceed, either bake the non-uniform WESML weights into ",
          "`input_data` via ", prepare_name, "(weights = ) / ", prepare_name,
          "(weights_col = ), or, if you deliberately want an unweighted fit, ",
          "strip the provenance with `attr(input_data, \"choice_sampling\") <- NULL`."
        )
        expect_identical(conditionMessage(case$result), expected)
        expect_null(conditionCall(case$result))
        expect_length(case$warnings, 0L)
        expect_identical(case$events, c("optimizer", "error"))
      }
    }
  }
})

test_that("almost equal weights remain nonuniform for provenance", {
  weights <- c(rep(1, 5), 1 + .Machine$double.eps)
  for (model in c("mnl", "mxl", "nl")) {
    case <- fit_provenance_case(model, "prepared", weights, provenance_metadata(),
                                "hessian")
    expect_s3_class(case$result, paste0("choicer_", model))
    expect_true(case$result$choice_sampling$weights_applied)
    expect_length(case$warnings, 1L)
    expect_identical(conditionMessage(case$warnings[[1]]), provenance_weight_warning)
    expect_null(conditionCall(case$warnings[[1]]))
  }
})
