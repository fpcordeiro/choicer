
#' Runs nested logit estimation
#'
#' Estimates a nested logit model via maximum likelihood.
#'
#' Two workflows are supported:
#' \describe{
#'   \item{Convenience}{Supply \code{data} and column names (including
#'     \code{nest_col}). Data preparation (\code{\link{prepare_nl_data}}) is
#'     handled automatically.}
#'   \item{Advanced}{Call \code{\link{prepare_nl_data}} (or build the input
#'     list manually) and pass it via \code{input_data}.}
#' }
#'
#' @param data Data frame containing choice data (convenience workflow).
#'   Mutually exclusive with \code{input_data}.
#' @param id_col Name of the column identifying choice situations.
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating chosen alternative (1/0).
#' @param covariate_cols Vector of column names for covariates. None may be
#'   named like a generated parameter: \code{Lambda_<k>}, or
#'   \code{ASC_<label>} with \code{use_asc = TRUE}.
#' @param nest_col Name of the column mapping each alternative to its nest
#'   (convenience workflow).
#' @param input_data List containing prepared input data for estimation
#'   (advanced workflow). Mutually exclusive with \code{data}.
#' @param use_asc Logical indicating whether to include alternative specific
#'   constants (ASCs).
#' @param theta_init Optional initial parameter vector in natural units,
#'   ordered as the coefficients, the nest parameters (\code{Lambda_<k>}, one
#'   per nest of two or more alternatives) and the constants. If \code{NULL}:
#'   0 for the coefficients and constants, 0.5 for each nest parameter.
#' @param param_names Optional vector of parameter names, which must be
#'   unique. If \code{NULL}, default names are generated.
#' @param optimizer Optimizer to use: \code{"nloptr"} (default), \code{"optim"},
#'   or a custom function. See \code{\link{run_mnlogit}} for details.
#' @param control List of optimizer-specific control parameters. With
#'   \code{scale_vars} other than \code{"none"}, settings in parameter units
#'   apply to the scaled coordinates.
#' @param weights Optional weight vector (convenience workflow). If \code{NULL},
#'   equal weights are used. All weights must be finite and strictly positive.
#' @param weights_col Optional name of a column in \code{data} holding per-row
#'   weights (convenience workflow only). The column must be constant within each
#'   \code{id_col} (one weight per choice situation) and is collapsed accordingly.
#'   Mutually exclusive with \code{weights}. All weights must be finite and strictly
#'   positive. Used for choice-based / WESML
#'   weighting; pair with \code{se_method = "sandwich"} for valid inference.
#' @param outside_opt_label Label for the outside option (convenience workflow).
#' @param include_outside_option Logical whether to include an outside option
#'   (convenience workflow).
#' @param keep_data Logical. If \code{TRUE} (default), stores prepared data in
#'   the returned object for post-estimation functions.
#' @param se_method Method for computing standard errors: \code{"hessian"}
#'   (default, analytical Hessian via \code{nl_loglik_hessian_parallel}),
#'   \code{"numeric"} (finite-difference oracle via
#'   \code{nl_loglik_numeric_hessian}), \code{"bhhh"} (outer product of
#'   gradients via \code{nl_bhhh_parallel}), \code{"sandwich"} (robust
#'   Huber--White / WESML variance \eqn{A^{-1} B A^{-1}}), or \code{"cluster"}
#'   (cluster-robust sandwich; requires \code{cluster_col} or a prepared
#'   \code{input_data} with a \code{cluster} field). Use \code{"sandwich"}
#'   under choice-based / WESML weighting. Any of these can also be recomputed
#'   post hoc via \code{vcov(fit, type = )}.
#' @param cluster_col Optional name of a column in \code{data} holding cluster
#'   labels for cluster-robust standard errors (e.g. a person id when the same
#'   decision maker contributes several choice situations). Must be constant
#'   within each \code{id_col}. Supplying \code{cluster_col} without an explicit
#'   \code{se_method} selects \code{se_method = "cluster"}.
#' @param nloptr_opts Deprecated. Use \code{optimizer} and \code{control}
#'   instead.
#' @param scale_vars How the optimizer's coordinates are scaled. The choice
#'   does not change the estimator; it changes how fast, and whether, the
#'   optimizer reaches the maximum. One of \code{"none"} (default), \code{"sd"}
#'   (sample standard deviation), \code{"mad"} (\code{stats::mad}),
#'   \code{"iqr"} (\code{stats::IQR(x) / 1.349}), or \code{"bhhh"}. With
#'   \code{"sd"}, \code{"mad"} or \code{"iqr"}, the optimizer works on the
#'   coefficients \code{X} would have with every column divided by the chosen
#'   scale (each coefficient times its column's scale), to improve the
#'   conditioning of its problem; the nest parameters and constants keep their
#'   scale, and the data are not divided. \code{"bhhh"} scales every parameter,
#'   nest parameters included, by the inverse square root of the diagonal of
#'   the BHHH (outer product of gradients) matrix at the start values, rounded
#'   to a power of two, from one extra gradient pass. Its scales are fixed at
#'   the start values, so when the nest parameters end far from their start
#'   (0.5) they can mislead the optimizer, which then converges more slowly
#'   than without them; hence \code{"none"} is the default here. Where the
#'   diagonal is unusable (no finite, positive entry), the parameters are left
#'   unscaled, with a message. With any choice but \code{"none"}, the
#'   optimizer, a custom one included, works in the scaled coordinates
#'   (\code{theta_init}, in natural units, and the nest parameters' lower bound
#'   of 1e-16 are mapped), so settings in parameter units in \code{control}
#'   (nloptr's \code{xtol_abs}; optim's \code{parscale}, and \code{pgtol},
#'   which then applies to the scaled gradient) apply to the scaled
#'   coordinates; \code{"none"} keeps the parameters' units. Coefficients are
#'   reported in natural units and the standard errors computed in them (for
#'   any choice but \code{"none"}, with the information matrix equilibrated
#'   before it is inverted), so, where the likelihood is well identified and
#'   every choice reaches its maximum, reported quantities do not depend on
#'   this choice beyond the optimizer's tolerance.
#' @returns A \code{choicer_nl} object (inherits from \code{choicer_fit}).
#'   Standard S3 methods available: \code{summary()}, \code{coef()},
#'   \code{vcov()}, \code{logLik()}, \code{AIC()}, \code{BIC()},
#'   \code{nobs()}.
#' @inheritSection prepare_mnl_data Column names
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 100; J <- 4
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, nest := ifelse(alt <= 2, "A", "B")]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#'
#' fit <- run_nestlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = c("x1", "x2"), nest_col = "nest"
#' )
#' summary(fit)
#' }
#' @importFrom nloptr nloptr
#' @export
run_nestlogit <- function(
    data = NULL,
    id_col = NULL,
    alt_col = NULL,
    choice_col = NULL,
    covariate_cols = NULL,
    nest_col = NULL,
    input_data = NULL,
    use_asc = TRUE,
    theta_init = NULL,
    param_names = NULL,
    optimizer = NULL,
    control = list(),
    weights = NULL,
    weights_col = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    keep_data = TRUE,
    se_method = c("hessian", "numeric", "bhhh", "sandwich", "cluster"),
    cluster_col = NULL,
    nloptr_opts = NULL,
    scale_vars = c("none", "sd", "mad", "iqr", "bhhh")
) {
  se_method_default <- missing(se_method)
  se_method <- match.arg(se_method)
  scale_vars <- match.arg(scale_vars)
  if (!is.null(cluster_col) && se_method_default) se_method <- "cluster"
  cl <- match.call()

  # Backward compatibility: nloptr_opts -> optimizer + control
  if (!is.null(nloptr_opts)) {
    message("'nloptr_opts' is deprecated. Use 'optimizer' and 'control' instead.")
    optimizer <- optimizer %||% "nloptr"
    control <- nloptr_opts
  }

  # --- Resolve input pathway --------------------------------------------------
  has_data <- !is.null(data)
  has_input <- !is.null(input_data)
  cs_meta <- if (has_data) attr(data, "choice_sampling") else attr(input_data, "choice_sampling")

  if (has_data && has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced), not both.")
  }
  if (!has_data && !has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced).")
  }
  if (has_input && !is.null(weights_col)) {
    stop("`weights_col` is only supported in the convenience (data) workflow. ",
         "Bake weights into `input_data` via prepare_nl_data(weights_col = ) ",
         "or supply `weights` to prepare_nl_data().")
  }
  if (has_input && !is.null(cluster_col)) {
    stop("`cluster_col` is only supported in the convenience (data) workflow. ",
         "Bake cluster labels into `input_data` via ",
         "prepare_nl_data(cluster_col = ).")
  }

  if (has_data) {
    # Convenience workflow: validate required column-name arguments
    if (is.null(id_col) || is.null(alt_col) || is.null(choice_col) ||
        is.null(covariate_cols) || is.null(nest_col)) {
      stop("Convenience workflow requires: id_col, alt_col, choice_col, ",
           "covariate_cols, and nest_col.")
    }
    # WESML provenance present but no weights supplied: auto-adopt the recorded
    # weight column, or error -- never silently fit unweighted under a WESML label.
    if (!is.null(cs_meta) && is.null(weights) && is.null(weights_col)) {
      wn <- cs_meta$weight_name
      if (!is.null(wn) && wn %in% names(data)) {
        weights_col <- wn
        message("Detected WESML choice-based-sampling provenance; applying attached ",
                "weights from column '", wn, "'.")
      } else {
        stop("Data carries WESML choice-based-sampling provenance but no weights were ",
             "supplied, and the recorded weight column (",
             if (is.null(wn)) "unknown" else paste0("'", wn, "'"),
             ") is not present in `data`. Pass `weights_col=` or `weights=` explicitly.",
             call. = FALSE)
      }
    }
    input_data <- prepare_nl_data(
      data = data,
      id_col = id_col,
      alt_col = alt_col,
      choice_col = choice_col,
      covariate_cols = covariate_cols,
      nest_col = nest_col,
      weights = weights,
      weights_col = weights_col,
      outside_opt_label = outside_opt_label,
      include_outside_option = include_outside_option,
      cluster_col = cluster_col
    )
  }

  if (se_method == "cluster" && is.null(input_data$cluster)) {
    stop("se_method = \"cluster\" needs cluster labels: pass `cluster_col=` ",
         "(convenience workflow) or prepare `input_data` with ",
         "prepare_nl_data(cluster_col = ).", call. = FALSE)
  }

  # Parameter dimensions
  J <- nrow(input_data$alt_mapping)
  K_x <- ncol(input_data$X)
  K_l <- sum(table(input_data$nest_idx) > 1)
  n_asc <- if (use_asc) J - 1 else 0
  n_params <- K_x + K_l + n_asc

  # Parameter names, checked for repeats before the optimizer runs
  names_supplied <- !is.null(param_names)
  if (!names_supplied) {
    beta_names <- colnames(input_data$X)
    if (is.null(beta_names)) beta_names <- paste0("X_", seq_len(K_x))
    lambda_names <- paste0("Lambda_", seq_len(K_l))
    alt_col <- names(input_data$alt_mapping)[2]
    asc_names <- if (use_asc) {
      paste0("ASC_", input_data$alt_mapping[[alt_col]][2:J])
    } else {
      character(0)
    }
    param_names <- c(beta_names, lambda_names, asc_names)
  }
  .check_param_names(list(param_names), c("Lambda_<k>", "ASC_<label>"),
                     supplied = names_supplied)

  # Initial parameter vector
  if (is.null(theta_init)) {
    theta_init <- c(rep(0, K_x), rep(0.5, K_l), rep(0, n_asc))
  }

  # Lower bounds: lambda must be > 0
  theta_lb <- c(rep(-Inf, K_x), rep(1e-16, K_l), rep(-Inf, n_asc))

  # Parameter index map
  param_map <- list(beta = seq_len(K_x))
  param_map$lambda <- K_x + seq_len(K_l)
  if (n_asc > 0) param_map$asc <- K_x + K_l + seq_len(n_asc)

  # The model's objective and gradient at theta; `...` takes opg_diag
  model_f <- function(theta, ...) {
    nl_loglik_gradient_parallel(
      theta = theta,
      X = input_data$X,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      nest_idx = input_data$nest_idx,
      M = input_data$M,
      weights = input_data$weights,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option,
      ...
    )
  }

  # --- The optimizer's coordinates (scale_vars) -------------------------------
  #   theta_natural = map$scale * theta_t + map$shift
  # "sd", "mad" and "iqr" take the map from the columns' scales (beta only);
  # "bhhh" from the BHHH diagonal at the start values, from one gradient pass.
  sX <- rep(1, K_x); names(sX) <- colnames(input_data$X)
  col_scaled <- scale_vars %in% c("sd", "mad", "iqr")
  if (col_scaled && K_x > 0) {
    sX <- .column_scales(input_data$X, scale_vars)
    .assert_scales_ok(sX, scale_vars, "fixed-coefficient")
  }
  if (scale_vars != "none") .check_theta_init(theta_init, n_params)
  t_bpass <- system.time(
    opg0 <- if (scale_vars == "bhhh") model_f(theta_init, opg_diag = TRUE)$opg_diag,
    gcFirst = FALSE
  )
  map <- .coordinate_map(scale_vars, param_map, n_params, sX = sX,
                         opg_diag = opg0)

  # Run the optimizer in the map's coordinates
  elapsed <- system.time({
    opt <- .optimize_in_coordinates(
      map = if (scale_vars != "none") map, optimizer = optimizer,
      theta_init = theta_init, eval_f = function(theta) model_f(theta),
      lower = theta_lb, control = control
    )
  })

  elapsed <- elapsed + t_bpass  # "bhhh"'s gradient pass counts
  message("Optimization run time ", convertTime(elapsed))

  # Estimates in natural units
  theta_hat <- opt$par
  names(theta_hat) <- param_names

  # Extract lambda values
  lambda <- theta_hat[param_map$lambda]

  # Choice-based-sampling provenance and a guardrail for weighted inference.
  choice_sampling <- .resolve_choice_sampling(
    weights = input_data$weights, se_method = se_method,
    cs_meta = cs_meta, has_input = has_input, prepare_fn = "prepare_nl_data"
  )

  # Compute vcov eagerly using the selected SE method. For "sandwich"
  # (robust / WESML) errors, form V = A^{-1} B A^{-1} with bread A = weighted
  # negated Hessian and meat B = weight-squared OPG (pass weights^2 to the
  # weight-free BHHH routine). For "cluster", the meat is the outer product of
  # within-cluster sums of weighted scores. At the natural estimates; a scaled
  # fit's matrices are inverted equilibrated, as post hoc.
  eq <- scale_vars != "none"
  if (se_method %in% c("sandwich", "cluster")) {
    A_bread <- nl_loglik_hessian_parallel(
      theta = theta_hat, X = input_data$X, alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx, nest_idx = input_data$nest_idx,
      M = input_data$M, weights = input_data$weights, use_asc = use_asc,
      include_outside_option = input_data$include_outside_option
    )
    B_meat <- if (se_method == "sandwich") {
      nl_bhhh_parallel(
        theta = theta_hat, X = input_data$X, alt_idx = input_data$alt_idx,
        choice_idx = input_data$choice_idx, nest_idx = input_data$nest_idx,
        M = input_data$M, weights = input_data$weights^2, use_asc = use_asc,
        include_outside_option = input_data$include_outside_option
      )
    } else {
      S_scores <- nl_scores_parallel(
        theta = theta_hat, X = input_data$X, alt_idx = input_data$alt_idx,
        choice_idx = input_data$choice_idx, nest_idx = input_data$nest_idx,
        M = input_data$M, use_asc = use_asc,
        include_outside_option = input_data$include_outside_option
      )
      .score_meat(S_scores, input_data$weights, "cluster", input_data$cluster)
    }
    vcov_result <- .sandwich_combine(A_bread, B_meat, equilibrate = eq)
  } else {
    hess <- switch(
      se_method,
      numeric = nl_loglik_numeric_hessian(
        theta = theta_hat, X = input_data$X, alt_idx = input_data$alt_idx,
        choice_idx = input_data$choice_idx, nest_idx = input_data$nest_idx,
        M = input_data$M, weights = input_data$weights, use_asc = use_asc,
        include_outside_option = input_data$include_outside_option
      ),
      bhhh = nl_bhhh_parallel(
        theta = theta_hat, X = input_data$X, alt_idx = input_data$alt_idx,
        choice_idx = input_data$choice_idx, nest_idx = input_data$nest_idx,
        M = input_data$M, weights = input_data$weights, use_asc = use_asc,
        include_outside_option = input_data$include_outside_option
      ),
      nl_loglik_hessian_parallel(
        theta = theta_hat, X = input_data$X, alt_idx = input_data$alt_idx,
        choice_idx = input_data$choice_idx, nest_idx = input_data$nest_idx,
        M = input_data$M, weights = input_data$weights, use_asc = use_asc,
        include_outside_option = input_data$include_outside_option
      )
    )
    vcov_result <- invert_hessian(hess, equilibrate = eq)
  }
  if (!is.null(vcov_result$vcov)) {
    rownames(vcov_result$vcov) <- param_names
    colnames(vcov_result$vcov) <- param_names
    names(vcov_result$se) <- param_names
  }

  # Build S3 object
  new_choicer_nl(
    call = cl,
    coefficients = theta_hat,
    loglik = -opt$value,
    nobs = input_data$N,
    n_params = n_params,
    convergence = opt$convergence,
    message = opt$message,
    data_spec = input_data$data_spec,
    alt_mapping = input_data$alt_mapping,
    param_map = param_map,
    use_asc = use_asc,
    include_outside_option = input_data$include_outside_option,
    optimizer = list(
      name = if (is.function(optimizer)) "custom" else (optimizer %||% "nloptr"),
      control = control,
      elapsed_time = elapsed[["elapsed"]],
      iterations = opt$iterations
    ),
    vcov = vcov_result$vcov,
    se = vcov_result$se,
    data = if (keep_data) {
      list(
        X = input_data$X,
        alt_idx = input_data$alt_idx,
        choice_idx = input_data$choice_idx,
        nest_idx = input_data$nest_idx,
        M = input_data$M,
        weights = input_data$weights,
        cluster = input_data$cluster,
        situation_ids = input_data$situation_ids
      )
    },
    lambda = lambda,
    nest_idx = input_data$nest_idx,
    se_method = se_method,
    choice_sampling = choice_sampling,
    scale_vars = scale_vars,
    sX = sX,
    param_scale = stats::setNames(map$scale, param_names),
    param_shift = stats::setNames(map$shift, param_names)
  )
}


#' Prepare inputs for nested logit estimation
#'
#' Validates inputs, builds design matrices, and constructs nest structure
#' for nested logit estimation. Calls \code{\link{prepare_mnl_data}} internally
#' for base data preparation, then adds nest-specific fields.
#'
#' @param data Data frame containing choice data.
#' @param id_col Name of the column identifying choice situations (individuals).
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating chosen alternative (1 = chosen, 0 = not chosen).
#' @param covariate_cols Vector of names of columns to be used as covariates.
#' @param nest_col Name of the column mapping each alternative to its nest.
#'   Every alternative must belong to exactly one nest.
#' @param weights Optional vector of weights for each choice situation. If \code{NULL}, equal weights are used. All weights must be finite and strictly positive.
#' @param weights_col Optional name of a column in \code{data} holding per-row
#'   weights. The column must be constant within each \code{id_col} (one weight
#'   per choice situation) and is collapsed accordingly. Mutually exclusive with
#'   \code{weights}. All weights must be finite and strictly positive.
#' @param outside_opt_label Label for the outside option (if any). If \code{NULL}, no outside option is assumed.
#' @param include_outside_option Logical indicating whether to include an outside option in the model.
#' @param cluster_col Optional name of a column in \code{data} holding cluster
#'   labels for cluster-robust standard errors. Must be constant within each
#'   \code{id_col}; collapsed to one label per choice situation and returned as
#'   \code{cluster}.
#' @returns A \code{choicer_data_nl} object (list) containing:
#'   \itemize{
#'     \item All fields from \code{\link{prepare_mnl_data}} (\code{X}, \code{alt_idx},
#'       \code{choice_idx}, \code{M}, \code{N}, \code{weights}, \code{cluster},
#'       \code{situation_ids}, \code{include_outside_option}, \code{alt_mapping},
#'       \code{dropped_cols}).
#'     \item \code{nest_idx}: Integer vector of length J mapping each alternative
#'       (in \code{alt_mapping} row order) to its nest.
#'     \item \code{data_spec}: List with column name metadata including \code{nest_col}.
#'   }
#' @inheritSection prepare_mnl_data Column names
#' @examples
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 4
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, nest := ifelse(alt <= 2, "A", "B")]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' input <- prepare_nl_data(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
#' input$nest_idx
#' input$alt_mapping
#' @export
prepare_nl_data <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    nest_col,
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    weights_col = NULL,
    cluster_col = NULL
) {
  # Names first: they need no data, and reading the data can stop for a
  # missing bit64. (prepare_mnl_data() below checks its own columns again.)
  .check_col_names(c(id_col, alt_col, choice_col, covariate_cols, weights_col,
                     cluster_col, nest_col), alt_col)
  # Only the two columns of the alternative-to-nest map; prepare_mnl_data()
  # below copies the columns it needs itself.
  dt <- .copy_cols(data, c(alt_col, nest_col))
  .load_bit64_for(dt)

  # Validate nest_col exists
  if (!nest_col %in% names(dt)) {
    stop("Missing column: ", nest_col)
  }

  # Extract unique alt -> nest mapping
  nest_map <- unique(dt[, c(alt_col, nest_col), with = FALSE])
  rm(dt)

  # Validate: each alternative belongs to exactly one nest
  if (anyDuplicated(nest_map[[alt_col]])) {
    bad_alts <- nest_map[[alt_col]][duplicated(nest_map[[alt_col]])]
    stop("Alternatives belong to multiple nests: ",
         paste(unique(bad_alts), collapse = ", "))
  }

  # Validate: at least 2 nests
  unique_nests <- unique(nest_map[[nest_col]])
  if (length(unique_nests) < 2) {
    stop("At least 2 nests are required; found ", length(unique_nests), ".")
  }

  # Validate: no missing nest assignments
  if (any(is.na(nest_map[[nest_col]]))) {
    stop("Missing nest assignments (NA) in column '", nest_col, "'.")
  }

  # Call prepare_mnl_data() for base data preparation
  result <- prepare_mnl_data(
    data = data,
    id_col = id_col,
    alt_col = alt_col,
    choice_col = choice_col,
    covariate_cols = covariate_cols,
    weights = weights,
    weights_col = weights_col,
    outside_opt_label = outside_opt_label,
    include_outside_option = include_outside_option,
    cluster_col = cluster_col
  )

  # Build nest_idx aligned with alt_mapping row order (inside alternatives only;

  # the outside option is handled implicitly in C++ when include_outside_option=TRUE)
  if (include_outside_option) {
    alt_labels <- result$alt_mapping[alt_int > 0][[alt_col]]
  } else {
    alt_labels <- result$alt_mapping[[alt_col]]
  }
  nest_labels <- nest_map[[nest_col]][match(alt_labels, nest_map[[alt_col]])]

  # Check all alternatives have a nest assignment
  if (any(is.na(nest_labels))) {
    missing_alts <- alt_labels[is.na(nest_labels)]
    stop("No nest assignment found for alternatives: ",
         paste(missing_alts, collapse = ", "))
  }

  # Convert nest labels to 1-based integers (sorted order)
  nest_levels <- sort(unique(nest_labels))
  nest_idx <- as.integer(factor(nest_labels, levels = nest_levels))

  # Validate: every nest has at least 1 alternative
  # (guaranteed by construction, but verify)
  if (length(unique(nest_idx)) != length(nest_levels)) {
    stop("Internal error: nest count mismatch after integer conversion.")
  }

  # Carry choice-based-sampling provenance from the MNL base preparation.
  cs_provenance <- attr(result, "choice_sampling")

  # Add NL-specific fields
  result$nest_idx <- nest_idx
  result$data_spec <- list(
    id_col = id_col,
    alt_col = alt_col,
    choice_col = choice_col,
    covariate_cols = covariate_cols,
    nest_col = nest_col,
    outside_opt_label = outside_opt_label
  )

  out <- structure(result, class = "choicer_data_nl")
  if (!is.null(cs_provenance)) {
    attr(out, "choice_sampling") <- cs_provenance
  }
  out
}
