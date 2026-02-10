
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
#' @param covariate_cols Vector of column names for covariates.
#' @param nest_col Name of the column mapping each alternative to its nest
#'   (convenience workflow).
#' @param input_data List containing prepared input data for estimation
#'   (advanced workflow). Mutually exclusive with \code{data}.
#' @param use_asc Logical indicating whether to include alternative specific
#'   constants (ASCs).
#' @param theta_init Optional initial parameter vector. If \code{NULL}, a
#'   default vector is used.
#' @param param_names Optional vector of parameter names. If \code{NULL},
#'   default names are generated.
#' @param optimizer Optimizer to use: \code{"nloptr"} (default), \code{"optim"},
#'   or a custom function. See \code{\link{run_mnlogit}} for details.
#' @param control List of optimizer-specific control parameters.
#' @param weights Optional weight vector (convenience workflow). If \code{NULL},
#'   equal weights are used.
#' @param outside_opt_label Label for the outside option (convenience workflow).
#' @param include_outside_option Logical whether to include an outside option
#'   (convenience workflow).
#' @param keep_data Logical. If \code{TRUE} (default), stores prepared data in
#'   the returned object for post-estimation functions.
#' @param nloptr_opts Deprecated. Use \code{optimizer} and \code{control}
#'   instead.
#' @returns A \code{choicer_nl} object (inherits from \code{choicer_fit}).
#'   Standard S3 methods available: \code{summary()}, \code{coef()},
#'   \code{vcov()}, \code{logLik()}, \code{AIC()}, \code{BIC()},
#'   \code{nobs()}.
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
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    keep_data = TRUE,
    nloptr_opts = NULL
) {
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

  if (has_data && has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced), not both.")
  }
  if (!has_data && !has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced).")
  }

  if (has_data) {
    # Convenience workflow: validate required column-name arguments
    if (is.null(id_col) || is.null(alt_col) || is.null(choice_col) ||
        is.null(covariate_cols) || is.null(nest_col)) {
      stop("Convenience workflow requires: id_col, alt_col, choice_col, ",
           "covariate_cols, and nest_col.")
    }
    input_data <- prepare_nl_data(
      data = data,
      id_col = id_col,
      alt_col = alt_col,
      choice_col = choice_col,
      covariate_cols = covariate_cols,
      nest_col = nest_col,
      weights = weights,
      outside_opt_label = outside_opt_label,
      include_outside_option = include_outside_option
    )
  }

  # Parameter dimensions
  J <- nrow(input_data$alt_mapping)
  K_x <- ncol(input_data$X)
  K_l <- sum(table(input_data$nest_idx) > 1)
  n_asc <- if (use_asc) J - 1 else 0
  n_params <- K_x + K_l + n_asc

  # Initial parameter vector
  if (is.null(theta_init)) {
    theta_init <- c(rep(0, K_x), rep(0.5, K_l), rep(0, n_asc))
  }

  # Lower bounds: lambda must be > 0
  theta_lb <- c(rep(-Inf, K_x), rep(1e-16, K_l), rep(-Inf, n_asc))

  # Build eval_f closure
  eval_f <- function(theta) {
    nl_loglik_gradient_parallel(
      theta = theta,
      X = input_data$X,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      nest_idx = input_data$nest_idx,
      M = input_data$M,
      weights = input_data$weights,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option
    )
  }

  # Run optimizer
  elapsed <- system.time({
    opt <- run_optimizer(
      optimizer = optimizer,
      theta_init = theta_init,
      eval_f = eval_f,
      lower = theta_lb,
      control = control
    )
  })

  message("Optimization run time ", convertTime(elapsed))

  # Parameter names and index map
  theta_hat <- opt$par

  if (is.null(param_names)) {
    beta_names <- colnames(input_data$X)
    if (is.null(beta_names)) beta_names <- paste0("X_", seq_len(K_x))
    lambda_names <- paste0("Lambda_", seq_len(K_l))
    alt_col <- names(input_data$alt_mapping)[2]
    asc_names <- if (use_asc) {
      paste0("ASC_", input_data$alt_mapping[2:J][[alt_col]])
    } else {
      character(0)
    }
    param_names <- c(beta_names, lambda_names, asc_names)
  }
  names(theta_hat) <- param_names

  # Parameter index map
  param_map <- list(beta = seq_len(K_x))
  param_map$lambda <- K_x + seq_len(K_l)
  if (n_asc > 0) param_map$asc <- K_x + K_l + seq_len(n_asc)

  # Extract lambda values
  lambda <- theta_hat[param_map$lambda]

  # Compute vcov eagerly
  hess <- nl_loglik_numeric_hessian(
    theta = theta_hat,
    X = input_data$X,
    alt_idx = input_data$alt_idx,
    choice_idx = input_data$choice_idx,
    nest_idx = input_data$nest_idx,
    M = input_data$M,
    weights = input_data$weights,
    use_asc = use_asc,
    include_outside_option = input_data$include_outside_option
  )
  vcov_result <- invert_hessian(hess)
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
        weights = input_data$weights
      )
    },
    lambda = lambda
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
#' @param weights Optional vector of weights for each choice situation. If \code{NULL}, equal weights are used.
#' @param outside_opt_label Label for the outside option (if any). If \code{NULL}, no outside option is assumed.
#' @param include_outside_option Logical indicating whether to include an outside option in the model.
#' @returns A \code{choicer_data_nl} object (list) containing:
#'   \itemize{
#'     \item All fields from \code{\link{prepare_mnl_data}} (\code{X}, \code{alt_idx},
#'       \code{choice_idx}, \code{M}, \code{N}, \code{weights}, \code{include_outside_option},
#'       \code{alt_mapping}, \code{dropped_cols}).
#'     \item \code{nest_idx}: Integer vector of length J mapping each alternative
#'       (in \code{alt_mapping} row order) to its nest.
#'     \item \code{data_spec}: List with column name metadata including \code{nest_col}.
#'   }
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
    include_outside_option = FALSE
) {
  dt <- as.data.table(data)[]

  # Validate nest_col exists
  if (!nest_col %in% names(dt)) {
    stop("Missing column: ", nest_col)
  }

  # Extract unique alt -> nest mapping
  nest_map <- unique(dt[, c(alt_col, nest_col), with = FALSE])

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
    outside_opt_label = outside_opt_label,
    include_outside_option = include_outside_option
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

  structure(result, class = "choicer_data_nl")
}
