
# S3 class constructors for choicer model objects
#
# Class hierarchy:
#   choicer_mnl -> choicer_fit
#   choicer_mxl -> choicer_fit
#   choicer_nl  -> choicer_fit

# --- Constructors -----------------------------------------------------------

#' Construct a choicer_mnl object
#' @param call The matched call from run_mnlogit()
#' @param coefficients Named numeric vector of point estimates
#' @param loglik Scalar log-likelihood at optimum (positive)
#' @param nobs Integer number of choice situations
#' @param n_params Integer number of estimated parameters
#' @param convergence Integer convergence code from optimizer
#' @param message Character string from optimizer
#' @param data_spec List with column name metadata (id_col, alt_col, etc.)
#' @param alt_mapping data.table mapping alternatives to summary statistics
#' @param param_map Named list of integer index vectors (beta, asc)
#' @param use_asc Logical whether ASCs were used
#' @param include_outside_option Logical whether outside option was included
#' @param optimizer List with optimizer metadata (name, control, elapsed_time)
#' @param vcov Named variance-covariance matrix (or NULL for lazy computation)
#' @param se Named numeric vector of standard errors (or NULL for lazy computation)
#' @param data List of prepared inputs (X, alt_idx, choice_idx, M, weights) or NULL
#' @returns A choicer_mnl object (S3 class)
#' @noRd
new_choicer_mnl <- function(call, coefficients, loglik,
                            nobs, n_params, convergence, message,
                            data_spec, alt_mapping, param_map,
                            use_asc, include_outside_option,
                            optimizer,
                            vcov = NULL, se = NULL, data = NULL) {
  structure(
    list(
      call = call,
      model = "mnl",
      coefficients = coefficients,
      vcov = vcov,
      se = se,
      loglik = loglik,
      nobs = nobs,
      n_params = n_params,
      convergence = convergence,
      message = message,
      data_spec = data_spec,
      alt_mapping = alt_mapping,
      param_map = param_map,
      use_asc = use_asc,
      include_outside_option = include_outside_option,
      optimizer = optimizer,
      data = data
    ),
    class = c("choicer_mnl", "choicer_fit")
  )
}

#' Construct a choicer_mxl object
#' @inheritParams new_choicer_mnl
#' @param draws_info List with Halton draw metadata (S, N, K_w) for regeneration
#' @param rc_dist Integer vector of random coefficient distributions
#' @param rc_correlation Logical whether correlated random coefficients
#' @param rc_mean Logical whether means were estimated
#' @param sigma Reconstructed covariance matrix of random coefficients (or NULL)
#' @returns A choicer_mxl object (S3 class)
#' @noRd
new_choicer_mxl <- function(call, coefficients, loglik,
                            nobs, n_params, convergence, message,
                            data_spec, alt_mapping, param_map,
                            use_asc, include_outside_option,
                            optimizer,
                            vcov = NULL, se = NULL, data = NULL,
                            draws_info = NULL,
                            rc_dist = NULL, rc_correlation = FALSE,
                            rc_mean = FALSE, sigma = NULL) {
  structure(
    list(
      call = call,
      model = "mxl",
      coefficients = coefficients,
      vcov = vcov,
      se = se,
      loglik = loglik,
      nobs = nobs,
      n_params = n_params,
      convergence = convergence,
      message = message,
      data_spec = data_spec,
      alt_mapping = alt_mapping,
      param_map = param_map,
      use_asc = use_asc,
      include_outside_option = include_outside_option,
      optimizer = optimizer,
      data = data,
      draws_info = draws_info,
      rc_dist = rc_dist,
      rc_correlation = rc_correlation,
      rc_mean = rc_mean,
      sigma = sigma
    ),
    class = c("choicer_mxl", "choicer_fit")
  )
}

#' Construct a choicer_nl object
#' @inheritParams new_choicer_mnl
#' @param lambda Named numeric vector of nest parameters
#' @returns A choicer_nl object (S3 class)
#' @noRd
new_choicer_nl <- function(call, coefficients, loglik,
                           nobs, n_params, convergence, message,
                           data_spec, alt_mapping, param_map,
                           use_asc, include_outside_option,
                           optimizer,
                           vcov = NULL, se = NULL, data = NULL,
                           lambda = NULL) {
  structure(
    list(
      call = call,
      model = "nl",
      coefficients = coefficients,
      vcov = vcov,
      se = se,
      loglik = loglik,
      nobs = nobs,
      n_params = n_params,
      convergence = convergence,
      message = message,
      data_spec = data_spec,
      alt_mapping = alt_mapping,
      param_map = param_map,
      use_asc = use_asc,
      include_outside_option = include_outside_option,
      optimizer = optimizer,
      data = data,
      lambda = lambda
    ),
    class = c("choicer_nl", "choicer_fit")
  )
}

# --- Lazy vcov computation --------------------------------------------------

#' Compute vcov and SE on demand (lazy evaluation)
#'
#' If `vcov` is already populated, returns the object unchanged.
#' Otherwise, uses `compute_hessian()` to recompute from stored data.
#'
#' @param object A choicer_fit object.
#' @returns The object with `vcov` and `se` populated.
#' @noRd
ensure_vcov <- function(object) {
  if (!is.null(object$vcov)) return(object)

  if (is.null(object$data)) {
    message("Cannot compute standard errors: no data stored. ",
            "Refit with keep_data = TRUE.")
    return(object)
  }

  hess <- compute_hessian(object)
  result <- invert_hessian(hess)

  object$vcov <- result$vcov
  object$se <- result$se

  if (!is.null(object$vcov)) {
    nms <- names(object$coefficients)
    rownames(object$vcov) <- nms
    colnames(object$vcov) <- nms
    names(object$se) <- nms
  }

  object
}

# --- Optimizer adapter -------------------------------------------------------

#' Run optimizer with standardized interface
#'
#' Dispatches to the appropriate optimizer and returns a normalized result.
#' Accepts either a string name ("nloptr", "optim") or a custom function.
#'
#' Custom function interface:
#' \code{my_optimizer(theta_init, eval_f, lower, upper, control)} where
#' \code{eval_f(theta)} returns \code{list(objective, gradient)}.
#' Must return a list with at least \code{par} (or \code{solution}) and
#' \code{value} (or \code{objective}). If the custom function accepts
#' \code{control} or \code{...}, the \code{control} list is forwarded.
#'
#' @param optimizer String name or function. If NULL, defaults to "nloptr".
#' @param theta_init Numeric vector of starting values.
#' @param eval_f Function(theta) returning list(objective, gradient).
#'   Data arguments are captured in the closure by the caller.
#' @param lower Numeric vector of lower bounds (default -Inf).
#' @param upper Numeric vector of upper bounds (default Inf).
#' @param control List of optimizer-specific control parameters.
#' @returns List with: par, value, convergence, message, iterations, raw
#' @noRd
run_optimizer <- function(optimizer, theta_init, eval_f,
                          lower = NULL, upper = NULL,
                          control = list()) {
  if (is.null(optimizer)) optimizer <- "nloptr"

  if (is.function(optimizer)) {
    fmls <- names(formals(optimizer))
    if ("control" %in% fmls || "..." %in% fmls) {
      result <- optimizer(theta_init, eval_f, lower = lower, upper = upper,
                          control = control)
    } else {
      if (length(control) > 0) {
        message("Custom optimizer does not accept 'control'; control arguments ignored.")
      }
      result <- optimizer(theta_init, eval_f, lower = lower, upper = upper)
    }
    return(normalize_optim_result(result))
  }

  if (!is.character(optimizer)) {
    stop("'optimizer' must be a string (e.g., \"nloptr\", \"optim\") or a function.")
  }

  switch(optimizer,
    nloptr = run_nloptr(theta_init, eval_f, lower, upper, control),
    optim = run_optim(theta_init, eval_f, lower, upper, control),
    stop("Unknown optimizer: '", optimizer, "'. ",
         "Supported: \"nloptr\", \"optim\", or pass a custom function.")
  )
}

#' nloptr adapter
#' @noRd
run_nloptr <- function(theta_init, eval_f, lower, upper, control) {
  opts <- list(
    algorithm = "NLOPT_LD_LBFGS",
    xtol_rel = 1.0e-8,
    maxeval = 1000L,
    print_level = 0L
  )
  opts[names(control)] <- control

  args <- list(x0 = theta_init, eval_f = eval_f, opts = opts)
  if (!is.null(lower)) args$lb <- lower
  if (!is.null(upper)) args$ub <- upper

  raw <- do.call(nloptr::nloptr, args)

  list(
    par = raw$solution,
    value = raw$objective,
    convergence = raw$status,
    message = raw$message,
    iterations = raw$iterations,
    raw = raw
  )
}

#' stats::optim adapter
#' @noRd
run_optim <- function(theta_init, eval_f, lower, upper, control) {
  fn <- function(theta) eval_f(theta)$objective
  gr <- function(theta) eval_f(theta)$gradient

  has_bounds <- !is.null(lower) || !is.null(upper)
  method <- if (has_bounds) "L-BFGS-B" else "BFGS"
  if (!is.null(control$method)) {
    method <- control$method
    control$method <- NULL
  }

  args <- list(par = theta_init, fn = fn, gr = gr, method = method,
               control = control)
  if (has_bounds) {
    if (!is.null(lower)) args$lower <- lower
    if (!is.null(upper)) args$upper <- upper
  }

  raw <- do.call(stats::optim, args)

  list(
    par = raw$par,
    value = raw$value,
    convergence = raw$convergence,
    message = raw$message %||% "",
    iterations = raw$counts[["function"]] %||% NA_integer_,
    raw = raw
  )
}

#' Normalize a custom optimizer result to standard fields
#' @noRd
normalize_optim_result <- function(result) {
  # Accept either nloptr-style (solution/objective) or optim-style (par/value)
  par <- result$par %||% result$solution
  value <- result$value %||% result$objective

  if (is.null(par) || is.null(value)) {
    stop("Custom optimizer must return a list with 'par'/'value' ",
         "(or 'solution'/'objective') elements.")
  }

  list(
    par = par,
    value = value,
    convergence = result$convergence %||% result$status %||% NA_integer_,
    message = result$message %||% "",
    iterations = result$iterations %||% NA_integer_,
    raw = result
  )
}

# --- Hessian / vcov helpers --------------------------------------------------

#' Compute Hessian from stored data (no closure needed)
#'
#' Dispatches to the correct C++ Hessian function based on \code{object$model}.
#' For MXL, regenerates Halton draws deterministically from \code{draws_info}.
#'
#' @param object A choicer_fit object with \code{keep_data = TRUE}.
#' @returns The Hessian matrix evaluated at \code{object$coefficients}.
#' @noRd
compute_hessian <- function(object) {
  if (is.null(object$data)) {
    stop("Cannot compute Hessian: no data stored. Refit with keep_data = TRUE.")
  }

  theta <- object$coefficients

  switch(object$model,
    mnl = mnl_loglik_hessian_parallel(
      theta = theta,
      X = object$data$X,
      alt_idx = object$data$alt_idx,
      choice_idx = object$data$choice_idx,
      M = object$data$M,
      weights = object$data$weights,
      use_asc = object$use_asc,
      include_outside_option = object$include_outside_option
    ),
    mxl = {
      eta_draws <- get_halton_normals(
        S = object$draws_info$S,
        N = object$draws_info$N,
        K_w = object$draws_info$K_w
      )
      mxl_hessian_parallel(
        theta = theta,
        X = object$data$X,
        W = object$data$W,
        alt_idx = object$data$alt_idx,
        choice_idx = object$data$choice_idx,
        M = object$data$M,
        weights = object$data$weights,
        eta_draws = eta_draws,
        rc_dist = object$rc_dist,
        rc_correlation = object$rc_correlation,
        rc_mean = object$rc_mean,
        use_asc = object$use_asc,
        include_outside_option = object$include_outside_option
      )
    },
    nl = nl_loglik_numeric_hessian(
      theta = theta,
      X = object$data$X,
      alt_idx = object$data$alt_idx,
      choice_idx = object$data$choice_idx,
      nest_idx = object$data$nest_idx,
      M = object$data$M,
      weights = object$data$weights,
      use_asc = object$use_asc,
      include_outside_option = object$include_outside_option
    ),
    stop("Unknown model type: '", object$model, "'.")
  )
}

#' Invert Hessian to get vcov matrix
#' @param hess Hessian matrix (negative second derivatives)
#' @returns List with vcov (matrix or NULL) and se (numeric vector)
#' @noRd
invert_hessian <- function(hess) {
  p_len <- nrow(hess)
  vcov_mat <- NULL
  se <- rep(NA_real_, p_len)

  singular_flag <- FALSE
  tryCatch({
    vcov_mat <- solve(hess)
  }, error = function(e) {
    singular_flag <<- TRUE
    message("Error computing vcov (likely singular Hessian): ", e$message)
  })

  if (!singular_flag && !is.null(vcov_mat)) {
    se <- sqrt(diag(vcov_mat))
  } else {
    message("Standard errors set to NA due to Hessian inversion failure.")
  }

  list(vcov = vcov_mat, se = se)
}
