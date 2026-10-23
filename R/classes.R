
# S3 class constructors for choicer model objects
#
# Class hierarchy:
#   choicer_mnl -> choicer_fit
#   choicer_mxl -> choicer_fit
#   choicer_nl  -> choicer_fit
#   choicer_mnp (standalone: posterior-draws object, no loglik / convergence /
#                lazy-Hessian contract, so it does not inherit choicer_fit)

# --- Constructors -----------------------------------------------------------

#' Construct a choicer_mnl object
#' @param call The matched call from run_mnlogit()
#' @param coefficients Named numeric vector of point estimates
#' @param loglik Scalar log-likelihood at optimum
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
#' @param scale_vars Character. How the optimizer's coordinates were scaled:
#'   \code{"none"} (default), \code{"sd"}, \code{"mad"}, \code{"iqr"}, or
#'   \code{"bhhh"}.
#' @param sX Named numeric vector of the column scales of X (\code{"sd"},
#'   \code{"mad"}, \code{"iqr"}); 1s otherwise.
#' @param se_method Character. Method used for standard errors: \code{"hessian"}
#'   (analytical Hessian, default), \code{"bhhh"} (outer product of gradients),
#'   or \code{"sandwich"} (robust Huber--White / WESML variance).
#' @param choice_sampling Optional list recording choice-based-sampling
#'   provenance (scheme, population/sample shares, meat type), or NULL.
#' @param param_scale,param_shift Named numeric vectors (m, c) of the map from
#'   the optimizer's coordinates to the parameters, theta = m * theta_t + c:
#'   1s and 0s without scaling.
#' @returns A choicer_mnl object (S3 class)
#' @noRd
new_choicer_mnl <- function(call, coefficients, loglik,
                            nobs, n_params, convergence, message,
                            data_spec, alt_mapping, param_map,
                            use_asc, include_outside_option,
                            optimizer,
                            vcov = NULL, se = NULL, data = NULL,
                            scale_vars = "none", sX = NULL,
                            se_method = "hessian",
                            choice_sampling = NULL,
                            param_scale = NULL, param_shift = NULL) {
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
      data = data,
      scale_vars = scale_vars,
      sX = sX,
      se_method = se_method,
      choice_sampling = choice_sampling,
      param_scale = param_scale,
      param_shift = param_shift
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
#' @param se_method Character. Method used for standard errors: "hessian"
#'   (analytical Hessian, default) or "bhhh" (outer product of gradients).
#' @param scale_vars Character. How the optimizer's coordinates were scaled:
#'   \code{"none"} (default), \code{"sd"}, \code{"mad"}, \code{"iqr"}, or
#'   \code{"bhhh"}.
#' @param sX Named numeric vector of the column scales of \code{X}
#'   (\code{scale_vars} \code{"sd"}, \code{"mad"} or \code{"iqr"}; e.g.
#'   \code{apply(X, 2, sd)}); 1s otherwise.
#' @param sW Named numeric vector of the column scales of \code{W}, as
#'   \code{sX}, with entries for log-normal columns (\code{rc_dist[k] == 1})
#'   carved out to 1; 1s otherwise.
#' @param choice_sampling Optional list recording choice-based-sampling
#'   provenance (scheme, population/sample shares, meat type), or NULL.
#' @param n_persons Number of decision makers of a panel fit (`person_col`),
#'   or NULL for a cross-sectional fit.
#' @param param_scale,param_shift Named numeric vectors (m, c) of the map from
#'   the optimizer's coordinates to the parameters, theta = m * theta_t + c:
#'   1s and 0s without scaling.
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
                            rc_mean = FALSE, sigma = NULL,
                            se_method = "hessian",
                            scale_vars = "none",
                            sX = NULL, sW = NULL,
                            choice_sampling = NULL,
                            n_persons = NULL,
                            param_scale = NULL, param_shift = NULL) {
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
      sigma = sigma,
      se_method = se_method,
      scale_vars = scale_vars,
      sX = sX,
      sW = sW,
      choice_sampling = choice_sampling,
      n_persons = n_persons,
      param_scale = param_scale,
      param_shift = param_shift
    ),
    class = c("choicer_mxl", "choicer_fit")
  )
}

#' Construct a choicer_nl object
#' @inheritParams new_choicer_mnl
#' @param lambda Named numeric vector of nest parameters
#' @param nest_idx Integer vector of length J mapping each (inside) alternative,
#'   in alt_mapping row order, to its 1-based nest index. Stored top-level so
#'   newdata prediction works even when data = NULL (keep_data = FALSE).
#' @param se_method Character. Method used for standard errors: \code{"hessian"}
#'   (analytical Hessian, default), \code{"numeric"} (finite-difference oracle),
#'   \code{"bhhh"} (outer product of gradients), or \code{"sandwich"} (robust
#'   Huber--White / WESML variance).
#' @param choice_sampling Optional list recording choice-based-sampling
#'   provenance (scheme, population/sample shares, meat type), or NULL.
#' @param scale_vars,sX As for \code{new_choicer_mnl()}.
#' @param param_scale,param_shift Named numeric vectors (m, c) of the map from
#'   the optimizer's coordinates to the parameters, theta = m * theta_t + c:
#'   1s and 0s without scaling.
#' @returns A choicer_nl object (S3 class)
#' @noRd
new_choicer_nl <- function(call, coefficients, loglik,
                           nobs, n_params, convergence, message,
                           data_spec, alt_mapping, param_map,
                           use_asc, include_outside_option,
                           optimizer,
                           vcov = NULL, se = NULL, data = NULL,
                           lambda = NULL, nest_idx = NULL,
                           se_method = "hessian",
                           choice_sampling = NULL,
                           scale_vars = "none", sX = NULL,
                           param_scale = NULL, param_shift = NULL) {
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
      lambda = lambda,
      nest_idx = nest_idx,
      se_method = se_method,
      choice_sampling = choice_sampling,
      scale_vars = scale_vars,
      sX = sX,
      param_scale = param_scale,
      param_shift = param_shift
    ),
    class = c("choicer_nl", "choicer_fit")
  )
}

#' Construct a choicer_mnp object
#' @param call The matched call from run_mnprobit()
#' @param coefficients Named numeric vector of posterior means of the
#'   identified coefficients (beta / sqrt(sigma_11))
#' @param se Named numeric vector of posterior standard deviations of the
#'   identified coefficient draws
#' @param vcov Posterior covariance matrix of the identified coefficient draws
#' @param sigma Posterior-mean identified covariance matrix of the utility
#'   differences (p x p)
#' @param draws List of draw matrices: beta / sigma (identified scale) and
#'   beta_raw / sigma_raw (unnormalized chain output)
#' @param prior List of resolved prior settings (beta_bar, A, nu, V)
#' @param mcmc List of resolved MCMC settings (R, burn, thin, seed, R_keep)
#' @param nobs Integer number of choice situations
#' @param n_params Integer number of coefficients
#' @param data_spec List with column name metadata (id_col, alt_col, etc.)
#' @param alt_mapping data.table mapping alternatives to summary statistics
#' @param base_alt Label of the base (differencing) alternative
#' @param param_map Named list of integer index vectors (beta, asc)
#' @param use_asc Logical whether ASCs were used
#' @param sampler List with sampler metadata (name, elapsed_time)
#' @param data List of prepared inputs (X, y, p) or NULL
#' @returns A choicer_mnp object (S3 class)
#' @noRd
new_choicer_mnp <- function(call, coefficients, se, vcov, sigma, draws,
                            prior, mcmc, nobs, n_params,
                            data_spec, alt_mapping, base_alt,
                            param_map, use_asc, sampler, data = NULL) {
  structure(
    list(
      call = call,
      model = "mnp",
      coefficients = coefficients,
      se = se,
      vcov = vcov,
      sigma = sigma,
      draws = draws,
      prior = prior,
      mcmc = mcmc,
      nobs = nobs,
      n_params = n_params,
      data_spec = data_spec,
      alt_mapping = alt_mapping,
      base_alt = base_alt,
      param_map = param_map,
      use_asc = use_asc,
      sampler = sampler,
      data = data
    ),
    # Intentionally not a choicer_fit: this is a posterior-draws object with
    # no loglik / convergence and an eagerly computed vcov, so none of the
    # choicer_fit methods (logLik, AIC, ensure_vcov, predict, ...) apply.
    class = "choicer_mnp"
  )
}

#' Construct a choicer_hmnl / choicer_hmnp object
#'
#' Internal constructors for the hierarchical Bayes fits. Both share the
#' `choicer_hb` parent class, which carries the print / summary / coef /
#' vcov / nobs methods (R/methods.R) and the post-estimation suite. Like
#' `choicer_mnp`, these are posterior-draws objects, intentionally NOT
#' `choicer_fit`.
#'
#' @param ... Named fields; see [run_hmnlogit()] for the layout.
#' @returns A classed list.
#' @noRd
new_choicer_hmnl <- function(call, coefficients, se, vcov, theta_summary,
                             sigma_d2_summary, W_mean, delta, xi, beta_i,
                             draws, accept, rhat, prior, mcmc, nobs,
                             n_persons, J, P, K_struct, param_map,
                             alt_mapping, data_spec, cf_active, sampler,
                             data = NULL, chains = NULL) {
  structure(
    list(
      call = call,
      model = "hmnl",
      coefficients = coefficients,
      se = se,
      vcov = vcov,
      theta_summary = theta_summary,
      sigma_d2_summary = sigma_d2_summary,
      W_mean = W_mean,
      delta = delta,
      xi = xi,
      beta_i = beta_i,
      draws = draws,
      accept = accept,
      rhat = rhat,
      prior = prior,
      mcmc = mcmc,
      nobs = nobs,
      n_persons = n_persons,
      J = J,
      P = P,
      K_struct = K_struct,
      param_map = param_map,
      alt_mapping = alt_mapping,
      data_spec = data_spec,
      cf_active = cf_active,
      sampler = sampler,
      data = data,
      chains = chains
    ),
    class = c("choicer_hmnl", "choicer_hb")
  )
}

#' @rdname new_choicer_hmnl
#' @noRd
new_choicer_hmnp <- function(call, coefficients, se, vcov, theta_summary,
                             sigma_d2_summary, W_mean, delta, xi, beta_i,
                             draws, rhat, prior, mcmc, nobs,
                             n_persons, J, P, K_struct, param_map,
                             alt_mapping, data_spec, cf_active, sampler,
                             data = NULL, chains = NULL) {
  structure(
    list(
      call = call,
      model = "hmnp",
      coefficients = coefficients,
      se = se,
      vcov = vcov,
      theta_summary = theta_summary,
      sigma_d2_summary = sigma_d2_summary,
      W_mean = W_mean,
      delta = delta,
      xi = xi,
      beta_i = beta_i,
      draws = draws,
      accept = NULL,   # fully conjugate: no Metropolis steps
      rhat = rhat,
      prior = prior,
      mcmc = mcmc,
      nobs = nobs,
      n_persons = n_persons,
      J = J,
      P = P,
      K_struct = K_struct,
      param_map = param_map,
      alt_mapping = alt_mapping,
      data_spec = data_spec,
      cf_active = cf_active,
      sampler = sampler,
      data = data,
      chains = chains
    ),
    class = c("choicer_hmnp", "choicer_hb")
  )
}

# --- Lazy vcov computation --------------------------------------------------

#' Compute vcov and SE on demand (lazy evaluation)
#'
#' If `vcov` is already populated, or the fit's variance failed at fit time
#' (its `variance_error`, see `.vcov_failed()`), returns the object
#' unchanged: the failed step is run again only on request, by
#' `vcov(fit, type = )`.
#' Otherwise, uses `compute_hessian()` to recompute from stored data.
#'
#' @param object A choicer_fit object.
#' @returns The object with `vcov` and `se` populated.
#' @noRd
ensure_vcov <- function(object) {
  if (!is.null(object$vcov) || !is.null(object[["variance_error"]])) {
    return(object)
  }

  if (is.null(object[["data"]])) {
    message("Cannot compute standard errors: no data stored. ",
            "Refit with keep_data = TRUE.")
    return(object)
  }

  se_method <- object$se_method %||% "hessian"
  result <- if (identical(se_method, "sandwich")) {
    compute_sandwich_vcov(object)
  } else if (identical(se_method, "cluster")) {
    # cluster = NULL routes to the stored, already-aligned fit-time labels.
    .assemble_score_vcov(object, type = "cluster")
  } else {
    invert_hessian(compute_hessian(object))
  }

  # Named before the fit holds it: naming a matrix that two objects share
  # copies it (one p x p matrix, 2.3 GB at 17,000 parameters).
  if (!is.null(result$vcov)) {
    nms <- names(object$coefficients)
    dimnames(result$vcov) <- list(nms, nms)
    names(result$se) <- nms
  }
  object$vcov <- result$vcov
  object$se <- result$se

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

#' Keep the last evaluation of the objective
#'
#' `stats::optim()` takes the objective (`fn`) and its gradient (`gr`) from
#' separate functions, and asks for the gradient at the point whose objective
#' it has just computed (L-BFGS-B's line search can also ask again for the
#' point it has just evaluated); `eval_f` computes both in one kernel pass.
#' The result answers a call at the last point evaluated from that
#' evaluation, and evaluates any other point afresh: BFGS's line search asks
#' for the objective alone at trial points, CG for the gradient back at its
#' line search's point once an interpolation step fails. Points match bit for
#' bit, names included (`identical(num.eq = FALSE, single.NA = FALSE)`: -0 is
#' not 0). A cached answer stands for a repeat of the call before it, which
#' leaves `.lift_sentinel()`'s running maximum as it was; a cache of earlier
#' points would not be exact, since the sentinel's value at a point that
#' overflows rises with that maximum.
#'
#' @param eval_f Function of theta returning `list(objective, gradient)`.
#' @returns The wrapped function.
#' @noRd
.one_point_cache <- function(eval_f) {
  force(eval_f)
  last_theta <- NULL
  last <- NULL
  function(theta) {
    if (!identical(theta, last_theta, num.eq = FALSE, single.NA = FALSE)) {
      res <- eval_f(theta)
      last_theta <<- theta
      last <<- res
    }
    last
  }
}

#' stats::optim adapter
#'
#' The objective and the gradient at the point optim has just evaluated share
#' one kernel pass (`.one_point_cache()`).
#' @noRd
run_optim <- function(theta_init, eval_f, lower, upper, control) {
  eval_f <- .one_point_cache(eval_f)
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

# --- The Newton phase (polish) -----------------------------------------------

#' Settings of the Newton phase
#'
#' `polish_control` of `run_mxlogit()`: a named list whose entries replace
#' the defaults, `max_iter` (20; a whole number, 0 certifies the start
#' without stepping), `tol_decrement` (1e-8, in log-likelihood units) and
#' `tol_step` (1e-12, relative). A NULL entry keeps the default.
#'
#' @param control A named list, or NULL.
#' @returns The complete settings.
#' @noRd
.polish_control <- function(control) {
  ctl <- list(max_iter = 20L, tol_decrement = 1e-8, tol_step = 1e-12)
  if (is.null(control)) return(ctl)
  if (!is.list(control)) stop("`polish_control` must be a list.", call. = FALSE)
  nms <- names(control)
  if (length(control) && (is.null(nms) || anyNA(nms) || !all(nzchar(nms)))) {
    stop("`polish_control` must be a named list.", call. = FALSE)
  }
  bad <- setdiff(nms, names(ctl))
  if (length(bad)) {
    stop("Unknown `polish_control` setting(s): ", paste(bad, collapse = ", "),
         ". Valid: ", paste(names(ctl), collapse = ", "), ".", call. = FALSE)
  }
  dup <- unique(nms[duplicated(nms)])
  if (length(dup)) {
    stop("Duplicate `polish_control` setting(s): ", paste(dup, collapse = ", "),
         ".", call. = FALSE)
  }
  control <- control[!vapply(control, is.null, TRUE)]
  ctl[names(control)] <- control
  m <- ctl$max_iter
  if (!is.numeric(m) || length(m) != 1L || !is.finite(m) || m < 0 ||
      m != round(m) || m > .Machine$integer.max) {
    stop("`polish_control$max_iter` must be a whole number >= 0.", call. = FALSE)
  }
  ctl$max_iter <- as.integer(m)
  for (nm in c("tol_decrement", "tol_step")) {
    v <- ctl[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop("`polish_control$", nm, "` must be a positive number.", call. = FALSE)
    }
  }
  ctl
}

#' Polish an optimizer's solution by damped Newton steps
#'
#' The Newton phase of `run_mxlogit(polish = )`. From the optimizer's
#' solution it takes damped Newton steps on the curvature `curv_f` returns,
#' the information (the Hessian of the negated log-likelihood) or the BHHH
#' matrix, until the Newton decrement, half of g' C^-1 g with g the gradient
#' and C the curvature, falls below `tol_decrement`. On a positive definite
#' Hessian the decrement estimates the log-likelihood the remaining Newton
#' step would gain, and sqrt(2 decrement) is the distance to the maximum of
#' the quadratic model in Hessian standard errors (Cauchy-Schwarz in the
#' C^-1 inner product bounds every linear combination's), so a converged
#' polish certifies the maximum the optimizer approached. BHHH is positive
#' semi-definite everywhere: where the Hessian H is positive definite its
#' sqrt(2 decrement) is that distance in the standard errors of the sandwich
#' H^-1 B H^-1, but on its own it certifies only a stationary point. The
#' estimator does not change.
#'
#' Each iteration works in natural units. Damped by Marquardt's diagonal the
#' Newton step does not depend on a diagonal rescaling of the parameters, to
#' within the factor of two of D's rounding; the cap, epsilon and `tol_step`
#' below are measured in the fit's coordinates (m, `scale`).
#' * The curvature at the current point, scaled by its Jacobi scale D
#'   (`.jacobi_scale()`, powers of two). The previous curvature is released
#'   before the next is computed.
#' * Bounds, by two-metric projection (Bertsekas 1982): a parameter within
#'   epsilon of a bound that its gradient points across is epsilon-active
#'   and moves to P(theta - D^2 g), onto the bound; epsilon =
#'   m min(1e-3, max |w / m|) over the parameters not held, with
#'   w = theta - P(theta - m^2 g) the projected gradient step, shrinks to
#'   zero at a stationary point. The other parameters are free, except those
#'   `fixed` and the flat ones, whose curvature column and gradient are zero
#'   (no information at this point, as where a log-Cholesky diagonal has gone
#'   so far toward -Inf that its derivatives underflow).
#' * A Newton step on the free set (`newton_step()`), (D C_FF D + lambda I) y
#'   = -D g_F, undamped (lambda = 0) when the free block's diagonal is
#'   positive and its Cholesky factorization succeeds with a reciprocal
#'   condition number above the machine epsilon. Otherwise, and when the step
#'   is not finite or its line search fails, Levenberg-Marquardt damping from
#'   max(1e-6, lambda_prev / 100) in factors of 10 up to 1e10; lambda I on the
#'   scaled block is Marquardt's lambda |diag C|, to within a factor of two.
#' * A free parameter whose step is longer than ten maximum steps (below),
#'   100 m, on a curvature that has all but vanished (|C_kk| m_k^2 below 1e-4
#'   times the free parameters' median) is held for the iteration (stalled)
#'   and the step solved again without it. BHHH's curvature vanishes so, like
#'   exp(4 L), on a log-Cholesky diagonal L whose variance goes to zero:
#'   shortened to the maximum step, its step would shorten every other
#'   parameter's by more than ten, freezing them. A step that is long because
#'   the start is far from the maximum is shortened instead.
#' * Before any step: converged when the undamped decrement is below
#'   `tol_decrement`, the epsilon-active parameters sit on their bounds and no
#'   parameter is flat or stalled (so the certificate is for the point
#'   returned; a flat or stalled parameter stops the polish unconverged,
#'   named in the message); then the iteration limit; then a step below
#'   `tol_step` relative to |theta| + m.
#' * A free step longer than 10 units of the fit's coordinates (|step| / m),
#'   the maximum step, is shortened to that length, direction kept (as in
#'   Dennis and Schnabel's methods and `nlm()`'s `stepmax`). Near a maximum
#'   the steps are far shorter.
#' * A line search on the projected path P(theta + alpha step), halving
#'   alpha up to 40 times (also where the projected step is not yet a descent
#'   direction, which it becomes for small alpha): Armijo's condition with
#'   1e-4; or, only where the predicted gain -g'd is within the objective's
#'   rounding allowance eps_f = 1e-10 max(|f|, 1), the approximate Armijo
#'   condition of Hager and Zhang (2005), the upper half of their approximate
#'   Wolfe conditions: f(alpha) <= f + eps_f with g(alpha)'d <= -0.8 g'd.
#'   Near the maximum of a likelihood summed over millions of units, rounding
#'   in f can exceed the gain of a correct step; the gradient carries no such
#'   noise. A trial point whose objective overflows (`overflow`), or whose
#'   objective or gradient is not finite, fails.
#'
#' An error from `eval_f`, `curv_f` or the step (a kernel's memory check, an
#' allocation) ends the polish at the last accepted point, with the error's
#' message in the result. A point worse than the start by more than eps_f
#' (gated steps of the second rule that add up) is not returned: the start
#' is, while `iterations` and the rest describe the path abandoned.
#'
#' Every stop is decided at the point whose curvature was evaluated last:
#' before a step is tried, or after its line search or damping failed, or on
#' an error in a trial evaluation, all at that point. So the last curvature is
#' at the point returned, unless the final guard returns the start after
#' steps, `curv_f` itself failed, the start could not be evaluated (an
#' error, or an objective or gradient that is not finite there), or `stop_f`
#' stopped the polish after a step, before the curvature there. With
#' `keep_curvature` the result then carries it, for a caller whose
#' variance needs that matrix at that point (the fit's Hessian or BHHH
#' matrix); without, the matrix is released with the polish.
#'
#' @param theta Start, the optimizer's solution, in natural units.
#' @param eval_f Function of theta returning `list(objective, gradient,
#'   overflow)`: the negated log-likelihood, its gradient, and whether the
#'   objective overflowed.
#' @param curv_f Function of theta returning the curvature matrix.
#' @param lower,upper Bounds, or NULL.
#' @param scale The fit's coordinate scales m (`param_scale`), or NULL (1s).
#' @param control From `.polish_control()`.
#' @param method The curvature's name for the result: "newton" or "bhhh".
#' @param fixed Positions of the parameters held at their values.
#' @param keep_curvature Whether to return the last curvature (see above).
#' @param stop_f NULL, or a function of theta and the gradient (the point
#'   reached and its gradient, after the start's evaluation and after every
#'   accepted step) returning NULL to go on, or the reason to stop converged
#'   before the curvature is evaluated there (`decrement` is then NA): a
#'   caller's own convergence rule beside the decrement (the profile
#'   likelihood's residual rule).
#' @returns List: `method`; `par` (named as `theta`) and `value`; `converged`;
#'   `iterations` (steps taken) and `damped` (those with lambda > 0);
#'   `decrement`, the undamped decrement at the last curvature evaluated, over
#'   the parameters not stalled (NA where the undamped block was not positive
#'   definite, or singular to working precision, and after a stop by
#'   `stop_f`), and `lambda`, the damping of
#'   the last step tried; `message`, why the polish stopped;
#'   `evaluations`, the objective's and the curvature's calls; `active`,
#'   `held` and `flat`, the names of the parameters on a bound, held, and
#'   left out of the last iteration's steps (flat or stalled); `improvement`,
#'   the objective's decrease; `error`, NULL or the message of the error that
#'   stopped it; `curvature_at_par`, with `keep_curvature` the last curvature
#'   evaluated when it was evaluated at `par` (bit for bit), otherwise NULL.
#' @references Bertsekas, D. P. (1982). Projected Newton methods for
#'   optimization problems with simple constraints. SIAM Journal on Control
#'   and Optimization 20(2), 221-246.
#'
#'   Dennis, J. E. and Schnabel, R. B. (1983). Numerical Methods for
#'   Unconstrained Optimization and Nonlinear Equations. Prentice-Hall.
#'
#'   Hager, W. W. and Zhang, H. (2005). A new conjugate gradient method with
#'   guaranteed descent and an efficient line search. SIAM Journal on
#'   Optimization 16(1), 170-192.
#' @noRd
run_polish <- function(theta, eval_f, curv_f, lower = NULL, upper = NULL,
                       scale = NULL, control = .polish_control(NULL),
                       method = "newton", fixed = integer(),
                       keep_curvature = FALSE, stop_f = NULL) {
  n <- length(theta)
  nm <- names(theta) %||% as.character(seq_len(n))
  lower <- lower %||% rep(-Inf, n)
  upper <- upper %||% rep(Inf, n)
  m <- scale %||% rep(1, n)
  clip <- function(x) pmin(pmax(x, lower), upper)
  rel <- function(d, x, mm = m) max(abs(d) / (abs(x) + mm))
  held <- seq_len(n) %in% fixed
  dname <- if (identical(method, "bhhh")) "BHHH decrement" else "Newton decrement"
  armijo <- 1e-4
  wolfe <- 0.1
  halvings <- 40L
  lambda0 <- 1e-6
  lambda_max <- 1e10
  max_step <- 10
  stall_ratio <- 10 * max_step
  collapse <- 1e-4

  n_f <- 0L
  n_c <- 0L
  evaluate <- function(x) {
    n_f <<- n_f + 1L
    eval_f(x)
  }
  theta_in <- theta
  theta <- clip(theta)
  theta0 <- theta
  f <- NA_real_
  f0 <- NA_real_
  g <- NULL
  eps_f <- NA_real_
  iter <- 0L
  damped <- 0L
  decrement <- NA_real_
  lambda <- NA_real_
  lambda_prev <- 0
  converged <- FALSE
  flat <- integer()
  stalled <- integer()
  # The last curvature and the point it was evaluated at (run() assigns
  # them here)
  C <- NULL
  theta_C <- NULL

  # A trial point P(theta + alpha step), alpha halved until a rule accepts
  line_search <- function(step) {
    alpha <- 1
    for (k in 0:halvings) {
      x <- clip(theta + alpha * step)
      d <- x - theta
      if (rel(d, theta) < control$tol_step) return(NULL)
      gd <- sum(g * d)
      if (isTRUE(gd < 0)) {
        e <- evaluate(x)
        if (!isTRUE(e$overflow) && is.finite(e$objective) &&
            all(is.finite(e$gradient)) &&
            (e$objective <= f + armijo * gd ||
             (-gd <= eps_f && e$objective <= f + eps_f &&
              sum(e$gradient * d) <= -(1 - 2 * wolfe) * gd))) {
          return(list(theta = x, f = e$objective, g = e$gradient))
        }
      }
      alpha <- alpha / 2
    }
    NULL
  }

  # Steps until a stop; returns why it stopped
  run <- function() {
    e <- evaluate(theta)
    if (isTRUE(e$overflow) || !is.finite(e$objective)) {
      return("the objective is not finite at the start")
    }
    if (!all(is.finite(e$gradient))) {
      return("the gradient is not finite at the start")
    }
    f0 <<- f <<- e$objective
    g <<- e$gradient
    eps_f <<- 1e-10 * max(abs(f0), 1)
    repeat {
      if (!is.null(stop_f)) {
        why <- stop_f(theta, g)
        if (!is.null(why)) {
          # no curvature at this point: no decrement either
          decrement <<- NA_real_
          converged <<- TRUE
          return(why)
        }
      }
      # the previous curvature is released before the next is computed
      C <<- NULL
      n_c <<- n_c + 1L
      C <<- curv_f(theta)
      theta_C <<- theta
      D <- .jacobi_scale(C)
      dC <- diag(C)
      w <- theta - clip(theta - m^2 * g)
      eps <- m * min(1e-3, max(c(0, abs(w[!held] / m[!held]))))
      act <- !held & ((theta - lower <= eps & g > 0) |
                        (upper - theta <= eps & g < 0))
      zero <- which(!held & dC == 0 & g == 0)
      flat <<- zero[vapply(zero, function(k) isTRUE(all(C[, k] == 0)), TRUE)]
      stalled <<- integer()
      decrement <<- NA_real_
      free <- which(!held & !act & !seq_len(n) %in% flat)
      lam <- if (isTRUE(all(dC[free] > 0))) 0 else max(lambda0, lambda_prev / 100)
      repeat {
        repeat {
          r <- newton_step(C, g, D, free, lam)
          if (r$info != 0L) break
          # a free parameter whose step runs past ten maximum steps on a
          # curvature that has all but vanished (in the fit's coordinates,
          # 1e4 times below the free parameters' median) is held, the step
          # re-solved
          kap <- abs(dC[free]) * m[free]^2
          off <- free[which(abs(r$step[free]) / m[free] > stall_ratio &
                              kap < collapse * stats::median(kap))]
          if (!length(off)) break
          stalled <<- c(stalled, off)
          free <- setdiff(free, off)
        }
        if (r$info < 0L) return("the curvature is not finite")
        ok <- r$info == 0L
        if (ok) {
          step <- r$step
          step[act] <- (clip(theta - D^2 * g) - theta)[act]
          dec <- -0.5 * sum(g[free] * step[free])
          ok <- all(is.finite(step)) && is.finite(dec)
        }
        if (ok) {
          lambda <<- lam
          if (lam == 0) decrement <<- dec
          if (lam == 0 && dec < control$tol_decrement &&
              (!any(act) ||
               rel(step[act], theta[act], m[act]) < control$tol_step)) {
            out <- c(flat, stalled)
            if (length(out)) {
              return(paste0(
                if (length(free)) {
                  paste("the", dname, "over the other parameters is below",
                        "tol_decrement, but ")
                },
                paste(nm[sort(out)], collapse = ", "),
                if (length(out) == 1L) " has" else " have",
                " no usable curvature at this point (",
                if (length(flat)) "no information",
                if (length(flat) && length(stalled)) "; ",
                if (length(stalled)) {
                  "a curvature all but vanished, a step more than 100 times its scale"
                },
                ")"))
            }
            converged <<- TRUE
            return(paste("the", dname, "is below tol_decrement"))
          }
          if (iter >= control$max_iter) {
            return("the iteration limit (max_iter) is reached")
          }
          if (rel(step, theta) < control$tol_step) {
            return("the step is below tol_step")
          }
          big <- max(c(0, abs(step[free]) / m[free]))
          if (big > max_step) step[free] <- step[free] * (max_step / big)
          ls <- line_search(step)
          if (!is.null(ls)) {
            theta <<- ls$theta
            f <<- ls$f
            g <<- ls$g
            iter <<- iter + 1L
            if (lam > 0) damped <<- damped + 1L
            lambda_prev <<- lam
            break
          }
        }
        lam <- if (lam == 0) max(lambda0, lambda_prev / 100) else 10 * lam
        if (lam > lambda_max) {
          return(if (ok) "the line search failed"
                 else "no damping makes the curvature positive definite")
        }
      }
    }
  }

  err <- NULL
  why <- tryCatch(run(), error = function(e) {
    err <<- conditionMessage(e)
    paste0("an error: ", err)
  })
  if (!is.finite(f)) {
    # nothing was accepted, not even the start
    theta <- theta_in
  } else if (f > f0 + eps_f) {
    theta <- theta0
    f <- f0
    converged <- FALSE
    why <- paste0(why, "; the last point was worse than the start, which ",
                  "is returned")
  }
  names(theta) <- names(theta_in)
  # The last curvature, when it was evaluated at the point returned
  curvature_at_par <- if (keep_curvature && !is.null(C) &&
                          identical(unname(theta_C), unname(theta), num.eq = FALSE)) C
  list(method = method, par = theta, value = f, converged = converged,
       iterations = iter, damped = damped, decrement = decrement,
       lambda = lambda, message = why,
       evaluations = c(objective = n_f, curvature = n_c),
       active = nm[theta <= lower | theta >= upper],
       held = nm[held], flat = nm[sort(c(flat, stalled))],
       improvement = f0 - f, error = err, curvature_at_par = curvature_at_par)
}

# --- The profile likelihood of the constants (concentrate_out_asc) -----------

#' Settings of the profile likelihood's inner solve
#'
#' `concentrate_control` of `run_mxlogit()`: a named list whose entries
#' replace the defaults, `max_iter` (25, a whole number >= 1: Newton steps per
#' inner solve; the first solve, from a cold start, may take four times as
#' many), `tol_residual` (1e-10: the solve stops when every constant's score,
#' and the reference's, is below it relative to its weighted chosen count plus
#' the mean weight), `tol_decrement` (1e-12: or when the Newton decrement over
#' the mean weight is) and `max_dense` (10000: the most constants whose block
#' is factored densely). A NULL entry keeps the default.
#'
#' @param control A named list, or NULL.
#' @returns The complete settings.
#' @noRd
.concentrate_control <- function(control) {
  ctl <- list(max_iter = 25L, tol_residual = 1e-10, tol_decrement = 1e-12,
              max_dense = 10000L)
  if (is.null(control)) return(ctl)
  if (!is.list(control)) stop("`concentrate_control` must be a list.", call. = FALSE)
  nms <- names(control)
  if (length(control) && (is.null(nms) || anyNA(nms) || !all(nzchar(nms)))) {
    stop("`concentrate_control` must be a named list.", call. = FALSE)
  }
  bad <- setdiff(nms, names(ctl))
  if (length(bad)) {
    stop("Unknown `concentrate_control` setting(s): ", paste(bad, collapse = ", "),
         ". Valid: ", paste(names(ctl), collapse = ", "), ".", call. = FALSE)
  }
  dup <- unique(nms[duplicated(nms)])
  if (length(dup)) {
    stop("Duplicate `concentrate_control` setting(s): ", paste(dup, collapse = ", "),
         ".", call. = FALSE)
  }
  control <- control[!vapply(control, is.null, TRUE)]
  ctl[names(control)] <- control
  for (nm in c("max_iter", "max_dense")) {
    m <- ctl[[nm]]
    if (!is.numeric(m) || length(m) != 1L || !is.finite(m) || m < 1 ||
        m != round(m) || m > .Machine$integer.max) {
      stop("`concentrate_control$", nm, "` must be a whole number >= 1.", call. = FALSE)
    }
    ctl[[nm]] <- as.integer(m)
  }
  for (nm in c("tol_residual", "tol_decrement")) {
    v <- ctl[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop("`concentrate_control$", nm, "` must be a positive number.", call. = FALSE)
    }
  }
  ctl
}

#' The profile objective for the outer optimizer
#'
#' The negated profile log-likelihood -l*(theta_c) = -max_delta l(theta_c,
#' delta) as a function of the continuous parameters, its gradient by the
#' envelope theorem (the gradient of the joint objective at (theta_c,
#' delta_hat(theta_c)), restricted to theta_c), and the bookkeeping that keeps
#' the inner solves safe. Generic: `solve(theta_c, delta0)` is the model's
#' inner solve.
#' * Warm starts come from the best outer point so far, the one with the
#'   lowest objective among those whose inner solve converged. The optimizer
#'   does not say which evaluations its line search accepted, but under a
#'   monotone line search the best point is the current iterate, or a better
#'   trial on its line; a rejected trial, being worse, never seeds the next
#'   solve (the audit's experiment, warm-started from every evaluated point,
#'   let the constants drift by 66.5 against their pinned reference).
#' * A theta_c bitwise equal to the best point returns its stored result
#'   without solving.
#' * A solve that does not converge returns the overflow sentinel, so that
#'   the optimizer backtracks; it updates neither the warm start nor the
#'   cache: at an unconverged delta the objective is only an upper bound, and
#'   the gradient not the profile's. Such solves are counted, in all and
#'   since the best point last improved: early ones are the line search's
#'   ordinary backtracking, late ones may have stopped the optimizer.
#'
#' @param solve Function(theta_c, delta0) returning `list(converged, value,
#'   delta, gradient, iterations, passes, residual, decrement)`: the negated
#'   joint log-likelihood at the solution, the constants there, and the full
#'   gradient there.
#' @param c_idx Positions of theta_c in the full parameter vector.
#' @param start The solved start: `list(theta_c, value, delta, gradient,
#'   residual, decrement)`.
#' @returns List: `f`, the objective of theta_c returning `list(objective,
#'   gradient, overflow)`; `state`, an environment holding `best` and the
#'   counts (`evaluations`, `solves`, `failed`, `failed_late` (since the best
#'   point last improved), `hits`, `inner_iterations`, `passes`,
#'   `max_residual`), and `error`, the last error that ended a solve since
#'   the best point last improved.
#' @noRd
.profile_objective <- function(solve, c_idx, start) {
  st <- new.env(parent = emptyenv())
  st$best <- start
  st$evaluations <- 0L
  st$solves <- 0L
  st$failed <- 0L
  st$failed_late <- 0L
  st$error <- NULL
  st$hits <- 0L
  st$inner_iterations <- 0L
  st$passes <- 0L
  st$max_residual <- 0
  f <- function(theta_c) {
    st$evaluations <- st$evaluations + 1L
    b <- st$best
    if (identical(theta_c, b$theta_c, num.eq = FALSE)) {
      st$hits <- st$hits + 1L
      return(list(objective = b$value, gradient = b$gradient[c_idx],
                  overflow = FALSE))
    }
    s <- solve(theta_c, b$delta)
    st$solves <- st$solves + 1L
    st$inner_iterations <- st$inner_iterations + s$iterations
    st$passes <- st$passes + s$passes
    if (!isTRUE(s$converged)) {
      st$failed <- st$failed + 1L
      st$failed_late <- st$failed_late + 1L
      if (!is.null(s$error)) st$error <- s$error
      return(list(objective = 1e10, gradient = rep(0, length(c_idx)),
                  overflow = TRUE))
    }
    st$max_residual <- max(st$max_residual, s$residual)
    if (s$value < b$value) {
      st$best <- list(theta_c = theta_c, value = s$value, delta = s$delta,
                      gradient = s$gradient, residual = s$residual,
                      decrement = s$decrement)
      st$failed_late <- 0L
      st$error <- NULL
    }
    list(objective = s$value, gradient = s$gradient[c_idx], overflow = FALSE)
  }
  list(f = f, state = st)
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
  if (is.null(object[["data"]])) {
    stop("Cannot compute Hessian: no data stored. Refit with keep_data = TRUE.")
  }

  theta <- object$coefficients

  switch(object$model,
    mnl = {
      se_method_mnl <- object$se_method %||% "hessian"
      if (identical(se_method_mnl, "bhhh")) {
        mnl_bhhh_parallel(
          theta = theta,
          X = object[["data"]]$X,
          alt_idx = object[["data"]]$alt_idx,
          choice_idx = object[["data"]]$choice_idx,
          M = object[["data"]]$M,
          weights = object[["data"]]$weights,
          use_asc = object$use_asc,
          include_outside_option = object$include_outside_option
        )
      } else {
        mnl_loglik_hessian_parallel(
          theta = theta,
          X = object[["data"]]$X,
          alt_idx = object[["data"]]$alt_idx,
          choice_idx = object[["data"]]$choice_idx,
          M = object[["data"]]$M,
          weights = object[["data"]]$weights,
          use_asc = object$use_asc,
          include_outside_option = object$include_outside_option
        )
      }
    },
    mxl = {
      # One draw block and one Ti entry per likelihood unit: in generate mode
      # a missing Ti would silently give cross-sectional derivatives.
      gp <- .mxl_gen_params(object$draws_info,
                            N = length(.unit_first(object[["data"]])))
      se_method <- object$se_method %||% "hessian"
      if (se_method == "bhhh") {
        mxl_bhhh_parallel(
          theta = theta,
          X = object[["data"]]$X,
          W = object[["data"]]$W,
          alt_idx = object[["data"]]$alt_idx,
          choice_idx = object[["data"]]$choice_idx,
          M = object[["data"]]$M,
          weights = object[["data"]]$weights,
          eta_draws = gp$eta_draws,
          rc_dist = object$rc_dist,
          rc_correlation = object$rc_correlation,
          rc_mean = object$rc_mean,
          use_asc = object$use_asc,
          include_outside_option = object$include_outside_option,
          gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S,
          Ti = object[["data"]]$Ti
        )
      } else {
        mxl_hessian_parallel(
          theta = theta,
          X = object[["data"]]$X,
          W = object[["data"]]$W,
          alt_idx = object[["data"]]$alt_idx,
          choice_idx = object[["data"]]$choice_idx,
          M = object[["data"]]$M,
          weights = object[["data"]]$weights,
          eta_draws = gp$eta_draws,
          rc_dist = object$rc_dist,
          rc_correlation = object$rc_correlation,
          rc_mean = object$rc_mean,
          use_asc = object$use_asc,
          include_outside_option = object$include_outside_option,
          gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S,
          Ti = object[["data"]]$Ti
        )
      }
    },
    nl = {
      se_method_nl <- object$se_method %||% "hessian"
      if (identical(se_method_nl, "numeric")) {
        nl_loglik_numeric_hessian(
          theta = theta,
          X = object[["data"]]$X,
          alt_idx = object[["data"]]$alt_idx,
          choice_idx = object[["data"]]$choice_idx,
          nest_idx = object[["data"]]$nest_idx,
          M = object[["data"]]$M,
          weights = object[["data"]]$weights,
          use_asc = object$use_asc,
          include_outside_option = object$include_outside_option
        )
      } else if (identical(se_method_nl, "bhhh")) {
        nl_bhhh_parallel(
          theta = theta,
          X = object[["data"]]$X,
          alt_idx = object[["data"]]$alt_idx,
          choice_idx = object[["data"]]$choice_idx,
          nest_idx = object[["data"]]$nest_idx,
          M = object[["data"]]$M,
          weights = object[["data"]]$weights,
          use_asc = object$use_asc,
          include_outside_option = object$include_outside_option
        )
      } else {
        nl_loglik_hessian_parallel(
          theta = theta,
          X = object[["data"]]$X,
          alt_idx = object[["data"]]$alt_idx,
          choice_idx = object[["data"]]$choice_idx,
          nest_idx = object[["data"]]$nest_idx,
          M = object[["data"]]$M,
          weights = object[["data"]]$weights,
          use_asc = object$use_asc,
          include_outside_option = object$include_outside_option
        )
      }
    },
    stop("Unknown model type: '", object$model, "'.")
  )
}

#' Invert an observed information matrix to get vcov
#'
#' Accepts either a negated-Hessian or a BHHH/OPG estimate of the observed
#' information matrix and returns its inverse plus standard errors.
#'
#' The matrix is inverted equilibrated: \eqn{V = D (D H D)^{-1} D}, with
#' \eqn{D} the Jacobi scale of `.jacobi_scale()` (powers of two, so the
#' products with \eqn{D} are exact and only the inversion rounds), by
#' `equilibrated_inverse()`, which scales and factors one working copy and
#' solves for the result, with `solve()`'s checks. An information matrix in
#' natural units can have curvatures spanning many orders of magnitude
#' (columns in different units, a rarely chosen alternative's constant);
#' equilibrated, it inverts as accurately as its well-scaled equivalent
#' instead of losing digits or failing as computationally singular.
#'
#' @param hess Observed information matrix (negated Hessian or BHHH/OPG).
#' @returns List with vcov (matrix, or NULL for a singular matrix) and se
#'   (numeric vector). Errors other than singularity (an allocation that
#'   fails, say) are passed on.
#' @noRd
invert_hessian <- function(hess) {
  p_len <- nrow(hess)
  vcov_mat <- NULL
  se <- rep(NA_real_, p_len)

  singular_flag <- FALSE
  tryCatch({
    vcov_mat <- equilibrated_inverse(hess, .jacobi_scale(hess))
  }, error = function(e) {
    # A singular information gives NA standard errors; any other error (an
    # allocation that fails, say) is the caller's
    if (!grepl("singular", conditionMessage(e), fixed = TRUE)) stop(e)
    singular_flag <<- TRUE
    message("Error inverting information matrix (likely singular): ", e$message)
  })

  if (!singular_flag && !is.null(vcov_mat)) {
    diag_vcov <- diag(vcov_mat)
    neg <- !is.na(diag_vcov) & diag_vcov < 0
    if (any(neg)) {
      message(
        "Information matrix is not positive definite; ",
        sum(neg), " variance(s) negative. Standard errors set to NA for ",
        "those parameters (optimum may not have been reached)."
      )
      diag_vcov[neg] <- NA_real_
    }
    se <- sqrt(diag_vcov)
  } else {
    message("Standard errors set to NA due to information matrix inversion failure.")
  }

  list(vcov = vcov_mat, se = se)
}

#' Combine bread and meat into a robust sandwich variance
#'
#' Forms \code{V = A^{-1} B A^{-1}} from the (weighted) negated Hessian
#' \code{A} (bread) and the (weight-squared) outer-product-of-gradients
#' \code{B} (meat), with the same singular / not-positive-definite guards as
#' \code{invert_hessian()}. Equilibrated by the Jacobi scale \eqn{D} of the
#' bread: \eqn{V = D [(DAD)^{-1} (DBD) (DAD)^{-1}] D}, the inner matrix
#' symmetrized before the exact products with \eqn{D}, by
#' `equilibrated_sandwich()` (two solves with the bread's LU factors, no
#' explicit inverse).
#'
#' @param A Bread matrix (observed information; weighted negated Hessian).
#' @param B Meat matrix (weighted outer product of per-individual scores).
#' @returns List with \code{vcov} (matrix, or NULL for a singular bread) and
#'   \code{se} (numeric). Errors other than singularity are passed on.
#' @noRd
.sandwich_combine <- function(A, B) {
  p_len <- nrow(A)
  vcov_mat <- NULL
  se <- rep(NA_real_, p_len)

  singular_flag <- FALSE
  tryCatch({
    vcov_mat <- equilibrated_sandwich(A, B, .jacobi_scale(A))
  }, error = function(e) {
    # As in invert_hessian(): only a singular bread is reported here
    if (!grepl("singular", conditionMessage(e), fixed = TRUE)) stop(e)
    singular_flag <<- TRUE
    message("Error inverting information (bread) matrix (likely singular): ",
            e$message)
  })

  if (!singular_flag && !is.null(vcov_mat)) {
    diag_vcov <- diag(vcov_mat)
    neg <- !is.na(diag_vcov) & diag_vcov < 0
    if (any(neg)) {
      message(
        "Sandwich covariance is not positive definite; ",
        sum(neg), " variance(s) negative. Standard errors set to NA for ",
        "those parameters."
      )
      diag_vcov[neg] <- NA_real_
    }
    se <- sqrt(diag_vcov)
  } else {
    message("Standard errors set to NA due to sandwich inversion failure.")
  }

  list(vcov = vcov_mat, se = se)
}

#' The vcov(type = ) that recomputes a fit's se_method
#'
#' \code{"robust"} for the sandwich; \code{"hessian"} for a nested logit's
#' \code{"numeric"}, whose bread \code{.compute_bread()} takes numerically for
#' such fits; the method's own name otherwise.
#'
#' @param se_method A fit's \code{se_method}.
#' @returns The \code{type} of \code{vcov.choicer_fit()}.
#' @noRd
.vcov_type <- function(se_method) {
  switch(se_method, sandwich = "robust", numeric = "hessian", se_method)
}

#' No variance for a fit whose variance step failed
#'
#' The error handler of the `run_*()` wrappers' variance step, which runs
#' after the optimizer: an error there (a kernel's memory check, an
#' allocation that fails) would otherwise
#' stop the call and discard the estimates. Warns, quoting the error and
#' naming the post-hoc route, and returns no variance: `vcov` NULL, `se` NA,
#' and the error's message, which the wrapper keeps as the fit's
#' `variance_error` (on such fits only; `$vcov` cannot partially match the
#' name), so that `ensure_vcov()` does not run the failed step again from
#' `summary()` and the other methods that call it, while `vcov(fit, type = )`
#' runs it on request. Nothing the size of the variance is allocated here:
#' where the p x p result itself did not fit, a matrix of NA would not fit
#' either. The mixed logit's memory check names in its own message the
#' threads that fit or how to lift the limit, so no advice on memory is added.
#'
#' @param e The error.
#' @param se_method The fit's `se_method`.
#' @param param_names The parameter names.
#' @param keep_data Whether the fit keeps its data.
#' @returns List with `vcov` (NULL), `se` (p NA, named) and `error` (the
#'   error's message).
#' @noRd
.vcov_failed <- function(e, se_method, param_names, keep_data) {
  retry <- sprintf("vcov(fit, type = \"%s\")", .vcov_type(se_method))
  err <- trimws(conditionMessage(e))
  if (!nzchar(err)) err <- paste("an error of class", class(e)[1L])
  if (!grepl("[.!?]$", err)) err <- paste0(err, ".")
  advice <- if (keep_data) {
    paste0("To compute it from the fit, call ", retry)
  } else {
    paste0("The fit keeps no data (keep_data = FALSE): refit with ",
           "keep_data = TRUE to compute it with ", retry)
  }
  warning(warningCondition(paste0(
    "The variance could not be computed (se_method = \"", se_method, "\"), ",
    "so the estimates are returned without it, with NA standard errors: ",
    err, " ", advice, "."),
    class = "choicer_vcov_failed"))
  list(vcov = NULL,
       se = stats::setNames(rep(NA_real_, length(param_names)), param_names),
       error = err)
}

#' Jacobi equilibration scale of an information matrix
#'
#' \eqn{D_k = 2^{\mathrm{round}(-\log_2 |A_{kk}| / 2)}} for a finite, nonzero
#' diagonal entry and 1 otherwise, within \eqn{[2^{-500}, 2^{500}]}: powers of
#' two, so \eqn{D A D} and the products with \eqn{D} are exact, and, unless
#' clamped, the scaled matrix's diagonal is within a factor of two of one in
#' magnitude. The magnitude, so that a matrix that is not positive definite
#' away from the optimum (a negative diagonal entry) is equilibrated too and
#' its negative variances are reported rather than a singular inversion. For
#' a positive semi-definite matrix \eqn{|A_{ij}| \le \sqrt{A_{ii} A_{jj}}},
#' so no entry of \eqn{D A D} exceeds about 2 in magnitude (about
#' \eqn{2^{24}} where a scale is clamped); an indefinite
#' one whose \eqn{|A_{ij}| / \sqrt{|A_{ii} A_{jj}|}} exceeds about 1e308
#' overflows, and its variances then come out NaN without a message.
#'
#' @param A Square matrix.
#' @returns Numeric vector of scales.
#' @noRd
.jacobi_scale <- function(A) {
  d <- diag(A)
  ok <- is.finite(d) & d != 0
  D <- rep(1, length(d))
  D[ok] <- 2^round(-0.5 * log2(abs(d[ok])))
  pmin(pmax(D, 2^-500), 2^500)
}

#' Robust (Huber-White) sandwich vcov for a fitted logit model
#'
#' Computes \code{V = A^{-1} B A^{-1}} from stored (natural-scale) data, where
#' \code{A = sum_i w_i (-H_i)} is the weighted negated Hessian and
#' \code{B = sum_i w_i^2 s_i s_i'} is the weight-squared outer product of
#' per-individual scores — the \code{"robust"} case of the shared
#' \code{.assemble_score_vcov()} path (each model's BHHH kernel with the
#' squared weights, without the score matrix). Valid under choice-based / WESML
#' weighting, where the plain inverse-Hessian is not. Supported for MNL, MXL
#' and NL fits.
#'
#' @param object A fitted \code{choicer_fit} object (MNL / MXL / NL) with
#'   \code{keep_data = TRUE}.
#' @returns List with \code{vcov} and \code{se}.
#' @noRd
compute_sandwich_vcov <- function(object) {
  .assemble_score_vcov(object, type = "robust")
}

#' Per-unit score matrix from stored data
#'
#' Dispatches to the internal \code{*_scores_parallel} C++ kernels and returns
#' the \code{U x p} matrix whose row \code{u} is the weight-free score of
#' likelihood unit \code{u} evaluated at the fitted coefficients. Units are
#' choice situations, or decision makers for a panel mixed logit
#' (\code{person_col}). For MXL, Halton draws are regenerated
#' deterministically from \code{draws_info} (mirroring
#' \code{compute_hessian()}).
#'
#' No variance route calls it: the kernels form the meats without the matrix
#' (\code{.mxl_score_meat()}, \code{.logit_meat()}).
#'
#' @param object A fitted \code{choicer_fit} object (MNL / MXL / NL) with
#'   \code{keep_data = TRUE}.
#' @returns Numeric matrix, \code{U x n_params} (\code{U = nobs} unless the
#'   fit is a panel mixed logit).
#' @noRd
compute_scores <- function(object) {
  if (is.null(object[["data"]])) {
    stop("Cannot compute scores: no data stored. Refit with keep_data = TRUE.")
  }
  if (!object$model %in% c("mnl", "mxl", "nl")) {
    stop("Score-based standard errors are implemented for multinomial (MNL), ",
         "mixed (MXL) and nested (NL) logit only.")
  }

  theta <- object$coefficients
  d <- object[["data"]]

  switch(object$model,
    mnl = mnl_scores_parallel(
      theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
      M = d$M, use_asc = object$use_asc,
      include_outside_option = object$include_outside_option
    ),
    nl = nl_scores_parallel(
      theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
      nest_idx = d$nest_idx, M = d$M, use_asc = object$use_asc,
      include_outside_option = object$include_outside_option
    ),
    mxl = {
      gp <- .mxl_gen_params(object$draws_info, N = length(.unit_first(d)))
      mxl_scores_parallel(
        theta = theta, X = d$X, W = d$W,
        alt_idx = d$alt_idx, choice_idx = d$choice_idx,
        M = d$M, eta_draws = gp$eta_draws,
        rc_dist = object$rc_dist, rc_correlation = object$rc_correlation,
        rc_mean = object$rc_mean, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option,
        gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S,
        Ti = d$Ti
      )
    }
  )
}

#' Bread matrix: weighted negated Hessian from stored data
#'
#' Always the analytical Hessian (numeric fallback for NL fits estimated with
#' \code{se_method = "numeric"}), regardless of \code{object$se_method} —
#' unlike \code{compute_hessian()}, which returns the BHHH matrix for
#' \code{se_method = "bhhh"} fits.
#'
#' @param gp For a mixed logit, the draw arguments of
#'   \code{.mxl_gen_params()} when the caller has built them already (a
#'   variance's meat and bread use the same draws); \code{NULL} builds them.
#' @noRd
.compute_bread <- function(object, gp = NULL) {
  theta <- object$coefficients
  d <- object[["data"]]
  w <- d$weights

  switch(object$model,
    mnl = mnl_loglik_hessian_parallel(
      theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
      M = d$M, weights = w, use_asc = object$use_asc,
      include_outside_option = object$include_outside_option
    ),
    nl = if (identical(object$se_method %||% "hessian", "numeric")) {
      nl_loglik_numeric_hessian(
        theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
        nest_idx = d$nest_idx, M = d$M, weights = w, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option
      )
    } else {
      nl_loglik_hessian_parallel(
        theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
        nest_idx = d$nest_idx, M = d$M, weights = w, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option
      )
    },
    mxl = {
      if (is.null(gp)) {
        gp <- .mxl_gen_params(object$draws_info, N = length(.unit_first(d)))
      }
      mxl_hessian_parallel(
        theta = theta, X = d$X, W = d$W,
        alt_idx = d$alt_idx, choice_idx = d$choice_idx,
        M = d$M, weights = w, eta_draws = gp$eta_draws,
        rc_dist = object$rc_dist, rc_correlation = object$rc_correlation,
        rc_mean = object$rc_mean, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option,
        gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S,
        Ti = d$Ti
      )
    }
  )
}

#' Meat matrix from a per-unit score matrix
#'
#' The meats of the score-based variance estimators, from the score matrix:
#' the R reference of the kernels' routes (\code{.mxl_score_meat()},
#' \code{.logit_meat()}), which form them without the matrix; no variance
#' route calls it.
#' \describe{
#'   \item{bhhh}{\code{crossprod(sqrt(w) * S)} — the ordinary weighted
#'     BHHH/OPG information \eqn{\sum_i w_i s_i s_i'}.}
#'   \item{robust}{\code{crossprod(w * S)} — the Huber-White / WESML sandwich
#'     meat \eqn{\sum_i w_i^2 s_i s_i'} (the \eqn{w^2} special case).}
#'   \item{cluster}{\code{crossprod(rowsum(w * S, cluster))} — within-cluster
#'     sums of weighted scores, then their outer product
#'     \eqn{\sum_g (\sum_{i \in g} w_i s_i)(\sum_{i \in g} w_i s_i)'}.
#'     Invariant to the ordering of cluster labels (\code{rowsum} groups by
#'     value). No small-sample (\eqn{G/(G-1)}) correction is applied.}
#' }
#'
#' @param S \code{N x p} score matrix (rows are weight-free scores of the
#'   likelihood units: choice situations, or decision makers for a panel
#'   mixed logit).
#' @param w Length-\code{N} weight vector, one weight per score row.
#' @param type One of \code{"bhhh"}, \code{"robust"}, \code{"cluster"}.
#' @param cluster Length-\code{N} cluster labels (required for
#'   \code{type = "cluster"}).
#' @returns \code{p x p} meat matrix.
#' @noRd
.score_meat <- function(S, w, type, cluster = NULL) {
  if (length(w) != nrow(S)) {
    stop("Internal error: ", length(w), " weights for ", nrow(S),
         " score rows (weights must be per likelihood unit).", call. = FALSE)
  }
  switch(type,
    bhhh = crossprod(sqrt(w) * S),
    robust = crossprod(w * S),
    cluster = {
      .check_cluster_labels(cluster, nrow(S))
      crossprod(rowsum(w * S, group = as.character(cluster)))
    },
    stop("Unknown meat type: '", type, "'.")
  )
}

#' Check cluster labels, one per likelihood unit
#'
#' The checks of the cluster meat, shared by \code{.score_meat()} and the
#' kernels' routes, which the \code{run_*()} wrappers also run on their
#' fit-time labels before the optimizer.
#'
#' @param cluster Cluster labels, one per likelihood unit (or \code{NULL}).
#' @param n_units Number of likelihood units.
#' @returns \code{cluster}, invisibly.
#' @noRd
.check_cluster_labels <- function(cluster, n_units) {
  if (is.null(cluster)) {
    stop("Clustered standard errors need `cluster`: a vector with one ",
         "cluster label per choice situation (or fit with `cluster_col=`).",
         call. = FALSE)
  }
  if (length(cluster) != n_units) {
    stop("`cluster` has length ", length(cluster), " but the model has ",
         n_units, " choice situations. Supply one cluster label per ",
         "choice situation, aligned with the prepared data ",
         "(situations sorted by id).", call. = FALSE)
  }
  if (anyNA(cluster)) {
    stop("`cluster` contains missing values.", call. = FALSE)
  }
  invisible(cluster)
}

#' Cluster labels coded 1..G in order of first appearance
#'
#' The codes the cluster meat kernels take:
#' \code{match(as.character(x), unique(as.character(x)))}, the grouping of
#' \code{rowsum()}, numbered as the labels first appear. Where the labels'
#' own values group them as their strings do, the codes come from the values,
#' without the strings (\code{as.character()} of 10^7 integer labels took
#' 9.5 s and 1.1 GB): a factor's integer codes, unless two of its levels are
#' the same string; plain integer, logical and character vectors as they are.
#' Every other vector (doubles, whose strings merge values that print alike,
#' dates, integer64 and other classes) goes through \code{as.character()}.
#' The partition and its numbering are the same either way.
#'
#' @param x Cluster labels, checked by \code{.check_cluster_labels()}.
#' @returns Integer codes, one per label.
#' @noRd
.cluster_codes <- function(x) {
  if (is.factor(x)) {
    x <- if (anyDuplicated(levels(x))) as.character(x) else as.integer(x)
  } else if (is.object(x) || !(is.integer(x) || is.logical(x) || is.character(x))) {
    x <- as.character(x)
  }
  match(x, unique(x))
}

#' Score-based meat of a mixed logit fit, from the C++ kernels
#'
#' The meats of \code{.score_meat()} without the \code{U x p} score matrix
#' (2 million decision makers by 2,620 parameters is 42 GB): \code{"bhhh"} is
#' \code{mxl_bhhh_parallel()} with the weights, \code{"robust"} the same
#' kernel with the squared weights, and \code{"cluster"}
#' \code{mxl_cluster_meat_parallel()} over the per-unit labels, checked as
#' \code{.score_meat()} checks them and grouped as \code{rowsum()} groups
#' them (\code{.cluster_codes()}).
#'
#' @param object A fitted \code{choicer_mxl} with stored data.
#' @param gp Draw arguments from \code{.mxl_gen_params()}.
#' @param type One of \code{"bhhh"}, \code{"robust"}, \code{"cluster"}.
#' @param cluster Cluster labels, one per likelihood unit
#'   (\code{type = "cluster"}).
#' @returns \code{p x p} meat matrix.
#' @noRd
.mxl_score_meat <- function(object, gp, type, cluster = NULL) {
  if (!type %in% c("bhhh", "robust", "cluster")) {
    stop("Unknown meat type: '", type, "'.")
  }
  d <- object[["data"]]
  theta <- object$coefficients
  if (identical(type, "cluster")) {
    .check_cluster_labels(cluster, length(.unit_first(d)))
    return(mxl_cluster_meat_parallel(
      theta = theta, X = d$X, W = d$W,
      alt_idx = d$alt_idx, choice_idx = d$choice_idx,
      M = d$M, weights = d$weights, cluster = .cluster_codes(cluster),
      eta_draws = gp$eta_draws,
      rc_dist = object$rc_dist, rc_correlation = object$rc_correlation,
      rc_mean = object$rc_mean, use_asc = object$use_asc,
      include_outside_option = object$include_outside_option,
      gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S,
      Ti = d$Ti
    ))
  }
  mxl_bhhh_parallel(
    theta = theta, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx,
    M = d$M, weights = if (identical(type, "robust")) d$weights^2 else d$weights,
    eta_draws = gp$eta_draws,
    rc_dist = object$rc_dist, rc_correlation = object$rc_correlation,
    rc_mean = object$rc_mean, use_asc = object$use_asc,
    include_outside_option = object$include_outside_option,
    gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S,
    Ti = d$Ti
  )
}

#' Score-based meat of a multinomial or nested logit fit, from the C++ kernels
#'
#' The meats of \code{.score_meat()} without the \code{N x p} score matrix, as
#' \code{.mxl_score_meat()} forms the mixed logit's: \code{"bhhh"} is the
#' BHHH kernel with the weights, \code{"robust"} the same kernel with the
#' squared weights, and \code{"cluster"} the cluster meat kernel over the
#' labels, checked as \code{.score_meat()} checks them and grouped as
#' \code{rowsum()} groups them (\code{.cluster_codes()}).
#'
#' @param object A fitted \code{choicer_mnl} or \code{choicer_nl} with stored
#'   data.
#' @param type One of \code{"bhhh"}, \code{"robust"}, \code{"cluster"}.
#' @param cluster Cluster labels, one per choice situation
#'   (\code{type = "cluster"}).
#' @returns \code{p x p} meat matrix.
#' @noRd
.logit_meat <- function(object, type, cluster = NULL) {
  if (!type %in% c("bhhh", "robust", "cluster")) {
    stop("Unknown meat type: '", type, "'.")
  }
  d <- object[["data"]]
  theta <- object$coefficients
  if (identical(type, "cluster")) {
    .check_cluster_labels(cluster, length(d$M))
    codes <- .cluster_codes(cluster)
    return(switch(object$model,
      mnl = mnl_cluster_meat_parallel(
        theta = theta, X = d$X, alt_idx = d$alt_idx,
        choice_idx = d$choice_idx, M = d$M, weights = d$weights,
        cluster = codes, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option
      ),
      nl = nl_cluster_meat_parallel(
        theta = theta, X = d$X, alt_idx = d$alt_idx,
        choice_idx = d$choice_idx, nest_idx = d$nest_idx, M = d$M,
        weights = d$weights, cluster = codes, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option
      )
    ))
  }
  w <- if (identical(type, "robust")) d$weights^2 else d$weights
  switch(object$model,
    mnl = mnl_bhhh_parallel(
      theta = theta, X = d$X, alt_idx = d$alt_idx,
      choice_idx = d$choice_idx, M = d$M, weights = w,
      use_asc = object$use_asc,
      include_outside_option = object$include_outside_option
    ),
    nl = nl_bhhh_parallel(
      theta = theta, X = d$X, alt_idx = d$alt_idx,
      choice_idx = d$choice_idx, nest_idx = d$nest_idx, M = d$M,
      weights = w, use_asc = object$use_asc,
      include_outside_option = object$include_outside_option
    )
  )
}

#' Resolve a user-supplied post-hoc cluster vector to prepared order
#'
#' Guards the alignment hazards of \code{vcov(fit, type = "cluster",
#' cluster = )}, where the supplied labels must line up with the prepared
#' (id-sorted) choice situations:
#' \itemize{
#'   \item A row-level vector (length \code{sum(M)}) is rejected outright — the
#'     most common mistake is passing raw per-alternative labels.
#'   \item A \emph{named} vector is realigned to the prepared order by matching
#'     its names against the stored choice-situation ids, so it is safe
#'     regardless of the order the user built it in. This is the recommended
#'     post-hoc form.
#'   \item An \emph{unnamed} vector of length \code{N} is taken to be in
#'     prepared order, with a warning flagging that assumption.
#' }
#' The fit-time \code{cluster_col=} path never reaches here: its labels are
#' collapsed and stored already aligned.
#'
#' @param object A fitted \code{choicer_fit} with stored data.
#' @param cluster User-supplied cluster labels (named or unnamed).
#' @returns A length-\code{N} vector aligned with the prepared choice
#'   situations.
#' @noRd
.resolve_cluster <- function(object, cluster) {
  d <- object[["data"]]
  N <- length(d$M)
  n_rows <- nrow(d$X)

  # Most common mistake: per-alternative (row-level) labels.
  if (n_rows != N && length(cluster) == n_rows) {
    stop("`cluster` has length ", n_rows, ", the number of stacked ",
         "alternative rows, but cluster-robust standard errors need one ",
         "label per choice situation (N = ", N, "). Collapse to one label ",
         "per choice situation, or set `cluster_col=` at fit time.",
         call. = FALSE)
  }

  nm <- names(cluster)
  ids <- d$situation_ids

  if (!is.null(nm)) {
    if (is.null(ids)) {
      stop("Cannot realign a named `cluster`: this fit does not carry ",
           "choice-situation ids (it predates id storage). Supply `cluster` ",
           "in the prepared, id-sorted order (unnamed), or set `cluster_col=` ",
           "at fit time.", call. = FALSE)
    }
    if (anyDuplicated(nm)) {
      stop("`cluster` has duplicate names; names must be the unique ",
           "choice-situation ids.", call. = FALSE)
    }
    idx <- match(as.character(ids), nm)
    if (anyNA(idx)) {
      miss <- as.character(ids)[is.na(idx)]
      stop("`cluster` names do not cover every choice situation (",
           length(miss), " of ", N, " missing, e.g. ",
           paste(utils::head(miss, 3L), collapse = ", "),
           "). Name `cluster` by choice-situation id.", call. = FALSE)
    }
    return(unname(cluster[idx]))
  }

  # Unnamed: must already be in prepared (id-sorted) order.
  if (length(cluster) != N) {
    stop("`cluster` has length ", length(cluster), " but the model has ", N,
         " choice situations. Supply one label per choice situation, name it ",
         "by choice-situation id, or set `cluster_col=` at fit time.",
         call. = FALSE)
  }
  # A panel fit orders situations by decision maker first; unless that is
  # also id order, unnamed labels (typically built in id order) would be
  # silently misaligned, so require names.
  if (!is.null(d$Ti) && !is.null(ids) && !.ids_sorted(ids)) {
    stop("Unnamed `cluster` is ambiguous for this panel fit: its choice ",
         "situations are ordered by decision maker (`person_col`), not by id. ",
         "Name `cluster` by choice-situation id (see `fit$data$situation_ids`) ",
         "or set `cluster_col=` at fit time.", call. = FALSE)
  }
  warning("Unnamed `cluster` is assumed to be in the prepared (id-sorted) ",
          "order. Name it by choice-situation id (see `fit$data$situation_ids`) ",
          "or set `cluster_col=` at fit time to guarantee alignment.",
          call. = FALSE)
  cluster
}

#' Score-based variance assembly (bhhh / robust / cluster) from stored data
#'
#' The single post-hoc entry point behind \code{vcov(fit, type = )},
#' \code{ensure_vcov()} (for \code{se_method = "cluster"} fits) and
#' \code{compute_sandwich_vcov()}, and every model's variance at fit time.
#' Operates entirely in natural space: the stored data and coefficients are
#' natural-scale even when the fit used \code{scale_vars}, so no
#' back-transform is needed; every matrix is inverted equilibrated
#' (\code{invert_hessian()}, \code{.sandwich_combine()}).
#'
#' @param object A fitted \code{choicer_fit} (MNL / MXL / NL) with
#'   \code{keep_data = TRUE}.
#' @param type One of \code{"hessian"}, \code{"bhhh"}, \code{"robust"},
#'   \code{"cluster"}.
#' @param cluster Cluster labels for \code{type = "cluster"}; defaults to the
#'   labels stored at fit time via \code{cluster_col}.
#' @param gp For a mixed logit, the draw arguments of \code{.mxl_gen_params()}
#'   when the caller holds them already (the fit itself, so that its variance
#'   uses its own draws without building a second store-mode cube);
#'   \code{NULL} builds them from \code{draws_info}.
#' @param what The cluster labels' name in messages (the fit itself says
#'   \code{cluster_col}).
#' @param bread,bhhh The weighted negated Hessian, and the BHHH matrix with
#'   the unsquared weights, at the fit's coefficients with its own data and
#'   draws, when the caller holds them already (a mixed logit polish's last
#'   curvature): used, unmodified, in place of \code{.compute_bread()} (the
#'   \code{"hessian"} type and the sandwiches' bread) and of the
#'   \code{"bhhh"} type's matrix. \code{NULL} computes them.
#' @returns List with \code{vcov} and \code{se}.
#' @noRd
.assemble_score_vcov <- function(object, type, cluster = NULL, gp = NULL,
                                 what = "`cluster`", bread = NULL,
                                 bhhh = NULL) {
  if (is.null(object[["data"]])) {
    stop("Cannot compute standard errors: no data stored. ",
         "Refit with keep_data = TRUE.")
  }
  if (!object$model %in% c("mnl", "mxl", "nl")) {
    stop("Score-based standard errors are implemented for multinomial (MNL), ",
         "mixed (MXL) and nested (NL) logit only.")
  }

  if (identical(type, "hessian")) {
    return(invert_hessian(bread %||% .compute_bread(object, gp)))
  }
  if (identical(type, "bhhh") && !is.null(bhhh)) return(invert_hessian(bhhh))

  # Scores, weights and cluster labels are per likelihood unit: choice
  # situations, or decision makers for a panel mixed logit.
  d <- object[["data"]]

  if (identical(object$model, "mxl")) {
    # The mixed logit's meats come from the C++ kernels, which sum the
    # scores without forming the U x p score matrix; the meat and the bread
    # share one set of draws, built once.
    if (identical(type, "cluster")) {
      cluster <- .unit_clusters(object, cluster, what)
      # before the draws, which in store mode can take gigabytes
      .check_cluster_labels(cluster, length(.unit_first(d)))
    }
    if (is.null(gp)) {
      gp <- .mxl_gen_params(object$draws_info, N = length(.unit_first(d)))
    }
    B <- .mxl_score_meat(object, gp, type, cluster)
    if (identical(type, "bhhh")) return(invert_hessian(B))
    return(.sandwich_combine(bread %||% .compute_bread(object, gp), B))
  }

  # The multinomial and nested logits' meats come from their kernels too,
  # without the N x p score matrix; the labels are checked before any kernel.
  if (identical(type, "cluster")) {
    cluster <- .unit_clusters(object, cluster, what)
    .check_cluster_labels(cluster, length(d$M))
  }
  B <- .logit_meat(object, type, cluster)
  if (identical(type, "bhhh")) return(invert_hessian(B))
  .sandwich_combine(bread %||% .compute_bread(object), B)
}

#' Cluster labels of a fit's likelihood units
#'
#' The labels \code{vcov(type = "cluster")} uses: the fit-time
#' \code{cluster_col} labels when \code{cluster} is \code{NULL}, else the
#' user's labels realigned to the prepared order (\code{.resolve_cluster()});
#' then one per likelihood unit (\code{.to_units()}: a panel fit's clusters
#' must nest decision makers).
#'
#' @param object A fitted \code{choicer_fit} with stored data.
#' @param cluster User-supplied labels, or \code{NULL}.
#' @param what The labels' name in messages.
#' @returns Cluster labels, one per likelihood unit.
#' @noRd
.unit_clusters <- function(object, cluster, what = "`cluster`") {
  d <- object[["data"]]
  if (is.null(cluster)) {
    # No explicit labels: use the fit-time cluster_col vector, already
    # collapsed and aligned to the prepared order.
    cluster <- d$cluster
    if (is.null(cluster)) {
      stop("Cluster-robust standard errors need cluster labels: pass ",
           "`cluster=` (named by choice-situation id) or fit with ",
           "`cluster_col=`.", call. = FALSE)
    }
  } else {
    # User-supplied labels: guard and realign to the prepared order.
    cluster <- .resolve_cluster(object, cluster)
  }
  # Per situation so far; panel fits need clusters that nest decision makers.
  .to_units(cluster, d, what)
}
