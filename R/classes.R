
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
#' @param scale_vars Character. Pre-estimation scaling applied to the design
#'   matrix: \code{"none"} (default), \code{"sd"}, \code{"mad"}, or \code{"iqr"}.
#' @param sX Named numeric vector of column scales used to standardize X during
#'   optimization. Defaults to a vector of 1s when scale_vars = 'none'.
#' @param se_method Character. Method used for standard errors: \code{"hessian"}
#'   (analytical Hessian, default), \code{"bhhh"} (outer product of gradients),
#'   or \code{"sandwich"} (robust Huber--White / WESML variance).
#' @param choice_sampling Optional list recording choice-based-sampling
#'   provenance (scheme, population/sample shares, meat type), or NULL.
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
                            choice_sampling = NULL) {
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
      choice_sampling = choice_sampling
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
#' @param scale_vars Character. Pre-estimation scaling applied to the design
#'   matrices: \code{"none"} (default) or \code{"sd"}.
#' @param sX Named numeric vector of column scales used to standardize \code{X}
#'   during optimization. Defaults to a vector of 1s (no scaling applied) when
#'   \code{scale_vars = "none"}; equals \code{apply(X, 2, sd)} when
#'   \code{scale_vars = "sd"}.
#' @param sW Named numeric vector of column scales used to standardize \code{W}
#'   during optimization. Defaults to a vector of 1s when
#'   \code{scale_vars = "none"}; equals \code{apply(W, 2, sd)} for normal
#'   random-coefficient columns when \code{scale_vars = "sd"}, with entries for
#'   log-normal columns (\code{rc_dist[k] == 1}) carved out to 1.
#' @param choice_sampling Optional list recording choice-based-sampling
#'   provenance (scheme, population/sample shares, meat type), or NULL.
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
                            choice_sampling = NULL) {
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
      choice_sampling = choice_sampling
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
                           choice_sampling = NULL) {
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
      choice_sampling = choice_sampling
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

  if (is.null(object[["data"]])) {
    message("Cannot compute standard errors: no data stored. ",
            "Refit with keep_data = TRUE.")
    return(object)
  }

  se_method <- object$se_method %||% "hessian"
  result <- if (identical(se_method, "sandwich")) {
    compute_sandwich_vcov(object)
  } else {
    invert_hessian(compute_hessian(object))
  }

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
      gp <- .mxl_gen_params(object$draws_info)
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
          gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S
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
          gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S
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
#' @param hess Observed information matrix (negated Hessian or BHHH/OPG).
#' @returns List with vcov (matrix or NULL) and se (numeric vector).
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
#' \code{invert_hessian()}.
#'
#' @param A Bread matrix (observed information; weighted negated Hessian).
#' @param B Meat matrix (weighted outer product of per-individual scores).
#' @returns List with \code{vcov} (matrix or NULL) and \code{se} (numeric).
#' @noRd
.sandwich_combine <- function(A, B) {
  p_len <- nrow(A)
  vcov_mat <- NULL
  se <- rep(NA_real_, p_len)

  singular_flag <- FALSE
  tryCatch({
    Ainv <- solve(A)
    vcov_mat <- Ainv %*% B %*% Ainv
  }, error = function(e) {
    singular_flag <<- TRUE
    message("Error inverting information (bread) matrix (likely singular): ",
            e$message)
  })

  if (!singular_flag && !is.null(vcov_mat)) {
    vcov_mat <- (vcov_mat + t(vcov_mat)) / 2   # symmetrize away FP asymmetry
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

#' Robust (Huber-White) sandwich vcov for a fitted logit model
#'
#' Computes \code{V = A^{-1} B A^{-1}} from stored (natural-scale) data, where
#' \code{A = sum_i w_i (-H_i)} is the weighted negated Hessian and
#' \code{B = sum_i w_i^2 s_i s_i'} is the weight-squared outer product of
#' per-individual scores. \code{B} is obtained with zero extra C++ by calling
#' the BHHH routine with squared weights (its per-individual score is
#' weight-free). Valid under choice-based / WESML weighting, where the plain
#' inverse-Hessian is not. Supported for MNL, MXL and NL fits.
#'
#' @param object A fitted \code{choicer_fit} object (MNL / MXL / NL) with
#'   \code{keep_data = TRUE}.
#' @returns List with \code{vcov} and \code{se}.
#' @noRd
compute_sandwich_vcov <- function(object) {
  if (is.null(object[["data"]])) {
    stop("Cannot compute sandwich vcov: no data stored. ",
         "Refit with keep_data = TRUE.")
  }
  if (!object$model %in% c("mnl", "mxl", "nl")) {
    stop("Sandwich standard errors are implemented for multinomial (MNL), ",
         "mixed (MXL) and nested (NL) logit only.")
  }

  # Recomputes A (bread) and B (meat) from the stored NATURAL-scale data and
  # natural-scale coefficients, so it needs no back-transform: by design it
  # operates entirely in natural space (unlike the eager run_mxlogit() path,
  # which works in scaled space and then back-transforms the result).
  theta <- object$coefficients
  d <- object[["data"]]
  w <- d$weights

  res <- switch(object$model,
    mnl = {
      A <- mnl_loglik_hessian_parallel(
        theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
        M = d$M, weights = w, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option
      )
      B <- mnl_bhhh_parallel(
        theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
        M = d$M, weights = w^2, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option
      )
      list(A = A, B = B)
    },
    nl = {
      A <- nl_loglik_hessian_parallel(
        theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
        nest_idx = d$nest_idx, M = d$M, weights = w, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option
      )
      B <- nl_bhhh_parallel(
        theta = theta, X = d$X, alt_idx = d$alt_idx, choice_idx = d$choice_idx,
        nest_idx = d$nest_idx, M = d$M, weights = w^2, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option
      )
      list(A = A, B = B)
    },
    mxl = {
      gp_sw <- .mxl_gen_params(object$draws_info)
      A <- mxl_hessian_parallel(
        theta = theta, X = d$X, W = d$W,
        alt_idx = d$alt_idx, choice_idx = d$choice_idx,
        M = d$M, weights = w, eta_draws = gp_sw$eta_draws,
        rc_dist = object$rc_dist, rc_correlation = object$rc_correlation,
        rc_mean = object$rc_mean, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option,
        gen_seed = gp_sw$gen_seed, gen_scramble = gp_sw$gen_scramble, gen_S = gp_sw$gen_S
      )
      B <- mxl_bhhh_parallel(
        theta = theta, X = d$X, W = d$W,
        alt_idx = d$alt_idx, choice_idx = d$choice_idx,
        M = d$M, weights = w^2, eta_draws = gp_sw$eta_draws,
        rc_dist = object$rc_dist, rc_correlation = object$rc_correlation,
        rc_mean = object$rc_mean, use_asc = object$use_asc,
        include_outside_option = object$include_outside_option,
        gen_seed = gp_sw$gen_seed, gen_scramble = gp_sw$gen_scramble, gen_S = gp_sw$gen_S
      )
      list(A = A, B = B)
    }
  )
  .sandwich_combine(res$A, res$B)
}
