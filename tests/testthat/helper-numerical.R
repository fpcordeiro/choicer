# Helper functions for numerical validation of gradients and Hessians

#' Validate gradient accuracy against numerical derivatives
#'
#' @param objective_fn Function that returns a scalar objective value
#' @param gradient_fn Function that returns a gradient vector
#' @param theta Parameter vector at which to evaluate
#' @param tol Tolerance for relative difference
#' @param method Method for numDeriv::grad (default: Richardson)
#' @return List with comparison results
validate_gradient <- function(
    objective_fn,
    gradient_fn,
    theta,
    tol = 1e-5,
    method = "Richardson"
) {
  analytic <- gradient_fn(theta)
  numeric <- numDeriv::grad(objective_fn, theta, method = method)

  max_abs_diff <- max(abs(analytic - numeric))

  # Use relative difference, scaling by gradient magnitude or 1
  scale <- pmax(abs(numeric), 1)
  max_rel_diff <- max(abs(analytic - numeric) / scale)

  list(
    analytic = analytic,
    numeric = numeric,
    max_abs_diff = max_abs_diff,
    max_rel_diff = max_rel_diff,
    passes = max_rel_diff < tol
  )
}

#' Validate Hessian accuracy against numerical derivatives
#'
#' @param objective_fn Function that returns a scalar objective value
#' @param hessian_fn Function that returns a Hessian matrix
#' @param theta Parameter vector at which to evaluate
#' @param tol Tolerance for absolute difference
#' @param method Method for numDeriv::hessian (default: Richardson)
#' @return List with comparison results
validate_hessian <- function(
    objective_fn,
    hessian_fn,
    theta,
    tol = 1e-4,
    method = "Richardson"
) {
  analytic <- hessian_fn(theta)
  numeric <- numDeriv::hessian(objective_fn, theta, method = method)

  max_abs_diff <- max(abs(analytic - numeric))

  # Check symmetry of analytical Hessian
  is_symmetric <- max(abs(analytic - t(analytic))) < 1e-12

  list(
    analytic = analytic,
    numeric = numeric,
    max_abs_diff = max_abs_diff,
    is_symmetric = is_symmetric,
    passes = max_abs_diff < tol
  )
}

#' Create wrapper functions for likelihood evaluation
#'
#' @param inputs List from prepare_*_data()
#' @param model_type One of "mnl", "mxl", "nl"
#' @param ... Additional arguments passed to likelihood function
#' @return List with objective_fn and gradient_fn
create_likelihood_wrappers <- function(inputs, model_type, ...) {
  extra_args <- list(...)

  if (model_type == "mnl") {
    objective_fn <- function(theta) {
      mnl_loglik_gradient_parallel(
        theta = theta,
        X = inputs$X,
        alt_idx = inputs$alt_idx,
        choice_idx = inputs$choice_idx,
        M = inputs$M,
        weights = inputs$weights,
        use_asc = extra_args$use_asc %||% TRUE,
        include_outside_option = inputs$include_outside_option %||% FALSE
      )$objective
    }

    gradient_fn <- function(theta) {
      mnl_loglik_gradient_parallel(
        theta = theta,
        X = inputs$X,
        alt_idx = inputs$alt_idx,
        choice_idx = inputs$choice_idx,
        M = inputs$M,
        weights = inputs$weights,
        use_asc = extra_args$use_asc %||% TRUE,
        include_outside_option = inputs$include_outside_option %||% FALSE
      )$gradient
    }
  } else if (model_type == "mxl") {
    objective_fn <- function(theta) {
      mxl_loglik_gradient_parallel(
        theta = theta,
        X = inputs$X,
        W = inputs$W,
        alt_idx = inputs$alt_idx,
        choice_idx = inputs$choice_idx,
        M = inputs$M,
        weights = inputs$weights,
        eta_draws = extra_args$eta_draws,
        rc_dist = extra_args$rc_dist,
        rc_correlation = inputs$rc_correlation,
        rc_mean = extra_args$rc_mean %||% FALSE,
        use_asc = extra_args$use_asc %||% TRUE,
        include_outside_option = inputs$include_outside_option %||% FALSE
      )$objective
    }

    gradient_fn <- function(theta) {
      mxl_loglik_gradient_parallel(
        theta = theta,
        X = inputs$X,
        W = inputs$W,
        alt_idx = inputs$alt_idx,
        choice_idx = inputs$choice_idx,
        M = inputs$M,
        weights = inputs$weights,
        eta_draws = extra_args$eta_draws,
        rc_dist = extra_args$rc_dist,
        rc_correlation = inputs$rc_correlation,
        rc_mean = extra_args$rc_mean %||% FALSE,
        use_asc = extra_args$use_asc %||% TRUE,
        include_outside_option = inputs$include_outside_option %||% FALSE
      )$gradient
    }
  } else if (model_type == "nl") {
    objective_fn <- function(theta) {
      nl_loglik_gradient_parallel(
        theta = theta,
        X = inputs$X,
        alt_idx = inputs$alt_idx,
        choice_idx = inputs$choice_idx,
        nest_idx = inputs$nest_idx,
        M = inputs$M,
        K = inputs$K,
        weights = inputs$weights,
        singleton_nests = inputs$singleton_nests,
        use_asc = extra_args$use_asc %||% TRUE,
        include_outside_option = inputs$include_outside_option %||% FALSE
      )$objective
    }

    gradient_fn <- function(theta) {
      nl_loglik_gradient_parallel(
        theta = theta,
        X = inputs$X,
        alt_idx = inputs$alt_idx,
        choice_idx = inputs$choice_idx,
        nest_idx = inputs$nest_idx,
        M = inputs$M,
        K = inputs$K,
        weights = inputs$weights,
        singleton_nests = inputs$singleton_nests,
        use_asc = extra_args$use_asc %||% TRUE,
        include_outside_option = inputs$include_outside_option %||% FALSE
      )$objective
    }
  }

  list(objective_fn = objective_fn, gradient_fn = gradient_fn)
}

# Null coalescing operator (if not available)
`%||%` <- function(x, y) if (is.null(x)) y else x
