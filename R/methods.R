
# S3 methods for choicer_fit objects

# Display names for model types
model_display_name <- function(model) {
  switch(model,
    mnl = "Multinomial Logit (MNL)",
    mxl = "Mixed Logit (MXL)",
    nl  = "Nested Logit (NL)",
    model
  )
}

# --- print -------------------------------------------------------------------

#' Print a choicer_fit object
#'
#' Prints a brief summary of the fitted model.
#'
#' @param x A choicer_fit object.
#' @param ... Additional arguments (ignored).
#' @returns The object invisibly.
#' @export
print.choicer_fit <- function(x, ...) {
  cat(model_display_name(x$model), "model\n")
  cat("  N obs:", x$nobs, " | Parameters:", x$n_params, "\n")
  cat("  Log-likelihood:", format(x$loglik, digits = 6), "\n")
  cat("  AIC:", format(-2 * x$loglik + 2 * x$n_params, digits = 6), "\n")
  if (!is.na(x$convergence)) {
    cat("  Convergence:", x$convergence, "(", x$message, ")\n")
  }
  invisible(x)
}

# --- coef --------------------------------------------------------------------

#' Extract coefficients from a choicer_fit object
#'
#' @param object A choicer_fit object.
#' @param ... Additional arguments (ignored).
#' @returns Named numeric vector of estimated coefficients.
#' @export
coef.choicer_fit <- function(object, ...) {
  object$coefficients
}

# --- vcov --------------------------------------------------------------------

#' Extract variance-covariance matrix from a choicer_fit object
#'
#' Triggers lazy Hessian computation if vcov has not been computed yet.
#'
#' @param object A choicer_fit object.
#' @param ... Additional arguments (ignored).
#' @returns Named variance-covariance matrix, or NULL if unavailable.
#' @export
vcov.choicer_fit <- function(object, ...) {
  object <- ensure_vcov(object)
  object$vcov
}

# --- logLik ------------------------------------------------------------------

#' Extract log-likelihood from a choicer_fit object
#'
#' Returns a logLik object, which enables AIC() and BIC() automatically.
#'
#' @param object A choicer_fit object.
#' @param ... Additional arguments (ignored).
#' @returns A logLik object with df and nobs attributes.
#' @export
logLik.choicer_fit <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- object$n_params
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}

# --- nobs --------------------------------------------------------------------

#' Extract number of observations from a choicer_fit object
#'
#' @param object A choicer_fit object.
#' @param ... Additional arguments (ignored).
#' @returns Integer number of choice situations.
#' @export
nobs.choicer_fit <- function(object, ...) {
  object$nobs
}

# --- summary: shared helpers -------------------------------------------------

#' Significance codes for p-values
#' @noRd
significance_code <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  ""
}

#' Build the standard coefficient table
#' @noRd
build_coef_table <- function(estimates, se, param_names) {
  zval <- estimates / se
  pval <- 2 * (1 - pnorm(abs(zval)))

  data.frame(
    Estimate  = estimates,
    Std_Error = se,
    z_value   = zval,
    Pr_z      = pval,
    Signif    = vapply(pval, significance_code, character(1)),
    row.names = param_names,
    stringsAsFactors = FALSE
  )
}

# --- summary: MNL ------------------------------------------------------------

#' Summary for multinomial logit model
#'
#' Computes and returns a coefficient summary table with standard errors,
#' z-values, p-values, and significance codes. Triggers lazy Hessian
#' computation if standard errors have not been computed yet.
#'
#' @param object A choicer_mnl object.
#' @param ... Additional arguments (ignored).
#' @returns A summary.choicer_mnl object (list with coefficients table and metadata).
#' @export
summary.choicer_mnl <- function(object, ...) {
  object <- ensure_vcov(object)

  coef_table <- build_coef_table(
    estimates = object$coefficients,
    se = object$se %||% rep(NA_real_, length(object$coefficients)),
    param_names = names(object$coefficients)
  )

  structure(
    list(
      model = object$model,
      coefficients = coef_table,
      loglik = object$loglik,
      nobs = object$nobs,
      n_params = object$n_params,
      convergence = object$convergence,
      message = object$message,
      elapsed_time = object$optimizer$elapsed_time
    ),
    class = "summary.choicer_mnl"
  )
}

#' Print summary for multinomial logit model
#' @param x A summary.choicer_mnl object.
#' @param ... Additional arguments (ignored).
#' @returns The object invisibly.
#' @export
print.summary.choicer_mnl <- function(x, ...) {
  cat(model_display_name(x$model), "model\n\n")
  print_coef_table(x$coefficients)
  cat("\n")
  print_footer(x)
  invisible(x)
}

# --- summary: MXL ------------------------------------------------------------

#' Summary for mixed logit model
#'
#' Computes coefficient summary with delta-method transformation for variance
#' parameters (Cholesky to covariance scale) and log-normal mean parameters.
#' Triggers lazy Hessian computation if standard errors have not been computed yet.
#'
#' @param object A choicer_mxl object.
#' @param ... Additional arguments (ignored).
#' @returns A summary.choicer_mxl object.
#' @export
summary.choicer_mxl <- function(object, ...) {
  object <- ensure_vcov(object)

  est <- object$coefficients
  se <- object$se %||% rep(NA_real_, length(est))

  # Apply delta method for log-normal mu and Cholesky -> Sigma
  if (!is.null(object$vcov) && !is.null(object$param_map)) {
    delta_result <- apply_mxl_delta_method(
      est_theta = est,
      se = se,
      vcov_mat = object$vcov,
      param_map = object$param_map,
      rc_dist = object$rc_dist,
      rc_correlation = object$rc_correlation,
      rc_mean = object$rc_mean
    )
    est <- delta_result$estimates
    se <- delta_result$se
  }

  # Build display names: L_ij -> Sigma_ij, Mu_x -> exp(Mu_x) for log-normal
  display_names <- names(object$coefficients)
  if (!is.null(object$param_map$sigma)) {
    idx_sigma <- object$param_map$sigma
    K_w <- length(object$rc_dist)
    if (object$rc_correlation) {
      sigma_display <- character(length(idx_sigma))
      k <- 1
      for (j in seq_len(K_w)) {
        for (i in j:K_w) {
          sigma_display[k] <- sprintf("Sigma_%d%d", i, j)
          k <- k + 1
        }
      }
    } else {
      sigma_display <- paste0("Sigma_", seq_len(K_w), seq_len(K_w))
    }
    display_names[idx_sigma] <- sigma_display
  }
  if (object$rc_mean && !is.null(object$param_map$mu)) {
    K_w <- length(object$rc_dist)
    for (k in seq_len(K_w)) {
      if (object$rc_dist[k] == 1) {
        idx <- object$param_map$mu[k]
        display_names[idx] <- paste0("exp(", display_names[idx], ")")
      }
    }
  }

  coef_table <- build_coef_table(
    estimates = est,
    se = se,
    param_names = display_names
  )

  structure(
    list(
      model = object$model,
      coefficients = coef_table,
      loglik = object$loglik,
      nobs = object$nobs,
      n_params = object$n_params,
      convergence = object$convergence,
      message = object$message,
      elapsed_time = object$optimizer$elapsed_time,
      sigma = object$sigma
    ),
    class = "summary.choicer_mxl"
  )
}

#' Print summary for mixed logit model
#' @param x A summary.choicer_mxl object.
#' @param ... Additional arguments (ignored).
#' @returns The object invisibly.
#' @export
print.summary.choicer_mxl <- function(x, ...) {
  cat(model_display_name(x$model), "model\n\n")
  print_coef_table(x$coefficients)
  cat("\n")
  if (!is.null(x$sigma)) {
    cat("Random coefficient covariance (Sigma):\n")
    print(x$sigma)
    cat("\n")
  }
  print_footer(x)
  invisible(x)
}

# --- summary: NL -------------------------------------------------------------

#' Summary for nested logit model
#'
#' Triggers lazy Hessian computation if standard errors have not been computed yet.
#'
#' @param object A choicer_nl object.
#' @param ... Additional arguments (ignored).
#' @returns A summary.choicer_nl object.
#' @export
summary.choicer_nl <- function(object, ...) {
  object <- ensure_vcov(object)

  coef_table <- build_coef_table(
    estimates = object$coefficients,
    se = object$se %||% rep(NA_real_, length(object$coefficients)),
    param_names = names(object$coefficients)
  )

  structure(
    list(
      model = object$model,
      coefficients = coef_table,
      loglik = object$loglik,
      nobs = object$nobs,
      n_params = object$n_params,
      convergence = object$convergence,
      message = object$message,
      elapsed_time = object$optimizer$elapsed_time
    ),
    class = "summary.choicer_nl"
  )
}

#' Print summary for nested logit model
#' @param x A summary.choicer_nl object.
#' @param ... Additional arguments (ignored).
#' @returns The object invisibly.
#' @export
print.summary.choicer_nl <- function(x, ...) {
  cat(model_display_name(x$model), "model\n\n")
  print_coef_table(x$coefficients)
  cat("\n")
  print_footer(x)
  invisible(x)
}

# --- Shared printing helpers -------------------------------------------------

#' Print a formatted coefficient table
#' @noRd
print_coef_table <- function(coef_table) {
  # Column widths
  param_width <- max(nchar(rownames(coef_table)), nchar("Parameter"))

  cat(sprintf(
    "%-*s  %10s %10s %8s %9s  %s\n",
    param_width, "Parameter", "Estimate", "Std.Error", "z-value", "Pr(>|z|)", ""
  ))

  for (i in seq_len(nrow(coef_table))) {
    cat(sprintf(
      "%-*s  %10.6f %10.6f %8.4f %9.2e  %s\n",
      param_width,
      rownames(coef_table)[i],
      coef_table$Estimate[i],
      coef_table$Std_Error[i],
      coef_table$z_value[i],
      coef_table$Pr_z[i],
      coef_table$Signif[i]
    ))
  }

  cat("---\nSignif. codes:  '***' 0.001 '**' 0.01 '*' 0.05\n")
}

#' Print model footer (log-likelihood, AIC, timing)
#' @noRd
print_footer <- function(x) {
  cat("Log-likelihood:", format(x$loglik, digits = 6), "\n")
  aic <- -2 * x$loglik + 2 * x$n_params
  bic <- -2 * x$loglik + log(x$nobs) * x$n_params
  cat("AIC:", format(aic, digits = 6), " | BIC:", format(bic, digits = 6), "\n")
  cat("N:", x$nobs, " | Parameters:", x$n_params, "\n")
  if (!is.null(x$elapsed_time)) {
    cat("Optimization time:", round(x$elapsed_time, 2), "s\n")
  }
  if (!is.na(x$convergence)) {
    cat("Convergence:", x$convergence, "(", x$message, ")\n")
  }
}

# --- predict: MNL ------------------------------------------------------------

#' Predict from a multinomial logit model
#'
#' Computes choice probabilities or aggregate market shares.
#'
#' @param object A choicer_mnl object.
#' @param type One of "probabilities" (individual-level choice probabilities)
#'   or "shares" (aggregate market shares).
#' @param ... Additional arguments (ignored).
#' @returns For "probabilities": a list with `choice_prob` and `utility` vectors.
#'   For "shares": a named numeric vector of market shares per alternative.
#' @export
predict.choicer_mnl <- function(object, type = c("probabilities", "shares"), ...) {
  type <- match.arg(type)

  if (is.null(object$data)) {
    stop("Prediction requires stored data. Refit with keep_data = TRUE.")
  }

  d <- object$data
  theta <- object$coefficients

  if (type == "probabilities") {
    mnl_predict(
      theta = theta,
      X = d$X,
      alt_idx = d$alt_idx,
      M = d$M,
      use_asc = object$use_asc,
      include_outside_option = object$include_outside_option
    )
  } else {
    mnl_predict_shares(
      theta = theta,
      X = d$X,
      alt_idx = d$alt_idx,
      M = d$M,
      weights = d$weights,
      use_asc = object$use_asc,
      include_outside_option = object$include_outside_option
    )
  }
}

# --- predict: MXL ------------------------------------------------------------

#' Predict from a mixed logit model
#'
#' Computes simulated choice probabilities using stored Halton draws.
#'
#' @param object A choicer_mxl object.
#' @param ... Additional arguments (ignored).
#' @returns A list with simulated choice probabilities.
#' @export
predict.choicer_mxl <- function(object, ...) {
  if (is.null(object$data) || is.null(object$draws_info)) {
    stop("Prediction requires stored data and draws. ",
         "Refit with keep_data = TRUE.")
  }
  stop("predict() for mixed logit is not yet implemented. ",
       "Use mxl_elasticities() for post-estimation analysis.")
}

# --- Delta method for MXL summary -------------------------------------------

#' Apply delta method for MXL variance parameters
#'
#' Transforms Cholesky parameters to covariance scale and
#' applies delta method for log-normal mu parameters.
#' Uses param_map for stable parameter indexing.
#'
#' @param est_theta Raw parameter estimates
#' @param se Raw standard errors
#' @param vcov_mat Variance-covariance matrix
#' @param param_map Named list with index vectors (beta, mu, sigma, asc)
#' @param rc_dist Integer vector of distribution types
#' @param rc_correlation Logical whether correlated
#' @param rc_mean Logical whether mu estimated
#' @returns List with transformed estimates and se
#' @noRd
apply_mxl_delta_method <- function(est_theta, se, vcov_mat,
                                   param_map, rc_dist, rc_correlation,
                                   rc_mean) {
  est <- est_theta
  se_out <- se

  K_w <- length(rc_dist)

  # Delta method for log-normal mu: exp(mu)
  if (rc_mean && !is.null(param_map$mu)) {
    idx_mu <- param_map$mu
    for (k in seq_len(K_w)) {
      if (rc_dist[k] == 1) {
        curr_idx <- idx_mu[k]
        mu_hat <- est[curr_idx]
        est[curr_idx] <- exp(mu_hat)
        se_out[curr_idx] <- exp(mu_hat) * se[curr_idx]
      }
    }
  }

  # Delta method for Cholesky -> Sigma
  if (!is.null(param_map$sigma)) {
    idx_sigma <- param_map$sigma
    L_params_hat <- est_theta[idx_sigma]
    Sigma_hat <- build_var_mat(L_params_hat, K_w, rc_correlation)

    if (rc_correlation) {
      est[idx_sigma] <- vech(Sigma_hat)
    } else {
      est[idx_sigma] <- diag(Sigma_hat)
    }

    J_mat <- jacobian_vech_Sigma(L_params_hat, K_w, rc_correlation)
    vcov_L <- vcov_mat[idx_sigma, idx_sigma]
    V_sigma <- J_mat %*% vcov_L %*% t(J_mat)
    se_out[idx_sigma] <- sqrt(diag(V_sigma))
  }

  list(estimates = est, se = se_out)
}

# --- Post-estimation generics ------------------------------------------------

#' Compute aggregate elasticities
#'
#' Computes a J x J matrix of aggregate elasticities. Entry (i, j) is the
#' percentage change in the probability of choosing alternative i when the
#' attribute of alternative j changes by 1\%.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @returns A J x J elasticity matrix with alternative labels.
#' @export
elasticities <- function(object, ...) UseMethod("elasticities")

#' Compute aggregate diversion ratios
#'
#' Computes a J x J matrix of diversion ratios. Entry (i, j) is the fraction
#' of demand lost by alternative j that is captured by alternative i when
#' alternative j becomes less attractive.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @returns A J x J diversion ratio matrix with alternative labels.
#' @export
diversion_ratios <- function(object, ...) UseMethod("diversion_ratios")

#' BLP contraction mapping
#'
#' Finds the ASC (delta) parameters such that predicted market shares match
#' target shares, using the contraction mapping of Berry, Levinsohn, and
#' Pakes (1995).
#'
#' @param object A fitted model object.
#' @param target_shares Numeric vector of target market shares (length J).
#' @param ... Additional arguments passed to methods.
#' @returns Converged delta (ASC) vector.
#' @export
blp <- function(object, target_shares, ...) UseMethod("blp")

# --- elasticities: MNL -------------------------------------------------------

#' Elasticities for multinomial logit model
#'
#' @param object A \code{choicer_mnl} object fitted with \code{keep_data = TRUE}.
#' @param elast_var Variable for elasticity computation: a column name (character)
#'   or 1-based index into the design matrix X.
#' @param ... Additional arguments (ignored).
#' @returns A J x J elasticity matrix with alternative labels.
#' @export
elasticities.choicer_mnl <- function(object, elast_var, ...) {
  if (is.null(object$data)) {
    stop("elasticities() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object$data
  idx <- resolve_var_index(elast_var, colnames(d$X))

  mat <- mnl_elasticities_parallel(
    theta = object$coefficients,
    X = d$X,
    alt_idx = d$alt_idx,
    choice_idx = d$choice_idx,
    M = d$M,
    weights = d$weights,
    elast_var_idx = idx,
    use_asc = object$use_asc,
    include_outside_option = object$include_outside_option
  )

  label_matrix(mat, object$alt_mapping)
}

# --- diversion_ratios: MNL ---------------------------------------------------

#' Diversion ratios for multinomial logit model
#'
#' @param object A \code{choicer_mnl} object fitted with \code{keep_data = TRUE}.
#' @param ... Additional arguments (ignored).
#' @returns A J x J diversion ratio matrix with alternative labels.
#' @export
diversion_ratios.choicer_mnl <- function(object, ...) {
  if (is.null(object$data)) {
    stop("diversion_ratios() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object$data

  mat <- mnl_diversion_ratios_parallel(
    theta = object$coefficients,
    X = d$X,
    alt_idx = d$alt_idx,
    M = d$M,
    weights = d$weights,
    use_asc = object$use_asc,
    include_outside_option = object$include_outside_option
  )

  label_matrix(mat, object$alt_mapping)
}

# --- blp: MNL ----------------------------------------------------------------

#' BLP contraction mapping for multinomial logit model
#'
#' @param object A \code{choicer_mnl} object fitted with \code{keep_data = TRUE}.
#' @param target_shares Numeric vector of target market shares (length J).
#' @param delta_init Initial guess for delta (ASC) values. If \code{NULL},
#'   uses the estimated ASCs from the fitted model.
#' @param tol Convergence tolerance (default 1e-8).
#' @param max_iter Maximum iterations (default 1000).
#' @param ... Additional arguments (ignored).
#' @returns Converged delta (ASC) vector.
#' @export
blp.choicer_mnl <- function(object, target_shares, delta_init = NULL,
                            tol = 1e-8, max_iter = 1000, ...) {
  if (is.null(object$data)) {
    stop("blp() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object$data
  pm <- object$param_map
  beta <- object$coefficients[pm$beta]

  J <- nrow(object$alt_mapping)
  if (is.null(delta_init)) {
    if (!is.null(pm$asc)) {
      # ASCs exclude the reference alternative (fixed to 0); prepend 0 for full J-vector
      delta_init <- c(0, object$coefficients[pm$asc])
    } else {
      delta_init <- rep(0, J)
    }
  }

  blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = d$X,
    beta = beta,
    alt_idx = d$alt_idx,
    M = d$M,
    weights = d$weights,
    include_outside_option = object$include_outside_option,
    tol = tol,
    max_iter = max_iter
  )
}

# --- elasticities: MXL -------------------------------------------------------

#' Elasticities for mixed logit model
#'
#' @param object A \code{choicer_mxl} object fitted with \code{keep_data = TRUE}.
#' @param elast_var Variable for elasticity computation: a column name (character)
#'   or 1-based index. Indexes into X columns for fixed coefficients, or W columns
#'   for random coefficients (when \code{is_random_coef = TRUE}).
#' @param is_random_coef Logical. \code{TRUE} if the variable has a random
#'   coefficient (is in W), \code{FALSE} if fixed (in X). Default \code{FALSE}.
#' @param ... Additional arguments (ignored).
#' @returns A J x J elasticity matrix with alternative labels.
#' @export
elasticities.choicer_mxl <- function(object, elast_var,
                                     is_random_coef = FALSE, ...) {
  if (is.null(object$data)) {
    stop("elasticities() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object$data

  col_names <- if (is_random_coef) colnames(d$W) else colnames(d$X)
  idx <- resolve_var_index(elast_var, col_names)

  eta_draws <- get_halton_normals(
    S = object$draws_info$S,
    N = object$draws_info$N,
    K_w = object$draws_info$K_w
  )

  mat <- mxl_elasticities_parallel(
    theta = object$coefficients,
    X = d$X,
    W = d$W,
    alt_idx = d$alt_idx,
    choice_idx = d$choice_idx,
    M = d$M,
    weights = d$weights,
    eta_draws = eta_draws,
    rc_dist = object$rc_dist,
    elast_var_idx = idx,
    is_random_coef = is_random_coef,
    rc_correlation = object$rc_correlation,
    rc_mean = object$rc_mean,
    use_asc = object$use_asc,
    include_outside_option = object$include_outside_option
  )

  label_matrix(mat, object$alt_mapping)
}

# --- blp: MXL ----------------------------------------------------------------

#' BLP contraction mapping for mixed logit model
#'
#' @param object A \code{choicer_mxl} object fitted with \code{keep_data = TRUE}.
#' @param target_shares Numeric vector of target market shares (length J).
#' @param delta_init Initial guess for delta (ASC) values. If \code{NULL},
#'   uses the estimated ASCs from the fitted model.
#' @param tol Convergence tolerance (default 1e-8).
#' @param max_iter Maximum iterations (default 1000).
#' @param ... Additional arguments (ignored).
#' @returns Converged delta (ASC) vector.
#' @export
blp.choicer_mxl <- function(object, target_shares, delta_init = NULL,
                            tol = 1e-8, max_iter = 1000, ...) {
  if (is.null(object$data)) {
    stop("blp() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object$data
  pm <- object$param_map

  beta <- object$coefficients[pm$beta]
  mu <- if (!is.null(pm$mu)) object$coefficients[pm$mu] else rep(0, object$draws_info$K_w)
  L_params <- object$coefficients[pm$sigma]

  J <- nrow(object$alt_mapping)
  if (is.null(delta_init)) {
    if (!is.null(pm$asc)) {
      delta_init <- c(0, object$coefficients[pm$asc])
    } else {
      delta_init <- rep(0, J)
    }
  }

  eta_draws <- get_halton_normals(
    S = object$draws_info$S,
    N = object$draws_info$N,
    K_w = object$draws_info$K_w
  )

  mxl_blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = d$X,
    W = d$W,
    beta = beta,
    mu = mu,
    L_params = L_params,
    alt_idx = d$alt_idx,
    M = d$M,
    weights = d$weights,
    eta_draws = eta_draws,
    rc_dist = object$rc_dist,
    rc_correlation = object$rc_correlation,
    rc_mean = object$rc_mean,
    include_outside_option = object$include_outside_option,
    tol = tol,
    max_iter = max_iter
  )
}
