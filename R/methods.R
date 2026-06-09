
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
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' print(fit)
#' }
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
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' coef(fit)
#' }
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
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' vcov(fit)
#' }
#' @export
vcov.choicer_fit <- function(object, ...) {
  object <- ensure_vcov(object)
  object$vcov
}

# --- wesml_vcov (robust / sandwich) -----------------------------------------

#' Robust (sandwich) variance for a weighted / choice-based mixed logit fit
#'
#' Recomputes the robust Huber-White sandwich variance
#' \eqn{V = A^{-1} B A^{-1}} for a fitted mixed logit, where the bread
#' \eqn{A = \sum_i w_i (-H_i)} is the weighted negated Hessian and the meat
#' \eqn{B = \sum_i w_i^2 s_i s_i'} is the weight-squared outer product of the
#' per-individual scores. This is the appropriate variance under choice-based
#' (endogenous stratified) / WESML weighting, where the inverse-Hessian and the
#' ordinary BHHH variance are invalid. It can be called on any fitted model
#' (e.g. one estimated with \code{se_method = "hessian"}) to obtain robust
#' standard errors post hoc, without refitting.
#'
#' If the stored weights are uniform (all equal), a warning is emitted: the
#' returned variance is then the ordinary robust (Huber-White) variance, not a
#' WESML-weighted variance. Refit with WESML weights for a choice-based-sampling
#' correction.
#'
#' @param object A fitted \code{choicer_mxl} object (requires
#'   \code{keep_data = TRUE}).
#' @param type Either \code{"vcov"} (default) to return the variance-covariance
#'   matrix or \code{"se"} to return the standard-error vector.
#' @param ... Unused.
#' @returns A variance-covariance matrix (\code{type = "vcov"}) or a named
#'   numeric vector of standard errors (\code{type = "se"}), in the raw
#'   parameter space.
#' @seealso \code{\link{wesml_weights}}, \code{\link{sample_by_choice}},
#'   \code{\link{run_mxlogit}}
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(1)
#' N <- 200L; J <- 3L
#' dt <- data.table(id = rep(seq_len(N), each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#' dt[, choice := as.integer(seq_len(.N) == sample.int(.N, 1L)), by = id]
#' fit <- run_mxlogit(dt, "id", "alt", "choice", "x1", "w1", S = 50L)
#' wesml_vcov(fit, "se")
#' }
#' @export
wesml_vcov <- function(object, ...) UseMethod("wesml_vcov")

#' @rdname wesml_vcov
#' @export
wesml_vcov.choicer_mxl <- function(object, type = c("vcov", "se"), ...) {
  type <- match.arg(type)
  if (is.null(object[["data"]])) {
    stop("wesml_vcov() needs the stored data; refit with keep_data = TRUE.")
  }
  w <- object[["data"]]$weights
  if (!is.null(w) && length(unique(w)) == 1L) {
    warning("Stored weights are uniform; wesml_vcov() returns the ordinary robust ",
            "(Huber-White) variance, not a WESML-weighted variance. Refit with WESML ",
            "weights for a choice-based-sampling correction.", call. = FALSE)
  }
  res <- compute_sandwich_vcov(object)
  nms <- names(object$coefficients)
  if (!is.null(res$vcov)) {
    rownames(res$vcov) <- nms
    colnames(res$vcov) <- nms
  }
  if (!is.null(res$se)) names(res$se) <- nms
  if (type == "se") res$se else res$vcov
}

# --- logLik ------------------------------------------------------------------

#' Extract log-likelihood from a choicer_fit object
#'
#' Returns a logLik object, which enables AIC() and BIC() automatically.
#'
#' @param object A choicer_fit object.
#' @param ... Additional arguments (ignored).
#' @returns A logLik object with df and nobs attributes.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
#' }
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
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' nobs(fit)
#' }
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
#' @param gof Logical; compute goodness-of-fit measures (McFadden R-squared,
#'   hit rate) for the summary footer. Involves an in-sample prediction pass
#'   (for mixed logit, a full simulation over draws); set to FALSE to skip.
#' @param ... Additional arguments (ignored).
#' @returns A summary.choicer_mnl object (list with coefficients table and
#'   metadata, including a `gof` element with goodness-of-fit measures from
#'   \code{\link{gof}}; its fields are NA when the model was fitted with
#'   \code{keep_data = FALSE}).
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' summary(fit)
#' }
#' @export
summary.choicer_mnl <- function(object, gof = TRUE, ...) {
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
      elapsed_time = object$optimizer$elapsed_time,
      gof = if (isTRUE(gof)) {
        # calling gof() here is safe: R skips the logical binding when
        # resolving a function call
        tryCatch(suppressMessages(gof(object)), error = function(e) NULL)
      }
    ),
    class = "summary.choicer_mnl"
  )
}

#' Print summary for multinomial logit model
#' @param x A summary.choicer_mnl object.
#' @param ... Additional arguments (ignored).
#' @returns The object invisibly.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' print(summary(fit))
#' }
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
#' @param gof Logical; compute goodness-of-fit measures (McFadden R-squared,
#'   hit rate) for the summary footer. Involves an in-sample prediction pass
#'   (for mixed logit, a full simulation over draws); set to FALSE to skip.
#' @param ... Additional arguments (ignored).
#' @returns A summary.choicer_mxl object (includes a `gof` element with
#'   goodness-of-fit measures from \code{\link{gof}}; its fields are NA when
#'   the model was fitted with \code{keep_data = FALSE}).
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mxlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = "x1", random_var_cols = "w1", S = 50L
#' )
#' summary(fit)
#' }
#' @export
summary.choicer_mxl <- function(object, gof = TRUE, ...) {
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
      for (i in seq_len(K_w)) {
        for (j in seq_len(i)) {
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
      sigma = object$sigma,
      se_method = object$se_method %||% "hessian",
      weighting = object$choice_sampling$scheme,
      weights_applied = object$choice_sampling$weights_applied,
      gof = if (isTRUE(gof)) {
        # calling gof() here is safe: R skips the logical binding when
        # resolving a function call
        tryCatch(suppressMessages(gof(object)), error = function(e) NULL)
      }
    ),
    class = "summary.choicer_mxl"
  )
}

#' Print summary for mixed logit model
#' @param x A summary.choicer_mxl object.
#' @param ... Additional arguments (ignored).
#' @returns The object invisibly.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mxlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = "x1", random_var_cols = "w1", S = 50L
#' )
#' print(summary(fit))
#' }
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
  cat("Std. Errors:", switch(
    x$se_method %||% "hessian",
    bhhh = "BHHH (OPG)",
    sandwich = "Sandwich (robust)",
    "Analytical Hessian"
  ), "\n")
  if (!is.null(x$weighting)) {
    # Backward-compat: when weights_applied is absent (older fits) treat as applied.
    applied <- !isFALSE(x$weights_applied)
    cat("Weighting:",
        if (identical(x$weighting, "wesml")) {
          if (applied) {
            "WESML choice-based"
          } else {
            "WESML provenance present but NOT applied (fit is unweighted)"
          }
        } else {
          "user-supplied"
        },
        "\n")
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
#' @param gof Logical; compute goodness-of-fit measures (McFadden R-squared,
#'   hit rate) for the summary footer. Involves an in-sample prediction pass
#'   (for mixed logit, a full simulation over draws); set to FALSE to skip.
#' @param ... Additional arguments (ignored).
#' @returns A summary.choicer_nl object (includes a `gof` element with
#'   goodness-of-fit measures from \code{\link{gof}}; its fields are NA when
#'   the model was fitted with \code{keep_data = FALSE}).
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 4
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, nest := ifelse(alt <= 2, "A", "B")]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_nestlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = c("x1", "x2"), nest_col = "nest"
#' )
#' summary(fit)
#' }
#' @export
summary.choicer_nl <- function(object, gof = TRUE, ...) {
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
      elapsed_time = object$optimizer$elapsed_time,
      gof = if (isTRUE(gof)) {
        # calling gof() here is safe: R skips the logical binding when
        # resolving a function call
        tryCatch(suppressMessages(gof(object)), error = function(e) NULL)
      }
    ),
    class = "summary.choicer_nl"
  )
}

#' Print summary for nested logit model
#' @param x A summary.choicer_nl object.
#' @param ... Additional arguments (ignored).
#' @returns The object invisibly.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 4
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, nest := ifelse(alt <= 2, "A", "B")]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_nestlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = c("x1", "x2"), nest_col = "nest"
#' )
#' print(summary(fit))
#' }
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
  print_gof_lines(x$gof)
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
#' Computes choice probabilities or aggregate market shares, either for the
#' data used at fit time (default) or for counterfactual `newdata`.
#'
#' @param object A choicer_mnl object.
#' @param type One of "probabilities" (individual-level choice probabilities)
#'   or "shares" (aggregate market shares).
#' @param newdata Optional data for counterfactual prediction. Either:
#'   * a data.frame in the same long format used at fit time (one row per
#'     id-alternative pair, with the fit-time id, alternative, and covariate
#'     columns; a choice column is not required). Alternative labels must have
#'     been seen at fit time; per-id subsets of alternatives are allowed.
#'   * a list with elements `X`, `alt_idx`, `M` (and optionally `weights`)
#'     matching the layout of `object$data` — the "modified design matrix"
#'     path for policy simulation (e.g., perturb a column of `object$data$X`).
#'     `alt_idx` must use the fit-time integer codes from `object$alt_mapping`.
#'
#'   When `NULL` (default), the data stored at fit time is used (requires
#'   `keep_data = TRUE`).
#' @param weights Optional numeric vector with one weight per choice situation,
#'   used for `type = "shares"` aggregation. For a data.frame `newdata`,
#'   supply one weight per id in order of first appearance in `newdata`
#'   (weights are realigned internally to the sorted row order). Defaults to
#'   equal weights. Ignored when `newdata` is `NULL` (the stored fit weights
#'   apply).
#' @param ... Additional arguments (ignored).
#' @returns For "probabilities": a list with `choice_prob` and `utility` vectors.
#'   For "shares": a named numeric vector of market shares per alternative.
#'   With a data.frame `newdata`, rows are ordered by id, then by fit-time
#'   alternative code (`alt_int` in `object$alt_mapping`).
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' predict(fit, type = "shares")
#' predict(fit, type = "probabilities")
#'
#' # Counterfactual: increase x1 for alternative 2
#' dt_cf <- copy(dt)[alt == 2, x1 := x1 + 1]
#' predict(fit, type = "shares", newdata = dt_cf)
#' }
#' @export
predict.choicer_mnl <- function(object, type = c("probabilities", "shares"),
                                newdata = NULL, weights = NULL, ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    if (is.null(object[["data"]])) {
      stop("Prediction requires stored data. Refit with keep_data = TRUE.")
    }
    d <- object[["data"]]
  } else {
    d <- resolve_predict_newdata(object, newdata, weights = weights)
  }

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
#' Computes simulated choice probabilities or aggregate market shares using
#' deterministic Halton draws, either for the data used at fit time (default)
#' or for counterfactual `newdata`.
#'
#' @param object A choicer_mxl object.
#' @param type Either "probabilities" (per-observation simulated choice
#'   probabilities) or "shares" (aggregate simulated market shares).
#' @param newdata Optional data for counterfactual prediction. Either:
#'   * a data.frame in the same long format used at fit time (one row per
#'     id-alternative pair, with the fit-time id, alternative, fixed-coefficient,
#'     and random-coefficient columns; a choice column is not required).
#'     Alternative labels must have been seen at fit time; per-id subsets of
#'     alternatives are allowed.
#'   * a list with elements `X`, `W`, `alt_idx`, `M` (and optionally
#'     `weights`) matching the layout of `object$data` — the "modified design
#'     matrix" path for policy simulation. `alt_idx` must use the fit-time
#'     integer codes from `object$alt_mapping`.
#'
#'   When `NULL` (default), the data stored at fit time is used (requires
#'   `keep_data = TRUE`). Halton draws are regenerated deterministically from
#'   `object$draws_info` with one block of draws per choice situation in
#'   `newdata`.
#' @param weights Optional numeric vector with one weight per choice situation,
#'   used for `type = "shares"` aggregation. For a data.frame `newdata`,
#'   supply one weight per id in order of first appearance in `newdata`
#'   (weights are realigned internally to the sorted row order). Defaults to
#'   equal weights. Ignored when `newdata` is `NULL` (the stored fit weights
#'   apply).
#' @param ... Additional arguments (ignored).
#' @returns For "probabilities": a list with `choice_prob` and `utility`
#'   vectors averaged across simulation draws. For "shares": a named numeric
#'   vector of simulated market shares per alternative. With a data.frame
#'   `newdata`, rows are ordered by id, then by fit-time alternative code
#'   (`alt_int` in `object$alt_mapping`).
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mxlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = "x1", random_var_cols = "w1", S = 50L
#' )
#' predict(fit, type = "shares")
#' predict(fit, type = "probabilities")
#' }
#' @export
predict.choicer_mxl <- function(object, type = c("probabilities", "shares"),
                                newdata = NULL, weights = NULL, ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    if (is.null(object[["data"]]) || is.null(object$draws_info)) {
      stop("Prediction requires stored data and draws. ",
           "Refit with keep_data = TRUE.")
    }
    d <- object[["data"]]
    N_draws <- object$draws_info$N
  } else {
    if (is.null(object$draws_info)) {
      stop("Prediction with newdata requires stored draw metadata ",
           "('draws_info'). Refit to enable newdata prediction.")
    }
    d <- resolve_predict_newdata(object, newdata, weights = weights)
    N_draws <- d$N
  }

  eta_draws <- get_halton_normals(
    S   = object$draws_info$S,
    N   = N_draws,
    K_w = object$draws_info$K_w
  )

  args <- list(
    theta                  = object$coefficients,
    X                      = d$X,
    W                      = d$W,
    alt_idx                = d$alt_idx,
    M                      = d$M,
    eta_draws              = eta_draws,
    rc_dist                = object$rc_dist,
    rc_correlation         = object$rc_correlation,
    rc_mean                = object$rc_mean,
    use_asc                = object$use_asc,
    include_outside_option = object$include_outside_option
  )

  if (type == "probabilities") {
    do.call(mxl_predict, args)
  } else {
    do.call(mxl_predict_shares, c(args, list(weights = d$weights)))
  }
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
      est[idx_sigma] <- vech_row(Sigma_hat)
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
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' elasticities(fit, "x1")
#' }
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
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' diversion_ratios(fit)
#' }
#' @export
diversion_ratios <- function(object, ...) UseMethod("diversion_ratios")

#' BLP contraction mapping
#'
#' Finds the ASC (delta) parameters such that predicted market shares match
#' target shares, using the contraction mapping of Berry, Levinsohn, and
#' Pakes (1995) \doi{10.2307/2171802}.
#'
#' @param object A fitted model object.
#' @param target_shares Numeric vector of target market shares (length J).
#' @param ... Additional arguments passed to methods.
#' @returns Converged delta (ASC) vector.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' blp(fit, target_shares = rep(1/J, J))
#' }
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
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' elasticities(fit, "x1")
#' }
#' @export
elasticities.choicer_mnl <- function(object, elast_var, ...) {
  if (is.null(object[["data"]])) {
    stop("elasticities() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object[["data"]]
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
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' diversion_ratios(fit)
#' }
#' @export
diversion_ratios.choicer_mnl <- function(object, ...) {
  if (is.null(object[["data"]])) {
    stop("diversion_ratios() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object[["data"]]

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
#' @param target_shares Numeric vector of target market shares.
#'   Length \code{J_inside} when no outside option, or \code{J_inside + 1}
#'   (with the outside option's share at index 1) when
#'   \code{include_outside_option = TRUE}.
#' @param delta_init Initial guess for delta (ASC) values. If \code{NULL},
#'   uses the estimated ASCs from the fitted model.
#' @param tol Convergence tolerance (default 1e-8).
#' @param max_iter Maximum iterations (default 1000).
#' @param ... Additional arguments (ignored).
#' @returns Converged delta (ASC) vector.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mnlogit(dt, "id", "alt", "choice", c("x1", "x2"))
#' blp(fit, target_shares = rep(1/J, J))
#' }
#' @export
blp.choicer_mnl <- function(object, target_shares, delta_init = NULL,
                            tol = 1e-8, max_iter = 1000, ...) {
  if (is.null(object[["data"]])) {
    stop("blp() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object[["data"]]
  pm <- object$param_map
  beta <- object$coefficients[pm$beta]

  J <- nrow(object$alt_mapping)
  if (is.null(delta_init)) {
    if (!is.null(pm$asc)) {
      delta_init <- if (object$include_outside_option) {
        object$coefficients[pm$asc]            # length J_inside (all ASCs free)
      } else {
        c(0, object$coefficients[pm$asc])      # length J_inside, baseline = 0
      }
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
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mxlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = "x1", random_var_cols = "w1", S = 50L
#' )
#' elasticities(fit, "x1")
#' elasticities(fit, "w1", is_random_coef = TRUE)
#' }
#' @export
elasticities.choicer_mxl <- function(object, elast_var,
                                     is_random_coef = FALSE, ...) {
  if (is.null(object[["data"]])) {
    stop("elasticities() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object[["data"]]

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

# --- diversion_ratios: MXL ---------------------------------------------------

#' Diversion ratios for mixed logit model
#'
#' Computes the attribute-based diversion ratio matrix. Entry (k, j) is the
#' fraction of demand lost by alternative j that is captured by alternative k
#' when a marginal change in alternative j's \code{wrt_var} attribute reduces
#' s_j.
#'
#' Unlike MNL, the MXL diversion ratio depends on which variable is perturbed:
#' the realised coefficient \eqn{\beta_{ik}^s} varies across individuals and
#' draws and does not cancel in the ratio. For a variable with a fixed
#' coefficient the result is independent of the variable (\eqn{\beta} cancels);
#' for a random-coefficient variable it is not.
#'
#' @param object A \code{choicer_mxl} object fitted with \code{keep_data = TRUE}.
#' @param wrt_var Variable used to perturb alternative j's utility: a column
#'   name (character) or 1-based index. Indexes into X columns for fixed
#'   coefficients, or W columns for random coefficients (when
#'   \code{is_random_coef = TRUE}).
#' @param is_random_coef Logical. \code{TRUE} if the variable has a random
#'   coefficient (is in W), \code{FALSE} if fixed (in X). Default \code{FALSE}.
#' @param ... Additional arguments (ignored).
#' @returns A J x J diversion ratio matrix with alternative labels.
#'   Cross-products are averaged across simulation draws inside the
#'   integration to avoid Jensen-style bias.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mxlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = "x1", random_var_cols = "w1", S = 50L
#' )
#' diversion_ratios(fit, "x1")
#' diversion_ratios(fit, "w1", is_random_coef = TRUE)
#' }
#' @export
diversion_ratios.choicer_mxl <- function(object, wrt_var,
                                         is_random_coef = FALSE, ...) {
  if (is.null(object[["data"]])) {
    stop("diversion_ratios() requires stored data. Refit with keep_data = TRUE.")
  }
  if (is.null(object$draws_info)) {
    stop("diversion_ratios() requires draws_info from a fitted MXL model.")
  }
  d <- object[["data"]]

  col_names <- if (is_random_coef) colnames(d$W) else colnames(d$X)
  idx <- resolve_var_index(wrt_var, col_names)

  eta_draws <- get_halton_normals(
    S   = object$draws_info$S,
    N   = object$draws_info$N,
    K_w = object$draws_info$K_w
  )

  mat <- mxl_diversion_ratios_parallel(
    theta                  = object$coefficients,
    X                      = d$X,
    W                      = d$W,
    alt_idx                = d$alt_idx,
    M                      = d$M,
    weights                = d$weights,
    eta_draws              = eta_draws,
    rc_dist                = object$rc_dist,
    elast_var_idx          = idx,
    is_random_coef         = is_random_coef,
    rc_correlation         = object$rc_correlation,
    rc_mean                = object$rc_mean,
    use_asc                = object$use_asc,
    include_outside_option = object$include_outside_option
  )

  label_matrix(mat, object$alt_mapping)
}

# --- blp: MXL ----------------------------------------------------------------

#' BLP contraction mapping for mixed logit model
#'
#' @param object A \code{choicer_mxl} object fitted with \code{keep_data = TRUE}.
#' @param target_shares Numeric vector of target market shares.
#'   Length \code{J_inside} when no outside option, or \code{J_inside + 1}
#'   (with the outside option's share at index 1) when
#'   \code{include_outside_option = TRUE}.
#' @param delta_init Initial guess for delta (ASC) values. If \code{NULL},
#'   uses the estimated ASCs from the fitted model.
#' @param tol Convergence tolerance (default 1e-8).
#' @param max_iter Maximum iterations (default 1000).
#' @param ... Additional arguments (ignored).
#' @returns Converged delta (ASC) vector.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_mxlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = "x1", random_var_cols = "w1", S = 50L
#' )
#' blp(fit, target_shares = rep(1/J, J))
#' }
#' @export
blp.choicer_mxl <- function(object, target_shares, delta_init = NULL,
                            tol = 1e-8, max_iter = 1000, ...) {
  if (is.null(object[["data"]])) {
    stop("blp() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object[["data"]]
  pm <- object$param_map

  beta <- object$coefficients[pm$beta]
  mu <- if (!is.null(pm$mu)) object$coefficients[pm$mu] else rep(0, object$draws_info$K_w)
  L_params <- object$coefficients[pm$sigma]

  J <- nrow(object$alt_mapping)
  if (is.null(delta_init)) {
    if (!is.null(pm$asc)) {
      delta_init <- if (object$include_outside_option) {
        object$coefficients[pm$asc]            # length J_inside (all ASCs free)
      } else {
        c(0, object$coefficients[pm$asc])      # length J_inside, baseline = 0
      }
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

# --- predict: NL -------------------------------------------------------------

#' Predict from a nested logit model
#'
#' Computes choice probabilities or aggregate market shares, either for the
#' data used at fit time (default) or for counterfactual `newdata`.
#'
#' @param object A choicer_nl object.
#' @param type One of "probabilities" (individual-level choice probabilities)
#'   or "shares" (aggregate market shares).
#' @param newdata Optional data for counterfactual prediction. Either:
#'   * a data.frame in the same long format used at fit time (one row per
#'     id-alternative pair, with the fit-time id, alternative, and covariate
#'     columns; a choice column is not required). Alternative labels must have
#'     been seen at fit time; per-id subsets of alternatives are allowed. The
#'     alternative-to-nest mapping always comes from the fitted object (it
#'     indexes the estimated `lambda` parameters), so a nest column in
#'     `newdata` is not required and is ignored if present.
#'   * a list with elements `X`, `alt_idx`, `M` (and optionally `weights`)
#'     matching the layout of `object$data` — the "modified design matrix"
#'     path for policy simulation. `alt_idx` must use the fit-time integer
#'     codes from `object$alt_mapping`.
#'
#'   When `NULL` (default), the data stored at fit time is used (requires
#'   `keep_data = TRUE`).
#' @param weights Optional numeric vector with one weight per choice situation,
#'   used for `type = "shares"` aggregation. For a data.frame `newdata`,
#'   supply one weight per id in order of first appearance in `newdata`
#'   (weights are realigned internally to the sorted row order). Defaults to
#'   equal weights. Ignored when `newdata` is `NULL` (the stored fit weights
#'   apply).
#' @param ... Additional arguments (ignored).
#' @returns For "probabilities": a list with `choice_prob` and `utility` vectors.
#'   For "shares": a named numeric vector of market shares per alternative.
#'   With a data.frame `newdata`, rows are ordered by id, then by fit-time
#'   alternative code (`alt_int` in `object$alt_mapping`).
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 4
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, nest := rep(c(1L, 1L, 2L, 2L), N)]
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
#' predict(fit, type = "shares")
#' predict(fit, type = "probabilities")
#' }
#' @export
predict.choicer_nl <- function(object, type = c("probabilities", "shares"),
                               newdata = NULL, weights = NULL, ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    if (is.null(object[["data"]])) {
      stop("Prediction requires stored data. Refit with keep_data = TRUE.")
    }
    d <- object[["data"]]
    nest_idx <- d$nest_idx
  } else {
    # The per-alternative (length J) nest mapping is authoritative from the
    # fit since it indexes the estimated lambda parameters. New fits store it
    # top-level; older fits only carry it inside $data (keep_data = TRUE).
    nest_idx <- object$nest_idx %||% object[["data"]]$nest_idx
    if (is.null(nest_idx)) {
      stop("Cannot resolve the alternative-to-nest mapping: this fit stores ",
           "no 'nest_idx'. Refit to enable newdata prediction.")
    }
    d <- resolve_predict_newdata(object, newdata, weights = weights)
  }

  theta <- object$coefficients

  if (type == "probabilities") {
    nl_predict(
      theta = theta,
      X = d$X,
      alt_idx = d$alt_idx,
      M = d$M,
      nest_idx = nest_idx,
      use_asc = object$use_asc,
      include_outside_option = object$include_outside_option
    )
  } else {
    nl_predict_shares(
      theta = theta,
      X = d$X,
      alt_idx = d$alt_idx,
      M = d$M,
      weights = d$weights,
      nest_idx = nest_idx,
      use_asc = object$use_asc,
      include_outside_option = object$include_outside_option
    )
  }
}

# --- elasticities: NL --------------------------------------------------------

#' Elasticities for nested logit model
#'
#' @param object A \code{choicer_nl} object fitted with \code{keep_data = TRUE}.
#' @param elast_var Variable for elasticity computation: a column name (character)
#'   or 1-based index into the design matrix X.
#' @param ... Additional arguments (ignored).
#' @returns A J x J elasticity matrix with alternative labels.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 4
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, nest := rep(c(1L, 1L, 2L, 2L), N)]
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
#' elasticities(fit, "x1")
#' }
#' @export
elasticities.choicer_nl <- function(object, elast_var, ...) {
  if (is.null(object[["data"]])) {
    stop("elasticities() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object[["data"]]
  idx <- resolve_var_index(elast_var, colnames(d$X))

  mat <- nl_elasticities_parallel(
    theta = object$coefficients,
    X = d$X,
    alt_idx = d$alt_idx,
    choice_idx = d$choice_idx,
    nest_idx = d$nest_idx,
    M = d$M,
    weights = d$weights,
    elast_var_idx = idx,
    use_asc = object$use_asc,
    include_outside_option = object$include_outside_option
  )

  label_matrix(mat, object$alt_mapping)
}

# --- diversion_ratios: NL ----------------------------------------------------

#' Diversion ratios for nested logit model
#'
#' @param object A \code{choicer_nl} object fitted with \code{keep_data = TRUE}.
#' @param ... Additional arguments (ignored).
#' @returns A J x J diversion ratio matrix with alternative labels.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 4
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, nest := rep(c(1L, 1L, 2L, 2L), N)]
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
#' diversion_ratios(fit)
#' }
#' @export
diversion_ratios.choicer_nl <- function(object, ...) {
  if (is.null(object[["data"]])) {
    stop("diversion_ratios() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object[["data"]]

  mat <- nl_diversion_ratios_parallel(
    theta = object$coefficients,
    X = d$X,
    alt_idx = d$alt_idx,
    nest_idx = d$nest_idx,
    M = d$M,
    weights = d$weights,
    use_asc = object$use_asc,
    include_outside_option = object$include_outside_option
  )

  label_matrix(mat, object$alt_mapping)
}

# --- blp: NL -----------------------------------------------------------------

#' BLP contraction mapping for nested logit model
#'
#' @param object A \code{choicer_nl} object fitted with \code{keep_data = TRUE}.
#' @param target_shares Numeric vector of target market shares.
#'   Length \code{J_inside} when no outside option, or \code{J_inside + 1}
#'   (with the outside option's share at index 1) when
#'   \code{include_outside_option = TRUE}.
#' @param delta_init Initial guess for delta (ASC) values. If \code{NULL},
#'   uses the estimated ASCs from the fitted model.
#' @param damping Contraction damping factor in (0, 1] (default 1).
#' @param tol Convergence tolerance (default 1e-8).
#' @param max_iter Maximum iterations (default 1000).
#' @param ... Additional arguments (ignored).
#' @returns Converged delta (ASC) vector.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 4
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, nest := rep(c(1L, 1L, 2L, 2L), N)]
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' fit <- run_nestlogit(dt, "id", "alt", "choice", c("x1", "x2"), "nest")
#' blp(fit, target_shares = rep(1/J, J))
#' }
#' @export
blp.choicer_nl <- function(object, target_shares, delta_init = NULL,
                           damping = 1, tol = 1e-8, max_iter = 1000, ...) {
  if (is.null(object[["data"]])) {
    stop("blp() requires stored data. Refit with keep_data = TRUE.")
  }
  d <- object[["data"]]
  pm <- object$param_map
  beta <- object$coefficients[pm$beta]

  # Build the full length-n_nests lambda vector expected by the kernel.
  # Singleton nests have lambda fixed to 1; non-singleton nests receive the
  # estimated lambdas, scattered by ascending nest index (matching the C++).
  n_nests <- max(d$nest_idx)
  lambda_full <- rep(1, n_nests)
  nest_counts <- tabulate(d$nest_idx, n_nests)
  non_singleton <- which(nest_counts > 1)
  lambda_full[non_singleton] <- object$coefficients[pm$lambda]

  J <- nrow(object$alt_mapping)
  if (is.null(delta_init)) {
    if (!is.null(pm$asc)) {
      delta_init <- if (object$include_outside_option) {
        object$coefficients[pm$asc]            # length J_inside (all ASCs free)
      } else {
        c(0, object$coefficients[pm$asc])      # length J_inside, baseline = 0
      }
    } else {
      delta_init <- rep(0, J)
    }
  }

  nl_blp_contraction(
    delta = delta_init,
    target_shares = target_shares,
    X = d$X,
    beta = beta,
    lambda = lambda_full,
    alt_idx = d$alt_idx,
    nest_idx = d$nest_idx,
    M = d$M,
    weights = d$weights,
    include_outside_option = object$include_outside_option,
    damping = damping,
    tol = tol,
    max_iter = max_iter
  )
}
