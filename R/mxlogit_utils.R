
# Normalize a bounds vector to nloptr's full-length representation. Accepts
# NULL (returns the default), a full-length unnamed numeric (returned as-is
# after a length check), or a named partial numeric whose names must be a
# subset of `param_names`. Used by run_mxlogit() for both `lower` and `upper`.
.normalize_bound <- function(b, param_names, default, which) {
  if (is.null(b) || length(b) == 0L) {
    return(rep(default, length(param_names)))
  }
  if (!is.numeric(b)) stop("`", which, "` must be numeric or NULL.")
  if (anyNA(b)) stop("`", which, "` must not contain NA/NaN.")
  if (is.null(names(b))) {
    if (length(b) != length(param_names)) {
      stop("Unnamed `", which, "` must have length n_params (",
           length(param_names), "); got ", length(b), ".")
    }
    return(b)
  }
  dups <- unique(names(b)[duplicated(names(b))])
  if (length(dups)) {
    stop("Duplicate name(s) in `", which, "`: ",
         paste(dups, collapse = ", "), ".")
  }
  bad <- setdiff(names(b), param_names)
  if (length(bad)) {
    stop("Unknown parameter name(s) in `", which, "`: ",
         paste(bad, collapse = ", "),
         ". Valid names: ", paste(param_names, collapse = ", "), ".")
  }
  out <- rep(default, length(param_names))
  names(out) <- param_names
  out[names(b)] <- b
  unname(out)
}

#' Runs mixed logit estimation
#'
#' Estimates a mixed logit model via simulated maximum likelihood.
#'
#' Two workflows are supported:
#' \describe{
#'   \item{Convenience}{Supply \code{data} and column names. Data preparation
#'     (\code{\link{prepare_mxl_data}}) and Halton draw generation
#'     (\code{\link{get_halton_normals}}) are handled automatically.}
#'   \item{Advanced}{Call \code{\link{prepare_mxl_data}} and
#'     \code{\link{get_halton_normals}} yourself, then pass the results via
#'     \code{input_data} and \code{eta_draws}.}
#' }
#'
#' @param data Data frame containing choice data (convenience workflow).
#'   Mutually exclusive with \code{input_data}.
#' @param id_col Name of the column identifying choice situations.
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating chosen alternative (1/0).
#' @param covariate_cols Vector of column names for fixed covariates.
#' @param random_var_cols Vector of column names for random coefficients.
#' @param input_data List output from \code{\link{prepare_mxl_data}} (advanced
#'   workflow). Mutually exclusive with \code{data}.
#' @param eta_draws Array of shape K_w x S x N with standard normal draws.
#'   Required for the advanced workflow; auto-generated from \code{S} in the
#'   convenience workflow.
#' @param S Integer number of Halton draws per individual (convenience workflow
#'   only). Default 100.
#' @param rc_dist Integer vector indicating distribution of random coefficients
#'   (0 = normal, 1 = log-normal). Default: all normal.
#' @param rc_mean Logical indicating whether to estimate means for random
#'   coefficients.
#' @param rc_correlation Logical indicating whether random coefficients are
#'   correlated (convenience workflow). Ignored when \code{input_data} is used
#'   (taken from the prepared data).
#' @param use_asc Logical indicating whether to include alternative-specific
#'   constants.
#' @param theta_init Initial parameter vector in natural-scale units. If
#'   \code{NULL}, defaults to zeros for the \eqn{\beta}, \eqn{\mu}, and ASC
#'   blocks, and \code{log(0.5)} on the Cholesky diagonal (so each diagonal
#'   factor \eqn{L_{pp} = 0.5}, i.e. a moderate random-coefficient variance of
#'   \code{0.25}). The zero-on-diagonal alternative corresponds to
#'   \eqn{L_{pp} = 1} (unit RC variance), which often lets the first L-BFGS
#'   step overshoot.
#' @param lower,upper Optional parameter bounds for the optimizer, in
#'   natural-scale units (forward-transformed internally to scaled space when
#'   \code{scale_vars != "none"}). Each accepts three forms:
#'   \describe{
#'     \item{\code{NULL}}{(default) Unbounded (\code{-Inf}/\code{Inf}).}
#'     \item{Unnamed numeric vector of length \code{n_params}}{Full-length
#'       vector ordered exactly like \code{theta_init} (the nloptr-native form).}
#'     \item{Named numeric vector}{Names must be a subset of the parameter
#'       names (\eqn{\beta} block: column names of \code{X};
#'       \eqn{\mu} block: \code{Mu_<col>} (if \code{rc_mean = TRUE});
#'       Cholesky block: \code{L_<i><j>} for \eqn{i \ge j}; ASC block:
#'       \code{ASC_<level>}). Unlisted parameters default to \eqn{\pm\infty}.
#'       This is the recommended form for typical use, e.g.
#'       \code{lower = c(L_11 = -5, L_22 = -5)} to clip Cholesky diagonals.}
#'   }
#' @param optimizer Optimizer to use: \code{"nloptr"} (default), \code{"optim"},
#'   or a custom function. See \code{\link{run_mnlogit}} for details.
#' @param control List of optimizer-specific control parameters.
#' @param se_method Method for computing standard errors. One of
#'   \code{"hessian"} (default) for the analytical Hessian of the simulated
#'   log-likelihood, \code{"bhhh"} for the BHHH/outer-product-of-gradients
#'   (OPG) estimator, or \code{"sandwich"} for the robust (Huber-White)
#'   variance \eqn{V = A^{-1} B A^{-1}} (bread \eqn{A} = weighted negated
#'   Hessian, meat \eqn{B} = weight-squared OPG). Use \code{"sandwich"} for
#'   valid inference under choice-based / WESML weighting, where the
#'   inverse-Hessian and ordinary BHHH are invalid; it reduces to the usual
#'   robust variance under uniform weights. BHHH scales better to large
#'   problems (many alternatives or simulation draws) but may underestimate
#'   standard errors in finite samples or away from the optimum.
#' @param scale_vars Pre-estimation column scaling for design matrices. One of
#'   \code{"none"} (default), \code{"sd"} (sample standard deviation),
#'   \code{"mad"} (\code{stats::mad}, i.e. 1.4826 \eqn{\times}
#'   median absolute deviation; SD-equivalent under normality), or
#'   \code{"iqr"} (\code{stats::IQR(x) / 1.349}; also SD-equivalent under
#'   normality). When not \code{"none"}, every column of \code{X} and \code{W}
#'   is divided by the chosen scale before optimization to improve Hessian
#'   conditioning. Robust scales (\code{"mad"}/\code{"iqr"}) better capture
#'   the bulk for heavy-tailed columns where SD is dominated by outliers, but
#'   \code{stats::mad} can return zero when more than half of a column's
#'   entries are identical (e.g., a sparse 0/1 dummy) and will then trigger
#'   the same near-constant-column error as \code{"sd"}. Coefficients and
#'   standard errors are back-transformed to the user's natural units via the
#'   delta method, so reported quantities are invariant to this choice.
#'   Columns of \code{W} associated with log-normal random coefficients
#'   (\code{rc_dist == 1}) are passed through unchanged, since the shifted
#'   log-normal parameterization does not admit a closed-form back-transform
#'   under multiplicative scaling.
#' @param weights Optional weight vector (convenience workflow). If \code{NULL},
#'   equal weights are used.
#' @param weights_col Optional name of a column in \code{data} holding a per-row
#'   weight (constant within each choice situation). Mutually exclusive with
#'   \code{weights}; the recommended way to pass WESML weights from
#'   \code{\link{sample_by_choice}} / \code{\link{wesml_weights}}, since
#'   alignment is by id rather than by position. Convenience workflow only.
#' @param outside_opt_label Label for the outside option (convenience workflow).
#' @param include_outside_option Logical whether to include an outside option
#'   (convenience workflow).
#' @param keep_data Logical. If \code{TRUE} (default), stores prepared data in
#'   the returned object for post-estimation functions.
#' @param nloptr_opts Deprecated. Use \code{optimizer} and \code{control}
#'   instead.
#' @returns A \code{choicer_mxl} object (inherits from \code{choicer_fit}).
#'   Standard S3 methods available: \code{summary()}, \code{coef()},
#'   \code{vcov()}, \code{logLik()}, \code{AIC()}, \code{BIC()},
#'   \code{nobs()}.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 100; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#'
#' fit <- run_mxlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = "x1", random_var_cols = c("w1", "w2"), S = 50L
#' )
#' summary(fit)
#' }
#' @importFrom nloptr nloptr
#' @export
run_mxlogit <- function(
    data = NULL,
    id_col = NULL,
    alt_col = NULL,
    choice_col = NULL,
    covariate_cols = NULL,
    random_var_cols = NULL,
    input_data = NULL,
    eta_draws = NULL,
    S = 100L,
    rc_dist = NULL,
    rc_mean = FALSE,
    rc_correlation = FALSE,
    use_asc = TRUE,
    theta_init = NULL,
    lower = NULL,
    upper = NULL,
    optimizer = NULL,
    control = list(),
    se_method = c("hessian", "bhhh", "sandwich"),
    scale_vars = c("none", "sd", "mad", "iqr"),
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    keep_data = TRUE,
    nloptr_opts = NULL,
    weights_col = NULL
) {
  cl <- match.call()

  se_method <- match.arg(se_method)
  scale_vars <- match.arg(scale_vars)

  # Backward compatibility: nloptr_opts -> optimizer + control
  if (!is.null(nloptr_opts)) {
    message("'nloptr_opts' is deprecated. Use 'optimizer' and 'control' instead.")
    optimizer <- optimizer %||% "nloptr"
    control <- nloptr_opts
  }

  # --- Resolve input pathway --------------------------------------------------
  has_data <- !is.null(data)
  has_input <- !is.null(input_data)
  cs_meta <- if (has_data) attr(data, "choice_sampling") else NULL

  if (has_data && has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced), not both.")
  }
  if (!has_data && !has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced).")
  }
  if (has_input && !is.null(weights_col)) {
    stop("`weights_col` is only supported in the convenience (data) workflow. ",
         "Bake weights into `input_data` via prepare_mxl_data(weights_col = ) ",
         "or supply `weights` to prepare_mxl_data().")
  }

  if (has_data) {
    # Convenience workflow: validate required column-name arguments
    if (is.null(id_col) || is.null(alt_col) || is.null(choice_col) ||
        is.null(covariate_cols) || is.null(random_var_cols)) {
      stop("Convenience workflow requires: id_col, alt_col, choice_col, ",
           "covariate_cols, and random_var_cols.")
    }
    input_data <- prepare_mxl_data(
      data = data,
      id_col = id_col,
      alt_col = alt_col,
      choice_col = choice_col,
      covariate_cols = covariate_cols,
      random_var_cols = random_var_cols,
      weights = weights,
      weights_col = weights_col,
      outside_opt_label = outside_opt_label,
      include_outside_option = include_outside_option,
      rc_correlation = rc_correlation
    )
    K_w <- ncol(input_data$W)
    eta_draws <- get_halton_normals(S, input_data$N, K_w)
  } else {
    # Advanced workflow
    if (is.null(eta_draws)) {
      stop("'eta_draws' is required when using 'input_data' (advanced workflow).")
    }
  }

  # Parameter dimensions
  J <- nrow(input_data$alt_mapping)
  K_x <- ncol(input_data$X)
  K_w <- ncol(input_data$W)
  rc_correlation <- input_data$rc_correlation
  L_size <- if (rc_correlation) K_w * (K_w + 1) / 2 else K_w
  mu_size <- if (rc_mean) K_w else 0
  n_asc <- J - 1
  n_params <- K_x + mu_size + L_size + n_asc

  if (is.null(rc_dist)) rc_dist <- rep(0L, K_w)

  # Parameter index map (built early so the scaling layer can address blocks)
  pos <- 0
  param_map <- list(beta = seq_len(K_x))
  pos <- K_x
  if (mu_size > 0) {
    param_map$mu <- pos + seq_len(mu_size)
    pos <- pos + mu_size
  }
  param_map$sigma <- pos + seq_len(L_size)
  pos <- pos + L_size
  param_map$asc <- pos + seq_len(n_asc)

  # Parameter names (built early; reused for theta_hat, vcov, se downstream)
  beta_names <- colnames(input_data$X)
  mu_names <- if (rc_mean) paste0("Mu_", colnames(input_data$W)) else character(0)
  if (rc_correlation) {
    sigma_names <- character(L_size)
    nm_idx <- 1L
    for (i in seq_len(K_w)) {
      for (j in seq_len(i)) {
        sigma_names[nm_idx] <- sprintf("L_%d%d", i, j)
        nm_idx <- nm_idx + 1L
      }
    }
  } else {
    sigma_names <- paste0("L_", seq_len(K_w), seq_len(K_w))
  }
  alt_col <- names(input_data$alt_mapping)[2]
  asc_names <- paste0("ASC_", input_data$alt_mapping[2:J][[alt_col]])
  param_names <- c(beta_names, mu_names, sigma_names, asc_names)

  # --- Variable scaling (optional) --------------------------------------------
  # Scale columns of X and W by their sample SD to improve Hessian conditioning.
  # Keep the natural-scale matrices for storage; theta_init is interpreted in
  # natural units and forward-transformed below; theta_hat and vcov are
  # back-transformed after optimization so reported quantities are in the
  # user's natural units. sX and sW are returned as 1s when scale_vars="none".
  natural_X <- input_data$X
  natural_W <- input_data$W
  sX <- rep(1, K_x); names(sX) <- colnames(input_data$X)
  sW <- rep(1, K_w); names(sW) <- colnames(input_data$W)
  if (scale_vars != "none") {
    if (K_x > 0) {
      sX_raw <- .column_scales(input_data$X, scale_vars)
      .assert_scales_ok(sX_raw, scale_vars, "fixed-coefficient")
      sX <- sX_raw
      input_data$X <- sweep(input_data$X, 2, sX, "/")
    }
    if (K_w > 0) {
      sW_raw <- .column_scales(input_data$W, scale_vars)
      normal_cols <- which(rc_dist == 0L)
      if (length(normal_cols) > 0L) {
        .assert_scales_ok(sW_raw, scale_vars, "normal random-coefficient",
                          idx = normal_cols)
      }
      # Preserve names from sW_raw; carve out log-normal columns (pass-through).
      sW <- sW_raw
      sW[rc_dist == 1L] <- 1
      input_data$W <- sweep(input_data$W, 2, sW, "/")
      n_lognormal <- sum(rc_dist == 1L)
      if (K_w > 0L && n_lognormal == K_w) {
        message("scale_vars='", scale_vars,
                "': all random-coefficient column(s) are log-normal; W not scaled.")
      } else if (n_lognormal > 0L) {
        message("scale_vars='", scale_vars,
                "': passing through log-normal random-coefficient column(s) ",
                "unchanged (no closed-form back-transform).")
      }
    }
  }

  # --- Natural <-> scaled Jacobian --------------------------------------------
  # Maps scaled-space parameters back to natural-scale units:
  #   theta_natural = bt_mult * theta_scaled + bt_shift
  # Inverse forward-transforms theta_init from natural to scaled space.
  # ASCs and any unset entries default to identity (mult=1, shift=0).
  bt_mult <- rep(1, n_params)
  bt_shift <- rep(0, n_params)
  if (scale_vars != "none") {
    if (K_x > 0) bt_mult[param_map$beta] <- 1 / sX
    if (mu_size > 0) bt_mult[param_map$mu] <- 1 / sW
    if (rc_correlation) {
      idx <- 1L
      for (i in seq_len(K_w)) {
        for (j in seq_len(i)) {
          pos <- param_map$sigma[idx]
          if (i == j) {
            bt_shift[pos] <- -log(sW[i])
          } else {
            bt_mult[pos] <- 1 / sW[i]
          }
          idx <- idx + 1L
        }
      }
    } else {
      for (i in seq_len(K_w)) {
        pos <- param_map$sigma[i]
        bt_shift[pos] <- -log(sW[i])
      }
    }
  }

  # Resolve theta_init (natural units); forward-transform to scaled space.
  # Default cold-start: zero on every block except the Cholesky diagonal,
  # which sits at log(0.5) so each diagonal factor L_pp = 0.5 (RC variance
  # 0.25). Starting at log(1) = 0 corresponds to L_pp = 1 (unit RC variance),
  # which is often too large for typical specs and lets the first L-BFGS step
  # push ell_pp far enough that L_pp underflows / overflows.
  if (is.null(theta_init)) {
    theta_init <- rep(0, n_params)
    if (K_w > 0L) {
      if (rc_correlation) {
        diag_idx <- param_map$sigma[cumsum(seq_len(K_w))]
      } else {
        diag_idx <- param_map$sigma
      }
      theta_init[diag_idx] <- log(0.5)
    }
  }
  if (scale_vars != "none") {
    theta_init <- (theta_init - bt_shift) / bt_mult
  }

  # Normalize lower/upper bounds (natural units in, scaled units out).
  lower <- .normalize_bound(lower, param_names, -Inf, "lower")
  upper <- .normalize_bound(upper, param_names,  Inf, "upper")
  if (scale_vars != "none") {
    lower <- (lower - bt_shift) / bt_mult
    upper <- (upper - bt_shift) / bt_mult
  }

  # Build eval_f closure
  eval_f <- function(theta) {
    mxl_loglik_gradient_parallel(
      theta = theta,
      X = input_data$X,
      W = input_data$W,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      M = input_data$M,
      weights = input_data$weights,
      rc_dist = rc_dist,
      rc_correlation = rc_correlation,
      rc_mean = rc_mean,
      eta_draws = eta_draws,
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
      lower = lower,
      upper = upper,
      control = control
    )
  })

  message("Optimization run time ", convertTime(elapsed))

  # Estimate at the optimum (in scaled space if scale_vars='sd')
  theta_hat <- opt$par
  names(theta_hat) <- param_names

  # Choice-based-sampling provenance and a guardrail for weighted inference.
  weights_nonuniform <- length(unique(input_data$weights)) > 1
  if (weights_nonuniform && se_method != "sandwich") {
    warning("Non-uniform weights detected. If these are sampling/WESML ",
            "weights, use se_method = 'sandwich' for valid inference.",
            call. = FALSE)
  }
  choice_sampling <- if (!is.null(cs_meta)) {
    utils::modifyList(as.list(cs_meta), list(se_method = se_method))
  } else if (weights_nonuniform) {
    list(scheme = "user", se_method = se_method)
  } else {
    NULL
  }

  # Compute vcov eagerly using the selected SE method.
  # For "sandwich" (robust / WESML) standard errors, form V = A^{-1} B A^{-1}
  # with bread A = weighted negated Hessian and meat B = weight-squared OPG
  # (pass weights^2 to the BHHH routine, whose per-individual score is
  # weight-free). Computed in scaled space; the back-transform below applies.
  if (se_method == "sandwich") {
    A_bread <- mxl_hessian_parallel(
      theta = theta_hat, X = input_data$X, W = input_data$W,
      alt_idx = input_data$alt_idx, choice_idx = input_data$choice_idx,
      M = input_data$M, weights = input_data$weights, eta_draws = eta_draws,
      rc_dist = rc_dist, rc_correlation = rc_correlation, rc_mean = rc_mean,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option
    )
    B_meat <- mxl_bhhh_parallel(
      theta = theta_hat, X = input_data$X, W = input_data$W,
      alt_idx = input_data$alt_idx, choice_idx = input_data$choice_idx,
      M = input_data$M, weights = input_data$weights^2, eta_draws = eta_draws,
      rc_dist = rc_dist, rc_correlation = rc_correlation, rc_mean = rc_mean,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option
    )
    vcov_result <- .sandwich_combine(A_bread, B_meat)
  } else {
  hess <- switch(
    se_method,
    hessian = mxl_hessian_parallel(
      theta = theta_hat,
      X = input_data$X,
      W = input_data$W,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      M = input_data$M,
      weights = input_data$weights,
      eta_draws = eta_draws,
      rc_dist = rc_dist,
      rc_correlation = rc_correlation,
      rc_mean = rc_mean,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option
    ),
    bhhh = mxl_bhhh_parallel(
      theta = theta_hat,
      X = input_data$X,
      W = input_data$W,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      M = input_data$M,
      weights = input_data$weights,
      eta_draws = eta_draws,
      rc_dist = rc_dist,
      rc_correlation = rc_correlation,
      rc_mean = rc_mean,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option
    )
  )
  vcov_result <- invert_hessian(hess)
  }
  if (!is.null(vcov_result$vcov)) {
    rownames(vcov_result$vcov) <- param_names
    colnames(vcov_result$vcov) <- param_names
    names(vcov_result$se) <- param_names
  }

  # --- Back-transform to natural scale ----------------------------------------
  # Uses the bt_mult / bt_shift map built before optimization:
  #   theta_natural = bt_mult * theta_scaled + bt_shift
  #   vcov_natural  = (bt_mult bt_mult') o vcov_scaled  (shifts don't enter)
  if (scale_vars != "none") {
    bt <- .backtransform_estimates(theta_hat, vcov_result, bt_mult, bt_shift, param_names)
    theta_hat <- bt$theta
    vcov_result <- bt$vcov_result
    input_data$X <- natural_X
    input_data$W <- natural_W
  }

  # Reconstruct Sigma for display (from back-transformed L params if scaled)
  L_params <- theta_hat[param_map$sigma]
  sigma_mat <- build_var_mat(L_params, K_w, rc_correlation)
  w_names <- colnames(input_data$W)
  if (!is.null(w_names)) {
    rownames(sigma_mat) <- w_names
    colnames(sigma_mat) <- w_names
  }

  # Draws info (metadata only, not the full array)
  draws_info <- list(
    S = dim(eta_draws)[2],
    N = dim(eta_draws)[3],
    K_w = K_w
  )

  # Build S3 object
  new_choicer_mxl(
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
        W = input_data$W,
        alt_idx = input_data$alt_idx,
        choice_idx = input_data$choice_idx,
        M = input_data$M,
        weights = input_data$weights
      )
    },
    draws_info = draws_info,
    rc_dist = rc_dist,
    rc_correlation = rc_correlation,
    rc_mean = rc_mean,
    sigma = sigma_mat,
    se_method = se_method,
    scale_vars = scale_vars,
    sX = sX,
    sW = sW,
    choice_sampling = choice_sampling
  )
}


#' Prepare inputs for `mxl_loglik_gradient_parallel()`
#'
#' Prepares and validates inputs for mixed logit estimation routine.
#'
#' @param data Data frame containing choice data
#' @param id_col Name of the column identifying choice situations (individuals)
#' @param alt_col Name of the column identifying alternatives
#' @param choice_col Name of the column indicating chosen alternative (1 = chosen, 0 = not chosen)
#' @param covariate_cols Vector of names of columns to be used as covariates
#' @param random_var_cols Vector of names of columns to be used as random variables
#' @param weights Optional vector of weights for each choice situation. If NULL, equal weights are used.
#' @param weights_col Optional name of a column in \code{data} holding a per-row weight (constant within each choice situation). Mutually exclusive with \code{weights}.
#' @param outside_opt_label Label for the outside option (if any). If NULL, no outside option is assumed.
#' @param include_outside_option Logical indicating whether to include an outside option in the model.
#' @param rc_correlation Logical indicating whether random coefficients are correlated. Default is FALSE.
#' @returns A `choicer_data_mxl` object (list) containing:
#'   \itemize{
#'     \item `X`: Fixed-coefficient design matrix (sum(M) x K_x).
#'     \item `W`: Random-coefficient design matrix (sum(M) x K_w).
#'     \item `alt_idx`: Integer vector of alternative indices.
#'     \item `choice_idx`: Integer vector of chosen alternative indices.
#'     \item `M`: Integer vector with number of alternatives per choice situation.
#'     \item `N`: Number of choice situations.
#'     \item `weights`: Vector of weights.
#'     \item `include_outside_option`: Logical flag.
#'     \item `rc_correlation`: Logical flag.
#'     \item `alt_mapping`: data.table mapping alternatives to summary statistics.
#'     \item `dropped_cols`: Names of columns dropped due to collinearity, if any.
#'     \item `data_spec`: List with column-name metadata.
#'   }
#' @examples
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' input <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("w1", "w2"))
#' str(input$X)
#' str(input$W)
#' @export
prepare_mxl_data <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    random_var_cols,
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    rc_correlation = FALSE,
    weights_col = NULL
) {

  ## Preliminary housekeeping --------------------------------------------------
  dt <- data.table::as.data.table(data)[]

  # Check if all relevant variables are available
  needed <- c(id_col, alt_col, choice_col, covariate_cols, random_var_cols)
  if (!is.null(weights) && !is.null(weights_col)) {
    stop("Supply only one of `weights` or `weights_col`.")
  }
  if (!is.null(weights_col)) needed <- c(needed, weights_col)
  if (!all(needed %in% names(dt)))
    stop("Missing columns: ",
         paste(setdiff(needed, names(dt)), collapse = ", "))

  # Drop non-relevant variables
  vars_to_drop <- setdiff(names(dt), needed)
  if (length(vars_to_drop) > 0) {
    dt[, (vars_to_drop) := NULL]
  }

  ## Drop ids with missing observations ----------------------------------------
  dt[, HAS_NA := rowSums(is.na(.SD)) > 0]
  ids_to_drop <- dt[HAS_NA==TRUE, get(id_col)] |> unique()
  if (length(ids_to_drop) > 0) {
    dt <- dt[!(get(id_col) %in% ids_to_drop)]
    warning("Removed ", length(ids_to_drop),
            " choice situations containing missing values.")
  }
  if (nrow(dt) == 0) {
    stop("All choice situations removed due to missing values.")
  }
  dt[, HAS_NA := NULL]

  ## Sanity checks ---------------------------------------------------------

  ## covariates must be numeric
  if (!all(vapply(dt[, ..covariate_cols], is.numeric, logical(1L))))
    stop("All covariates must be numeric.")
  if (!all(vapply(dt[, ..random_var_cols], is.numeric, logical(1L))))
    stop("All covariates must be numeric.")

  ## choice column must be 0 or 1
  bad_choice <- dt[[choice_col]] %in% c(0, 1) == FALSE
  if (any(bad_choice))
    stop("`", choice_col, "` must contain only 0 and 1.")

  ## Exactly one '1' per choice situation
  by_id <- dt[, .(chosen = sum(get(choice_col))), by = id_col]
  if (include_outside_option == FALSE && any(by_id$chosen != 1)) {
    stop("Each ", id_col, " must have exactly one chosen alternative (one '1' in ",
         choice_col, ").")
  }
  if (include_outside_option && any(by_id$chosen > 1)) {
    stop("Each ", id_col, " must have at most one chosen alternative (one '1' in ",
         choice_col, "). An id with no explicit choice is assumed to be outside option.")
  }

  ## Create integer alternative codes ------------------------------------------

  if (!is.null(outside_opt_label) && include_outside_option==FALSE) {
    levels <- c(outside_opt_label, sort(setdiff(unique(dt[[alt_col]]), outside_opt_label)))
  } else {
    levels <- sort(unique(dt[[alt_col]]))
  }

  dt[, alt_int := as.integer(factor(get(alt_col), levels = levels))]

  ## Order rows ----------------------------------------------------------------
  ##   within each id: ascending alternative id
  ##   between ids   : ascending id
  data.table::setorderv(dt, c(id_col, "alt_int"))

  ## index of each row within its choice set
  dt[, idx_in_group := seq_len(.N), by = id_col]

  ## Build objects -------------------------------------------------------------
  ## design matrix
  X <- as.matrix(dt[, ..covariate_cols])
  X_res <- check_collinearity(X)
  X <- X_res$mat
  if (!is.null(X_res$dropped)) dropped_vars <- X_res$dropped # accumulate dropped vars if we had multiple checks


  W <- as.matrix(dt[, ..random_var_cols])
  W_res <- check_collinearity(W)
  W <- W_res$mat
  if (!is.null(W_res$dropped)) {
     if(exists("dropped_vars")) dropped_vars <- c(dropped_vars, W_res$dropped)
     else dropped_vars <- W_res$dropped
  }

  cols_to_drop <- union(covariate_cols, random_var_cols)
  dt[, (cols_to_drop) := NULL]

  ## alternative ids used for delta coefficients
  alt_idx <- as.integer(dt$alt_int)                             # length == sum(M)

  ## M[i] - # alternatives per choice situation
  M <- dt[, .N, by = id_col][["N"]]                             # length N

  ## N: number of individuals / choice situations
  ids <- dt[, get(id_col)][!duplicated(dt[[id_col]])]  # vector of ids in *current* order
  N   <- length(ids)

  ## Collapse a row-level weight column to one weight per choice situation.
  ## Done AFTER ordering/filtering so alignment is by id, never by position.
  if (!is.null(weights_col)) {
    if (!is.numeric(dt[[weights_col]])) {
      stop("`", weights_col, "` must be numeric.")
    }
    nuniq <- dt[, data.table::uniqueN(get(weights_col)), by = id_col][["V1"]]
    if (any(nuniq != 1L)) {
      stop("`", weights_col, "` must be constant within each '", id_col,
           "' (one weight per choice situation).")
    }
    wmap <- dt[, get(weights_col)[1L], by = id_col]
    weights <- wmap[["V1"]][match(ids, wmap[[id_col]])]
    if (any(!is.finite(weights))) {
      stop("`", weights_col, "` produced non-finite weights.")
    }
  }

  ## choice_idx[i] - 1-based index *within* the choice set data
  ## 0 == outside option (only if chosen = 0 for all inside options & include_outside_option == TRUE)
  if (include_outside_option) {
    # start with all-zero (everyone assumed to pick the outside good)
    choice_idx <- integer(N)
    chosen_dt <- dt[get(choice_col) == 1, .(pos = idx_in_group), by = id_col]

    # match chosen ids back to the master index vector
    data.table::setkeyv(chosen_dt, id_col)
    choice_idx[match(chosen_dt[[id_col]], ids)] <- chosen_dt$pos
  } else {
    # exactly one explicit choice per id
    choice_idx <- dt[get(choice_col) == 1, idx_in_group]
  }

  # Weights default = 1
  if (is.null(weights)) weights <- rep(1, N)

  ## Alternative summary -------------------------------------------------------

  if (include_outside_option) {
    inside_alt_mapping <- dt[
      , .(N_OBS = .N, N_CHOICES = sum(get(choice_col))),
      keyby = c("alt_int", alt_col)
    ]
    outside_alt_mapping <- data.table::data.table(alt_int=0L, N_OBS = N, N_CHOICES = sum(choice_idx == 0L))
    outside_alt_mapping[[alt_col]] <- outside_opt_label
    alt_mapping <- list(outside_alt_mapping, inside_alt_mapping) |>
      data.table::rbindlist(use.names = TRUE, fill = TRUE)
    data.table::setcolorder(alt_mapping, c("alt_int", alt_col, "N_OBS", "N_CHOICES"))
  } else {
    alt_mapping <- dt[
      , .(N_OBS = .N, N_CHOICES = sum(get(choice_col))),
      keyby = c("alt_int", alt_col)
    ]
  }

  alt_mapping[, `:=`(
    TAKE_RATE = N_CHOICES / N_OBS,
    MKT_SHARE = N_CHOICES / sum(N_CHOICES)
  )]

  ## Final validity checks -----------------------------------------------------
  stopifnot(
    length(alt_idx)    == nrow(X),
    length(choice_idx) == N,
    length(M)          == N,
    length(weights)    == N,
    all(is.finite(X)),
    all(is.finite(W))
  )

  ## Return output -------------------------------------------------------------
  structure(
    list(
      X           = X,
      W           = W,
      alt_idx     = alt_idx,
      choice_idx  = as.integer(choice_idx),
      M           = M,
      N           = N,
      weights     = weights,
      include_outside_option = include_outside_option,
      rc_correlation = rc_correlation,
      alt_mapping = alt_mapping[],
      dropped_cols = if(exists("dropped_vars")) dropped_vars else NULL,
      data_spec = list(
        id_col = id_col,
        alt_col = alt_col,
        choice_col = choice_col,
        covariate_cols = covariate_cols,
        random_var_cols = random_var_cols,
        outside_opt_label = outside_opt_label
      )
    ),
    class = "choicer_data_mxl"
  )
}

#' Halton draws for mixed logit
#'
#' Create halton normal draws in appropriate format for mixed logit estimation
#'
#' @param S Number of draws for each choice situation
#' @param N number of choice situations
#' @param K_w dimension of random coefficients (number of columns in W matrix)
#' @returns K_w x S x N array with halton standard normal draws
#' @examples
#' draws <- get_halton_normals(S = 50, N = 10, K_w = 2)
#' dim(draws)  # 2 x 50 x 10
#' @importFrom randtoolbox halton
#' @export
get_halton_normals <- function(S, N, K_w) {
  # Generate all needed Halton draws at once
  # We need S * N draws for each of K_w dimensions
  total_draws <- S * N

  # Generate Halton sequence
  # randtoolbox::halton returns a matrix of size total_draws x K_w
  # (but drops to vector when K_w = 1, so ensure matrix)
  halton_seq <- randtoolbox::halton(n = total_draws, dim = K_w, normal = TRUE)
  if (!is.matrix(halton_seq)) halton_seq <- matrix(halton_seq, ncol = 1)

  # Initialize the eta_draws array
  eta_draws <- array(0, dim = c(K_w, S, N))
  
  # Fill the array
  # The original code used: start_index = (i - 1) * S + 1 for each individual
  # This corresponds to taking chunks of S rows from the halton sequence
  for (i in 1:N) {
     start_row <- (i - 1) * S + 1
     end_row   <- i * S
     
     # halton_seq[start:end, ] is S x K_x
     # we want K_x x S for eta_draws[, , i]
     eta_draws[, , i] <- t(halton_seq[start_row:end_row, , drop=FALSE])
  }

  return(eta_draws)
}
