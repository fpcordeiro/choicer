
#' Runs Bayesian multinomial probit estimation
#'
#' Estimates a multinomial probit model by Gibbs sampling with data
#' augmentation (Albert & Chib 1993; McCulloch & Rossi 1994). The model is
#' specified in utility differences against a base alternative: for choice
#' situation \eqn{i} with \eqn{J} alternatives, \eqn{w_i = X_i \beta +
#' \epsilon_i} with \eqn{\epsilon_i \sim N_{J-1}(0, \Sigma)}.
#'
#' Two workflows are supported:
#' \describe{
#'   \item{Convenience (default)}{Supply \code{data} and column names. Data
#'     preparation (\code{\link{prepare_mnp_data}}) is handled automatically.}
#'   \item{Advanced}{Call \code{\link{prepare_mnp_data}} yourself and pass the
#'     result via \code{input_data}.}
#' }
#'
#' \strong{Identification.} The multinomial probit likelihood is invariant to
#' a common rescaling \eqn{(\beta, \Sigma) \to (c\beta, c^2\Sigma)}. The
#' sampler runs on the non-identified parameterization (unrestricted
#' \eqn{\Sigma} with an inverse-Wishart prior) and identified quantities are
#' computed by normalizing each kept draw by \eqn{\sigma_{11}}:
#' \eqn{\beta / \sqrt{\sigma_{11}}} and \eqn{\Sigma / \sigma_{11}}. This is
#' the McCulloch & Rossi (1994) default, which keeps all Gibbs conditionals
#' conjugate and mixes better than the fully identified sampler of McCulloch,
#' Polson & Rossi (2000). Reported coefficients, standard deviations, and
#' credible intervals are posterior summaries of the identified draws.
#'
#' \strong{Reproducibility.} The sampler uses its own thread-safe RNG with
#' one stream per (iteration, observation), so results are reproducible
#' independent of the number of OpenMP threads (see
#' \code{set_num_threads()}). When \code{mcmc$seed} is not supplied, a
#' master seed is drawn from R's RNG, so \code{set.seed()} controls the run.
#'
#' \strong{Scope.} Balanced choice sets are required: every choice situation
#' must contain the same \eqn{J} alternatives. To model an outside option,
#' include it as explicit rows with zero covariates and set \code{base_alt}
#' to its label.
#'
#' @param data Data frame containing choice data (convenience workflow).
#'   Mutually exclusive with \code{input_data}.
#' @param id_col Name of the column identifying choice situations (individuals).
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating chosen alternative (1 = chosen, 0 = not chosen).
#' @param covariate_cols Vector of names of columns to be used as covariates.
#' @param input_data List output from \code{\link{prepare_mnp_data}} (advanced
#'   workflow). Mutually exclusive with \code{data}.
#' @param base_alt Label of the base (reference) alternative used for utility
#'   differencing. If \code{NULL} (default), the first alternative in sort
#'   order is used.
#' @param use_asc Logical indicating whether to include alternative-specific
#'   constants (one intercept per non-base alternative in the differenced
#'   utilities).
#' @param prior Named list of prior settings, merged over defaults:
#'   \describe{
#'     \item{\code{beta_bar}}{Prior mean of \eqn{\beta} (default \code{rep(0, K)}).}
#'     \item{\code{A}}{Prior precision of \eqn{\beta} (default \code{0.01 * diag(K)}).}
#'     \item{\code{nu}}{Inverse-Wishart degrees of freedom (default \code{p + 3}).}
#'     \item{\code{V}}{Inverse-Wishart scale matrix (default \code{nu * diag(p)}).}
#'   }
#' @param mcmc Named list of MCMC settings, merged over defaults:
#'   \describe{
#'     \item{\code{R}}{Total Gibbs iterations (default \code{10000}).}
#'     \item{\code{burn}}{Burn-in iterations discarded (default \code{floor(R / 5)}).}
#'     \item{\code{thin}}{Keep every \code{thin}-th post-burn-in draw (default \code{1}).}
#'     \item{\code{seed}}{Master RNG seed (default: drawn from R's RNG).}
#'     \item{\code{trace}}{Print progress every \code{trace} iterations (default \code{0}, silent).}
#'   }
#' @param keep_data Logical. If \code{TRUE} (default), stores prepared data in
#'   the returned object.
#' @returns A \code{choicer_mnp} object. S3 methods available:
#'   \code{summary()}, \code{coef()} (posterior means of identified
#'   coefficients), \code{vcov()} (posterior covariance of identified
#'   coefficient draws), \code{nobs()}. Posterior draws are stored in
#'   \code{$draws} (\code{beta} / \code{sigma} on the identified scale,
#'   \code{beta_raw} / \code{sigma_raw} unnormalized).
#' @references
#' Albert, J. H., & Chib, S. (1993). Bayesian Analysis of Binary and
#' Polychotomous Response Data. \emph{Journal of the American Statistical
#' Association}, 88(422), 669-679.
#'
#' McCulloch, R., & Rossi, P. E. (1994). An exact likelihood analysis of the
#' multinomial probit model. \emph{Journal of Econometrics}, 64(1-2), 207-240.
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 200; J <- 3; beta_true <- c(1.0, -0.5)
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, U := drop(as.matrix(.SD) %*% beta_true) + rnorm(.N), .SDcols = c("x1", "x2")]
#' dt[, choice := as.integer(U == max(U)), by = id]
#'
#' fit <- run_mnprobit(dt, "id", "alt", "choice", c("x1", "x2"),
#'                     mcmc = list(R = 500, burn = 100))
#' summary(fit)
#' coef(fit)
#' }
#' @export
run_mnprobit <- function(
    data = NULL,
    id_col = NULL,
    alt_col = NULL,
    choice_col = NULL,
    covariate_cols = NULL,
    input_data = NULL,
    base_alt = NULL,
    use_asc = TRUE,
    prior = list(),
    mcmc = list(),
    keep_data = TRUE
) {
  cl <- match.call()

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
    if (is.null(id_col) || is.null(alt_col) || is.null(choice_col) ||
        is.null(covariate_cols)) {
      stop("Convenience workflow requires: id_col, alt_col, choice_col, ",
           "and covariate_cols.")
    }
    input_list <- prepare_mnp_data(
      data, id_col, alt_col, choice_col, covariate_cols,
      base_alt = base_alt, use_asc = use_asc
    )
  } else {
    input_list <- input_data
  }

  # Resolve alt_col for naming (needed by both workflows)
  if (is.null(alt_col)) alt_col <- names(input_list$alt_mapping)[2]

  X <- input_list$X
  K <- ncol(X)
  p <- input_list$p
  N <- input_list$N

  # --- Priors (bayesm-like defaults) -------------------------------------------
  beta_bar <- prior$beta_bar %||% rep(0, K)
  A <- prior$A %||% (0.01 * diag(K))
  nu <- prior$nu %||% (p + 3)
  V <- prior$V %||% (nu * diag(p))

  if (length(beta_bar) != K) {
    stop("prior$beta_bar must have length ", K, ".")
  }
  if (!is.matrix(A) || any(dim(A) != K)) {
    stop("prior$A must be a ", K, " x ", K, " matrix.")
  }
  if (nu < p) {
    stop("prior$nu must be >= p = ", p, ".")
  }
  if (!is.matrix(V) || any(dim(V) != p)) {
    stop("prior$V must be a ", p, " x ", p, " matrix.")
  }

  # --- MCMC settings ------------------------------------------------------------
  R <- as.integer(mcmc$R %||% 10000L)
  burn <- as.integer(mcmc$burn %||% (R %/% 5L))
  thin <- as.integer(mcmc$thin %||% 1L)
  trace <- as.integer(mcmc$trace %||% 0L)
  seed <- as.numeric(mcmc$seed %||% sample.int(.Machine$integer.max, 1L))

  if (R <= burn) stop("mcmc$R must be greater than mcmc$burn.")
  if (thin < 1L) stop("mcmc$thin must be >= 1.")

  # --- Run the Gibbs sampler ------------------------------------------------------
  elapsed <- system.time({
    out <- mnp_gibbs(
      X = X, y = input_list$y, p = p,
      beta_bar = beta_bar, A = A, nu = nu, V = V,
      R = R, burn = burn, thin = thin,
      seed = seed, trace = trace
    )
  })

  message("MCMC run time ", convertTime(elapsed))

  # --- Post-process: identified quantities, normalized per draw -----------------
  # E[beta / sqrt(sigma_11)] != E[beta] / sqrt(E[sigma_11]); all identified
  # summaries are computed on the per-draw normalized draws.
  param_names <- colnames(X)
  sigma_names <- character(p * (p + 1L) / 2L)
  k <- 1L
  for (i in seq_len(p)) {
    for (j in seq_len(i)) {
      sigma_names[k] <- sprintf("Sigma_%d%d", i, j)
      k <- k + 1L
    }
  }

  betadraw <- out$betadraw
  sigmadraw <- out$sigmadraw
  sigma11 <- sigmadraw[, 1L]
  beta_id <- betadraw / sqrt(sigma11)
  sigma_id <- sigmadraw / sigma11
  colnames(betadraw) <- param_names
  colnames(beta_id) <- param_names
  colnames(sigmadraw) <- sigma_names
  colnames(sigma_id) <- sigma_names

  coefficients <- colMeans(beta_id)
  se <- apply(beta_id, 2, stats::sd)
  vcov <- stats::cov(beta_id)

  # Posterior-mean identified Sigma reconstructed as a p x p matrix
  sigma_mean <- colMeans(sigma_id)
  Sigma_mat <- matrix(0, p, p)
  k <- 1L
  for (i in seq_len(p)) {
    for (j in seq_len(i)) {
      Sigma_mat[i, j] <- sigma_mean[k]
      Sigma_mat[j, i] <- sigma_mean[k]
      k <- k + 1L
    }
  }
  w_names <- paste0("w_", input_list$alt_mapping[2:input_list$J][[alt_col]])
  dimnames(Sigma_mat) <- list(w_names, w_names)

  # --- Build S3 object -----------------------------------------------------------
  new_choicer_mnp(
    call = cl,
    coefficients = coefficients,
    se = se,
    vcov = vcov,
    sigma = Sigma_mat,
    draws = list(
      beta = beta_id,
      sigma = sigma_id,
      beta_raw = betadraw,
      sigma_raw = sigmadraw
    ),
    prior = list(beta_bar = beta_bar, A = A, nu = nu, V = V),
    mcmc = list(R = R, burn = burn, thin = thin, seed = seed,
                R_keep = out$R_keep),
    nobs = N,
    n_params = K,
    data_spec = input_list$data_spec,
    alt_mapping = input_list$alt_mapping,
    base_alt = input_list$base_alt,
    param_map = input_list$param_map,
    use_asc = input_list$use_asc,
    sampler = list(
      name = "gibbs_mr1994",
      elapsed_time = elapsed[["elapsed"]]
    ),
    data = if (keep_data) {
      list(X = X, y = input_list$y, p = p)
    }
  )
}

#' Prepare inputs for `mnp_gibbs()`
#'
#' Prepares and validates inputs for Bayesian multinomial probit estimation.
#' Covariates are differenced against the base alternative, so the design
#' matrix has one row per (choice situation, non-base alternative) pair.
#' Balanced choice sets are required: every choice situation must contain
#' the same \eqn{J} alternatives.
#'
#' @param data Data frame containing choice data.
#' @param id_col Name of the column identifying choice situations (individuals).
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating chosen alternative (1 = chosen, 0 = not chosen).
#' @param covariate_cols Vector of names of columns to be used as covariates.
#' @param base_alt Label of the base (reference) alternative used for utility
#'   differencing. If \code{NULL} (default), the first alternative in sort
#'   order is used.
#' @param use_asc Logical indicating whether to include alternative-specific
#'   constants (one intercept per non-base alternative).
#' @returns A list containing:
#'   \itemize{
#'     \item `X`: Stacked differenced design matrix ((N * p) x K), covariate
#'       columns first, then ASC columns when `use_asc = TRUE`.
#'     \item `y`: Integer vector of choices (0 = base alternative, j in 1..p
#'       for the j-th non-base alternative), one per choice situation.
#'     \item `p`: Number of utility differences (J - 1).
#'     \item `J`: Number of alternatives.
#'     \item `N`: Number of choice situations.
#'     \item `K`: Number of columns of `X`.
#'     \item `alt_mapping`: Data.table mapping alternatives to summary statistics
#'       (the base alternative is `alt_int = 1`).
#'     \item `base_alt`: Resolved label of the base alternative.
#'     \item `param_map`: Named list of integer index vectors (beta, asc).
#'     \item `use_asc`: Logical flag.
#'     \item `dropped_cols`: Names of columns dropped due to collinearity, if any.
#'     \item `data_spec`: List with column name metadata.
#'   }
#' @examples
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), x2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' input <- prepare_mnp_data(dt, "id", "alt", "choice", c("x1", "x2"))
#' str(input$X)
#' input$alt_mapping
#' @export
prepare_mnp_data <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    base_alt = NULL,
    use_asc = TRUE
) {
  ## Preliminary housekeeping --------------------------------------------------
  dt <- data.table::as.data.table(data)[]

  # Check if all relevant variables are available
  needed <- c(id_col, alt_col, choice_col, covariate_cols)
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
  ids_to_drop <- dt[HAS_NA == TRUE, get(id_col)] |> unique()
  if (length(ids_to_drop) > 0) {
    dt <- dt[!(get(id_col) %in% ids_to_drop)]
    warning("Removed ", length(ids_to_drop),
            " choice situations containing missing values.")
  }
  if (nrow(dt) == 0) {
    stop("All choice situations removed due to missing values.")
  }
  dt[, HAS_NA := NULL]

  ## Sanity checks -------------------------------------------------------------

  ## Covariates must be numeric
  if (!all(vapply(dt[, ..covariate_cols], is.numeric, logical(1L))))
    stop("All covariates must be numeric.")

  ## choice column must be 0/1 and exactly one '1' per choice situation
  bad_choice <- dt[[choice_col]] %in% c(0, 1) == FALSE
  if (any(bad_choice))
    stop("`", choice_col, "` must contain only 0 and 1.")

  by_id <- dt[, .(chosen = sum(get(choice_col))), by = id_col]
  if (any(by_id$chosen != 1)) {
    stop("Each ", id_col, " must have exactly one chosen alternative (one '1' in ",
         choice_col, ").")
  }

  ## Balanced choice sets (required by the MNP differencing) -------------------
  alts <- unique(dt[[alt_col]])
  J <- length(alts)
  if (J < 2) stop("Need at least 2 alternatives.")
  counts <- dt[, .N, by = id_col][["N"]]
  n_pairs <- nrow(unique(dt[, c(id_col, alt_col), with = FALSE]))
  if (any(counts != J) || n_pairs != nrow(dt)) {
    stop("MNP requires balanced choice sets: every choice situation must ",
         "contain the same ", J, " alternatives.")
  }

  ## Create integer alternative codes (base alternative first) -----------------
  if (!is.null(base_alt)) {
    if (!(base_alt %in% alts)) {
      stop("base_alt '", base_alt, "' is not one of the alternatives.")
    }
    levels <- c(base_alt, sort(setdiff(alts, base_alt)))
  } else {
    levels <- sort(alts)
  }

  dt[, alt_int := as.integer(factor(get(alt_col), levels = levels))]

  ## Order rows ----------------------------------------------------------------
  ##   within each id: ascending alternative id (base first)
  ##   between ids   : ascending id
  data.table::setorderv(dt, c(id_col, "alt_int"))

  N <- as.integer(nrow(dt) / J)
  p <- J - 1L

  ## Build objects -------------------------------------------------------------
  ## Differenced design matrix: row (i, j) is X_ij - X_i,base
  diff_dt <- dt[, lapply(.SD, function(v) v[-1L] - v[1L]),
                by = id_col, .SDcols = covariate_cols]
  diff_dt[, (id_col) := NULL]
  X_diff <- as.matrix(diff_dt)

  ## Alternatives summary (base alternative is alt_int = 1)
  alt_mapping <- dt[
    , .(N_OBS = .N, N_CHOICES = sum(get(choice_col))),
    keyby = c("alt_int", alt_col)
  ]
  alt_mapping[, `:=`(
    TAKE_RATE = N_CHOICES / N_OBS,
    MKT_SHARE = N_CHOICES / sum(N_CHOICES)
  )]

  ## ASC block: in differenced utilities the ASCs are J - 1 intercepts, one
  ## per non-base alternative
  if (use_asc) {
    X_asc <- diag(1, p)[rep(seq_len(p), N), , drop = FALSE]
    colnames(X_asc) <- paste0("ASC_", alt_mapping[2:J][[alt_col]])
    X <- cbind(X_diff, X_asc)
  } else {
    X <- X_diff
  }

  ## Collinearity check on the combined design: alternative-invariant
  ## covariates difference into the span of the ASC columns
  X_res <- check_collinearity(X)
  X <- X_res$mat

  ## Parameter index map (robust to collinearity drops)
  cov_cols_kept <- intersect(colnames(X), colnames(X_diff))
  asc_cols_kept <- setdiff(colnames(X), cov_cols_kept)
  param_map <- list(beta = match(cov_cols_kept, colnames(X)))
  if (length(asc_cols_kept) > 0) {
    param_map$asc <- match(asc_cols_kept, colnames(X))
  }

  ## y: 0 = base alternative, j in 1..p for the j-th non-base alternative
  y <- as.integer(dt[get(choice_col) == 1, alt_int] - 1L)

  ## Final validity checks -----------------------------------------------------
  stopifnot(
    nrow(X) == N * p,
    length(y) == N,
    all(y >= 0L & y <= p),
    all(is.finite(X))
  )

  ## return output -------------------------------------------------------------
  structure(
    list(
      X = X,
      y = y,
      p = p,
      J = as.integer(J),
      N = N,
      K = ncol(X),
      alt_mapping = alt_mapping,
      base_alt = alt_mapping[1][[alt_col]],
      param_map = param_map,
      use_asc = use_asc,
      dropped_cols = if (length(X_res$dropped) > 0) X_res$dropped else NULL,
      data_spec = list(
        id_col = id_col,
        alt_col = alt_col,
        choice_col = choice_col,
        covariate_cols = covariate_cols,
        base_alt = alt_mapping[1][[alt_col]]
      )
    ),
    class = "choicer_data_mnp"
  )
}
