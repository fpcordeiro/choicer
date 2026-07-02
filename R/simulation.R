# Data-generating processes (DGPs) for discrete choice models.
# Returns classed `choicer_sim` objects consumed by recovery_table() and
# monte_carlo() in R/recovery.R and benchmark_fit() in _benchmarks/bench_helpers.R.

# Internal helpers =============================================================

.lse <- function(x) {
  a <- max(x)
  if (a == Inf) return(Inf)
  if (a == -Inf) return(-Inf)
  a + log(sum(exp(x - a)))
}

.draw_choice_set <- function(N, J, vary = TRUE) {
  if (!vary) {
    sizes <- rep(J, N)
    idx <- replicate(N, seq_len(J), simplify = FALSE)
    return(list(sizes = sizes, idx = idx))
  }
  sizes <- sample(2L:J, N, replace = TRUE)
  idx <- vector("list", N)
  for (i in seq_len(N)) idx[[i]] <- sort(sample(J, sizes[i]))
  list(sizes = sizes, idx = idx)
}

# choicer_sim S3 class =========================================================

#' Construct a `choicer_sim` object
#'
#' Wraps simulated data, true parameter values, and DGP settings into a
#' classed list. Returned by [simulate_mnl_data()], [simulate_mxl_data()],
#' and [simulate_nl_data()], and consumed by [recovery_table()].
#'
#' @param data A `data.table` of simulated choice observations.
#' @param true_params Named list of true DGP parameters
#'   (e.g. `beta`, `delta`, `Sigma`, `mu`, `lambdas`).
#' @param settings Named list of DGP settings (e.g. `N`, `J`, `K_x`).
#' @param model Character scalar: `"mnl"`, `"mxl"`, `"nl"`, `"mnp"`,
#'   `"hmnl"`, or `"hmnp"`.
#' @returns A list of class `choicer_sim`.
#' @export
new_choicer_sim <- function(data, true_params, settings, model) {
  if (!(is.data.frame(data) || data.table::is.data.table(data))) {
    stop("`data` must be a data.frame or data.table.")
  }
  if (!is.list(true_params)) stop("`true_params` must be a named list.")
  if (!is.list(settings))    stop("`settings` must be a named list.")
  if (!is.character(model) || length(model) != 1L ||
      !(model %in% c("mnl", "mxl", "nl", "mnp", "hmnl", "hmnp"))) {
    stop("`model` must be one of \"mnl\", \"mxl\", \"nl\", \"mnp\", ",
         "\"hmnl\", or \"hmnp\".")
  }
  structure(
    list(data = data, true_params = true_params, settings = settings, model = model),
    class = "choicer_sim"
  )
}

#' @export
print.choicer_sim <- function(x, ...) {
  cat("<choicer_sim: ", x$model, ">\n", sep = "")
  cat("  settings:\n")
  for (nm in names(x$settings)) {
    v <- x$settings[[nm]]
    if (is.list(v)) {
      parts <- vapply(
        v,
        function(z) paste0("(", paste(z, collapse = ","), ")"),
        character(1)
      )
      cat("    ", nm, " = ", paste(parts, collapse = ", "), "\n", sep = "")
    } else {
      cat("    ", nm, " = ", paste(v, collapse = " "), "\n", sep = "")
    }
  }
  cat("  rows in $data: ", nrow(x$data), "\n", sep = "")
  cat("  true_params: ", paste(names(x$true_params), collapse = ", "), "\n", sep = "")
  invisible(x)
}

# MNL DGP ======================================================================

#' Simulate multinomial logit data
#'
#' Generates synthetic choice data with i.i.d. Gumbel errors, optionally with
#' varying choice-set sizes and an outside option (alt = 0). Choices are
#' determined by argmax of utility; covariates are drawn as Uniform(-1, 1).
#'
#' @param N Number of choice situations.
#' @param J Number of inside alternatives.
#' @param beta Fixed coefficients for `x1..x{K_x}` (length `K_x = length(beta)`).
#' @param delta Alternative-specific constants for inside alternatives
#'   (length `J`). Defaults to an alternating pattern of `c(0.5, -0.5)`.
#' @param seed Random seed. Pass `NULL` to skip `set.seed()` (useful inside
#'   `monte_carlo()` where the caller manages RNG).
#' @param outside_option Logical; if `TRUE` (default) an outside option with
#'   `alt = 0` and zero covariates is added to every choice set.
#' @param vary_choice_set Logical; if `TRUE` (default) choice set size is
#'   sampled uniformly from `2:J`; if `FALSE` every individual faces all `J`
#'   inside alternatives.
#' @returns A `choicer_sim` object.
#' @examples
#' \donttest{
#' sim <- simulate_mnl_data(N = 1000, J = 5, seed = 123)
#' print(sim)
#' }
#' @export
simulate_mnl_data <- function(N = 5000,
                              J = 5,
                              beta = c(0.8, -0.6),
                              delta = NULL,
                              seed = 123,
                              outside_option = TRUE,
                              vary_choice_set = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  K_x <- length(beta)
  if (is.null(delta)) delta <- rep(c(0.5, -0.5), length.out = J)
  x_names <- paste0("x", seq_len(K_x))

  cs <- .draw_choice_set(N, J, vary = vary_choice_set)
  alt_indices <- cs$idx
  choice_set_sizes <- cs$sizes

  # Inside alternatives. Column insertion order: id, x1..x{K_x}, alt, delta_val.
  dt <- data.table::data.table(id = rep(seq_len(N), times = choice_set_sizes))
  for (k in seq_len(K_x)) dt[, (x_names[k]) := stats::runif(.N, -1, 1)]
  dt[, `:=`(
    alt       = alt_indices[[id]],
    delta_val = delta[alt_indices[[id]]]
  ), by = id]

  if (outside_option) {
    dt_outside <- dt[, .(alt = 0L), keyby = id]
    dt <- data.table::rbindlist(list(dt_outside, dt), use.names = TRUE, fill = TRUE)
    fill_cols <- setdiff(names(dt), c("id", "alt"))
    dt[alt == 0L, (fill_cols) := 0]
  }
  data.table::setkey(dt, id, alt)

  dt[, epsilon := -log(-log(stats::runif(.N)))]
  xb_expr <- paste(sprintf("%s * beta[%d]", x_names, seq_len(K_x)), collapse = " + ")
  util_expr <- parse(text = paste("delta_val +", xb_expr, "+ epsilon"))
  dt[, utility := eval(util_expr)]
  dt[, choice := data.table::fifelse(seq_len(.N) == which.max(utility), 1L, 0L), by = id]
  dt[, c("delta_val", "epsilon", "utility") := NULL]

  new_choicer_sim(
    data        = dt,
    true_params = list(beta = beta, delta = delta),
    settings    = list(N = N, J = J, K_x = K_x,
                       outside_option = outside_option,
                       vary_choice_set = vary_choice_set),
    model       = "mnl"
  )
}

# MXL DGP ======================================================================

#' Simulate mixed logit data
#'
#' Generates synthetic choice data with random coefficients drawn from a
#' multivariate normal (optionally log-normal per dimension) and an additional
#' mean shifter `mu`. Random coefficients are parameterized via the lower
#' Cholesky factor of `Sigma`. Covariates are Uniform(-1, 1) by default;
#' columns named in `price_cols` are drawn as `-Uniform(0.1, 3)` to mimic
#' strictly-negative price variables.
#'
#' @param N Number of choice situations.
#' @param J Number of inside alternatives.
#' @param beta Fixed coefficients for `x1..x{K_x}` (length `K_x = length(beta)`).
#' @param delta ASCs for inside alternatives (length `J`). Defaults to an
#'   alternating pattern of `c(0.5, -0.5)`.
#' @param mu Mean shifter for random coefficients (length `K_w = ncol(Sigma)`).
#'   Defaults to a zero vector.
#' @param Sigma Covariance matrix of random coefficients (square, `K_w x K_w`).
#' @param rc_dist Integer vector (length `K_w`): `0L` for normal, `1L` for
#'   log-normal. Default `NULL` is treated as all-normal.
#' @param rc_correlation Logical; if `NULL` (default) it is auto-detected from
#'   the off-diagonal entries of `Sigma`.
#' @param price_cols Character vector of `w*` column names to draw as
#'   `-Uniform(0.1, 3)` instead of `Uniform(-1, 1)`. Default `NULL`.
#' @param seed Random seed (`NULL` skips `set.seed()`).
#' @param outside_option Logical; include outside option with `alt = 0`.
#' @param vary_choice_set Logical; if `TRUE` (default) choice set size is
#'   sampled uniformly from `2:J`.
#' @returns A `choicer_sim` object. `true_params` includes `beta`, `delta`,
#'   `Sigma`, `L_params` (packed Cholesky parameters), `mu`, `rc_dist`,
#'   `rc_correlation`.
#' @details Random coefficients are constructed to match the estimator's
#'   parameterization in `src/mxlogit.cpp`. For every dimension the raw draw
#'   is `L %*% eta` where `eta ~ N(0, I)`. A normal random coefficient
#'   (`rc_dist = 0`) is then `gamma_k = mu_k + (L %*% eta)_k`. A log-normal
#'   random coefficient (`rc_dist = 1`) follows the shifted log-normal
#'   `beta_k = exp(mu_k) + exp((L %*% eta)_k)` -- not the textbook
#'   `exp(mu_k + sigma_k * eta)` -- so `mu_k` in `true_params$mu` is on the
#'   same scale the estimator recovers and `recovery_table()` can compare
#'   like-for-like.
#' @examples
#' \donttest{
#' sim <- simulate_mxl_data(N = 1000, J = 4, seed = 123)
#' print(sim)
#' }
#' @export
simulate_mxl_data <- function(N = 5000,
                              J = 4,
                              beta = c(0.8, -0.6),
                              delta = NULL,
                              mu = NULL,
                              Sigma = matrix(c(1.0, 0.5, 0.5, 1.5), nrow = 2),
                              rc_dist = NULL,
                              rc_correlation = NULL,
                              price_cols = NULL,
                              seed = 123,
                              outside_option = TRUE,
                              vary_choice_set = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  K_x <- length(beta)
  K_w <- ncol(Sigma)
  if (is.null(delta)) delta <- rep(c(0.5, -0.5), length.out = J)
  if (is.null(mu)) mu <- rep(0, K_w)
  if (is.null(rc_dist)) rc_dist <- rep(0L, K_w)
  if (is.null(rc_correlation)) {
    off <- Sigma[lower.tri(Sigma)]
    rc_correlation <- any(off != 0)
  }
  stopifnot(
    length(mu) == K_w,
    length(rc_dist) == K_w,
    all(rc_dist %in% c(0L, 1L))
  )
  x_names <- paste0("x", seq_len(K_x))
  w_names <- paste0("w", seq_len(K_w))
  if (!is.null(price_cols) && !all(price_cols %in% w_names)) {
    stop("`price_cols` must be a subset of ", paste(w_names, collapse = ", "))
  }

  # Random coefficients matching the estimator's parameterization in
  # src/mxlogit.cpp (see @details). For log-normal dimensions the DGP must
  # add exp(mu) (not mu) so recovery_table() compares like-for-like against
  # the recovered mu parameter.
  L_true <- t(chol(Sigma))
  gamma_i <- t(L_true %*% matrix(stats::rnorm(N * K_w), nrow = K_w, ncol = N))  # N x K_w
  for (r in seq_len(K_w)) {
    if (rc_dist[r] == 1L) gamma_i[, r] <- exp(gamma_i[, r])
  }
  if (any(mu != 0)) {
    for (r in seq_len(K_w)) {
      shift <- if (rc_dist[r] == 1L) exp(mu[r]) else mu[r]
      gamma_i[, r] <- gamma_i[, r] + shift
    }
  }

  # Packed Cholesky parameters matching build_var_mat() in C++:
  #   rc_correlation = TRUE  -> log-diag and free off-diag elements (row-major lower)
  #   rc_correlation = FALSE -> log-diag only
  L_params <- if (rc_correlation) {
    L_size <- K_w * (K_w + 1) / 2
    vals <- numeric(L_size)
    pos <- 1L
    for (row in seq_len(K_w)) {
      for (col in seq_len(row)) {
        vals[pos] <- if (row == col) log(L_true[row, col]) else L_true[row, col]
        pos <- pos + 1L
      }
    }
    vals
  } else {
    log(diag(L_true))
  }

  cs <- .draw_choice_set(N, J, vary = vary_choice_set)
  alt_indices <- cs$idx
  choice_set_sizes <- cs$sizes

  # Inside alternatives. Insertion order: id, w1..wKw, x1..xKx, gamma1..gammaKw, alt, delta_val.
  dt <- data.table::data.table(id = rep(seq_len(N), times = choice_set_sizes))
  for (k in seq_len(K_w)) {
    wk <- w_names[k]
    if (!is.null(price_cols) && wk %in% price_cols) {
      dt[, (wk) := -stats::runif(.N, 0.1, 3)]
    } else {
      dt[, (wk) := stats::runif(.N, -1, 1)]
    }
  }
  for (k in seq_len(K_x)) dt[, (x_names[k]) := stats::runif(.N, -1, 1)]
  for (k in seq_len(K_w)) dt[, (paste0("gamma", k)) := gamma_i[id, k]]
  dt[, `:=`(
    alt       = alt_indices[[id]],
    delta_val = delta[alt_indices[[id]]]
  ), by = id]

  if (outside_option) {
    dt_outside <- dt[, .(alt = 0L), keyby = id]
    dt <- data.table::rbindlist(list(dt_outside, dt), use.names = TRUE, fill = TRUE)
    fill_cols <- setdiff(names(dt), c("id", "alt"))
    dt[alt == 0L, (fill_cols) := 0]
  }
  data.table::setkey(dt, id, alt)

  dt[, epsilon := -log(-log(stats::runif(.N)))]
  xb_expr  <- paste(sprintf("%s * beta[%d]", x_names, seq_len(K_x)), collapse = " + ")
  wg_expr  <- paste(sprintf("%s * gamma%d", w_names, seq_len(K_w)), collapse = " + ")
  util_expr <- parse(text = paste("delta_val +", xb_expr, "+", wg_expr, "+ epsilon"))
  dt[, utility := eval(util_expr)]
  dt[, choice := data.table::fifelse(seq_len(.N) == which.max(utility), 1L, 0L), by = id]

  drop_cols <- c("delta_val", paste0("gamma", seq_len(K_w)), "epsilon", "utility")
  dt[, (drop_cols) := NULL]

  new_choicer_sim(
    data        = dt,
    true_params = list(
      beta = beta, delta = delta, Sigma = Sigma, L_params = L_params,
      mu = mu, rc_dist = rc_dist, rc_correlation = rc_correlation
    ),
    settings    = list(
      N = N, J = J, K_x = K_x, K_w = K_w,
      outside_option = outside_option,
      vary_choice_set = vary_choice_set
    ),
    model       = "mxl"
  )
}

# MNP DGP ======================================================================

#' Simulate multinomial probit data
#'
#' Generates synthetic choice data from the MNP data-generating process
#' estimated by [run_mnprobit()]: latent utility differences against the base
#' alternative (alternative 1),
#' \deqn{w_i = X_i \beta + \delta + \varepsilon_i, \qquad
#'       \varepsilon_i \sim N_{J-1}(0, \Sigma),}
#' with alternative \eqn{j > 1} chosen iff
#' \eqn{w_{ij} > \max(0, \max_{k \neq j} w_{ik})} and the base chosen iff all
#' \eqn{w_{ij} < 0}. Covariates are Uniform(-1, 1). Choice sets are balanced
#' (every individual faces all `J` alternatives), as the MNP estimator
#' requires; there is no outside-option flag — model an outside good as a
#' zero-covariate base alternative instead.
#'
#' The MNP likelihood only identifies parameters up to scale, so
#' `true_params` is reported on the *identified* scale (normalized by
#' \eqn{\sigma_{11}}): `beta` \eqn{= \beta / \sqrt{\sigma_{11}}}, `delta`
#' \eqn{= \delta / \sqrt{\sigma_{11}}}, and `Sigma` \eqn{= \Sigma /
#' \sigma_{11}} — the scale on which [run_mnprobit()] reports its posterior.
#' With the default `Sigma` (\eqn{\sigma_{11} = 1}) the DGP scale and the
#' identified scale coincide.
#'
#' @param N Number of choice situations.
#' @param J Number of alternatives (alternative 1 is the base).
#' @param beta Fixed coefficients for `x1..x{K_x}` (length `K_x = length(beta)`).
#' @param delta ASCs of the differenced utilities, one per non-base
#'   alternative (length `J - 1`). Defaults to an alternating pattern of
#'   `c(0.5, -0.5)`.
#' @param Sigma Covariance matrix of the differenced errors
#'   (`(J-1) x (J-1)`).
#' @param seed Random seed (`NULL` skips `set.seed()`).
#' @returns A `choicer_sim` object. `true_params` contains `beta`, `delta`,
#'   and `Sigma` on the identified scale (see Details).
#' @examples
#' \donttest{
#' sim <- simulate_mnp_data(N = 1000, J = 3, seed = 123)
#' print(sim)
#' }
#' @export
simulate_mnp_data <- function(N = 5000,
                              J = 3,
                              beta = c(0.8, -0.6),
                              delta = NULL,
                              Sigma = matrix(c(1.0, 0.5, 0.5, 1.5), nrow = 2),
                              seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  p <- J - 1L
  K_x <- length(beta)
  Sigma <- as.matrix(Sigma)
  if (nrow(Sigma) != p || ncol(Sigma) != p) {
    stop("`Sigma` must be a ", p, " x ", p, " matrix (J - 1 utility differences).")
  }
  if (is.null(delta)) delta <- rep(c(0.5, -0.5), length.out = p)
  if (length(delta) != p) stop("`delta` must have length J - 1 = ", p, ".")
  x_names <- paste0("x", seq_len(K_x))

  # Balanced choice sets: ids 1..N, alternatives 1..J.
  dt <- data.table::data.table(
    id  = rep(seq_len(N), each = J),
    alt = rep(seq_len(J), N)
  )
  for (k in seq_len(K_x)) dt[, (x_names[k]) := stats::runif(.N, -1, 1)]

  # Latent utility differences vs the base (alt 1): rows of W are the J - 1
  # non-base components; delta recycles down the rows (one ASC per component).
  V <- matrix(as.matrix(dt[, ..x_names]) %*% beta, nrow = J)        # J x N
  L <- t(chol(Sigma))
  eps <- L %*% matrix(stats::rnorm(p * N), p, N)
  W <- (V[-1L, , drop = FALSE] - rep(1, p) %o% V[1L, ]) + delta + eps

  chosen <- apply(W, 2, function(w) if (max(w) < 0) 1L else which.max(w) + 1L)
  dt[, choice := as.integer(alt == chosen[id])]

  s11 <- Sigma[1, 1]
  new_choicer_sim(
    data        = dt,
    true_params = list(
      beta  = beta / sqrt(s11),
      delta = delta / sqrt(s11),
      Sigma = Sigma / s11
    ),
    settings    = list(N = N, J = J, K_x = K_x, base_alt = 1L),
    model       = "mnp"
  )
}

# NL DGP =======================================================================

#' Simulate nested logit data
#'
#' Generates synthetic choice data with nested logit probabilities computed
#' analytically (log-sum-exp over inclusive values), then samples choices from
#' the implied multinomial. The outside option (`j = 0`) sits in a singleton
#' nest with `lambda = 1`.
#'
#' @param N Number of choice situations.
#' @param beta Fixed coefficients for covariates `X, W` (length 2 by default).
#' @param delta Named numeric vector of ASCs for inside alternatives.
#' @param nests List of integer vectors defining nest membership for inside
#'   alternatives.
#' @param lambdas Numeric vector of dissimilarity parameters, one per nest.
#' @param seed Random seed (`NULL` skips `set.seed()`).
#' @returns A `choicer_sim` object. `true_params` includes `beta`, `delta`,
#'   `lambdas`; `settings` includes the `nest_structure`. The returned
#'   `data` retains a `nest` column (integer, with `0L` for the outside
#'   option) for convenient use with [run_nestlogit()].
#' @note Unlike [simulate_mnl_data()] and [simulate_mxl_data()], this
#'   function does not expose `outside_option` or `vary_choice_set` flags.
#'   The outside option (`j = 0`) is always present as a singleton nest with
#'   `lambda = 1`, and every individual faces the full set of inside
#'   alternatives. Add these flags if downstream use cases need them.
#' @examples
#' \donttest{
#' sim <- simulate_nl_data(N = 2000, seed = 123)
#' print(sim)
#' }
#' @export
simulate_nl_data <- function(N       = 10000,
                             beta    = c(1.5, -0.8),
                             delta   = c("1" = 0.5, "2" = 0.3, "3" = -0.2,
                                         "4" = -0.5, "5" = 0.4),
                             nests   = list(c(1, 2), c(3, 4, 5)),
                             lambdas = c(0.8, 0.2),
                             seed    = 123) {
  if (!is.null(seed)) set.seed(seed)
  J_inside <- length(delta)

  dt <- data.table::CJ(id = 1:N, j = 0:J_inside)

  # Nest assignments: outside option sits in nest 0 with lambda = 1.
  dt[, nest := 0L]
  for (g in seq_along(nests)) {
    dt[j %in% nests[[g]], nest := as.integer(g)]
  }
  lambda_all <- c(1.0, lambdas)
  dt[, lambda := lambda_all[nest + 1]]

  delta_all <- c(0, delta)
  dt[, delta_val := delta_all[j + 1]]

  dt[, X := 0.0][j > 0, X := stats::rnorm(.N, mean = 1, sd = 1)]
  dt[, W := 0.0][j > 0, W := stats::runif(.N, -2, 2)]

  dt[, V := delta_val + beta[1] * X + beta[2] * W]
  dt[, V_over_lambda := V / lambda]

  dt_iv <- dt[, .(IV = .lse(V_over_lambda)), by = .(id, nest, lambda)]
  dt_iv <- unique(dt_iv, by = c("id", "nest"))
  dt_iv[, log_denom := .lse(lambda * IV), by = id]
  dt_iv[, nest_prob := exp(lambda * IV - log_denom)]

  dt[dt_iv[, .(id, nest, nest_prob, IV)],
     on = c("id", "nest"),
     `:=`(nest_prob = i.nest_prob, IV = i.IV)]
  dt[, cond_prob := exp(V_over_lambda - IV)]
  dt[, P_j := cond_prob * nest_prob]
  dt[j == 0, P_j := nest_prob]  # singleton nest: P(j|k) = 1

  dt[, choice := data.table::fifelse(j == sample(j, size = 1, prob = P_j), 1L, 0L),
     by = id]

  # Keep `nest` in the returned data (run_nestlogit() needs it); drop the
  # other working columns used only by the DGP.
  dt[, c("lambda", "delta_val", "V", "V_over_lambda",
         "IV", "nest_prob", "cond_prob", "P_j") := NULL]

  new_choicer_sim(
    data        = dt,
    true_params = list(beta = beta, delta = delta, lambdas = lambdas),
    settings    = list(N = N, J_inside = J_inside, nest_structure = nests),
    model       = "nl"
  )
}

# HB DGPs (HMNL / HMNP) ========================================================

# Shared panel DGP behind simulate_hmnl_data() / simulate_hmnp_data(). The
# two models share the utility structure
#   U_ijt = x_ijt' gamma_i + delta_j + eps,   U_iot = eps  (outside)
# with delta_j = z_j' theta + xi_j, xi_j ~ N(0, sigma_d^2) and
# beta_i ~ MVN(beta, W); they differ only in the shock (EV1 vs N(0, sigma^2))
# and in rc_dist (HMNL only: log-normal coordinates enter utility as
# exp(beta_ik)). The outside option carries its OWN utility shock — it is
# stochastic with systematic utility 0, matching both estimators' models.
.simulate_hb_panel <- function(N, T, J, beta, W, theta, sigma_d, Z, rc_dist,
                               include_outside, seed, vary_choice_set,
                               shock = c("ev1", "normal"), sigma = 1) {
  shock <- match.arg(shock)
  if (!is.null(seed)) set.seed(seed)
  stopifnot(N >= 1, T >= 1, J >= 1)
  if (vary_choice_set && J < 2) {
    stop("`vary_choice_set = TRUE` requires J >= 2.")
  }
  K_x <- length(beta)
  if (is.null(W)) W <- diag(0.5, K_x)
  W <- as.matrix(W)
  if (nrow(W) != K_x || ncol(W) != K_x) {
    stop("`W` must be a ", K_x, " x ", K_x, " matrix.")
  }
  if (is.null(rc_dist)) rc_dist <- rep(0L, K_x)
  rc_dist <- as.integer(rc_dist)
  stopifnot(length(rc_dist) == K_x, all(rc_dist %in% c(0L, 1L)))
  if (!is.finite(sigma_d) || sigma_d < 0) {
    stop("`sigma_d` must be a non-negative number.")
  }

  # Mean-function design: intercept always present (matching the prep's Z
  # convention); `Z` supplies the J x (P - 1) NON-intercept columns and
  # defaults to Uniform(-1, 1) draws so theta recovery is non-trivial.
  P <- length(theta)
  if (P < 1) stop("`theta` must at least contain the intercept entry.")
  if (is.null(Z)) {
    Z_extra <- if (P > 1) {
      matrix(stats::runif(J * (P - 1), -1, 1), J, P - 1)
    } else {
      matrix(0, J, 0)
    }
  } else {
    Z_extra <- as.matrix(Z)
    if (nrow(Z_extra) != J || ncol(Z_extra) != P - 1) {
      stop("`Z` must be a ", J, " x ", P - 1, " matrix ",
           "(alternative-level covariates excluding the intercept).")
    }
  }
  z_names <- if (P > 1) paste0("z", seq_len(P - 1)) else character(0)
  colnames(Z_extra) <- z_names
  Z_full <- cbind("(Intercept)" = rep(1, J), Z_extra)

  # Alternative-level effects: delta_j = z_j' theta + xi_j (realized once,
  # global across all respondents and tasks).
  xi <- stats::rnorm(J, 0, sigma_d)
  delta <- drop(Z_full %*% theta) + xi

  # Respondent-level tastes: beta_i ~ MVN(beta, W); log-normal coordinates
  # enter utility as exp(beta_ik) (hierarchy normal on the log/chain scale).
  L_true <- t(chol(W))
  beta_i <- t(L_true %*% matrix(stats::rnorm(N * K_x), nrow = K_x, ncol = N)) +
    matrix(beta, N, K_x, byrow = TRUE)                          # N x K_x
  gamma_i <- beta_i
  for (k in seq_len(K_x)) {
    if (rc_dist[k] == 1L) gamma_i[, k] <- exp(gamma_i[, k])
  }

  x_names <- paste0("x", seq_len(K_x))
  n_tasks <- N * T

  cs <- .draw_choice_set(n_tasks, J, vary = vary_choice_set)
  alt_indices <- cs$idx
  choice_set_sizes <- cs$sizes

  # Inside alternatives. Insertion order: pid, task, x1..xK, alt, delta_val,
  # z1..z{P-1}. Task ids are globally unique so id_col = "task" works with
  # or without person_col = "pid".
  dt <- data.table::data.table(
    pid  = rep((seq_len(n_tasks) - 1L) %/% T + 1L, times = choice_set_sizes),
    task = rep(seq_len(n_tasks), times = choice_set_sizes)
  )
  for (k in seq_len(K_x)) dt[, (x_names[k]) := stats::runif(.N, -1, 1)]
  dt[, `:=`(
    alt       = alt_indices[[task]],
    delta_val = delta[alt_indices[[task]]]
  ), by = task]
  for (p in seq_along(z_names)) dt[, (z_names[p]) := Z_extra[alt, p]]

  if (include_outside) {
    dt_outside <- dt[, .(pid = pid[1L], alt = 0L), keyby = task]
    dt <- data.table::rbindlist(list(dt_outside, dt),
                                use.names = TRUE, fill = TRUE)
    fill_cols <- setdiff(names(dt), c("pid", "task", "alt"))
    dt[alt == 0L, (fill_cols) := 0]
  }
  data.table::setkey(dt, pid, task, alt)

  # The outside rows (x = 0, delta = 0) receive their own shock: the outside
  # good is stochastic, not a deterministic utility-0 bound.
  if (shock == "ev1") {
    dt[, epsilon := -log(-log(stats::runif(.N)))]
  } else {
    dt[, epsilon := stats::rnorm(.N, 0, sigma)]
  }
  xmat <- as.matrix(dt[, ..x_names])
  gmat <- gamma_i[dt$pid, , drop = FALSE]
  dt[, utility := rowSums(xmat * gmat) + delta_val + epsilon]
  dt[, choice := data.table::fifelse(seq_len(.N) == which.max(utility), 1L, 0L),
     by = task]
  dt[, c("delta_val", "epsilon", "utility") := NULL]

  list(data = dt, delta = delta, xi = xi, Z = Z_full, W = W,
       rc_dist = rc_dist, K_x = K_x, P = P)
}

#' Simulate hierarchical multinomial logit data
#'
#' Generates synthetic panel choice data from the hierarchical (random
#' coefficients + alternative-level random effects) logit DGP: respondents
#' `i = 1..N` face `T` choice situations each, with utilities
#' \deqn{U_{ijt} = x_{ijt}'\gamma_i + \delta_j + \epsilon_{ijt}, \qquad
#'       U_{iot} = \epsilon_{iot},}
#' i.i.d. Gumbel shocks (including a shock on the outside option, whose
#' systematic utility is 0), \eqn{\beta_i \sim N(\beta, W)} with
#' \eqn{\gamma_{ik} = \beta_{ik}} or \eqn{\exp(\beta_{ik})} per `rc_dist`,
#' and \eqn{\delta_j = z_j'\theta + \xi_j} with
#' \eqn{\xi_j \sim N(0, \sigma_d^2)}. Covariates are Uniform(-1, 1); the
#' alternative-level covariates `z*` are constant within each alternative.
#'
#' Log-normal coordinates are reported on the *chain* (log) scale in
#' `true_params$beta` — the scale on which the estimator's hierarchy
#' operates — while entering utility as `exp(beta_ik)`.
#'
#' @param N Number of respondents.
#' @param T Number of choice situations per respondent.
#' @param J Number of inside alternatives.
#' @param beta Population means of the structural random coefficients
#'   (length `K_x = length(beta)`; chain scale for log-normal coordinates).
#' @param W Covariance of the random coefficients (`K_x` x `K_x`).
#'   Defaults to `diag(0.5, K_x)`.
#' @param theta Mean-function coefficients for
#'   \eqn{\delta_j = z_j'\theta + \xi_j}; the first entry is the intercept,
#'   entries `2..P` load on the alternative-level covariates.
#' @param sigma_d Standard deviation of the alternative-level effects
#'   \eqn{\xi_j}.
#' @param Z Optional `J x (length(theta) - 1)` matrix of alternative-level
#'   covariates (excluding the intercept). Default `NULL` draws them
#'   Uniform(-1, 1).
#' @param rc_dist Integer vector (length `K_x`): `0L` for normal, `1L` for
#'   log-normal coordinates. Default `NULL` is all-normal.
#' @param include_outside Logical; if `TRUE` (default) an outside option
#'   with `alt = 0`, zero covariates, and its own Gumbel shock is added to
#'   every choice set.
#' @param seed Random seed (`NULL` skips `set.seed()`).
#' @param vary_choice_set Logical; if `TRUE` choice-set size is sampled
#'   uniformly from `2:J` per task. Default `FALSE`.
#' @returns A `choicer_sim` object. `true_params` contains `beta`, `W`,
#'   `theta`, `sigma_d`, the realized `delta` and `xi` vectors, the full
#'   mean-function design `Z` (intercept first), and `rc_dist`.
#' @examples
#' \donttest{
#' sim <- simulate_hmnl_data(N = 100, T = 4, J = 4, seed = 123)
#' print(sim)
#' sim$true_params$delta
#' }
#' @export
simulate_hmnl_data <- function(N = 500,
                               T = 10,
                               J = 4,
                               beta = c(0.8, -0.6),
                               W = NULL,
                               theta = c(0.5, -0.4),
                               sigma_d = 0.5,
                               Z = NULL,
                               rc_dist = NULL,
                               include_outside = TRUE,
                               seed = 123,
                               vary_choice_set = FALSE) {
  sim <- .simulate_hb_panel(
    N = N, T = T, J = J, beta = beta, W = W, theta = theta,
    sigma_d = sigma_d, Z = Z, rc_dist = rc_dist,
    include_outside = include_outside, seed = seed,
    vary_choice_set = vary_choice_set, shock = "ev1"
  )
  new_choicer_sim(
    data        = sim$data,
    true_params = list(
      beta = beta, W = sim$W, theta = theta, sigma_d = sigma_d,
      delta = sim$delta, xi = sim$xi, Z = sim$Z, rc_dist = sim$rc_dist
    ),
    settings    = list(
      N = N, T = T, J = J, K_x = sim$K_x, P = sim$P,
      include_outside = include_outside,
      vary_choice_set = vary_choice_set
    ),
    model       = "hmnl"
  )
}

#' Simulate hierarchical multinomial probit data
#'
#' Generates synthetic panel choice data from the hierarchical probit DGP
#' with iid normal utility shocks:
#' \deqn{U_{ijt} = x_{ijt}'\beta_i + \delta_j + \epsilon_{ijt}, \qquad
#'       U_{iot} = \epsilon_{iot}, \qquad
#'       \epsilon \sim N(0, \sigma^2),}
#' choice by argmax within the task. The outside option is stochastic — it
#' carries its own \eqn{N(0, \sigma^2)} shock on top of systematic utility
#' 0, exactly as in the estimator. \eqn{\beta_i \sim N(\beta, W)} (normal
#' only) and \eqn{\delta_j = z_j'\theta + \xi_j},
#' \eqn{\xi_j \sim N(0, \sigma_d^2)}, as in [simulate_hmnl_data()].
#'
#' The iid-probit likelihood identifies parameters only up to the common
#' scale \eqn{\sigma}, so `true_params` is reported on the *identified*
#' scale: `beta` \eqn{= \beta/\sigma}, `W` \eqn{= W/\sigma^2}, `theta`
#' \eqn{= \theta/\sigma}, `sigma_d` \eqn{= \sigma_d/\sigma}, `delta`
#' \eqn{= \delta/\sigma}, `xi` \eqn{= \xi/\sigma}. With the default
#' `sigma = 1` the DGP scale and the identified scale coincide.
#'
#' @inheritParams simulate_hmnl_data
#' @param beta Population means of the structural random coefficients
#'   (length `K_x = length(beta)`); all coordinates are normal.
#' @param include_outside Logical; if `TRUE` (default) an outside option
#'   with `alt = 0`, zero covariates, and its own normal shock is added to
#'   every choice set.
#' @param sigma Standard deviation of the iid utility shocks (DGP scale).
#' @returns A `choicer_sim` object. `true_params` contains `beta`, `W`,
#'   `theta`, `sigma_d`, the realized `delta` and `xi`, and the full
#'   mean-function design `Z` — all on the identified scale (see Details).
#' @examples
#' \donttest{
#' sim <- simulate_hmnp_data(N = 100, T = 4, J = 4, seed = 123)
#' print(sim)
#' sim$true_params$delta
#' }
#' @export
simulate_hmnp_data <- function(N = 500,
                               T = 10,
                               J = 4,
                               beta = c(0.8, -0.6),
                               W = NULL,
                               theta = c(0.5, -0.4),
                               sigma_d = 0.5,
                               Z = NULL,
                               include_outside = TRUE,
                               seed = 123,
                               vary_choice_set = FALSE,
                               sigma = 1) {
  if (!is.finite(sigma) || sigma <= 0) {
    stop("`sigma` must be a positive number.")
  }
  sim <- .simulate_hb_panel(
    N = N, T = T, J = J, beta = beta, W = W, theta = theta,
    sigma_d = sigma_d, Z = Z, rc_dist = NULL,
    include_outside = include_outside, seed = seed,
    vary_choice_set = vary_choice_set, shock = "normal", sigma = sigma
  )
  new_choicer_sim(
    data        = sim$data,
    true_params = list(
      beta = beta / sigma, W = sim$W / sigma^2, theta = theta / sigma,
      sigma_d = sigma_d / sigma, delta = sim$delta / sigma,
      xi = sim$xi / sigma, Z = sim$Z
    ),
    settings    = list(
      N = N, T = T, J = J, K_x = sim$K_x, P = sim$P,
      include_outside = include_outside,
      vary_choice_set = vary_choice_set, sigma = sigma
    ),
    model       = "hmnp"
  )
}
