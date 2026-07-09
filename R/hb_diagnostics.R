# Convergence diagnostics for the hierarchical Bayes models (choicer_hmnl /
# choicer_hmnp): split-R-hat (rhat), rank-normalized ESS/MCSE (ess, mcse),
# a consolidated diagnostic table (.hb_diagnostic_table), traceplots
# (traceplot), and the posterior-predictive share check (ppc_shares).

# --- Shared rank-normalization / FFT helpers --------------------------------

#' Blom's rank-normal transform of a pooled column matrix
#'
#' @param pooled Numeric matrix (rows = pooled draws across all chains,
#'   columns = parameters).
#' @returns A matrix of the same shape, each column rank-normalized
#'   (average-rank ties, Blom's transform).
#' @noRd
.rank_normalize_pooled <- function(pooled) {
  n_rows <- nrow(pooled)
  n_cols <- ncol(pooled)
  z <- matrix(NA_real_, n_rows, n_cols, dimnames = dimnames(pooled))
  for (p in seq_len(n_cols)) {
    r <- rank(pooled[, p], ties.method = "average")
    z[, p] <- stats::qnorm((r - 3 / 8) / (n_rows - 2 * 3 / 8 + 1))
  }
  z
}

#' Split a pooled draw matrix back into per-chain matrices
#'
#' @param pooled Numeric matrix (rows = pooled draws, in chain-major
#'   `rbind()` order), columns = parameters.
#' @param n_iter Number of draws per chain.
#' @param n_chains Number of chains.
#' @returns A list of `n_chains` matrices, each `n_iter` rows.
#' @noRd
.split_pooled_by_chain <- function(pooled, n_iter, n_chains) {
  out <- vector("list", n_chains)
  idx <- 1L
  for (cc in seq_len(n_chains)) {
    out[[cc]] <- pooled[idx:(idx + n_iter - 1L), , drop = FALSE]
    idx <- idx + n_iter
  }
  out
}

#' Autocovariance of a single (half-)chain via FFT
#'
#' Computes the classical biased autocovariance
#' `gamma_hat[t] = (1/n) sum_{i=1}^{n-t} yc[i] yc[i+t]` for `t = 0, ..., n-1`
#' via zero-padded FFT (Wiener-Khinchin), matching the direct double-sum
#' formula to floating-point precision.
#'
#' @param y Numeric vector (one half-chain).
#' @returns Numeric vector of length `length(y)`, `gamma_hat[0..n-1]`
#'   (1-indexed in R: `out[1]` is lag 0).
#' @noRd
.fft_autocov <- function(y) {
  n <- length(y)
  yc <- y - mean(y)
  L <- stats::nextn(2L * n)
  padded <- c(yc, rep(0, L - n))
  Ff <- stats::fft(padded)
  S <- Ff * Conj(Ff)
  ac <- Re(stats::fft(S, inverse = TRUE))
  ac[seq_len(n)] / (L * n)
}

#' Rank-normalized ESS (Geyer initial monotone sequence) for one parameter
#'
#' Implements Steps A-E of the rank-normalized ESS algorithm (Vehtari,
#' Gelman, Simpson, Carpenter & Buerkner 2021) for a single parameter's
#' pooled draws, split-half across chains. Used both for `ess_bulk` (on the
#' raw draws) and internally by `ess_tail` (on 0/1 quantile-indicator
#' sequences) -- the rank-normalization in Step A is always applied fresh to
#' whatever is passed in.
#'
#' @param chain_vecs List of `C` numeric vectors (one chain's draws each),
#'   identical length `n_iter`.
#' @returns A single numeric ESS value, or `NA_real_` if any intermediate
#'   quantity is non-finite or the pooled variance is non-positive.
#' @noRd
.ess_scalar <- function(chain_vecs) {
  n_chains <- length(chain_vecs)
  n_iter <- length(chain_vecs[[1L]])
  pooled <- unlist(chain_vecs, use.names = FALSE)
  n_pooled <- length(pooled)

  r <- rank(pooled, ties.method = "average")
  z <- stats::qnorm((r - 3 / 8) / (n_pooled - 2 * 3 / 8 + 1))
  z_chains <- vector("list", n_chains)
  idx <- 1L
  for (cc in seq_len(n_chains)) {
    z_chains[[cc]] <- z[idx:(idx + n_iter - 1L)]
    idx <- idx + n_iter
  }

  half <- n_iter %/% 2L
  if (half < 2L) return(NA_real_)
  halves <- unlist(lapply(z_chains, function(v) {
    list(v[seq_len(half)], v[(length(v) - half + 1L):length(v)])
  }), recursive = FALSE)
  m_halves <- length(halves)
  n <- half

  gamma_mat <- vapply(halves, .fft_autocov, numeric(n))
  if (is.null(dim(gamma_mat))) gamma_mat <- matrix(gamma_mat, nrow = n)
  s2 <- gamma_mat[1L, ] * n / (n - 1L)
  W <- mean(s2)
  chain_means <- vapply(halves, mean, numeric(1L))
  B <- n * stats::var(chain_means)
  var_plus <- (n - 1L) / n * W + B / n
  if (!is.finite(var_plus) || var_plus <= 0 || !is.finite(W) || W <= 0) {
    return(NA_real_)
  }

  t_max <- n - 1L
  if (t_max < 1L) {
    ess_val <- m_halves * n
    return(if (is.finite(ess_val)) ess_val else NA_real_)
  }
  rho_hat <- numeric(t_max)
  for (t in seq_len(t_max)) {
    rho_hat[t] <- 1 - (W - mean(gamma_mat[t + 1L, ])) / var_plus
  }

  # Geyer's initial monotone sequence, INCLUDING the implicit lag-0
  # autocorrelation rho_0 = 1: pairs are Gamma_m = rho_{2m} + rho_{2m+1} for
  # m = 0, 1, 2, ... (0-indexed), so Gamma_0 = rho_0 + rho_1 = 1 + rho_1.
  # Omitting rho_0 here (pairing rho_1+rho_2, rho_3+rho_4, ...) silently
  # drops 1 from tau_hat and roughly doubles the reported ESS -- this is the
  # standard Stan/Vehtari-et-al. convention (tau = -1 + 2*sum(Gamma^mono),
  # which reduces to tau = 1 + 2*sum_{t>=1} rho_t once rho_0=1 is folded in).
  rho_full <- c(1, rho_hat)   # rho_full[i] is lag (i - 1), i = 1..(t_max + 1)
  n_lags <- length(rho_full)
  k_max <- n_lags %/% 2L
  if (k_max < 1L) {
    tau_hat <- 1
  } else {
    Pk <- numeric(k_max)
    for (m in seq_len(k_max)) {
      i1 <- 2L * (m - 1L) + 1L
      i2 <- 2L * (m - 1L) + 2L
      Pk[m] <- rho_full[i1] + rho_full[i2]
    }
    k0 <- which(Pk < 0)[1L]
    if (is.na(k0)) k0 <- k_max + 1L
    if (k0 == 1L) {
      tau_hat <- -1
    } else {
      P_keep <- Pk[seq_len(k0 - 1L)]
      P_mono <- cummin(P_keep)
      tau_hat <- -1 + 2 * sum(P_mono)
    }
  }
  ess_raw <- if (is.finite(tau_hat) && tau_hat > 0) (m_halves * n) / tau_hat else (m_halves * n)
  ess_val <- min(ess_raw, m_halves * n)
  if (!is.finite(ess_val)) return(NA_real_)
  ess_val
}

# --- Exported diagnostics ----------------------------------------------------

#' Split-\eqn{\widehat{R}} convergence diagnostic
#'
#' Computes the split-\eqn{\widehat{R}} (potential scale reduction factor) of
#' Gelman et al. for each column of a matrix of posterior draws. Every chain
#' is split in half, so the diagnostic detects non-stationarity within a
#' single chain as well as disagreement across chains; values near 1 indicate
#' convergence, and values above roughly 1.05 warrant a longer run.
#'
#' @param draws A matrix of posterior draws (rows = iterations, columns =
#'   parameters) for a single chain, or a list of such matrices (one per
#'   chain, identical dimensions).
#' @param rank Logical; if `FALSE` (default) computes the classic split
#'   \eqn{\widehat{R}} (unchanged from prior releases, bit-for-bit). If
#'   `TRUE`, computes the rank-normalized, folded \eqn{\widehat{R}} of
#'   Vehtari, Gelman, Simpson, Carpenter & Buerkner (2021, *Bayesian
#'   Analysis*): the max of a bulk (rank-normalized draws) and a fold (rank-
#'   normalized absolute deviation from the pooled median) split-\eqn{R-hat},
#'   which is more robust to heavy tails and scale differences across chains.
#' @returns Named numeric vector with one \eqn{\widehat{R}} per parameter
#'   (`NA` for parameters with zero variance).
#' @examples
#' set.seed(42)
#' draws <- matrix(rnorm(2000), ncol = 2,
#'                 dimnames = list(NULL, c("a", "b")))
#' rhat(draws)          # ~1: white noise is stationary
#' drifting <- cbind(a = cumsum(rnorm(1000)))
#' rhat(drifting)       # >> 1: a random walk is not
#' rhat(draws, rank = TRUE)   # rank-normalized variant
#' @export
rhat <- function(draws, rank = FALSE) {
  if (is.matrix(draws)) draws <- list(draws)
  if (!is.list(draws) || !all(vapply(draws, is.matrix, logical(1L)))) {
    stop("`draws` must be a matrix or a list of matrices.")
  }
  dims <- unique(lapply(draws, dim))
  if (length(dims) != 1L) {
    stop("All chains must have identical dimensions.")
  }
  n_iter <- nrow(draws[[1L]])
  if (n_iter < 4L) {
    stop("Need at least 4 draws per chain to compute split-R-hat.")
  }

  if (!isTRUE(rank)) {
    ## Split every chain in half; each half is one 'chain' in the classic
    ## between/within variance decomposition.
    half <- n_iter %/% 2L
    splits <- unlist(lapply(draws, function(m) {
      list(m[seq_len(half), , drop = FALSE],
           m[(nrow(m) - half + 1L):nrow(m), , drop = FALSE])
    }), recursive = FALSE)

    m_chains <- length(splits)
    chain_means <- vapply(splits, colMeans, numeric(ncol(draws[[1L]])))
    chain_vars <- vapply(splits, function(m) apply(m, 2L, stats::var),
                         numeric(ncol(draws[[1L]])))
    if (is.null(dim(chain_means))) {   # single-parameter edge: keep matrix form
      chain_means <- matrix(chain_means, nrow = 1L)
      chain_vars <- matrix(chain_vars, nrow = 1L)
    }

    B <- half * apply(chain_means, 1L, stats::var)   # between-half variance
    W <- rowMeans(chain_vars)                        # within-half variance
    var_plus <- (half - 1L) / half * W + B / half
    out <- sqrt(var_plus / W)
    out[!is.finite(out)] <- NA_real_
    names(out) <- colnames(draws[[1L]])
    return(out)
  }

  ## rank = TRUE: rank-normalized, folded R-hat (Vehtari et al. 2021).
  .split_rhat_raw <- function(chain_list) {
    half <- n_iter %/% 2L
    splits <- unlist(lapply(chain_list, function(m) {
      list(m[seq_len(half), , drop = FALSE],
           m[(nrow(m) - half + 1L):nrow(m), , drop = FALSE])
    }), recursive = FALSE)
    chain_means <- vapply(splits, colMeans, numeric(ncol(chain_list[[1L]])))
    chain_vars <- vapply(splits, function(m) apply(m, 2L, stats::var),
                         numeric(ncol(chain_list[[1L]])))
    if (is.null(dim(chain_means))) {
      chain_means <- matrix(chain_means, nrow = 1L)
      chain_vars <- matrix(chain_vars, nrow = 1L)
    }
    B <- half * apply(chain_means, 1L, stats::var)
    W <- rowMeans(chain_vars)
    var_plus <- (half - 1L) / half * W + B / half
    r <- sqrt(var_plus / W)
    r[!is.finite(r)] <- NA_real_
    r
  }

  pooled <- do.call(rbind, draws)
  z_pooled <- .rank_normalize_pooled(pooled)
  z_chains <- .split_pooled_by_chain(z_pooled, n_iter, length(draws))
  rhat_bulk <- .split_rhat_raw(z_chains)

  med <- apply(pooled, 2L, stats::median)
  zeta <- lapply(draws, function(m) abs(sweep(m, 2L, med, "-")))
  zeta_pooled <- do.call(rbind, zeta)
  z_zeta_pooled <- .rank_normalize_pooled(zeta_pooled)
  z_zeta_chains <- .split_pooled_by_chain(z_zeta_pooled, n_iter, length(draws))
  rhat_fold <- .split_rhat_raw(z_zeta_chains)

  out <- pmax(rhat_bulk, rhat_fold)
  names(out) <- colnames(draws[[1L]])
  out
}

#' Rank-normalized effective sample size (bulk and tail)
#'
#' Computes the rank-normalized bulk and tail effective sample size (ESS) of
#' Vehtari, Gelman, Simpson, Carpenter & Buerkner (2021, *Bayesian
#' Analysis*) for each column of a matrix of posterior draws. `ess_bulk`
#' measures the effective number of independent draws available for
#' estimating the posterior mean/median; `ess_tail` is the minimum of the
#' 5th- and 95th-percentile indicator ESS, relevant for tail-quantile /
#' credible-interval precision. Autocovariance is computed via
#' \code{\link[stats]{fft}} (Geyer's initial monotone sequence estimator of
#' the integrated autocorrelation time).
#'
#' @param draws A matrix of posterior draws (rows = iterations, columns =
#'   parameters) for a single chain, or a list of such matrices (one per
#'   chain, identical dimensions).
#' @returns A numeric matrix, one row per parameter, two columns
#'   (`"bulk"`, `"tail"`); `NA` for parameters with zero variance or a
#'   non-finite result.
#' @examples
#' set.seed(42)
#' draws <- matrix(rnorm(2000), ncol = 2,
#'                 dimnames = list(NULL, c("a", "b")))
#' ess(draws)
#' @export
ess <- function(draws) {
  if (is.matrix(draws)) draws <- list(draws)
  if (!is.list(draws) || !all(vapply(draws, is.matrix, logical(1L)))) {
    stop("`draws` must be a matrix or a list of matrices.")
  }
  dims <- unique(lapply(draws, dim))
  if (length(dims) != 1L) {
    stop("All chains must have identical dimensions.")
  }
  n_iter <- nrow(draws[[1L]])
  if (n_iter < 4L) {
    stop("Need at least 4 draws per chain to compute ESS.")
  }

  p <- ncol(draws[[1L]])
  cn <- colnames(draws[[1L]])
  bulk <- numeric(p)
  tail_ess <- numeric(p)
  for (j in seq_len(p)) {
    chain_vecs <- lapply(draws, function(m) m[, j])
    bulk[j] <- .ess_scalar(chain_vecs)

    x <- unlist(chain_vecs, use.names = FALSE)
    q_lo <- stats::quantile(x, 0.05, names = FALSE, type = 7)
    q_hi <- stats::quantile(x, 0.95, names = FALSE, type = 7)
    i_lo <- lapply(chain_vecs, function(v) as.numeric(v <= q_lo))
    i_hi <- lapply(chain_vecs, function(v) as.numeric(v <= q_hi))
    ess_lo <- .ess_scalar(i_lo)
    ess_hi <- .ess_scalar(i_hi)
    tail_ess[j] <- if (is.na(ess_lo) || is.na(ess_hi)) NA_real_ else min(ess_lo, ess_hi)
  }
  out <- cbind(bulk = bulk, tail = tail_ess)
  rownames(out) <- cn
  out
}

#' Monte Carlo standard error of posterior summaries
#'
#' Monte Carlo standard error (MCSE) of the posterior mean or median, using
#' the rank-normalized bulk/tail ESS from \code{\link{ess}}.
#' \code{kind = "mean"}: `MCSE = SD / sqrt(ess_bulk)`. \code{kind =
#' "median"}: `MCSE = sqrt(pi/2) * SD / sqrt(ess_tail)` (the closed-form
#' normal asymptotic-efficiency approximation; `2/pi` is the asymptotic
#' relative efficiency of the sample median vs. the mean under normality).
#'
#' @inheritParams ess
#' @param kind `"mean"` (default) or `"median"`.
#' @returns Named numeric vector, one value per parameter (`NA` when the
#'   underlying ESS or pooled SD is undefined).
#' @examples
#' set.seed(42)
#' draws <- matrix(rnorm(2000), ncol = 2,
#'                 dimnames = list(NULL, c("a", "b")))
#' mcse(draws)
#' mcse(draws, kind = "median")
#' @export
mcse <- function(draws, kind = c("mean", "median")) {
  kind <- match.arg(kind)
  if (is.matrix(draws)) draws <- list(draws)
  if (!is.list(draws) || !all(vapply(draws, is.matrix, logical(1L)))) {
    stop("`draws` must be a matrix or a list of matrices.")
  }
  dims <- unique(lapply(draws, dim))
  if (length(dims) != 1L) {
    stop("All chains must have identical dimensions.")
  }
  n_iter <- nrow(draws[[1L]])
  if (n_iter < 4L) {
    stop("Need at least 4 draws per chain to compute ESS.")
  }

  ess_mat <- ess(draws)
  p <- ncol(draws[[1L]])
  cn <- colnames(draws[[1L]])
  out <- rep(NA_real_, p)
  for (j in seq_len(p)) {
    x <- unlist(lapply(draws, function(m) m[, j]), use.names = FALSE)
    sd_pooled <- stats::sd(x)
    if (!is.finite(sd_pooled) || sd_pooled == 0) next
    e <- if (kind == "mean") ess_mat[j, "bulk"] else ess_mat[j, "tail"]
    if (is.na(e) || e <= 0) next
    out[j] <- if (kind == "mean") sd_pooled / sqrt(e) else sqrt(pi / 2) * sd_pooled / sqrt(e)
  }
  names(out) <- cn
  out
}

#' Build the R-hat table for a hierarchical Bayes fit
#'
#' Internal: assembles rank-normalized split-R-hat over the b, theta, and
#' sigma_d^2 draws of one or more chains into a named vector, used by
#' [run_hmnlogit()] (and `run_hmnprobit()`) for the convergence warning and
#' the `object$rhat` field.
#'
#' @param chain_list List with one element per chain, each a list holding
#'   matrices `b`, `theta` and vector `sigma_d2`.
#' @returns Named numeric vector of rank-normalized split-R-hat values.
#' @noRd
.hb_rhat_table <- function(chain_list) {
  bind_block <- function(name) {
    lapply(chain_list, function(ch) {
      m <- ch[[name]]
      if (is.null(dim(m))) m <- matrix(m, ncol = 1L,
                                       dimnames = list(NULL, name))
      m
    })
  }
  c(
    rhat(bind_block("b"), rank = TRUE),
    rhat(bind_block("theta"), rank = TRUE),
    rhat(bind_block("sigma_d2"), rank = TRUE)
  )
}

#' Broadened convergence-warning offender list
#'
#' Internal: computes rank-R-hat and ESS_bulk over all columns of the given
#' small blocks (b, theta, sigma_d2) plus every column of the delta block,
#' and returns the names/labels failing the Vehtari thresholds (rank-R-hat >
#' 1.01 or ESS_bulk < 400). Used by the [run_hmnlogit()] / `run_hmnprobit()`
#' fit-time convergence warning (this replaces the old bare `rhat_tab >
#' 1.05` check).
#'
#' `sigma2` (HMNP's raw, non-identified parameter-expansion scale) is
#' deliberately EXCLUDED from this gate even when present in
#' `chain_blocks_small`: by design it is not expected to converge in R-hat/
#' ESS terms (the chain is only meaningful up to scale, via the per-draw
#' normalization applied to `b`/`theta`/`delta`/`sigma_d2`), so flagging it
#' would be a false positive on every HMNP fit.
#'
#' @param chain_blocks_small List with one element per chain, each a named
#'   list of matrices/vectors for the small blocks (arbitrary block names,
#'   e.g. `b`, `theta`, `sigma_d2`, `sigma2`).
#' @param chain_list_delta List with one element per chain, each the
#'   J-column `delta` matrix.
#' @returns Character vector of unique offending parameter/alternative names
#'   (possibly empty).
#' @noRd
.hb_convergence_offenders <- function(chain_blocks_small, chain_list_delta) {
  block_names <- setdiff(names(chain_blocks_small[[1L]]), "sigma2")
  offenders <- character(0)
  for (bn in block_names) {
    mats <- lapply(chain_blocks_small, function(ch) {
      m <- ch[[bn]]
      if (is.null(dim(m))) m <- matrix(m, ncol = 1L, dimnames = list(NULL, bn))
      m
    })
    r <- tryCatch(rhat(mats, rank = TRUE), error = function(e) NULL)
    e <- tryCatch(ess(mats), error = function(e) NULL)
    if (!is.null(r)) offenders <- c(offenders, names(r)[r > 1.01])
    if (!is.null(e)) offenders <- c(offenders, rownames(e)[e[, "bulk"] < 400])
  }
  r_delta <- tryCatch(rhat(chain_list_delta, rank = TRUE), error = function(e) NULL)
  e_delta <- tryCatch(ess(chain_list_delta), error = function(e) NULL)
  if (!is.null(r_delta)) offenders <- c(offenders, names(r_delta)[r_delta > 1.01])
  if (!is.null(e_delta)) offenders <- c(offenders, rownames(e_delta)[e_delta[, "bulk"] < 400])
  unique(stats::na.omit(offenders))
}

#' Representative subset of `delta` alternatives for traceplot/diagnostics
#'
#' Deterministic auto-selection: the alternative(s) with the 3 highest
#' rank-R-hat plus the 3 lowest ESS_bulk (deduplicated, index order as
#' tiebreak), capped at 6 total.
#'
#' @param delta_list List of `C` `delta` matrices (J columns), one per chain.
#' @param top_n Number of R-hat / ESS extremes to take from each side
#'   (default 3).
#' @param cap Maximum number of columns returned (default 6).
#' @returns Integer vector of column indices into `delta_list[[1]]`.
#' @noRd
.hb_delta_representative <- function(delta_list, top_n = 3L, cap = 6L) {
  rhat_delta <- rhat(delta_list, rank = TRUE)
  ess_delta <- ess(delta_list)
  J <- length(rhat_delta)
  ord_rhat <- order(-rhat_delta, seq_len(J))
  ord_ess <- order(ess_delta[, "bulk"], seq_len(J))
  top_rhat <- utils::head(ord_rhat, top_n)
  bottom_ess <- utils::head(ord_ess, top_n)
  sel <- sort(unique(c(top_rhat, bottom_ess)))
  utils::head(sel, cap)
}

#' Build the consolidated diagnostic table for a hierarchical Bayes fit
#'
#' Internal: assembles one row per `b`/`theta` parameter, a single row for
#' `sigma_d^2` (and, for HMNP, `sigma^2 (raw)`), and a single worst-case
#' summary row for the (potentially ~200-column) `delta` block, using all
#' retained chains (`object$chains`) when available.
#'
#' @param object A `choicer_hmnl` or `choicer_hmnp` fit.
#' @returns A list: `table` (data.frame), `delta_worst_label` (character
#'   scalar), `chains` (integer), `R_keep` (integer).
#' @noRd
.hb_diagnostic_table <- function(object) {
  chain_list <- object$chains
  if (is.null(chain_list) || length(chain_list) < 2L) {
    chain_list <- list(object$draws)
  }
  n_chains <- length(chain_list)
  R_keep <- nrow(chain_list[[1L]]$b)

  as_matrix <- function(m, fallback_name) {
    if (is.null(dim(m))) m <- matrix(m, ncol = 1L,
                                     dimnames = list(NULL, fallback_name))
    m
  }

  build_rows <- function(mat_list, prefix) {
    Rhat <- rhat(mat_list, rank = TRUE)
    ess_mat <- ess(mat_list)
    mcse_mean <- mcse(mat_list, kind = "mean")
    cn <- colnames(mat_list[[1L]])
    data.frame(
      Block = paste0(prefix, "[", cn, "]"),
      Rhat = as.numeric(Rhat),
      ESS_bulk = as.numeric(ess_mat[, "bulk"]),
      ESS_tail = as.numeric(ess_mat[, "tail"]),
      MCSE_mean = as.numeric(mcse_mean),
      stringsAsFactors = FALSE
    )
  }

  build_single_row <- function(mat_list, label) {
    Rhat <- rhat(mat_list, rank = TRUE)
    ess_mat <- ess(mat_list)
    mcse_mean <- mcse(mat_list, kind = "mean")
    data.frame(
      Block = label,
      Rhat = as.numeric(Rhat),
      ESS_bulk = as.numeric(ess_mat[, "bulk"]),
      ESS_tail = as.numeric(ess_mat[, "tail"]),
      MCSE_mean = as.numeric(mcse_mean),
      stringsAsFactors = FALSE
    )
  }

  b_list <- lapply(chain_list, function(ch) ch$b)
  theta_list <- lapply(chain_list, function(ch) ch$theta)
  sigma_d2_list <- lapply(chain_list, function(ch) as_matrix(ch$sigma_d2, "sigma_d^2"))

  rows <- rbind(
    build_rows(b_list, "b"),
    build_rows(theta_list, "theta"),
    build_single_row(sigma_d2_list, "sigma_d^2")
  )

  if (!is.null(chain_list[[1L]]$sigma2)) {
    sigma2_list <- lapply(chain_list, function(ch) as_matrix(ch$sigma2, "sigma^2 (raw)"))
    rows <- rbind(rows, build_single_row(sigma2_list, "sigma^2 (raw)"))
  }

  delta_list <- lapply(chain_list, function(ch) ch$delta)
  J <- ncol(delta_list[[1L]])
  rhat_delta <- rhat(delta_list, rank = TRUE)
  ess_delta <- ess(delta_list)
  worst_idx <- which.max(rhat_delta)
  # which.max() on an all-NA vector (every delta column zero-variance /
  # non-finite) returns integer(0); fall back to the first column so the
  # delta summary row is never silently dropped by the rbind() below.
  if (length(worst_idx) == 0L) worst_idx <- 1L
  delta_worst_label <- colnames(delta_list[[1L]])[worst_idx]
  row_delta <- data.frame(
    Block = sprintf("delta (J=%d)", J),
    Rhat = as.numeric(rhat_delta[worst_idx]),
    ESS_bulk = as.numeric(ess_delta[worst_idx, "bulk"]),
    ESS_tail = as.numeric(ess_delta[worst_idx, "tail"]),
    MCSE_mean = NA_real_,
    stringsAsFactors = FALSE
  )
  rows <- rbind(rows, row_delta)
  rownames(rows) <- NULL

  list(
    table = rows,
    delta_worst_label = delta_worst_label,
    chains = n_chains,
    R_keep = R_keep
  )
}

#' Print the consolidated diagnostic table
#'
#' @param dt The list returned by `.hb_diagnostic_table()`.
#' @returns Invisible `NULL`.
#' @noRd
.print_hb_diagnostic_table <- function(dt) {
  tbl <- dt$table
  is_delta_row <- grepl("^delta \\(J=", tbl$Block)
  is_sigma2_row <- tbl$Block == "sigma^2 (raw)"

  fmt_rhat <- sprintf("%.3f", tbl$Rhat)
  fmt_bulk <- sprintf("%d", as.integer(round(tbl$ESS_bulk)))
  fmt_tail <- sprintf("%d", as.integer(round(tbl$ESS_tail)))
  fmt_mcse <- ifelse(is.na(tbl$MCSE_mean), "\u2014", sprintf("%.4f", tbl$MCSE_mean))
  fmt_rhat[is_delta_row] <- paste0(fmt_rhat[is_delta_row], "*")
  fmt_bulk[is_delta_row] <- paste0(fmt_bulk[is_delta_row], "*")
  fmt_tail[is_delta_row] <- paste0(fmt_tail[is_delta_row], "*")

  block_disp <- tbl$Block
  block_disp[is_sigma2_row] <- paste0(block_disp[is_sigma2_row], "^")

  block_w <- max(nchar(block_disp), nchar("Block"))
  rhat_w <- max(nchar(fmt_rhat), nchar("R-hat"))
  bulk_w <- max(nchar(fmt_bulk), nchar("ESS_bulk"))
  tail_w <- max(nchar(fmt_tail), nchar("ESS_tail"))
  mcse_w <- max(nchar(fmt_mcse), nchar("MCSE(mean)"))

  cat(sprintf("%-*s  %*s  %*s  %*s  %*s\n",
              block_w, "Block", rhat_w, "R-hat",
              bulk_w, "ESS_bulk", tail_w, "ESS_tail",
              mcse_w, "MCSE(mean)"))
  for (i in seq_len(nrow(tbl))) {
    cat(sprintf("%-*s  %*s  %*s  %*s  %*s\n",
                block_w, block_disp[i], rhat_w, fmt_rhat[i],
                bulk_w, fmt_bulk[i], tail_w, fmt_tail[i],
                mcse_w, fmt_mcse[i]))
  }
  if (any(is_delta_row)) {
    cat(sprintf("*worst: delta[%s]\n", dt$delta_worst_label))
  }
  if (any(is_sigma2_row)) {
    cat("^sigma^2 (raw) is the non-identified parameter-expansion scale ",
        "(expected to not converge by design; excluded from the ",
        "convergence-failure check).\n", sep = "")
  }
  invisible(NULL)
}

#' Traceplot for a hierarchical Bayes fit
#'
#' Generic dispatching on the fit's class. See
#' [traceplot.choicer_hb()] for the `choicer_hmnl` / `choicer_hmnp` method.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @returns The object, invisibly.
#' @export
traceplot <- function(object, ...) UseMethod("traceplot")

#' Traceplot method for hierarchical Bayes fits
#'
#' Overlaid per-chain traceplots (base `graphics`, no new dependency) for a
#' `choicer_hmnl` / `choicer_hmnp` fit's population coefficients (`b`), delta
#' mean-function coefficients (`theta`), alternative-effect variance
#' (`sigma_d2`), and, opt-in, a representative subset of the (potentially
#' ~200-column) alternative effects (`delta`).
#'
#' @param object A `choicer_hmnl` or `choicer_hmnp` fit.
#' @param block Character vector, any non-empty subset of `c("b", "theta",
#'   "sigma_d2", "delta")` (HMNP fits additionally accept `"sigma2"`, the raw
#'   non-identified shock-variance chain). Default plots `b`, `theta`,
#'   `sigma_d2` (not `delta`, which can have ~200 columns).
#' @param which Only consulted when `"delta"` is in `block`. `NULL`
#'   (default) auto-selects a representative subset (the 3 highest rank-R-hat
#'   plus 3 lowest ESS_bulk alternatives, deduplicated, capped at 6). If
#'   supplied, a character vector of alternative labels or an integer vector
#'   of column indices.
#' @param ... Additional arguments (ignored).
#' @returns `object`, invisibly.
#' @examples
#' \donttest{
#' sim <- simulate_hmnl_data(N = 60, T = 2, J = 4, seed = 1)
#' fit <- suppressWarnings(run_hmnlogit(sim$data, "task", "alt", "choice", c("x1", "x2"),
#'                     person_col = "pid",
#'                     mcmc = list(R = 300, burn = 100), chains = 2))
#' traceplot(fit)
#' }
#' @export
traceplot.choicer_hb <- function(object, block = c("b", "theta", "sigma_d2"),
                                  which = NULL, ...) {
  valid_blocks <- c("b", "theta", "sigma_d2", "delta")
  if (identical(object$model, "hmnp")) valid_blocks <- c(valid_blocks, "sigma2")
  if (!is.character(block) || length(block) < 1L || !all(block %in% valid_blocks)) {
    stop("`block` must be one or more of ", paste(valid_blocks, collapse = ", "), ".")
  }

  chain_list <- object$chains
  n_chains <- length(chain_list)
  if (is.null(chain_list) || n_chains < 1L) {
    message("`object$chains` is unavailable (fit predates this feature); ",
            "plotting chain 1 (object$draws) only.")
    chain_list <- list(object$draws)
    n_chains <- 1L
  }
  R_keep <- nrow(chain_list[[1L]]$b)

  panels <- list()
  for (blk in block) {
    if (blk == "delta") {
      delta_list <- lapply(chain_list, function(ch) ch$delta)
      sel <- which
      if (is.null(sel)) {
        sel <- .hb_delta_representative(delta_list)
      } else if (is.character(sel)) {
        idx <- match(sel, colnames(delta_list[[1L]]))
        if (any(is.na(idx))) {
          stop("`which` contains alternative label(s) not found in delta columns.")
        }
        sel <- idx
      } else if (is.numeric(sel)) {
        sel <- as.integer(sel)
        if (any(sel < 1L | sel > ncol(delta_list[[1L]]))) {
          stop("`which` integer index out of range for delta columns.")
        }
      } else {
        stop("`which` must be a character vector of alternative labels or ",
             "an integer vector of column indices.")
      }
      for (j in sel) {
        panels[[length(panels) + 1L]] <- list(
          label = paste0("delta[", colnames(delta_list[[1L]])[j], "]"),
          series = lapply(delta_list, function(m) m[, j])
        )
      }
    } else {
      mat_list <- lapply(chain_list, function(ch) {
        m <- ch[[blk]]
        if (is.null(dim(m))) m <- matrix(m, ncol = 1L)
        m
      })
      cn <- colnames(mat_list[[1L]])
      for (j in seq_len(ncol(mat_list[[1L]]))) {
        lab <- if (blk == "sigma_d2") {
          "sigma_d^2"
        } else if (blk == "sigma2") {
          "sigma^2 (raw)"
        } else if (!is.null(cn)) {
          paste0(blk, "[", cn[j], "]")
        } else {
          blk
        }
        panels[[length(panels) + 1L]] <- list(
          label = lab,
          series = lapply(mat_list, function(m) m[, j])
        )
      }
    }
  }

  n_panels <- length(panels)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  ncol_p <- ceiling(sqrt(n_panels))
  nrow_p <- ceiling(n_panels / ncol_p)
  graphics::par(mfrow = c(nrow_p, ncol_p))
  for (pnl in panels) {
    yr <- range(unlist(pnl$series), na.rm = TRUE)
    graphics::plot(seq_len(R_keep), pnl$series[[1L]], type = "l",
                   col = 1L, xlab = "iteration", ylab = pnl$label,
                   main = pnl$label, ylim = yr)
    if (n_chains > 1L) {
      for (cc in seq_len(n_chains)[-1L]) {
        graphics::lines(seq_len(R_keep), pnl$series[[cc]], col = cc)
      }
      graphics::legend("topright", legend = paste("chain", seq_len(n_chains)),
                        col = seq_len(n_chains), lty = 1, cex = 0.7, bty = "n")
    }
  }
  invisible(object)
}

#' Posterior-predictive share check for hierarchical Bayes fits
#'
#' Compares each alternative's observed take rate (share of choice
#' situations in which it was chosen, including the outside option) with its
#' posterior-predictive share from [predict.choicer_hb()]. Large systematic
#' gaps indicate model misfit — e.g. a missing covariate or an
#' outside-option share the delta level cannot rationalize.
#'
#' @param object A `choicer_hmnl` or `choicer_hmnp` fit (with
#'   `keep_data = TRUE`).
#' @param n_draws Posterior draws to integrate over (default 200).
#' @returns A `data.table` with columns `alternative`, `observed`,
#'   `predicted`, `lower`, `upper` (95% posterior-predictive interval), and
#'   `covered` (is the observed share inside the interval).
#' @examples
#' \donttest{
#' sim <- simulate_hmnl_data(N = 100, T = 3, J = 4, seed = 42)
#' fit <- suppressWarnings(run_hmnlogit(sim$data, "task", "alt", "choice", c("x1", "x2"),
#'                     person_col = "pid",
#'                     mcmc = list(R = 500, burn = 200)))
#' ppc_shares(fit)
#' }
#' @export
ppc_shares <- function(object, n_draws = 200L) {
  if (!inherits(object, "choicer_hb")) {
    stop("`object` must be a choicer_hmnl or choicer_hmnp fit.")
  }
  pred <- stats::predict(object, n_draws = n_draws)

  am <- object$alt_mapping
  spec <- object$data_spec
  obs_lab <- ifelse(am$alt_int == 0, "(outside)",
                    as.character(am[[spec$alt_col]]))
  obs <- stats::setNames(am$N_CHOICES / sum(am$N_CHOICES), obs_lab)

  out <- data.table::data.table(
    alternative = pred$alternative,
    observed = as.numeric(obs[pred$alternative]),
    predicted = pred$share,
    lower = pred$lower,
    upper = pred$upper
  )
  out[, covered := observed >= lower & observed <= upper]
  out[]
}
