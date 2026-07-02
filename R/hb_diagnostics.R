# Convergence diagnostics for the hierarchical Bayes models (choicer_hmnl /
# choicer_hmnp). Phase 1 ships split-R-hat; the posterior-predictive share
# check arrives with the post-estimation phase of
# _plans/hierarchical_bayes_plan.md.

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
#' @returns Named numeric vector with one \eqn{\widehat{R}} per parameter
#'   (`NA` for parameters with zero variance).
#' @examples
#' set.seed(42)
#' draws <- matrix(rnorm(2000), ncol = 2,
#'                 dimnames = list(NULL, c("a", "b")))
#' rhat(draws)          # ~1: white noise is stationary
#' drifting <- cbind(a = cumsum(rnorm(1000)))
#' rhat(drifting)       # >> 1: a random walk is not
#' @export
rhat <- function(draws) {
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
  out
}

#' Build the R-hat table for a hierarchical Bayes fit
#'
#' Internal: assembles split-R-hat over the b, theta, and sigma_d^2 draws of
#' one or more chains into a named vector, used by [run_hmnlogit()] (and
#' `run_hmnprobit()`) for the convergence warning and the summary display.
#'
#' @param chain_list List with one element per chain, each a list holding
#'   matrices `b`, `theta` and vector `sigma_d2`.
#' @returns Named numeric vector of split-R-hat values.
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
    rhat(bind_block("b")),
    rhat(bind_block("theta")),
    rhat(bind_block("sigma_d2"))
  )
}
