# ============================================================================
# Simulation-based calibration (plan §2, claim C7). If theta ~ prior,
# y | theta, and L posterior draws are taken given y, the rank of any f(theta)
# among the draws' f is uniform on 0..L (Cook, Gelman & Rubin 2006; Talts et
# al. 2018; Modrák et al. 2023). Definitions only. The test half (ranks, the
# ECDF test and its plot) comes first; the loop half (thinning, ESS and each
# replication's diag$sbc) follows it. The prior draws are the DGPs'
# (lib/dgp_bayes.R).
# ============================================================================

# --- SBC: ranks and the uniformity test -------------------------------------

# The ranks of `truth` among posterior draws: the count of draws below it,
# ties split uniformly at random (the caller's RNG, drawn only where there
# are ties), an integer in 0..L. `draws` is a vector of L draws for a scalar
# truth, or an L x P matrix whose column p holds truth[p]'s draws. A truth
# or a draw that is NA gives NA.
sbc_rank <- function(truth, draws) {
  draws <- as.matrix(draws)
  if (ncol(draws) != length(truth)) {
    stop("sbc_rank(): `draws` must be L x length(truth) (one column per quantity)")
  }
  L <- nrow(draws)
  tm <- matrix(truth, L, length(truth), byrow = TRUE)
  below <- colSums(draws < tm)
  ties <- colSums(draws == tm)
  u <- numeric(length(truth))
  k <- which(ties > 0)
  if (length(k)) u[k] <- stats::runif(length(k))
  as.integer(below + floor(u * (ties + 1)))
}

.sbc_cache <- new.env(parent = emptyenv())

# The two-sided binomial tail probability of every count c = 0..R at every
# evaluation point z_i = i/(L+1): 2 min(P(X <= c), P(X >= c)) capped at 1,
# X ~ Binomial(R, z_i); an (R+1) x L matrix.
.sbc_tails <- function(R, L) {
  key <- paste("tails", R, L)
  if (!is.null(.sbc_cache[[key]])) return(.sbc_cache[[key]])
  z <- rep(seq_len(L) / (L + 1), each = R + 1)
  cnt <- rep(0:R, L)
  tl <- pmin(1, 2 * pmin(stats::pbinom(cnt, R, z),
                         stats::pbinom(cnt - 1, R, z, lower.tail = FALSE)))
  # (point i, count c) and (point L + 1 - i, count R - c) are equal in exact
  # arithmetic (z_{L+1-i} = 1 - z_i), but pbinom's two tails differ in the
  # last bits; made bitwise equal, so that `tails > s` cannot split a tie and
  # a rank vector and its mirror image (rank -> L - rank) get the same p.
  tl <- matrix(tl, R + 1, L)
  .sbc_cache[[key]] <- pmin(tl, tl[(R + 1):1, L:1, drop = FALSE])
}

# The ECDF test's statistic for count vectors: the smallest tail probability
# over the evaluation points, where cum[i] = #{rank < i} = R F(z_i). `counts`
# is an (L+1) x n matrix of rank counts (one column per sample).
.sbc_stat <- function(counts, tails) {
  L <- ncol(tails)
  cum <- 0
  stat <- rep(Inf, ncol(counts))
  for (i in seq_len(L)) {
    cum <- cum + counts[i, ]
    stat <- pmin(stat, tails[cbind(cum + 1, i)])
  }
  stat
}

# Point i's acceptance set at level s, {c : tail probability of c > s}, as an
# interval [lo, hi] of counts (the tail probability is unimodal in c, so the
# set is an interval; checked), or NULL when empty.
.sbc_accept_set <- function(tails, i, s) {
  ok <- which(tails[, i] > s) - 1L
  if (!length(ok)) return(NULL)
  if (any(diff(ok) != 1L)) stop("internal: acceptance set at point ", i, " is not an interval")
  c(ok[1], ok[length(ok)])
}

# The exact null probability P(min_i tail_i > s) that the ECDF stays inside
# every point's acceptance set, for R ranks uniform on 0..L. The rank counts
# n_0..n_L are multinomial, i.e. L + 1 independent Poisson(R / (L + 1))
# counts conditioned on their sum being R; so the cumulative count c_i =
# #{rank < i} is a Poisson walk, c_i = c_{i-1} + n_{i-1}, whose mass is
# carried over the acceptance intervals one point at a time (a convolution
# with the Poisson pmf restricted to the interval), then multiplied by the
# last bin's pmf at R - c_L and divided by P(sum = R) = dpois(R, R).
# O(L w^2) for intervals of width w (at most R + 1).
.sbc_inside <- function(R, L, s) {
  tails <- .sbc_tails(R, L)
  pois <- stats::dpois(0:R, R / (L + 1))
  lo <- 0L
  hi <- 0L
  m <- 1
  for (i in seq_len(L)) {
    a <- .sbc_accept_set(tails, i, s)
    if (is.null(a)) return(0)
    nlo <- max(a[1], lo)                  # the walk never decreases
    nhi <- a[2]
    if (nlo > nhi) return(0)
    d <- c(hi - lo + 1L, nhi - nlo + 1L)
    k <- .col(d) - .row(d) + (nlo - lo)   # increment from count lo + r - 1 to nlo + s - 1
    B <- matrix(0, d[1], d[2])
    B[k >= 0] <- pois[k[k >= 0] + 1L]
    m <- drop(m %*% B)
    lo <- nlo
    hi <- nhi
  }
  sum(m * pois[R - (lo:hi) + 1L]) / stats::dpois(R, R)
}

# The level gamma and the simultaneous band of the test at `alpha`, cached by
# (R, L, alpha): gamma is the largest attainable tail probability s with
# P(min_i tail_i <= s) <= alpha (binary search over the table's distinct
# values above alpha / L, where the union bound already gives <= alpha), and
# the band for F(z_i) is point i's acceptance set at gamma over R. So the
# ECDF leaves the band exactly when stat <= gamma, i.e. when p <= alpha.
.sbc_band <- function(R, L, alpha) {
  key <- paste("band", R, L, alpha)
  if (!is.null(.sbc_cache[[key]])) return(.sbc_cache[[key]])
  tails <- .sbc_tails(R, L)
  u <- sort(unique(as.vector(tails)))
  base <- max(c(0, u[u <= alpha / L]))
  u <- u[u > base]
  lo <- 0L                                # u[lo] (or base) has P(min <= s) <= alpha
  hi <- length(u) + 1L                    # u[hi] (or beyond) has it > alpha
  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    if (1 - .sbc_inside(R, L, u[mid]) <= alpha) lo <- mid else hi <- mid
  }
  gamma <- if (lo) u[lo] else base
  band <- vapply(seq_len(L), function(i) .sbc_accept_set(tails, i, gamma) %||% c(NA, NA),
                 numeric(2))
  .sbc_cache[[key]] <- list(gamma = gamma, lower = band[1, ] / R, upper = band[2, ] / R)
}

# The ECDF uniformity test of SBC ranks (Säilynoja, Bürkner & Vehtari 2022).
# At z_i = i/(L+1), i = 1..L, R F(z_i) = #{rank < i} is Binomial(R, z_i) under
# uniformity; `stat` is the smallest two-sided binomial tail probability over
# i; `p` = P(min_i tail_i <= stat) under the exact discrete null (ranks
# uniform on 0..L; .sbc_inside()), with no Monte Carlo floor; `gamma` and
# `lower`/`upper` are the level and the simultaneous 1 - alpha band for
# F(z_i) (.sbc_band()), so that the ECDF leaves the band exactly when
# p <= alpha. NA ranks are dropped; R is the count used.
sbc_ecdf_test <- function(ranks, L, alpha = 0.05) {
  if (!is.numeric(L) || length(L) != 1L || L < 1 || L != floor(L)) {
    stop("sbc_ecdf_test(): `L` must be a positive whole number")
  }
  ranks <- ranks[!is.na(ranks)]
  if (any(ranks < 0 | ranks > L | ranks != floor(ranks))) {
    stop("sbc_ecdf_test(): ranks must be whole numbers in 0..L = 0..", L)
  }
  R <- length(ranks)
  z <- seq_len(L) / (L + 1)
  if (R < 1) {
    return(list(stat = NA_real_, p = NA_real_, gamma = NA_real_, z = z,
                lower = rep(NA_real_, L), upper = rep(NA_real_, L), R = 0L))
  }
  counts <- matrix(tabulate(ranks + 1, L + 1), ncol = 1)
  stat <- .sbc_stat(counts, .sbc_tails(R, L))
  band <- .sbc_band(R, L, alpha)
  list(stat = stat, p = max(0, 1 - .sbc_inside(R, L, stat)), gamma = band$gamma, z = z,
       lower = band$lower, upper = band$upper, R = R)
}

# The ECDF difference F(z) - z of SBC ranks with the test's simultaneous
# bands (base graphics, on the current device); returns the test invisibly.
sbc_ecdf_plot <- function(ranks, L, test = NULL, main = "") {
  if (is.null(test)) test <- sbc_ecdf_test(ranks, L)
  ranks <- ranks[!is.na(ranks)]
  x <- c(0, test$z, 1)
  if (!test$R) {
    graphics::plot(0:1, c(-0.1, 0.1), type = "n", xlab = "z", ylab = "ECDF difference", main = main)
    graphics::text(0.5, 0, "no ranks")
    return(invisible(test))
  }
  Fz <- cumsum(tabulate(ranks + 1, L + 1))[seq_len(L)] / test$R
  y <- c(0, Fz - test$z, 0)
  lo <- c(0, test$lower - test$z, 0)
  hi <- c(0, test$upper - test$z, 0)
  graphics::plot(range(x), range(y, lo, hi), type = "n", xlab = "z = i / (L + 1)",
                 ylab = "ECDF difference", main = main)
  graphics::polygon(c(x, rev(x)), c(lo, rev(hi)), col = grDevices::grey(0.85), border = NA)
  graphics::abline(h = 0, col = grDevices::grey(0.5))
  graphics::lines(x, y, type = "s")
  graphics::mtext(sprintf("R = %d, L = %d, p = %.3g", test$R, L, test$p), side = 3, line = 0.2,
                  cex = 0.8)
  invisible(test)
}

# --- SBC: thinning, ESS and each replication's record (WP8) -----------------

# A replication of SBC is one unit of the runner (run_unit()): the family's
# simulate() draws the truth from the prior and the data given it (F9:
# db_simulate(prior = ), lib/dgp_bayes.R), and its method runs the sampler
# with the same prior and passes the chain to sbc_diag(), whose result is
# the record's diag$sbc (metrics_cell()'s sbc_ecdf rows, report.R's plots).

# The indices of L draws equally spaced over a chain of n kept draws, the
# last included: floor(l n / L), l = 1..L.
sbc_thin_index <- function(n, L = 99L) {
  if (n < L) stop("sbc_thin_index(): the chain holds ", n, " draws, fewer than L = ", L)
  floor(seq_len(L) * n / L)
}

# The bulk effective sample size of each column of one chain's draws (n x
# q), choicer::ess()'s rank-normalized estimate (Vehtari et al. 2021); NA for
# a constant column.
sbc_ess <- function(draws) {
  draws <- as.matrix(draws)
  stats::setNames(choicer::ess(draws)[, "bulk"], colnames(draws))
}

# One replication's diag$sbc from the truth (a named numeric: the
# quantities) and one chain's kept draws (n x q, columns named as truth):
# L draws equally spaced (sbc_thin_index()) and the rank of every quantity
# among them (sbc_rank(); ties split by the caller's RNG); with `loglik` (the
# log-likelihood at the truth), also the rank of the log-likelihood
# (Modrák et al. 2023), `loglik_draws` holding it at every kept draw (length
# n, then thinned like the draws) or at the L thinned draws. Plan §2: the
# chain should hold an ESS of at least L for every quantity, so that the L
# draws are nearly independent; a replication below `min_ess` is flagged,
# not dropped. Returns list(L, ranks, ess (per quantity, the chain's; the
# log-likelihood's when given at every draw), ess_min, flagged, n).
sbc_diag <- function(truth, draws, L = 99L, loglik = NULL, loglik_draws = NULL, min_ess = L) {
  draws <- as.matrix(draws)
  q <- names(truth)
  if (is.null(q) || !all(q %in% colnames(draws))) {
    stop("sbc_diag(): `truth` must be named by columns of `draws`")
  }
  draws <- draws[, q, drop = FALSE]
  n <- nrow(draws)
  idx <- sbc_thin_index(n, L)
  ess <- sbc_ess(draws)
  ranks <- sbc_rank(unname(truth), draws[idx, , drop = FALSE])
  names(ranks) <- q
  if (!is.null(loglik)) {
    ll <- as.numeric(loglik_draws)
    if (length(ll) == n) {
      ess <- c(ess, loglik = sbc_ess(cbind(loglik = ll))[[1]])
      ll <- ll[idx]
    } else if (length(ll) != L) {
      stop("sbc_diag(): `loglik_draws` must hold the log-likelihood at the n = ", n,
           " kept draws or at the L = ", L, " thinned ones")
    }
    ranks <- c(ranks, loglik = sbc_rank(loglik, ll))
  }
  ess_min <- suppressWarnings(min(ess, na.rm = TRUE))
  list(L = as.integer(L), ranks = ranks, ess = ess, ess_min = ess_min,
       flagged = !is.finite(ess_min) || ess_min < min_ess, n = n)
}
