# The suite's own probability engines: choicer's MNL, NL and MXL in choicer's
# exact parameterizations, written independently of choicer's kernels, so that
# a test of an estimate against the truth tests choicer and not a shared
# kernel. Definitions only (no side effects when sourced); base R + data.table.
#
# Usage: source("_validation/lib/simulate.R"); the self-check is
# `_validation/checks/wp0_simulate.R`.
#
# ---------------------------------------------------------------------------
# choicer's conventions (learned from fits of 0.2.1.9000, 2026-11-02; the
# self-check re-verifies every one of them against names(coef(fit)) and
# predict()):
#
# Data layout. Long rows, one per offered alternative, a situation's rows
# contiguous. choicer orders rows by (person_col, id_col, alternative code),
# so a DGP that writes rows in that order gets row-aligned predictions. The
# engines read a design list `d` (see design()): `sit` (integer 1..N per row,
# contiguous), `alt` (choicer's alternative code per row), `X`, `W` (numeric
# matrices, one row per row), and for an MXL panel `unit` (decision maker
# 1..U per situation, contiguous).
#
# Alternative codes. choicer codes the alternatives present in the data as
# their positions in sort(unique(labels)): numeric labels sort numerically
# (1, 2, 10), character labels in the session's collation. alt_mapping lists
# them in that order (alt_int 1..J); with include_outside_option = TRUE its
# first row is the outside option (alt_int 0, label NA).
#
# Outside option. No physical rows: the DGP writes inside alternatives only
# and the fit sets include_outside_option = TRUE. Its utility is 0 in every
# model; a situation with no chosen row chose it.
#
# Constants. Named ASC_<label>, in code order, after every other block.
# Without an outside option the first code's (the smallest label) is the
# reference, fixed at 0 and absent from theta; with one, every inside
# alternative has a constant. use_asc = FALSE: no constants, no `asc` block.
#
# Parameter order and names (fit$param_map gives the same blocks):
#   MNL  beta (X columns, covariate_cols order), asc
#   NL   beta, lambda (Lambda_1..Lambda_K), asc
#   MXL  beta, mu (Mu_<W column>, only with rc_mean = TRUE),
#        sigma (L_<i><j>), asc
#
# Nested logit (RU2, mlogit's un.nest.el = TRUE). Nests are coded by sorted
# nest label. Lambda_<k> numbers only the nests with two or more alternatives
# in the data, in nest-code order (k counts those nests, it is not the nest
# code); a singleton nest has lambda = 1 and no parameter. Within nest k the
# utility is divided by lambda_k:
#   log I_k = log sum_{j in k} exp(V_j / lambda_k),
#   P_j = exp(V_j / lambda_k - log I_k) * exp(lambda_k log I_k) / D,
#   D = sum_l exp(lambda_l log I_l) [+ 1 for the outside option, a degenerate
#   nest with V = 0 and lambda = 1], over the nests the situation offers.
#
# Mixed logit. Coefficient of random column k at draw eta (K_w-vector),
# g = L eta:
#   normal (rc_dist 0):     mu_k + g_k            (mu_k = 0 without rc_mean)
#   log-normal (rc_dist 1): exp(mu_k) + exp(g_k)  with rc_mean (shifted),
#                           exp(g_k)              without rc_mean.
# Utility: X beta + W b + delta[alt]. L is lower triangular, Sigma = L L' is
# the covariance of g. Correlated: `sigma` holds L's lower triangle row by
# row (L_11, L_21, L_22, L_31, ...), the diagonal as log L_ii, off-diagonal
# entries as they are. Uncorrelated: `sigma` holds log L_kk (the log standard
# deviations of g), named L_11, L_22, ... Draws: choicer's `eta_draws` cube
# K_w x S x U, one block per likelihood unit: situations in id order in the
# cross-section, decision makers in person order with person_col (their
# situations contiguous, person-major). get_halton_normals(S, U, K_w) is the
# store-mode cube. predict() never reads a supplied eta_draws cube: it forms
# store-mode points for one block per situation (a panel fit's too), equal to
# get_halton_normals(S, N, K_w) up to an ulp, so mxl_probs() on that cube
# (without d$unit) reproduces it. With every parameter fixed by named lower =
# upper, run_mxlogit(input_data =, eta_draws =)'s fit$loglik is choicer's
# simulated log-likelihood on the cube; there, prepare_mxl_data() needs the
# same rc_correlation as run_mxlogit() (the names follow input_data's).
#
# Weights. One per likelihood unit (situation; decision maker in a panel,
# where choicer reads the weight of the person's first situation). The
# log-likelihood is sum_u w_u log L_u with the weights as given.
#
# Seeds. Functions with a `seed` set it for their own draws and restore the
# caller's RNG state afterwards; seed = NULL draws from the current stream.
# ---------------------------------------------------------------------------

# --- Specifications and choicer's parameter vector ---------------------------

# A model specification, with choicer's defaults filled in. Fields: beta (X
# column names), alts (inside alternative labels in the data), outside
# (include_outside_option), use_asc; NL: nests (nest label of each of
# `alts`); MXL: random (W column names), rc_dist (0 normal, 1 log-normal per
# random column), rc_mean, correlated (rc_correlation).
.spec <- function(model, spec) {
  model <- match.arg(model, c("mnl", "nl", "mxl"))
  s <- utils::modifyList(list(beta = character(0), alts = NULL, outside = FALSE,
                              use_asc = TRUE, nests = NULL,
                              random = character(0), rc_dist = NULL,
                              rc_mean = FALSE, correlated = FALSE), spec)
  if (anyDuplicated(s$alts)) stop("`spec$alts` repeats a label.", call. = FALSE)
  if (s$use_asc && !length(s$alts)) {
    stop("`spec$alts` is needed with use_asc = TRUE.", call. = FALSE)
  }
  if (model == "nl" && length(s$nests) != length(s$alts)) {
    stop("`spec$nests` needs one nest label per element of `spec$alts`.",
         call. = FALSE)
  }
  if (model == "mxl") {
    if (!length(s$random)) stop("`spec$random` is empty.", call. = FALSE)
    if (is.null(s$rc_dist)) s$rc_dist <- rep(0L, length(s$random))
    if (length(s$rc_dist) != length(s$random) || !all(s$rc_dist %in% 0:1)) {
      stop("`spec$rc_dist` needs a 0 or 1 per random column.", call. = FALSE)
    }
    s$rc_dist <- as.integer(s$rc_dist)
  }
  s$model <- model
  s
}

# choicer's alternative codes: the labels' positions in this vector.
.alt_levels <- function(spec) sort(unique(spec$alts))

# The labels of the free constants, in theta's order.
.asc_labels <- function(spec) {
  if (!spec$use_asc) return(character(0))
  lev <- .alt_levels(spec)
  if (spec$outside) lev else lev[-1L]
}

# NL: nest levels, each alternative code's nest code, and the nests with a
# dissimilarity parameter (two or more alternatives).
.nest_info <- function(spec) {
  lev <- .alt_levels(spec)
  nest_lev <- sort(unique(spec$nests))
  nest_of_alt <- match(spec$nests[match(lev, spec$alts)], nest_lev)
  multi <- which(tabulate(nest_of_alt, length(nest_lev)) > 1L)
  list(levels = nest_lev, nest_of_alt = nest_of_alt, multi = multi)
}

# Lower-triangle entries of a K x K matrix, row by row.
.n_sigma <- function(K_w, correlated) {
  K_w <- as.integer(K_w)
  if (correlated) (K_w * (K_w + 1L)) %/% 2L else K_w
}

# Positions of each block in choicer's theta (fit$param_map's blocks; empty
# mu and asc blocks are absent, as there).
theta_blocks <- function(model, spec) {
  spec <- .spec(model, spec)
  out <- list(beta = seq_along(spec$beta))
  pos <- length(spec$beta)
  add <- function(name, n) {
    if (n > 0L) out[[name]] <<- pos + seq_len(n)
    pos <<- pos + n
  }
  if (spec$model == "nl") add("lambda", length(.nest_info(spec)$multi))
  if (spec$model == "mxl") {
    K_w <- length(spec$random)
    if (spec$rc_mean) add("mu", K_w)
    add("sigma", .n_sigma(K_w, spec$correlated))
  }
  add("asc", length(.asc_labels(spec)))
  out
}

# choicer's coefficient names, in its order: names(coef(fit)).
theta_names <- function(model, spec) {
  spec <- .spec(model, spec)
  K_w <- length(spec$random)
  lambda <- if (spec$model == "nl") {
    paste0("Lambda_", seq_along(.nest_info(spec)$multi), recycle0 = TRUE)
  }
  mu <- if (spec$model == "mxl" && spec$rc_mean) paste0("Mu_", spec$random)
  sigma <- if (spec$model == "mxl") {
    if (spec$correlated) {
      i <- rep(seq_len(K_w), seq_len(K_w))
      sprintf("L_%d%d", i, sequence(seq_len(K_w)))
    } else {
      sprintf("L_%d%d", seq_len(K_w), seq_len(K_w))
    }
  }
  c(spec$beta, lambda, mu, sigma,
    paste0("ASC_", .asc_labels(spec), recycle0 = TRUE))
}

# L from choicer's sigma block (the log-diagonal parameterization above).
L_from_theta <- function(sigma, K_w, correlated) {
  if (length(sigma) != .n_sigma(K_w, correlated)) {
    stop("`sigma` has the wrong length for K_w = ", K_w, ".", call. = FALSE)
  }
  if (!correlated) return(diag(exp(sigma), nrow = K_w))
  Lt <- matrix(0, K_w, K_w)
  Lt[upper.tri(Lt, diag = TRUE)] <- sigma   # column-major upper = row-major lower
  L <- t(Lt)
  diag(L) <- exp(diag(L))
  L
}

# choicer's sigma block from L (lower triangular, positive diagonal).
theta_from_L <- function(L, correlated) {
  if (any(L[upper.tri(L)] != 0) || any(diag(L) <= 0)) {
    stop("`L` must be lower triangular with a positive diagonal.", call. = FALSE)
  }
  if (!correlated) {
    if (any(L[lower.tri(L)] != 0)) {
      stop("An uncorrelated model needs a diagonal `L`.", call. = FALSE)
    }
    return(log(diag(L)))
  }
  Lt <- t(L)
  diag(Lt) <- log(diag(Lt))
  Lt[upper.tri(Lt, diag = TRUE)]
}

# A block's values in the order of `labels` (by name when named).
.block_values <- function(v, labels, what) {
  if (length(v) != length(labels) && is.null(names(v))) {
    stop("`blocks$", what, "` needs ", length(labels), " values.", call. = FALSE)
  }
  if (!is.null(names(v))) {
    miss <- setdiff(labels, names(v))
    if (length(miss)) {
      stop("`blocks$", what, "` lacks: ", paste(miss, collapse = ", "), ".",
           call. = FALSE)
    }
    v <- v[labels]
  }
  unname(as.numeric(v))
}

# choicer's theta from named blocks: beta; NL lambda (one per nest with a
# parameter, in nest order, or named by nest label, singletons at 1); MXL mu
# (with rc_mean), and sigma (choicer's block) or L (the matrix); asc (the
# free constants in theta's order, or named by alternative label; a full
# named set is normalized to the reference, delta_j - delta_ref, without an
# outside option). Returns the named vector in choicer's order.
pack_theta <- function(model, spec, blocks) {
  spec <- .spec(model, spec)
  unknown <- setdiff(names(blocks), c("beta", "mu", "sigma", "L", "lambda", "asc"))
  if (length(unknown)) {
    stop("Unknown blocks: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  pos <- theta_blocks(spec$model, spec)
  theta <- rep(NA_real_, sum(lengths(pos)))
  theta[pos$beta] <- .block_values(blocks[["beta"]], spec$beta, "beta")
  if (spec$model == "nl" && length(pos$lambda)) {
    ni <- .nest_info(spec)
    lam <- blocks[["lambda"]]
    if (!is.null(names(lam))) {
      single <- intersect(names(lam), ni$levels[-ni$multi])
      if (any(lam[single] != 1)) {
        stop("Singleton nests have lambda = 1: ", paste(single, collapse = ", "),
             ".", call. = FALSE)
      }
    }
    theta[pos$lambda] <- .block_values(lam, ni$levels[ni$multi], "lambda")
  }
  if (spec$model == "mxl") {
    K_w <- length(spec$random)
    if (spec$rc_mean) theta[pos$mu] <- .block_values(blocks[["mu"]], spec$random, "mu")
    sigma <- if (!is.null(blocks[["L"]])) theta_from_L(blocks[["L"]], spec$correlated) else blocks[["sigma"]]
    theta[pos$sigma] <- .block_values(unname(sigma),   # positional: choicer's layout
                                      seq_len(.n_sigma(K_w, spec$correlated)), "sigma")
  }
  if (length(pos$asc)) {
    asc <- blocks[["asc"]]
    ref <- .alt_levels(spec)[1L]
    if (!spec$outside && !is.null(names(asc)) && as.character(ref) %in% names(asc)) {
      asc <- asc - asc[[as.character(ref)]]
    }
    theta[pos$asc] <- .block_values(asc, as.character(.asc_labels(spec)), "asc")
  }
  if (anyNA(theta)) stop("Some blocks are missing or NA.", call. = FALSE)
  stats::setNames(theta, theta_names(spec$model, spec))
}

# The engines' parameter list from choicer's theta (named in any order, or
# unnamed in choicer's): beta, delta (per alternative code, the reference at
# 0; NULL without constants), the outside flag; NL lambda (per nest code,
# singletons at 1) and nest_of_alt; MXL mu (zeros without rc_mean), sigma, L,
# rc_dist, rc_mean, correlated.
unpack_theta <- function(model, spec, theta) {
  spec <- .spec(model, spec)
  nm <- theta_names(spec$model, spec)
  if (length(theta) != length(nm)) {
    stop("`theta` has ", length(theta), " values; the model has ", length(nm),
         ".", call. = FALSE)
  }
  if (!is.null(names(theta))) {
    if (!setequal(names(theta), nm)) {
      stop("`theta`'s names differ from choicer's: ",
           paste(setdiff(union(names(theta), nm), intersect(names(theta), nm)),
                 collapse = ", "), ".", call. = FALSE)
    }
    theta <- theta[nm]
  }
  theta <- stats::setNames(as.numeric(theta), nm)
  pos <- theta_blocks(spec$model, spec)
  lev <- .alt_levels(spec)
  par <- list(model = spec$model, beta = theta[pos$beta], outside = spec$outside,
              levels = lev, delta = NULL)
  if (spec$use_asc) {
    delta <- stats::setNames(numeric(length(lev)), lev)
    delta[as.character(.asc_labels(spec))] <- theta[pos$asc]
    par$delta <- delta
  }
  if (spec$model == "nl") {
    ni <- .nest_info(spec)
    par$lambda <- stats::setNames(rep(1, length(ni$levels)), ni$levels)
    par$lambda[ni$multi] <- theta[pos$lambda]
    par$nest_of_alt <- ni$nest_of_alt
  }
  if (spec$model == "mxl") {
    K_w <- length(spec$random)
    par$mu <- stats::setNames(if (spec$rc_mean) unname(theta[pos$mu]) else numeric(K_w),
                              spec$random)
    par$sigma <- theta[pos$sigma]
    par$L <- L_from_theta(par$sigma, K_w, spec$correlated)
    dimnames(par$L) <- list(spec$random, spec$random)
    par$rc_dist <- spec$rc_dist
    par$rc_mean <- spec$rc_mean
    par$correlated <- spec$correlated
  }
  par
}

# --- Data layout ---------------------------------------------------------------

# The columns `cols` of `dt` as a numeric matrix (zero columns allowed).
.cols_matrix <- function(dt, cols) {
  m <- matrix(0, nrow(dt), length(cols), dimnames = list(NULL, cols))
  for (k in seq_along(cols)) m[, k] <- dt[[cols[k]]]
  m
}

# The engines' design list from a long table whose situations' rows are
# contiguous and sorted by alternative code (choicer's row order): sit, alt
# (choicer's codes), X, W, and with `person` the decision maker of each
# situation (`unit`, contiguous).
design <- function(dt, spec, x = spec$beta, w = spec$random, id = "id",
                   alt = "alt", person = NULL) {
  sit <- data.table::rleid(dt[[id]])
  N <- sit[length(sit)]
  if (N != data.table::uniqueN(dt[[id]])) {
    stop("The rows of each situation must be contiguous.", call. = FALSE)
  }
  a <- match(dt[[alt]], .alt_levels(spec))
  if (anyNA(a)) stop("Alternatives outside `spec$alts`.", call. = FALSE)
  if (any(diff(a) <= 0L & diff(sit) == 0L)) {
    stop("Within a situation, rows must be sorted by alternative code, without ",
         "repeats.", call. = FALSE)
  }
  d <- list(sit = sit, alt = a, X = .cols_matrix(dt, x))
  if (length(w)) d$W <- .cols_matrix(dt, w)
  if (!is.null(person)) {
    first <- which(c(TRUE, diff(sit) != 0L))
    d$unit <- data.table::rleid(dt[[person]][first])
    if (d$unit[N] != data.table::uniqueN(dt[[person]])) {
      stop("Each decision maker's situations must be contiguous.", call. = FALSE)
    }
  }
  d
}

# The chosen row of each situation (0 for the outside option) from a 0/1
# choice column.
chosen_rows <- function(choice, sit) {
  chosen <- integer(sit[length(sit)])
  r <- which(choice == 1)
  chosen[sit[r]] <- r
  chosen
}

# Situation count; stops unless sit is 1..N with contiguous runs.
.n_sit <- function(sit) {
  n <- length(sit)
  if (!n || sit[1L] != 1L || is.unsorted(sit) || any(diff(sit) > 1L)) {
    stop("`sit` must run 1..N, each situation's rows contiguous.", call. = FALSE)
  }
  sit[n]
}

# --- Group operations (vectorized over situations) -------------------------------

# Per-group max and sum of v over integer groups g in 1..G (-Inf / 0 for a
# group without rows): data.table's GForce, no R loop over groups.
.group_max <- function(v, g, G) {
  r <- data.table::data.table(g = g, v = v)[, .(m = max(v)), by = g]
  out <- rep(-Inf, G)
  out[r$g] <- r$m
  out
}
.group_sum <- function(v, g, G) {
  r <- data.table::data.table(g = g, v = v)[, .(s = sum(v)), by = g]
  out <- numeric(G)
  out[r$g] <- r$s
  out
}

# Per-group log-sum-exp of v; `zero` adds a term exp(0) (the outside option).
.group_lse <- function(v, g, G, zero = FALSE) {
  m <- .group_max(v, g, G)
  if (zero) m <- pmax(m, 0)
  m[!is.finite(m)] <- 0
  s <- .group_sum(exp(v - m[g]), g, G)
  if (zero) s <- s + exp(-m)
  m + log(s)
}

# --- Utilities and probabilities ----------------------------------------------

# Systematic utility without random coefficients: X beta + delta[alt].
.utility <- function(d, par) {
  V <- if (length(par$beta)) drop(d$X %*% unname(par$beta)) else numeric(length(d$sit))
  if (!is.null(par$delta)) V <- V + unname(par$delta)[d$alt]
  V
}

# Logit log-probabilities of V within situations: logp per row, logp0 (the
# outside option's, -Inf without one) per situation.
.logit_logp <- function(V, sit, N, outside) {
  D <- .group_lse(V, sit, N, zero = outside)
  list(logp = V - D[sit], logp0 = if (outside) -D else rep(-Inf, N))
}

# Probabilities from log-probabilities: p per row, p0 per situation.
.probs <- function(lp) list(p = exp(lp$logp), p0 = exp(lp$logp0))

# MNL utilities X beta + delta[alt].
mnl_utility <- function(d, par) .utility(d, par)

# MNL probabilities: p (inside, per row) and p0 (outside, per situation; 0
# without an outside option).
mnl_probs <- function(d, par) {
  .probs(.logit_logp(.utility(d, par), d$sit, .n_sit(d$sit), par$outside))
}

# NL log-probabilities (RU2; see the header): logp per row, logp0 per
# situation.
.nl_logp <- function(d, par, N) {
  V <- .utility(d, par)
  lam_nest <- unname(par$lambda)
  K <- length(lam_nest)
  nest <- par$nest_of_alt[d$alt]
  key <- (as.numeric(d$sit) - 1) * K + (nest - 1)   # one group per (situation, nest)
  ukey <- unique(key)
  g <- match(key, ukey)
  G <- length(ukey)
  vl <- V / lam_nest[nest]
  log_I <- .group_lse(vl, g, G)
  term <- lam_nest[ukey %% K + 1] * log_I             # lambda_k log I_k
  D <- .group_lse(term, ukey %/% K + 1, N, zero = par$outside)
  list(logp = vl - log_I[g] + term[g] - D[d$sit],
       logp0 = if (par$outside) -D else rep(-Inf, N))
}

# NL probabilities: p (per row), p0 (per situation).
nl_probs <- function(d, par) .probs(.nl_logp(d, par, .n_sit(d$sit)))

# MXL realized coefficients at draws eta (K_w x m): a K_w x m matrix, choicer's
# convention (header).
rc_coefs <- function(mu, L, eta, rc_dist, rc_mean) {
  eta <- matrix(eta, nrow(L))
  b <- L %*% eta
  ln <- rc_dist == 1L
  if (any(ln)) b[ln, ] <- exp(b[ln, , drop = FALSE])
  if (rc_mean) b <- b + ifelse(ln, exp(mu), mu)
  b
}

# Each decision maker's (or situation's) realized coefficients from R's RNG:
# a U x K_w matrix, eta ~ N(0, I) drawn unit by unit.
draw_tastes <- function(U, mu, L, rc_dist, rc_mean, seed = NULL) {
  K_w <- nrow(L)
  eta <- .with_seed(seed, matrix(stats::rnorm(K_w * U), K_w, U))
  t(rc_coefs(mu, L, eta, rc_dist, rc_mean))
}

# Each row's likelihood unit: its decision maker (panel) or its situation.
.unit_of_row <- function(d) if (is.null(d$unit)) d$sit else d$unit[d$sit]

# W b for coefficients B (K_w x U, one column per unit).
.random_part <- function(W, B, unit_row) {
  v <- 0
  for (k in seq_len(ncol(W))) v <- v + W[, k] * B[k, unit_row]
  v
}

# MXL probabilities given each unit's tastes (U x K_w, as draw_tastes()
# returns): the DGP's conditional logit. p per row, p0 per situation.
mxl_cond_probs <- function(d, par, tastes) {
  N <- .n_sit(d$sit)
  V <- .utility(d, par) + .random_part(d$W, t(tastes), .unit_of_row(d))
  .probs(.logit_logp(V, d$sit, N, par$outside))
}

# Checks a draws cube against the design: K_w x S x U, U units.
.check_cube <- function(d, eta) {
  U <- if (is.null(d$unit)) d$sit[length(d$sit)] else d$unit[length(d$unit)]
  if (length(dim(eta)) != 3L || dim(eta)[1L] != ncol(d$W) || dim(eta)[3L] != U) {
    stop("`eta` must be a K_w x S x U cube (U = ", U, " units).", call. = FALSE)
  }
  U
}

# Per-draw log-probabilities, handed to f(lp, s) draw by draw.
.mxl_draws <- function(d, par, eta, f) {
  N <- .n_sit(d$sit)
  U <- .check_cube(d, eta)
  K_w <- dim(eta)[1L]
  base <- .utility(d, par)
  unit_row <- .unit_of_row(d)
  for (s in seq_len(dim(eta)[2L])) {
    B <- rc_coefs(par$mu, par$L, matrix(eta[, s, ], K_w, U), par$rc_dist, par$rc_mean)
    f(.logit_logp(base + .random_part(d$W, B, unit_row), d$sit, N, par$outside), s)
  }
}

# MXL simulated probabilities on a draws cube (one block per unit: situation,
# or decision maker with d$unit): the mean over draws of the logit
# probabilities. p per row, p0 per situation.
mxl_probs <- function(d, par, eta) {
  p <- 0
  p0 <- 0
  .mxl_draws(d, par, eta, function(lp, s) {
    p <<- p + exp(lp$logp)
    p0 <<- p0 + exp(lp$logp0)
  })
  list(p = p / dim(eta)[2L], p0 = p0 / dim(eta)[2L])
}

# --- Choices ---------------------------------------------------------------------

# One choice per situation by inverse CDF, one uniform per situation (the
# outside option first, then the rows in order): `u` when given (drawn by
# the caller), else runif(N) under `seed`. Returns `chosen` (the chosen row
# per situation, 0 for the outside option) and `choice` (0/1 per row).
draw_choice <- function(P_inside, P_outside = NULL, sit, seed = NULL, u = NULL) {
  N <- .n_sit(sit)
  if (!is.null(u) && !is.null(seed)) stop("Give `seed` or `u`, not both.", call. = FALSE)
  if (!is.null(u) && length(u) != N) stop("`u` needs one uniform per situation.", call. = FALSE)
  M <- tabulate(sit, N)
  start <- cumsum(c(1L, M[-N]))
  p0 <- if (is.null(P_outside)) numeric(N) else P_outside
  # Each situation's cumulative inside probabilities, by position in the set
  cum <- numeric(length(sit))
  last <- integer(N)                                 # last position with p > 0
  ord <- order(M, decreasing = TRUE)
  n_ge <- rev(cumsum(rev(tabulate(M))))              # situations with >= k rows
  for (k in seq_len(max(M))) {
    s <- ord[seq_len(n_ge[k])]
    r <- start[s] + (k - 1L)
    cum[r] <- if (k == 1L) P_inside[r] else cum[r - 1L] + P_inside[r]
    pos <- P_inside[r] > 0
    last[s[pos]] <- k
  }
  total <- p0 + cum[start + M - 1L]
  if (any(!(total > 0))) stop("A situation has no positive probability.", call. = FALSE)
  thr <- (if (is.null(u)) .with_seed(seed, stats::runif(N)) else u) * total - p0
  k <- pmin(tabulate(sit[cum <= thr[sit]], N) + 1L, last)
  chosen <- ifelse(thr < 0, 0L, start + k - 1L)
  choice <- integer(length(sit))
  choice[chosen[chosen > 0L]] <- 1L
  list(chosen = chosen, choice = choice)
}

# --- Log-likelihoods ---------------------------------------------------------------

# Each situation's log-probability of its chosen row (0: the outside option).
.chosen_logp <- function(lp, chosen) {
  out <- lp$logp0
  i <- chosen > 0L
  out[i] <- lp$logp[chosen[i]]
  out
}

# sum_u w_u ll_u, or the units' contributions w_u ll_u.
.weighted <- function(ll, weights, by_unit) {
  if (!is.null(weights)) {
    if (length(weights) != length(ll)) {
      stop("`weights` needs one value per likelihood unit (", length(ll), ").",
           call. = FALSE)
    }
    ll <- weights * ll
  }
  if (by_unit) ll else sum(ll)
}

# MNL log-likelihood at par for chosen rows (chosen_rows()), weighted.
mnl_loglik <- function(d, par, chosen, weights = NULL, by_unit = FALSE) {
  lp <- .logit_logp(.utility(d, par), d$sit, .n_sit(d$sit), par$outside)
  .weighted(.chosen_logp(lp, chosen), weights, by_unit)
}

# NL log-likelihood at par, weighted.
nl_loglik <- function(d, par, chosen, weights = NULL, by_unit = FALSE) {
  .weighted(.chosen_logp(.nl_logp(d, par, .n_sit(d$sit)), chosen), weights, by_unit)
}

# MXL simulated log-likelihood on a draws cube: sum_u w_u log mean_s
# prod_{t in u} P_ts (the product over a decision maker's situations with
# d$unit, one situation per unit otherwise), log-mean-exp over draws.
mxl_loglik <- function(d, par, eta, chosen, weights = NULL, by_unit = FALSE) {
  U <- .check_cube(d, eta)
  m <- rep(-Inf, U)                                  # running max over draws
  acc <- numeric(U)                                  # sum exp(a_s - m)
  .mxl_draws(d, par, eta, function(lp, s) {
    a <- .chosen_logp(lp, chosen)
    if (!is.null(d$unit)) a <- .group_sum(a, d$unit, U)
    up <- a > m
    e <- exp(a - m)
    e[a == -Inf] <- 0                                # a zero-probability draw
    acc <<- ifelse(up, acc * exp(m - a) + 1, acc + e)
    m <<- pmax(m, a)
  })
  .weighted(m + log(acc / dim(eta)[2L]), weights, by_unit)
}

# --- RNG ---------------------------------------------------------------------------

# Evaluates expr under set.seed(seed), restoring the caller's RNG state; the
# current stream when seed is NULL.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  env <- globalenv()
  old <- if (exists(".Random.seed", envir = env, inherits = FALSE)) {
    get(".Random.seed", envir = env, inherits = FALSE)
  }
  on.exit(if (is.null(old)) {
    rm(".Random.seed", envir = env)
  } else {
    assign(".Random.seed", old, envir = env)
  })
  set.seed(seed)
  expr
}
