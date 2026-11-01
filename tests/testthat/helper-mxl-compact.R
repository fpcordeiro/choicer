# Many-alternative designs with sparse choice sets for the mixed logit's
# derivative kernels (test-mxl-compact-derivatives.R) and the share
# Jacobian (test-mxl-blp-jacobian.R): the normalized reference
# alternative, the outside option, panels whose situations offer different
# alternatives, a situation listing an alternative twice, the reference
# never offered, extra trailing ASCs.

# A prepared many-alternative design in the stacked layout, with the fields
# mxlp_call() reads. U decision makers with 1-3 situations each (panel) or one
# (cross-section); each situation offers `m_range` distinct alternatives drawn
# from 1..J_codes (from 2..J_codes when `absent_ref`), plus one repeated row
# when `dup`; weights vary across decision makers.
cdx_fixture <- function(name, seed, panel = TRUE,
                        include_outside_option = FALSE, J = 30L,
                        m_range = c(5L, 10L), U = if (panel) 10L else 24L,
                        rc_dist = c(0L, 1L), rc_correlation = TRUE,
                        rc_mean = TRUE, use_asc = TRUE, K_x = 2L, S = 12L,
                        dup = FALSE, absent_ref = FALSE, J_codes = J,
                        ref_only = FALSE, shuffle = FALSE) {
  set.seed(seed)
  rc_dist <- as.integer(rc_dist)
  K_w <- length(rc_dist)
  Ti <- if (panel) sample(1:3, U, replace = TRUE) else rep(1L, U)
  N <- sum(Ti)
  unit <- rep(seq_len(U), Ti)
  pool <- if (absent_ref) 2:J_codes else seq_len(J_codes)
  alts <- lapply(seq_len(N), function(t) {
    # `ref_only`: every situation of unit 1, and about a third of the others,
    # offer only the reference alternative (no free ASC)
    if (ref_only && (unit[t] == 1L || stats::runif(1L) < 1 / 3)) return(1L)
    m <- sample(m_range[1L]:m_range[2L], 1L)
    a <- sort(pool[sample.int(length(pool), min(m, length(pool)))])
    if (dup) a <- sort(c(a, a[sample.int(length(a), 1L)]))
    # `shuffle`: rows in no particular order (repeats not adjacent)
    if (shuffle) a <- a[sample.int(length(a))]
    a
  })
  M <- lengths(alts)
  n_rows <- sum(M)
  lo <- if (include_outside_option) 0L else 1L
  choice_idx <- vapply(M, function(m) sample(lo:m, 1L), integer(1L))
  w_u <- stats::runif(U, 0.5, 2)
  L_size <- if (rc_correlation) (K_w * (K_w + 1L)) %/% 2L else K_w
  n_delta <- if (!use_asc) 0L else if (include_outside_option) J else J - 1L
  theta <- c(stats::rnorm(K_x, sd = 0.4), if (rc_mean) stats::rnorm(K_w, sd = 0.3),
             stats::rnorm(L_size, sd = 0.3) - 0.5, stats::rnorm(n_delta, sd = 0.4))
  list(
    name = name, X = matrix(stats::rnorm(n_rows * K_x), n_rows, K_x),
    W = matrix(stats::rnorm(n_rows * K_w), n_rows, K_w),
    alt_idx = as.integer(unlist(alts)), choice_idx = as.integer(choice_idx),
    M = as.integer(M), weights = rep(w_u, Ti),
    Ti = if (panel) as.integer(Ti) else NULL,
    eta = get_halton_normals(S, U, K_w), theta = theta, rc_dist = rc_dist,
    rc_correlation = rc_correlation, rc_mean = rc_mean, use_asc = use_asc,
    include_outside_option = include_outside_option, K_w = K_w, K_x = K_x,
    J = J, S = S, U = U, N = N
  )
}

# Index of the first ASC in theta, and the global parameter indices of the
# compact block of decision maker u: the continuous parameters and the free
# ASCs of the alternatives in u's rows.
cdx_delta_start <- function(fx) {
  L_size <- if (fx$rc_correlation) (fx$K_w * (fx$K_w + 1L)) %/% 2L else fx$K_w
  fx$K_x + (if (fx$rc_mean) fx$K_w else 0L) + L_size + 1L
}
cdx_block <- function(fx, u) {
  Ti <- if (is.null(fx$Ti)) rep(1L, fx$N) else fx$Ti
  sit <- which(rep(seq_along(Ti), Ti) == u)
  row_off <- c(0L, cumsum(fx$M))
  codes <- fx$alt_idx[unlist(lapply(sit, function(t) (row_off[t] + 1L):row_off[t + 1L]))]
  free <- if (fx$include_outside_option) codes else codes[codes > 1L] - 1L
  d0 <- cdx_delta_start(fx)
  c(seq_len(d0 - 1L), sort(unique(free)) + d0 - 1L)
}

cdx_cells <- function() {
  list(
    cdx_fixture("xsec", 1001, panel = FALSE),
    cdx_fixture("xsec, outside option", 1002, panel = FALSE,
                include_outside_option = TRUE),
    cdx_fixture("panel", 1003),
    cdx_fixture("panel, outside option", 1004, include_outside_option = TRUE),
    cdx_fixture("panel, repeated alternative", 1005, dup = TRUE),
    cdx_fixture("panel, reference never offered", 1006, absent_ref = TRUE),
    cdx_fixture("panel, trailing ASCs", 1007, J_codes = 24L),
    cdx_fixture("panel, normal, uncorrelated, no means", 1008,
                rc_dist = c(0L, 0L), rc_correlation = FALSE, rc_mean = FALSE),
    cdx_fixture("panel, no ASCs", 1009, use_asc = FALSE),
    cdx_fixture("one decision maker", 1010, U = 1L),
    cdx_fixture("panel, situations with the reference only", 1020,
                ref_only = TRUE),
    cdx_fixture("xsec, shuffled rows, repeated alternatives", 1021,
                panel = FALSE, dup = TRUE, shuffle = TRUE)
  )
}
