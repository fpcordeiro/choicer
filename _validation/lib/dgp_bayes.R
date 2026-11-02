# ============================================================================
# DB: the Bayesian samplers' data (F9; plan §1 C7, §5 F9, Q1.4, Q4). Prior
# draws on the samplers' own (raw) scales, data from a draw, the identified
# quantities SBC ranks, their log-likelihoods, the samplers' fits, and the
# quick tier's Q1.4 (choicer's MNP against bayesm) and Q4 (mini-SBC) checks.
# Definitions only, except the registration of DGP "DB" (dgp_register()).
#
# Variants (the model each one's data follow):
#   "default"  MNP (run_mnprobit()) on D2's design: travellers resampled from
#              mode_choice (get_dgp("D2")), J = 4 all offered, covariates wait,
#              travel, vcost and income_air divided by DB_X_SCALE (half-hours,
#              4 hours, $25, $20k: each differenced column has an SD near 1).
#              Levels U_j = x_j' beta + a_j + e_j, e ~ N(0, Omega); in
#              differences against car (base_alt = "car"), w_i = Xd_i beta +
#              eps_i, eps_i ~ N_3(0, Sigma) over air, bus, train (choicer's
#              order), with constants ASC_air, ASC_bus, ASC_train; y_i = car if
#              max_j w_ij < 0, else the argmax. Identified: beta / sqrt(Sigma_11),
#              Sigma / Sigma_11 (Sigma_11 the variance of e_air - e_car). Why car,
#              not choicer's default base air: the normalization then falls on
#              bus - air, bus the least chosen mode (4-6%), and the identified
#              draws mixed 2-10 times slower (bulk ESS per 10,000 draws at
#              N = 500, three datasets: minimum 6-53 against 45-57; 2026-11-02).
#   "hmnl"     HMNL (run_hmnlogit()) on D1-L (get_dgp("D1"), variant
#              "labeled", shock "logit"): b, W, theta, sigma_d2 and xi through
#              WP4's `theta` override; the persons' tastes beta_i ~ N(b, W).
#   "hmnp"     HMNP (run_hmnprobit()) on D1-L with shock "probit", WP4's
#              model at sigma = 1. The sampler's raw scale has a free shock
#              variance sigma^2; a raw draw (b, W, theta, sigma_d2, xi,
#              sigma^2) enters the override identified, as b / sigma,
#              W / sigma^2, theta / sigma, sigma_d2 / sigma^2, xi / sigma, which
#              gives the data the same distribution.
#
# theta0 (DB_THETA0; BvM's fixed truths, identified, choicer's draw names).
# "default", 2026-11-02 (db_mnp_theta0()): the probit counterpart of D2's MNL
# theta0 (its WESML fit of mode_choice, D2_THETA0$default): level
# coefficients beta_MNL * scale / (pi / sqrt(6)) (Gumbel errors have SD
# pi / sqrt(6)); Omega with air's error independent and the ground modes'
# correlated 0.5, all of variance 1, so that Sigma_11 = 2, Sigma / Sigma_11 =
# (1, 0.25, 0.5, 0.25, 0.25, 0.5); the constants then re-solved so that the
# population's expected shares (D2's 210 traveller types and weights, GHK
# with 4,000 draws) equal Q = (air .14, bus .09, car .64, train .13), as D2's
# theta0 does (db_mnp_calibrate(): 39 iterations, log-share gaps below 1e-10;
# GHK's own error in the shares, another seed's draws, 6e-5; without the
# step the shares were .152/.069/.680/.099). "hmnl": D1-L's theta0 (SPL_TRUE, WP4,
# 2026-11-02), copied. "hmnp": the same divided by the probit scale
# pi / sqrt(6) (W and sigma_d2 by its square).
#
# Priors (DB_PRIOR; the samplers' arguments, passed to them explicitly).
# Proper and moderately informative; sigma_d2 ~ IG(c0, d0) (half_cauchy =
# FALSE: the plan's choice, conjugate, and with a light tail). `sbc`:
# centered at theta0, and as tight as sane data need: of 300 datasets from
# prior draws (db_prior_predictive()), the share of those with a share below
# 0.02 or above 0.9 fell from 0.23 / 0.05 (MNP, SD 0.5), 0.08 (HMNL) and 0.20
# (HMNP; theta SD 0.5, sigma_d2 ~ IG(4, 0.75)) to 0.10 / 0.007, 0.04 and
# 0.15 with the values below (2026-11-02). `bvm`: centered at 0 (theta0 is
# not the prior's center), SD 2 on every coefficient.
#   MNP   beta ~ N(beta_bar, A^-1); Sigma ~ IW(nu, V), E Sigma = V / (nu - 4).
#         sbc: beta_bar = theta0's, SD 0.35; nu = 10, V = 6 Sigma0 (E =
#         Sigma0, theta0's). bvm: beta_bar = 0, SD 2; nu = 10, V = 6 Sigma_iid
#         (iid level errors: 1 on the diagonal, 0.5 off it).
#   HMNL  b ~ N(b_bar, A^-1); W ~ IW(nu, V), E W = V / (nu - 3); theta ~
#         N(theta_bar, A_theta^-1); sigma_d2 ~ IG(c0, d0). sbc: b_bar and
#         theta_bar theta0's, SDs 0.3; nu = 7, V = 4 W0 (E W = W0, theta0's);
#         IG(6, 1.25) (mean 0.25, SD 0.125). bvm: b_bar = theta_bar = 0, SD
#         2; nu = 5, V = 0.5 I (E W = 0.25 I); IG(4, 0.75) (mean 0.25).
#   HMNP  sbc: HMNL's divided by the probit scale (raw sigma^2 near 1), and
#         sigma^2 ~ IG(10, 9) (mean 1, SD 0.35).
# The draws follow choicer's conventions (src/bayes_samplers.h,
# src/hb_internal.h): A and A_theta are precisions; IW(nu, V) is the inverse
# of a Wishart(nu, V^-1); IG(shape, scale) is scale / Gamma(shape, 1).
#
# Randomness, from the caller's `seed`: the design (D2's or D1-L's simulate()
# at `seed`), the probit shocks (MNP: subseed 1), the prior draw (subseed 2),
# the sampler (subseed 3), GHK's uniforms (subseed 4); db_subseed().
# ============================================================================

DB_VARIANTS <- c("default", "hmnl", "hmnp")
DB_MODEL <- c(default = "mnp", hmnl = "hmnl", hmnp = "hmnp")
DB_ALTS <- c("air", "bus", "car", "train")              # D2's labels and row order
DB_BASE <- "car"                                        # run_mnprobit()'s base_alt
DB_NONBASE <- setdiff(DB_ALTS, DB_BASE)                 # w's order: air, bus, train
DB_X <- c("wait", "travel", "vcost", "income_air")
DB_X_SCALE <- c(wait = 30, travel = 240, vcost = 25, income_air = 20)
DB_P <- 3L                                              # utility differences
DB_BETA <- c(DB_X, paste0("ASC_", DB_NONBASE))
DB_Q <- c(air = 0.14, bus = 0.09, car = 0.64, train = 0.13)   # D2's population shares
DB_HB_X <- c("price", "feature")
DB_HB_ALTS <- sprintf("P%d", 1:6)
DB_T <- 8L                                              # D1-L tasks per person
DB_PROBIT <- pi / sqrt(6)                               # Gumbel SD over the probit's
DB_GHK_DRAWS <- 64L
DB_GH_NODES <- 24L

# Sigma's lower triangle, row by row (Sigma_11, Sigma_21, Sigma_22, ...), as
# run_mnprobit() names its draws, and back.
.db_vech <- function(S) {
  p <- nrow(S)
  i <- rep(seq_len(p), seq_len(p))
  j <- sequence(seq_len(p))
  stats::setNames(S[cbind(i, j)], sprintf("Sigma_%d%d", i, j))
}
.db_unvech <- function(v, p = DB_P) {
  S <- matrix(0, p, p)
  S[cbind(rep(seq_len(p), seq_len(p)), sequence(seq_len(p)))] <- v
  S + t(S) - diag(diag(S), p)
}
.DB_SIGMA_NAMES <- names(.db_vech(diag(DB_P)))

# theta0's level errors (DB_ALTS' order; ground modes correlated 0.5), and
# the iid errors' Sigma / Sigma_11 (the BvM prior's center).
.DB_OMEGA <- matrix(c(1, 0, 0, 0, 0, 1, 0.5, 0.5, 0, 0.5, 1, 0.5, 0, 0.5, 0.5, 1), 4L)
.DB_SIGMA_IID <- .db_unvech(c(1, 0.5, 1, 0.5, 0.5, 1))

DB_THETA0 <- list(
  default = c(wait = -2.1884803716763224, travel = -0.42155687351273091,
              vcost = -0.15500591024560662, income_air = 0.010435638674514278,
              ASC_air = 3.0269231188888872, ASC_bus = 2.0965172480977001,
              ASC_train = 2.2184203567622491, Sigma_11 = 1, Sigma_21 = 0.25, Sigma_22 = 0.5,
              Sigma_31 = 0.25, Sigma_32 = 0.25, Sigma_33 = 0.5),
  hmnl = c(price = -0.8, feature = 0.6, "W[price,price]" = 0.16, "W[feature,price]" = -0.06,
           "W[feature,feature]" = 0.25, "(Intercept)" = 1.2599756978744163, z = 0.6,
           sigma_d2 = 0.25, delta_P1 = 0.20095099512410092, delta_P2 = 1.2699138033474233,
           delta_P3 = 1.576155721719076, delta_P4 = 1.3628114951425077,
           delta_P5 = 2.2567632004099405, delta_P6 = 1.3125785151523788))
DB_THETA0$hmnp <- DB_THETA0$hmnl /
  DB_PROBIT^ifelse(grepl("^W\\[|^sigma_d2$", names(DB_THETA0$hmnl)), 2, 1)

# The HB priors in the samplers' argument names, from the logit scale's
# centers and SDs; s = DB_PROBIT gives HMNP's raw scale.
.db_hb_prior <- function(b_bar, b_sd, nu, V, theta_bar, theta_sd, c0, d0, s = 1) {
  K <- length(b_bar)
  P <- length(theta_bar)
  list(b_bar = b_bar / s, A = diag(s^2 / b_sd^2, K), nu = nu, V = V / s^2,
       theta_bar = theta_bar / s, A_theta = diag(s^2 / theta_sd^2, P),
       sd_prior = list(half_cauchy = FALSE, s_d = 1, c0 = c0, d0 = d0 / s^2))
}
.DB_W0 <- matrix(c(0.16, -0.06, -0.06, 0.25), 2L)
.DB_THETA_HB <- c(1.2599756978744163, 0.6)
DB_PRIOR <- list(
  sbc = list(
    mnp = list(beta_bar = unname(DB_THETA0$default[DB_BETA]), A = diag(1 / 0.35^2, 7L),
               nu = 10, V = 6 * .db_unvech(DB_THETA0$default[.DB_SIGMA_NAMES])),
    hmnl = .db_hb_prior(c(-0.8, 0.6), 0.3, 7, 4 * .DB_W0, .DB_THETA_HB, 0.3, 6, 1.25),
    hmnp = c(.db_hb_prior(c(-0.8, 0.6), 0.3, 7, 4 * .DB_W0, .DB_THETA_HB, 0.3, 6, 1.25,
                          s = DB_PROBIT), list(a0 = 10, s0 = 9))),
  bvm = list(
    mnp = list(beta_bar = rep(0, 7L), A = diag(1 / 4, 7L), nu = 10, V = 6 * .DB_SIGMA_IID),
    hmnl = .db_hb_prior(c(0, 0), 2, 5, 0.5 * diag(2L), c(0, 0), 2, 4, 0.75)))

# A seed derived from `seed` for stream k (integers in [1, 2^31 - 2]).
db_subseed <- function(seed, k) 1 + (as.numeric(seed) + k * 1000003) %% (2^31 - 2)

# --- Prior draws ------------------------------------------------------------

# N(mu, A^-1), A a precision.
.db_rmvn_prec <- function(mu, A) drop(mu + backsolve(chol(A), stats::rnorm(length(mu))))
# IW(nu, V): the inverse of a Wishart(nu, V^-1) draw.
.db_riwish <- function(nu, V) {
  S <- solve(stats::rWishart(1L, nu, solve(V))[, , 1L])
  (S + t(S)) / 2
}
# IG(shape, scale): scale / Gamma(shape, 1).
.db_rinvgamma <- function(shape, scale) scale / stats::rgamma(1L, shape)
# sigma_d2 from an HB sampler's sd_prior: IG(c0, d0), or sigma_d ~
# half-Cauchy(0, s_d), the Makalic-Schmidt mixture's marginal.
.db_rsigma_d2 <- function(sdp) {
  if (isTRUE(sdp$half_cauchy)) (sdp$s_d * abs(stats::rcauchy(1L)))^2
  else .db_rinvgamma(sdp$c0, sdp$d0)
}

# One draw of the model's raw parameters from `prior` (the caller's RNG):
# MNP list(beta, Sigma); HMNL list(b, W, theta, sigma_d2, xi); HMNP the same
# and sigma2. RNG order: beta (b), Sigma (W), theta, sigma_d2, xi, sigma2.
db_prior_draw <- function(model, prior) {
  if (model == "mnp") {
    return(list(beta = stats::setNames(.db_rmvn_prec(prior$beta_bar, prior$A), DB_BETA),
                Sigma = .db_riwish(prior$nu, prior$V)))
  }
  b <- stats::setNames(.db_rmvn_prec(prior$b_bar, prior$A), DB_HB_X)
  W <- .db_riwish(prior$nu, prior$V)
  theta <- .db_rmvn_prec(prior$theta_bar, prior$A_theta)
  sigma_d2 <- .db_rsigma_d2(prior$sd_prior)
  xi <- stats::setNames(stats::rnorm(length(DB_HB_ALTS), 0, sqrt(sigma_d2)), DB_HB_ALTS)
  out <- list(b = b, W = W, theta = theta, sigma_d2 = sigma_d2, xi = xi)
  if (model == "hmnp") out$sigma2 <- .db_rinvgamma(prior$a0, prior$s0)
  out
}

# MNP's identified quantities of a raw (beta, Sigma): beta / sqrt(Sigma_11)
# and the vech of Sigma / Sigma_11, named as run_mnprobit()'s draws.
db_mnp_identified <- function(beta, Sigma) {
  c(beta[DB_BETA] / sqrt(Sigma[1, 1]), .db_vech(Sigma / Sigma[1, 1]))
}

# D1-L's override (WP4) of an HB raw draw: identified (raw / sigma) for HMNP.
.db_hb_override <- function(raw) {
  s <- sqrt(raw$sigma2 %||% 1)
  list(b = raw$b / s, W = raw$W / s^2, theta = raw$theta / s, sigma_d2 = raw$sigma_d2 / s^2,
       xi = raw$xi / s)
}

# --- Data -------------------------------------------------------------------

# MNP's differenced design, choicer's layout: rows (situation i, non-base
# alternative j), i-major, j in air, bus, train; columns DB_X (x_ij -
# x_i,car), then one constant per non-base alternative. `data` holds the 4
# modes of each situation in DB_ALTS' order.
db_mnp_design <- function(data) {
  J <- length(DB_ALTS)
  N <- nrow(data) / J
  if (N != floor(N) || !identical(as.character(data$alt), rep(DB_ALTS, N))) {
    stop("DB: MNP data must hold the 4 modes of every situation in label order")
  }
  nb <- match(DB_NONBASE, DB_ALTS)
  b <- match(DB_BASE, DB_ALTS)
  Xd <- vapply(DB_X, function(x) {
    m <- matrix(data[[x]], N, J, byrow = TRUE)
    as.vector(t(m[, nb, drop = FALSE] - m[, b]))
  }, numeric(N * DB_P))
  Xd <- matrix(Xd, ncol = length(DB_X), dimnames = list(NULL, DB_X))
  asc <- diag(DB_P)[rep(seq_len(DB_P), N), , drop = FALSE]
  colnames(asc) <- paste0("ASC_", DB_NONBASE)
  cbind(Xd, asc)
}

# The alternative chosen at each y (0 = the base), as DB_ALTS' positions.
.db_mnp_alt_of_y <- function(y) match(c(DB_BASE, DB_NONBASE), DB_ALTS)[y + 1L]

# MNP data: D2's travellers at `seed` (their MNL choices discarded), the
# covariates scaled, choices from w = Xd beta + eps, eps ~ N(0, Sigma)
# (subseed 1). Returns list(data, Xd, y): y = 0 (car), 1..3 (air, bus,
# train), run_mnprobit()'s coding at base_alt = "car".
db_mnp_data <- function(N, seed, beta, Sigma) {
  d2 <- get_dgp("D2")$simulate(N, seed, "default")$data
  data <- d2[, c("id", "alt", DB_X, "type"), with = FALSE]
  for (x in DB_X) data.table::set(data, j = x, value = data[[x]] / DB_X_SCALE[[x]])
  Xd <- db_mnp_design(data)
  eps <- vs_with_seed(db_subseed(seed, 1L), matrix(stats::rnorm(N * DB_P), N) %*% chol(Sigma))
  w <- matrix(drop(Xd %*% beta[colnames(Xd)]), N, DB_P, byrow = TRUE) + eps
  y <- ifelse(apply(w, 1L, max) < 0, 0L, max.col(w, ties.method = "first"))
  ch <- as.integer(match(data$alt, DB_ALTS) == rep(.db_mnp_alt_of_y(y), each = length(DB_ALTS)))
  data.table::set(data, j = "choice", value = ch)
  data.table::setcolorder(data, c("id", "alt", "choice"))
  list(data = data, Xd = Xd, y = y)
}

# A dataset of `variant` with N situations (MNP) or persons (HMNL/HMNP, T
# tasks each): list(data, truth, meta). The truth is drawn from `prior` (a
# DB_PRIOR entry; subseed 2) when given, else it is `theta` (MNP: identified,
# named as DB_THETA0$default; HB: WP4's override list) or theta0 (MNP:
# DB_THETA0's, raw at Sigma_11 = 1; HMNL: WP4's; HMNP: WP4's over the probit
# scale). `truth` is identified and named as choicer's draws; meta: variant,
# model, N, T, raw (the raw draw, or NULL), prior, the realized shares, and
# MNP's Xd and y, or the HB's true tastes (N x 2, identified) and delta.
db_simulate <- function(N, seed, variant = "default", prior = NULL, theta = NULL, T = DB_T) {
  variant <- match.arg(variant, DB_VARIANTS)
  model <- DB_MODEL[[variant]]
  N <- as.integer(N)
  raw <- if (!is.null(prior)) vs_with_seed(db_subseed(seed, 2L), db_prior_draw(model, prior))
  if (model == "mnp") {
    if (is.null(raw)) {
      th <- theta %||% DB_THETA0$default
      raw <- list(beta = th[DB_BETA], Sigma = .db_unvech(th[.DB_SIGMA_NAMES]))
    }
    s <- db_mnp_data(N, seed, raw$beta, raw$Sigma)
    sh <- tabulate(.db_mnp_alt_of_y(s$y), length(DB_ALTS)) / N
    return(list(data = s$data, truth = db_mnp_identified(raw$beta, raw$Sigma),
                meta = list(variant = variant, model = model, N = N, T = 1L, raw = raw,
                            prior = prior, shares = stats::setNames(sh, DB_ALTS), Xd = s$Xd,
                            y = s$y)))
  }
  ov <- if (!is.null(raw)) .db_hb_override(raw)
        else theta %||% if (model == "hmnp") .db_hb_override(list(
          b = SPL_TRUE$b, W = SPL_TRUE$W, theta = SPL_TRUE$theta, sigma_d2 = SPL_TRUE$sigma_d2,
          xi = sqrt(SPL_TRUE$sigma_d2) * SPL_XI_STD, sigma2 = DB_PROBIT^2))
  sim <- get_dgp("D1")$simulate(N, seed, "labeled", T = as.integer(T), theta = ov,
                                shock = if (model == "hmnl") "logit" else "probit")
  ch <- sim$data$alt[sim$data$choice == 1L]
  sh <- tabulate(match(ch, DB_HB_ALTS), length(DB_HB_ALTS)) / (N * T)
  list(data = sim$data, truth = sim$truth,
       meta = list(variant = variant, model = model, N = N, T = as.integer(T), raw = raw,
                   prior = prior, shares = c(outside = 1 - sum(sh), stats::setNames(sh, DB_HB_ALTS)),
                   tastes = sim$meta$tastes, delta = sim$meta$hb$delta))
}

# --- Probabilities and log-likelihoods (SBC's test quantity, Modrák et al.
# 2023) --------------------------------------------------------------------------

# MNP: P(y_i | Xd_i, beta, Sigma) by GHK with the uniforms U (N x R x 3;
# fixed, so that the result is a fixed function of (theta, y): SBC's ranks
# need the same U at the truth and at every draw). For y = c the event is
# M w < 0, M = I (c = 0) or rows -e_c and e_k - e_c (k != c); with u =
# M (w - mu) ~ N(0, M Sigma M') and bounds -M mu, GHK draws u's components in
# turn from their truncated conditionals.
db_mnp_prob <- function(beta, Sigma, Xd, y, U) {
  p <- ncol(Sigma)
  N <- length(y)
  R <- dim(U)[2L]
  mu <- matrix(drop(Xd %*% beta[colnames(Xd)]), N, p, byrow = TRUE)
  out <- numeric(N)
  for (c in 0:p) {
    i <- which(y == c)
    if (!length(i)) next
    M <- diag(p)
    if (c > 0) M[, c] <- -1
    C <- t(chol(M %*% Sigma %*% t(M)))
    b <- -mu[i, , drop = FALSE] %*% t(M)                 # n_c x p bounds
    prob <- matrix(1, length(i), R)
    eta <- vector("list", p)
    for (k in seq_len(p)) {
      m <- b[, k]
      for (l in seq_len(k - 1L)) m <- m - C[k, l] * eta[[l]]
      Pk <- stats::pnorm(m / C[k, k])
      prob <- prob * Pk
      if (k < p) eta[[k]] <- stats::qnorm(pmax(matrix(U[i, , k], length(i), R) * Pk, 1e-300))
    }
    out[i] <- rowMeans(prob)
  }
  out
}

db_mnp_loglik <- function(beta, Sigma, Xd, y, U) {
  sum(log(pmax(db_mnp_prob(beta, Sigma, Xd, y, U), 1e-300)))
}

# The HB datasets' layout for the log-likelihoods: X (rows x 2), each row's
# person (1..N in data order) and alternative (1..6), each task's chosen
# alternative (0 = outside). D1-L offers all 6 products in every task.
db_hb_layout <- function(data) {
  J <- length(DB_HB_ALTS)
  n_t <- nrow(data) / J
  if (n_t != floor(n_t) || !identical(as.character(data$alt), rep(DB_HB_ALTS, n_t))) {
    stop("DB: HB data must hold the 6 products of every task in label order")
  }
  ch <- matrix(data$choice, n_t, J, byrow = TRUE)
  list(X = as.matrix(data[, DB_HB_X, with = FALSE]),
       person = match(data$person, unique(data$person)), persons = unique(data$person),
       alt = rep(seq_len(J), n_t), chosen = ifelse(rowSums(ch) == 0, 0L, max.col(ch, "first")),
       J = J, n_tasks = n_t)
}

# Systematic utilities, tasks x 6, given the persons' tastes B (N x 2, rows
# in data order) and delta (6, in DB_HB_ALTS' order).
.db_hb_utility <- function(lay, B, delta) {
  v <- rowSums(lay$X * B[lay$person, , drop = FALSE]) + unname(delta)[lay$alt]
  matrix(v, lay$n_tasks, lay$J, byrow = TRUE)
}

# HMNL: each task's log P(chosen | beta_i, delta), logit with the outside
# option's utility 0; their sum is run_hmnlogit()'s loglik trace.
db_hmnl_logp <- function(lay, B, delta) {
  V <- .db_hb_utility(lay, B, delta)
  m <- pmax(0, apply(V, 1L, max))
  lse <- m + log(exp(-m) + rowSums(exp(V - m)))
  ifelse(lay$chosen > 0, V[cbind(seq_len(lay$n_tasks), pmax(lay$chosen, 1L))], 0) - lse
}
db_hmnl_loglik <- function(lay, B, delta) sum(db_hmnl_logp(lay, B, delta))

# Gauss-Hermite rule for E f(z), z ~ N(0, 1) (Golub-Welsch).
.db_gh <- function(n) {
  jac <- matrix(0, n, n)
  jac[cbind(1:(n - 1L), 2:n)] <- jac[cbind(2:n, 1:(n - 1L))] <- sqrt(1:(n - 1L) / 2)
  eg <- eigen(jac, symmetric = TRUE)
  list(z = sqrt(2) * eg$values, w = eg$vectors[1L, ]^2)
}

# HMNP (identified scale, unit shocks): each task's log P(chosen), P(c) =
# E_e prod_{k != c} Phi(V_c - V_k + e), e ~ N(0, 1), over the 6 products and
# the outside option (V = 0), by Gauss-Hermite; and their sum.
db_hmnp_logp <- function(lay, B, delta, nodes = DB_GH_NODES) {
  gh <- .db_gh(nodes)
  V <- cbind(.db_hb_utility(lay, B, delta), 0)
  ch <- ifelse(lay$chosen > 0, lay$chosen, ncol(V))
  a <- V[cbind(seq_len(lay$n_tasks), ch)] - V                # V_c - V_k
  s <- vapply(gh$z, function(e) {
    rowSums(stats::pnorm(a + e, log.p = TRUE)) - stats::pnorm(e, log.p = TRUE)
  }, numeric(lay$n_tasks))
  s <- s + rep(log(gh$w), each = lay$n_tasks)
  m <- apply(s, 1L, max)
  m + log(rowSums(exp(s - m)))
}
db_hmnp_loglik <- function(lay, B, delta, nodes = DB_GH_NODES) {
  sum(db_hmnp_logp(lay, B, delta, nodes))
}

# --- Fits and draws ---------------------------------------------------------

# The model's sampler on `data` with `prior` and `mcmc` (R, burn, thin,
# seed), one chain.
db_fit <- function(model, data, prior, mcmc, keep_beta_i = "none") {
  switch(model,
         mnp = choicer::run_mnprobit(data, "id", "alt", "choice", DB_X, base_alt = DB_BASE,
                                     prior = prior, mcmc = mcmc, keep_data = FALSE),
         hmnl = choicer::run_hmnlogit(data, "id", "alt", "choice", DB_HB_X,
                                      person_col = "person", alt_covariate_cols = "z",
                                      prior = prior, mcmc = mcmc, keep_beta_i = keep_beta_i),
         hmnp = choicer::run_hmnprobit(data, "id", "alt", "choice", DB_HB_X,
                                       person_col = "person", alt_covariate_cols = "z",
                                       prior = prior, mcmc = mcmc, keep_beta_i = keep_beta_i))
}

# The kept draws of the identified quantities, named as db_simulate()'s
# truth: MNP beta and Sigma without Sigma_11 (1 by normalization); HB b, W,
# theta, sigma_d2 and delta (choicer normalizes HMNP's chains).
db_draws <- function(fit, model) {
  if (model == "mnp") return(cbind(fit$draws$beta, fit$draws$sigma[, -1L, drop = FALSE]))
  ch <- fit$chains[[1L]]
  cbind(ch$b, ch$w_vech, ch$theta, sigma_d2 = ch$sigma_d2,
        `colnames<-`(ch$delta, paste0("delta_", colnames(ch$delta))))
}

# The identified quantities' blocks (the rep record's est$block):
# coefficients "beta", constants and deltas "asc", covariances "sigma",
# theta and sigma_d2 "other".
db_block <- function(nm) {
  ifelse(nm %in% c(DB_X, DB_HB_X), "beta",
         ifelse(grepl("^ASC_|^delta_", nm), "asc",
                ifelse(grepl("^Sigma_|^W\\[", nm), "sigma", "other")))
}

# The persons' taste draws (identified), N x 2 x length(idx), persons in data
# order, of an HB fit run with keep_beta_i = "draws", at kept draws `idx`.
db_beta_i_draws <- function(fit, lay, idx) {
  i <- match(as.character(lay$persons), as.character(fit$data$person_ids))
  if (anyNA(i)) stop("DB: the fit's persons do not match the data's")
  aperm(fit$beta_i$draws[, i, idx, drop = FALSE], c(2L, 1L, 3L))
}

# --- Moments, targets and the prior predictive --------------------------------

# Calibration moments: the realized shares (HB: and the outside option's),
# their minimum and maximum.
db_moments <- function(sim) {
  sh <- sim$meta$shares
  c(stats::setNames(sh, paste0("share_", names(sh))), min_share = min(sh), max_share = max(sh))
}

DB_TARGETS <- data.frame(moment = c("min_share", "max_share"), lo = c(0.02, 0),
                         hi = c(1, 0.9), stringsAsFactors = FALSE)

# The prior predictive spread of the shares: `n` datasets of `variant` from
# DB_PRIOR$sbc (seeds seed + 1..n), each's realized shares. Returns the
# shares (n x alternatives), their 5/50/95% quantiles, and the rates of
# datasets with a share below 0.02 or above 0.9.
db_prior_predictive <- function(variant = "default", n = 200L, N = NULL, T = DB_T, seed = 1L) {
  model <- DB_MODEL[[variant]]
  N <- N %||% if (model == "mnp") 300L else 100L
  sh <- t(vapply(seq_len(n), function(r) {
    db_simulate(N, seed + r, variant, prior = DB_PRIOR$sbc[[model]], T = T)$meta$shares
  }, numeric(if (model == "mnp") 4L else 7L)))
  list(shares = sh, quantiles = apply(sh, 2L, stats::quantile, c(0.05, 0.5, 0.95)),
       extreme_rate = c(below_0.02 = mean(apply(sh, 1L, min) < 0.02),
                        above_0.9 = mean(apply(sh, 1L, max) > 0.9)))
}

# --- Calibration of theta0 (never run at source time) -------------------------

# The population's expected MNP shares at identified theta (DB_ALTS' order):
# D2's 210 traveller types weighted by pi, each type's probabilities by GHK
# with R draws (uniforms from `seed`).
db_mnp_pop_shares <- function(theta, R = 4000L, seed = 1L) {
  ty <- d2_types()
  long <- data.table::copy(ty$long)
  for (x in DB_X) data.table::set(long, j = x, value = long[[x]] / DB_X_SCALE[[x]])
  Xd <- db_mnp_design(long)
  U <- vs_with_seed(seed, array(stats::runif(ty$n * R * DB_P), c(ty$n, R, DB_P)))
  Sigma <- .db_unvech(theta[.DB_SIGMA_NAMES])
  P <- vapply(0:DB_P, function(c) db_mnp_prob(theta[DB_BETA], Sigma, Xd, rep(c, ty$n), U),
              numeric(ty$n))
  s <- drop(ty$pi %*% P)
  stats::setNames(s[match(DB_ALTS, c(DB_BASE, DB_NONBASE))], DB_ALTS)
}

# DB_THETA0$default's construction from D2's MNL theta0 (header): level
# coefficients and constants over pi / sqrt(6), Sigma from Omega, then the
# constants re-solved by ASC_j += (log Q_j - log s_j) - (log Q_car - log
# s_car) until the log-share gaps are below `tol` (common uniforms: a fixed
# point of a smooth map). Returns the theta and the path's last gaps.
db_mnp_calibrate <- function(R = 4000L, seed = 1L, tol = 1e-10, max_iter = 200L) {
  th <- D2_THETA0$default
  a <- c(air = 0, th[paste0("ASC_", c("bus", "car", "train"))])
  names(a) <- DB_ALTS
  M <- diag(4L)[match(DB_NONBASE, DB_ALTS), , drop = FALSE]
  M[, match(DB_BASE, DB_ALTS)] <- -1
  Sigma <- M %*% .DB_OMEGA %*% t(M)
  s <- sqrt(Sigma[1, 1]) * DB_PROBIT                     # level probit scale, then Sigma_11
  theta <- c(th[DB_X] * DB_X_SCALE[DB_X] / s,
             stats::setNames((a[DB_NONBASE] - a[[DB_BASE]]) / s, paste0("ASC_", DB_NONBASE)),
             .db_vech(Sigma / Sigma[1, 1]))
  asc <- paste0("ASC_", DB_NONBASE)
  for (it in seq_len(max_iter)) {
    r <- log(DB_Q[DB_ALTS]) - log(db_mnp_pop_shares(theta, R, seed))
    gap <- r[DB_NONBASE] - r[[DB_BASE]]
    if (max(abs(gap)) < tol) break
    theta[asc] <- theta[asc] + gap
  }
  list(theta = theta, gap = gap, iterations = it)
}

# --- Q1.4: choicer's MNP against bayesm::rmnpGibbs ---------------------------

# Q1.4 (plan §5): choicer's run_mnprobit() and bayesm's rmnpGibbs() on the
# same MNP data (DB at theta0, N situations) and prior (DB_PRIOR$bvm$mnp),
# each `burn` + `kept` iterations. bayesm differences against its last
# alternative, so its alternatives are ordered air, bus, train, car (car,
# choicer's base_alt here, last), and its design is built by createX() from
# the wide covariates, independently of choicer's prep; its draws are
# normalized as choicer's (beta / sqrt(Sigma_11), Sigma / Sigma_11). Check
# rows (test, check, value, threshold, pass): the two designs equal; per
# identified quantity, |mean_c - mean_b| / sqrt(mcse_c^2 + mcse_b^2) < 4
# (MCSE = SD / sqrt(bulk ESS), choicer::mcse()) and |sd_c / sd_b - 1| < 0.10;
# a last row (pass NA) with the seconds and the smallest ESS.
ref_q14 <- function(N = 300L, burn = 5000L, kept = 200000L, seed = 20261102L) {
  vs_require("bayesm")
  sim <- db_simulate(N, seed, "default")
  prior <- DB_PRIOR$bvm$mnp
  t0 <- proc.time()[["elapsed"]]
  fc <- suppressMessages(db_fit("mnp", sim$data, prior,
                                list(R = burn + kept, burn = burn, seed = db_subseed(seed, 3L))))
  t_c <- proc.time()[["elapsed"]] - t0
  ord <- c(DB_NONBASE, DB_BASE)
  wide <- do.call(cbind, lapply(DB_X, function(x) {
    matrix(sim$data[[x]], N, length(DB_ALTS), byrow = TRUE)[, match(ord, DB_ALTS)]
  }))
  Xb <- bayesm::createX(p = 4, na = length(DB_X), nd = NULL, Xa = wide, Xd = NULL, INT = TRUE,
                        DIFF = TRUE, base = 4)
  Xb <- Xb[, c(seq_along(DB_X) + DB_P, seq_len(DB_P))]  # to choicer's order: DB_X, constants
  dx <- max(abs(Xb - sim$meta$Xd))
  yb <- ifelse(sim$meta$y == 0L, 4L, sim$meta$y)
  pb <- list(betabar = prior$beta_bar, A = prior$A, nu = prior$nu, V = prior$V)
  t0 <- proc.time()[["elapsed"]]
  utils::capture.output(fb <- vs_with_seed(db_subseed(seed, 5L), bayesm::rmnpGibbs(
    Data = list(p = 4, y = yb, X = Xb), Prior = pb,
    Mcmc = list(R = burn + kept, keep = 1, nprint = 0))))
  t_b <- proc.time()[["elapsed"]] - t0
  keep <- burn + seq_len(kept)
  s11 <- fb$sigmadraw[keep, 1L]
  vi <- rep(seq_len(DB_P), seq_len(DB_P))
  vj <- sequence(seq_len(DB_P))
  Sb <- fb$sigmadraw[keep, (vj - 1L) * DB_P + vi, drop = FALSE] / s11
  db <- cbind(`colnames<-`(fb$betadraw[keep, , drop = FALSE] / sqrt(s11), DB_BETA),
              `colnames<-`(Sb, .DB_SIGMA_NAMES))
  dc <- db_draws(fc, "mnp")
  db <- db[, colnames(dc)]
  # MCSEs from the chains thinned to at most 60,000 draws (conservative for
  # the full chains' means): choicer::ess() returns NA past about 65,000
  # draws (an integer overflow in .fft_autocov()).
  th <- seq(1L, kept, by = ceiling(kept / 60000))
  mc <- choicer::mcse(dc[th, , drop = FALSE])
  mb <- choicer::mcse(db[th, , drop = FALSE])
  z <- abs(colMeans(dc) - colMeans(db)) / sqrt(mc^2 + mb^2)
  r <- abs(apply(dc, 2L, stats::sd) / apply(db, 2L, stats::sd) - 1)
  q <- colnames(dc)
  e <- min(choicer::ess(dc[th, , drop = FALSE])[, "bulk"],
           choicer::ess(db[th, , drop = FALSE])[, "bulk"])
  rbind(
    data.frame(test = "Q1.4", check = "differenced design: bayesm::createX() vs the DGP's",
               value = dx, threshold = 1e-12, pass = dx < 1e-12, stringsAsFactors = FALSE),
    data.frame(test = "Q1.4", check = paste("mean", q, "|diff| / combined MCSE"),
               value = unname(z), threshold = 4, pass = unname(!is.na(z) & z < 4),
               stringsAsFactors = FALSE),
    data.frame(test = "Q1.4", check = paste("sd", q, "|choicer / bayesm - 1|"),
               value = unname(r), threshold = 0.10, pass = unname(!is.na(r) & r < 0.10),
               stringsAsFactors = FALSE),
    data.frame(test = "Q1.4",
               check = sprintf("seconds: choicer %.1f, bayesm %.1f; smallest bulk ESS %.0f",
                               t_c, t_b, e),
               value = t_c + t_b, threshold = NA_real_, pass = NA, stringsAsFactors = FALSE))
}

# --- Q4: the mini-SBC ----------------------------------------------------------

# Q4's check rows from the records of F9's quick cells (q_sbc_mnp,
# q_sbc_hmnl, q_sbc_hmnp; R$quick reps each, through run_rep()): their
# sbc_ecdf rows, gated by F9's gates_quick and Holm-adjusted together at
# alpha = 0.001 (plan Q4: Holm p > 0.001). One row per quantity, value the
# Holm-adjusted p. Returns list(rows, tab).
q4_check <- function(fam, records, alpha = 0.001) {
  key <- vapply(records, function(r) paste(r$cell, r$method, sep = "\r"), "")
  tab <- do.call(rbind, lapply(unname(split(records, key)), metrics_cell,
                               gates = fam$gates_quick, tol_set = "quick"))
  tab <- gate_family(tab, alpha = alpha)
  g <- tab[tab$gated & tab$metric == "sbc_ecdf", , drop = FALSE]
  rows <- data.frame(test = "Q4", check = paste0(g$cell, " ", g$estimand, " (R = ", g$R, ")"),
                     value = g$p_holm, threshold = alpha, pass = g$verdict == "pass",
                     stringsAsFactors = FALSE)
  list(rows = rows, tab = tab)
}

dgp_register(list(
  id = "DB",
  name = paste("Bayesian samplers' data: MNP on D2's design (default), HMNL and HMNP on D1-L;",
               "prior draws for SBC"),
  theta0 = DB_THETA0,
  simulate = function(N, seed, variant = "default", ...) db_simulate(N, seed, variant, ...),
  moments = db_moments, targets = DB_TARGETS, variants = DB_VARIANTS))
