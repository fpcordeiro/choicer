# ============================================================================
# The geographic generator's R reference (plan §4: D3 hospital and D4
# school share one generator; WP6 ports it to C++, lib/geo_gen.cpp, which
# must give identical data at small sizes). This file holds D4, school choice
# from a census (Abdulkadiroglu, Agarwal & Pathak 2017), on two designs:
# D4-40 (F2 and F3) and D4-small; WP6 adds D3. Definitions only, except the
# registration of DGPs "D4" and "D4-40" (dgp_register()), which builds each
# design's fixed world of schools once (milliseconds).
#
# Structure (what the port reproduces):
#   1. geo_draw(): all randomness, as uniforms, in one documented order.
#      The world (schools) is drawn from the design's fixed world_seed, so
#      that every replication faces the same schools and the same truth; the
#      sample (students) from the replication's seed. World, in this order:
#      cluster centres (x, y per cluster, market-major), schools (component,
#      x, y per school, market-major), the sector keys, the quality noise and
#      the constants' noise xi (one per school each). Sample, in this order:
#      students' locations (component, x, y per student), SES noise, tastes
#      (K_w = 2 per student, drawn in every variant so that variants share a
#      seed's students), then the uniforms that pick the choices (one per
#      student). Normals are qnorm() of these uniforms.
#   2. geo_build(): deterministic construction: cluster centres, locations
#      (a mixture of a uniform background and truncated normals around the
#      market's clusters, by inversion), sectors (a permutation of the
#      market's sector counts by the keys' ranks), school quality, student
#      SES, and the long table: each student's K nearest schools of their
#      market (ties by school index), rows in choicer's order (student id,
#      then school label), with distance and the interactions.
#   3. Probabilities from lib/simulate.R's engines at theta0 and choices by
#      inverse CDF (draw_choice()) from the sample's choice uniforms.
# Port notes: distances are sqrt(dx * dx + dy * dy) in double and the
# interactions single products, so a C++ port must not contract them into
# fused multiply-adds; qnorm()/pnorm() are R's (Rmath).
#
# D4 (school). M markets, each a square of side `side` km whose population
# clusters around `clusters` centres (a share `background` spread
# uniformly); schools are located by the same density. Each market's schools
# are split into sectors (public, charter, private) by fixed counts. Student
# i's choice set is the K nearest schools of their market plus an outside
# option (no rows; utility 0). SES_i = rho_ses g(x_i) + sqrt(1 - rho_ses^2)
# z_i, with g the market's standardized east-west gradient, so that SES
# sorts spatially. Quality q_j = rho_q g(x_j) + sector_q + tau_q z_j. Utility
# of school j for student i:
#   MNL / NL: delta_j + b_d d_ij + b_ds d_ij SES_i + b_qs q_j SES_i
#   MXL:      the same, plus exp(s_d eta_i1) (-d_ij) + s_q eta_i2 q_j,
# with constants delta_j = a q_j + sector_delta + sigma_xi xi_j (skewed
# shares). Covariates: dist, dist_ses, q_ses (X); MXL also ndist = -dist
# (log-normal) and quality (normal) in W. NL: nests by sector (RU2), the
# outside option its own degenerate nest.
#
# The MXL parameterization (resolving choicer's global rc_mean, which gives
# every random column a mean or none, R/mxlogit_utils.R:1253). With a
# constant for every school, a mean on quality, a school-level attribute, is
# collinear with the constants, so rc_mean = TRUE is not identified. The
# design uses rc_mean = FALSE, uncorrelated:
#   - quality's normal term is mean-zero, s_q eta_i2 q_j; quality's mean
#     effect is in the constants (a q_j);
#   - distance's coefficient is b_d - exp(s_d eta_i1): the fixed b_d (dist in
#     X) plus choicer's log-normal exp((L eta)_1) on ndist = -dist in W. It
#     is the shifted log-normal of rc_mean = TRUE with its location free in
#     sign (median b_d - 1 per km; exp(g) has median 1).
# Exactly choicer's model (rc_dist = c(1, 0), rc_mean = FALSE,
# rc_correlation = FALSE) and identified: b_d by the mean distance gradient,
# s_d by its dispersion, s_q by how substitution among schools of different
# quality varies with the composition of the nearest-K sets.
#
# Designs. D4-40 (DGP "D4-40", F2's and F3's): one market of 40 schools,
# K = 15, so that choice sets still vary by location (each student's 15
# nearest of 40) while J = 40 constants keep J/sqrt(N) at 1.26, 0.63, 0.32
# on the 1k/4k/16k rungs (S1, 2026-11-02: with D4-small's 160 the distance
# coefficient's bias/SD was -0.163 J/sqrt(N), the incidental constants'
# bias C6, which the gates are not meant to test). One market of the same
# geography rather than four of 10 schools: with 10 schools and K = 15 every
# student would face the whole market, and sigma_q, identified by how the
# nearest-K sets' composition varies, would not be. D4-small (DGP "D4", 4
# markets x 40, 160 constants) stays for WP6's scaled designs.
#
# theta0. The structural values below are literals chosen 2026-11-02 and
# shared by the designs. The tastes' dispersion (s_d = 0.8, s_q = 1.75) and
# quality's spatial sorting (rho_q = 0.8, tau_q = 0.6) were raised from s_d =
# 0.5, s_q = 0.8, rho_q = 0.4 so that sigma_q / SE(sigma_q) >= 5 at N = 1000
# (plan §4) on D4-small; SES's sorting was then lowered from rho_ses = 0.5
# to 0.45 so that D4-40's NL variant meets the quality-gradient target.
# sigma_q / SE(sigma_q) at N = 1000 (checks/wp5_geo.R --only=ident; a 16k
# fit's SE times 4, and the median of five 1k fits): D4-40 6.7 and 6.6,
# D4-small 5.5 and 5.1. Each design's sector_delta (per variant) was
# calibrated the same day by d4_calibrate(), the expected shares at N =
# 40000, seed 1, set to outside .15, public .65, charter .09, private .11:
# plan §4's outside share, and charter and private shares somewhat above the
# national split in NCES's Digest of Education Statistics (about 7% of
# public school students in charter schools and 9% of K-12 students in
# private schools, fall 2021), as in the urban choice districts these
# designs stand for; the outside option takes students out of the market
# (other districts, home schooling). The constants follow from the fixed
# world and these literals (d4_theta0()).
#
# Alternatives never chosen in a replication (plan §4): drop_unbounded_alts()
# drops them, and the situations offering an alternative chosen wherever
# offered, until stable (ML's limit as delta_j -> -Inf or +Inf); families
# record the counts, and geo_truth() restates theta0 on the alternatives
# left.
# ============================================================================

# --- Designs and the truth -----------------------------------------------------

# A design: the geography (markets of `schools` schools, K nearest in the
# student's market, a square of side `side` km, `clusters` population
# centres), the sector counts per market, the world's seed, and the sector
# constants calibrated for it, per variant (d4_calibrate()).
GEO_D4_SMALL <- list(
  kind = "school", name = "D4-small", markets = 4L, schools = 40L, K = 15L, side = 10,
  clusters = 3L, cluster_sd = 1.2, background = 0.25,
  sectors = c(public = 24L, charter = 8L, private = 8L), world_seed = 20261102L,
  sector_delta = list(default = c(public = 2.0613, charter = 0.8518, private = 0.9528),
                      nl = c(public = 2.2155, charter = 0.6877, private = 0.7675),
                      mnl = c(public = 1.4733, charter = 0.5430, private = 0.6660)))

GEO_D4_40 <- list(
  kind = "school", name = "D4-40", markets = 1L, schools = 40L, K = 15L, side = 10,
  clusters = 3L, cluster_sd = 1.2, background = 0.25,
  sectors = c(public = 24L, charter = 8L, private = 8L), world_seed = 40261102L,
  sector_delta = list(default = c(public = 2.0587, charter = 0.9423, private = 1.7015),
                      nl = c(public = 2.2358, charter = 0.6381, private = 1.2110),
                      mnl = c(public = 1.4976, charter = 0.5082, private = 1.1507)))

D4_X <- c("dist", "dist_ses", "q_ses")
D4_W <- c("ndist", "quality")
D4_VARIANTS <- c(default = "mxl", nl = "nl", mnl = "mnl")

# The structural truth, shared by the designs; d4_par() adds a design's
# sector constants.
D4_TRUE <- list(
  rho_ses = 0.45, rho_q = 0.8, tau_q = 0.6,
  sector_q = c(public = 0, charter = 0.2, private = 0.4),
  a = 0.5, sigma_xi = 0.6,
  beta = list(default = c(dist = -0.3, dist_ses = 0.2, q_ses = 0.4),
              nl = c(dist = -1.3, dist_ses = 0.2, q_ses = 0.4),
              mnl = c(dist = -1.3, dist_ses = 0.2, q_ses = 0.4)),
  s = c(ndist = 0.8, quality = 1.75),          # sd of g = L eta (L diagonal)
  lambda = c(public = 0.5, charter = 0.7, private = 0.85))

d4_par <- function(design) c(D4_TRUE, list(sector_delta = design$sector_delta))

# The engines' specification of a D4 variant over the alternatives `alts`
# (labels; their sectors from the schools table).
d4_spec <- function(variant, schools, alts = schools$alt) {
  model <- D4_VARIANTS[[variant]]
  spec <- list(beta = D4_X, alts = alts, outside = TRUE)
  if (model == "nl") spec$nests <- schools$sector[match(alts, schools$alt)]
  if (model == "mxl") {
    spec <- c(spec, list(random = D4_W, rc_dist = c(1L, 0L), rc_mean = FALSE,
                         correlated = FALSE))
  }
  spec
}

# A variant's named blocks (pack_theta()'s): beta, NL lambda by nest label,
# MXL L, and the constants by school label, delta_j = a q_j + sector_delta +
# sigma_xi xi_j.
d4_blocks <- function(variant, schools, par) {
  model <- D4_VARIANTS[[variant]]
  delta <- par$a * schools$quality + par$sector_delta[[variant]][schools$sector] +
    par$sigma_xi * schools$xi
  b <- list(beta = par$beta[[variant]], asc = stats::setNames(delta, schools$alt))
  if (model == "nl") b$lambda <- par$lambda
  if (model == "mxl") b$L <- diag(unname(par$s[D4_W]), 2L)
  b
}

# theta0 of a variant in choicer's names and order.
d4_theta0 <- function(variant, schools, par) {
  pack_theta(D4_VARIANTS[[variant]], d4_spec(variant, schools), d4_blocks(variant, schools, par))
}

# --- 1. Draws ------------------------------------------------------------------------

# Every uniform a dataset uses, in the order of the header: `world` from
# design$world_seed, `sample` from `seed` (N students, K_w tastes each).
geo_draw <- function(design, N, seed, K_w = 2L) {
  N <- as.integer(N)
  J <- design$markets * design$schools
  world <- vs_with_seed(design$world_seed, {
    centers <- matrix(stats::runif(2L * design$clusters * design$markets), 2L)
    school <- matrix(stats::runif(3L * J), 3L)
    sector <- stats::runif(J)
    quality <- stats::runif(J)
    xi <- stats::runif(J)
    list(centers = centers, school = school, sector = sector, quality = quality, xi = xi)
  })
  sample <- vs_with_seed(seed, {
    loc <- matrix(stats::runif(3L * N), 3L)
    ses <- stats::runif(N)
    eta <- matrix(stats::runif(K_w * N), K_w)
    choice <- stats::runif(N)
    list(loc = loc, ses = ses, eta = eta, choice = choice)
  })
  list(world = world, sample = sample, N = N, K_w = as.integer(K_w))
}

# --- 2. Construction -------------------------------------------------------------------

# Students per market: N split as evenly as possible, the first markets
# taking the remainder.
.geo_market_sizes <- function(N, M) (N %/% M) + (seq_len(M) <= N %% M)

# A truncated normal on [lo, hi] by inversion of the uniform u.
.geo_tnorm <- function(u, m, s, lo, hi) {
  a <- stats::pnorm((lo - m) / s)
  b <- stats::pnorm((hi - m) / s)
  m + s * stats::qnorm(a + u * (b - a))
}

# Locations from (component, x, y) uniforms (3 x n) in markets `market`:
# uniform on the square with probability `background`, else a truncated
# normal around one of the market's cluster centres (`ctr`, 2 x C*M).
.geo_locate <- function(u, market, ctr, design) {
  C <- design$clusters
  p0 <- design$background
  bg <- u[1L, ] < p0
  k <- pmin(C, 1L + floor((u[1L, ] - p0) / (1 - p0) * C))
  col <- (market - 1L) * C + ifelse(bg, 1L, k)
  side <- design$side
  sd <- design$cluster_sd
  x <- ifelse(bg, side * u[2L, ], .geo_tnorm(u[2L, ], ctr[1L, col], sd, 0, side))
  y <- ifelse(bg, side * u[3L, ], .geo_tnorm(u[3L, ], ctr[2L, col], sd, 0, side))
  list(x = x, y = y)
}

# The market's standardized east-west gradient (sd 1 for a uniform x).
.geo_gradient <- function(x, side) (x - side / 2) / (side / sqrt(12))

# School labels that sort the same in every locale: "S001", ...
.geo_labels <- function(J) sprintf("S%0*d", max(3L, nchar(J)), seq_len(J))

# Each student's K nearest schools of their market, ties by school index:
# a data.table (id, j, dist) in (id, j) order, built market by market.
geo_sets <- function(students, schools, K) {
  out <- lapply(sort(unique(students$market)), function(m) {
    st <- students[students$market == m, ]
    sc <- schools[schools$market == m, ]
    n <- nrow(st)
    J <- nrow(sc)
    if (J < K) stop("market ", m, " has ", J, " schools, fewer than K = ", K, call. = FALSE)
    dx <- rep(st$x, each = J) - rep(sc$x, n)
    dy <- rep(st$y, each = J) - rep(sc$y, n)
    long <- data.table::data.table(id = rep(st$id, each = J), j = rep(sc$j, n),
                                   dist = sqrt(dx * dx + dy * dy))
    data.table::setorderv(long, c("id", "dist", "j"))
    long <- long[data.table::rowidv(long, "id") <= K]
    data.table::setorderv(long, c("id", "j"))
    long
  })
  data.table::rbindlist(out)
}

# Cluster centres (2 x C*M), in the middle 60% of each market's square.
.geo_centers <- function(design, world) design$side * (0.2 + 0.6 * world$centers)

# The world's schools from its draws: j, alt (label), market, sector, x, y,
# quality and xi (the constants' noise, standard normal).
geo_schools <- function(design, world, par = D4_TRUE) {
  if (design$kind != "school") stop("geo_schools(): only D4 (kind \"school\") so far", call. = FALSE)
  M <- design$markets
  J <- M * design$schools
  sm <- rep(seq_len(M), each = design$schools)
  loc <- .geo_locate(world$school, sm, .geo_centers(design, world), design)
  key_rank <- stats::ave(world$sector, sm, FUN = rank)
  sector <- rep(names(design$sectors), design$sectors)[key_rank]
  quality <- par$rho_q * .geo_gradient(loc$x, design$side) + par$sector_q[sector] +
    par$tau_q * stats::qnorm(world$quality)
  data.table::data.table(j = seq_len(J), alt = .geo_labels(J), market = sm, sector = sector,
                         x = loc$x, y = loc$y, quality = unname(quality),
                         xi = stats::qnorm(world$xi))
}

# The world and the sample from the draws: `schools` (j, alt, market,
# sector, x, y, quality, xi), `students` (id, market, x, y, ses; ids 1..N,
# market by market), `long` (id, alt, market, sector, dist, dist_ses, q_ses,
# ndist, quality) in choicer's row order, and `eta`, the tastes' normals
# (K_w x N).
geo_build <- function(design, draws, par = D4_TRUE) {
  s <- draws$sample
  M <- design$markets
  side <- design$side
  schools <- geo_schools(design, draws$world, par)
  # Students
  N <- draws$N
  stm <- rep(seq_len(M), .geo_market_sizes(N, M))
  loc <- .geo_locate(s$loc, stm, .geo_centers(design, draws$world), design)
  ses <- par$rho_ses * .geo_gradient(loc$x, side) + sqrt(1 - par$rho_ses^2) * stats::qnorm(s$ses)
  students <- data.table::data.table(id = seq_len(N), market = stm, x = loc$x, y = loc$y,
                                     ses = ses)
  # Choice sets and covariates
  long <- geo_sets(students, schools, design$K)
  ses_r <- students$ses[long$id]
  q_r <- schools$quality[long$j]
  long[, `:=`(alt = schools$alt[j], market = schools$market[j], sector = schools$sector[j],
              dist_ses = dist * ses_r, q_ses = q_r * ses_r, ndist = -dist, quality = q_r)]
  long[, j := NULL]
  data.table::setcolorder(long, c("id", "alt", "market", "sector", "dist", "dist_ses", "q_ses",
                                  "ndist", "quality"))
  list(schools = schools, students = students, long = long,
       eta = matrix(stats::qnorm(s$eta), nrow(s$eta)))
}

# --- 3. Probabilities and choices --------------------------------------------------------

# The probabilities choices are drawn from, at theta: MNL/NL the model's;
# MXL the logit given each student's tastes (rc_coefs() at the drawn eta).
.geo_probs <- function(model, d, par, eta) {
  switch(model,
         mnl = mnl_probs(d, par),
         nl = nl_probs(d, par),
         mxl = mxl_cond_probs(d, par, t(rc_coefs(par$mu, par$L, eta, par$rc_dist, par$rc_mean))))
}

# Choices by inverse CDF from the sample's choice uniforms.
.geo_choose <- function(pr, sit, sample) draw_choice(pr$p, pr$p0, sit, u = sample$choice)

# A D4 dataset: list(data, truth, meta). `data` the long table with
# `choice` (0/1); `truth` theta0 over every school; `meta`: variant, model,
# spec and blocks (for geo_truth()), N, the design's summary, the expected
# shares (mean probability per school label and of the outside option) and
# counts, the schools and students tables, the tastes' normals `eta` (K_w x
# N; MXL), and with keep_probs the probabilities the choices were drawn from
# (`probs`: p per row, p0 per student).
d4_simulate <- function(N, seed, variant = "default", design = GEO_D4_SMALL, par = d4_par(design),
                        keep_probs = FALSE) {
  if (!variant %in% names(D4_VARIANTS)) {
    stop("D4 has no variant '", variant, "'; variants: ", paste(names(D4_VARIANTS), collapse = ", "),
         call. = FALSE)
  }
  model <- D4_VARIANTS[[variant]]
  draws <- geo_draw(design, N, seed)
  b <- geo_build(design, draws, par)
  spec <- d4_spec(variant, b$schools)
  blocks <- d4_blocks(variant, b$schools, par)
  truth <- pack_theta(model, spec, blocks)
  d <- design(b$long, spec)
  pr <- .geo_probs(model, d, unpack_theta(model, spec, truth), b$eta)
  data <- b$long
  data[, choice := .geo_choose(pr, d$sit, draws$sample)$choice]
  expected <- stats::setNames(numeric(nrow(b$schools)), b$schools$alt)
  e <- rowsum(pr$p, data$alt)
  expected[rownames(e)] <- e[, 1L]
  meta <- list(variant = variant, model = model, spec = spec, blocks = blocks, N = draws$N,
               design = design[c("name", "markets", "schools", "K", "side")],
               J = nrow(b$schools), expected = expected, outside_expected = sum(pr$p0),
               shares = expected / draws$N, outside_share = sum(pr$p0) / draws$N,
               schools = b$schools, students = b$students,
               eta = if (model == "mxl") b$eta, probs = if (keep_probs) pr)
  list(data = data, truth = truth, meta = meta)
}

# theta0 restated on the alternatives `alts` left after a drop: the
# constants of those alternatives, and the NL's dissimilarities of the nests
# still holding two or more of them (a singleton nest's lambda leaves the
# model: lambda log I = V there).
geo_truth <- function(meta, alts) {
  spec <- meta$spec
  keep <- spec$alts %in% alts
  spec$alts <- spec$alts[keep]
  if (!is.null(spec$nests)) spec$nests <- spec$nests[keep]
  b <- meta$blocks
  if (!is.null(b$lambda)) {
    n <- table(spec$nests)
    b$lambda <- b$lambda[names(b$lambda) %in% names(n)[n > 1L]]
  }
  b$asc <- b$asc[as.character(spec$alts)]
  pack_theta(meta$model, spec, b)
}

# --- Alternatives without a finite constant -------------------------------------------------

# Drops the alternatives never chosen and the situations offering an
# alternative chosen wherever offered, repeating until neither is left (each
# drop can create the other). With an outside option, a situation left
# without rows chose it and is gone too. Generic: any long table with a
# situation id, an alternative label and a 0/1 choice column. Returns
# list(data (rows in their order), dropped_alts (labels), dropped_ids
# (situations), counts = c(alts, situations, rows, rounds)).
drop_unbounded_alts <- function(data, id_col = "id", alt_col = "alt", choice_col = "choice") {
  dt <- data.table::as.data.table(data)
  n0 <- nrow(dt)
  ids0 <- unique(dt[[id_col]])
  alts0 <- unique(dt[[alt_col]])
  rounds <- 0L
  repeat {
    st <- data.table::data.table(a = dt[[alt_col]], y = dt[[choice_col]])[
      , .(n = .N, ch = sum(y)), by = a]
    never <- st$a[st$ch == 0]
    always <- st$a[st$ch == st$n & st$ch > 0]
    if (!length(never) && !length(always)) break
    rounds <- rounds + 1L
    drop_id <- unique(dt[[id_col]][dt[[alt_col]] %in% always])
    dt <- dt[!(dt[[alt_col]] %in% never) & !(dt[[id_col]] %in% drop_id)]
  }
  dropped_alts <- setdiff(alts0, unique(dt[[alt_col]]))
  dropped_ids <- setdiff(ids0, unique(dt[[id_col]]))
  list(data = dt, dropped_alts = dropped_alts, dropped_ids = dropped_ids,
       counts = c(alts = length(dropped_alts), situations = length(dropped_ids),
                  rows = n0 - nrow(dt), rounds = rounds))
}

# A simulated dataset (list(data, truth, meta)) after drop_unbounded_alts():
# the data, theta0 restated on the alternatives left (geo_truth()), and in
# meta the counts dropped (`drop`) and the dropped labels.
geo_drop <- function(sim) {
  dr <- drop_unbounded_alts(sim$data)
  sim$data <- dr$data
  sim$truth <- geo_truth(sim$meta, unique(dr$data$alt))
  sim$meta$drop <- dr$counts
  sim$meta$dropped_alts <- dr$dropped_alts
  sim
}

# --- Calibration moments -----------------------------------------------------------------

# D4's moments on a simulated dataset: outside and sector shares (of all
# students), the quality gradient (OLS slope of the chosen school's quality
# on SES, inside choosers), and untargeted: the share choosing their nearest
# school and the mean distance travelled (inside choosers), the skew of
# school shares (Gini, max/median), the expected count per school per 1,000
# students (min, 10th percentile, median), the number of schools J and J/sqrt(1000).
d4_moments <- function(sim) {
  dt <- sim$data
  N <- sim$meta$N
  ch <- dt[choice == 1L]
  sec <- table(factor(ch$sector, levels = c("public", "charter", "private")))
  ses <- sim$meta$students$ses[ch$id]
  nearest <- dt[, .(dmin = min(dist)), by = id]
  e1k <- sim$meta$expected / N * 1000
  sh <- sort(sim$meta$expected / sum(sim$meta$expected))
  J <- length(sh)
  c(outside_share = 1 - nrow(ch) / N,
    share_public = sec[["public"]] / N, share_charter = sec[["charter"]] / N,
    share_private = sec[["private"]] / N,
    quality_gradient = stats::cov(ses, ch$quality) / stats::var(ses),
    nearest_share = mean(ch$dist == nearest$dmin[match(ch$id, nearest$id)]),
    mean_dist = mean(ch$dist),
    share_gini = sum((2 * seq_len(J) - J - 1) * sh) / J,
    share_max_med = max(sh) / stats::median(sh),
    count_min_1k = min(e1k), count_p10_1k = unname(stats::quantile(e1k, 0.1)),
    count_med_1k = stats::median(e1k),
    J = J, J_over_sqrt_1k = J / sqrt(1000))
}

D4_TARGETS <- data.frame(
  moment = c("outside_share", "share_public", "share_charter", "share_private",
             "quality_gradient"),
  lo = c(0.10, 0.55, 0.05, 0.07, 0.15),
  hi = c(0.20, 0.75, 0.15, 0.15, 0.40),
  stringsAsFactors = FALSE)

# Calibrates sector_delta for a variant: expected shares (mean probabilities
# over N students of one seed, common draws across iterations) to `target`
# (outside, public, charter, private), by the contraction sector_delta_s +=
# log(target_s) - log(share_s), the outside option's utility fixed at 0.
# Not run when the file is sourced; its result is frozen in the design's
# sector_delta.
d4_calibrate <- function(variant, N = 40000L, seed = 1L,
                         target = c(outside = 0.15, public = 0.65, charter = 0.09,
                                    private = 0.11),
                         design = GEO_D4_SMALL, par = d4_par(design), tol = 1e-10,
                         max_iter = 400L) {
  model <- D4_VARIANTS[[variant]]
  draws <- geo_draw(design, N, seed)
  b <- geo_build(design, draws, par)
  spec <- d4_spec(variant, b$schools)
  d <- design(b$long, spec)
  sec <- b$long$sector
  for (it in seq_len(max_iter)) {
    th <- pack_theta(model, spec, d4_blocks(variant, b$schools, par))
    pr <- .geo_probs(model, d, unpack_theta(model, spec, th), b$eta)
    share <- c(outside = sum(pr$p0), tapply(pr$p, sec, sum)) / N
    step <- log(target[names(par$sector_delta[[variant]])]) -
      log(share[names(par$sector_delta[[variant]])])
    par$sector_delta[[variant]] <- par$sector_delta[[variant]] + step
    if (max(abs(step)) < tol) break
  }
  list(sector_delta = par$sector_delta[[variant]], shares = share, iterations = it)
}

# --- Registration --------------------------------------------------------------------------

# Registers a D4 design as a DGP: theta0 per variant from the design's
# fixed world (milliseconds), simulate() on the design.
d4_register <- function(id, design, name) {
  world <- geo_schools(design, geo_draw(design, 0L, 1L)$world)
  par <- d4_par(design)
  dgp_register(list(
    id = id, name = name, design = design,
    theta0 = lapply(stats::setNames(nm = names(D4_VARIANTS)), d4_theta0, schools = world,
                    par = par),
    simulate = function(N, seed, variant = "default", ...) {
      d4_simulate(N, seed, variant, design = design, ...)
    },
    moments = d4_moments, targets = D4_TARGETS, variants = names(D4_VARIANTS)))
}

d4_register("D4", GEO_D4_SMALL,
            "school choice, D4-small (4 markets x 40 schools, K = 15, outside option)")
d4_register("D4-40", GEO_D4_40,
            "school choice, D4-40 (1 market of 40 schools, K = 15, outside option)")
