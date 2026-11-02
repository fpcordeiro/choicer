# ============================================================================
# The geographic generator's R reference (plan §4: D3 hospital and D4
# school share one generator). This file holds D4, school choice from a
# census (Abdulkadiroglu, Agarwal & Pathak 2017), on three designs (D4-40
# for F2 and F3, D4-small, and the census of F7), and D3, hospital choice
# from claims (Capps, Dranove & Satterthwaite 2003; Ho 2006; Gowrisankaran,
# Nevo & Town 2015), on two (F6's 200 hospitals and F7's claims).
# lib/geo_gen.cpp ports the parts whose cost grows with the long table (the
# choice sets and their rows, the logit probabilities and the choices) to
# C++ for population scale: simulate(engine = "cpp") uses it, engine = "r"
# (the default) this file alone, and the two give identical data
# (checks/wp6_geo.R). Definitions only, except the registration of DGPs "D4",
# "D4-40", "D4-census", "D3" and "D3-claims" (dgp_register()), which builds
# each design's fixed world once (milliseconds).
#
# Structure (what the port reproduces):
#   1. geo_draw(): all randomness, as uniforms, in one documented order.
#      The world (schools, hospitals) is drawn from the design's fixed
#      world_seed, so that every replication faces the same world and the
#      same truth; the sample (students, patients) from the replication's
#      seed. D4 world, in this order: cluster centres (x, y per cluster,
#      market-major), schools (component, x, y per school, market-major), the
#      sector keys, the quality noise and the constants' noise xi (one per
#      school each). D4 sample: students' locations (component, x, y per
#      student), SES noise, tastes (K_w = 2 per student, drawn in every
#      variant so that variants share a seed's students), then the uniforms
#      that pick the choices (one per student). D3 world: town centres (x, y
#      per town), hospitals (component, x, y per hospital), the beds noise,
#      the teaching keys and xi (one per hospital each). D3 sample: patients'
#      locations (component, x, y per patient), age keys (one per patient),
#      with a panel the admission counts (one per patient), then severities
#      (one per admission), tastes (K_w = 2 per patient) and the choice
#      uniforms (one per admission). Normals are qnorm() of these uniforms.
#   2. Construction: locations (a mixture of a uniform background and
#      truncated normals around the clusters, by inversion), the world's
#      attributes and the units' (deterministic, in R in both engines), then
#      the choice sets and the long table (geo_sets() and geo_sets_radius()
#      here, geo_knn() and geo_radius() in the port): rows in choicer's
#      order (situation id, then alternative label; a patient's admissions
#      contiguous), with distance and the interactions (formed in R in both
#      engines).
#   3. Probabilities at theta0 (lib/simulate.R's engines; geo_choose() in
#      the port, MNL and MXL) and choices by inverse CDF (draw_choice()) from
#      the sample's choice uniforms.
# Port notes: distances are sqrt(dx * dx + dy * dy) in double and the
# interactions single products, so a C++ port must not contract them into
# fused multiply-adds; qnorm()/pnorm() are R's (Rmath) and stay in R.
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
# D4 designs. D4-40 (DGP "D4-40", F2's and F3's): one market of 40 schools,
# K = 15, so that choice sets still vary by location (each student's 15
# nearest of 40) while J = 40 constants keep J/sqrt(N) at 1.26, 0.63, 0.32
# on the 1k/4k/16k rungs (S1, 2026-11-02: with D4-small's 160 the distance
# coefficient's bias/SD was -0.163 J/sqrt(N), the incidental constants'
# bias C6, which the gates are not meant to test). One market of the same
# geography rather than four of 10 schools: with 10 schools and K = 15 every
# student would face the whole market, and sigma_q, identified by how the
# nearest-K sets' composition varies, would not be. D4-small (DGP "D4", 4
# markets x 40, 160 constants) stays for the scaled designs. D4-census (DGP
# "D4-census", F7's): 150 markets of 20 schools (12 public, 4 charter, 4
# private), K = 15, J = 3,000 constants, 10^6 students in the replications
# (about 6,700 per market, 1.5 x 10^7 rows).
#
# D4 theta0. The structural values below are literals chosen 2026-11-02 and
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
# D3 (hospital). One square region of side `side`, in units of 10 km
# (choicer's shifted log-normal has a median of at least 1 per unit of its
# column, too steep per km; at theta0 the median distance coefficient is
# -3.9 per 10 km), whose population clusters around `clusters` towns of
# Zipf sizes (weight w_c proportional to c^-zipf, spread cluster_sd sqrt(C
# w_c), centres in the middle 80% of the square); a share `background` of
# patients lives in the countryside, uniformly (rural_i = 1), and a share
# `hosp_background` of the hospitals. Hospital j has log beds_j = log 150 +
# z_j + 0.5 urban_j, teaches with probability plogis(-3 + 1.5 (log beds_j -
# log 150)), and has the constant delta_j = a (log beds_j - mean) +
# sigma_xi xi_j (skewed shares). Hospitals are numbered by beds,
# largest first, and labelled by that integer: the reference alternative
# (the first label; no outside option) is the largest hospital. Patient i is
# 65 or older with probability 0.4 and, in a panel design, has T_i = min(T_max,
# floor(U^(-1 / T_alpha))) admissions (P(T >= t) = t^-T_alpha for t <= T_max:
# a discrete Pareto tail, capped), else one; each admission has a severity
# sev ~ N(0, 1). The choice set, the same at each of the patient's
# admissions: the hospitals within `radius`, at least the K_min and at most
# the K_max nearest (ties by index); sets cross every boundary and connect
# the constants. Utility of hospital j at admission a of patient i, d = d_ij:
#   MNL: delta_j + b_1 d + b_2 d^2 + b_3 d age65_i + b_4 d rural_i
#        + b_5 teach_j sev_a
#   MXL: delta_j + b_2 d^2 + b_3 d age65_i + b_4 d rural_i
#        - (exp(mu_d) + exp(s_d eta_i1)) d + (mu_t + s_t eta_i2) teach_j sev_a,
# the tastes per patient, shared by the patient's admissions (the panel of
# Revelt & Train 1998). Covariates dist, dist2, dist_age65, dist_rural,
# teach_sev (MNL's X); MXL: X = dist2, dist_age65, dist_rural and W = (ndist =
# -dist, teach_sev), rc_dist = c(1, 0), rc_mean = TRUE, uncorrelated: both
# random columns vary by patient within a hospital, so their means are
# identified beside the constants (unlike D4's quality), and dist leaves X
# (its coefficient is the shift exp(mu_d)).
#
# D3 designs. D3 (DGP "D3", F6's): 200 hospitals over 200 x 200 km (10
# towns), a cross-section (T = 1), radius 3 (30 km), sets of 3 to 30
# hospitals (mean 21.5; a third at the cap). D3-claims (DGP "D3-claims",
# F7's): 2,000 hospitals over ten times the area (100 towns), a panel of
# admissions (T_alpha = 2, T_max = 20: three in four patients admitted once,
# mean 1.6), radius 4, sets of 3 to 60 hospitals (mean 39); N counts
# patients: N = 625,000 gives about 10^6 admissions and 3.9 x 10^7 rows.
#
# D3 theta0: literals chosen 2026-11-02 (D3_TRUE; the beds' spread and a
# raised from 0.6 and 1 so that the shares' skew reaches its target), the
# distance location (MXL mu_d, MNL b_1) calibrated by d3_calibrate() on D3
# the same day so that the expected share of admissions to the nearest
# hospital is 0.40 (N = 40000, seed 1, the port: mu_d = 1.074716, b_1 =
# -4.142944). Targets (D3_TARGETS): plausible ranges for US inpatient
# choice, not estimates from the papers cited: a third to a half of patients
# choose the nearest hospital, the distance elasticity of a hospital's
# choice probability about -1 to -3, shares skewed by size (J sum s_j^2 of 2
# to 6), and in the panel three quarters admitted once and a few percent five
# times or more (about one in five Medicare discharges is readmitted within
# 30 days: Jencks, Williams & Coleman 2009). At N = 200,000 (seed 1) D3's
# MXL variant gives 0.398, -1.97 and 2.21, D3-claims' 0.398, -1.99, 2.31 and
# 0.750, 0.041, 1.60 (checks/wp6_geo.R --only=calib).
#
# Alternatives never chosen in a replication (plan §4): drop_unbounded_alts()
# drops them, and the situations offering an alternative chosen wherever
# offered, until stable (ML's limit as delta_j -> -Inf or +Inf); families
# record the counts, and geo_truth() restates theta0 on the alternatives
# left.
# ============================================================================

# --- D4: designs and the truth --------------------------------------------------

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

GEO_D4_CENSUS <- list(
  kind = "school", name = "D4-census", markets = 150L, schools = 20L, K = 15L, side = 10,
  clusters = 3L, cluster_sd = 1.2, background = 0.25,
  sectors = c(public = 12L, charter = 4L, private = 4L), world_seed = 50261102L,
  sector_delta = list(default = c(public = 3.0054, charter = 1.7382, private = 1.7566),
                      nl = c(public = 2.8849, charter = 1.2893, private = 1.2819),
                      mnl = c(public = 2.3338, charter = 1.1699, private = 1.2352)))

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

# --- D3: designs and the truth --------------------------------------------------

# A design: the region (side, in 10 km), the towns (count, Zipf exponent,
# base spread), the rural shares of patients and hospitals, the hospitals,
# the choice sets (radius, K_min, K_max), the admissions (panel or not, the
# Pareto tail and its cap), the world's seed.
GEO_D3 <- list(
  kind = "hospital", name = "D3", side = 20, clusters = 10L, zipf = 1, cluster_sd = 1.3,
  background = 0.25, hosp_background = 0.10, hospitals = 200L, radius = 3, K_min = 3L,
  K_max = 30L, panel = FALSE, T_alpha = 2, T_max = 20L, world_seed = 30261102L)

GEO_D3_CLAIMS <- list(
  kind = "hospital", name = "D3-claims", side = 63.25, clusters = 100L, zipf = 1,
  cluster_sd = 1.3, background = 0.25, hosp_background = 0.10, hospitals = 2000L, radius = 4,
  K_min = 3L, K_max = 60L, panel = TRUE, T_alpha = 2, T_max = 20L, world_seed = 31261102L)

D3_X <- list(default = c("dist2", "dist_age65", "dist_rural"),
             mnl = c("dist", "dist2", "dist_age65", "dist_rural", "teach_sev"))
D3_W <- c("ndist", "teach_sev")
D3_VARIANTS <- c(default = "mxl", mnl = "mnl")

# The structural truth, shared by the designs (distances in 10 km): the
# patients' and hospitals' attributes, the constants' a and sigma_xi, the
# MNL's beta and the MXL's (beta, mu = (mu_d, mu_t), s = sd of g = L eta).
D3_TRUE <- list(
  age65 = 0.4, beds_mu = log(150), beds_sd = 1.0, beds_urban = 0.5, teach = c(-3, 1.5),
  a = 1.5, sigma_xi = 0.5,
  beta = list(default = c(dist2 = 0.15, dist_age65 = -0.3, dist_rural = 0.8),
              mnl = c(dist = -4.1429, dist2 = 0.15, dist_age65 = -0.3, dist_rural = 0.8,
                      teach_sev = 0.5)),
  mu = c(ndist = 1.0747, teach_sev = 0.5),
  s = c(ndist = 0.8, teach_sev = 0.6))

# The engines' specification of a D3 variant over the hospitals `alts`.
d3_spec <- function(variant, hospitals, alts = hospitals$alt) {
  model <- D3_VARIANTS[[variant]]
  spec <- list(beta = D3_X[[variant]], alts = alts, outside = FALSE)
  if (model == "mxl") {
    spec <- c(spec, list(random = D3_W, rc_dist = c(1L, 0L), rc_mean = TRUE,
                         correlated = FALSE))
  }
  spec
}

# A variant's named blocks: beta, MXL mu and L, and the constants by label,
# delta_j = a (log beds_j - mean) + sigma_xi xi_j (pack_theta() takes them
# relative to the reference, the largest hospital).
d3_blocks <- function(variant, hospitals, par) {
  model <- D3_VARIANTS[[variant]]
  delta <- par$a * (hospitals$log_beds - mean(hospitals$log_beds)) + par$sigma_xi * hospitals$xi
  b <- list(beta = par$beta[[variant]], asc = stats::setNames(delta, hospitals$alt))
  if (model == "mxl") {
    b$mu <- par$mu[D3_W]
    b$L <- diag(unname(par$s[D3_W]), 2L)
  }
  b
}

d3_theta0 <- function(variant, hospitals, par) {
  pack_theta(D3_VARIANTS[[variant]], d3_spec(variant, hospitals), d3_blocks(variant, hospitals, par))
}

# --- 1. Draws ------------------------------------------------------------------------

# Every uniform a dataset uses, in the order of the header: `world` from
# design$world_seed, `sample` from `seed` (N students or patients, K_w tastes
# each). D3 also returns T, the admissions per patient (all 1 without a
# panel), which sets how many severities and choice uniforms follow.
geo_draw <- function(design, N, seed, K_w = 2L) {
  N <- as.integer(N)
  if (design$kind == "hospital") return(.d3_draw(design, N, seed, K_w))
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

# Admissions per patient from their count uniforms: min(T_max,
# floor(u^(-1/T_alpha))) in a panel design, else 1.
.d3_admissions <- function(u, design) {
  if (!isTRUE(design$panel)) return(rep(1L, length(u)))
  as.integer(pmin(design$T_max, floor(u^(-1 / design$T_alpha))))
}

.d3_draw <- function(design, N, seed, K_w) {
  J <- design$hospitals
  world <- vs_with_seed(design$world_seed, {
    centers <- matrix(stats::runif(2L * design$clusters), 2L)
    hosp <- matrix(stats::runif(3L * J), 3L)
    beds <- stats::runif(J)
    teach <- stats::runif(J)
    xi <- stats::runif(J)
    list(centers = centers, hosp = hosp, beds = beds, teach = teach, xi = xi)
  })
  sample <- vs_with_seed(seed, {
    loc <- matrix(stats::runif(3L * N), 3L)
    age <- stats::runif(N)
    adm <- if (isTRUE(design$panel)) .d3_admissions(stats::runif(N), design) else rep(1L, N)
    sev <- stats::runif(sum(adm))
    eta <- matrix(stats::runif(K_w * N), K_w)
    choice <- stats::runif(sum(adm))
    list(loc = loc, age = age, T = adm, sev = sev, eta = eta, choice = choice)
  })
  list(world = world, sample = sample, N = N, A = sum(sample$T), T = sample$T,
       K_w = as.integer(K_w))
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
  if (design$kind != "school") stop("geo_schools(): a D4 design (kind \"school\")", call. = FALSE)
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
# ndist, quality) in choicer's row order, `code` (each row's school index,
# choicer's alternative code over every school) and `eta`, the tastes'
# normals (K_w x N). The sets come from geo_sets() (engine "r") or the
# port's geo_knn() (engine "cpp"); the rest is shared.
geo_build <- function(design, draws, par = D4_TRUE, engine = "r") {
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
  long <- if (.geo_engine(engine) == "cpp") {
    data.table::setDT(.geo_gen()$geo_knn(students$x, students$y, students$market, schools$x,
                                         schools$y, schools$market, design$K))
  } else {
    geo_sets(students, schools, design$K)
  }
  ses_r <- students$ses[long$id]
  q_r <- schools$quality[long$j]
  long[, `:=`(alt = schools$alt[j], market = schools$market[j], sector = schools$sector[j],
              dist_ses = dist * ses_r, q_ses = q_r * ses_r, ndist = -dist, quality = q_r)]
  code <- long$j
  long[, j := NULL]
  data.table::setcolorder(long, c("id", "alt", "market", "sector", "dist", "dist_ses", "q_ses",
                                  "ndist", "quality"))
  list(schools = schools, students = students, long = long, code = code,
       eta = matrix(stats::qnorm(s$eta), nrow(s$eta)))
}

# D3's town weights (Zipf) and spreads.
.d3_towns <- function(design) {
  w <- seq_len(design$clusters)^(-design$zipf)
  w <- w / sum(w)
  list(w = w, sd = design$cluster_sd * sqrt(design$clusters * w))
}

# D3 locations from (component, x, y) uniforms (3 x n): rural (uniform on
# the square) with probability p0, else a truncated normal around town k,
# drawn with the town weights. Returns x, y and rural (0/1).
.d3_locate <- function(u, ctr, design, p0) {
  tw <- .d3_towns(design)
  side <- design$side
  bg <- u[1L, ] < p0
  k <- pmin(design$clusters, findInterval((u[1L, ] - p0) / (1 - p0), cumsum(tw$w)) + 1L)
  k[bg] <- 1L
  x <- ifelse(bg, side * u[2L, ], .geo_tnorm(u[2L, ], ctr[1L, k], tw$sd[k], 0, side))
  y <- ifelse(bg, side * u[3L, ], .geo_tnorm(u[3L, ], ctr[2L, k], tw$sd[k], 0, side))
  list(x = x, y = y, rural = as.integer(bg))
}

# Town centres (2 x C), in the middle 80% of the square.
.d3_centers <- function(design, world) design$side * (0.1 + 0.8 * world$centers)

# The world's hospitals from its draws, numbered by beds (largest first):
# j, alt (= j), x, y, urban, log_beds, teaching (0/1) and xi.
geo_hospitals <- function(design, world, par = D3_TRUE) {
  if (design$kind != "hospital") stop("geo_hospitals(): a D3 design (kind \"hospital\")",
                                      call. = FALSE)
  loc <- .d3_locate(world$hosp, .d3_centers(design, world), design, design$hosp_background)
  urban <- 1L - loc$rural
  log_beds <- par$beds_mu + par$beds_sd * stats::qnorm(world$beds) + par$beds_urban * urban
  zb <- (log_beds - par$beds_mu) / par$beds_sd
  teaching <- as.integer(world$teach < stats::plogis(par$teach[1L] + par$teach[2L] * zb))
  o <- order(log_beds, decreasing = TRUE)
  J <- length(o)
  data.table::data.table(j = seq_len(J), alt = seq_len(J), x = loc$x[o], y = loc$y[o],
                         urban = urban[o], log_beds = log_beds[o], teaching = teaching[o],
                         xi = stats::qnorm(world$xi)[o])
}

# Each patient's choice set: the hospitals within design$radius, at least
# the K_min and at most the K_max nearest (ties by index), as a data.table
# (patient, j, dist) in (patient, j) order; brute force over every hospital,
# `chunk` patients at a time.
geo_sets_radius <- function(patients, hospitals, design, chunk = 2000L) {
  J <- nrow(hospitals)
  n <- nrow(patients)
  k_min <- min(design$K_min, J)
  out <- lapply(split(seq_len(n), (seq_len(n) - 1L) %/% chunk), function(ix) {
    m <- length(ix)
    dx <- rep(patients$x[ix], each = J) - rep(hospitals$x, m)
    dy <- rep(patients$y[ix], each = J) - rep(hospitals$y, m)
    long <- data.table::data.table(patient = rep(patients$patient[ix], each = J),
                                   j = rep(hospitals$j, m), dist = sqrt(dx * dx + dy * dy))
    data.table::setorderv(long, c("patient", "dist", "j"))
    rk <- data.table::rowidv(long, "patient")
    n_r <- .group_sum(as.numeric(long$dist <= design$radius), data.table::rleid(long$patient), m)
    take <- pmin(design$K_max, pmax(k_min, n_r))
    long <- long[rk <= rep(take, each = J)]
    data.table::setorderv(long, c("patient", "j"))
    long
  })
  data.table::rbindlist(out)
}

# The patients from the draws: patient (1..N), x, y, rural, age65 (0/1), T.
d3_patients <- function(design, draws, par = D3_TRUE) {
  s <- draws$sample
  loc <- .d3_locate(s$loc, .d3_centers(design, draws$world), design, design$background)
  data.table::data.table(patient = seq_len(draws$N), x = loc$x, y = loc$y, rural = loc$rural,
                         age65 = as.integer(s$age < par$age65), T = s$T)
}

# The world and the sample from the draws: `hospitals`, `patients`, `sev`
# (one per admission), `long` (id = admission, patient, alt, dist, dist2,
# dist_age65, dist_rural, teach_sev, ndist) in choicer's row order (a
# patient's admissions contiguous, numbered 1..A patient by patient), and
# `eta`, the tastes' normals (K_w x N, one column per patient). The sets
# come from geo_sets_radius() (engine "r") or the port's geo_radius()
# (engine "cpp"); the covariates are formed here in both.
d3_build <- function(design, draws, par = D3_TRUE, engine = "r") {
  s <- draws$sample
  hospitals <- geo_hospitals(design, draws$world, par)
  patients <- d3_patients(design, draws, par)
  if (.geo_engine(engine) == "cpp") {
    long <- data.table::setDT(.geo_gen()$geo_radius(
      patients$x, patients$y, patients$T, hospitals$x, hospitals$y, design$side, design$radius,
      design$K_min, design$K_max))
  } else {
    sets <- geo_sets_radius(patients, hospitals, design)
    # Each patient's rows repeated over its admissions, admissions numbered
    # patient by patient.
    len <- tabulate(sets$patient, nrow(patients))
    first <- cumsum(c(1L, len[-length(len)]))
    adm_len <- rep(len, patients$T)
    rows <- sequence(adm_len, from = rep(first, patients$T))
    long <- data.table::data.table(id = rep(seq_along(adm_len), adm_len),
                                   patient = sets$patient[rows], j = sets$j[rows],
                                   dist = sets$dist[rows])
  }
  sev <- stats::qnorm(s$sev)
  data.table::setnames(long, "j", "alt")
  long[, `:=`(dist2 = dist * dist, dist_age65 = dist * patients$age65[patient],
              dist_rural = dist * patients$rural[patient],
              teach_sev = hospitals$teaching[alt] * sev[id], ndist = -dist)]
  list(hospitals = hospitals, patients = patients, sev = sev, long = long,
       eta = matrix(stats::qnorm(s$eta), nrow(s$eta)))
}

# --- The port (lib/geo_gen.cpp) --------------------------------------------------------

# The compiled port's functions: geo_gen_compile() fills this environment.
.geo_cpp <- new.env(parent = baseenv())

# Compiles lib/geo_gen.cpp with Rcpp::sourceCpp() into `cache` (default:
# the environment variable VS_GEO_CACHE, else _validation/output/
# geo_gen_cache, through the output guard) with PKG_CXXFLAGS =
# -ffp-contract=off, loads it and checks it (.geo_selftest()); returns the
# cache directory, invisibly. sourceCpp() reuses a build cached there for
# the same source, so processes sharing a cache compile once: mc.R --warm
# calls this before shards start, and every shard then loads the cached
# build. The cache is locked while a process compiles or loads
# (.geo_lock(): a directory created atomically, waited for up to `wait`
# seconds), so shards started without the warm pass do not build one cache
# at once.
geo_gen_compile <- function(cache = NULL, rebuild = FALSE, verbose = FALSE, wait = 600) {
  if (is.null(cache)) cache <- Sys.getenv("VS_GEO_CACHE", unset = "")
  if (!nzchar(cache)) cache <- file.path(vs_root(), "output", "geo_gen_cache")
  cache <- vs_out_dir(out = cache)
  lock <- .geo_lock(cache, wait)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  old <- Sys.getenv("PKG_CXXFLAGS", unset = NA)
  Sys.setenv(PKG_CXXFLAGS = trimws(paste(if (!is.na(old)) old, "-ffp-contract=off")))
  on.exit(if (is.na(old)) Sys.unsetenv("PKG_CXXFLAGS") else Sys.setenv(PKG_CXXFLAGS = old),
          add = TRUE)
  Rcpp::sourceCpp(file.path(vs_root(), "lib", "geo_gen.cpp"), cacheDir = cache, env = .geo_cpp,
                  rebuild = rebuild, verbose = verbose)
  .geo_selftest(cache)
  invisible(cache)
}

# Takes the cache's lock, `geo_gen.lock`, a directory (dir.create() is
# atomic), waiting for another process's for up to `wait` seconds; returns
# its path. A lock left by a killed process is reported, not broken.
.geo_lock <- function(cache, wait) {
  lock <- file.path(cache, "geo_gen.lock")
  t0 <- proc.time()[["elapsed"]]
  while (!dir.create(lock, showWarnings = FALSE)) {
    if (proc.time()[["elapsed"]] - t0 > wait) {
      stop("geo_gen_compile(): the build cache ", cache, " has been locked by another process ",
           "for ", wait, " s (", lock, "); if no process is compiling there, remove that ",
           "directory", call. = FALSE)
    }
    Sys.sleep(0.5)
  }
  lock
}

# The loaded port's distances against R's, bit for bit: geo_knn() on fixed
# points (400 units x 50 alternatives in one market, every alternative in
# each set) against sqrt(dx * dx + dy * dy) evaluated by R. A build that
# contracts multiply-adds differs in a few percent of them; then the loaded
# functions are removed and the build refused, naming the flags.
.geo_selftest <- function(cache) {
  pts <- vs_with_seed(20261102L, list(u = matrix(stats::runif(800L, 0, 20), 2L),
                                      a = matrix(stats::runif(100L, 0, 20), 2L)))
  n <- ncol(pts$u)
  J <- ncol(pts$a)
  res <- .geo_cpp$geo_knn(pts$u[1L, ], pts$u[2L, ], rep(1L, n), pts$a[1L, ], pts$a[2L, ],
                          rep(1L, J), J)
  dx <- rep(pts$u[1L, ], each = J) - rep(pts$a[1L, ], n)
  dy <- rep(pts$u[2L, ], each = J) - rep(pts$a[2L, ], n)
  ref <- sqrt(dx * dx + dy * dy)
  bad <- if (identical(res$j, rep(seq_len(J), n))) sum(res$dist != ref) else n * J
  if (bad > 0L) {
    rm(list = ls(.geo_cpp, all.names = TRUE), envir = .geo_cpp)
    stop("geo_gen_compile(): the build in ", cache, " computes ", bad, " of ", n * J,
         " distances differently from R (floating-point contraction into fused multiply-adds). ",
         "It is compiled with PKG_CXXFLAGS=-ffp-contract=off; a CXXFLAGS with -ffp-contract=fast ",
         "or -ffast-math in ~/.R/Makevars (or the file R_MAKEVARS_USER names) comes after it ",
         "and wins. Remove that flag, then rebuild: geo_gen_compile(rebuild = TRUE).",
         call. = FALSE)
  }
  invisible(TRUE)
}

# The port's functions, compiled (or loaded from the cache) on first use.
.geo_gen <- function() {
  if (!exists("geo_choose", envir = .geo_cpp, inherits = FALSE)) geo_gen_compile()
  .geo_cpp
}

.geo_engine <- function(engine) {
  if (!identical(engine, "r") && !identical(engine, "cpp")) {
    stop("engine must be \"r\" or \"cpp\"", call. = FALSE)
  }
  engine
}

# --- 3. Probabilities and choices --------------------------------------------------------

# The probabilities choices are drawn from, at theta: MNL/NL the model's;
# MXL the logit given each unit's tastes (rc_coefs() at the drawn eta; the
# unit is the situation, or the decision maker with d$unit).
.geo_probs <- function(model, d, par, eta) {
  switch(model,
         mnl = mnl_probs(d, par),
         nl = nl_probs(d, par),
         mxl = mxl_cond_probs(d, par, t(rc_coefs(par$mu, par$L, eta, par$rc_dist, par$rc_mean))))
}

# A long table's probabilities at theta (par, unpack_theta()'s) and its
# choices from the uniforms u (one per situation): list(p, p0, choice,
# expected), `expected` the sum of p per alternative (named by spec$alts).
# Engine "r": lib/simulate.R's engines and draw_choice(); "cpp" (MNL and
# MXL): geo_choose(), with `code` each row's alternative code (spec$alts
# must be in code order) and `unit` each situation's decision maker (MXL;
# NULL: the situation), returning p per row only with keep_p (8 bytes a
# row; the engines' p is there anyway). The NL always takes the R engines.
.geo_outcomes <- function(model, long, spec, par, eta, u, code, unit = NULL, person = NULL,
                          engine = "r", keep_p = FALSE) {
  expected <- stats::setNames(numeric(length(spec$alts)), spec$alts)
  if (.geo_engine(engine) == "cpp" && model != "nl") {
    if (!identical(.alt_levels(spec), spec$alts)) stop("spec$alts must be in code order")
    mxl <- model == "mxl"
    tastes <- if (mxl) t(rc_coefs(par$mu, par$L, eta, par$rc_dist, par$rc_mean))
              else matrix(0, 0L, 0L)
    res <- .geo_gen()$geo_choose(
      long$id, code, lapply(spec$beta, function(k) long[[k]]), unname(par$beta),
      unname(par$delta), if (mxl) lapply(spec$random, function(k) long[[k]]) else list(),
      tastes, if (mxl) as.integer(unit %||% seq_along(u)) else integer(0), spec$outside, u,
      keep_p)
    expected[] <- res$expected
    return(list(p = if (keep_p) res$p, p0 = res$p0, choice = res$choice, expected = expected))
  }
  d <- design(long, spec, person = person)
  pr <- .geo_probs(model, d, par, eta)
  choice <- draw_choice(pr$p, pr$p0, d$sit, u = u)$choice
  e <- rowsum(pr$p, long$alt)
  expected[rownames(e)] <- e[, 1L]
  list(p = pr$p, p0 = pr$p0, choice = choice, expected = expected)
}

# A D4 dataset: list(data, truth, meta). `data` the long table with
# `choice` (0/1); `truth` theta0 over every school; `meta`: variant, model,
# spec and blocks (for geo_truth()), N, the design's summary, the expected
# shares (mean probability per school label and of the outside option) and
# counts, the schools and students tables, the tastes' normals `eta` (K_w x
# N; MXL), the engine, and with keep_probs the probabilities the choices
# were drawn from (`probs`: p per row, p0 per student). engine "cpp" builds
# the sets, probabilities and choices with the port (identical data).
d4_simulate <- function(N, seed, variant = "default", design = GEO_D4_SMALL, par = d4_par(design),
                        keep_probs = FALSE, engine = "r") {
  if (!variant %in% names(D4_VARIANTS)) {
    stop("D4 has no variant '", variant, "'; variants: ", paste(names(D4_VARIANTS), collapse = ", "),
         call. = FALSE)
  }
  model <- D4_VARIANTS[[variant]]
  draws <- geo_draw(design, N, seed)
  b <- geo_build(design, draws, par, engine)
  spec <- d4_spec(variant, b$schools)
  blocks <- d4_blocks(variant, b$schools, par)
  truth <- pack_theta(model, spec, blocks)
  out <- .geo_outcomes(model, b$long, spec, unpack_theta(model, spec, truth), b$eta,
                       draws$sample$choice, b$code, engine = engine, keep_p = keep_probs)
  data <- b$long
  data[, choice := out$choice]
  meta <- list(variant = variant, model = model, spec = spec, blocks = blocks, N = draws$N,
               design = design[c("name", "markets", "schools", "K", "side")],
               J = nrow(b$schools), expected = out$expected, outside_expected = sum(out$p0),
               shares = out$expected / draws$N, outside_share = sum(out$p0) / draws$N,
               schools = b$schools, students = b$students,
               eta = if (model == "mxl") b$eta, engine = engine,
               probs = if (keep_probs) out[c("p", "p0")])
  list(data = data, truth = truth, meta = meta)
}

# A D3 dataset: list(data, truth, meta). `data` the long table (id =
# admission, patient, alt, dist, dist2, dist_age65, dist_rural, teach_sev,
# ndist, choice); `truth` theta0 over every hospital; `meta`: variant,
# model, spec, blocks, N (patients), A (admissions), the design's summary,
# J, the expected counts and shares per hospital, the hospitals and patients
# tables, the severities, the tastes' normals `eta` (K_w x N; MXL), the
# engine, and with keep_probs `probs` (p per row). The MXL's unit is the
# patient: its tastes are shared by the patient's admissions.
d3_simulate <- function(N, seed, variant = "default", design = GEO_D3, par = D3_TRUE,
                        keep_probs = FALSE, engine = "r") {
  if (!variant %in% names(D3_VARIANTS)) {
    stop("D3 has no variant '", variant, "'; variants: ", paste(names(D3_VARIANTS), collapse = ", "),
         call. = FALSE)
  }
  model <- D3_VARIANTS[[variant]]
  draws <- geo_draw(design, N, seed)
  b <- d3_build(design, draws, par, engine)
  spec <- d3_spec(variant, b$hospitals)
  blocks <- d3_blocks(variant, b$hospitals, par)
  truth <- pack_theta(model, spec, blocks)
  unit <- rep(b$patients$patient, b$patients$T)
  out <- .geo_outcomes(model, b$long, spec, unpack_theta(model, spec, truth), b$eta,
                       draws$sample$choice, b$long$alt, unit = unit,
                       person = if (isTRUE(design$panel)) "patient", engine = engine,
                       keep_p = keep_probs)
  data <- b$long
  data[, choice := out$choice]
  meta <- list(variant = variant, model = model, spec = spec, blocks = blocks, N = draws$N,
               A = draws$A, design = design[c("name", "hospitals", "side", "radius", "K_min",
                                              "K_max", "panel")],
               J = nrow(b$hospitals), expected = out$expected, shares = out$expected / draws$A,
               hospitals = b$hospitals, patients = b$patients, sev = b$sev,
               eta = if (model == "mxl") b$eta, engine = engine,
               probs = if (keep_probs) out[c("p", "p0")])
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
# (situations), counts = c(alts, situations, rows, rounds)). A data.table
# is not copied when nothing is dropped (rows are dropped into a new table:
# the input is never modified).
drop_unbounded_alts <- function(data, id_col = "id", alt_col = "alt", choice_col = "choice") {
  dt <- if (data.table::is.data.table(data)) data else data.table::as.data.table(data)
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

# The connected components of the alternatives linked by shared choice
# sets (rows: situation id and alternative label): the number of
# components; one means every constant is identified against the
# reference. Label propagation: each alternative takes the smallest label
# among the alternatives it shares a set with, until stable.
.geo_components <- function(id, alt) {
  sit <- data.table::rleid(id)
  a <- match(alt, sort(unique(alt)))
  J <- max(a)
  comp <- seq_len(J)
  repeat {
    m <- -.group_max(-comp[a], sit, sit[length(sit)])
    new <- pmin(comp, -.group_max(-m[sit], a, J))
    if (identical(new, comp)) break
    comp <- new
  }
  length(unique(comp))
}

# A D3 dataset's probabilities (p per row), from meta$probs or recomputed
# with the dataset's engine.
.d3_probs <- function(sim) {
  if (!is.null(sim$meta$probs)) return(sim$meta$probs$p)
  m <- sim$meta
  par <- unpack_theta(m$model, m$spec, sim$truth)
  .geo_outcomes(m$model, sim$data, m$spec, par, m$eta, rep(0.5, m$A), sim$data$alt,
                unit = rep(m$patients$patient, m$patients$T),
                person = if (isTRUE(m$design$panel)) "patient", engine = m$engine,
                keep_p = TRUE)$p
}

# D3's moments on a simulated dataset (before any drop): the share of
# admissions to the nearest hospital of their set (realized), the distance
# elasticity (the probability-weighted mean over admissions of each
# hospital's own elasticity d log P / d log d = (dV/dd) d (1 - P), at the
# patients' realized tastes), the size skew J sum s_j^2 of the expected
# shares, and with a panel the admissions per patient (share admitted once,
# share five times or more, mean); untargeted: the expected nearest share,
# the mean distance travelled (km), the teaching share of admissions, the
# set sizes (mean, 10th and 90th percentiles; shares at K_min and K_max),
# the expected count per hospital (min, 10th percentile, median), J, and
# the connected components of the sets (one: every constant identified).
d3_moments <- function(sim) {
  dt <- sim$data
  m <- sim$meta
  A <- m$A
  p <- .d3_probs(sim)
  sit <- data.table::rleid(dt$id)
  dmin <- -.group_max(-dt$dist, sit, A)
  near <- dt$dist == dmin[sit]
  ch <- dt$choice == 1L
  b <- unpack_theta(m$model, m$spec, sim$truth)$beta
  slope <- b[["dist2"]] * 2 * dt$dist2 + b[["dist_age65"]] * dt$dist_age65 +
    b[["dist_rural"]] * dt$dist_rural
  if (m$model == "mnl") {
    slope <- slope + b[["dist"]] * dt$dist
  } else {
    par <- unpack_theta(m$model, m$spec, sim$truth)
    tastes <- t(rc_coefs(par$mu, par$L, m$eta, par$rc_dist, par$rc_mean))
    slope <- slope + tastes[dt$patient, 1L] * dt$ndist
  }
  size <- tabulate(sit, A)
  sh <- m$expected / sum(m$expected)
  e <- m$expected
  T <- m$patients$T
  c(nearest_share = sum(ch & near) / A,
    dist_elasticity = sum(p * slope * (1 - p)) / A,
    share_hhi_rel = length(sh) * sum(sh^2),
    share_T1 = if (isTRUE(m$design$panel)) mean(T == 1L) else NA_real_,
    share_T5plus = if (isTRUE(m$design$panel)) mean(T >= 5L) else NA_real_,
    mean_T = mean(T),
    nearest_expected = sum(p[near]) / A,
    mean_dist_km = 10 * sum(dt$dist[ch]) / A,
    teaching_share = sum(m$hospitals$teaching[dt$alt[ch]]) / A,
    set_mean = mean(size), set_p10 = unname(stats::quantile(size, 0.1)),
    set_p90 = unname(stats::quantile(size, 0.9)),
    set_at_K_min = mean(size == m$design$K_min), set_at_K_max = mean(size == m$design$K_max),
    count_min = min(e), count_p10 = unname(stats::quantile(e, 0.1)), count_med = stats::median(e),
    J = length(e), components = .geo_components(dt$id, dt$alt))
}

D3_TARGETS <- data.frame(
  moment = c("nearest_share", "dist_elasticity", "share_hhi_rel"),
  lo = c(0.30, -3.0, 2), hi = c(0.55, -1.0, 6), stringsAsFactors = FALSE)

D3_PANEL_TARGETS <- rbind(D3_TARGETS, data.frame(
  moment = c("share_T1", "share_T5plus", "mean_T"),
  lo = c(0.65, 0.01, 1.3), hi = c(0.85, 0.06, 1.8), stringsAsFactors = FALSE))

# Calibrates the distance location of a variant (MXL mu_d, MNL b_1) by
# bisection on [lo, hi] so that the expected share of admissions to the
# nearest hospital (mean probability, N patients of one seed, common draws)
# is `target`; more distance aversion raises it. Not run when the file is
# sourced; its result is frozen in D3_TRUE.
d3_calibrate <- function(variant, N = 40000L, seed = 1L, target = 0.40, design = GEO_D3,
                         par = D3_TRUE, lo = -3, hi = 3, tol = 1e-6, engine = "r") {
  model <- D3_VARIANTS[[variant]]
  draws <- geo_draw(design, N, seed)
  b <- d3_build(design, draws, par, engine)
  spec <- d3_spec(variant, b$hospitals)
  unit <- rep(b$patients$patient, b$patients$T)
  sit <- data.table::rleid(b$long$id)
  near <- b$long$dist == -.group_max(-b$long$dist, sit, draws$A)[sit]
  share_at <- function(x) {
    if (model == "mxl") par$mu[["ndist"]] <- x else par$beta$mnl[["dist"]] <- -exp(x)
    th <- pack_theta(model, spec, d3_blocks(variant, b$hospitals, par))
    p <- .geo_outcomes(model, b$long, spec, unpack_theta(model, spec, th), b$eta,
                       draws$sample$choice, b$long$alt, unit = unit,
                       person = if (isTRUE(design$panel)) "patient", engine = engine,
                       keep_p = TRUE)$p
    sum(p[near]) / draws$A
  }
  # The MNL's b_1 = -exp(x), so that both locations are searched on a log scale.
  f_lo <- share_at(lo) - target
  for (it in seq_len(200L)) {
    mid <- (lo + hi) / 2
    f_mid <- share_at(mid) - target
    if (sign(f_mid) == sign(f_lo)) { lo <- mid; f_lo <- f_mid } else hi <- mid
    if (hi - lo < tol) break
  }
  x <- (lo + hi) / 2
  list(value = if (model == "mxl") c(mu_d = x) else c(b_dist = -exp(x)),
       nearest_expected = share_at(x), iterations = it)
}

# --- The families' counterfactuals (F6, F7) --------------------------------------------------

# The alternatives that take h's choosers when h closes, at theta0 under the
# DGP's own probabilities (lib/simulate.R's engines; the MXL at the units'
# realized tastes) on a simulated dataset before any drop: the diversion
# ratios D_hk = (S_k without h - S_k) / S_h over the situations offering h,
# named by label ("outside" for the outside option) in decreasing order.
# The families name their derived estimands by its first labels.
geo_recipients <- function(sim, h) {
  m <- sim$meta
  par <- unpack_theta(m$model, m$spec, sim$truth)
  ids <- unique(sim$data$id[sim$data$alt == h])
  aff <- sim$data[sim$data$id %in% ids]
  panel <- isTRUE(m$design$panel)
  shares <- function(dt) {
    units <- if (is.null(dt$patient)) unique(dt$id) else unique(dt$patient)
    eta <- if (m$model == "mxl") m$eta[, units, drop = FALSE]
    out <- .geo_outcomes(m$model, dt, m$spec, par, eta, rep(0.5, data.table::uniqueN(dt$id)),
                         person = if (panel) "patient")
    c(out$expected, outside = sum(out$p0))
  }
  s0 <- shares(aff)
  s1 <- shares(aff[aff$alt != h])
  d <- (s1 - s0) / s0[[as.character(h)]]
  d <- d[names(d) != as.character(h) & (names(d) != "outside" | m$spec$outside)]
  sort(d, decreasing = TRUE)
}

# choicer's post-estimation for the derived estimands (never used to
# simulate). `fit` is a choicer fit, or its at_truth() copy for the truth;
# `data` a long table with the fit's columns.

# Predicted shares of `newdata` (predict(type = "shares"): the mean over its
# situations), named by label, "outside" for the outside option.
geo_shares <- function(fit, newdata) {
  s <- stats::predict(fit, newdata = newdata, type = "shares")[, 1L]
  am <- fit$alt_mapping
  stats::setNames(s, ifelse(am$alt_int == 0L, "outside",
                            as.character(am[[fit$data_spec$alt_col]])))
}

# Closing alternative h, over the situations of `data` that offer it: the
# diversion ratios D_hk = (S_k without h - S_k) / S_h to the labels `to`
# ("outside" allowed), and the change in consumer surplus per situation of
# the sample (n_total situations), in utils (the designs have no money
# price): the summed change in logsum() over the situations offering h,
# over n_total. Each situation keeps its draws in both scenarios (one block
# per situation, by position).
geo_closure <- function(fit, data, h, to, n_total) {
  ids <- unique(data$id[data$alt == h])
  aff <- data[data$id %in% ids]
  cut <- aff[aff$alt != h]
  s0 <- geo_shares(fit, aff)
  s1 <- geo_shares(fit, cut)
  list(div = (s1[to] - s0[to]) / s0[[as.character(h)]],
       dcs = sum(logsum(fit, newdata = cut) - logsum(fit, newdata = aff)) / n_total)
}

# Raising school j's quality by `shift` (D4: quality and q_ses = quality x
# SES on its rows; its constant held, as a fit with constants must), over
# the students offered j: the change in j's share of all n_total students,
# and the change in consumer surplus per student (utils). `ses` is SES by
# student id.
geo_quality <- function(fit, data, j, shift, ses, n_total) {
  ids <- unique(data$id[data$alt == j])
  aff <- data[data$id %in% ids]
  up <- data.table::copy(aff)
  r <- which(up$alt == j)
  q_new <- up$quality[r] + shift
  data.table::set(up, r, "quality", q_new)
  data.table::set(up, r, "q_ses", q_new * ses[up$id[r]])
  s0 <- geo_shares(fit, aff)
  s1 <- geo_shares(fit, up)
  list(dshare = (s1[[j]] - s0[[j]]) * length(ids) / n_total,
       dcs = sum(logsum(fit, newdata = up) - logsum(fit, newdata = aff)) / n_total)
}

# Each constant's expected count at theta0 in a rep (the sum of the DGP's
# probabilities, by ASC name, the fit's constants only): F6 and F7 put it in
# diag$asc_expected.
geo_asc_expected <- function(fit, sim) {
  asc <- names(stats::coef(fit))[fit$param_map$asc]
  e <- sim$meta$expected
  stats::setNames(unname(e[sub("^ASC_", "", asc)]), asc)
}

# Expected-count terciles of the constants over a cell's records (each
# constant's mean diag$asc_expected): a named vector ASC name -> "low",
# "mid", "high", the `groups` of pooled_tails() (plan §2, F6 and F7).
geo_count_groups <- function(records) {
  e <- lapply(Filter(function(r) isTRUE(r$ok), records), function(r) r$diag$asc_expected)
  e <- Filter(Negate(is.null), e)
  if (!length(e)) return(NULL)
  nm <- unique(unlist(lapply(e, names)))
  m <- rowMeans(do.call(cbind, lapply(e, function(x) x[nm])), na.rm = TRUE)
  g <- cut(m, unique(stats::quantile(m, c(0, 1 / 3, 2 / 3, 1))), include.lowest = TRUE,
           labels = FALSE)
  stats::setNames(factor(c("low", "mid", "high")[g], levels = c("low", "mid", "high")), nm)
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

# Registers a D3 design as a DGP (N counts patients).
d3_register <- function(id, design, name) {
  world <- geo_hospitals(design, geo_draw(design, 0L, 1L)$world)
  dgp_register(list(
    id = id, name = name, design = design,
    theta0 = lapply(stats::setNames(nm = names(D3_VARIANTS)), d3_theta0, hospitals = world,
                    par = D3_TRUE),
    simulate = function(N, seed, variant = "default", ...) {
      d3_simulate(N, seed, variant, design = design, ...)
    },
    moments = d3_moments, targets = if (isTRUE(design$panel)) D3_PANEL_TARGETS else D3_TARGETS,
    variants = names(D3_VARIANTS)))
}

d4_register("D4", GEO_D4_SMALL,
            "school choice, D4-small (4 markets x 40 schools, K = 15, outside option)")
d4_register("D4-40", GEO_D4_40,
            "school choice, D4-40 (1 market of 40 schools, K = 15, outside option)")
d4_register("D4-census", GEO_D4_CENSUS,
            "school choice, census (150 markets x 20 schools, K = 15, outside option)")
d3_register("D3", GEO_D3,
            "hospital choice, D3 (200 hospitals, sets within 30 km, cross-section)")
d3_register("D3-claims", GEO_D3_CLAIMS,
            "hospital choice, claims (2,000 hospitals, sets within 40 km, admissions panel)")
