#!/usr/bin/env Rscript
# ============================================================================
# Golden identity harness: every compiled entry point of the MNL, MXL and NL
# kernels (and the preparations, the HB/MNP samplers and fits through the
# run_*() wrappers) over grids of flags, draw modes, weights, panels, edge
# inputs and error cases. `baseline` saves the outputs of one build; `check`
# compares another build's outputs with them, leaf by leaf, with identical()
# (bit for bit), or, for the tolerance group (the MXL estimation kernels'
# summation-order outputs; --tol-slots), within 1e-10 of max(1, max|base|).
#
# Usage (the build under test is the FIRST entry of R_LIBS; run at one thread):
#   R_LIBS=<lib-A> Rscript _validation/golden/refactor_identity_check.R baseline --ref=A [--force]
#   R_LIBS=<lib-B> Rscript _validation/golden/refactor_identity_check.R check --ref=A \
#     [--tol-slots=grad,hess,bhhh,scores,tastes | none] [--verbose]
# Options: --ref=<name> names the baseline (refactor_baseline_<name>.rds);
#   --tol-slots replaces the tolerance group's slots (none: everything
#   identical); --out=<dir> is where baselines, reports and call logs go, by
#   default _validation/output/golden (ignored by git); --force overwrites an
#   existing baseline; --verbose prints every entry.
#
# Outputs never go into a package source tree's tracked files: an output
# directory inside a package tree (an ancestor holds a DESCRIPTION) is used
# only where git ignores it (`git check-ignore`), and the harness refuses
# otherwise, or when it cannot tell.
#
# Batteries: inline mnl, mxl, nl, errors, mxl_gen, mxl_wt, mxl_panel,
# mxl_edge, mxl_gen_pred, mxl_gen_wide, mxl_gen_long; then, from the files
# here, in this order (later files use earlier files' helpers):
# prep_batteries.R (prep, cube), p2_batteries.R (mnl_x, nl_x, idx_dbl,
# err_p2, hb, fits), p3_batteries.R (pred_x, pred_panel, pred_dbl, err_p3,
# pred_chunk), p4_batteries.R (mxl_asc, mxl_meat, mxl_vcov, mxl_edge_asc),
# p6b_batteries.R (scale_kern, scale_fit, scale_cols, scale_mnl),
# p6d_batteries.R (p6d_diag, p6d_coord, p6d_vcov, p6d_new),
# start_batteries.R (start_unit), optim_batteries.R (optim_ro, optim_fit),
# uasc_batteries.R (uasc), shares_batteries.R (pred_shares),
# p5a_batteries.R (ucon, polish), p5b_batteries.R (pdelta, profile),
# p8_batteries.R (p8_mnl, p8_nl, p8_edge, p8_fit), p8b_batteries.R
# (p8b_mnl, p8b_nl, p8b_edge, p8b_mxl, p8b_cl) and p8c_batteries.R
# (p8c_grad, p8c_polish, p8c_mxl).
# blp_batteries.R (blp), sourced first: every BLP call of the batteries goes
# through its golden_blp(), whose leaves (class "golden_blp") form the solver
# tier: the share inversion may move within its tolerance (BLP engine plan,
# user decision U1), so such a leaf passes if it is identical to the
# baseline's or meets the tier's rules (gb_tier(): the residual recomputed
# through the untouched predict kernels within tol, delta within the
# conditioning-scaled bound of a reference solution, each component's
# reference held, no regression where the baseline met its tol, an error
# exactly where the call's inputs fail one of the engine's checks or the call
# site expects one), and the run prints the tier's report.
# `check --bitwise` compares every leaf with identical(num.eq = FALSE,
# single.NA = FALSE), which also tells +0 from -0 and NaN payloads apart
# (identical()'s defaults do not); the tolerance group is unchanged.
# A new battery: a file here that defines run_<name>(), sourced and added to
# `results` below. A build that lacks what a battery calls records the error
# as that battery's output, so a battery can be added before the change it
# covers. Results at more than one OpenMP thread are not bit-stable (partial
# sums combine in varying order), so the harness sets one thread.
# ============================================================================

.lib <- strsplit(Sys.getenv("R_LIBS"), .Platform$path.sep, fixed = TRUE)[[1]]
.lib <- .lib[nzchar(.lib)]
if (!length(.lib)) {
  stop("Set R_LIBS=<lib> so the build under test is explicit, e.g. ",
       "R_LIBS=<lib> Rscript _validation/golden/refactor_identity_check.R check --ref=main")
}
.lib <- normalizePath(.lib[1], mustWork = TRUE)
suppressMessages(library(choicer, lib.loc = .lib))
.pkg_path <- getNamespaceInfo("choicer", "path")
if (normalizePath(dirname(.pkg_path)) != .lib) {
  stop("choicer was loaded from ", .pkg_path, ", not from R_LIBS[1] = ", .lib)
}
cat(sprintf("choicer %s loaded from %s\n",
            as.character(utils::packageVersion("choicer", lib.loc = .lib)),
            .pkg_path))

# The identity harness deliberately exercises implementation kernels.  Keep
# the installed-package lookup explicit now that those kernels are internal.
mnl_loglik_gradient_parallel <- choicer:::mnl_loglik_gradient_parallel
mnl_predict <- choicer:::mnl_predict
mnl_predict_shares <- choicer:::mnl_predict_shares
mnl_loglik_hessian_parallel <- choicer:::mnl_loglik_hessian_parallel
mnl_elasticities_parallel <- choicer:::mnl_elasticities_parallel
mnl_diversion_ratios_parallel <- choicer:::mnl_diversion_ratios_parallel
mxl_loglik_gradient_parallel <- choicer:::mxl_loglik_gradient_parallel
mxl_hessian_parallel <- choicer:::mxl_hessian_parallel
mxl_bhhh_parallel <- choicer:::mxl_bhhh_parallel
mxl_scores_parallel <- choicer:::mxl_scores_parallel
mxl_conditional_tastes_parallel <- choicer:::mxl_conditional_tastes_parallel
mxl_predict <- choicer:::mxl_predict
mxl_logsum <- choicer:::mxl_logsum
mxl_predict_shares <- choicer:::mxl_predict_shares
mxl_elasticities_parallel <- choicer:::mxl_elasticities_parallel
mxl_diversion_ratios_parallel <- choicer:::mxl_diversion_ratios_parallel
build_var_mat <- choicer:::build_var_mat
jacobian_vech_Sigma <- choicer:::jacobian_vech_Sigma
nl_loglik_gradient_parallel <- choicer:::nl_loglik_gradient_parallel
nl_predict <- choicer:::nl_predict
nl_predict_shares <- choicer:::nl_predict_shares
nl_elasticities_parallel <- choicer:::nl_elasticities_parallel
nl_diversion_ratios_parallel <- choicer:::nl_diversion_ratios_parallel
nl_loglik_numeric_hessian <- choicer:::nl_loglik_numeric_hessian

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else ""
if (!mode %in% c("baseline", "check")) {
  stop("Usage: R_LIBS=<lib> Rscript _validation/golden/refactor_identity_check.R ",
       "baseline|check [--ref=<name>] [--tol-slots=<slots>|none] [--out=<dir>] ",
       "[--bitwise] [--force] [--verbose]")
}
force_overwrite <- "--force" %in% args
verbose <- "--verbose" %in% args
bitwise <- "--bitwise" %in% args

# The batteries are sourced from this script's directory.
.script_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) != 1L) stop("Run this file with Rscript.")
  dirname(normalizePath(f, mustWork = TRUE))
})
# The output directory, absolute whether or not it exists yet.
.abs_path <- function(p) {
  p <- path.expand(p)
  if (!grepl("^/", p)) p <- file.path(getwd(), p)
  rest <- character()
  while (!dir.exists(p)) { rest <- c(basename(p), rest); p <- dirname(p) }
  do.call(file.path, as.list(c(normalizePath(p), rest)))
}
.out_arg <- sub("^--out=", "", grep("^--out=", args, value = TRUE))
.out_dir <- .abs_path(if (length(.out_arg)) .out_arg[1] else
  file.path(dirname(.script_dir), "output", "golden"))
ref_arg <- sub("^--ref=", "", grep("^--ref=", args, value = TRUE))
if (length(ref_arg) && !grepl("^[A-Za-z0-9._-]+$", ref_arg[1])) {
  stop("--ref names a file: use letters, digits, '.', '_' and '-' only")
}
ref_name <- if (length(ref_arg)) paste0("_", ref_arg[1]) else ""
# Never write into a package source tree's tracked files: inside a package
# tree, only where git ignores each file the harness writes (git check-ignore
# exits 0 for an ignored path, 1 for one not ignored, 128 outside a
# repository; given several paths it exits 0 if any is ignored, so each is
# checked on its own).
local({
  d <- .out_dir
  repeat {
    if (file.exists(file.path(d, "DESCRIPTION"))) break
    if (dirname(d) == d) return(invisible())
    d <- dirname(d)
  }
  anc <- .out_dir
  while (!dir.exists(anc)) anc <- dirname(anc)
  outs <- c(paste0("refactor_baseline", ref_name, ".rds"),
            paste0("refactor_check_report", ref_name, ".csv"),
            sprintf("optim_calls%s_%s_%s.rds", ref_name, basename(.lib), mode))
  for (o in outs) {
    st <- suppressWarnings(tryCatch(
      system2("git", c("-C", shQuote(anc), "check-ignore", "-q", "--",
                       shQuote(file.path(.out_dir, o))),
              stdout = FALSE, stderr = FALSE),
      error = function(e) -1L))
    if (!identical(as.integer(st), 0L)) {
      stop("Refusing to write ", o, " into ", .out_dir, ": it lies inside the package tree ",
           d, " and git does not ignore it (git check-ignore exit ", st, "). ",
           "Pass --out= with an ignored directory or one outside the tree.")
    }
  }
})
dir.create(.out_dir, showWarnings = FALSE, recursive = TRUE)
cat("Outputs in", .out_dir, "\n")
baseline_path <- file.path(.out_dir, paste0("refactor_baseline", ref_name, ".rds"))
report_path <- file.path(.out_dir, paste0("refactor_check_report", ref_name, ".csv"))

getFromNamespace("set_num_threads", "choicer")(1)
cat("OpenMP max threads inside choicer:",
    choicer::thread_info()$omp_get_max_threads, "\n")

get_halton_normals <- getFromNamespace("get_halton_normals", "choicer")

# ---------------------------------------------------------------------------
# Data generators (all RNG under fixed seeds)
# ---------------------------------------------------------------------------

# Unbalanced choice sets over J global alternatives. Individual 1 always sees
# the full set so every alternative has positive predicted share (needed for
# the BLP contraction targets).
make_choice_data <- function(seed, N, J, K, include_outside_option) {
  set.seed(seed)
  M <- integer(N)
  alt_list <- vector("list", N)
  for (i in seq_len(N)) {
    if (i == 1) {
      alts <- seq_len(J)
    } else {
      alts <- sort(sample.int(J, sample(2:J, 1)))
    }
    M[i] <- length(alts)
    alt_list[[i]] <- alts
  }
  alt_idx <- unlist(alt_list)
  total <- sum(M)
  X <- matrix(rnorm(total * K), total, K)
  choice_idx <- integer(N)
  for (i in seq_len(N)) {
    lo <- if (include_outside_option) 0L else 1L
    choice_idx[i] <- sample(lo:M[i], 1L)
  }
  weights <- runif(N, 0.5, 2)
  list(X = X, alt_idx = alt_idx, choice_idx = choice_idx, M = M,
       weights = weights, J = J, K = K, N = N, total = total)
}

n_delta_free <- function(J, use_asc, ioo) {
  if (!use_asc) 0L else if (ioo) J else J - 1L
}

# ---------------------------------------------------------------------------
# MNL battery
# ---------------------------------------------------------------------------

run_mnl <- function() {
  out <- list()
  d <- make_choice_data(seed = 101, N = 40, J = 5, K = 3,
                        include_outside_option = FALSE)
  d_oo <- make_choice_data(seed = 102, N = 40, J = 5, K = 3,
                           include_outside_option = TRUE)
  for (use_asc in c(TRUE, FALSE)) {
    for (ioo in c(TRUE, FALSE)) {
      dat <- if (ioo) d_oo else d
      n_par <- dat$K + n_delta_free(dat$J, use_asc, ioo)
      for (rep in 1:2) {
        set.seed(1000 + 10 * rep + ioo)
        theta <- rnorm(n_par, sd = 0.4)
        key <- sprintf("mnl|asc=%d|oo=%d|rep=%d", use_asc, ioo, rep)
        res <- list()
        res$grad <- mnl_loglik_gradient_parallel(
          theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, dat$weights,
          use_asc, ioo)
        res$pred <- mnl_predict(theta, dat$X, dat$alt_idx, dat$M, use_asc, ioo)
        res$shares <- mnl_predict_shares(
          theta, dat$X, dat$alt_idx, dat$M, dat$weights, use_asc, ioo)
        res$hess <- mnl_loglik_hessian_parallel(
          theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, dat$weights,
          use_asc, ioo)
        res$elas <- mnl_elasticities_parallel(
          theta, dat$X, dat$alt_idx, dat$choice_idx, dat$M, dat$weights,
          elast_var_idx = 1L, use_asc, ioo)
        res$dr <- mnl_diversion_ratios_parallel(
          theta, dat$X, dat$alt_idx, dat$M, dat$weights, use_asc, ioo)
        if (use_asc && rep == 1) {
          beta <- theta[seq_len(dat$K)]
          res$blp <- golden_blp("blp_contraction", list(
            rep(0, dat$J), res$shares, dat$X, beta, dat$alt_idx, dat$M,
            dat$weights, ioo), ref = gb_delta_star(theta, dat$K, dat$J, TRUE, ioo))
        }
        out[[key]] <- res
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# MXL battery
# ---------------------------------------------------------------------------

run_mxl <- function() {
  out <- list()
  N <- 30; J <- 4; K_x <- 2; K_w <- 2; S <- 8
  d <- make_choice_data(seed = 201, N = N, J = J, K = K_x,
                        include_outside_option = FALSE)
  d_oo <- make_choice_data(seed = 202, N = N, J = J, K = K_x,
                           include_outside_option = TRUE)
  set.seed(203)
  eta <- get_halton_normals(S, N, K_w)
  rc_dist_mixes <- list(normal = c(0L, 0L), mixed = c(0L, 1L),
                        lognormal = c(1L, 1L))

  for (use_asc in c(TRUE, FALSE)) {
    for (ioo in c(TRUE, FALSE)) {
      dat <- if (ioo) d_oo else d
      for (rc_corr in c(TRUE, FALSE)) {
        for (rc_mean in c(TRUE, FALSE)) {
          for (dist_name in names(rc_dist_mixes)) {
            rc_dist <- rc_dist_mixes[[dist_name]]
            for (w_type in c("row", "alt")) {
              set.seed(300 + 7 * ioo + 13 * rc_corr)
              W <- if (w_type == "row") {
                matrix(rnorm(dat$total * K_w), dat$total, K_w)
              } else {
                matrix(rnorm(J * K_w), J, K_w)
              }
              L_size <- if (rc_corr) K_w * (K_w + 1) / 2 else K_w
              n_par <- K_x + (if (rc_mean) K_w else 0L) + L_size +
                n_delta_free(J, use_asc, ioo)
              set.seed(400 + 3 * use_asc + 5 * rc_mean)
              theta <- rnorm(n_par, sd = 0.3)
              key <- sprintf(
                "mxl|asc=%d|oo=%d|corr=%d|mean=%d|dist=%s|W=%s",
                use_asc, ioo, rc_corr, rc_mean, dist_name, w_type)
              res <- list()
              res$grad <- mxl_loglik_gradient_parallel(
                theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                dat$weights, eta, rc_dist, rc_corr, rc_mean, use_asc, ioo)
              res$hess <- mxl_hessian_parallel(
                theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                dat$weights, eta, rc_dist, rc_corr, rc_mean, use_asc, ioo)
              res$bhhh <- mxl_bhhh_parallel(
                theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                dat$weights, eta, rc_dist, rc_corr, rc_mean, use_asc, ioo)
              # S0: weight-free per-situation scores (no weights argument).
              res$scores <- mxl_scores_parallel(
                theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                eta, rc_dist, rc_corr, rc_mean, use_asc, ioo)
              res$pred <- mxl_predict(
                theta, dat$X, W, dat$alt_idx, dat$M, eta, rc_dist,
                rc_corr, rc_mean, use_asc, ioo)
              res$logsum <- mxl_logsum(
                theta, dat$X, W, dat$alt_idx, dat$M, eta, rc_dist,
                rc_corr, rc_mean, use_asc, ioo)
              res$shares <- mxl_predict_shares(
                theta, dat$X, W, dat$alt_idx, dat$M, dat$weights, eta,
                rc_dist, rc_corr, rc_mean, use_asc, ioo)
              res$elas_fixed <- mxl_elasticities_parallel(
                theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                dat$weights, eta, rc_dist, elast_var_idx = 1L,
                is_random_coef = FALSE, rc_corr, rc_mean, use_asc, ioo)
              res$elas_random <- mxl_elasticities_parallel(
                theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                dat$weights, eta, rc_dist, elast_var_idx = 1L,
                is_random_coef = TRUE, rc_corr, rc_mean, use_asc, ioo)
              res$dr_fixed <- mxl_diversion_ratios_parallel(
                theta, dat$X, W, dat$alt_idx, dat$M, dat$weights, eta,
                rc_dist, elast_var_idx = 1L, is_random_coef = FALSE,
                rc_corr, rc_mean, use_asc, ioo)
              res$dr_random <- mxl_diversion_ratios_parallel(
                theta, dat$X, W, dat$alt_idx, dat$M, dat$weights, eta,
                rc_dist, elast_var_idx = 2L, is_random_coef = TRUE,
                rc_corr, rc_mean, use_asc, ioo)
              if (dist_name == "normal") {
                beta <- theta[seq_len(K_x)]
                mu <- if (rc_mean) theta[K_x + seq_len(K_w)] else rep(0, K_w)
                L_params <- theta[K_x + (if (rc_mean) K_w else 0L) +
                                    seq_len(L_size)]
                res$blp <- golden_blp("mxl_blp_contraction", list(
                  rep(0, J), res$shares, dat$X, W, beta, mu, L_params,
                  dat$alt_idx, dat$M, dat$weights, eta, rc_dist,
                  rc_corr, rc_mean, ioo), ref = gb_delta_star(
                    theta, K_x + (if (rc_mean) K_w else 0L) + L_size, J,
                    use_asc, ioo))
              }
              out[[key]] <- res
            }
          }
        }
      }
    }
  }

  # Cholesky utilities
  set.seed(500)
  for (rc_corr in c(TRUE, FALSE)) {
    L_size <- if (rc_corr) 6L else 3L
    L_params <- rnorm(L_size, sd = 0.5)
    key <- sprintf("mxl_chol|corr=%d", rc_corr)
    out[[key]] <- list(
      L = getFromNamespace("build_L_mat", "choicer")(L_params, 3L, rc_corr),
      Sigma = build_var_mat(L_params, 3L, rc_corr),
      jac = jacobian_vech_Sigma(L_params, 3L, rc_corr)
    )
  }
  out
}

# ---------------------------------------------------------------------------
# S0: shared cell builder for the extra MXL batteries. Same RNG recipe as the
# loop body of run_mxl(), so W and theta equal those of the matching
# store-mode cell there.
# ---------------------------------------------------------------------------

mxl_cell_inputs <- function(dat, J, K_x, K_w, use_asc, ioo, rc_corr, rc_mean,
                            w_type) {
  set.seed(300 + 7 * ioo + 13 * rc_corr)
  W <- if (w_type == "row") {
    matrix(rnorm(dat$total * K_w), dat$total, K_w)
  } else {
    matrix(rnorm(J * K_w), J, K_w)
  }
  L_size <- if (rc_corr) K_w * (K_w + 1) / 2 else K_w
  n_par <- K_x + (if (rc_mean) K_w else 0L) + L_size +
    n_delta_free(J, use_asc, ioo)
  set.seed(400 + 3 * use_asc + 5 * rc_mean)
  theta <- rnorm(n_par, sd = 0.3)
  list(W = W, theta = theta)
}

# ---------------------------------------------------------------------------
# S0: MXL estimation kernels in generate mode (on-the-fly Halton draws).
# Empty cube + gen_seed/gen_S/gen_scramble; full flag grid x scramble {0, 1}.
# ---------------------------------------------------------------------------

run_mxl_gen <- function() {
  out <- list()
  N <- 30; J <- 4; K_x <- 2; K_w <- 2; S <- 8
  d <- make_choice_data(seed = 201, N = N, J = J, K = K_x,
                        include_outside_option = FALSE)
  d_oo <- make_choice_data(seed = 202, N = N, J = J, K = K_x,
                           include_outside_option = TRUE)
  eta_empty <- array(0, dim = c(K_w, 0L, 0L))
  gen_seed <- 7L
  rc_dist_mixes <- list(normal = c(0L, 0L), mixed = c(0L, 1L),
                        lognormal = c(1L, 1L))

  for (scr in c(0L, 1L)) {
    for (use_asc in c(TRUE, FALSE)) {
      for (ioo in c(TRUE, FALSE)) {
        dat <- if (ioo) d_oo else d
        for (rc_corr in c(TRUE, FALSE)) {
          for (rc_mean in c(TRUE, FALSE)) {
            for (dist_name in names(rc_dist_mixes)) {
              rc_dist <- rc_dist_mixes[[dist_name]]
              for (w_type in c("row", "alt")) {
                inp <- mxl_cell_inputs(dat, J, K_x, K_w, use_asc, ioo,
                                       rc_corr, rc_mean, w_type)
                theta <- inp$theta; W <- inp$W
                key <- sprintf(
                  "mxl_gen|scr=%d|asc=%d|oo=%d|corr=%d|mean=%d|dist=%s|W=%s",
                  scr, use_asc, ioo, rc_corr, rc_mean, dist_name, w_type)
                res <- list()
                res$grad <- mxl_loglik_gradient_parallel(
                  theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                  dat$weights, eta_empty, rc_dist, rc_corr, rc_mean, use_asc,
                  ioo, gen_seed = gen_seed, gen_scramble = scr, gen_S = S)
                res$hess <- mxl_hessian_parallel(
                  theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                  dat$weights, eta_empty, rc_dist, rc_corr, rc_mean, use_asc,
                  ioo, gen_seed = gen_seed, gen_scramble = scr, gen_S = S)
                res$bhhh <- mxl_bhhh_parallel(
                  theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                  dat$weights, eta_empty, rc_dist, rc_corr, rc_mean, use_asc,
                  ioo, gen_seed = gen_seed, gen_scramble = scr, gen_S = S)
                res$scores <- mxl_scores_parallel(
                  theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                  eta_empty, rc_dist, rc_corr, rc_mean, use_asc, ioo,
                  gen_seed = gen_seed, gen_scramble = scr, gen_S = S)
                out[[key]] <- res
              }
            }
          }
        }
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# S0: gradient / Hessian / BHHH under a second, wider-range non-uniform weight
# vector (store mode, same draws as run_mxl) on a few flag combinations.
# ---------------------------------------------------------------------------

run_mxl_wt <- function() {
  out <- list()
  N <- 30; J <- 4; K_x <- 2; K_w <- 2; S <- 8
  d <- make_choice_data(seed = 201, N = N, J = J, K = K_x,
                        include_outside_option = FALSE)
  d_oo <- make_choice_data(seed = 202, N = N, J = J, K = K_x,
                           include_outside_option = TRUE)
  set.seed(203)
  eta <- get_halton_normals(S, N, K_w)
  set.seed(4242)
  w_alt <- runif(N, 0.2, 5)
  rc_dist_mixes <- list(normal = c(0L, 0L), mixed = c(0L, 1L),
                        lognormal = c(1L, 1L))
  cells <- list(
    list(asc = TRUE,  ioo = FALSE, corr = TRUE,  mean = TRUE,
         dist = "mixed",     W = "row"),
    list(asc = TRUE,  ioo = TRUE,  corr = TRUE,  mean = FALSE,
         dist = "normal",    W = "alt"),
    list(asc = FALSE, ioo = TRUE,  corr = FALSE, mean = TRUE,
         dist = "lognormal", W = "row"),
    list(asc = TRUE,  ioo = FALSE, corr = FALSE, mean = FALSE,
         dist = "normal",    W = "alt"),
    list(asc = TRUE,  ioo = TRUE,  corr = TRUE,  mean = TRUE,
         dist = "lognormal", W = "alt")
  )
  for (cc in cells) {
    dat <- if (cc$ioo) d_oo else d
    rc_dist <- rc_dist_mixes[[cc$dist]]
    inp <- mxl_cell_inputs(dat, J, K_x, K_w, cc$asc, cc$ioo, cc$corr,
                           cc$mean, cc$W)
    theta <- inp$theta; W <- inp$W
    key <- sprintf("mxl_wt|asc=%d|oo=%d|corr=%d|mean=%d|dist=%s|W=%s",
                   cc$asc, cc$ioo, cc$corr, cc$mean, cc$dist, cc$W)
    res <- list()
    res$grad <- mxl_loglik_gradient_parallel(
      theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M, w_alt, eta,
      rc_dist, cc$corr, cc$mean, cc$asc, cc$ioo)
    res$hess <- mxl_hessian_parallel(
      theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M, w_alt, eta,
      rc_dist, cc$corr, cc$mean, cc$asc, cc$ioo)
    res$bhhh <- mxl_bhhh_parallel(
      theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M, w_alt, eta,
      rc_dist, cc$corr, cc$mean, cc$asc, cc$ioo)
    out[[key]] <- res
  }
  out
}

# ---------------------------------------------------------------------------
# NL battery
# ---------------------------------------------------------------------------

run_nl <- function() {
  out <- list()
  N <- 40; J <- 5; K <- 3
  # nests: {1,2} -> 1, {3,4} -> 2, {5} -> 3 (singleton)
  nest_idx <- c(1L, 1L, 2L, 2L, 3L)
  n_nonsing <- 2L
  d <- make_choice_data(seed = 301, N = N, J = J, K = K,
                        include_outside_option = FALSE)
  d_oo <- make_choice_data(seed = 302, N = N, J = J, K = K,
                           include_outside_option = TRUE)
  for (use_asc in c(TRUE, FALSE)) {
    for (ioo in c(TRUE, FALSE)) {
      dat <- if (ioo) d_oo else d
      n_par <- K + n_nonsing + n_delta_free(J, use_asc, ioo)
      for (rep in 1:2) {
        set.seed(2000 + 10 * rep + ioo)
        theta <- rnorm(n_par, sd = 0.4)
        theta[K + seq_len(n_nonsing)] <- runif(n_nonsing, 0.6, 1.0)
        key <- sprintf("nl|asc=%d|oo=%d|rep=%d", use_asc, ioo, rep)
        res <- list()
        res$grad <- nl_loglik_gradient_parallel(
          theta, dat$X, dat$alt_idx, dat$choice_idx, nest_idx, dat$M,
          dat$weights, use_asc, ioo)
        res$pred <- nl_predict(
          theta, dat$X, dat$alt_idx, dat$M, nest_idx, use_asc, ioo)
        res$shares <- nl_predict_shares(
          theta, dat$X, dat$alt_idx, dat$M, dat$weights, nest_idx,
          use_asc, ioo)
        res$elas <- nl_elasticities_parallel(
          theta, dat$X, dat$alt_idx, dat$choice_idx, nest_idx, dat$M,
          dat$weights, elast_var_idx = 1L, use_asc, ioo)
        res$dr <- nl_diversion_ratios_parallel(
          theta, dat$X, dat$alt_idx, nest_idx, dat$M, dat$weights,
          use_asc, ioo)
        if (rep == 1) {
          res$num_hess <- nl_loglik_numeric_hessian(
            theta, dat$X, dat$alt_idx, dat$choice_idx, nest_idx, dat$M,
            dat$weights, use_asc, ioo)
        }
        if (use_asc && rep == 1) {
          beta <- theta[seq_len(K)]
          lambda_full <- rep(1, 3)
          lambda_full[1:2] <- theta[K + 1:2]
          res$blp <- golden_blp("nl_blp_contraction", list(
            rep(0, J), res$shares, dat$X, beta, lambda_full, dat$alt_idx,
            nest_idx, dat$M, dat$weights, ioo),
            ref = gb_delta_star(theta, K + 2L, J, TRUE, ioo))
        }
        out[[key]] <- res
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Error-message subharness
# ---------------------------------------------------------------------------

run_errors <- function() {
  grab <- function(expr) tryCatch({ expr; "<<no error>>" },
                                  error = function(e) conditionMessage(e))
  out <- list()
  d <- make_choice_data(seed = 901, N = 10, J = 4, K = 2,
                        include_outside_option = FALSE)
  nest_idx <- c(1L, 1L, 2L, 2L)
  set.seed(902)
  eta <- get_halton_normals(4, d$N, 2)
  W <- matrix(rnorm(d$total * 2), d$total, 2)

  # MNL: ASC expected but theta has only beta
  out$mnl_missing_asc <- grab(mnl_loglik_gradient_parallel(
    rnorm(2), d$X, d$alt_idx, d$choice_idx, d$M, d$weights, TRUE, FALSE))
  # MXL: theta too short for L block (gradient path)
  out$mxl_short_L <- grab(mxl_loglik_gradient_parallel(
    rnorm(2), d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta,
    c(0L, 0L), FALSE, FALSE, FALSE, FALSE))
  # MXL: theta too short for delta block (gradient path, use_asc)
  out$mxl_short_delta <- grab(mxl_loglik_gradient_parallel(
    rnorm(4), d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta,
    c(0L, 0L), FALSE, FALSE, TRUE, FALSE))
  # MXL: delta stop on a predict-family path
  out$mxl_predict_short_delta <- grab(mxl_predict(
    rnorm(4), d$X, W, d$alt_idx, d$M, eta, c(0L, 0L),
    FALSE, FALSE, TRUE, FALSE))
  # MXL: bad rc_dist length
  out$mxl_bad_rcdist <- grab(mxl_loglik_gradient_parallel(
    rnorm(7), d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta,
    c(0L), FALSE, FALSE, TRUE, FALSE))
  # MXL: mismatched eta draws (wrong N)
  set.seed(903)
  eta_bad <- get_halton_normals(4, d$N + 1, 2)
  out$mxl_bad_eta <- grab(mxl_loglik_gradient_parallel(
    rnorm(7), d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta_bad,
    c(0L, 0L), FALSE, FALSE, TRUE, FALSE))
  # NL: non-positive lambda
  out$nl_bad_lambda <- grab(nl_loglik_gradient_parallel(
    c(rnorm(2), -0.5, 0.8, rnorm(3)), d$X, d$alt_idx, d$choice_idx,
    nest_idx, d$M, d$weights, TRUE, FALSE))
  # NL: theta too short for lambda block
  out$nl_short <- grab(nl_loglik_gradient_parallel(
    rnorm(3), d$X, d$alt_idx, d$choice_idx, nest_idx, d$M, d$weights,
    FALSE, FALSE))

  # --- Always-on validation added with the parser/validator consolidation ---
  set.seed(904)
  # MNL: trailing parameters not parsed when use_asc = FALSE
  out$mnl_too_long <- grab(mnl_loglik_gradient_parallel(
    rnorm(4), d$X, d$alt_idx, d$choice_idx, d$M, d$weights, FALSE, FALSE))
  # MNL: weights length mismatch
  out$mnl_weights_len <- grab(mnl_loglik_gradient_parallel(
    rnorm(5), d$X, d$alt_idx, d$choice_idx, d$M, c(d$weights, 1), TRUE, FALSE))
  # MNL: choice_idx length mismatch
  out$mnl_choice_len <- grab(mnl_loglik_gradient_parallel(
    rnorm(5), d$X, d$alt_idx, d$choice_idx[-1], d$M, d$weights, TRUE, FALSE))
  # MNL: sum(M) does not match nrow(X)
  M_bad <- d$M; M_bad[1] <- M_bad[1] + 1L
  out$mnl_m_rows <- grab(mnl_loglik_gradient_parallel(
    rnorm(5), d$X, d$alt_idx, d$choice_idx, M_bad, d$weights, TRUE, FALSE))
  # MNL: alt_idx must be 1-based
  alt_bad <- d$alt_idx; alt_bad[1] <- 0L
  out$mnl_alt_zero <- grab(mnl_predict(
    rnorm(5), d$X, alt_bad, d$M, TRUE, FALSE))
  # MNL: delta block does not cover all referenced alternatives
  out$mnl_delta_coverage <- grab(mnl_predict(
    rnorm(4), d$X, d$alt_idx, d$M, TRUE, FALSE))  # delta covers 3 of 4 alts
  # MXL: trailing parameters not parsed when use_asc = FALSE
  out$mxl_too_long <- grab(mxl_predict(
    rnorm(6), d$X, W, d$alt_idx, d$M, eta, c(0L, 0L),
    FALSE, FALSE, FALSE, FALSE))
  # MXL: predict family now rejects mismatched eta draws too
  set.seed(905)
  eta_bad2 <- get_halton_normals(4, d$N + 1, 2)
  out$mxl_predict_bad_eta <- grab(mxl_predict(
    rnorm(4), d$X, W, d$alt_idx, d$M, eta_bad2, c(0L, 0L),
    FALSE, FALSE, FALSE, FALSE))
  # MXL: alt-level W with too few rows for the referenced alternatives
  W_bad <- matrix(rnorm(4), 2, 2)
  out$mxl_w_rows <- grab(mxl_predict(
    rnorm(4), d$X, W_bad, d$alt_idx, d$M, eta, c(0L, 0L),
    FALSE, FALSE, FALSE, FALSE))
  # NL: nest_idx shorter than the number of referenced alternatives
  out$nl_nest_short <- grab(nl_loglik_gradient_parallel(
    c(rnorm(2), 0.8), d$X, d$alt_idx, d$choice_idx, c(1L, 1L, 2L),
    d$M, d$weights, FALSE, FALSE))
  # NL: trailing parameters not parsed when use_asc = FALSE
  out$nl_too_long <- grab(nl_loglik_gradient_parallel(
    c(rnorm(2), 0.8, 0.9, 1.1), d$X, d$alt_idx, d$choice_idx, nest_idx,
    d$M, d$weights, FALSE, FALSE))

  # --- S0 additions: the four MXL estimation kernels (cross-section) ---
  # Valid theta for (K_x = 2, K_w = 2, rc_corr = FALSE, rc_mean = FALSE,
  # use_asc = TRUE, J = 4): 2 + 2 + 3 = 7 parameters.
  set.seed(906)
  th7 <- rnorm(7)
  rcd <- c(0L, 0L)
  eta_empty <- array(0, dim = c(2L, 0L, 0L))
  choice_bad <- d$choice_idx; choice_bad[4] <- 99L
  # Cube slices != N (message must be kept for the cross-section)
  out$mxl_hess_bad_eta <- grab(mxl_hessian_parallel(
    th7, d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta_bad,
    rcd, FALSE, FALSE, TRUE, FALSE))
  out$mxl_bhhh_bad_eta <- grab(mxl_bhhh_parallel(
    th7, d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta_bad,
    rcd, FALSE, FALSE, TRUE, FALSE))
  out$mxl_scores_bad_eta <- grab(mxl_scores_parallel(
    th7, d$X, W, d$alt_idx, d$choice_idx, d$M, eta_bad,
    rcd, FALSE, FALSE, TRUE, FALSE))
  # Generate mode with gen_S = 0
  out$mxl_grad_gen_S0 <- grab(mxl_loglik_gradient_parallel(
    th7, d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta_empty,
    rcd, FALSE, FALSE, TRUE, FALSE, gen_seed = 7L, gen_S = 0L))
  out$mxl_hess_gen_S0 <- grab(mxl_hessian_parallel(
    th7, d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta_empty,
    rcd, FALSE, FALSE, TRUE, FALSE, gen_seed = 7L, gen_S = 0L))
  out$mxl_bhhh_gen_S0 <- grab(mxl_bhhh_parallel(
    th7, d$X, W, d$alt_idx, d$choice_idx, d$M, d$weights, eta_empty,
    rcd, FALSE, FALSE, TRUE, FALSE, gen_seed = 7L, gen_S = 0L))
  out$mxl_scores_gen_S0 <- grab(mxl_scores_parallel(
    th7, d$X, W, d$alt_idx, d$choice_idx, d$M, eta_empty,
    rcd, FALSE, FALSE, TRUE, FALSE, gen_seed = 7L, gen_S = 0L))
  # Chosen index out of range
  out$mxl_grad_bad_choice <- grab(mxl_loglik_gradient_parallel(
    th7, d$X, W, d$alt_idx, choice_bad, d$M, d$weights, eta,
    rcd, FALSE, FALSE, TRUE, FALSE))
  out$mxl_hess_bad_choice <- grab(mxl_hessian_parallel(
    th7, d$X, W, d$alt_idx, choice_bad, d$M, d$weights, eta,
    rcd, FALSE, FALSE, TRUE, FALSE))
  out$mxl_bhhh_bad_choice <- grab(mxl_bhhh_parallel(
    th7, d$X, W, d$alt_idx, choice_bad, d$M, d$weights, eta,
    rcd, FALSE, FALSE, TRUE, FALSE))
  out$mxl_scores_bad_choice <- grab(mxl_scores_parallel(
    th7, d$X, W, d$alt_idx, choice_bad, d$M, eta,
    rcd, FALSE, FALSE, TRUE, FALSE))
  # Weights length mismatch
  out$mxl_grad_weights_len <- grab(mxl_loglik_gradient_parallel(
    th7, d$X, W, d$alt_idx, d$choice_idx, d$M, c(d$weights, 1), eta,
    rcd, FALSE, FALSE, TRUE, FALSE))
  out$mxl_hess_weights_len <- grab(mxl_hessian_parallel(
    th7, d$X, W, d$alt_idx, d$choice_idx, d$M, c(d$weights, 1), eta,
    rcd, FALSE, FALSE, TRUE, FALSE))
  out$mxl_bhhh_weights_len <- grab(mxl_bhhh_parallel(
    th7, d$X, W, d$alt_idx, d$choice_idx, d$M, c(d$weights, 1), eta,
    rcd, FALSE, FALSE, TRUE, FALSE))
  # Theta too short for the L block (scores path)
  out$mxl_scores_short_L <- grab(mxl_scores_parallel(
    rnorm(2), d$X, W, d$alt_idx, d$choice_idx, d$M, eta,
    rcd, FALSE, FALSE, FALSE, FALSE))
  out
}

# ---------------------------------------------------------------------------
# P1 (kernel-efficiency session): panel cells. The five unit kernels
# (gradient, Hessian, BHHH, scores, conditional tastes) with a decision-maker
# structure Ti (T_u in 1..5, unbalanced choice sets, weights constant within
# each decision maker) in store mode and in generate mode with both scramble
# settings, over the full asc x oo x corr x mean x dist x W grid. Also the
# conditional tastes of the cross-sectional cells (Ti = NULL).
# ---------------------------------------------------------------------------

make_panel_data <- function(seed, U, J, K_x, ioo) {
  set.seed(seed)
  Ti <- sample(c(1:5, sample.int(5L, U - 5L, replace = TRUE)))
  N <- sum(Ti)
  alts <- lapply(seq_len(N), function(t) {
    if (t == 1L) seq_len(J) else sort(sample.int(J, sample(2:J, 1L)))
  })
  M <- lengths(alts)
  alt_idx <- unlist(alts)
  total <- sum(M)
  X <- matrix(rnorm(total * K_x), total, K_x)
  lo <- if (ioo) 0L else 1L
  choice_idx <- vapply(M, function(m) sample(lo:m, 1L), integer(1))
  w_u <- runif(U, 0.5, 2)
  list(X = X, alt_idx = alt_idx, choice_idx = as.integer(choice_idx), M = M,
       weights = rep(w_u, Ti), Ti = as.integer(Ti), U = U, N = N,
       total = total, J = J, K = K_x)
}

run_mxl_panel <- function() {
  out <- list()
  U <- 12; J <- 4; K_x <- 2; K_w <- 2; S <- 8
  d <- make_panel_data(seed = 701, U = U, J = J, K_x = K_x, ioo = FALSE)
  d_oo <- make_panel_data(seed = 702, U = U, J = J, K_x = K_x, ioo = TRUE)
  set.seed(703)
  eta <- get_halton_normals(S, U, K_w)
  eta_empty <- array(0, dim = c(K_w, 0L, 0L))
  rc_dist_mixes <- list(normal = c(0L, 0L), mixed = c(0L, 1L),
                        lognormal = c(1L, 1L))
  modes <- list(store = list(eta = eta, gs = -1L, sc = 1L, gS = 0L),
                gen0 = list(eta = eta_empty, gs = 11L, sc = 0L, gS = S),
                gen1 = list(eta = eta_empty, gs = 11L, sc = 1L, gS = S))
  for (mode_name in names(modes)) {
    md <- modes[[mode_name]]
    for (use_asc in c(TRUE, FALSE)) {
      for (ioo in c(TRUE, FALSE)) {
        dat <- if (ioo) d_oo else d
        for (rc_corr in c(TRUE, FALSE)) {
          for (rc_mean in c(TRUE, FALSE)) {
            for (dist_name in names(rc_dist_mixes)) {
              rc_dist <- rc_dist_mixes[[dist_name]]
              for (w_type in c("row", "alt")) {
                inp <- mxl_cell_inputs(dat, J, K_x, K_w, use_asc, ioo,
                                       rc_corr, rc_mean, w_type)
                theta <- inp$theta; W <- inp$W
                key <- sprintf(
                  "mxl_panel|%s|asc=%d|oo=%d|corr=%d|mean=%d|dist=%s|W=%s",
                  mode_name, use_asc, ioo, rc_corr, rc_mean, dist_name, w_type)
                res <- list()
                res$grad <- mxl_loglik_gradient_parallel(
                  theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                  dat$weights, md$eta, rc_dist, rc_corr, rc_mean, use_asc, ioo,
                  gen_seed = md$gs, gen_scramble = md$sc, gen_S = md$gS,
                  Ti = dat$Ti)
                res$hess <- mxl_hessian_parallel(
                  theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                  dat$weights, md$eta, rc_dist, rc_corr, rc_mean, use_asc, ioo,
                  gen_seed = md$gs, gen_scramble = md$sc, gen_S = md$gS,
                  Ti = dat$Ti)
                res$bhhh <- mxl_bhhh_parallel(
                  theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                  dat$weights, md$eta, rc_dist, rc_corr, rc_mean, use_asc, ioo,
                  gen_seed = md$gs, gen_scramble = md$sc, gen_S = md$gS,
                  Ti = dat$Ti)
                res$scores <- mxl_scores_parallel(
                  theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                  md$eta, rc_dist, rc_corr, rc_mean, use_asc, ioo,
                  gen_seed = md$gs, gen_scramble = md$sc, gen_S = md$gS,
                  Ti = dat$Ti)
                res$tastes <- mxl_conditional_tastes_parallel(
                  theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                  md$eta, rc_dist, rc_corr, rc_mean, use_asc, ioo,
                  gen_seed = md$gs, gen_scramble = md$sc, gen_S = md$gS,
                  Ti = dat$Ti)
                out[[key]] <- res
              }
            }
          }
        }
      }
    }
  }

  # Cross-sectional conditional tastes on the run_mxl() cells (store mode).
  N <- 30
  dx <- make_choice_data(seed = 201, N = N, J = J, K = K_x,
                         include_outside_option = FALSE)
  dx_oo <- make_choice_data(seed = 202, N = N, J = J, K = K_x,
                            include_outside_option = TRUE)
  set.seed(203)
  eta_x <- get_halton_normals(S, N, K_w)
  for (use_asc in c(TRUE, FALSE)) {
    for (ioo in c(TRUE, FALSE)) {
      dat <- if (ioo) dx_oo else dx
      for (rc_corr in c(TRUE, FALSE)) {
        for (rc_mean in c(TRUE, FALSE)) {
          for (dist_name in names(rc_dist_mixes)) {
            rc_dist <- rc_dist_mixes[[dist_name]]
            for (w_type in c("row", "alt")) {
              inp <- mxl_cell_inputs(dat, J, K_x, K_w, use_asc, ioo,
                                     rc_corr, rc_mean, w_type)
              key <- sprintf(
                "mxl_xtastes|asc=%d|oo=%d|corr=%d|mean=%d|dist=%s|W=%s",
                use_asc, ioo, rc_corr, rc_mean, dist_name, w_type)
              out[[key]] <- list(tastes = mxl_conditional_tastes_parallel(
                inp$theta, dat$X, inp$W, dat$alt_idx, dat$choice_idx, dat$M,
                eta_x, rc_dist, rc_corr, rc_mean, use_asc, ioo))
            }
          }
        }
      }
    }
  }
  out
}

# ---------------------------------------------------------------------------
# P1: extreme-utility cells (item 1 of the kernel-efficiency session). One
# chosen alternative sits `gap` utils below its competitors (through X[, 1]
# with beta_1 = 0.5), in a cross-section and inside a multi-situation decision
# maker, plus a fat-tailed cell (Cauchy covariates). Reported per kernel. The
# log(P(chosen)) primitive loses precision once P is subnormal (gap < -708)
# and gives -Inf past about -745, so these entries are expected to move when
# the kernels switch to the stable log probability.
# ---------------------------------------------------------------------------

run_mxl_edge <- function() {
  out <- list()
  U <- 10; J <- 4; K_x <- 2; K_w <- 2; S <- 8
  rc_dist <- c(0L, 1L)
  base <- make_panel_data(seed = 801, U = U, J = J, K_x = K_x, ioo = TRUE)
  set.seed(802)
  W <- matrix(rnorm(base$total * K_w), base$total, K_w)
  # beta (2), mu (2), L (3, corr), delta (J = 4 with the outside option)
  theta <- c(0.5, -0.3, 0.2, -0.4, log(0.6), 0.2, log(0.5), rnorm(J, sd = 0.3))
  row_off <- c(0L, cumsum(base$M))
  u_edge <- which(base$Ti >= 2L)[2L]
  t_edge <- sum(base$Ti[seq_len(u_edge - 1L)]) + 1L
  for (gap in c(-36, -720, -740, -1e5)) {
    for (panel in c(TRUE, FALSE)) {
      dat <- base
      if (dat$choice_idx[t_edge] == 0L) dat$choice_idx[t_edge] <- 1L
      rows <- (row_off[t_edge] + 1L):row_off[t_edge + 1L]
      dat$X[rows, 1L] <- 0
      dat$X[rows[dat$choice_idx[t_edge]], 1L] <- gap / theta[1L]
      Ti <- if (panel) dat$Ti else NULL
      n_units <- if (panel) U else dat$N
      eta <- get_halton_normals(S, n_units, K_w)
      w <- if (panel) dat$weights else runif(dat$N, 0.5, 2)
      key <- sprintf("mxl_edge|gap=%g|%s", gap, if (panel) "panel" else "xsec")
      out[[key]] <- mxl_edge_cell(theta, dat, W, w, eta, rc_dist, Ti)
    }
  }
  # Fat tails: Cauchy covariates, so some choice probabilities underflow.
  set.seed(803)
  dat <- base
  dat$X[, 1L] <- rcauchy(dat$total, scale = 20)
  eta <- get_halton_normals(S, U, K_w)
  out[["mxl_edge|cauchy|panel"]] <-
    mxl_edge_cell(theta, dat, W, dat$weights, eta, rc_dist, dat$Ti)
  out
}

mxl_edge_cell <- function(theta, dat, W, w, eta, rc_dist, Ti) {
  res <- list()
  res$grad <- mxl_loglik_gradient_parallel(
    theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M, w, eta, rc_dist,
    TRUE, TRUE, TRUE, TRUE, Ti = Ti)
  res$hess <- mxl_hessian_parallel(
    theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M, w, eta, rc_dist,
    TRUE, TRUE, TRUE, TRUE, Ti = Ti)
  res$bhhh <- mxl_bhhh_parallel(
    theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M, w, eta, rc_dist,
    TRUE, TRUE, TRUE, TRUE, Ti = Ti)
  res$scores <- mxl_scores_parallel(
    theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M, eta, rc_dist,
    TRUE, TRUE, TRUE, TRUE, Ti = Ti)
  res$tastes <- mxl_conditional_tastes_parallel(
    theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M, eta, rc_dist,
    TRUE, TRUE, TRUE, TRUE, Ti = Ti)
  res
}

# ---------------------------------------------------------------------------
# Scale audit, Phase 1 (Halton hybrid): generate-mode draws through every
# kernel that consumes them.
#   mxl_gen_pred: the prediction family in generate mode (the S0 battery
#     covers only the estimation kernels), both scramble settings, over
#     asc x oo x dist x W with correlated random coefficients and means.
#   mxl_gen_wide: K_w = 13, so dimensions 12 and 13 use bases 37 and 41,
#     past the compile-time-constant bases (2..31) of the block generator;
#     estimation and prediction kernels.
#   mxl_gen_long: S = 1000 draws, global Halton indices up to 3e4 (15 binary
#     digits of odometer carries), cross-section and panel.
# All draws come from HaltonGen; the estimation slots join the tolerance
# group for later class-T phases, everything else must stay identical().
# ---------------------------------------------------------------------------

gen_pred_outputs <- function(theta, dat, W, rc_dist, rc_corr, rc_mean, use_asc,
                             ioo, gen_seed, scr, S, K_x, K_w, L_size, J,
                             with_blp = TRUE) {
  eta_empty <- array(0, dim = c(K_w, 0L, 0L))
  g <- list(gen_seed = gen_seed, gen_scramble = scr, gen_S = S)
  res <- list()
  res$pred <- do.call(mxl_predict, c(list(
    theta, dat$X, W, dat$alt_idx, dat$M, eta_empty, rc_dist, rc_corr,
    rc_mean, use_asc, ioo), g))
  res$logsum <- do.call(mxl_logsum, c(list(
    theta, dat$X, W, dat$alt_idx, dat$M, eta_empty, rc_dist, rc_corr,
    rc_mean, use_asc, ioo), g))
  res$shares <- do.call(mxl_predict_shares, c(list(
    theta, dat$X, W, dat$alt_idx, dat$M, dat$weights, eta_empty, rc_dist,
    rc_corr, rc_mean, use_asc, ioo), g))
  res$elas_fixed <- do.call(mxl_elasticities_parallel, c(list(
    theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M, dat$weights,
    eta_empty, rc_dist, elast_var_idx = 1L, is_random_coef = FALSE, rc_corr,
    rc_mean, use_asc, ioo), g))
  res$elas_random <- do.call(mxl_elasticities_parallel, c(list(
    theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M, dat$weights,
    eta_empty, rc_dist, elast_var_idx = 1L, is_random_coef = TRUE, rc_corr,
    rc_mean, use_asc, ioo), g))
  res$dr_fixed <- do.call(mxl_diversion_ratios_parallel, c(list(
    theta, dat$X, W, dat$alt_idx, dat$M, dat$weights, eta_empty, rc_dist,
    elast_var_idx = 1L, is_random_coef = FALSE, rc_corr, rc_mean, use_asc,
    ioo), g))
  res$dr_random <- do.call(mxl_diversion_ratios_parallel, c(list(
    theta, dat$X, W, dat$alt_idx, dat$M, dat$weights, eta_empty, rc_dist,
    elast_var_idx = 2L, is_random_coef = TRUE, rc_corr, rc_mean, use_asc,
    ioo), g))
  if (with_blp && all(rc_dist == 0L)) {
    beta <- theta[seq_len(K_x)]
    mu <- if (rc_mean) theta[K_x + seq_len(K_w)] else rep(0, K_w)
    L_params <- theta[K_x + (if (rc_mean) K_w else 0L) + seq_len(L_size)]
    res$blp <- golden_blp("mxl_blp_contraction", c(list(
      rep(0, J), res$shares, dat$X, W, beta, mu, L_params, dat$alt_idx,
      dat$M, dat$weights, eta_empty, rc_dist, rc_corr, rc_mean, ioo), g),
      ref = gb_delta_star(theta, K_x + (if (rc_mean) K_w else 0L) + L_size, J,
                          use_asc, ioo))
  }
  res
}

gen_est_outputs <- function(theta, dat, W, rc_dist, rc_corr, rc_mean, use_asc,
                            ioo, gen_seed, scr, S, K_w, Ti = NULL,
                            hess = TRUE) {
  eta_empty <- array(0, dim = c(K_w, 0L, 0L))
  g <- list(gen_seed = gen_seed, gen_scramble = scr, gen_S = S, Ti = Ti)
  args_w <- list(theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                 dat$weights, eta_empty, rc_dist, rc_corr, rc_mean, use_asc,
                 ioo)
  args_nw <- list(theta, dat$X, W, dat$alt_idx, dat$choice_idx, dat$M,
                  eta_empty, rc_dist, rc_corr, rc_mean, use_asc, ioo)
  res <- list()
  res$grad <- do.call(mxl_loglik_gradient_parallel, c(args_w, g))
  if (hess) {
    res$hess <- do.call(mxl_hessian_parallel, c(args_w, g))
    res$bhhh <- do.call(mxl_bhhh_parallel, c(args_w, g))
  }
  res$scores <- do.call(mxl_scores_parallel, c(args_nw, g))
  res$tastes <- do.call(mxl_conditional_tastes_parallel, c(args_nw, g))
  res
}

run_mxl_gen_pred <- function() {
  out <- list()
  N <- 30; J <- 4; K_x <- 2; K_w <- 2; S <- 8
  d <- make_choice_data(seed = 201, N = N, J = J, K = K_x,
                        include_outside_option = FALSE)
  d_oo <- make_choice_data(seed = 202, N = N, J = J, K = K_x,
                           include_outside_option = TRUE)
  rc_dist_mixes <- list(normal = c(0L, 0L), mixed = c(0L, 1L))
  rc_corr <- TRUE; rc_mean <- TRUE
  L_size <- K_w * (K_w + 1) / 2
  for (scr in c(0L, 1L)) {
    for (use_asc in c(TRUE, FALSE)) {
      for (ioo in c(TRUE, FALSE)) {
        dat <- if (ioo) d_oo else d
        for (dist_name in names(rc_dist_mixes)) {
          for (w_type in c("row", "alt")) {
            inp <- mxl_cell_inputs(dat, J, K_x, K_w, use_asc, ioo, rc_corr,
                                   rc_mean, w_type)
            key <- sprintf("mxl_gen_pred|scr=%d|asc=%d|oo=%d|dist=%s|W=%s",
                           scr, use_asc, ioo, dist_name, w_type)
            out[[key]] <- gen_pred_outputs(
              inp$theta, dat, inp$W, rc_dist_mixes[[dist_name]], rc_corr,
              rc_mean, use_asc, ioo, gen_seed = 7L, scr = scr, S = S,
              K_x = K_x, K_w = K_w, L_size = L_size, J = J)
          }
        }
      }
    }
  }
  out
}

run_mxl_gen_wide <- function() {
  out <- list()
  N <- 30; J <- 4; K_x <- 2; K_w <- 13L; S <- 16
  d <- make_choice_data(seed = 211, N = N, J = J, K = K_x,
                        include_outside_option = FALSE)
  d_oo <- make_choice_data(seed = 212, N = N, J = J, K = K_x,
                           include_outside_option = TRUE)
  rc_dist_mixes <- list(normal = rep(0L, K_w),
                        mixed = rep(c(0L, 1L), length.out = K_w))
  use_asc <- TRUE; rc_mean <- TRUE
  for (scr in c(0L, 1L)) {
    for (ioo in c(TRUE, FALSE)) {
      dat <- if (ioo) d_oo else d
      for (rc_corr in c(TRUE, FALSE)) {
        for (dist_name in names(rc_dist_mixes)) {
          rc_dist <- rc_dist_mixes[[dist_name]]
          L_size <- if (rc_corr) K_w * (K_w + 1) / 2 else K_w
          set.seed(900 + 7 * ioo + 13 * rc_corr)
          W <- matrix(rnorm(dat$total * K_w, sd = 0.5), dat$total, K_w)
          n_par <- K_x + K_w + L_size + n_delta_free(J, use_asc, ioo)
          theta <- rnorm(n_par, sd = 0.2)
          # Moderate random-coefficient scale: L diagonal log(0.3) +- noise
          diag_pos <- if (rc_corr) cumsum(seq_len(K_w)) else seq_len(K_w)
          theta[K_x + K_w + diag_pos] <- log(0.3) + theta[K_x + K_w + diag_pos]
          key <- sprintf("mxl_gen_wide|scr=%d|oo=%d|corr=%d|dist=%s",
                         scr, ioo, rc_corr, dist_name)
          out[[key]] <- c(
            gen_est_outputs(theta, dat, W, rc_dist, rc_corr, rc_mean,
                            use_asc, ioo, gen_seed = 23L, scr = scr, S = S,
                            K_w = K_w),
            gen_pred_outputs(theta, dat, W, rc_dist, rc_corr, rc_mean,
                             use_asc, ioo, gen_seed = 23L, scr = scr, S = S,
                             K_x = K_x, K_w = K_w, L_size = L_size, J = J,
                             with_blp = FALSE))
        }
      }
    }
  }
  out
}

run_mxl_gen_long <- function() {
  out <- list()
  N <- 30; J <- 4; K_x <- 2; K_w <- 3L; S <- 1000
  d <- make_choice_data(seed = 221, N = N, J = J, K = K_x,
                        include_outside_option = FALSE)
  d_oo <- make_choice_data(seed = 222, N = N, J = J, K = K_x,
                           include_outside_option = TRUE)
  p_oo <- make_panel_data(seed = 223, U = 12, J = J, K_x = K_x, ioo = TRUE)
  rc_dist <- c(0L, 1L, 0L)
  use_asc <- TRUE; rc_corr <- TRUE; rc_mean <- TRUE
  L_size <- K_w * (K_w + 1) / 2
  cells <- list(xsec = list(dat = d, ioo = FALSE, Ti = NULL),
                xsec_oo = list(dat = d_oo, ioo = TRUE, Ti = NULL),
                panel_oo = list(dat = p_oo, ioo = TRUE, Ti = p_oo$Ti))
  for (scr in c(0L, 1L)) {
    for (cell in names(cells)) {
      cc <- cells[[cell]]
      dat <- cc$dat
      set.seed(950 + nchar(cell))
      W <- matrix(rnorm(dat$total * K_w), dat$total, K_w)
      n_par <- K_x + K_w + L_size + n_delta_free(J, use_asc, cc$ioo)
      theta <- rnorm(n_par, sd = 0.3)
      key <- sprintf("mxl_gen_long|scr=%d|%s", scr, cell)
      res <- gen_est_outputs(theta, dat, W, rc_dist, rc_corr, rc_mean,
                             use_asc, cc$ioo, gen_seed = 31L, scr = scr,
                             S = S, K_w = K_w, Ti = cc$Ti, hess = FALSE)
      if (is.null(cc$Ti)) {
        res <- c(res, gen_pred_outputs(
          theta, dat, W, rc_dist, rc_corr, rc_mean, use_asc, cc$ioo,
          gen_seed = 31L, scr = scr, S = S, K_x = K_x, K_w = K_w,
          L_size = L_size, J = J, with_blp = FALSE)[c("pred", "shares",
                                                      "logsum")])
      }
      out[[key]] <- res
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Run, save or compare
# ---------------------------------------------------------------------------

# The solver tier's recorder (golden_blp()), used by the batteries below.
source(file.path(.script_dir, "blp_batteries.R"))

# Original batteries first, in their original order; S0 batteries last.
cat("Running MNL battery...\n");   res_mnl <- run_mnl()
cat("Running MXL battery...\n");   res_mxl <- run_mxl()
cat("Running NL battery...\n");    res_nl <- run_nl()
cat("Running error battery...\n"); res_err <- run_errors()
cat("Running MXL generate-mode battery (S0)...\n"); res_mxl_gen <- run_mxl_gen()
cat("Running MXL weights battery (S0)...\n");       res_mxl_wt <- run_mxl_wt()
cat("Running MXL panel battery (P1)...\n");         res_mxl_panel <- run_mxl_panel()
cat("Running MXL extreme-utility battery (P1)...\n"); res_mxl_edge <- run_mxl_edge()
cat("Running MXL generate-mode prediction battery (H1)...\n"); res_mxl_gen_pred <- run_mxl_gen_pred()
cat("Running MXL generate-mode K_w = 13 battery (H1)...\n");   res_mxl_gen_wide <- run_mxl_gen_wide()
cat("Running MXL generate-mode S = 1000 battery (H1)...\n");   res_mxl_gen_long <- run_mxl_gen_long()
source(file.path(.script_dir, "prep_batteries.R"))
cat("Running data-preparation battery (6A)...\n");  res_prep <- run_prep()
cat("Running store-mode cube battery (6A)...\n");   res_cube <- run_cube()
source(file.path(.script_dir, "p2_batteries.R"))
cat("Running MNL extras battery (P2)...\n");        res_mnl_x <- run_mnl_x()
cat("Running NL extras battery (P2)...\n");         res_nl_x <- run_nl_x()
cat("Running double-index battery (P2)...\n");      res_idx_dbl <- run_idx_dbl()
cat("Running MNL/NL error battery (P2)...\n");      res_err_p2 <- run_err_p2()
cat("Running HB/MNP battery (P2)...\n");            res_hb <- run_hb()
cat("Running fits battery (P2)...\n");              res_fits <- run_fits()
source(file.path(.script_dir, "p3_batteries.R"))
cat("Running MXL prediction shapes battery (P3)...\n");  res_pred_x <- run_pred_x()
cat("Running MXL panel post-estimation battery (P3)...\n"); res_pred_panel <- run_pred_panel()
cat("Running MXL double-index battery (P3)...\n");      res_pred_dbl <- run_pred_dbl()
cat("Running MXL prediction error battery (P3)...\n");  res_err_p3 <- run_err_p3()
cat("Running MXL chunked-draws battery (P3)...\n");     res_pred_chunk <- run_pred_chunk()
source(file.path(.script_dir, "p4_batteries.R"))
cat("Running MXL many-alternative derivative battery (P4)...\n"); res_mxl_asc <- run_mxl_asc()
cat("Running MXL meat battery (P4)...\n");               res_mxl_meat <- run_mxl_meat()
cat("Running MXL variance-route fits battery (P4)...\n"); res_mxl_vcov <- run_mxl_vcov()
cat("Running MXL non-finite battery (P4)...\n");         res_mxl_edge_asc <- run_mxl_edge_asc()
source(file.path(.script_dir, "p6b_batteries.R"))
cat("Running scale-on-load kernel battery (P6B)...\n"); res_scale_kern <- run_scale_kern()
cat("Running scaled MXL fits battery (P6B)...\n");      res_scale_fit <- run_scale_fit()
cat("Running column-scales battery (P6B)...\n");        res_scale_cols <- run_scale_cols()
cat("Running scaled MNL fits battery (P6B)...\n");      res_scale_mnl <- run_scale_mnl()
source(file.path(.script_dir, "p6d_batteries.R"))
cat("Running OPG-diagonal battery (P6D)...\n");         res_p6d_diag <- run_p6d_diag()
cat("Running coordinate-map battery (P6D)...\n");       res_p6d_coord <- run_p6d_coord()
cat("Running variance-route fits battery (P6D)...\n");  res_p6d_vcov <- run_p6d_vcov()
cat("Running new-option fits battery (P6D)...\n");      res_p6d_new <- run_p6d_new()
source(file.path(.script_dir, "start_batteries.R"))
cat("Running default-start battery (6D-f4)...\n");     res_start_unit <- run_start_unit()
source(file.path(.script_dir, "optim_batteries.R"))
cat("Running optim adapter battery (6D-f2)...\n");      res_optim_ro <- run_optim_ro()
cat("Running optim fits battery (6D-f2)...\n");         res_optim_fit <- run_optim_fit()
source(file.path(.script_dir, "uasc_batteries.R"))
cat("Running use_asc battery (6D-f6)...\n");            res_uasc <- run_uasc()
source(file.path(.script_dir, "shares_batteries.R"))
cat("Running predicted-shares battery...\n");          res_pred_shares <- run_pred_shares()
source(file.path(.script_dir, "p5a_batteries.R"))
cat("Running unbounded-constants battery (5A)...\n");   res_ucon <- run_ucon()
cat("Running polish battery (5A)...\n");                res_polish <- run_polish_bat()
source(file.path(.script_dir, "p5b_batteries.R"))
cat("Running constants' block battery (5B)...\n");     res_pdelta <- run_pdelta()
cat("Running profile battery (5B)...\n");              res_profile <- run_profile()
source(file.path(.script_dir, "p8_batteries.R"))
cat("Running MNL derivative battery (P8)...\n");      res_p8_mnl <- run_p8_mnl()
cat("Running NL derivative battery (P8)...\n");       res_p8_nl <- run_p8_nl()
cat("Running MNL/NL non-finite battery (P8)...\n");   res_p8_edge <- run_p8_edge()
cat("Running MNL/NL fits battery (P8)...\n");         res_p8_fit <- run_p8_fit()
source(file.path(.script_dir, "p8b_batteries.R"))
cat("Running MNL kernels battery (P8B)...\n");        res_p8b_mnl <- run_p8b_mnl()
cat("Running NL kernels battery (P8B)...\n");         res_p8b_nl <- run_p8b_nl()
cat("Running MNL/NL edge battery (P8B)...\n");        res_p8b_edge <- run_p8b_edge()
cat("Running MXL elasticity battery (P8B)...\n");     res_p8b_mxl <- run_p8b_mxl()
cat("Running cluster-label battery (P8B)...\n");      res_p8b_cl <- run_p8b_cl()
source(file.path(.script_dir, "p8c_batteries.R"))
cat("Running gradient overflow battery (P8C)...\n");  res_p8c_grad <- run_p8c_grad()
cat("Running MNL/NL polish battery (P8C)...\n");      res_p8c_polish <- run_p8c_polish()
cat("Running MXL polish battery (P8C)...\n");         res_p8c_mxl <- run_p8c_mxl()
cat("Running share-inversion battery (BLP A1)...\n");  res_blp <- run_blp()
oc_save_calls()

results <- list(mnl = res_mnl, mxl = res_mxl, nl = res_nl, errors = res_err,
                mxl_gen = res_mxl_gen, mxl_wt = res_mxl_wt,
                mxl_panel = res_mxl_panel, mxl_edge = res_mxl_edge,
                mxl_gen_pred = res_mxl_gen_pred,
                mxl_gen_wide = res_mxl_gen_wide,
                mxl_gen_long = res_mxl_gen_long,
                prep = res_prep, cube = res_cube,
                mnl_x = res_mnl_x, nl_x = res_nl_x, idx_dbl = res_idx_dbl,
                err_p2 = res_err_p2, hb = res_hb, fits = res_fits,
                pred_x = res_pred_x, pred_panel = res_pred_panel,
                pred_dbl = res_pred_dbl, err_p3 = res_err_p3,
                pred_chunk = res_pred_chunk,
                mxl_asc = res_mxl_asc, mxl_meat = res_mxl_meat,
                mxl_vcov = res_mxl_vcov, mxl_edge_asc = res_mxl_edge_asc,
                scale_kern = res_scale_kern, scale_fit = res_scale_fit,
                scale_cols = res_scale_cols, scale_mnl = res_scale_mnl,
                p6d_diag = res_p6d_diag, p6d_coord = res_p6d_coord,
                p6d_vcov = res_p6d_vcov, p6d_new = res_p6d_new,
                start_unit = res_start_unit,
                optim_ro = res_optim_ro, optim_fit = res_optim_fit,
                uasc = res_uasc, pred_shares = res_pred_shares,
                ucon = res_ucon, polish = res_polish,
                pdelta = res_pdelta, profile = res_profile,
                p8_mnl = res_p8_mnl, p8_nl = res_p8_nl,
                p8_edge = res_p8_edge, p8_fit = res_p8_fit,
                p8b_mnl = res_p8b_mnl, p8b_nl = res_p8b_nl,
                p8b_edge = res_p8b_edge, p8b_mxl = res_p8b_mxl,
                p8b_cl = res_p8b_cl,
                p8c_grad = res_p8c_grad, p8c_polish = res_p8c_polish,
                p8c_mxl = res_p8c_mxl, blp = res_blp)

flatten <- function(x, prefix = "") {
  if (is.list(x) && !isTRUE(attr(x, "golden_leaf"))) {
    nms <- names(x)
    if (is.null(nms)) nms <- as.character(seq_along(x))
    out <- list()
    for (i in seq_along(x)) {
      out <- c(out, flatten(x[[i]], paste0(prefix, "$", nms[i])))
    }
    out
  } else {
    stats::setNames(list(x), prefix)
  }
}

# Tolerance group = outputs of the four MXL estimation kernels. Everything
# else (MNL, NL, MXL prediction family, Cholesky utilities, error messages)
# must stay identical().
TOL_REL <- 1e-10
TOL_BATTERIES <- c("mxl", "mxl_gen", "mxl_wt", "mxl_panel", "mxl_edge",
                   "mxl_gen_wide", "mxl_gen_long",
                   "mxl_asc", "mxl_meat", "mxl_vcov", "mxl_edge_asc")
TOL_SLOTS <- c("grad", "hess", "bhhh", "scores", "tastes")
# --tol-slots=a,b,c replaces the tolerance slots (P4: hess,bhhh,meat,bread,
# vcov,se,posthoc, so that grad, scores and tastes must stay identical).
.tol_arg <- sub("^--tol-slots=", "", grep("^--tol-slots=", args, value = TRUE))
if (length(.tol_arg)) TOL_SLOTS <- strsplit(.tol_arg[1], ",", fixed = TRUE)[[1]]
cat("Tolerance slots:", paste(TOL_SLOTS, collapse = ", "), "\n")
GROUP_ORDER <- c("mnl", "mxl:predict-family", "mxl:chol",
                 paste0("mxl:", TOL_SLOTS), paste0("mxl_gen:", TOL_SLOTS),
                 paste0("mxl_wt:", TOL_SLOTS), paste0("mxl_panel:", TOL_SLOTS),
                 paste0("mxl_edge:", TOL_SLOTS), "nl", "errors",
                 "mxl_gen_pred", "mxl_gen_wide",
                 paste0("mxl_gen_wide:", TOL_SLOTS), "mxl_gen_long",
                 paste0("mxl_gen_long:", TOL_SLOTS))

# Flattened names look like "$<battery>$<config>$<slot>[$<sub>]".
entry_info <- function(keys) {
  parts <- strsplit(keys, "$", fixed = TRUE)
  battery <- vapply(parts, function(p) p[2], "")
  config <- vapply(parts, function(p) if (length(p) >= 3) p[3] else NA_character_, "")
  slot <- vapply(parts, function(p) if (length(p) >= 4) p[4] else NA_character_, "")
  tol <- battery %in% TOL_BATTERIES & !is.na(slot) & slot %in% TOL_SLOTS
  group <- ifelse(tol, paste0(battery, ":", slot),
                  ifelse(battery == "mxl" & startsWith(config, "mxl_chol"),
                         "mxl:chol",
                         ifelse(battery == "mxl", "mxl:predict-family",
                                battery)))
  data.frame(entry = keys, battery = battery, slot = slot, group = group,
             tol = tol, stringsAsFactors = FALSE)
}

# max|a - b| / max(1, max|b|), b = baseline. Shape changes or a different
# non-finite pattern count as an infinite deviation.
rel_dev <- function(a, b) {
  if (!is.numeric(a) || !is.numeric(b)) return(Inf)
  if (length(a) != length(b) || !identical(dim(a), dim(b))) return(Inf)
  if (!length(b)) return(0)
  fa <- is.finite(a); fb <- is.finite(b)
  if (!identical(fa, fb)) return(Inf)
  if (!identical(a[!fb], b[!fb])) return(Inf)
  if (!any(fb)) return(0)
  max(abs(a[fb] - b[fb])) / max(1, max(abs(b[fb])))
}
abs_dev <- function(a, b) {
  if (!is.numeric(a) || !is.numeric(b) || length(a) != length(b)) {
    return(NA_real_)
  }
  d <- abs(as.numeric(a) - as.numeric(b))
  if (all(is.na(d))) NA_real_ else max(d, na.rm = TRUE)
}
fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))

meta <- list(
  lib = .lib, pkg_path = .pkg_path,
  version = as.character(utils::packageVersion("choicer", lib.loc = .lib)),
  built = utils::packageDescription("choicer", lib.loc = .lib)$Built,
  grad_formals = names(formals(mxl_loglik_gradient_parallel)),
  has_tastes_kernel = exists("mxl_conditional_tastes_parallel",
                             envir = asNamespace("choicer"), inherits = FALSE),
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  R = R.version.string, blas = unname(extSoftVersion()["BLAS"]))

new_flat <- flatten(results)
new_info <- entry_info(names(new_flat))

if (mode == "baseline") {
  if (file.exists(baseline_path) && !force_overwrite) {
    stop("Baseline already exists at ", baseline_path,
         "; pass --force to overwrite it.")
  }
  attr(results, "meta") <- meta
  saveRDS(results, baseline_path)
  cat(sprintf("Baseline saved to %s\n", baseline_path))
  cat(sprintf(paste0("  configs: %d MNL, %d MXL (incl. %d Cholesky), %d MXL-gen, ",
                     "%d MXL-wt, %d NL; %d error cases\n"),
              length(res_mnl), length(res_mxl),
              sum(startsWith(names(res_mxl), "mxl_chol")), length(res_mxl_gen),
              length(res_mxl_wt), length(res_nl), length(res_err)))
  n_solver <- sum(vapply(new_flat, inherits, TRUE, what = "golden_blp"))
  cat(sprintf(paste0("  leaf entries: %d = %d identical-group + %d tolerance-group ",
                     "+ %d solver-tier (BLP)\n"),
              nrow(new_info), sum(!new_info$tol) - n_solver, sum(new_info$tol),
              n_solver))
  # a reference a call site gives must reproduce its targets
  noted <- names(new_flat)[vapply(new_flat, function(x) {
    inherits(x, "golden_blp") && !is.na(x$ref_note)
  }, TRUE)]
  if (length(noted)) {
    cat(sprintf("  WARNING: %d solver-tier leaves whose given reference misses its targets:\n",
                length(noted)))
    cat(sprintf("    %s\n", head(noted, 20)), sep = "")
  }
  tab <- table(factor(new_info$group, levels = GROUP_ORDER))
  print(tab[tab > 0])
  cat("  md5:", unname(tools::md5sum(baseline_path)), "\n")
  quit(status = 0)
}

# mode == "check"
if (!file.exists(baseline_path)) {
  stop("No baseline found at ", baseline_path,
       " — run with 'baseline' on the pre-refactor build first.")
}
baseline <- readRDS(baseline_path)
bmeta <- attr(baseline, "meta")
if (!is.null(bmeta)) {
  cat(sprintf("Baseline: choicer %s from %s (built %s; created %s)\n",
              bmeta$version, bmeta$lib, bmeta$built, bmeta$created))
}
cat(sprintf("Current : choicer %s from %s (built %s); Ti arg: %s; tastes kernel: %s\n",
            meta$version, meta$lib, meta$built,
            "Ti" %in% meta$grad_formals, meta$has_tastes_kernel))
cat("Baseline md5:", unname(tools::md5sum(baseline_path)), "\n\n")

base_flat <- flatten(baseline)

fail <- FALSE
only_base <- setdiff(names(base_flat), names(new_flat))
only_new <- setdiff(names(new_flat), names(base_flat))
if (length(only_base) || length(only_new)) {
  fail <- TRUE
  show <- function(x) {
    if (length(x) > 20) c(head(x, 20), sprintf("... (%d more)", length(x) - 20)) else x
  }
  cat("MISMATCH in output structure:\n")
  cat(sprintf("  only in baseline (%d):\n", length(only_base)))
  if (length(only_base)) cat(paste0("    ", show(only_base), "\n"), sep = "")
  cat(sprintf("  only in new (%d):\n", length(only_new)))
  if (length(only_new)) cat(paste0("    ", show(only_new), "\n"), sep = "")
}

keys <- intersect(names(base_flat), names(new_flat))
rep_df <- entry_info(keys)
# The solver tier: BLP leaves (golden_blp), each battery's in a group of its own
rep_df$solver <- vapply(keys, function(k) {
  inherits(base_flat[[k]], "golden_blp") || inherits(new_flat[[k]], "golden_blp")
}, TRUE, USE.NAMES = FALSE)
rep_df$tol[rep_df$solver] <- FALSE
rep_df$group[rep_df$solver] <- paste0(rep_df$battery[rep_df$solver], ":blp")
rep_df$criterion <- ifelse(rep_df$solver, "solver",
                           ifelse(rep_df$tol, sprintf("rel<=%.0e", TOL_REL),
                                  if (bitwise) "bitwise" else "identical"))
rep_df$identical <- NA
rep_df$rel_dev <- NA_real_
rep_df$abs_dev <- NA_real_
rep_df$pass <- NA
rep_df$note <- NA_character_
solver_res <- list()
for (i in seq_along(keys)) {
  b <- base_flat[[keys[i]]]
  n <- new_flat[[keys[i]]]
  ident <- if (bitwise) {
    identical(gb_core(b), gb_core(n), num.eq = FALSE, single.NA = FALSE)
  } else {
    identical(gb_core(b), gb_core(n))
  }
  rep_df$identical[i] <- ident
  if (rep_df$solver[i]) {
    tr <- if (ident) {
      list(pass = TRUE, why = "identical", res_ratio = NA_real_,
           bound_ratio = NA_real_, vs_main = NA_real_, vs_main_res = NA_real_)
    } else {
      tryCatch(gb_tier(b, n), error = function(e) {
        list(pass = FALSE, why = paste("the tier failed:", conditionMessage(e)),
             res_ratio = NA_real_, bound_ratio = NA_real_, vs_main = NA_real_,
             vs_main_res = NA_real_)
      })
    }
    tr$identical <- ident
    solver_res[[keys[i]]] <- tr
    rep_df$pass[i] <- tr$pass
    rep_df$note[i] <- paste(tr$why, collapse = "; ")
    next
  }
  rd <- if (is.numeric(b)) (if (ident) 0 else rel_dev(n, b)) else NA_real_
  rep_df$rel_dev[i] <- rd
  rep_df$abs_dev[i] <- if (is.numeric(b)) (if (ident) 0 else abs_dev(n, b)) else NA_real_
  rep_df$pass[i] <- if (rep_df$tol[i]) ident || (is.finite(rd) && rd <= TOL_REL) else ident
}
utils::write.csv(rep_df, report_path, row.names = FALSE)

# ---- Per-group summary ----
grp_levels <- c(intersect(GROUP_ORDER, rep_df$group),
                setdiff(unique(rep_df$group), GROUP_ORDER))
summ <- do.call(rbind, lapply(grp_levels, function(g) {
  x <- rep_df[rep_df$group == g, , drop = FALSE]
  rd <- replace(x$rel_dev, is.na(x$rel_dev), -Inf)
  w <- which.max(rd)
  data.frame(group = g, criterion = x$criterion[1], n = nrow(x),
             bit_identical = sum(x$identical), fail = sum(!x$pass),
             worst_rel_dev = if (length(w) && rd[w] > -Inf) fmt(rd[w]) else "NA",
             stringsAsFactors = FALSE)
}))
cat("Per-group summary (worst_rel_dev = max|new - base| / max(1, max|base|)):\n")
print(summ, row.names = FALSE, right = FALSE)

if (verbose) {
  cat("\nAll entries:\n")
  cat(sprintf("  %-4s %-10s rel=%-9s %s\n", ifelse(rep_df$pass, "ok", "FAIL"),
              rep_df$criterion, fmt(rep_df$rel_dev), rep_df$entry), sep = "")
}

bad <- rep_df[!rep_df$pass, , drop = FALSE]
if (nrow(bad)) {
  fail <- TRUE
  show_n <- min(nrow(bad), 60L)
  cat(sprintf("\nFAILING ENTRIES (%d; full per-entry table in %s):\n",
              nrow(bad), report_path))
  for (i in seq_len(show_n)) {
    cat(sprintf("  %-10s %s  rel=%s abs=%s\n", bad$criterion[i], bad$entry[i],
                fmt(bad$rel_dev[i]), fmt(bad$abs_dev[i])))
    b <- base_flat[[bad$entry[i]]]
    n <- new_flat[[bad$entry[i]]]
    if (isTRUE(bad$solver[i])) {
      cat("      solver tier:", bad$note[i], "\n")
      ne <- gb_field(n, "error")
      if (length(ne) == 1L && !is.na(ne)) cat("      new error:", ne, "\n")
    } else if (is.character(b) || is.character(n)) {
      cat("      base: ", paste(format(b), collapse = " | "), "\n",
          "      new:  ", paste(format(n), collapse = " | "), "\n", sep = "")
    } else if (!(is.numeric(b) && is.numeric(n) && length(b) == length(n))) {
      cat("      all.equal:", paste(all.equal(b, n), collapse = "; "), "\n")
    }
  }
  if (nrow(bad) > show_n) cat(sprintf("  ... %d more\n", nrow(bad) - show_n))
}

tl <- rep_df[rep_df$tol, , drop = FALSE]
tl_rd <- replace(tl$rel_dev, is.na(tl$rel_dev), Inf)
worst_i <- which.max(tl_rd)
worst_txt <- if (length(worst_i)) {
  sprintf("%s (%s)", fmt(tl_rd[worst_i]), tl$entry[worst_i])
} else "NA"
cat(sprintf(paste0("\nTolerance group (MXL estimation kernels): %d entries, ",
                   "%d bit-identical; worst rel dev = %s\n"),
            nrow(tl), sum(tl$identical), worst_txt))
nz <- tl[!tl$identical, , drop = FALSE]
if (nrow(nz)) {
  nz <- head(nz[order(-replace(nz$rel_dev, is.na(nz$rel_dev), Inf)), ], 10)
  cat("  largest deviations:\n")
  cat(sprintf("    rel=%s  %s\n", fmt(nz$rel_dev), nz$entry), sep = "")
}

sv <- rep_df[rep_df$solver, , drop = FALSE]
if (nrow(sv)) {
  cat(sprintf(paste0("\nSolver tier (BLP leaves): %d entries, %d identical to the ",
                     "baseline, %d failing\n"),
              nrow(sv), sum(sv$identical), sum(!sv$pass)))
  rpt <- gb_tier_summary(sv$entry, base_flat, new_flat, solver_res[sv$entry])
  cat(sprintf(paste0("  converged: baseline %d, new %d; errors: baseline %d, new %d; ",
                     "no reference: %d; ill-conditioned (2 kappa tol > 1e-4): %d; ",
                     "changed and unverified (no residual recomputable): %d\n"),
              sum(rpt$main_conv), sum(rpt$new_conv), sum(rpt$main_error),
              sum(rpt$new_error), sum(!rpt$has_ref), sum(rpt$illcond),
              sum(rpt$unverified)))
  worst <- function(x) if (any(is.finite(x))) fmt(max(x[is.finite(x)])) else "NA"
  cat(sprintf(paste0("  worst residual / tol: %s; worst error / bound: %s ",
                     "(changed leaves that converged)\n"),
              worst(rpt$res_ratio), worst(rpt$bound_ratio)))
  less <- rpt[is.finite(rpt$vs_main) & rpt$vs_main > 1, , drop = FALSE]
  cat(sprintf(paste0("  less accurate than the baseline (reported, not gated): %d ",
                     "entries%s\n"),
              nrow(less), if (nrow(less)) sprintf(", worst x%s", fmt(max(less$vs_main))) else ""))
  if (nrow(less)) {
    less <- head(less[order(-less$vs_main), ], 10)
    cat(sprintf("    x%s  %s\n", fmt(less$vs_main), less$entry), sep = "")
  }
  worse <- rpt[is.finite(rpt$vs_main_res) & rpt$vs_main_res > 1, , drop = FALSE]
  cat(sprintf(paste0("  unconverged on both sides with a larger residual than the ",
                     "baseline's (reported): %d entries%s\n"),
              nrow(worse), if (nrow(worse)) sprintf(", worst x%s", fmt(max(worse$vs_main_res))) else ""))
  if (nrow(worse)) {
    worse <- head(worse[order(-worse$vs_main_res), ], 10)
    cat(sprintf("    x%s  %s\n", fmt(worse$vs_main_res), worse$entry), sep = "")
  }
  base_miss <- vapply(sv$entry, function(k) {
    b <- base_flat[[k]]
    isTRUE(gb_field(b, "converged")) && is.finite(b$residual) && is.finite(b$tol) &&
      b$residual > b$tol
  }, TRUE)
  cat(sprintf(paste0("  baseline leaves that reported convergence with a residual ",
                     "above their tol: %d\n"), sum(base_miss)))
  if (any(rpt$ref_note)) {
    cat(sprintf("  WARNING: %d leaves whose given reference misses its targets\n",
                sum(rpt$ref_note)))
  }
  mo <- rpt$entry[rpt$main_only]
  cat(sprintf(paste0("  converged by the baseline's own flag, not by the new build's ",
                     "(reported): %d entries\n"), length(mo)))
  if (length(mo)) cat(sprintf("    %s\n", utils::head(mo, 10)), sep = "")
  uv <- rpt$entry[rpt$unverified]
  if (length(uv)) {
    cat(sprintf("  WARNING: %d changed leaves whose residual cannot be recomputed:\n",
                length(uv)))
    cat(sprintf("    %s\n", utils::head(uv, 10)), sep = "")
  }
  if (verbose) print(rpt, row.names = FALSE)
}

n_idg <- sum(!rep_df$tol & !rep_df$solver)
if (fail) {
  if (nrow(bad) && all(bad$battery == "mxl_edge") && !length(only_base) && !length(only_new)) {
    cat("\nAll failing entries are in the mxl_edge (extreme-utility) battery.\n")
  }
  cat(sprintf(paste0("\nFAIL: %d of %d compared entries failed%s ",
                     "(identical-group failures: %d; tolerance-group failures: %d; ",
                     "solver-tier failures: %d). Worst tolerance-group rel dev: %s\n"),
              nrow(bad), nrow(rep_df),
              if (length(only_base) || length(only_new)) " + output-structure mismatch" else "",
              sum(!bad$tol & !bad$solver), sum(bad$tol), sum(bad$solver), worst_txt))
  quit(status = 1)
}
cat(sprintf(paste0("\nPASS: all %d entries OK — %d identical-group entries ",
                   "bit-identical; %d tolerance-group entries within rel %.0e ",
                   "(%d bit-identical); %d solver-tier entries (%d identical). ",
                   "Worst tolerance-group rel dev: %s\n"),
            nrow(rep_df), n_idg, nrow(tl), TOL_REL, sum(tl$identical),
            sum(rep_df$solver), sum(rep_df$solver & rep_df$identical), worst_txt))
quit(status = 0)
