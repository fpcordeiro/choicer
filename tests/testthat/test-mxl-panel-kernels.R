# Kernel tests for the panel mixed logit likelihood (Revelt & Train 1998).
#
# The four MXL estimation kernels and mxl_conditional_tastes_parallel() take a
# decision-maker structure `Ti` (the number of consecutive choice situations of
# each decision maker). Every check compares the kernels with a computation
# that shares no code with the C++:
#   * the brute-force oracle in helper-mxl-panel.R, derived from the estimand
#     (objective, per-decision-maker log-likelihood, conditional tastes);
#   * numerical derivatives of the kernel objective and of the oracle;
#   * exact identities: BHHH and the gradient as weighted sums of the unit
#     score rows, Hessian symmetry and linearity in the weights, Ti = rep(1, N)
#     as the cross-section, generate mode as store mode, and the collapse of
#     the panel likelihood to a cross-section on shared draws when S = 1.
# The cells of mxlp_cell_configs() are small unbalanced panels (T_u in 1..4,
# 1 to 4 inside alternatives per situation) that jointly cover every model
# option; each failure message names its cell.

mxlp_cells <- mxlp_build_cells()

# --- Objective, gradient and Hessian -----------------------------------------

test_that("panel objective equals the brute-force oracle", {
  for (fx in mxlp_cells) {
    res <- mxlp_call("gradient", fx)
    orc <- mxlp_oracle(fx$theta, fx)
    expect_length(res$gradient, length(fx$theta))
    mxlp_expect_close(res$objective, -orc$loglik, 1e-10,
                      sprintf("[%s] objective vs oracle", fx$name))
  }
})

test_that("panel gradient matches numDeriv of kernel objective and oracle", {
  skip_if_not_installed("numDeriv")
  for (fx in mxlp_cells) {
    grad <- mxlp_call("gradient", fx)$gradient
    obj_kernel <- function(th) mxlp_call("gradient", fx, theta = th)$objective
    obj_oracle <- function(th) -mxlp_oracle(th, fx)$loglik
    mxlp_expect_close(grad, numDeriv::grad(obj_kernel, fx$theta), TOL_GRAD,
                      sprintf("[%s] gradient vs numDeriv(kernel objective)",
                              fx$name))
    mxlp_expect_close(grad, numDeriv::grad(obj_oracle, fx$theta), TOL_GRAD,
                      sprintf("[%s] gradient vs numDeriv(oracle)", fx$name))
  }
})

test_that("panel Hessian matches numDeriv of the kernel objective", {
  skip_if_not_installed("numDeriv")
  for (fx in mxlp_cells) {
    H <- mxlp_call("hessian", fx)
    p <- length(fx$theta)
    expect_equal(dim(H), c(p, p), label = sprintf("[%s] dim(H)", fx$name))
    obj <- function(th) mxlp_call("gradient", fx, theta = th)$objective
    mxlp_expect_close(H, numDeriv::hessian(obj, fx$theta), TOL_HESS,
                      sprintf("[%s] Hessian vs numDeriv", fx$name))
  }
})

test_that("panel Hessian is symmetric and linear in the weights", {
  for (fx in mxlp_cells) {
    H <- mxlp_call("hessian", fx)
    mxlp_expect_close(H, t(H), 1e-10,
                      sprintf("[%s] Hessian symmetry", fx$name))
    H2 <- mxlp_call("hessian", fx, weights = 2 * fx$weights)
    mxlp_expect_close(H2, 2 * H, 1e-10,
                      sprintf("[%s] H(2w) vs 2 H(w)", fx$name))
  }
})

# --- Unit scores, BHHH and conditional tastes --------------------------------

test_that("BHHH and the gradient are weighted sums over the unit score rows", {
  for (fx in mxlp_cells) {
    S_u <- mxlp_call("scores", fx)
    w_u <- mxlp_unit_weights(fx)
    expect_equal(dim(S_u), c(fx$U, length(fx$theta)),
                 label = sprintf("[%s] dim(scores)", fx$name))
    mxlp_expect_close(mxlp_call("bhhh", fx), crossprod(sqrt(w_u) * S_u),
                      1e-10, sprintf("[%s] BHHH vs sum_u w_u s_u s_u'",
                                     fx$name))
    # The kernel returns the gradient of the negated log-likelihood.
    mxlp_expect_close(colSums(w_u * S_u), -mxlp_call("gradient", fx)$gradient,
                      1e-10, sprintf("[%s] sum_u w_u s_u vs -gradient",
                                     fx$name))
  }
})

test_that("score rows match numDeriv of each decision maker's log-likelihood", {
  skip_if_not_installed("numDeriv")
  for (fx in mxlp_cells) {
    S_u <- mxlp_call("scores", fx)
    unit_of <- rep(seq_len(fx$U), fx$Ti)
    for (u in fx$probe_units) {
      what <- sprintf("[%s] score row %d (T_u = %d)", fx$name, u, fx$Ti[u])
      g_oracle <- numDeriv::grad(mxlp_oracle_unit_loglik, fx$theta,
                                 fx = fx, u = u)
      mxlp_expect_close(S_u[u, ], g_oracle, 1e-6,
                        paste(what, "vs numDeriv(oracle l_u)"))
      # Indicator weights, constant within every unit, select decision
      # maker u in the kernel's own objective.
      w_ind <- as.numeric(unit_of == u)
      loglik_u <- function(th) {
        -mxlp_call("gradient", fx, theta = th, weights = w_ind)$objective
      }
      mxlp_expect_close(S_u[u, ], numDeriv::grad(loglik_u, fx$theta), 1e-6,
                        paste(what, "vs numDeriv(indicator-weighted kernel)"))
    }
  }
})

test_that("conditional tastes equal the oracle's posterior draw moments", {
  for (fx in mxlp_cells) {
    ct <- mxlp_call("tastes", fx)
    orc <- mxlp_oracle(fx$theta, fx)
    expect_named(ct, c("mean", "sd"), ignore.order = TRUE)
    expect_equal(dim(ct$mean), c(fx$K_w, fx$U),
                 label = sprintf("[%s] dim(tastes$mean)", fx$name))
    expect_equal(dim(ct$sd), c(fx$K_w, fx$U),
                 label = sprintf("[%s] dim(tastes$sd)", fx$name))
    mxlp_expect_close(ct$mean, orc$mean, 1e-10,
                      sprintf("[%s] conditional mean vs oracle", fx$name))
    mxlp_expect_close(ct$sd, orc$sd, 1e-10,
                      sprintf("[%s] conditional SD vs oracle", fx$name))
  }
})

# --- Equivalences -----------------------------------------------------------

test_that("Ti = rep(1, N) reproduces the cross-sectional kernels", {
  for (fx in mxlp_cells) {
    xs <- mxlp_cross_section(fx)
    ones <- rep(1L, xs$N)
    for (k in mxlp_kernels) {
      res_null <- mxlp_call(k, xs, Ti = NULL)
      res_ones <- mxlp_call(k, xs, Ti = ones)
      mxlp_expect_close(res_ones, res_null, 1e-12,
                        sprintf("[%s] %s, Ti = rep(1, N) vs NULL",
                                fx$name, k))
    }
    orc <- mxlp_oracle(xs$theta, xs)
    mxlp_expect_close(mxlp_call("gradient", xs)$objective, -orc$loglik, 1e-10,
                      sprintf("[%s] cross-section objective vs oracle",
                              fx$name))
    ct <- mxlp_call("tastes", xs)
    expect_equal(dim(ct$mean), c(xs$K_w, xs$N),
                 label = sprintf("[%s] cross-section dim(tastes$mean)",
                                 fx$name))
    mxlp_expect_close(ct$mean, orc$mean, 1e-10,
                      sprintf("[%s] cross-section conditional mean vs oracle",
                              fx$name))
    mxlp_expect_close(ct$sd, orc$sd, 1e-10,
                      sprintf("[%s] cross-section conditional SD vs oracle",
                              fx$name))
  }
})

test_that("generate mode reproduces store-mode panel draws", {
  # gen_scramble = 0 generates plain Halton normals: decision maker u gets
  # the same block as slice u of get_halton_normals(S, U, K_w).
  for (fx in mxlp_cells) {
    for (k in mxlp_kernels) {
      mxlp_expect_close(mxlp_call(k, fx, generate = TRUE), mxlp_call(k, fx),
                        1e-12, sprintf("[%s] %s, generate vs store",
                                       fx$name, k))
    }
  }
})

test_that("with S = 1 the panel kernels are a cross-section on shared draws", {
  # With S = 1, l_u = sum_t log P_t(c_u), so the panel objective, gradient and
  # Hessian equal the cross-sectional ones when each situation uses its
  # decision maker's draw; unit scores are within-unit sums of situation
  # scores; and the conditional tastes are that draw, with zero SD.
  for (fx in mxlp_cells) {
    f1 <- fx
    f1$S <- 1L
    f1$eta <- get_halton_normals(1L, fx$U, fx$K_w)
    unit_of <- rep(seq_len(fx$U), fx$Ti)
    xs <- f1
    xs[["Ti"]] <- NULL
    xs$eta <- f1$eta[, , unit_of, drop = FALSE]

    pan <- mxlp_call("gradient", f1)
    crs <- mxlp_call("gradient", xs)
    mxlp_expect_close(pan$objective, crs$objective, 1e-10,
                      sprintf("[%s] S = 1 objective", fx$name))
    mxlp_expect_close(pan$gradient, crs$gradient, 1e-10,
                      sprintf("[%s] S = 1 gradient", fx$name))
    mxlp_expect_close(mxlp_call("hessian", f1), mxlp_call("hessian", xs),
                      1e-10, sprintf("[%s] S = 1 Hessian", fx$name))
    mxlp_expect_close(mxlp_call("scores", f1),
                      rowsum(mxlp_call("scores", xs), unit_of), 1e-10,
                      sprintf("[%s] S = 1 unit scores", fx$name))
    ct <- mxlp_call("tastes", f1)
    mxlp_expect_close(ct$mean, mxlp_oracle(f1$theta, f1)$mean, 1e-10,
                      sprintf("[%s] S = 1 conditional mean", fx$name))
    mxlp_expect_close(ct$sd, matrix(0, fx$K_w, fx$U), 1e-10,
                      sprintf("[%s] S = 1 conditional SD", fx$name))
  }
})

# --- Numerical stability -----------------------------------------------------

test_that("panel kernels stay accurate when a unit's likelihood underflows", {
  # Decision maker 2 makes 150 low-probability choices, so exp(lambda_us) is
  # exactly 0 in double precision for every draw. The panel likelihood and
  # the posterior draw weights must be computed in log space.
  skip_if_not_installed("numDeriv")
  fx <- mxlp_long_panel_fixture()
  orc <- mxlp_oracle(fx$theta, fx)
  expect_true(all(exp(orc$lambda[2L, ]) == 0))  # the fixture does underflow

  res <- mxlp_call("gradient", fx)
  mxlp_expect_close(res$objective, -orc$loglik, 1e-10,
                    "[long panel] objective vs oracle")
  obj <- function(th) mxlp_call("gradient", fx, theta = th)$objective
  mxlp_expect_close(res$gradient, numDeriv::grad(obj, fx$theta), TOL_GRAD,
                    "[long panel] gradient vs numDeriv")
  mxlp_expect_close(mxlp_call("hessian", fx), numDeriv::hessian(obj, fx$theta),
                    TOL_HESS, "[long panel] Hessian vs numDeriv")

  S_u <- mxlp_call("scores", fx)
  mxlp_expect_close(colSums(S_u), -res$gradient, 1e-10,
                    "[long panel] sum_u s_u vs -gradient")
  mxlp_expect_close(mxlp_call("bhhh", fx), crossprod(S_u), 1e-10,
                    "[long panel] BHHH vs sum_u s_u s_u'")
  ct <- mxlp_call("tastes", fx)
  mxlp_expect_close(ct$mean, orc$mean, 1e-10,
                    "[long panel] conditional mean vs oracle")
  mxlp_expect_close(ct$sd, orc$sd, 1e-10,
                    "[long panel] conditional SD vs oracle")
})

# Cells with choice sets of 2 to 24 alternatives, straddling Armadillo's
# 16-element local storage (the kernels keep grow-only utility and probability
# buffers across situations of different sizes).
mxlp_wide_cells <- list(
  mxlp_fixture("wide|n+l|corr|mean|W=row", seed = 71, J = 24L,
               rc_dist = c(0L, 1L), rc_correlation = TRUE, rc_mean = TRUE),
  mxlp_fixture("wide|n+n|diag|oo|W=alt", seed = 72, J = 24L,
               include_outside_option = TRUE, W_layout = "alt",
               weight_type = "person")
)

test_that("choice sets of more than 16 alternatives match the oracle", {
  skip_if_not_installed("numDeriv")
  for (fx in mxlp_wide_cells) {
    orc <- mxlp_oracle(fx$theta, fx)
    res <- mxlp_call("gradient", fx)
    mxlp_expect_close(res$objective, -orc$loglik, 1e-10,
                      sprintf("[%s] objective vs oracle", fx$name))
    obj <- function(th) mxlp_call("gradient", fx, theta = th)$objective
    mxlp_expect_close(res$gradient, numDeriv::grad(obj, fx$theta), TOL_GRAD,
                      sprintf("[%s] gradient vs numDeriv", fx$name))
    mxlp_expect_close(mxlp_call("hessian", fx),
                      numDeriv::hessian(obj, fx$theta), TOL_HESS,
                      sprintf("[%s] Hessian vs numDeriv", fx$name))
    ct <- mxlp_call("tastes", fx)
    mxlp_expect_close(ct$mean, orc$mean, 1e-10,
                      sprintf("[%s] conditional mean vs oracle", fx$name))
    mxlp_expect_close(ct$sd, orc$sd, 1e-10,
                      sprintf("[%s] conditional SD vs oracle", fx$name))
  }
})

test_that("draw batches reproduce the single-batch kernels", {
  # A decision maker whose R x S matrices would exceed the per-thread budget
  # has its draws processed a few at a time, the score folded in batch by
  # batch with a streaming log-sum-exp; the Hessian instead forms each
  # situation's draws in turn. Forcing batches of one and of at most five
  # draws must reproduce the one-batch results up to rounding, in store and
  # generate
  # mode, with wide choice sets, and on the long panel whose likelihood
  # underflows (where the streaming reference moves the most).
  cells <- c(mxlp_cells, list(mxlp_long_panel_fixture()), mxlp_wide_cells)
  for (fx in cells) {
    for (k in mxlp_kernels) {
      ref <- mxlp_call(k, fx)
      for (b in c(1L, 5L)) {
        mxlp_expect_close(mxlp_call(k, fx, draw_batch = b), ref, 1e-10,
                          sprintf("[%s] %s, batches of at most %d draws",
                                  fx$name, k, b))
      }
      mxlp_expect_close(mxlp_call(k, fx, generate = TRUE, draw_batch = 3L),
                        mxlp_call(k, fx, generate = TRUE), 1e-10,
                        sprintf(paste("[%s] %s, generate mode, batches of at",
                                      "most 3 draws"), fx$name, k))
    }
  }
})

test_that("the batched path is taken, and chosen automatically past the budget", {
  # One thread, so that a repeated call is bitwise reproducible and any
  # difference comes from the batches. At S = 601 the long panel's second
  # decision maker stacks 450 rows, and 450 x 601 exceeds the 2^18
  # row-draws that fit the score kernels' budget, so they split its draws on
  # their own, into batches of 300 and 301 draws; one forced batch of all S
  # draws is the reference.
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  fx <- mxlp_long_panel_fixture()
  expect_false(identical(mxlp_call("gradient", fx, draw_batch = 1L),
                         mxlp_call("gradient", fx)))
  fx <- mxlp_long_panel_fixture(S = 601L)
  for (k in mxlp_kernels) {
    mxlp_expect_close(mxlp_call(k, fx), mxlp_call(k, fx, draw_batch = fx$S),
                      1e-10, sprintf("[S = 601] %s, automatic batches", k))
  }
  expect_false(identical(mxlp_call("gradient", fx),
                         mxlp_call("gradient", fx, draw_batch = fx$S)))
})

test_that("draw batches with no positive weight are folded exactly", {
  # Unit u's first five draws give its choices probability exactly zero: its
  # one random coefficient enters only the chosen alternative of one
  # situation (W = 10 there, 0 elsewhere in the unit), and those draws are
  # -1e308, so the chosen utility overflows to -Inf. lambda_us is -Inf for
  # s <= 5 and finite after, which batches of at most 1, 2 and 5 draws meet
  # at every step of the streaming fold (all -Inf batches, a mixed batch, the
  # first batch with weight). Those draws carry weight exactly zero, so the
  # objective, scores and BHHH stay finite; the Hessian and the tastes, which
  # multiply the -Inf draws by their zero weight, are NaN either way.
  fx <- mxlp_fixture("neg-inf draws", seed = 91, rc_dist = 0L)
  u <- fx$probe_units[2L]
  off <- mxlp_unit_offsets(fx$Ti, fx$N)
  sits <- (off[u] + 1L):off[u + 1L]
  row0 <- c(0L, cumsum(fx$M))
  fx$W[(row0[sits[1L]] + 1L):row0[sits[length(sits)] + 1L], 1L] <- 0
  fx$W[row0[sits[1L]] + fx$choice_idx[sits[1L]], 1L] <- 10
  fx$eta[1L, 1:5, u] <- -1e308
  lam <- mxlp_oracle(fx$theta, fx)$lambda[u, ]
  expect_true(all(lam[1:5] == -Inf) && all(is.finite(lam[-(1:5)])),
              label = "only draws 1-5 give unit u's choices probability 0")
  ref <- lapply(c("gradient", "scores", "bhhh"), mxlp_call, fx = fx)
  expect_true(all(is.finite(unlist(ref))))
  for (b in c(1L, 2L, 5L)) {
    for (i in 1:3) {
      k <- c("gradient", "scores", "bhhh")[i]
      mxlp_expect_close(mxlp_call(k, fx, draw_batch = b), ref[[i]], 1e-10,
                        sprintf("%s, batches of at most %d draws", k, b))
    }
  }
})

test_that("overflowing utilities give the same sentinel, skip and NaN in draw batches", {
  for (panel in c(FALSE, TRUE)) {
    fx <- mxlp_edge_fixture(-36, panel)
    u <- fx$edge_unit
    t <- fx$edge_situation
    fx$theta[1L] <- 2
    fx$X[sum(fx$M[seq_len(t - 1L)]) + fx$choice_idx[t], 1L] <- -1e308
    ref_scores <- mxlp_call("scores", fx)
    ref_hess <- mxlp_call("hessian", fx)
    for (b in c(1L, 5L)) {
      what <- sprintf("[%s] batches of at most %d draws", fx$name, b)
      res <- mxlp_call("gradient", fx, draw_batch = b)
      expect_equal(res$objective, 1e10, label = paste(what, "sentinel"))
      expect_true(all(res$gradient == 0), label = paste(what, "zero gradient"))
      sc <- mxlp_call("scores", fx, draw_batch = b)
      expect_true(any(is.nan(sc[u, ])), label = paste(what, "NaN score row"))
      expect_identical(is.nan(sc[u, ]), is.nan(ref_scores[u, ]),
                       label = paste(what, "NaN pattern of the score row"))
      mxlp_expect_close(sc[-u, ], ref_scores[-u, ], 1e-10,
                        paste(what, "other score rows"))
      mxlp_expect_close(mxlp_call("hessian", fx, draw_batch = b), ref_hess,
                        1e-10, paste(what, "Hessian without the unit"))
      ct <- mxlp_call("tastes", fx, draw_batch = b)
      expect_true(all(is.na(ct$mean[, u])) && all(is.na(ct$sd[, u])),
                  label = paste(what, "NA tastes"))
    }
  }
})

test_that("the raw-array softmax and log-sum-exp match Armadillo's bit for bit", {
  # The draw loops use stable_softmax_n(), log_sum_exp_n() and
  # max_shifted_lse_n() on reused buffers; they must reproduce
  # stable_softmax(), logSumExp() and the log-sum's max + log(accu(exp()))
  # exactly, signed zeros, infinities and NaN included.
  set.seed(5)
  cases <- c(lapply(1:40, function(n) stats::rnorm(n, sd = 3)),
             list(0, -0, c(0, -0), c(-0, 0), c(1, 1), c(-Inf, 0), c(Inf, 1),
                  c(Inf, Inf), c(NaN, 1), c(1, NaN, 2), c(-Inf, -Inf),
                  c(1e300, -1e300), c(-800, -700, -745),
                  c(700, 710, 709.5, -1e308)))
  for (v in cases) {
    r <- choicer:::test_softmax_n(v)
    what <- paste0("v = c(", paste(format(v), collapse = ", "), ")")
    expect_true(identical(r$v_raw, r$v_arma, num.eq = FALSE),
                label = paste(what, ": shifted utilities"))
    expect_true(identical(r$p_raw, r$p_arma, num.eq = FALSE),
                label = paste(what, ": probabilities"))
    expect_true(identical(r$log_denom_raw, r$log_denom_arma, num.eq = FALSE),
                label = paste(what, ": log denominator"))
    expect_true(identical(r$lse_raw, r$lse_arma, num.eq = FALSE),
                label = paste(what, ": log-sum-exp"))
    expect_true(identical(r$lse_shift_raw, r$lse_shift_arma, num.eq = FALSE),
                label = paste(what, ": max-shifted log-sum"))
  }
})

# The edge fixtures put one chosen alternative `gap` utils below its
# competitors, once in a cross-section and once inside a multi-situation
# decision maker (see mxlp_edge_fixture()).

test_that("Hessian keeps a unit whose simulated likelihood is below 1e-12", {
  # A gap of -36 puts P(chosen) near 1e-16 at every draw: tiny but finite, so
  # the unit belongs in the Hessian. The pre-panel kernel dropped every unit
  # whose simulated probability was below 1e-12.
  skip_if_not_installed("numDeriv")
  for (panel in c(FALSE, TRUE)) {
    fx <- mxlp_edge_fixture(-36, panel)
    lik <- sum(exp(mxlp_oracle(fx$theta, fx)$lambda[fx$edge_unit, ]))
    expect_true(lik > 0 && lik < 1e-12,
                label = sprintf("[%s] summed unit likelihood in (0, 1e-12)",
                                fx$name))
    grad <- function(th) drop(mxlp_call("gradient", fx, theta = th)$gradient)
    mxlp_expect_close(mxlp_call("hessian", fx),
                      numDeriv::jacobian(grad, fx$theta), TOL_HESS,
                      sprintf("[%s] Hessian vs numDeriv::jacobian(gradient)",
                              fx$name))
  }
})

test_that("subnormal choice probabilities keep the objective and scores exact", {
  # P(chosen) is subnormal (positive, below .Machine$double.xmin) at every
  # draw. The kernels form log P(chosen) as the chosen utility minus the
  # log-sum-exp denominator, never as the log of the subnormal probability,
  # so the objective and the unit's score stay exact even at a gap of -740,
  # where P(chosen) itself keeps only about 6 significant bits. (Earlier
  # kernels took log(P(chosen)) and lost that precision; before the panel
  # work they also returned a zero gradient here.)
  skip_if_not_installed("numDeriv")
  for (gap in c(-720, -740)) {
    for (panel in c(FALSE, TRUE)) {
      fx <- mxlp_edge_fixture(gap, panel)
      u <- fx$edge_unit
      p <- exp(mxlp_situation_logp(fx$theta, fx, fx$edge_situation))
      expect_true(all(p > 0 & p < .Machine$double.xmin),
                  label = sprintf("[%s] P(chosen) subnormal at every draw",
                                  fx$name))
      res <- mxlp_call("gradient", fx)
      S_u <- mxlp_call("scores", fx)
      B <- mxlp_call("bhhh", fx)
      mxlp_expect_close(res$objective, -mxlp_oracle(fx$theta, fx)$loglik,
                        1e-10, sprintf("[%s] objective vs oracle", fx$name))
      expect_true(all(is.finite(S_u)) && all(is.finite(B)),
                  label = sprintf("[%s] finite scores and BHHH", fx$name))
      w_u <- mxlp_unit_weights(fx)
      mxlp_expect_close(colSums(w_u * S_u), -res$gradient, 1e-10,
                        sprintf("[%s] sum_u w_u s_u vs -gradient", fx$name))
      mxlp_expect_close(B, crossprod(sqrt(w_u) * S_u), 1e-10,
                        sprintf("[%s] BHHH vs sum_u w_u s_u s_u'", fx$name))
      g_oracle <- numDeriv::grad(mxlp_oracle_unit_loglik, fx$theta,
                                 fx = fx, u = u)
      mxlp_expect_close(S_u[u, ], g_oracle, 1e-6,
                        sprintf("[%s] score row %d vs numDeriv(oracle l_u)",
                                fx$name, u))
    }
  }
})

test_that("a choice probability that underflows to zero keeps an exact log-likelihood", {
  # Gaps of -1000 and -1e5 make P(chosen) exactly 0 in double precision at
  # every draw. In log space the choice keeps a finite, exact log
  # probability, so the unit's log-likelihood and its derivatives stay
  # exact: no optimizer sentinel, the Hessian keeps the unit, and its
  # conditional tastes are defined. With 10^7 choice situations some sit this
  # far in the tail, and one such unit used to turn the whole objective into
  # the sentinel. Checked for an inside choice with and without the outside
  # option, and for the outside option chosen against inside alternatives
  # far above it.
  skip_if_not_installed("numDeriv")
  cases <- expand.grid(layout = c("inside", "no_outside", "outside"),
                       panel = c(FALSE, TRUE), stringsAsFactors = FALSE)
  for (k in seq_len(nrow(cases))) {
    G <- H <- list()
    for (gap in c(-1000, -1e5)) {
      fx <- mxlp_edge_fixture(gap, cases$panel[k], cases$layout[k])
      p <- exp(mxlp_situation_logp(fx$theta, fx, fx$edge_situation))
      expect_true(all(p == 0),
                  label = sprintf("[%s] P(chosen) = 0 at every draw", fx$name))

      orc <- mxlp_oracle(fx$theta, fx)
      res <- mxlp_call("gradient", fx)
      mxlp_expect_close(res$objective, -orc$loglik, 1e-10,
                        sprintf("[%s] objective vs oracle", fx$name))
      mxlp_expect_close(colSums(mxlp_unit_weights(fx) *
                                  mxlp_call("scores", fx)),
                        -res$gradient, 1e-10,
                        sprintf("[%s] sum_u w_u s_u vs -gradient", fx$name))
      ct <- mxlp_call("tastes", fx)
      mxlp_expect_close(ct$mean, orc$mean, 1e-10,
                        sprintf("[%s] conditional mean vs oracle", fx$name))
      mxlp_expect_close(ct$sd, orc$sd, 1e-10,
                        sprintf("[%s] conditional SD vs oracle", fx$name))
      G[[length(G) + 1L]] <- res$gradient
      H[[length(H) + 1L]] <- mxlp_call("hessian", fx)
      if (gap == -1000) {
        obj <- function(th) -mxlp_oracle(th, fx)$loglik
        mxlp_expect_close(res$gradient, numDeriv::grad(obj, fx$theta),
                          TOL_GRAD,
                          sprintf("[%s] gradient vs numDeriv(oracle)", fx$name))
        grad <- function(th) drop(mxlp_call("gradient", fx, theta = th)$gradient)
        J <- numDeriv::jacobian(grad, fx$theta,
                                method.args = list(d = 1e-3, r = 6))
        mxlp_expect_close(H[[1L]], J, TOL_HESS,
                          sprintf("[%s] Hessian vs numDeriv::jacobian(gradient)",
                                  fx$name))
      }
    }
    # Once P(chosen) is 0 at every draw, the gap enters the unit's
    # log-likelihood only linearly, through beta_1 * x_1, so moving it from
    # -1000 to -1e5 shifts the gradient only along beta_1, by
    # w_u (1e5 - 1000) / beta_1, and leaves the Hessian unchanged. (At -1e5
    # the objective is about 1e5, too large for finite differences at these
    # tolerances. An uncentered Louis identity would miss the Hessian check by
    # about 0.15: the draw weights sum to one only to ulp(lambda), about
    # 1e-11, times g^2 = 4e10.)
    w_u <- mxlp_unit_weights(fx)[fx$edge_unit]
    shift <- replace(numeric(length(fx$theta)), 1L,
                     w_u * (1e5 - 1000) / fx$theta[1L])
    mxlp_expect_close(G[[2L]] - G[[1L]], shift, 1e-9,
                      sprintf("[%s] gradient shift from gap -1000 to -1e5",
                              fx$name))
    mxlp_expect_close(H[[2L]], H[[1L]], 1e-4,
                      sprintf("[%s] Hessian at gap -1e5 vs gap -1000",
                              fx$name))
  }
})

test_that("extreme finite log-likelihoods retain normalized posterior moments", {
  # One person, two binary tasks. Task 1 chooses an extremely unlikely
  # alternative; task 2 retains ordinary curvature. Analytical binary-logit
  # derivatives avoid finite differences of an objective of order 1e9.
  fx <- list(
    theta = c(1, log(0.3)), X = matrix(c(-1000, 0, 0.1, 0), 4L),
    W = matrix(c(0.2, -0.2, 0.3, -0.3), 4L), alt_idx = rep(1:2, 2L),
    choice_idx = c(1L, 2L), M = c(2L, 2L), Ti = 2L, weights = c(1, 1),
    K_w = 1L, S = 20L, eta = get_halton_normals(20L, 1L, 1L),
    rc_dist = 0L, rc_correlation = FALSE, rc_mean = FALSE, use_asc = FALSE,
    include_outside_option = FALSE)
  z <- 0.3 * drop(fx$eta)
  v2 <- (0.1 + 0.3 * z) - (-0.3 * z)
  p2 <- plogis(v2)
  for (gap in c(1000, 1e8, 1e9)) {
    fx$X[1L, 1L] <- -gap
    v1 <- (-gap + 0.2 * z) - (-0.2 * z)
    # At these gaps P_1 is zero to machine precision, log(P_1) = v1.
    logp2 <- -pmax(v2, 0) - log1p(exp(-abs(v2)))
    lambda <- v1 + logp2
    w <- exp(lambda - max(lambda))
    omega <- w / sum(w)
    G <- rbind(-gap - 0.1 * p2, 0.4 * z - 0.6 * z * p2)
    score <- drop(G %*% omega)
    # Covariance is translation invariant: remove the common score first.
    D <- G - G[, 1L]
    D <- D - drop(D %*% omega)
    H <- matrix(0, 2L, 2L)
    for (s in seq_len(fx$S)) {
      q <- c(0.1, 0.6 * z[s])
      H_s <- -p2[s] * (1 - p2[s]) * tcrossprod(q)
      H_s[2L, 2L] <- H_s[2L, 2L] + G[2L, s]
      H <- H - omega[s] * (H_s + tcrossprod(D[, s]))
    }
    mean_z <- sum(omega * z)
    sd_z <- sqrt(sum(omega * (z - mean_z)^2))
    for (gen in c(FALSE, TRUE)) {
      for (b in c(1L, 7L, fx$S)) {
        label <- sprintf("gap %g, generate %s, batch %d", gap, gen, b)
        res <- mxlp_call("gradient", fx, generate = gen, draw_batch = b)
        expect_false(res$overflow)
        mxlp_expect_close(res$gradient, -score, 1e-8, label)
        mxlp_expect_close(mxlp_call("scores", fx, generate = gen, draw_batch = b),
                          score, 1e-8, label)
        mxlp_expect_close(mxlp_call("hessian", fx, generate = gen, draw_batch = b),
                          H, 1e-8, label)
        mxlp_expect_close(mxlp_call("bhhh", fx, generate = gen, draw_batch = b),
                          tcrossprod(score), 1e-8, label)
        ct <- mxlp_call("tastes", fx, generate = gen, draw_batch = b)
        mxlp_expect_close(ct$mean, mean_z, 1e-8, label)
        mxlp_expect_close(ct$sd, sd_z, 1e-8, label)
      }
    }
  }
  # With W = 0 the draws are irrelevant to the likelihood. This checks the
  # weight sum directly through the score, including the streaming path,
  # and ensures centering does not create curvature from identical scores.
  fx$W[,] <- 0
  fx$eta[] <- 1
  for (gap in c(1e9, 1e16)) {
    fx$X[1L, 1L] <- -gap
    score <- c(-gap - 0.1 * plogis(0.1), 0)
    H <- diag(c(0.01 * plogis(0.1) * plogis(-0.1), 0))
    for (b in c(1L, 7L, fx$S)) {
      res <- mxlp_call("gradient", fx, draw_batch = b)
      mxlp_expect_close(res$gradient, -score, 1e-14, "draw-independent score")
      mxlp_expect_close(mxlp_call("hessian", fx, draw_batch = b),
                        H, 1e-14, "draw-independent Hessian")
      ct <- mxlp_call("tastes", fx, draw_batch = b)
      mxlp_expect_close(ct$mean, 0.3, 1e-14, "identical tastes")
      mxlp_expect_close(ct$sd, 0, 1e-14, "zero taste variance")
    }
  }
})

test_that("a zero-weight outlier does not determine Hessian centering", {
  skip_if_not_installed("numDeriv")
  # The first log-normal draw has a huge finite score but zero posterior
  # weight. Using it as a centering reference erases variation among the
  # other draws; permuting the same draws must leave curvature unchanged.
  eta <- c(1, 0, 0.02, -0.02)
  args <- list(
    theta = c(0.1, log(50)), X = matrix(c(1, 0), 2L),
    W = matrix(c(-1, 0), 2L), alt_idx = 1:2, choice_idx = 1L, M = 2L,
    weights = 1, eta_draws = array(eta, c(1L, 4L, 1L)), rc_dist = 1L,
    rc_correlation = FALSE, rc_mean = FALSE, use_asc = FALSE,
    include_outside_option = FALSE)
  grad <- function(theta) {
    args$theta <- theta
    drop(do.call(mxl_loglik_gradient_parallel, args)$gradient)
  }
  expected <- numDeriv::jacobian(grad, args$theta)
  for (order in list(1:4, c(2L, 3L, 4L, 1L))) {
    args$eta_draws[] <- eta[order]
    for (b in c(1L, 4L)) {
      args$draw_batch <- b
      mxlp_expect_close(do.call(mxl_hessian_parallel, args), expected, 1e-6,
                        "Hessian independent of the zero-weight draw's position")
    }
  }
})

test_that("overflowing utilities trigger the sentinel, the Hessian skip and NA tastes", {
  # With log-space probabilities a unit's log-likelihood is non-finite only
  # when its utilities are. Here the chosen alternative's fixed utility
  # overflows to -Inf (x = -1e308, beta_1 = 2), so every draw gives the
  # unit's choices zero probability. The objective returns the optimizer
  # sentinel, the Hessian leaves the unit out, and the unit's conditional
  # tastes are NA while every other unit's stay exact.
  for (panel in c(FALSE, TRUE)) {
    fx <- mxlp_edge_fixture(-36, panel)
    u <- fx$edge_unit
    t <- fx$edge_situation
    fx$theta[1L] <- 2
    fx$X[sum(fx$M[seq_len(t - 1L)]) + fx$choice_idx[t], 1L] <- -1e308
    expect_true(all(mxlp_situation_logp(fx$theta, fx, t) == -Inf),
                label = sprintf("[%s] chosen utility is -Inf at every draw",
                                fx$name))

    res <- mxlp_call("gradient", fx)
    expect_equal(res$objective, 1e10,
                 label = sprintf("[%s] sentinel objective", fx$name))
    expect_true(all(res$gradient == 0),
                label = sprintf("[%s] sentinel gradient is zero", fx$name))

    mxlp_expect_close(mxlp_call("hessian", fx),
                      mxlp_call("hessian", mxlp_drop_unit(fx, u)), 1e-10,
                      sprintf("[%s] Hessian vs Hessian without unit %d",
                              fx$name, u))

    ct <- mxlp_call("tastes", fx)
    orc <- mxlp_oracle(fx$theta, fx)
    expect_true(all(is.na(ct$mean[, u])) && all(is.na(ct$sd[, u])),
                label = sprintf("[%s] conditional tastes of unit %d are NA",
                                fx$name, u))
    mxlp_expect_close(ct$mean[, -u], orc$mean[, -u], 1e-10,
                      sprintf("[%s] other units' conditional mean vs oracle",
                              fx$name))
    mxlp_expect_close(ct$sd[, -u], orc$sd[, -u], 1e-10,
                      sprintf("[%s] other units' conditional SD vs oracle",
                              fx$name))
  }
})

# --- Input validation --------------------------------------------------------

test_that("the kernels accept double-typed index vectors", {
  # The integer inputs are read in place when they are integer (the output of
  # prepare_mxl_data()) and coerced once per call when they are not.
  on.exit(set_num_threads(2L), add = TRUE)
  set_num_threads(1L)
  fx <- mxlp_long_panel_fixture()
  fd <- fx
  fd$alt_idx <- as.double(fx$alt_idx)
  fd$choice_idx <- as.double(fx$choice_idx)
  for (k in mxlp_kernels) {
    expect_identical(mxlp_call(k, fd), mxlp_call(k, fx),
                     label = sprintf("%s with double-typed indices", k))
  }
})

test_that("input errors give the same message in every kernel and draw mode", {
  # One broken input per case, in the checks' order: situations and design
  # height, index lengths, alternative codes, choices, the delta block, the
  # draws, an alternative-level W. The last three cases break two checks at
  # once and pin which one reports first.
  fx <- mxlp_cells[[1L]]   # row-aligned W, no outside option, ASCs
  fo <- mxlp_cells[[2L]]   # outside option
  fa <- mxlp_cells[[3L]]   # alternative-level W, ASCs
  n <- sum(fx$M)
  N <- fx$N
  kname <- c(gradient = "mxl_loglik_gradient_parallel",
             hessian = "mxl_hessian_parallel", bhhh = "mxl_bhhh_parallel",
             scores = "mxl_scores_parallel",
             tastes = "mxl_conditional_tastes_parallel")
  err <- function(k, f, ...) {
    tryCatch({
      mxlp_call(k, f, ...)
      "<no error>"
    }, error = conditionMessage)
  }
  with_field <- function(f, name, value) {
    f[[name]] <- value
    f
  }
  bad_choice <- function(k) {
    sprintf("Invalid chosen alternative index for individual 3 (%s)", kname[[k]])
  }
  delta_msg <- "Theta's delta (ASC) block implies %d alternatives but alt_idx references alternative %d."
  cases <- list(
    list("M = 0", with_field(fx, "M", replace(fx$M, 3L, 0L)), list(),
         "M must be positive for every individual (M[3] = 0)."),
    list("M = NA", with_field(fx, "M", replace(fx$M, 3L, NA_integer_)), list(),
         "M must be positive for every individual (M[3] = -2147483648)."),
    list("X one row short", with_field(fx, "X", fx$X[-1L, , drop = FALSE]),
         list(), sprintf("X has %d rows but sum(M) is %d.", n - 1L, n)),
    list("alt_idx one short", with_field(fx, "alt_idx", fx$alt_idx[-1L]),
         list(), sprintf(paste("alt_idx length (%d) does not match the number",
                               "of rows of X (%d)."), n - 1L, n)),
    list("choice_idx one short",
         with_field(fx, "choice_idx", fx$choice_idx[-1L]), list(),
         sprintf("choice_idx length (%d) does not match N (%d)", N - 1L, N)),
    list("alternative code 0", with_field(fx, "alt_idx", replace(fx$alt_idx, 5L, 0L)),
         list(), "alt_idx must use 1-based alternative indices (found 0)."),
    list("alternative code NA",
         with_field(fx, "alt_idx", replace(fx$alt_idx, 5L, NA_integer_)), list(),
         "alt_idx must use 1-based alternative indices (found NA)."),
    list("alternative code -1",
         with_field(fx, "alt_idx", replace(fx$alt_idx, 5L, -1L)), list(),
         "alt_idx must use 1-based alternative indices (found -1)."),
    list("choice past the set",
         with_field(fx, "choice_idx", replace(fx$choice_idx, 4L, fx$M[4L] + 1L)),
         list(), bad_choice),
    list("choice NA",
         with_field(fx, "choice_idx", replace(fx$choice_idx, 4L, NA_integer_)),
         list(), bad_choice),
    list("choice past the set, outside option",
         with_field(fo, "choice_idx", replace(fo$choice_idx, 4L, fo$M[4L] + 1L)),
         list(), bad_choice),
    list("choice NA, outside option",
         with_field(fo, "choice_idx", replace(fo$choice_idx, 4L, NA_integer_)),
         list(), bad_choice),
    list("delta block one short", fx, list(theta = fx$theta[-length(fx$theta)]),
         sprintf(delta_msg, fx$J - 1L, fx$J)),
    list("alternative-level W one row short",
         with_field(fa, "W", fa$W[-fa$J, , drop = FALSE]), list(),
         sprintf(paste("W must be row-aligned with X (%d rows) or contain one",
                       "row per global alternative (at least %d rows); got %d",
                       "rows."), sum(fa$M), fa$J, fa$J - 1L)),
    list("a code past both the delta block and W: delta first",
         with_field(fa, "alt_idx", replace(fa$alt_idx, 5L, fa$J + 1L)), list(),
         sprintf(delta_msg, fa$J, fa$J + 1L)),
    list("M one situation too long: Ti first",
         with_field(fx, "M", c(fx$M, 2L)), list(),
         sprintf(paste("sum(Ti) (%d) does not match the number of choice",
                       "situations (%d)."), N, N + 1L)),
    list(paste("two situations merged in M (cross-section): the weights'",
               "length first, then choice_idx's"),
         with_field(fx, "M", c(fx$M[1L] + fx$M[2L], fx$M[-(1:2)])),
         list(Ti = NULL),
         function(k) {
           if (k %in% c("gradient", "hessian", "bhhh")) {
             sprintf("weights length (%d) does not match N (%d)", N, N - 1L)
           } else {
             sprintf("choice_idx length (%d) does not match N (%d)", N, N - 1L)
           }
         })
  )
  for (cs in cases) {
    for (gen in c(FALSE, TRUE)) {
      for (k in mxlp_kernels) {
        expected <- if (is.function(cs[[4L]])) cs[[4L]](k) else cs[[4L]]
        msg <- do.call(err, c(list(k, cs[[2L]], generate = gen), cs[[3L]]))
        expect_identical(msg, expected, label = sprintf(
          "%s: %s (%s)", cs[[1L]], k, if (gen) "generate" else "store"))
      }
    }
  }

  # Weighted kernels: the weights' length.
  for (gen in c(FALSE, TRUE)) {
    for (k in c("gradient", "hessian", "bhhh")) {
      expect_identical(err(k, fx, weights = fx$weights[-1L], generate = gen),
                       sprintf("weights length (%d) does not match N (%d)",
                               N - 1L, N),
                       label = paste(k, "with one weight short"))
    }
  }
  # Store mode: the draw cube.
  for (k in mxlp_kernels) {
    expect_identical(
      err(k, fx, eta = fx$eta[, , -1L, drop = FALSE]),
      sprintf(paste("eta_draws 3rd dimension (%d) does not match the number",
                    "of decision makers (%d)"), fx$U - 1L, fx$U),
      label = paste(k, "with one slice short"))
    expect_identical(
      err(k, fx, Ti = NULL, eta = fx$eta),
      sprintf("eta_draws 3rd dimension (%d) does not match N (%d)", fx$U, N),
      label = paste(k, "cross-section with one slice per decision maker"))
    expect_identical(
      err(k, fx, eta = array(0, c(fx$K_w + 1L, fx$S, fx$U))),
      sprintf("eta_draws 1st dimension (%d) does not match K_w (%d)",
              fx$K_w + 1L, fx$K_w),
      label = paste(k, "with a draw row too many"))
  }
  # Generate mode: the draw count and the primes table.
  f0 <- with_field(fx, "S", 0L)
  fw <- fx
  fw$W <- cbind(fx$W, matrix(0, nrow(fx$W), 129L - fx$K_w))
  fw$K_w <- 129L
  fw$rc_dist <- rep(0L, 129L)
  fw$theta <- c(fx$theta[1:2], rep(log(0.5), 129L), fx$theta[-(1:(2 + fx$K_w))])
  for (k in mxlp_kernels) {
    expect_identical(err(k, f0, generate = TRUE),
                     "gen_S must be positive when gen_seed >= 0",
                     label = paste(k, "with gen_S = 0"))
    expect_identical(err(k, fw, generate = TRUE),
                     paste("K_w exceeds the primes table size (128); reduce",
                           "K_w or extend the primes table."),
                     label = paste(k, "with 129 random coefficients"))
  }
})

test_that("alternative codes are checked in parallel above 10^6 rows", {
  # Past 10^6 stacked rows the scan for the smallest and largest code runs
  # in an OpenMP reduction; a bad code in the last row must still be found.
  N <- 250001L
  M <- rep(4L, N)
  n <- sum(M)
  X <- matrix(0, n, 1L)
  alt <- rep(1:4, N)
  call_grad <- function(alt_idx, W = matrix(0, n, 1L), use_asc = TRUE,
                        theta = c(0, log(0.5), 0, 0, 0)) {
    tryCatch({
      mxl_loglik_gradient_parallel(
        theta, X, W, alt_idx, rep(1L, N), M, rep(1, N),
        array(0, c(1L, 0L, 0L)), 0L, rc_correlation = FALSE,
        rc_mean = FALSE, use_asc = use_asc, include_outside_option = FALSE,
        gen_seed = 0L, gen_scramble = 0L, gen_S = 1L)
      "<no error>"
    }, error = conditionMessage)
  }
  expect_identical(call_grad(replace(alt, n, NA_integer_)),
                   "alt_idx must use 1-based alternative indices (found NA).")
  expect_identical(call_grad(replace(alt, n, 0L)),
                   "alt_idx must use 1-based alternative indices (found 0).")
  expect_identical(call_grad(replace(alt, n, 5L)),
                   paste("Theta's delta (ASC) block implies 4 alternatives but",
                         "alt_idx references alternative 5."))
  expect_identical(
    call_grad(replace(alt, n, 5L), W = matrix(0, 4L, 1L), use_asc = FALSE,
              theta = c(0, log(0.5))),
    sprintf(paste("W must be row-aligned with X (%d rows) or contain one row",
                  "per global alternative (at least 5 rows); got 4 rows."), n))
})

test_that("panel kernels reject non-positive or missing Ti entries", {
  fx <- mxlp_cells[[2L]]
  Ti <- fx$Ti
  # Same length and sum as Ti, so only positivity is violated. Uniform
  # weights keep the regrouped units weight-constant.
  Ti_zero <- Ti
  Ti_zero[2L] <- Ti[1L] + Ti[2L]
  Ti_zero[1L] <- 0L
  Ti_neg <- Ti
  Ti_neg[2L] <- Ti[1L] + Ti[2L] + 1L
  Ti_neg[1L] <- -1L
  Ti_na <- Ti
  Ti_na[2L] <- NA_integer_
  w1 <- rep(1, fx$N)
  for (k in mxlp_kernels) {
    expect_error(mxlp_call(k, fx, Ti = Ti_zero, weights = w1),
                 "Ti must be positive", label = paste(k, "with Ti[1] = 0"))
    expect_error(mxlp_call(k, fx, Ti = Ti_neg, weights = w1),
                 "Ti must be positive", label = paste(k, "with Ti[1] = -1"))
    expect_error(mxlp_call(k, fx, Ti = Ti_na, weights = w1),
                 "Ti must be positive", label = paste(k, "with Ti[2] = NA"))
  }
})

test_that("panel kernels reject an empty Ti (needs the empty-Ti guard)", {
  # Ti = integer(0) describes no decision maker at all; it must not fall
  # back to the cross-section. The draw cube is empty too, so no other check
  # can fire first. Each kernel's message is captured so that every kernel
  # reports its own result.
  fx <- mxlp_cells[[2L]]
  eta_empty <- fx$eta[, , integer(0), drop = FALSE]
  for (k in mxlp_kernels) {
    msg <- tryCatch({
      mxlp_call(k, fx, Ti = integer(0), weights = rep(1, fx$N),
                eta = eta_empty)
      "<no error>"
    }, error = conditionMessage)
    expect_match(msg, "at least one respondent",
                 label = paste(k, "error message with Ti = integer(0)"))
  }
})

test_that("panel kernels reject Ti not summing to the number of situations", {
  fx <- mxlp_cells[[2L]]
  Ti_long <- fx$Ti
  Ti_long[fx$U] <- Ti_long[fx$U] + 1L
  u <- fx$probe_units[2L]  # an interior unit with T_u >= 2
  Ti_short <- fx$Ti
  Ti_short[u] <- Ti_short[u] - 1L
  w1 <- rep(1, fx$N)
  for (k in mxlp_kernels) {
    expect_error(mxlp_call(k, fx, Ti = Ti_long, weights = w1),
                 "does not match the number of choice situations",
                 label = paste(k, "with sum(Ti) = N + 1"))
    expect_error(mxlp_call(k, fx, Ti = Ti_short, weights = w1),
                 "does not match the number of choice situations",
                 label = paste(k, "with sum(Ti) = N - 1"))
  }
})

test_that("panel kernels reject weights that vary within a decision maker", {
  fx <- mxlp_cells[[2L]]
  u <- fx$probe_units[2L]  # an interior unit with T_u >= 2
  first <- sum(fx$Ti[seq_len(u - 1L)]) + 1L
  w_bad <- fx$weights
  w_bad[first + 1L] <- 1.5 * w_bad[first + 1L]
  for (k in c("gradient", "hessian", "bhhh")) {
    expect_error(mxlp_call(k, fx, weights = w_bad),
                 "constant within each decision maker",
                 label = paste(k, "with a within-unit weight change"))
  }
})

test_that("store-mode draws need one slice per decision maker", {
  fx <- mxlp_cells[[2L]]
  eta_situations <- get_halton_normals(fx$S, fx$N, fx$K_w)  # N != U slices
  eta_short <- fx$eta[, , -1L, drop = FALSE]                  # U - 1 slices
  for (k in mxlp_kernels) {
    expect_error(mxlp_call(k, fx, eta = eta_situations),
                 "eta_draws 3rd dimension",
                 label = paste(k, "with one slice per situation"))
    expect_error(mxlp_call(k, fx, eta = eta_short),
                 "eta_draws 3rd dimension",
                 label = paste(k, "with U - 1 slices"))
  }
  # The cross-section keeps its check: one slice per situation.
  expect_error(mxlp_call("gradient", fx, Ti = NULL, eta = fx$eta),
               "eta_draws 3rd dimension")
})
