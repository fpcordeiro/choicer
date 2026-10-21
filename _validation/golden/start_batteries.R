# 6D follow-up 4 (unit-aware default start): run_mxlogit()'s column steps
# (design_column_step()) and its default start (.mxl_default_start()) on fixed
# matrices and parameter layouts. Absent on main, where every entry records
# the error (deliberate, class D).
run_start_unit <- function() {
  out <- list()
  set.seed(2026)
  n <- 500L
  mats <- list(
    dense = cbind(rnorm(n), 1e3 * runif(n), 1e8 + rnorm(n)),
    dummy = cbind(as.numeric(seq_len(n) %% 37L == 0L),
                  1 - as.numeric(seq_len(n) %% 3L == 0L)),
    sparse = cbind(ifelse(seq_len(n) %% 4L == 0L, 2 + rexp(n), 0)),
    balanced = cbind(rep(c(0, 1), length.out = n)),
    edge = cbind(rep(3, n), c(rep(3, n - 1), 3 * (1 + 4 * .Machine$double.eps)),
                 c(Inf, rnorm(n - 1)), c(NA, rnorm(n - 1))),
    rows0 = matrix(numeric(), 0, 2),
    rows1 = matrix(c(1, 2), 1, 2))
  for (nm in names(mats)) {
    out[[paste0("step|", nm)]] <- p2_val(choicer:::design_column_step(mats[[nm]]))
  }
  W <- cbind(mats$dense[, 1:2], mats$dummy[, 1])
  layouts <- list(
    diag = list(n = 8L, pm = list(beta = 1:3, sigma = 4:6, asc = 7:8),
                corr = FALSE, rcd = c(0L, 0L, 0L)),
    corr_mu_ln = list(n = 14L, pm = list(beta = 1:3, mu = 4:6, sigma = 7:12,
                                         asc = 13:14),
                      corr = TRUE, rcd = c(0L, 1L, 0L)))
  for (nm in names(layouts)) {
    L <- layouts[[nm]]
    out[[paste0("start|", nm)]] <- p2_val(
      choicer:::.mxl_default_start(L$n, L$pm, W, L$rcd, L$corr))
  }
  out
}
