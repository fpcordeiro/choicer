# The default start of choicer 0.2.x (theta_init = NULL there): zeros, with
# L_pp = 0.5 on the Cholesky diagonal of `fit`'s parameter layout. Tests whose
# unscaled ("none") fits must reach the maximum of a badly scaled design start
# them here: from the unit-aware default the unscaled optimizer can drive a
# random coefficient's variance onto its zero plateau, where the gradient of
# log L_pp vanishes and the information is singular.
mxl_start_0_2 <- function(fit) {
  theta <- rep(0, length(coef(fit)))
  s <- fit$param_map$sigma
  K_w <- length(fit$rc_dist)
  theta[if (fit$rc_correlation) s[cumsum(seq_len(K_w))] else s] <- log(0.5)
  theta
}
