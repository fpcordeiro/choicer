simulate_mnlogit_benchmark_data <- function(N, J, beta = c(0.8, -0.6),
                                            delta = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- as.integer(N)
  J <- as.integer(J)
  K <- length(beta)
  if (N <= 0L || J < 2L) stop("N must be positive and J must be at least 2.")

  if (is.null(delta)) {
    delta <- c(0, rep(c(0.4, -0.3, 0.2, -0.15), length.out = J - 1L))
  }
  if (length(delta) != J) stop("delta must have length J.")
  delta <- delta - delta[[1L]]

  x_cols <- paste0("x", seq_len(K))
  dt <- data.table::data.table(
    id = rep(seq_len(N), each = J),
    alt = rep(seq_len(J), times = N)
  )
  for (k in seq_len(K)) {
    dt[, (x_cols[[k]]) := stats::runif(.N, -1, 1)]
  }

  X <- as.matrix(dt[, x_cols, with = FALSE])
  dt[, deterministic_utility := as.numeric(X %*% beta) + delta[alt]]
  dt[, epsilon := -log(-log(stats::runif(.N)))]
  dt[, utility := deterministic_utility + epsilon]
  dt[, choice := as.integer(seq_len(.N) == which.max(utility)), by = id]
  dt[, c("deterministic_utility", "epsilon", "utility") := NULL]
  data.table::setcolorder(dt, c("id", "alt", "choice", x_cols))
  data.table::setkey(dt, id, alt)

  list(
    data = dt,
    x_cols = x_cols,
    true_params = list(beta = stats::setNames(beta, x_cols), delta = delta),
    settings = list(N = N, J = J, K_x = K)
  )
}

mnlogit_wide_data <- function(dt, x_cols, J) {
  wide <- data.table::dcast(
    dt,
    id ~ alt,
    value.var = x_cols,
    sep = "_"
  )
  choices <- dt[choice == 1L, .(id, choice_alt = alt)]
  wide <- merge(wide, choices, by = "id", sort = FALSE)
  wide[, ID := id]
  wide[, CHOICE := as.integer(choice_alt)]
  for (j in seq_len(J)) {
    wide[, paste0("avail_", j) := 1L]
  }
  data.table::setorder(wide, id)
  wide[]
}
