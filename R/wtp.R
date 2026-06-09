
# Willingness-to-pay (WTP) with delta-method standard errors

# --- Internal helpers ---------------------------------------------------------

#' Delta-method estimate and SE for a WTP ratio
#'
#' Computes the WTP ratio and its delta-method standard error from the
#' variance block \code{V} of the parameters involved. Three numerator kinds
#' are supported, matching the package's random-coefficient parameterization
#' (see \code{src/mxlogit.cpp}: log-normal coefficients are the shifted
#' log-normal \code{beta_k = exp(mu_k) + exp((L eta)_k)}):
#' \itemize{
#'   \item \code{"linear"}: \code{g = -num/den}; gradients
#'     \code{dg/dnum = -1/den}, \code{dg/dden = num/den^2}.
#'   \item \code{"logshift"} (shifted log-normal median, numerator
#'     \code{exp(num) + 1}): \code{g = -(exp(num) + 1)/den}; gradients
#'     \code{dg/dnum = -exp(num)/den},
#'     \code{dg/dden = (exp(num) + 1)/den^2}.
#'   \item \code{"unit"} (log-normal with \code{rc_mean = FALSE}: median
#'     numerator is exactly 1, no free parameter): \code{g = -1/den};
#'     gradient \code{dg/dden = 1/den^2}; \code{V} is the 1x1 price block.
#' }
#'
#' @param num Numerator coefficient (attribute coefficient, or mu for the
#'   shifted log-normal median case; ignored for \code{kind = "unit"}).
#' @param den Denominator coefficient (price coefficient).
#' @param V Variance-covariance block: 2x2 of \code{c(num, den)} for
#'   \code{"linear"}/\code{"logshift"}, 1x1 of \code{den} for \code{"unit"}.
#' @param kind One of \code{"linear"}, \code{"logshift"}, \code{"unit"}.
#' @returns List with \code{estimate} and \code{se} (NA when the variance is
#'   unavailable or negative).
#' @noRd
.wtp_delta_ratio <- function(num, den, V, kind = "linear") {
  if (kind == "unit") {
    estimate <- -1 / den
    grad <- 1 / den^2
  } else if (kind == "logshift") {
    e <- exp(num)
    estimate <- -(e + 1) / den
    grad <- c(-e / den, (e + 1) / den^2)
  } else {
    estimate <- -num / den
    grad <- c(-1 / den, num / den^2)
  }
  var_g <- as.numeric(t(grad) %*% V %*% grad)
  se <- if (is.finite(var_g) && var_g >= 0) sqrt(var_g) else NA_real_

  list(estimate = estimate, se = se)
}

#' Validate the price variable and return its coefficient index
#' @noRd
.wtp_price_index <- function(object, price_var) {
  if (!is.character(price_var) || length(price_var) != 1L || is.na(price_var)) {
    stop("'price_var' must be a single variable name (character).")
  }
  cf_names <- names(object$coefficients)
  beta_idx <- object$param_map$beta
  pos <- match(price_var, cf_names[beta_idx])
  if (is.na(pos)) {
    stop("Price variable '", price_var, "' not found among fixed-coefficient ",
         "variables. Available: ",
         paste(cf_names[beta_idx], collapse = ", "))
  }
  beta_idx[pos]
}

#' Resolve attribute variables into WTP row specifications
#'
#' Each row is a list with \code{label} (display name), \code{idx} (coefficient
#' index into \code{coef(object)}, or NA for \code{kind = "unit"}), and
#' \code{kind} (\code{"linear"}, \code{"logshift"} for shifted log-normal
#' medians with an estimated location, or \code{"unit"} for log-normal
#' coefficients without one, whose median numerator is exactly 1).
#'
#' @param object A choicer_fit object.
#' @param price_var Validated price variable name.
#' @param attr_vars User-requested attributes or NULL (defaults).
#' @param w_names Random-coefficient (W) column names; empty for MNL/NL.
#' @returns List of row specifications.
#' @noRd
.wtp_resolve_rows <- function(object, price_var, attr_vars,
                              w_names = character(0)) {
  cf_names <- names(object$coefficients)
  pm <- object$param_map
  beta_names <- cf_names[pm$beta]
  asc_names <- if (!is.null(pm$asc)) cf_names[pm$asc] else character(0)
  mu_names <- if (!is.null(pm$mu)) cf_names[pm$mu] else character(0)
  has_mu <- isTRUE(object$rc_mean) && !is.null(pm$mu)

  make_row <- function(label, idx, kind) {
    list(label = label, idx = idx, kind = kind)
  }
  # Random-coefficient row for W variable k. Shifted log-normal
  # (beta_k = exp(mu_k) + exp((L eta)_k)): median numerator exp(mu_k) + 1
  # with rc_mean = TRUE, exactly 1 without.
  rc_row <- function(k) {
    if (object$rc_dist[k] == 1) {
      if (has_mu) {
        make_row(w_names[k], pm$mu[k], "logshift")
      } else {
        make_row(w_names[k], NA_integer_, "unit")
      }
    } else {
      if (!has_mu) return(NULL)  # normal RC mean is 0 by construction
      make_row(cf_names[pm$mu[k]], pm$mu[k], "linear")
    }
  }

  rows <- list()
  add_row <- function(row) {
    if (!is.null(row)) rows[[length(rows) + 1L]] <<- row
  }

  if (is.null(attr_vars)) {
    for (v in setdiff(beta_names, price_var)) {
      add_row(make_row(v, pm$beta[match(v, beta_names)], "linear"))
    }
    for (k in seq_along(w_names)) {
      add_row(rc_row(k))
    }
    return(rows)
  }

  if (!is.character(attr_vars)) {
    stop("'attr_vars' must be a character vector of variable names or NULL.")
  }
  for (v in attr_vars) {
    if (identical(v, price_var)) {
      stop("'attr_vars' must not include the price variable '", price_var, "'.")
    }
    if (v %in% beta_names) {
      add_row(make_row(v, pm$beta[match(v, beta_names)], "linear"))
    } else if (v %in% asc_names) {
      add_row(make_row(v, pm$asc[match(v, asc_names)], "linear"))
    } else if (v %in% w_names || v %in% mu_names) {
      k <- if (v %in% w_names) match(v, w_names) else match(v, mu_names)
      row <- rc_row(k)
      if (is.null(row)) {
        stop("Variable '", v, "' is a normal random coefficient without an ",
             "estimated mean (rc_mean = FALSE), so its mean WTP is 0 by ",
             "construction. Refit with rc_mean = TRUE to obtain its WTP.")
      }
      add_row(row)
    } else {
      stop("Variable '", v, "' not found. Available: ",
           paste(c(setdiff(beta_names, price_var), asc_names,
                   if (length(w_names) > 0) w_names else mu_names),
                 collapse = ", "))
    }
  }
  rows
}

#' Assemble the WTP table from resolved rows
#' @noRd
.wtp_build_table <- function(object, price_var, rows, level) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.")
  }

  object <- ensure_vcov(object)
  cf <- object$coefficients
  V <- object$vcov
  price_idx <- .wtp_price_index(object, price_var)
  theta_p <- unname(cf[price_idx])

  n <- length(rows)
  labels <- vapply(rows, function(r) r$label, character(1))
  is_median <- vapply(rows, function(r) r$kind %in% c("logshift", "unit"),
                      logical(1))
  est <- rep(NA_real_, n)
  se <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    r <- rows[[i]]
    par_idx <- if (r$kind == "unit") price_idx else c(r$idx, price_idx)
    Vr <- if (!is.null(V)) {
      V[par_idx, par_idx, drop = FALSE]
    } else {
      matrix(NA_real_, length(par_idx), length(par_idx))
    }
    num <- if (r$kind == "unit") NA_real_ else unname(cf[r$idx])
    res <- .wtp_delta_ratio(num, theta_p, Vr, kind = r$kind)
    est[i] <- res$estimate
    se[i] <- res$se
  }

  q <- stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(
    Estimate  = est,
    Std_Error = se,
    z_value   = est / se,
    CI_lower  = est - q * se,
    CI_upper  = est + q * se,
    row.names = labels,
    stringsAsFactors = FALSE
  )
  attr(out, "price_var") <- price_var
  attr(out, "level") <- level
  attr(out, "median_rows") <- labels[is_median]
  class(out) <- c("choicer_wtp", "data.frame")
  out
}

# --- Generic and methods ------------------------------------------------------

#' Compute willingness to pay
#'
#' Computes willingness-to-pay (WTP) ratios with delta-method standard errors
#' from a fitted choice model. For an attribute coefficient
#' \eqn{\theta_k} and a price coefficient \eqn{\theta_p}, the WTP is
#' \deqn{WTP_k = -\theta_k / \theta_p,}
#' the marginal rate of substitution between the attribute and price. Standard
#' errors use the delta method with analytic gradients
#' \eqn{\partial g/\partial \theta_k = -1/\theta_p} and
#' \eqn{\partial g/\partial \theta_p = \theta_k/\theta_p^2}, applied to the
#' corresponding 2x2 block of \code{vcov(object)}.
#'
#' For mixed logit models, random coefficients are included via their
#' estimated location parameters. The package's log-normal random coefficient
#' is the \emph{shifted} log-normal
#' \eqn{\beta_k = \exp(\mu_k) + \exp((L\eta)_k)} (see
#' \code{run_mxlogit()}), so:
#' \itemize{
#'   \item Normal random coefficient \eqn{k} (\code{rc_mean = TRUE}): mean WTP
#'     \eqn{-\mu_k / \theta_p}, labeled \code{Mu_x}.
#'   \item Log-normal random coefficient \eqn{k} (\code{rc_mean = TRUE}):
#'     \strong{median} WTP \eqn{-(\exp(\mu_k) + 1) / \theta_p}, since the
#'     median of \eqn{\exp((L\eta)_k)} is 1. (The mean,
#'     \eqn{\exp(\mu_k) + \exp(\sigma_k^2/2)}, is highly sensitive to the
#'     estimated variance; the median is the more robust summary.) These rows
#'     are labeled by the attribute name and flagged as medians when printed.
#'   \item Log-normal random coefficient with \code{rc_mean = FALSE}:
#'     \eqn{\beta_k = \exp((L\eta)_k)} has median 1, so the median WTP is
#'     \eqn{-1/\theta_p} with uncertainty driven solely by \eqn{\theta_p}.
#' }
#' Normal random coefficients with \code{rc_mean = FALSE} have mean 0 by
#' construction and are excluded from the table.
#'
#' The price variable must have a \emph{fixed} coefficient. A random price
#' coefficient is rejected: the ratio of two random coefficients generally has
#' no finite moments (the denominator has positive density at 0), so mean or
#' median WTP computed from location parameters would be meaningless. Use a
#' fixed price coefficient, or estimate the model in WTP space.
#'
#' @param object A fitted model object (\code{choicer_mnl}, \code{choicer_mxl},
#'   or \code{choicer_nl}).
#' @param price_var Name of the price variable. Must be a fixed-coefficient
#'   variable (a column of the design matrix \code{X}).
#' @param attr_vars Character vector of attributes to report. Defaults to all
#'   fixed-coefficient variables other than \code{price_var} (plus, for mixed
#'   logit with \code{rc_mean = TRUE}, all random coefficients). ASC names
#'   (e.g. \code{"ASC_2"}) may also be supplied; the WTP of an ASC is
#'   \eqn{-ASC_j / \theta_p}.
#' @param level Confidence level for the normal-approximation interval
#'   \eqn{Estimate \pm z_{1-(1-level)/2} \times SE}. Default 0.95.
#' @param ... Additional arguments passed to methods.
#' @returns A data.frame of class \code{choicer_wtp} with one row per
#'   attribute and columns \code{Estimate}, \code{Std_Error}, \code{z_value},
#'   \code{CI_lower}, \code{CI_upper}. Attributes \code{price_var} and
#'   \code{level} record the inputs; \code{median_rows} lists rows that are
#'   median (rather than mean) WTP. Standard errors are NA when the
#'   variance-covariance matrix is unavailable.
#' @examples
#' \donttest{
#' library(data.table)
#' sim <- simulate_mnl_data(N = 1000, J = 4, beta = c(0.8, -0.6), seed = 123,
#'                          outside_option = FALSE, vary_choice_set = FALSE)
#' fit <- run_mnlogit(sim$data, "id", "alt", "choice", c("x1", "x2"))
#' # treat x2 as the price variable
#' wtp(fit, price_var = "x2")
#' wtp(fit, price_var = "x2", attr_vars = c("x1", "ASC_2"), level = 0.90)
#' }
#' @export
wtp <- function(object, price_var, attr_vars = NULL, level = 0.95, ...) {
  UseMethod("wtp")
}

#' @rdname wtp
#' @export
wtp.choicer_fit <- function(object, price_var, attr_vars = NULL,
                            level = 0.95, ...) {
  .wtp_price_index(object, price_var)
  rows <- .wtp_resolve_rows(object, price_var, attr_vars)
  .wtp_build_table(object, price_var, rows, level)
}

#' @rdname wtp
#' @export
wtp.choicer_mxl <- function(object, price_var, attr_vars = NULL,
                            level = 0.95, ...) {
  w_names <- names(object$sW) %||% colnames(object[["data"]]$W) %||% character(0)
  if (is.character(price_var) && length(price_var) == 1L &&
      price_var %in% w_names) {
    stop("Random price coefficients are not supported: the WTP ratio of two ",
         "random coefficients generally has no finite moments. Use a fixed ",
         "price coefficient (a 'covariate_cols' variable) or estimate the ",
         "model in WTP space.")
  }
  .wtp_price_index(object, price_var)
  rows <- .wtp_resolve_rows(object, price_var, attr_vars, w_names = w_names)
  out <- .wtp_build_table(object, price_var, rows, level)
  if (!isTRUE(object$rc_mean) && any(object$rc_dist == 0) &&
      is.null(attr_vars)) {
    attr(out, "rc_note") <- paste0(
      "Normal random-coefficient means are 0 by construction ",
      "(rc_mean = FALSE); those coefficients are excluded.")
  }
  out
}

# --- print --------------------------------------------------------------------

#' Print a WTP table
#'
#' @param x A \code{choicer_wtp} object.
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments passed to \code{print.data.frame}.
#' @returns The object invisibly.
#' @export
print.choicer_wtp <- function(x, digits = 4, ...) {
  level <- attr(x, "level") %||% 0.95
  cat(sprintf("Willingness to pay (WTP), price variable: '%s' (%g%% CI)\n",
              attr(x, "price_var"), 100 * level))
  if (nrow(x) == 0) {
    cat("No attributes to report.\n")
  } else {
    print(structure(x, class = "data.frame"), digits = digits, ...)
  }
  median_rows <- attr(x, "median_rows")
  if (length(median_rows) > 0) {
    cat("Note: median WTP (log-normal random coefficient): ",
        paste(median_rows, collapse = ", "), "\n", sep = "")
  }
  rc_note <- attr(x, "rc_note")
  if (!is.null(rc_note)) {
    cat("Note:", rc_note, "\n")
  }
  invisible(x)
}
