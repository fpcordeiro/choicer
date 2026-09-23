# Conditional (individual-level) tastes of a mixed logit fit: Revelt & Train
# (2000); Train (2009, ch. 11). The frequentist analogue of the hierarchical
# Bayes beta_i summaries of run_hmnlogit().

#' Conditional (individual-level) tastes from a mixed logit fit
#'
#' For each decision maker, the mean and standard deviation of the random
#' coefficients conditional on the choices they were observed to make,
#' \eqn{E[\beta_n \mid y_n, x_n; \hat\theta]} (Revelt and Train 2000; Train
#' 2009, ch. 11): the frequentist counterpart of the individual-level
#' \code{beta_i} summaries of \code{\link{run_hmnlogit}}.
#'
#' \strong{Estimand.} The fit estimates the population distribution of tastes,
#' \eqn{f(\beta \mid \theta)}. Among decision makers who face \eqn{x_n} and
#' make the choices \eqn{y_n}, tastes follow the conditional distribution
#' \deqn{h(\beta \mid y_n, x_n; \theta) =
#'   \frac{P(y_n \mid x_n, \beta) \, f(\beta \mid \theta)}
#'        {\int P(y_n \mid x_n, \beta) \, f(\beta \mid \theta) \, d\beta},}
#' where \eqn{P(y_n \mid x_n, \beta) = \prod_t P_{nt}(j_{nt} \mid \beta)} is
#' the product of the decision maker's logit probabilities. Its mean and
#' standard deviation are simulated at \eqn{\hat\theta} with the estimation
#' draws, as weighted averages over draws \eqn{\beta_{ns}} with weights
#' \eqn{\omega_{ns} \propto P(y_n \mid x_n, \beta_{ns})}. Coefficients are on
#' the scale they enter utility: \eqn{\mu_k + (L\eta)_k} for a normal
#' coefficient and \eqn{\exp(\mu_k) + \exp((L\eta)_k)} for the shifted
#' log-normal, with the \eqn{\mu} term present only when
#' \code{rc_mean = TRUE}. If a variable's mean is carried instead by the same
#' variable in \code{covariate_cols} (with \code{rc_mean = FALSE}), only its
#' random part is reported. This is the utility scale; the \code{beta_i}
#' summaries of \code{run_hmnlogit()} are on the scale of the underlying
#' normal, i.e. the log of the coefficient for log-normal coordinates.
#'
#' \strong{What it is not.} The conditional mean is not the decision maker's
#' own \eqn{\beta_n}. It is the average taste among people who would make
#' \eqn{n}'s choices, so it shrinks toward the population mean, the more so
#' the fewer choice situations \eqn{n} contributes, and its dispersion across
#' decision makers understates the population dispersion. By the law of
#' iterated expectations the conditional means average, in large samples, to
#' the population mean (\code{mean_cond} vs \code{model_mean} in
#' \code{population}), and by the law of total variance
#' \eqn{Var(\beta) = Var(E[\beta \mid y]) + E[Var(\beta \mid y)]}: the share
#' of taste variance revealed by choices (\code{revealed}) and the share left
#' unresolved (\code{residual}) should sum to about one. These are numerical
#' checks, not specification tests: for normal random coefficients with free
#' means, both identities are the likelihood's first-order conditions for
#' \eqn{\mu} and \eqn{\Sigma}, so they hold at a converged maximum whatever
#' the true taste distribution. A departure points to too few draws
#' \code{S} or an unconverged \eqn{\hat\theta}; when a random coefficient's
#' mean is fixed at zero (\code{rc_mean = FALSE} and no fixed coefficient on
#' the same variable), a \code{mean_cond} far from zero can also flag that
#' restriction. For log-normal coefficients the variance identity holds on
#' the underlying normal scale, not the utility scale reported here, so
#' expect departures of a few percent in either direction at moderate
#' \code{S}. The \code{sd} is the dispersion of the conditional distribution
#' of tastes, not a standard error: it ignores the sampling error in
#' \eqn{\hat\theta}.
#'
#' \strong{Cross-section.} Without \code{person_col} each choice situation is
#' its own decision maker. A single choice reveals little about tastes, so
#' the conditional means stay close to the population mean (\code{revealed}
#' small); fit with \code{person_col} to learn about tastes from repeated
#' choices.
#'
#' Units whose choices have zero simulated probability at every draw get
#' \code{NA} and are left out of the \code{population} summaries.
#'
#' @param object A \code{choicer_mxl} fit with stored data
#'   (\code{keep_data = TRUE}).
#' @param ... Unused.
#' @returns An object of class \code{choicer_tastes}: a list with
#'   \describe{
#'     \item{\code{mean}, \code{sd}}{\eqn{K_w \times U} matrices of
#'       conditional means and standard deviations; rows are the random
#'       coefficients (\code{random_var_cols}), columns the \eqn{U} likelihood
#'       units, named by decision-maker id (\code{person_col}) or by
#'       choice-situation id.}
#'     \item{\code{unit}}{\code{"decision maker"} or \code{"choice situation"}.}
#'     \item{\code{n_units}}{\eqn{U}.}
#'     \item{\code{S}}{Simulation draws per unit.}
#'     \item{\code{weights}}{Unit weights used in \code{population}.}
#'     \item{\code{population}}{A data.frame with one row per random
#'       coefficient: \code{model_mean} and \code{model_sd}, the moments of
#'       \eqn{f(\beta \mid \hat\theta)}; \code{mean_cond} and \code{sd_cond},
#'       the weighted mean and standard deviation of the conditional means
#'       across units; \code{revealed}
#'       \eqn{= Var(E[\beta \mid y]) / Var(\beta)}; and \code{residual}
#'       \eqn{= E[Var(\beta \mid y)] / Var(\beta)}.}
#'   }
#' @references Revelt, D. and Train, K. (2000). Customer-specific taste
#'   parameters and mixed logit: households' choice of electricity supplier.
#'   Working Paper E00-274, Department of Economics, University of California,
#'   Berkeley.
#'
#'   Train, K. (2009). \emph{Discrete Choice Methods with Simulation}, 2nd
#'   ed., Ch. 11. Cambridge University Press.
#' @seealso \code{\link{run_mxlogit}} (\code{person_col}),
#'   \code{\link{simulate_mxl_data}} (\code{true_params$gamma_i})
#' @examples
#' \donttest{
#' sim <- simulate_mxl_data(N = 150, T = 5, J = 3, seed = 1)
#' fit <- run_mxlogit(sim$data, "id", "alt", "choice", c("x1", "x2"),
#'                    c("w1", "w2"), outside_opt_label = 0L,
#'                    rc_correlation = TRUE, person_col = "pid", S = 50L)
#' ct <- conditional_tastes(fit)
#' ct
#' # conditional means track the realized tastes, with shrinkage
#' cor(ct$mean["w1", ], sim$true_params$gamma_i["w1", colnames(ct$mean)])
#' }
#' @export
conditional_tastes <- function(object, ...) UseMethod("conditional_tastes")

#' @rdname conditional_tastes
#' @export
conditional_tastes.default <- function(object, ...) {
  stop("conditional_tastes() is implemented for mixed logit fits from ",
       "run_mxlogit(); for hierarchical Bayes fits (run_hmnlogit(), ",
       "run_hmnprobit()) see the individual-level summaries in `object$beta_i`.",
       call. = FALSE)
}

#' @rdname conditional_tastes
#' @export
conditional_tastes.choicer_mxl <- function(object, ...) {
  d <- object[["data"]]
  if (is.null(d) || is.null(object$draws_info)) {
    stop("conditional_tastes() needs the stored data and draws; refit with ",
         "keep_data = TRUE.", call. = FALSE)
  }
  first <- .unit_first(d)
  n_units <- length(first)
  gp <- .mxl_gen_params(object$draws_info, N = n_units)
  ct <- mxl_conditional_tastes_parallel(
    theta = object$coefficients, X = d$X, W = d$W,
    alt_idx = d$alt_idx, choice_idx = d$choice_idx, M = d$M,
    eta_draws = gp$eta_draws, rc_dist = object$rc_dist,
    rc_correlation = object$rc_correlation, rc_mean = object$rc_mean,
    use_asc = object$use_asc,
    include_outside_option = object$include_outside_option,
    gen_seed = gp$gen_seed, gen_scramble = gp$gen_scramble, gen_S = gp$gen_S,
    Ti = d$Ti
  )

  panel <- !is.null(d$Ti)
  ids <- if (panel) d$person_ids else d$situation_ids
  if (is.null(ids)) ids <- seq_len(n_units)
  dn <- list(colnames(d$W), as.character(ids))
  cmean <- ct$mean
  csd <- ct$sd
  dimnames(cmean) <- dn
  dimnames(csd) <- dn
  w <- d$weights[first]

  structure(
    list(
      mean = cmean,
      sd = csd,
      unit = if (panel) "decision maker" else "choice situation",
      n_units = n_units,
      S = object$draws_info$S,
      weights = w,
      population = .taste_population(object, cmean, csd, w)
    ),
    class = "choicer_tastes"
  )
}

#' Population moments of the random coefficients vs their conditional means
#'
#' `model_mean` / `model_sd` are the moments of f(beta | theta-hat) on the
#' utility scale (normal: mu_k or 0 and sqrt(Sigma_kk); shifted log-normal:
#' exp(mu_k) 1{rc_mean} + exp(Sigma_kk / 2) and
#' sqrt((exp(Sigma_kk) - 1) exp(Sigma_kk))). The rest are unit-weighted
#' moments of the conditional means and variances; units with NA tastes are
#' dropped.
#'
#' @param object A `choicer_mxl` fit.
#' @param cmean,csd K_w x U conditional means and standard deviations.
#' @param w Length-U unit weights.
#' @returns A data.frame with one row per random coefficient.
#' @noRd
.taste_population <- function(object, cmean, csd, w) {
  K_w <- nrow(cmean)
  lognormal <- object$rc_dist == 1L
  shift <- rep(0, K_w)  # mu_final: mu (normal) or exp(mu) (log-normal)
  if (isTRUE(object$rc_mean)) {
    mu <- object$coefficients[object$param_map$mu]
    shift <- ifelse(lognormal, exp(mu), mu)
  }
  s2 <- diag(as.matrix(object$sigma))  # variance of (L eta)_k
  model_mean <- shift + ifelse(lognormal, exp(s2 / 2), 0)
  model_var <- ifelse(lognormal, (exp(s2) - 1) * exp(s2), s2)

  ok <- colSums(is.na(cmean) | is.na(csd)) == 0L
  wn <- w[ok] / sum(w[ok])
  m <- cmean[, ok, drop = FALSE]
  mean_cond <- drop(m %*% wn)
  var_cond <- drop((m - mean_cond)^2 %*% wn)             # Var_w(E[beta | y])
  var_resid <- drop(csd[, ok, drop = FALSE]^2 %*% wn)    # E_w[Var(beta | y)]

  data.frame(
    model_mean = model_mean,
    model_sd = sqrt(model_var),
    mean_cond = mean_cond,
    sd_cond = sqrt(var_cond),
    revealed = var_cond / model_var,
    residual = var_resid / model_var,
    row.names = rownames(cmean)
  )
}

#' Print conditional tastes
#'
#' @param x A \code{choicer_tastes} object from
#'   \code{\link{conditional_tastes}}.
#' @param digits Number of decimal places for the population table.
#' @param ... Additional arguments (ignored).
#' @returns The object invisibly.
#' @export
print.choicer_tastes <- function(x, digits = 3, ...) {
  cat("Conditional tastes E[beta | observed choices] for ", x$n_units, " ",
      x$unit, "s (S = ", x$S, " draws each)\n", sep = "")
  n_na <- sum(colSums(is.na(x$mean)) > 0L)
  cat("Unit-level means and SDs in $mean and $sd",
      if (n_na > 0L) paste0("; ", n_na, " ", x$unit, "(s) with NA tastes ",
                            "excluded below"),
      ".\n", sep = "")
  cat("Population distribution vs conditional means",
      " (revealed + residual should be close to 1):\n", sep = "")
  print(round(x$population, digits))
  invisible(x)
}
