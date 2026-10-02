
# Normalize a bounds vector to nloptr's full-length representation. Accepts
# NULL (returns the default), a full-length unnamed numeric (returned as-is
# after a length check), or a named partial numeric whose names must be a
# subset of `param_names`. Used by run_mxlogit() for both `lower` and `upper`.
.normalize_bound <- function(b, param_names, default, which) {
  if (is.null(b) || length(b) == 0L) {
    return(rep(default, length(param_names)))
  }
  if (!is.numeric(b)) stop("`", which, "` must be numeric or NULL.")
  if (anyNA(b)) stop("`", which, "` must not contain NA/NaN.")
  if (is.null(names(b))) {
    if (length(b) != length(param_names)) {
      stop("Unnamed `", which, "` must have length n_params (",
           length(param_names), "); got ", length(b), ".")
    }
    return(b)
  }
  dups <- unique(names(b)[duplicated(names(b))])
  if (length(dups)) {
    stop("Duplicate name(s) in `", which, "`: ",
         paste(dups, collapse = ", "), ".")
  }
  bad <- setdiff(names(b), param_names)
  if (length(bad)) {
    stop("Unknown parameter name(s) in `", which, "`: ",
         paste(bad, collapse = ", "),
         ". Valid names: ", paste(param_names, collapse = ", "), ".")
  }
  out <- rep(default, length(param_names))
  names(out) <- param_names
  out[names(b)] <- b
  unname(out)
}

# Labels of the random-coefficient covariance block, in the order of its
# Cholesky parameters: <prefix>_<i><j> over the lower triangle by rows, or
# the diagonal <prefix>_<k><k> without correlation. "L" names the Cholesky
# factor in the parameter vector; summary() shows "Sigma", the covariance
# L L' it implies.
.mxl_cov_names <- function(prefix, K_w, rc_correlation) {
  if (!rc_correlation) {
    return(paste0(prefix, "_", seq_len(K_w), seq_len(K_w)))
  }
  i <- rep(seq_len(K_w), seq_len(K_w))
  sprintf("%s_%d%d", prefix, i, sequence(seq_len(K_w)))
}

# Labels summary() prints for a mixed logit's parameters: the covariance
# block as Sigma_<i><j>, and the mean of a log-normal coefficient as
# exp(Mu_<variable>). run_mxlogit() checks them for repeats before the
# kernels have validated rc_dist, so only its first K_w entries are read.
.mxl_summary_labels <- function(param_names, param_map, K_w, rc_dist,
                                rc_correlation, rc_mean) {
  labels <- param_names
  if (!is.null(param_map$sigma)) {
    labels[param_map$sigma] <- .mxl_cov_names("Sigma", K_w, rc_correlation)
  }
  if (rc_mean && !is.null(param_map$mu)) {
    idx <- param_map$mu[rc_dist[seq_len(K_w)] %in% 1L]
    labels[idx] <- paste0("exp(", labels[idx], ")")
  }
  labels
}

#' Runs mixed logit estimation
#'
#' Estimates a mixed logit model via simulated maximum likelihood.
#'
#' Two workflows are supported:
#' \describe{
#'   \item{Convenience}{Supply \code{data} and column names. Data preparation
#'     (\code{\link{prepare_mxl_data}}) and Halton draw generation
#'     (\code{\link{get_halton_normals}}) are handled automatically.}
#'   \item{Advanced}{Call \code{\link{prepare_mxl_data}} and
#'     \code{\link{get_halton_normals}} yourself, then pass the results via
#'     \code{input_data} and \code{eta_draws}.}
#' }
#'
#' \strong{Cross-section vs panel.} With \code{person_col}, each decision
#' maker \eqn{n} draws one taste vector \eqn{\beta_n} from the mixing
#' distribution \eqn{f(\beta \mid \theta)} and keeps it in all of their choice
#' situations \eqn{t = 1, \dots, T_n}. The likelihood integrates the product
#' of their conditional logit probabilities (Revelt and Train 1998),
#' \deqn{L_n(\theta) = \int \prod_{t=1}^{T_n} P_{nt}(j_{nt} \mid \beta)
#'   f(\beta \mid \theta) \, d\beta,}
#' simulated with \code{S} draws per decision maker, so repeated choices
#' speak to each decision maker's tastes, and
#' \code{\link{conditional_tastes}} recovers what they reveal. Without
#' \code{person_col} (the default), each choice situation is an independent
#' draw from the mixing distribution: \eqn{L_i(\theta) = \int P_i(j_i \mid
#' \beta) f(\beta \mid \theta) \, d\beta}. That is the right model when tastes
#' are redrawn in every situation (or each decision maker chooses once). On
#' panel data with stable tastes it is a composite (pseudo-) likelihood: it
#' still targets the population mixing distribution, but it ignores the
#' within-person correlation of choices, so it is less efficient, its
#' standard errors must be clustered by decision maker (\code{cluster_col}),
#' and its log-likelihood, AIC and BIC are not comparable with those of the
#' panel fit.
#'
#' @param data Data frame containing choice data (convenience workflow).
#'   Mutually exclusive with \code{input_data}.
#' @param id_col Name of the column identifying choice situations.
#' @param alt_col Name of the column identifying alternatives.
#' @param choice_col Name of the column indicating chosen alternative (1/0).
#' @param covariate_cols Vector of column names for fixed covariates. None
#'   may take a name the model gives its parameters, or prints for them in
#'   \code{summary()}: \code{ASC_<label>}, \code{L_<i><j>} and
#'   \code{Sigma_<i><j>}, and with \code{rc_mean = TRUE}
#'   \code{Mu_<variable>} and, for a log-normal coefficient,
#'   \code{exp(Mu_<variable>)}.
#' @param random_var_cols Vector of column names for random coefficients.
#' @param input_data List output from \code{\link{prepare_mxl_data}} (advanced
#'   workflow). Mutually exclusive with \code{data}.
#' @param eta_draws Array of shape K_w x S x U with standard normal draws, one
#'   K_w x S block per likelihood unit: U is the number of decision makers
#'   when \code{input_data} was prepared with \code{person_col}
#'   (\code{length(input_data$Ti)}), and the number of choice situations
#'   otherwise. Required for the advanced workflow; auto-generated from
#'   \code{S} in the convenience workflow. Post-hoc methods
#'   (\code{vcov(fit, type = )}, \code{\link{conditional_tastes}}, prediction)
#'   regenerate Halton draws with \code{\link{get_halton_normals}} from the
#'   recorded draw count, so build \code{eta_draws} with it for them to
#'   reproduce the estimation draws.
#' @param S Integer number of Halton draws per decision maker, or per choice
#'   situation without \code{person_col} (convenience workflow only). Default
#'   100.
#' @param rc_dist Integer vector indicating distribution of random coefficients
#'   (0 = normal, 1 = log-normal). Default: all normal.
#' @param rc_mean Logical indicating whether to estimate means for random
#'   coefficients.
#' @param rc_correlation Logical indicating whether random coefficients are
#'   correlated (convenience workflow). Ignored when \code{input_data} is used
#'   (taken from the prepared data).
#' @param use_asc Logical indicating whether to include alternative-specific
#'   constants.
#' @param theta_init Initial parameter vector in natural-scale units. If
#'   \code{NULL}, defaults to zeros for the \eqn{\beta}, \eqn{\mu}, and ASC
#'   blocks, and \code{log(0.5)} on the Cholesky diagonal (so each diagonal
#'   factor \eqn{L_{pp} = 0.5}, i.e. a moderate random-coefficient variance of
#'   \code{0.25}). The zero-on-diagonal alternative corresponds to
#'   \eqn{L_{pp} = 1} (unit RC variance), which often lets the first L-BFGS
#'   step overshoot.
#' @param lower,upper Optional parameter bounds for the optimizer, in
#'   natural-scale units (forward-transformed internally to scaled space when
#'   \code{scale_vars != "none"}). Each accepts three forms:
#'   \describe{
#'     \item{\code{NULL}}{(default) Unbounded (\code{-Inf}/\code{Inf}).}
#'     \item{Unnamed numeric vector of length \code{n_params}}{Full-length
#'       vector ordered exactly like \code{theta_init} (the nloptr-native form).}
#'     \item{Named numeric vector}{Names must be a subset of the parameter
#'       names (\eqn{\beta} block: column names of \code{X};
#'       \eqn{\mu} block: \code{Mu_<col>} (if \code{rc_mean = TRUE});
#'       Cholesky block: \code{L_<i><j>} for \eqn{i \ge j}; ASC block:
#'       \code{ASC_<level>}). Unlisted parameters default to \eqn{\pm\infty}.
#'       This is the recommended form for typical use, e.g.
#'       \code{lower = c(L_11 = -5, L_22 = -5)} to clip Cholesky diagonals.}
#'   }
#' @param optimizer Optimizer to use: \code{"nloptr"} (default), \code{"optim"},
#'   or a custom function. See \code{\link{run_mnlogit}} for details.
#' @param control List of optimizer-specific control parameters.
#' @param se_method Method for computing standard errors. One of
#'   \code{"hessian"} (default) for the analytical Hessian of the simulated
#'   log-likelihood, \code{"bhhh"} for the BHHH/outer-product-of-gradients
#'   (OPG) estimator, \code{"sandwich"} for the robust (Huber-White)
#'   variance \eqn{V = A^{-1} B A^{-1}} (bread \eqn{A} = weighted negated
#'   Hessian, meat \eqn{B} = weight-squared OPG), or \code{"cluster"} for the
#'   cluster-robust sandwich (requires \code{cluster_col} or a prepared
#'   \code{input_data} with a \code{cluster} field). Use \code{"sandwich"} for
#'   valid inference under choice-based / WESML weighting, where the
#'   inverse-Hessian and ordinary BHHH are invalid; it reduces to the usual
#'   robust variance under uniform weights. BHHH scales better to large
#'   problems (many alternatives or simulation draws) but may underestimate
#'   standard errors in finite samples or away from the optimum. Any of these
#'   can also be recomputed post hoc via \code{vcov(fit, type = )}. Without
#'   \code{person_col}, clustering repairs the inference, not the likelihood:
#'   the cross-sectional likelihood still treats each choice situation as an
#'   independent draw from the mixing distribution, and on panel data is a
#'   less efficient composite likelihood for the same taste distribution.
#'   With \code{person_col} the likelihood unit is the decision maker, so
#'   scores are per decision maker and
#'   \code{"sandwich"} is already robust to dependence across a decision
#'   maker's choice situations; \code{"cluster"} is then needed only for
#'   coarser groups (e.g. households or markets) that nest decision makers.
#' @param cluster_col Optional name of a column in \code{data} holding cluster
#'   labels for cluster-robust standard errors (e.g. a person id when the same
#'   decision maker contributes several choice situations to a
#'   cross-sectional fit). Must be constant within each \code{id_col}, and
#'   within each decision maker when \code{person_col} is used (clusters must
#'   nest decision makers). Supplying \code{cluster_col} without an explicit
#'   \code{se_method} selects \code{se_method = "cluster"}.
#' @param scale_vars Pre-estimation column scaling for design matrices. One of
#'   \code{"none"} (default), \code{"sd"} (sample standard deviation),
#'   \code{"mad"} (\code{stats::mad}, i.e. 1.4826 \eqn{\times}
#'   median absolute deviation; SD-equivalent under normality), or
#'   \code{"iqr"} (\code{stats::IQR(x) / 1.349}; also SD-equivalent under
#'   normality). When not \code{"none"}, every column of \code{X} and \code{W}
#'   is divided by the chosen scale before optimization to improve Hessian
#'   conditioning. Robust scales (\code{"mad"}/\code{"iqr"}) better capture
#'   the bulk for heavy-tailed columns where SD is dominated by outliers, but
#'   \code{stats::mad} can return zero when more than half of a column's
#'   entries are identical (e.g., a sparse 0/1 dummy) and will then trigger
#'   the same near-constant-column error as \code{"sd"}. Coefficients and
#'   standard errors are back-transformed to the user's natural units via the
#'   delta method, so reported quantities are invariant to this choice.
#'   Columns of \code{W} associated with log-normal random coefficients
#'   (\code{rc_dist == 1}) are passed through unchanged, since the shifted
#'   log-normal parameterization does not admit a closed-form back-transform
#'   under multiplicative scaling.
#' @param weights Optional weight vector (convenience workflow), one weight per
#'   choice situation in ascending-id order (see
#'   \code{\link{prepare_mxl_data}}). With \code{person_col}, situations are
#'   reordered by decision maker, so a positional vector is accepted only when
#'   that leaves the id order unchanged; otherwise supply
#'   \code{weights_col}. If \code{NULL},
#'   equal weights are used. All weights must be finite and strictly positive.
#'   With \code{person_col} they are decision-maker weights, constant within
#'   each decision maker: the objective is \eqn{\sum_n w_n \log L_n}.
#' @param weights_col Optional name of a column in \code{data} holding a per-row
#'   weight (constant within each choice situation, finite and strictly positive).
#'   Mutually exclusive with
#'   \code{weights}; the recommended way to pass WESML weights from
#'   \code{\link{sample_by_choice}} / \code{\link{wesml_weights}}, since
#'   alignment is by id rather than by position. Convenience workflow only. If
#'   \code{data} carries choice-based-sampling provenance (a
#'   \code{"choice_sampling"} attribute, as attached by
#'   \code{\link{sample_by_choice}} / \code{\link{wesml_weights}}) and neither
#'   \code{weights} nor \code{weights_col} is supplied, the recorded weight
#'   column is auto-detected and applied (with a message); if that column is
#'   absent the call errors rather than silently fitting an unweighted model
#'   under a WESML label. With \code{person_col} the weight must also be
#'   constant within each decision maker, and WESML is not supported:
#'   choice-based weights vary with the chosen alternative across a decision
#'   maker's situations, so they cannot weight the panel likelihood.
#' @param outside_opt_label Label for the outside option (convenience workflow).
#' @param include_outside_option Logical whether to include an outside option
#'   (convenience workflow).
#' @param draws Draw storage mode. One of \code{"store"} (default) or \code{"generate"}.
#'   \code{"store"} pre-materializes the full \eqn{K_w \times S \times U} Halton cube, one
#'   block per likelihood unit (\code{U} decision makers with \code{person_col}, choice
#'   situations otherwise; existing behavior, exact reproducibility). It supports at
#'   most \eqn{2^{31} - 1} points (\eqn{S \times U}) and \eqn{2^{32} - 1} values
#'   (\eqn{K_w \times S \times U}); \code{predict()}, \code{elasticities()},
#'   \code{diversion_ratios()}, \code{blp()}, \code{logsum()},
#'   \code{consumer_surplus()} and \code{gof()} (hence \code{summary()}'s fit
#'   statistics) regenerate one block per choice situation, so for them the
#'   number of choice situations replaces \eqn{U} (see
#'   \code{\link{get_halton_normals}}). \code{"generate"}
#'   computes each unit's draws on-the-fly in C++ from a stored seed, eliminating the O(U)
#'   cube; recommended for memory-constrained or large-N settings. With
#'   \code{scramble = "permuted"}, each base-\eqn{b} digit position in each
#'   dimension receives a seeded permutation shared across sequence indices.
#'   This is not Owen's nested-uniform scramble and does not carry standard
#'   randomized-QMC unbiasedness or error-estimation guarantees. Only supported
#'   in the convenience workflow.
#' @param seed Integer master seed for the on-the-fly generator. Used only when
#'   \code{draws = "generate"}. If \code{NULL} (default), a seed is drawn from R's
#'   RNG at call time (so \code{set.seed()} governs reproducibility). Ignored when
#'   \code{draws = "store"}.
#' @param scramble Scrambling mode for on-the-fly Halton draws. One of
#'   \code{"permuted"} (default) for seeded position-wise digit permutations or
#'   \code{"none"} for plain Halton (identity permutations). The historical value
#'   \code{"owen"} is accepted with a deprecation warning as an alias for
#'   \code{"permuted"}; the implementation is not Owen's nested-uniform scramble.
#'   \code{"none"} reproduces the randtoolbox sequence exactly. Simulation-draw
#'   sensitivity should be assessed by increasing \code{S} and, for
#'   \code{"permuted"}, varying \code{seed}. Used only when
#'   \code{draws = "generate"}.
#' @param keep_data Logical. If \code{TRUE} (default), stores prepared data in
#'   the returned object for post-estimation functions.
#' @param nloptr_opts Deprecated. Use \code{optimizer} and \code{control}
#'   instead.
#' @param person_col Optional name of the column in \code{data} identifying
#'   decision makers (respondents). When supplied, all choice situations of a
#'   decision maker share one draw of the random coefficients (panel
#'   likelihood); \code{id_col} must still identify choice situations uniquely
#'   across the data set. \code{NULL} (default) makes each choice situation its
#'   own decision maker (cross-sectional likelihood). Convenience workflow
#'   only; in the advanced workflow pass \code{person_col} to
#'   \code{\link{prepare_mxl_data}}.
#' @returns A \code{choicer_mxl} object (inherits from \code{choicer_fit}).
#'   Standard S3 methods available: \code{summary()}, \code{coef()},
#'   \code{vcov()}, \code{logLik()}, \code{AIC()}, \code{BIC()},
#'   \code{nobs()}. A panel fit (\code{person_col}) also carries
#'   \code{n_persons}, the number of decision makers; \code{nobs()} remains
#'   the number of choice situations, which is also the \eqn{n} in BIC's
#'   \eqn{\log n} penalty.
#' @references Revelt, D. and Train, K. (1998). Mixed logit with repeated
#'   choices: households' choices of appliance efficiency level.
#'   \emph{Review of Economics and Statistics} 80(4), 647-657.
#' @seealso \code{\link{conditional_tastes}}, \code{\link{run_hmnlogit}} (the
#'   hierarchical Bayes counterpart)
#' @inheritSection prepare_mnl_data Column names
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(42)
#' N <- 100; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#'
#' fit <- run_mxlogit(
#'   data = dt, id_col = "id", alt_col = "alt", choice_col = "choice",
#'   covariate_cols = "x1", random_var_cols = c("w1", "w2"), S = 50L
#' )
#' summary(fit)
#' }
#' @importFrom nloptr nloptr
#' @export
run_mxlogit <- function(
    data = NULL,
    id_col = NULL,
    alt_col = NULL,
    choice_col = NULL,
    covariate_cols = NULL,
    random_var_cols = NULL,
    input_data = NULL,
    eta_draws = NULL,
    S = 100L,
    rc_dist = NULL,
    rc_mean = FALSE,
    rc_correlation = FALSE,
    use_asc = TRUE,
    theta_init = NULL,
    lower = NULL,
    upper = NULL,
    optimizer = NULL,
    control = list(),
    se_method = c("hessian", "bhhh", "sandwich", "cluster"),
    scale_vars = c("none", "sd", "mad", "iqr"),
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    draws       = c("store", "generate"),
    seed        = NULL,
    scramble    = c("permuted", "none", "owen"),
    keep_data = TRUE,
    nloptr_opts = NULL,
    weights_col = NULL,
    cluster_col = NULL,
    person_col = NULL
) {
  cl <- match.call()

  se_method_default <- missing(se_method)
  se_method <- match.arg(se_method)
  if (!is.null(cluster_col) && se_method_default) se_method <- "cluster"
  scale_vars <- match.arg(scale_vars)
  draws   <- match.arg(draws)
  scramble <- match.arg(scramble)
  if (identical(scramble, "owen")) {
    warning("scramble = \"owen\" is a deprecated alias for ",
            "scramble = \"permuted\". The implemented position-wise digit ",
            "permutation is not Owen's nested-uniform scramble.",
            call. = FALSE)
    scramble <- "permuted"
  }

  # Validate seed parameter
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0) {
      stop("'seed' must be NULL or a single non-negative integer.")
    }
    seed <- as.integer(seed)
    if (draws != "generate") {
      message("'seed' is ignored when draws = 'store'.")
    }
  }

  # Backward compatibility: nloptr_opts -> optimizer + control
  if (!is.null(nloptr_opts)) {
    message("'nloptr_opts' is deprecated. Use 'optimizer' and 'control' instead.")
    optimizer <- optimizer %||% "nloptr"
    control <- nloptr_opts
  }

  # --- Resolve input pathway --------------------------------------------------
  has_data <- !is.null(data)
  has_input <- !is.null(input_data)
  cs_meta <- if (has_data) attr(data, "choice_sampling") else attr(input_data, "choice_sampling")

  if (has_data && has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced), not both.")
  }
  if (!has_data && !has_input) {
    stop("Supply either 'data' (convenience) or 'input_data' (advanced).")
  }
  if (has_input && !is.null(weights_col)) {
    stop("`weights_col` is only supported in the convenience (data) workflow. ",
         "Bake weights into `input_data` via prepare_mxl_data(weights_col = ) ",
         "or supply `weights` to prepare_mxl_data().")
  }
  if (has_input && !is.null(cluster_col)) {
    stop("`cluster_col` is only supported in the convenience (data) workflow. ",
         "Bake cluster labels into `input_data` via ",
         "prepare_mxl_data(cluster_col = ).")
  }
  if (has_input && !is.null(person_col)) {
    stop("`person_col` is only supported in the convenience (data) workflow. ",
         "Bake the panel structure into `input_data` via ",
         "prepare_mxl_data(person_col = ).")
  }

  if (has_data) {
    # Convenience workflow: validate required column-name arguments
    if (is.null(id_col) || is.null(alt_col) || is.null(choice_col) ||
        is.null(covariate_cols) || is.null(random_var_cols)) {
      stop("Convenience workflow requires: id_col, alt_col, choice_col, ",
           "covariate_cols, and random_var_cols.")
    }
    # WESML provenance present but no weights supplied: auto-adopt the recorded
    # weight column, or error -- never silently fit unweighted under a WESML label.
    # (With person_col, prepare_mxl_data() rejects WESML provenance outright.)
    if (has_data && !is.null(cs_meta) && is.null(person_col) &&
        is.null(weights) && is.null(weights_col)) {
      wn <- cs_meta$weight_name
      if (!is.null(wn) && wn %in% names(data)) {
        weights_col <- wn
        message("Detected WESML choice-based-sampling provenance; applying attached ",
                "weights from column '", wn, "'.")
      } else {
        stop("Data carries WESML choice-based-sampling provenance but no weights were ",
             "supplied, and the recorded weight column (",
             if (is.null(wn)) "unknown" else paste0("'", wn, "'"),
             ") is not present in `data`. Pass `weights_col=` or `weights=` explicitly.",
             call. = FALSE)
      }
    }
    input_data <- prepare_mxl_data(
      data = data,
      id_col = id_col,
      alt_col = alt_col,
      choice_col = choice_col,
      covariate_cols = covariate_cols,
      random_var_cols = random_var_cols,
      weights = weights,
      weights_col = weights_col,
      outside_opt_label = outside_opt_label,
      include_outside_option = include_outside_option,
      rc_correlation = rc_correlation,
      cluster_col = cluster_col,
      person_col = person_col
    )
  } else {
    # Advanced workflow
    if (draws == "generate") {
      stop("draws=generate is only supported in the convenience workflow. ",
           "In the advanced workflow, supply 'eta_draws' directly.")
    }
    if (is.null(eta_draws)) {
      stop("'eta_draws' is required when using 'input_data' (advanced workflow).")
    }
    # `S` is a convenience-workflow argument; record the draws actually used so
    # post-hoc methods regenerate the same Halton blocks.
    S <- dim(eta_draws)[2L]
  }

  if (se_method == "cluster" && is.null(input_data$cluster)) {
    stop("se_method = \"cluster\" needs cluster labels: pass `cluster_col=` ",
         "(convenience workflow) or prepare `input_data` with ",
         "prepare_mxl_data(cluster_col = ).", call. = FALSE)
  }

  # Parameter dimensions
  J <- nrow(input_data$alt_mapping)
  K_x <- ncol(input_data$X)
  K_w <- ncol(input_data$W)
  rc_correlation <- input_data$rc_correlation
  L_size <- if (rc_correlation) K_w * (K_w + 1) / 2 else K_w
  mu_size <- if (rc_mean) K_w else 0
  n_asc <- J - 1
  n_params <- K_x + mu_size + L_size + n_asc

  if (is.null(rc_dist)) rc_dist <- rep(0L, K_w)

  # Parameter index map (built early so the scaling layer can address blocks)
  pos <- 0
  param_map <- list(beta = seq_len(K_x))
  pos <- K_x
  if (mu_size > 0) {
    param_map$mu <- pos + seq_len(mu_size)
    pos <- pos + mu_size
  }
  param_map$sigma <- pos + seq_len(L_size)
  pos <- pos + L_size
  param_map$asc <- pos + seq_len(n_asc)

  # Parameter names (built early; reused for theta_hat, vcov, se downstream)
  beta_names <- colnames(input_data$X)
  mu_names <- if (rc_mean) paste0("Mu_", colnames(input_data$W)) else character(0)
  sigma_names <- .mxl_cov_names("L", K_w, rc_correlation)
  alt_col <- names(input_data$alt_mapping)[2]
  asc_names <- paste0("ASC_", input_data$alt_mapping[[alt_col]][2:J])
  param_names <- c(beta_names, mu_names, sigma_names, asc_names)
  # summary() relabels some parameters: its labels may not repeat either
  .check_param_names(
    list(param_names,
         .mxl_summary_labels(param_names, param_map, K_w, rc_dist,
                             rc_correlation, rc_mean)),
    c("ASC_<label>", "Mu_<variable>", "L_<i><j>", "Sigma_<i><j>",
      "exp(Mu_<variable>)")
  )

  # Draws for the convenience workflow, built once the inputs have passed
  # the checks above: a store-mode cube can take gigabytes.
  if (has_data) {
    if (draws == "store") {
      # One K_w x S draw block per likelihood unit: per decision maker with
      # person_col, per choice situation otherwise.
      n_units <- length(.unit_first(input_data))
      eta_draws <- get_halton_normals(S, n_units, K_w)
    } else {
      # generate mode: no cube ever materialized; empty placeholder
      eta_draws <- array(0, dim = c(K_w, 0L, 0L))
      # Draw seed from R RNG when not supplied (like run_mnprobit)
      if (is.null(seed)) {
        seed <- sample.int(.Machine$integer.max, 1L)
      }
    }
  }

  # --- Variable scaling (optional) --------------------------------------------
  # Scale columns of X and W by their sample SD to improve Hessian conditioning.
  # Keep the natural-scale matrices for storage; theta_init is interpreted in
  # natural units and forward-transformed below; theta_hat and vcov are
  # back-transformed after optimization so reported quantities are in the
  # user's natural units. sX and sW are returned as 1s when scale_vars="none".
  natural_X <- input_data$X
  natural_W <- input_data$W
  sX <- rep(1, K_x); names(sX) <- colnames(input_data$X)
  sW <- rep(1, K_w); names(sW) <- colnames(input_data$W)
  if (scale_vars != "none") {
    if (K_x > 0) {
      sX_raw <- .column_scales(input_data$X, scale_vars)
      .assert_scales_ok(sX_raw, scale_vars, "fixed-coefficient")
      sX <- sX_raw
      input_data$X <- sweep(input_data$X, 2, sX, "/")
    }
    if (K_w > 0) {
      sW_raw <- .column_scales(input_data$W, scale_vars)
      normal_cols <- which(rc_dist == 0L)
      if (length(normal_cols) > 0L) {
        .assert_scales_ok(sW_raw, scale_vars, "normal random-coefficient",
                          idx = normal_cols)
      }
      # Preserve names from sW_raw; carve out log-normal columns (pass-through).
      sW <- sW_raw
      sW[rc_dist == 1L] <- 1
      input_data$W <- sweep(input_data$W, 2, sW, "/")
      n_lognormal <- sum(rc_dist == 1L)
      if (K_w > 0L && n_lognormal == K_w) {
        message("scale_vars='", scale_vars,
                "': all random-coefficient column(s) are log-normal; W not scaled.")
      } else if (n_lognormal > 0L) {
        message("scale_vars='", scale_vars,
                "': passing through log-normal random-coefficient column(s) ",
                "unchanged (no closed-form back-transform).")
      }
    }
  }

  # --- Natural <-> scaled Jacobian --------------------------------------------
  # Maps scaled-space parameters back to natural-scale units:
  #   theta_natural = bt_mult * theta_scaled + bt_shift
  # Inverse forward-transforms theta_init from natural to scaled space.
  # ASCs and any unset entries default to identity (mult=1, shift=0).
  bt_mult <- rep(1, n_params)
  bt_shift <- rep(0, n_params)
  if (scale_vars != "none") {
    if (K_x > 0) bt_mult[param_map$beta] <- 1 / sX
    if (mu_size > 0) bt_mult[param_map$mu] <- 1 / sW
    if (rc_correlation) {
      idx <- 1L
      for (i in seq_len(K_w)) {
        for (j in seq_len(i)) {
          pos <- param_map$sigma[idx]
          if (i == j) {
            bt_shift[pos] <- -log(sW[i])
          } else {
            bt_mult[pos] <- 1 / sW[i]
          }
          idx <- idx + 1L
        }
      }
    } else {
      for (i in seq_len(K_w)) {
        pos <- param_map$sigma[i]
        bt_shift[pos] <- -log(sW[i])
      }
    }
  }

  # Resolve theta_init (natural units); forward-transform to scaled space.
  # Default cold-start: zero on every block except the Cholesky diagonal,
  # which sits at log(0.5) so each diagonal factor L_pp = 0.5 (RC variance
  # 0.25). Starting at log(1) = 0 corresponds to L_pp = 1 (unit RC variance),
  # which is often too large for typical specs and lets the first L-BFGS step
  # push ell_pp far enough that L_pp underflows / overflows.
  if (is.null(theta_init)) {
    theta_init <- rep(0, n_params)
    if (K_w > 0L) {
      if (rc_correlation) {
        diag_idx <- param_map$sigma[cumsum(seq_len(K_w))]
      } else {
        diag_idx <- param_map$sigma
      }
      theta_init[diag_idx] <- log(0.5)
    }
  }
  if (scale_vars != "none") {
    theta_init <- (theta_init - bt_shift) / bt_mult
  }

  # Normalize lower/upper bounds (natural units in, scaled units out).
  lower <- .normalize_bound(lower, param_names, -Inf, "lower")
  upper <- .normalize_bound(upper, param_names,  Inf, "upper")
  if (scale_vars != "none") {
    lower <- (lower - bt_shift) / bt_mult
    upper <- (upper - bt_shift) / bt_mult
  }

  # Resolve generate-mode parameters for C++ kernels.
  # In store mode (draws="store"): gen_seed_cpp = -1L triggers cube path (unchanged behavior).
  gen_seed_cpp     <- if (draws == "generate") seed else -1L
  gen_scramble_cpp <- if (draws == "generate") (if (scramble == "permuted") 1L else 0L) else 1L
  gen_S_cpp        <- if (draws == "generate") S else 0L

  # Build eval_f closure (the kernel's overflow sentinel kept above the path)
  eval_f <- .lift_sentinel(function(theta) {
    mxl_loglik_gradient_parallel(
      theta = theta,
      X = input_data$X,
      W = input_data$W,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      M = input_data$M,
      weights = input_data$weights,
      rc_dist = rc_dist,
      rc_correlation = rc_correlation,
      rc_mean = rc_mean,
      eta_draws = eta_draws,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option,
      gen_seed = gen_seed_cpp,
      gen_scramble = gen_scramble_cpp,
      gen_S = gen_S_cpp,
      Ti = input_data$Ti
    )
  })

  # Run optimizer
  elapsed <- system.time({
    opt <- run_optimizer(
      optimizer = optimizer,
      theta_init = theta_init,
      eval_f = eval_f,
      lower = lower,
      upper = upper,
      control = control
    )
  })

  message("Optimization run time ", convertTime(elapsed))

  # Estimate at the optimum (in scaled space if scale_vars='sd')
  theta_hat <- opt$par
  names(theta_hat) <- param_names

  # Choice-based-sampling provenance and a guardrail for weighted inference.
  choice_sampling <- .resolve_choice_sampling(
    weights = input_data$weights, se_method = se_method,
    cs_meta = cs_meta, has_input = has_input, prepare_fn = "prepare_mxl_data"
  )

  # Compute vcov eagerly using the selected SE method.
  # For "sandwich" (robust / WESML) standard errors, form V = A^{-1} B A^{-1}
  # with bread A = weighted negated Hessian and meat B = weight-squared OPG
  # (pass weights^2 to the BHHH routine, whose per-unit score is
  # weight-free). For "cluster", the meat is the outer product of
  # within-cluster sums of weighted scores. Scores, weights and cluster labels
  # are per likelihood unit (decision maker with person_col). Computed in
  # scaled space; the back-transform below applies.
  if (se_method %in% c("sandwich", "cluster")) {
    A_bread <- mxl_hessian_parallel(
      theta = theta_hat, X = input_data$X, W = input_data$W,
      alt_idx = input_data$alt_idx, choice_idx = input_data$choice_idx,
      M = input_data$M, weights = input_data$weights, eta_draws = eta_draws,
      rc_dist = rc_dist, rc_correlation = rc_correlation, rc_mean = rc_mean,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option,
      gen_seed = gen_seed_cpp, gen_scramble = gen_scramble_cpp, gen_S = gen_S_cpp,
      Ti = input_data$Ti
    )
    B_meat <- if (se_method == "sandwich") {
      mxl_bhhh_parallel(
        theta = theta_hat, X = input_data$X, W = input_data$W,
        alt_idx = input_data$alt_idx, choice_idx = input_data$choice_idx,
        M = input_data$M, weights = input_data$weights^2, eta_draws = eta_draws,
        rc_dist = rc_dist, rc_correlation = rc_correlation, rc_mean = rc_mean,
        use_asc = use_asc,
        include_outside_option = input_data$include_outside_option,
        gen_seed = gen_seed_cpp, gen_scramble = gen_scramble_cpp, gen_S = gen_S_cpp,
        Ti = input_data$Ti
      )
    } else {
      S_scores <- mxl_scores_parallel(
        theta = theta_hat, X = input_data$X, W = input_data$W,
        alt_idx = input_data$alt_idx, choice_idx = input_data$choice_idx,
        M = input_data$M, eta_draws = eta_draws,
        rc_dist = rc_dist, rc_correlation = rc_correlation, rc_mean = rc_mean,
        use_asc = use_asc,
        include_outside_option = input_data$include_outside_option,
        gen_seed = gen_seed_cpp, gen_scramble = gen_scramble_cpp, gen_S = gen_S_cpp,
        Ti = input_data$Ti
      )
      u <- .unit_first(input_data)
      .score_meat(S_scores, input_data$weights[u], "cluster",
                  .to_units(input_data$cluster, input_data, "`cluster_col`"))
    }
    vcov_result <- .sandwich_combine(A_bread, B_meat)
  } else {
  hess <- switch(
    se_method,
    hessian = mxl_hessian_parallel(
      theta = theta_hat,
      X = input_data$X,
      W = input_data$W,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      M = input_data$M,
      weights = input_data$weights,
      eta_draws = eta_draws,
      rc_dist = rc_dist,
      rc_correlation = rc_correlation,
      rc_mean = rc_mean,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option,
      gen_seed = gen_seed_cpp, gen_scramble = gen_scramble_cpp, gen_S = gen_S_cpp,
      Ti = input_data$Ti
    ),
    bhhh = mxl_bhhh_parallel(
      theta = theta_hat,
      X = input_data$X,
      W = input_data$W,
      alt_idx = input_data$alt_idx,
      choice_idx = input_data$choice_idx,
      M = input_data$M,
      weights = input_data$weights,
      eta_draws = eta_draws,
      rc_dist = rc_dist,
      rc_correlation = rc_correlation,
      rc_mean = rc_mean,
      use_asc = use_asc,
      include_outside_option = input_data$include_outside_option,
      gen_seed = gen_seed_cpp, gen_scramble = gen_scramble_cpp, gen_S = gen_S_cpp,
      Ti = input_data$Ti
    )
  )
  vcov_result <- invert_hessian(hess)
  }
  if (!is.null(vcov_result$vcov)) {
    rownames(vcov_result$vcov) <- param_names
    colnames(vcov_result$vcov) <- param_names
    names(vcov_result$se) <- param_names
  }

  # --- Back-transform to natural scale ----------------------------------------
  # Uses the bt_mult / bt_shift map built before optimization:
  #   theta_natural = bt_mult * theta_scaled + bt_shift
  #   vcov_natural  = (bt_mult bt_mult') o vcov_scaled  (shifts don't enter)
  if (scale_vars != "none") {
    bt <- .backtransform_estimates(theta_hat, vcov_result, bt_mult, bt_shift, param_names)
    theta_hat <- bt$theta
    vcov_result <- bt$vcov_result
    input_data$X <- natural_X
    input_data$W <- natural_W
  }

  # Reconstruct Sigma for display (from back-transformed L params if scaled)
  L_params <- theta_hat[param_map$sigma]
  sigma_mat <- build_var_mat(L_params, K_w, rc_correlation)
  w_names <- colnames(input_data$W)
  if (!is.null(w_names)) {
    rownames(sigma_mat) <- w_names
    colnames(sigma_mat) <- w_names
  }

  # Draws info (metadata only, not the full array). N is the number of choice
  # situations, the draw-block count of the per-situation prediction sites;
  # estimation-type sites use one block per likelihood unit (.unit_first()).
  draws_info <- list(
    S       = S,
    N       = input_data$N,
    K_w     = K_w,
    mode    = draws,
    seed    = if (draws == "generate") seed else NULL,
    scramble = if (draws == "generate") scramble else NULL
  )

  # Build S3 object
  new_choicer_mxl(
    call = cl,
    coefficients = theta_hat,
    loglik = -opt$value,
    nobs = input_data$N,
    n_params = n_params,
    convergence = opt$convergence,
    message = opt$message,
    data_spec = input_data$data_spec,
    alt_mapping = input_data$alt_mapping,
    param_map = param_map,
    use_asc = use_asc,
    include_outside_option = input_data$include_outside_option,
    optimizer = list(
      name = if (is.function(optimizer)) "custom" else (optimizer %||% "nloptr"),
      control = control,
      elapsed_time = elapsed[["elapsed"]],
      iterations = opt$iterations
    ),
    vcov = vcov_result$vcov,
    se = vcov_result$se,
    data = if (keep_data) {
      list(
        X = input_data$X,
        W = input_data$W,
        alt_idx = input_data$alt_idx,
        choice_idx = input_data$choice_idx,
        M = input_data$M,
        weights = input_data$weights,
        cluster = input_data$cluster,
        situation_ids = input_data$situation_ids,
        Ti = input_data$Ti,
        person_ids = input_data$person_ids
      )
    },
    draws_info = draws_info,
    rc_dist = rc_dist,
    rc_correlation = rc_correlation,
    rc_mean = rc_mean,
    sigma = sigma_mat,
    se_method = se_method,
    scale_vars = scale_vars,
    sX = sX,
    sW = sW,
    choice_sampling = choice_sampling,
    n_persons = if (!is.null(input_data$Ti)) length(input_data$Ti)
  )
}


#' Prepare inputs for mixed logit estimation
#'
#' Prepares and validates inputs for mixed logit estimation routine.
#'
#' Rows are ordered by choice-situation id and, within a situation, by
#' alternative. With \code{person_col}, rows are ordered by decision maker
#' first, then by situation id and alternative, so that each decision maker's
#' situations are contiguous. A positional \code{weights} vector is read in
#' the prepared order, which is ascending id in the cross-section; with
#' \code{person_col} it is accepted only when the prepared order is also
#' ascending id, and an error otherwise points to \code{weights_col}, which is
#' aligned by id and is the safer interface.
#'
#' @param data Data frame containing choice data
#' @param id_col Name of the column identifying choice situations
#' @param alt_col Name of the column identifying alternatives
#' @param choice_col Name of the column indicating chosen alternative (1 = chosen, 0 = not chosen)
#' @param covariate_cols Vector of names of columns to be used as covariates
#' @param random_var_cols Vector of names of columns to be used as random variables
#' @param weights Optional vector of weights, one per choice situation in the
#'   prepared order (ascending id; see Details). If NULL, equal weights are
#'   used. All weights
#'   must be finite and strictly positive. With \code{person_col} they are
#'   decision-maker weights and must be constant within each decision maker.
#' @param weights_col Optional name of a column in \code{data} holding a
#'   per-row weight (constant within each choice situation, and within each
#'   decision maker when \code{person_col} is used; finite and strictly
#'   positive). Mutually exclusive with \code{weights}.
#' @param outside_opt_label Label for the outside option (if any). If NULL, no outside option is assumed.
#' @param include_outside_option Logical indicating whether to include an outside option in the model.
#' @param rc_correlation Logical indicating whether random coefficients are correlated. Default is FALSE.
#' @param cluster_col Optional name of a column in \code{data} holding cluster
#'   labels for cluster-robust standard errors. Must be constant within each
#'   \code{id_col}, and within each decision maker when \code{person_col} is
#'   used (clusters must nest decision makers); collapsed to one label per
#'   choice situation and returned as \code{cluster}.
#' @param person_col Optional name of the column identifying decision makers
#'   (respondents). When supplied, all choice situations of a decision maker
#'   share one draw of the random coefficients (panel likelihood);
#'   \code{id_col} must still identify choice situations uniquely across the
#'   data set. \code{NULL} (default) makes each choice situation its own
#'   decision maker (cross-sectional likelihood). Choice-based (WESML) samples
#'   are not supported with \code{person_col}: their weights vary with the
#'   chosen alternative, not by decision maker.
#' @returns A `choicer_data_mxl` object (list) containing:
#'   \itemize{
#'     \item `X`: Fixed-coefficient design matrix (sum(M) x K_x, double).
#'     \item `W`: Random-coefficient design matrix (sum(M) x K_w, double).
#'     \item `alt_idx`: Integer vector of alternative indices.
#'     \item `choice_idx`: Integer vector of chosen alternative indices.
#'     \item `M`: Integer vector with number of alternatives per choice situation.
#'     \item `N`: Number of choice situations.
#'     \item `weights`: Vector of weights.
#'     \item `cluster`: Vector of cluster labels (or `NULL`).
#'     \item `situation_ids`: Choice-situation ids in prepared (sorted) order.
#'     \item `Ti`: Number of choice situations of each decision maker, in
#'       prepared order (`NULL` without `person_col`).
#'     \item `person_ids`: Decision-maker ids in prepared order (`NULL`
#'       without `person_col`).
#'     \item `include_outside_option`: Logical flag.
#'     \item `rc_correlation`: Logical flag.
#'     \item `alt_mapping`: data.table mapping alternatives to summary statistics.
#'     \item `dropped_cols`: Names of columns dropped due to collinearity, if any.
#'     \item `data_spec`: List with column-name metadata (incl. `person_col`).
#'   }
#' @inheritSection prepare_mnl_data Column names
#' @examples
#' library(data.table)
#' set.seed(42)
#' N <- 50; J <- 3
#' dt <- data.table(id = rep(1:N, each = J), alt = rep(1:J, N))
#' dt[, `:=`(x1 = rnorm(.N), w1 = rnorm(.N), w2 = rnorm(.N))]
#' dt[, choice := 0L]
#' dt[, choice := sample(c(1L, rep(0L, J - 1))), by = id]
#' input <- prepare_mxl_data(dt, "id", "alt", "choice", "x1", c("w1", "w2"))
#' str(input$X)
#' str(input$W)
#' @export
prepare_mxl_data <- function(
    data,
    id_col,
    alt_col,
    choice_col,
    covariate_cols,
    random_var_cols,
    weights = NULL,
    outside_opt_label = NULL,
    include_outside_option = FALSE,
    rc_correlation = FALSE,
    weights_col = NULL,
    cluster_col = NULL,
    person_col = NULL
) {

  ## Preliminary housekeeping --------------------------------------------------
  # Capture any choice-based-sampling provenance before column drops / coercion,
  # so it can be carried onto the returned object for the advanced pathway.
  cs_provenance <- attr(data, "choice_sampling")
  if (!is.null(person_col) && !is.null(cs_provenance)) {
    stop("`person_col` cannot be used with a choice-based (WESML) sample: ",
         "WESML weights are defined per choice situation (they vary with the ",
         "chosen alternative), so they cannot weight the panel likelihood ",
         "sum_n w_n log L_n, which has one weight per decision maker. Fit the ",
         "cross-sectional likelihood with the WESML weights instead, clustering ",
         "on the decision maker via `cluster_col`.", call. = FALSE)
  }

  # Check if all relevant variables are available
  needed <- c(id_col, alt_col, choice_col, covariate_cols, random_var_cols)
  if (!is.null(weights) && !is.null(weights_col)) {
    stop("Supply only one of `weights` or `weights_col`.")
  }
  if (!is.null(weights_col)) needed <- c(needed, weights_col)
  if (!is.null(cluster_col)) needed <- c(needed, cluster_col)
  if (!is.null(person_col)) needed <- c(needed, person_col)
  # Names first: they need no data, and reading the data can stop for a
  # missing bit64.
  .check_col_names(needed, alt_col)
  # The covariates are read in place from `src` (the caller's data, never
  # modified, when possible); only the index columns are copied into `dt`,
  # and X and W are gathered by source row at the end.
  prep_src <- .prep_source(data, needed)
  src <- prep_src$src
  if (!all(needed %in% names(src)))
    stop("Missing columns: ",
         paste(setdiff(needed, names(src)), collapse = ", "))
  dt <- .copy_cols(src, c(id_col, alt_col, choice_col, weights_col,
                          cluster_col, person_col))
  dt[, .choicer_row := seq_len(.N)]
  .int64_cols_to_double(dt, weights_col)

  ## Drop ids with missing observations ----------------------------------------
  has_na <- .rows_with_na(src, prep_src$scan)
  ids_to_drop <- unique(dt[[id_col]][has_na[dt$.choicer_row]])
  if (length(ids_to_drop) > 0) {
    # computed outside dt[...], where a column named ids_to_drop would mask it
    keep <- !(dt[[id_col]] %in% ids_to_drop)
    dt <- dt[keep]
    warning("Removed ", length(ids_to_drop),
            " choice situations containing missing values.")
  }
  if (nrow(dt) == 0) {
    stop("All choice situations removed due to missing values.")
  }

  ## Sanity checks ---------------------------------------------------------

  ## covariates must be numeric
  if (!all(vapply(match(covariate_cols, names(src)),
                  function(j) is.numeric(.subset2(src, j)), NA)))
    stop("All covariates must be numeric.")
  if (!all(vapply(match(random_var_cols, names(src)),
                  function(j) is.numeric(.subset2(src, j)), NA)))
    stop("All covariates must be numeric.")

  ## choice column must be 0 or 1
  if (!.is_zero_one(dt[[choice_col]]))
    stop("`", choice_col, "` must contain only 0 and 1.")

  ## Panel: choice situations nest in decision makers. Checked before the
  ## per-id choice counts, which would misfire if ids restart within persons.
  if (!is.null(person_col)) {
    n_persons_per_id <- .n_distinct_by(dt, person_col, id_col)
    if (any(n_persons_per_id != 1L)) {
      stop("Each '", id_col, "' must belong to exactly one '", person_col,
           "': `id_col` must identify choice situations uniquely across ",
           "decision makers, but ", sum(n_persons_per_id != 1L), " id(s) ",
           "appear under several. If situation ids restart within each ",
           "decision maker, build a unique id, e.g. paste(", person_col, ", ",
           id_col, ").", call. = FALSE)
    }
  }

  ## Exactly one '1' per choice situation
  n_chosen <- .n_chosen_by(dt, choice_col, id_col)
  if (include_outside_option == FALSE && any(n_chosen != 1)) {
    stop("Each ", id_col, " must have exactly one chosen alternative (one '1' in ",
         choice_col, ").")
  }
  if (include_outside_option && any(n_chosen > 1)) {
    stop("Each ", id_col, " must have at most one chosen alternative (one '1' in ",
         choice_col, "). An id with no explicit choice is assumed to be outside option.")
  }

  ## Create integer alternative codes ------------------------------------------

  if (!is.null(outside_opt_label) && include_outside_option==FALSE) {
    levels <- c(outside_opt_label, sort(setdiff(unique(dt[[alt_col]]), outside_opt_label)))
  } else {
    levels <- sort(unique(dt[[alt_col]]))
  }

  dt <- .code_alternatives(dt, alt_col, levels)

  ## Order rows ----------------------------------------------------------------
  ##   within each id: ascending alternative id
  ##   between ids   : ascending id
  ##   with person_col, first by decision maker, so that each decision maker's
  ##   situations are contiguous (the kernels' unit layout)
  data.table::setorderv(dt, c(person_col, id_col, ".choicer_alt_int"))

  ## index of each row within its choice set (computed outside dt[...], where
  ## a column named id_col would take the local's place)
  .choicer_idx <- data.table::rowidv(dt, cols = id_col)
  dt[, .choicer_idx_in_group := .choicer_idx]
  rm(.choicer_idx)

  ## Build objects -------------------------------------------------------------
  ## design matrix
  .check_design_size(nrow(dt), random_var_cols,
                     "The random-coefficient design matrix W")
  warn_once <- .int64_warn_once()  # a column can be in both X and W
  X <- warn_once(.gather_matrix(src, covariate_cols, dt$.choicer_row,
                                "The design matrix X"))
  X_res <- check_collinearity(X)
  X <- X_res$mat
  if (!is.null(X_res$dropped)) dropped_vars <- X_res$dropped # accumulate dropped vars if we had multiple checks


  W <- warn_once(.gather_matrix(src, random_var_cols, dt$.choicer_row,
                                "The random-coefficient design matrix W"))
  W_res <- check_collinearity(W)
  W <- W_res$mat
  if (!is.null(W_res$dropped)) {
     if(exists("dropped_vars")) dropped_vars <- c(dropped_vars, W_res$dropped)
     else dropped_vars <- W_res$dropped
  }

  dt[, .choicer_row := NULL]
  .check_index_covariates(c(covariate_cols, random_var_cols),
                          c(id_col, alt_col, choice_col, weights_col,
                            cluster_col, person_col))

  ## alternative ids used for delta coefficients
  alt_idx <- as.integer(dt$.choicer_alt_int)                    # length == sum(M)

  ## M[i] - # alternatives per choice situation
  M <- .n_rows_by(dt, id_col)                                   # length N

  ## N: number of individuals / choice situations
  ids <- dt[[id_col]][!duplicated(dt[[id_col]])]  # vector of ids in *current* order
  N   <- length(ids)

  ## Panel layout: Ti[u] consecutive situations for decision maker u, in
  ## sorted person order; NULL in the cross-section.
  Ti <- person_ids <- NULL
  if (!is.null(person_col)) {
    person <- .collapse_situation_col(dt, person_col, id_col, ids)
    person_ids <- unique(person)
    person_pos <- match(person, person_ids)
    Ti <- tabulate(person_pos, length(person_ids))
    stopifnot(!is.unsorted(person_pos))  # contiguous by construction (sort)
    # A positional `weights` vector is read in the prepared order, which is
    # now decision maker first; unless that is also ascending-id order, a
    # vector built in id order would be silently misaligned.
    if (!is.null(weights) && !.ids_sorted(ids)) {
      stop("With `person_col`, choice situations are reordered by decision ",
           "maker, so a positional `weights` vector is ambiguous. Supply the ",
           "weights as a column via `weights_col` (aligned by id).",
           call. = FALSE)
    }
  }

  ## Collapse a row-level weight column to one weight per choice situation.
  ## Done AFTER ordering/filtering so alignment is by id, never by position.
  if (!is.null(weights_col)) {
    if (!is.numeric(dt[[weights_col]])) {
      stop("`", weights_col, "` must be numeric.")
    }
    nuniq <- .n_distinct_by(dt, weights_col, id_col)
    if (any(nuniq != 1L)) {
      stop("`", weights_col, "` must be constant within each '", id_col,
           "' (one weight per choice situation).")
    }
    weights <- .first_by(dt, weights_col, id_col, ids)
    if (any(!is.finite(weights))) {
      stop("`", weights_col, "` produced non-finite weights.")
    }
  }

  ## Collapse a row-level cluster column to one label per choice situation
  ## (same alignment discipline as weights_col).
  cluster <- if (!is.null(cluster_col)) {
    .collapse_situation_col(dt, cluster_col, id_col, ids)
  }

  ## choice_idx[i] - 1-based index *within* the choice set data
  ## 0 == outside option (only if chosen = 0 for all inside options & include_outside_option == TRUE)
  chosen <- dt[[choice_col]] == 1
  if (include_outside_option) {
    # start with all-zero (everyone assumed to pick the outside good), then
    # match the chosen rows' ids (at most one per id) to the master index
    choice_idx <- integer(N)
    choice_idx[match(dt[[id_col]][chosen], ids)] <-
      dt$.choicer_idx_in_group[chosen]
  } else {
    # exactly one explicit choice per id
    choice_idx <- dt$.choicer_idx_in_group[chosen]
  }
  rm(chosen)

  # Weights default = 1
  if (is.null(weights)) weights <- rep(1, N)
  weights <- .int64_to_double(weights, "`weights`")

  ## Weights must be finite and strictly positive. Zero/negative weights would
  ## silently invalidate weighted and WESML sandwich inference (w in the bread,
  ## w^2 in the meat). Validated here so every resolution path (weights=,
  ## weights_col=, and the uniform default) is covered.
  if (any(!is.finite(weights))) {
    stop("Weights must be finite, but non-finite values (NA/NaN/Inf) were found.",
         call. = FALSE)
  }
  if (any(weights <= 0)) {
    stop("Weights must be strictly positive, but values <= 0 were found.",
         call. = FALSE)
  }

  ## Panel: the likelihood has one term per decision maker, so weights and
  ## cluster labels must be constant within each (identity otherwise).
  units <- list(M = M, Ti = Ti)
  .to_units(weights, units, "Weights")
  .to_units(cluster, units, "`cluster_col`")

  ## Alternative summary -------------------------------------------------------
  ## One inside-alternative aggregation; the outside branch only prepends its
  ## synthetic alt_int = 0 row.
  alt_mapping <- .alt_counts(dt, alt_col, choice_col)
  if (include_outside_option) {
    outside_alt_mapping <- data.table::data.table(alt_int=0L, N_OBS = N, N_CHOICES = sum(choice_idx == 0L))
    outside_alt_mapping[[alt_col]] <- outside_opt_label
    alt_mapping <- list(outside_alt_mapping, alt_mapping) |>
      data.table::rbindlist(use.names = TRUE, fill = TRUE)
    data.table::setcolorder(alt_mapping, c("alt_int", alt_col, "N_OBS", "N_CHOICES"))
  }

  alt_mapping[, `:=`(
    TAKE_RATE = N_CHOICES / N_OBS,
    MKT_SHARE = N_CHOICES / sum(N_CHOICES)
  )]

  ## Final validity checks -----------------------------------------------------
  stopifnot(
    length(alt_idx)    == nrow(X),
    length(choice_idx) == N,
    length(M)          == N,
    length(weights)    == N,
    is.null(Ti) || sum(Ti) == N,
    all(is.finite(X)),
    all(is.finite(W))
  )

  ## Return output -------------------------------------------------------------
  out <- structure(
    list(
      X           = X,
      W           = W,
      alt_idx     = alt_idx,
      choice_idx  = as.integer(choice_idx),
      M           = M,
      N           = N,
      weights     = weights,
      cluster     = cluster,
      situation_ids = ids,
      Ti          = Ti,
      person_ids  = person_ids,
      include_outside_option = include_outside_option,
      rc_correlation = rc_correlation,
      alt_mapping = alt_mapping[],
      dropped_cols = if(exists("dropped_vars")) dropped_vars else NULL,
      data_spec = list(
        id_col = id_col,
        alt_col = alt_col,
        choice_col = choice_col,
        covariate_cols = covariate_cols,
        random_var_cols = random_var_cols,
        outside_opt_label = outside_opt_label,
        person_col = person_col
      )
    ),
    class = "choicer_data_mxl"
  )
  if (!is.null(cs_provenance)) {
    attr(out, "choice_sampling") <- cs_provenance
  }
  out
}

#' Halton draws for mixed logit
#'
#' Create halton normal draws in appropriate format for mixed logit estimation
#'
#' Draw unit \eqn{i} receives the \eqn{S} consecutive points
#' \eqn{(i - 1)S + 1, \ldots, iS} of the \eqn{K_w}-dimensional Halton sequence
#' of \code{randtoolbox::halton()}, mapped to standard normals. The cube is
#' filled a block of units at a time, so it never holds a second copy of the
#' sequence. Store mode supports at most \eqn{2^{31} - 1} points
#' (\eqn{S \times N}), the largest starting index \code{halton(start = )}
#' accepts, and \eqn{2^{32} - 1} values (\eqn{K_w \times S \times N}), the
#' largest cube choicer supports;
#' \code{draws = "generate"} in \code{\link{run_mxlogit}} has neither limit and
#' never materializes the cube.
#'
#' @param S Number of draws per draw unit
#' @param N number of draw units: choice situations, or decision makers for a
#'   panel fit (\code{person_col} in \code{\link{run_mxlogit}})
#' @param K_w dimension of random coefficients (number of columns in W matrix)
#' @returns K_w x S x N array with halton standard normal draws
#' @examples
#' draws <- get_halton_normals(S = 50, N = 10, K_w = 2)
#' dim(draws)  # 2 x 50 x 10
#' @importFrom randtoolbox halton
#' @export
get_halton_normals <- function(S, N, K_w) {
  bad <- function(x) {
    !is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != round(x)
  }
  if (bad(S)) stop("`S` must be a single positive whole number.")
  if (bad(N)) stop("`N` must be a single positive whole number.")
  if (bad(K_w)) stop("`K_w` must be a single positive whole number.")
  n_points <- as.numeric(S) * N
  if (n_points > .Machine$integer.max) {
    stop("A store-mode draw cube needs S * N = ",
         format(n_points, big.mark = ",", scientific = FALSE),
         " points of the Halton sequence, more than 2^31 - 1, the largest ",
         "starting index randtoolbox::halton(start = ) accepts. Refit with ",
         "run_mxlogit(draws = \"generate\").")
  }
  # The kernels addressed the cube with Armadillo's 32-bit indices, which
  # misread a larger one; choicer now defines ARMA_64BIT_WORD, and this stop
  # stays until cubes that large have been run end to end.
  if (n_points * K_w > 2^32 - 1) {
    stop("A store-mode draw cube needs K_w * S * N = ",
         format(n_points * K_w, big.mark = ",", scientific = FALSE),
         " values, more than 2^32 - 1, the largest cube choicer supports. ",
         "Refit with run_mxlogit(draws = \"generate\").")
  }
  .halton_cube(S, N, K_w)
}

#' Fill the store-mode Halton cube a block of draw units at a time
#'
#' Unit i's draws are points (i - 1) * S + 1, ..., i * S of the sequence.
#' `randtoolbox::halton(start = )` (randtoolbox >= 1.31.0) reproduces any slice
#' of the full sequence bit for bit, so the blocks give exactly the draws of a
#' single call over all S * N points. They avoid that call's copies of the
#' whole sequence and its 32-bit offsets into the output, which overflow past
#' 2^31 - 1 values.
#'
#' @param S,N,K_w As in [get_halton_normals()], already validated.
#' @param block Target number of values generated per `halton()` call; a
#'   block holds at least one unit.
#' @returns K_w x S x N array.
#' @noRd
.halton_cube <- function(S, N, K_w, block = 2^22) {
  eta <- array(0, dim = c(K_w, S, N))
  units_per_block <- max(1, floor(block / (as.numeric(S) * K_w)))
  i0 <- 1
  while (i0 <= N) {
    i1 <- min(N, i0 + units_per_block - 1)
    # (S * (i1 - i0 + 1)) x K_w, or a vector when K_w = 1; its transpose is
    # eta[, , i0:i1] in memory order either way.
    h <- randtoolbox::halton(n = S * (i1 - i0 + 1), dim = K_w, normal = TRUE,
                             start = (i0 - 1) * S + 1)
    eta[, , i0:i1] <- t(h)
    i0 <- i1 + 1
  }
  eta
}


#' Keep the likelihood kernel's overflow sentinel above the optimizer's path
#'
#' `mxl_loglik_gradient_parallel()` returns `objective = 1e10` with a zero
#' gradient and `overflow = TRUE` where the utilities overflow, so that a
#' line search backtracks. A finite objective can also equal 1e10, so only
#' the explicit flag identifies the sentinel.
#' Its objective is otherwise finite however poor the fit, and unbounded, so
#' after a start above 1e10 (a badly scaled warm start at population scale)
#' the fixed sentinel would look like an improvement. Wraps `eval_f` to report
#' the sentinel as ten times the largest objective seen so far, and never
#' below 1e10; objectives are negated log-likelihoods, hence nonnegative.
#'
#' @param eval_f Function of theta returning `list(objective, gradient, overflow)`.
#' @returns The wrapped function.
#' @noRd
.lift_sentinel <- function(eval_f) {
  f_max <- 0
  function(theta) {
    res <- eval_f(theta)
    if (isTRUE(res$overflow)) {
      res$objective <- max(1e10, min(10 * f_max, .Machine$double.xmax))
    } else if (is.finite(res$objective)) {
      f_max <<- max(f_max, res$objective)
    }
    res$overflow <- NULL
    res
  }
}


#' Resolve draw parameters for post-estimation regeneration sites
#'
#' When the fitted object used generate mode, returns an empty placeholder cube
#' plus the three gen_* integers. When in store mode, materialises the Halton
#' cube from the stored metadata.
#'
#' @param draws_info List from a fitted choicer_mxl object.
#' @param N Number of draw blocks. Prediction sites keep the default, one block
#'   per choice situation; estimation-type sites (Hessian, scores, conditional
#'   tastes) pass the number of likelihood units, `length(.unit_first(d))`,
#'   which is the number of decision makers for a panel fit.
#' @noRd
.mxl_gen_params <- function(draws_info, N = draws_info$N) {
  mode <- draws_info$mode %||% "store"
  if (mode == "generate") {
    list(
      eta_draws    = array(0, dim = c(draws_info$K_w, 0L, 0L)),
      gen_seed     = as.integer(draws_info$seed),
      # "owen" is the legacy serialized label for the same position-wise
      # permutation. Read it silently so old fitted objects keep working.
      gen_scramble = if (draws_info$scramble %in% c("permuted", "owen")) 1L else 0L,
      gen_S        = as.integer(draws_info$S)
    )
  } else {
    list(
      eta_draws    = get_halton_normals(draws_info$S, N, draws_info$K_w),
      gen_seed     = -1L,
      gen_scramble = 1L,
      gen_S        = 0L
    )
  }
}
