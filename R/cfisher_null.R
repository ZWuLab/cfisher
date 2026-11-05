# ==============================================================================
# File: cfisher_null.R
# Purpose: The null distribution of the count T of exposed cases. For subjects with
#   exposure x, outcome y (1 = case) and fitted null probabilities mu_i of being a case
#   (a logistic model of y on the confounders without x), T is the number of cases among
#   the exposed subjects. Its null distribution is either the conditional Bernoulli
#   distribution of T given the exposed set and the case total r (Fisher's exact test
#   controlling for confounders), or the Poisson-binomial distribution of T given the
#   exposed set only (the unconditional test).
#
# EXPORTED FUNCTIONS:
#   - cfisher.null()         The null distribution of T as a data frame (support, prob)
#
# INTERNAL HELPERS:
#   - .prepare()             Validate the inputs and compute the quantities every test needs
#   - .null_distribution()   Conditional or unconditional null distribution of T
#   - .pb_pmf()              Poisson-binomial pmf of T, in log scale
#   - .cb_pmf()              Conditional Bernoulli pmf of T given the case total, in log scale
#   - .pb_log()              Log Poisson-binomial pmf: probability scale when an underflow
#                            bound allows it, log scale otherwise (calls the C++ code)
#   - .pb_scaled()           Scaled Poisson-binomial recursion in probability scale
#   - .lse()                 Log of a sum of exponentials, overflow-safe
#   - .center_total()        Shift the log-odds so that the fitted probabilities sum to r
#   - .fitted_null()         Fitted probabilities and log-odds of the logistic null model
#   - .design_matrix()       Design matrix of the confounders, with an intercept
#   - .check_binary()        Validation of 0/1 vectors
#
# Dependencies: base R, stats and the compiled recursion in src/pb_pmf.cpp (Rcpp).

#################### EXPORTED MAIN FUNCTIONS ####################

#' Null distribution of the number of exposed cases
#'
#' Returns, as a data frame, the null distribution of \eqn{T}, the number of cases among
#' the exposed subjects, that \code{\link{cfisher.test}} uses. Other packages can take it
#' as a discrete null distribution (a support and its probabilities).
#'
#' @inheritParams cfisher.test
#'
#' @return A data frame with columns \code{support} (the possible values of \eqn{T}, in
#'   increasing order) and \code{prob} (their probabilities, summing to one). A
#'   probability below the range of double precision (about \eqn{10^{-308}}) is shown as
#'   0 or with reduced precision; \code{\link{cfisher.test}} works with the logarithms,
#'   which keep such probabilities.
#'
#' @details
#' Let \eqn{m} be the number of exposed subjects, \eqn{r} the number of cases and
#' \eqn{\mu_i} the fitted null probability that subject \eqn{i} is a case, with odds
#' \eqn{w_i = \mu_i / (1 - \mu_i)}.
#'
#' With \code{conditional = TRUE} the distribution is the conditional Bernoulli
#' distribution of \eqn{T} given the exposed set and \eqn{r}:
#' \deqn{P(T = t) = e_t(w_E) e_{r-t}(w_{E^c}) / e_r(w),}
#' where \eqn{E} is the exposed set and \eqn{e_k} is the elementary symmetric function of
#' order \eqn{k} of the odds in a set. It is computed as the ratio of Poisson-binomial
#' probabilities
#' \deqn{P(T = t) = f_E(t) f_{E^c}(r - t) / \sum_s f_E(s) f_{E^c}(r - s),}
#' where \eqn{f_A} is the probability mass function of the number of cases in the set
#' \eqn{A}, by the forward recursion over subjects, in \eqn{O(m^2 + (n - m) r)} time, so
#' that no elementary symmetric function is formed. The distribution does not change when
#' every odds is multiplied by a constant, so before the recursion the log-odds are shifted by the constant that makes
#' \eqn{\sum_i \mu_i = r}. This keeps the recursion at the center of the distributions;
#' for a fitted model with an intercept the shift is zero. If all odds are equal the
#' distribution is hypergeometric, and if the odds take one value on the exposed and
#' another on the unexposed it is Fisher's noncentral hypergeometric distribution.
#'
#' The recursion runs in probability scale, which is fast, as long as a bound on its error
#' from underflow is below \eqn{10^{-12}} of every probability it has to deliver (roughly,
#' of those above \eqn{10^{-290}}). Otherwise it runs in log scale, at about 20 times the
#' cost, so that no probability is lost to underflow, even far below the range of double
#' precision. In both scales the probabilities carry the rounding error of
#' double-precision arithmetic, a relative error that grows with the number of subjects.
#'
#' With \code{conditional = FALSE} it is the Poisson-binomial distribution of \eqn{T}
#' given the exposed set, the distribution of a sum of \eqn{m} independent Bernoulli
#' variables with probabilities \eqn{\mu_i}, \eqn{i \in E}, computed in \eqn{O(m^2)} time.
#' The number of cases is left random.
#'
#' Assumptions: given the confounders, the outcomes are independent Bernoulli variables
#' with the null probabilities \eqn{\mu_i} (the outcome does not depend on the exposure
#' given the confounders). The fitted probabilities are plugged in as if known.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' z <- rnorm(n)
#' x <- rbinom(n, 1, plogis(-2 + z))
#' y <- rbinom(n, 1, plogis(-1 + z))
#' cfisher.null(x, y, z)
#' cfisher.null(x, y, z, conditional = FALSE)
#'
#' @references
#' Wu, Z. and Contador, G. (2026). Fisher's exact test controlling for confounders. Draft
#' manuscript.
#'
#' Chen, S. X. and Liu, J. S. (1997). Statistical applications of the Poisson-binomial and
#' conditional Bernoulli distributions. \emph{Statistica Sinica}, 7(4), 875-892.
#'
#' @seealso \code{\link{cfisher.test}}
#' @export
cfisher.null <- function(x, y, z = NULL, conditional = TRUE, fit = NULL) {
  d <- .prepare(x, y, z, fit)
  nd <- .null_distribution(d, conditional)
  data.frame(support = nd$support, prob = nd$prob)
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' Validate the inputs and compute what every test needs
#'
#' Returns a list with the validated 0/1 vectors x and y, the null probabilities mu, the
#' exposed set (indices), the exposed count m, the case total r, the sample size n and
#' the observed count t of exposed cases; eta is the vector of null log-odds (from the
#' linear predictor of a logit fit, see .fitted_null()).
#' @keywords internal
#' @noRd
.prepare <- function(x, y, z, fit) {
  x <- .check_binary(x, "x")
  y <- .check_binary(y, "y")
  n <- length(y)
  if (length(x) != n) stop("'x' and 'y' must have the same length.")
  if (all(x == 0L) || all(x == 1L))
    stop("The exposure 'x' is constant (all ", x[1], "); the count of exposed cases has ",
         "no null distribution to test. Drop it before calling cfisher.test().")
  if (all(y == 0L) || all(y == 1L))
    stop("The outcome 'y' is constant (all ", y[1], "); there is nothing to test.")
  nf <- .fitted_null(y, z, fit)
  exposed <- which(x == 1L)
  list(x = x, y = y, mu = nf$mu, eta = nf$eta, exposed = exposed,
       m = length(exposed), r = sum(y), n = n, t = sum(y[exposed]))
}

#' The null distribution of T as list(support, prob), from the output of .prepare()
#' @keywords internal
#' @noRd
.null_distribution <- function(d, conditional) {
  if (!is.logical(conditional) || length(conditional) != 1L || is.na(conditional))
    stop("'conditional' must be TRUE or FALSE.")
  if (conditional) return(.cb_pmf(d$eta, d$exposed, d$r))
  .pb_pmf(d$eta[d$exposed])
}

#' The Poisson-binomial distribution of a sum of Bernoulli(plogis(eta_i)) variables as
#' list(support, prob, logprob), from the log-odds eta; no probability is lost to
#' underflow (exact = FALSE: see .pb_log())
#' @keywords internal
#' @noRd
.pb_pmf <- function(eta, exact = TRUE) {
  lp <- .pb_log(eta, exact = exact)
  list(support = 0:length(eta), prob = exp(lp), logprob = lp)
}

#' Conditional Bernoulli pmf of a subset count given the total
#'
#' The pmf of T = sum over the exposed set E of Y_i, given that the sum of all Y_i is r,
#' for independent Bernoulli(mu_i) variables with log-odds eta_i = qlogis(mu_i):
#' P(T = t | N = r) = e_t(w_E) e_{r-t}(w_Ec) / e_r(w), computed as
#' f_E(t) f_Ec(r - t) / sum_s f_E(s) f_Ec(r - s) with f the Poisson-binomial pmfs,
#' combined in log scale. Returns list(support, prob, logprob), where logprob keeps the
#' far tails that prob underflows (exact = FALSE: see .pb_log(); the probabilities are
#' then accurate in absolute terms). Takes the log-odds, so that no probability near 1
#' loses its complement. Assumes every eta_i finite and E a non-empty proper subset. For
#' r = 0 and r = n the pmf is a point mass; otherwise the log-odds are first centered
#' (see .center_total()).
#' @keywords internal
#' @noRd
.cb_pmf <- function(eta, exposed, r, exact = TRUE) {
  n <- length(eta)
  exposed <- as.integer(exposed)
  m <- length(exposed)
  if (m == 0L || m == n) stop("'exposed' must be a non-empty proper subset of the subjects.")
  if (any(!is.finite(eta)))
    stop("The null log-odds must be finite (probabilities strictly inside (0, 1)).")
  r <- as.integer(r)
  if (r < 0L || r > n) stop("The case total must be between 0 and the number of subjects.")
  if (r == 0L) return(list(support = 0L, prob = 1, logprob = 0))   # no cases: T = 0
  if (r == n)  return(list(support = m, prob = 1, logprob = 0))    # all cases: T = m
  # Center the total: the conditional pmf is homogeneous in the odds, so shifting every
  # log-odds by a constant leaves it unchanged. The constant that makes sum(mu) = r puts
  # both recursions at the center of their distributions, whatever r is. The shift stays
  # on the log-odds scale (D2 of the review round 2026-11-04).
  eta <- .center_total(eta, r)
  t_min <- max(0L, r - (n - m))
  t_max <- min(m, r)
  support <- t_min:t_max
  # log P(sum over E = t) at t = 0..m and log P(sum over E^c = s) at s = 0..min(r, n - m),
  # each accurate at the indices the support needs
  lf_e <- .pb_log(eta[exposed], need = support, exact = exact)
  lf_u <- .pb_log(eta[-exposed], kmax = r, need = r - support, exact = exact)
  lp <- lf_e[support + 1L] + lf_u[r - support + 1L]  # log of f_E(t) f_Ec(r - t)
  ltot <- .lse(lp)                                   # log of the normalizing sum
  if (!is.finite(ltot))
    stop("The normalizing sum is not finite and positive; this should not happen. ",
         "r = ", r, ", sum(mu) = ", signif(sum(stats::plogis(eta)), 4), ".")
  lp <- lp - ltot
  list(support = support, prob = exp(lp), logprob = lp)
}

# Bound on the relative error from underflow that .pb_log() accepts for the probabilities it
# computes in probability scale; above it, they are recomputed in log scale.
.PB_REL_TOL <- 1e-12

#' Log Poisson-binomial pmf at 0..kmax, from the log-odds eta, with no needed index lost
#' to underflow
#'
#' Runs the probability-scale recursion and keeps it when a bound on its error from
#' underflow is below .PB_REL_TOL times every needed probability; otherwise reruns the
#' recursion in log scale. The bound covers underflow only. Each update of an entry is
#' q_j f(k) + mu_j f(k - 1), and its inputs, its products and their sum can lose at most a
#' few multiples of 2^-1075 to underflow in the working scale, where every entry is at most
#' 1 and which is never smaller than probability scale (rescaling only enlarges it). Every
#' later update is a convex combination, so such an error reaches any final entry with a
#' weight of at most 1. With at most n * (kmax + 1) updates, the absolute error from
#' underflow of a final probability is below n (kmax + 1) 2^-1072. Besides it, both
#' recursions carry the ordinary rounding error of double precision, a relative error that
#' grows with n (and, in log scale, with the size of the log probability); the bound does
#' not cover it. The inputs are plogis(eta) and plogis(-eta), so that a probability near 1
#' keeps its complement, and the log-scale recursion takes their logs from plogis(...,
#' log.p = TRUE), which never underflows. With exact = FALSE the probability-scale result
#' is returned as it is, which suffices for tail sums and means; its far tail can be -Inf.
#' @keywords internal
#' @noRd
.pb_log <- function(eta, kmax = length(eta), need = 0:kmax, exact = TRUE) {
  eta <- as.numeric(eta)
  kmax <- as.integer(min(kmax, length(eta)))
  s <- .pb_scaled(eta, kmax)
  lf <- log(s$f) - s$logscale                        # log of the probability-scale result
  log_floor <- log(length(eta)) + log(kmax + 1) - 1072 * log(2) - log(.PB_REL_TOL)
  if (!exact || all(lf[need + 1L] >= log_floor)) return(lf)
  # some needed probability is too small for the bound: the log-scale recursion
  .pb_log_cpp(stats::plogis(eta, log.p = TRUE), stats::plogis(-eta, log.p = TRUE), kmax)
}

#' The scaled Poisson-binomial recursion from the log-odds eta: list(f, logscale) with
#' pmf = f * exp(-logscale)
#' @keywords internal
#' @noRd
.pb_scaled <- function(eta, kmax = length(eta)) {
  eta <- as.numeric(eta)
  kmax <- as.integer(min(kmax, length(eta)))
  .pb_scaled_cpp(stats::plogis(eta), stats::plogis(-eta), kmax)
}

#' Log of the sum of exp(v), shifted by the maximum so that it neither overflows nor
#' underflows; -Inf when every entry is -Inf
#' @keywords internal
#' @noRd
.lse <- function(v) {
  mx <- max(v)
  if (!is.finite(mx)) return(mx)
  mx + log(sum(exp(v - mx)))
}

#' Shift the log-odds eta by a common constant c so that sum(plogis(eta + c)) = r
#' (0 < r < n); returns the shifted log-odds eta + c
#' @keywords internal
#' @noRd
.center_total <- function(eta, r) {
  g <- function(c) sum(stats::plogis(eta + c)) - r
  if (abs(g(0)) < 1e-10 * length(eta)) return(eta)
  # Newton's method on g(c) = sum(plogis(eta + c)) - r, which is increasing and smooth
  # (g' = sum of mu (1 - mu) > 0). An undamped step can saturate every probability when
  # the target is far from sum(mu) with few subjects, so a zero or non-finite derivative
  # or value hands over to the bracketing fallback, which always converges.
  c0 <- 0
  for (it in 1:50) {
    m_c <- stats::plogis(eta + c0); gc <- sum(m_c) - r
    if (!is.finite(gc)) break
    if (abs(gc) < 1e-9) return(eta + c0)
    dg <- sum(m_c * stats::plogis(-(eta + c0)))      # sum of mu (1 - mu), without 1 - mu
    if (!is.finite(dg) || dg <= 0) break
    c0 <- c0 - gc / dg
  }
  lo <- -1; hi <- 1
  while (g(lo) > 0) lo <- 2 * lo
  while (g(hi) < 0) hi <- 2 * hi
  eta + stats::uniroot(g, c(lo, hi), tol = 1e-12)$root
}

#' Fitted probabilities and log-odds of the logistic null model of y on the confounders z
#'
#' fit takes precedence: a binomial glm (its fitted values) or a vector of probabilities.
#' Otherwise a logistic model with intercept is fitted by glm.fit; z NULL gives the
#' intercept-only model, whose fitted probability is mean(y). Returns list(mu, eta), with
#' eta the log-odds taken from the linear predictor of a logit fit, so that a fitted value
#' near 1 keeps its complement, and qlogis(mu) otherwise.
#' @keywords internal
#' @noRd
.fitted_null <- function(y, z, fit) {
  n <- length(y)
  if (!is.null(fit)) {
    eta <- NULL                                         # set from a logit glm below
    if (inherits(fit, "glm")) {
      if (!identical(fit$family$family, "binomial"))
        stop("'fit' must be a glm of the binomial family.")
      mu <- stats::fitted(fit)
      eta <- if (identical(fit$family$link, "logit")) as.numeric(fit$linear.predictors)
    } else {
      mu <- fit
      if (!is.numeric(mu)) stop("'fit' must be a binomial glm or a numeric vector of probabilities.")
    }
    mu <- as.numeric(mu)
    if (length(mu) != n) stop("'fit' must give one fitted probability per subject.")
    if (any(!is.finite(mu)) || any(mu <= 0) || any(mu >= 1))
      stop("The fitted probabilities must lie strictly inside (0, 1).")
    if (is.null(eta) || length(eta) != n) eta <- stats::qlogis(mu)
    return(list(mu = .warn_separation(mu), eta = eta))
  }
  if (is.null(z)) {
    mu <- rep(mean(y), n)
    return(list(mu = mu, eta = stats::qlogis(mu)))
  }
  zm <- .design_matrix(z, n)
  g <- stats::glm.fit(zm, y, family = stats::binomial())
  mu <- as.numeric(g$fitted.values)
  if (any(mu <= 0) || any(mu >= 1))
    stop("The null model fits some probabilities at 0 or 1 (separation); ",
         "the test is not defined for them.")
  list(mu = .warn_separation(mu), eta = as.numeric(g$linear.predictors))
}

#' Warn when a fitted null probability is within 1e-8 of 0 or 1; returns mu
#'
#' Such values signal (quasi-)separation of the outcome by the confounders, or an
#' extreme supplied probability. The null distribution is then computed, but the odds
#' ratio estimate and interval can be unreliable (a flat likelihood).
#' @keywords internal
#' @noRd
.warn_separation <- function(mu) {
  if (any(mu < 1e-8 | mu > 1 - 1e-8))
    warning("Some fitted null probabilities are within 1e-8 of 0 or 1, which suggests ",
            "(quasi-)separation of the outcome by the confounders; those subjects are ",
            "nearly certain cases or controls and carry almost no information on the ",
            "exposure.", call. = FALSE)
  mu
}

#' Design matrix of the confounders z with an intercept (an intercept alone for NULL)
#'
#' z is a vector, a matrix or a data frame (factors and characters become dummy columns).
#' @keywords internal
#' @noRd
.design_matrix <- function(z, n) {
  if (is.null(z)) return(matrix(1, n, 1L))
  if (anyNA(z)) stop("'z' must not contain missing values.")
  if (is.factor(z) || is.character(z)) z <- data.frame(z = z)
  zm <- if (is.data.frame(z)) stats::model.matrix(~ ., data = z) else cbind(1, as.matrix(z))
  if (nrow(zm) != n) stop("'z' must have one row per subject.")
  zm
}

#' Validate a 0/1 vector and return it as integer
#' @keywords internal
#' @noRd
.check_binary <- function(v, name) {
  if (is.logical(v)) v <- as.integer(v)
  if (!is.numeric(v) || anyNA(v) || !all(v %in% c(0, 1)))
    stop("'", name, "' must be a 0/1 (or logical) vector without missing values.")
  as.integer(v)
}
