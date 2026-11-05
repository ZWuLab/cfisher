# ==============================================================================
# File: cfisher_test.R
# Purpose: Fisher's exact test of a binary exposure against a binary outcome while
#   controlling for confounders. The statistic is the count T of exposed cases; its null
#   distribution (R/cfisher_null.R) is the conditional Bernoulli distribution given the
#   exposed set and the case total (the default), or the Poisson-binomial distribution
#   given the exposed set. This file holds the test, its p-values, and the estimate and
#   confidence interval of the odds ratio, which come from inverting the test over
#   the hypothesized odds ratio, with the confounder slopes refitted under each hypothesized
#   value (or, when only null probabilities are supplied, from the exponential tilt of the
#   null distribution).
#
# EXPORTED FUNCTIONS:
#   - cfisher.test()         The test; returns an object of class "htest"
#
# INTERNAL HELPERS:
#   - .discrete_pvalue()     Left, right and two-sided (mid-)p-values on a pmf in log scale
#   - .odds_ratio_inference() Estimate and interval of the odds ratio by inverting the test
#   - .tilt()                Exponential tilt of a pmf by exp(beta * t), from log probabilities
#   - .solve_increasing()    Root of an increasing function by an expanding bracket
#   - .logistic_refit()      Logistic MLE with an offset by damped Newton (the refit)
#   - .full_mle()            MLE of the exposure's coefficient (unconditional estimate)
#   - .second_crossing()     Whether a solved function crosses its target a second time
#   - .nuisance_design()     Design matrix of the null model for the refit (z or a glm)
#   - .unidentified()        Whether the exposure lies in the span of the null model's design
#   - .in_span()             Whether a vector lies in the column space of a matrix
#
# Dependencies: base R and stats only (plus the compiled recursion via cfisher_null.R).

#################### EXPORTED MAIN FUNCTIONS ####################

#' Fisher's exact test controlling for confounders
#'
#' Tests whether a binary exposure \code{x} is associated with a binary outcome \code{y}
#' while controlling for confounders \code{z}, in the manner of Fisher's exact test. The
#' statistic is \eqn{T}, the number of cases among the exposed subjects. Its null
#' distribution is exact, given a logistic model of the outcome on the confounders that
#' does not contain the exposure. Without confounders the test reproduces
#' \code{\link[stats]{fisher.test}} on the \eqn{2 \times 2} table of \code{x} and \code{y}.
#'
#' @param x Integer or logical vector of 0/1 exposure values, one per subject.
#' @param y Integer or logical vector of 0/1 outcomes (1 = case), the same length as
#'   \code{x}.
#' @param z Confounders: \code{NULL} (the default; the null model has an intercept only,
#'   so the test is Fisher's exact test), or a numeric vector, matrix or data frame with
#'   one row per subject. An intercept is always added. Ignored when \code{fit} is given.
#' @param alternative Character. \code{"two.sided"} (the default), \code{"greater"} (the
#'   exposed are enriched among the cases, p-value \eqn{P(T \ge t)}) or \code{"less"}
#'   (\eqn{P(T \le t)}).
#' @param conf.int Logical. Whether to compute the estimate and the confidence interval of
#'   the odds ratio. If \code{FALSE} both are omitted, which is faster when many tests are
#'   run. (Unlike \code{\link[stats]{fisher.test}}, which still returns its estimate, the
#'   omission of the estimate here is intentional, for speed in scans.)
#' @param conf.level Confidence level of the interval, strictly between 0 and 1, such that
#'   each tail probability of the interval lies between \eqn{10^{-12}} and
#'   \eqn{1 - 10^{-12}}.
#' @param midp Logical. If \code{TRUE} the observed count carries half its probability in
#'   the p-value (the mid-p value). Only the p-value changes; the estimate and the
#'   confidence interval are those of the exact test.
#' @param conditional Logical. \code{TRUE} (the default): the null distribution is the
#'   conditional Bernoulli distribution of \eqn{T} given the exposed subjects and the
#'   number of cases, which is Fisher's exact test controlling for confounders. \code{FALSE}:
#'   the null distribution is the Poisson-binomial distribution of \eqn{T} given the
#'   exposed subjects only, the unconditional test (the null of the BinomiRare test and of
#'   the conditional randomization test).
#' @param fit Optional null probabilities, which take precedence over \code{z}: a fitted
#'   \code{\link[stats]{glm}} object of the binomial family whose fitted values are the
#'   probabilities that each subject is a case under the null hypothesis, or a numeric
#'   vector of such probabilities, strictly between 0 and 1, one per subject. Passing
#'   \code{fit} lets a scan over many exposures fit the null model once. The model must not
#'   contain the exposure. For the odds ratio, a \code{glm} must have the logit link and
#'   no weights or offset, and its model matrix is used to refit the slopes; a vector of
#'   probabilities allows no refit (see Details).
#'
#' @return An object of class \code{"htest"} with components
#' \item{statistic}{the observed count \eqn{T = t} of exposed cases.}
#' \item{parameter}{\code{exposed} (the number \eqn{m} of exposed subjects), \code{cases}
#'   (the number \eqn{r} of cases) and \code{n} (the number of subjects).}
#' \item{p.value}{the p-value.}
#' \item{conf.int}{a confidence interval for the odds ratio, with attribute
#'   \code{conf.level}. For the two-sided alternative each tail carries half of one minus
#'   \code{conf.level}, as in \code{\link[stats]{fisher.test}}; for the one-sided
#'   alternatives the interval is one-sided. A limit is \code{0} or \code{Inf} when the
#'   observed count is at an end of the support or when the limit lies beyond the search
#'   range (see Details), and \code{NA} (with a warning) when its root search fails or the
#'   effect is not identified. Omitted if \code{conf.int = FALSE}.}
#' \item{estimate}{the estimate of the odds ratio, \eqn{\exp(\hat\beta)}, where
#'   \eqn{\hat\beta} solves \eqn{E_\beta(T) = t}; \code{0} or \code{Inf} when the observed
#'   count is at an end of the support or when \eqn{E_\beta(T)} reaches \eqn{t} only in
#'   the limit (see Details), and \code{NA} (with a warning) when its root search fails,
#'   when it lies beyond the search range under the conditional distribution, when the
#'   observed count is certain at every odds ratio or when the effect is not identified.
#'   With \code{conditional = FALSE} it is named \code{"odds ratio"}, because
#'   it is then the maximum likelihood estimate under the unconditional distribution.
#'   Omitted if \code{conf.int = FALSE}.}
#' \item{null.value}{the odds ratio under the null hypothesis, 1.}
#' \item{alternative}{the alternative.}
#' \item{method}{a string naming the test, whether it is conditional and mid-p, and how
#'   the odds ratio was obtained (when it is computed), which is either\cr
#'   \code{"slopes refitted under each hypothesized odds ratio"} or\cr
#'   \code{"tilt at the supplied odds"}.}
#' \item{data.name}{the names of the data.}
#' \item{null.distribution}{a data frame with columns \code{support} and \code{prob}, the
#'   null distribution of \eqn{T} that the test used.}
#' \item{expected}{the mean of \eqn{T} under that null distribution.}
#'
#' @details
#' Let \eqn{E} be the exposed subjects, \eqn{m = |E|}, \eqn{r} the number of cases and
#' \eqn{\mu_i} the null probability that subject \eqn{i} is a case, the fitted value of the
#' logistic regression of \code{y} on \code{z} (with intercept, without \code{x}), with
#' odds \eqn{w_i = \mu_i / (1 - \mu_i)}.
#'
#' \emph{Null distribution.} Given the confounders and the exposure, the null hypothesis
#' makes the outcomes independent Bernoulli(\eqn{\mu_i}) variables. The exposed set is held
#' fixed. With \code{conditional = TRUE} the number of cases \eqn{r} is held fixed too, and
#' \eqn{T} has the conditional Bernoulli distribution
#' \deqn{q(t) = e_t(w_E) e_{r-t}(w_{E^c}) / e_r(w),}
#' with \eqn{e_k} the elementary symmetric function of order \eqn{k} of the odds in a set.
#' It is the hypergeometric distribution when all odds are equal, and Fisher's noncentral
#' hypergeometric distribution when the odds take one value on the exposed and another on
#' the unexposed. It is computed as a ratio of Poisson-binomial probabilities by a forward
#' recursion over subjects in \eqn{O(m^2 + (n - m) r)} time (see \code{\link{cfisher.null}}).
#' With \code{conditional = FALSE}, \eqn{T} has the Poisson-binomial distribution of the
#' exposed subjects' probabilities.
#'
#' \emph{P-values.} The one-sided p-values are \eqn{P(T \ge t)} and \eqn{P(T \le t)}. The
#' two-sided p-value is the sum of the probabilities of the support points that are no more
#' probable than the observed one, the convention of \code{\link[stats]{fisher.test}}; two
#' probabilities are treated as tied when they agree to a relative tolerance of
#' \eqn{10^{-7}}, the tolerance of \code{fisher.test}. The mid-p value subtracts half the
#' probability of the observed count.
#'
#' \emph{Odds ratio.} The exposure enters the logistic model with coefficient \eqn{\beta},
#' so that \eqn{e^\beta} is the odds ratio of the exposure given the confounders. For each
#' hypothesized \eqn{\beta} the slopes of the confounders are refitted by maximum likelihood
#' with \eqn{\beta x_i} as an offset, the null distribution of \eqn{T} is computed at the
#' resulting odds (conditional or Poisson-binomial, as chosen), and the test is inverted. The
#' estimate solves \eqn{E_\beta(T) = t}. The lower limit solves \eqn{P_\beta(T \ge t) =
#' \alpha/2} and the upper limit \eqn{P_\beta(T \le t) = \alpha/2}, where \eqn{\alpha} is
#' one minus \code{conf.level}; a one-sided interval puts all of \eqn{\alpha} in one tail. The
#' equations are solved for \eqn{\beta} in \eqn{[-40, 40]} by root finding with an
#' expanding bracket. A limit beyond that range is reported as \code{0} or \code{Inf}. The
#' estimate is \code{0} or \code{Inf} when \eqn{E_\beta(T) - t} is within \eqn{10^{-9}} of
#' zero at an end of that range, which is taken to mean that \eqn{E_\beta(T)} reaches
#' \eqn{t} only in the limit, as when the logistic model that contains the exposure has no
#' finite maximum likelihood estimate. With \code{conditional = FALSE} and
#' the slopes refitted, any other estimate beyond the range is found as the maximum
#' likelihood estimate of \eqn{\beta} in the logistic model that contains the exposure,
#' because \eqn{E_\beta(T) = t} is then the profile likelihood equation of that model for
#' \eqn{\beta}. Otherwise such an estimate is \code{NA} with a warning. Each refit starts
#' from the nearest refit already computed, and a refit that does not converge makes the
#' quantity being solved \code{NA}, with a warning. The refit runs on the columns of the
#' design that the null fit keeps. The refitted tail probabilities need
#' not be monotone in \eqn{\beta}, so the set of odds ratios that the test does not reject
#' need not be an interval; the function warns when the values it evaluated show more than
#' one crossing.
#' Without confounders, the conditional test (the default) gives the estimate and the
#' interval of \code{\link[stats]{fisher.test}}. Refitting is needed because the null model,
#' fitted without the exposure, absorbs part of the exposure effect when the exposure is
#' common, so tilting its distribution estimates the wrong family under an alternative.
#' The refit is a plug-in: the p-value and the interval treat the fitted slopes as known, and
#' they are exact at the true slopes. When the exposure is a linear combination of the
#' intercept and the confounders, its effect is not identified: the function warns, and
#' the estimate and the interval are \code{NA}. The coverage of the refitted interval was checked in a
#' simulation reported in the accompanying paper (Wu and Contador, 2026, a draft; see
#' References). When \code{fit} is a numeric vector, no
#' refit is possible, and the distribution at \eqn{\beta} is the exponential tilt
#' \eqn{q_\beta(t) \propto q(t) e^{\beta t}} of the distribution at those probabilities.
#' This is exact when they are the baseline probabilities under the alternative (the
#' probabilities at \eqn{\beta = 0} of the model that contains the exposure), and
#' otherwise it is not recommended.
#'
#' \emph{Numerical range.} No probability of the null distribution is lost to underflow, at
#' any point of the support \eqn{\max(0, r - n + m) \le t \le \min(m, r)}, including
#' those far below the range of double precision (about \eqn{10^{-308}}). The recursion
#' runs in probability scale when a bound on its error from underflow is below
#' \eqn{10^{-12}} of every probability needed, and in log scale otherwise (see
#' \code{\link{cfisher.null}}). Apart from underflow, the probabilities carry the rounding
#' error of double-precision arithmetic, a relative error that grows with the number of
#' subjects. The p-values and the tilt are computed from the log probabilities, so a tiny
#' p-value keeps that precision down to about \eqn{10^{-308}}, below which it is returned
#' as 0 or with reduced precision. The distributions that the refit computes at each
#' hypothesized odds ratio stay in probability scale, which is fast; the inversion uses
#' only their tail probabilities at the level \eqn{\alpha/2} and their means, for which
#' the far tail that probability scale can lose is negligible.
#'
#' \emph{Assumptions.} Given the confounders, the outcomes are independent, and under the
#' null hypothesis the outcome does not depend on the exposure. The fitted null
#' probabilities are plugged in as if known, so with confounders the test is exact at the
#' true probabilities and approximate at the fitted ones. The conditional distribution does
#' not change when every odds is multiplied by a constant, so only the slopes of the null
#' model matter to it, not its intercept. The exposure must take both values and the outcome
#' too; otherwise the function stops with a message. A supplied probability of exactly 0
#' or 1 is an error. The function warns when some fitted probability is within
#' \eqn{10^{-8}} of 0 or 1, which suggests separation of the outcome by the confounders, as
#' when a stratum has no cases or no controls. The subjects concerned are then nearly
#' certain cases or controls and carry almost no information on the exposure. Under
#' complete separation every subject is such, the null distribution of \eqn{T} is
#' concentrated at the observed count, and the p-value carries no information. For the odds
#' ratio, under separation the refitted model has no maximum likelihood estimate at any
#' odds ratio, so the estimate and the interval are \code{NA} with a warning. A stratum
#' without cases or without controls carries no information on the exposure, and the test
#' on the data without such strata gives the odds ratio.
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' z <- rnorm(n)                                # a confounder
#' x <- rbinom(n, 1, plogis(-1 + 1.2 * z))      # exposure depends on z
#' y <- rbinom(n, 1, plogis(-0.5 + z))          # outcome depends on z, not on x
#' fisher.test(table(x, y))$p.value             # ignores z
#' cfisher.test(x, y, z)                        # controls for z
#'
#' # fit the null model once and test many exposures
#' fit <- glm(y ~ z, family = binomial)
#' cfisher.test(x, y, fit = fit, alternative = "greater", conf.int = FALSE)$p.value
#'
#' @references
#' Wu, Z. and Contador, G. (2026). Fisher's exact test controlling for confounders. Draft
#' manuscript.
#'
#' Chen, S. X. and Liu, J. S. (1997). Statistical applications of the Poisson-binomial and
#' conditional Bernoulli distributions. \emph{Statistica Sinica}, 7(4), 875-892.
#'
#' Sofer, T. (2017). BinomiRare: A robust test of the association of a rare variant with a
#' disease for pooled analysis and meta-analysis, with application to the HCHS/SOL.
#' \emph{Genetic Epidemiology}, 41(5), 388-395. \doi{10.1002/gepi.22044}
#'
#' Berrett, T. B., Wang, Y., Barber, R. F. and Samworth, R. J. (2020). The conditional
#' permutation test for independence while controlling for confounders. \emph{Journal of
#' the Royal Statistical Society Series B}, 82(1), 175-197. \doi{10.1111/rssb.12340}
#'
#' Candes, E., Fan, Y., Janson, L. and Lv, J. (2018). Panning for gold:
#' \sQuote{model-X} knockoffs for high dimensional controlled variable selection. \emph{Journal of
#' the Royal Statistical Society Series B}, 80(3), 551-577. \doi{10.1111/rssb.12265}
#'
#' @seealso \code{\link{cfisher.null}}, \code{\link[stats]{fisher.test}}
#' @export
cfisher.test <- function(x, y, z = NULL, alternative = c("two.sided", "greater", "less"),
                         conf.int = TRUE, conf.level = 0.95, midp = FALSE,
                         conditional = TRUE, fit = NULL) {
  # names for the printout, taken before the arguments are evaluated further
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  if (!is.null(fit)) {
    dname <- paste0(dname, ", controlling for the supplied null fit")
  } else if (!is.null(z)) {
    dname <- paste0(dname, ", controlling for ", deparse1(substitute(z)))
  }
  alternative <- match.arg(alternative)
  if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int))
    stop("'conf.int' must be TRUE or FALSE.")
  if (!is.logical(midp) || length(midp) != 1L || is.na(midp))
    stop("'midp' must be TRUE or FALSE.")
  if (conf.int && (!is.numeric(conf.level) || length(conf.level) != 1L ||
                   is.na(conf.level) || conf.level <= 0 || conf.level >= 1))
    stop("'conf.level' must be a single number strictly between 0 and 1.")
  # the root search resolves tail probabilities down to about 1e-12 (code review 2)
  if (conf.int) {
    a_tail <- if (alternative == "two.sided") (1 - conf.level) / 2 else 1 - conf.level
    if (min(a_tail, 1 - a_tail) < 1e-12)
      stop("'conf.level' must leave each tail probability of the interval between 1e-12 ",
           "and 1 - 1e-12.")
  }

  d <- .prepare(x, y, z, fit)
  nd <- .null_distribution(d, conditional)
  p <- .discrete_pvalue(d$t, nd$support, nd$logprob, alternative, midp)
  # An exposure in the span of the null model's design is absorbed by the null fit, so
  # its effect is not identified (review round 2026-11-04, D4). A vector fit has no design.
  unidentified <- .unidentified(d, z, fit)
  if (unidentified)
    warning("The exposure 'x' is a linear combination of the intercept and the confounders, ",
            "so its effect is not identified and the p-value carries no information about it",
            if (conf.int) "; the estimate and the interval are NA" else "", ".", call. = FALSE)

  res <- list(statistic = c(T = d$t),
              parameter = c(exposed = d$m, cases = d$r, n = d$n),
              p.value = p)
  if (conf.int) {
    oi <- .odds_ratio_inference(d, nd, z, fit, conditional, alternative, conf.level,
                                unidentified)
    res$conf.int <- oi$conf.int
    res$estimate <- c(oi$estimate)
    names(res$estimate) <- if (conditional) "conditional odds ratio" else "odds ratio"
  }
  res$null.value <- c("odds ratio" = 1)
  res$alternative <- alternative
  res$method <- paste0(if (conditional) "Fisher's exact test controlling for confounders"
                       else "Exact unconditional test controlling for confounders",
                       if (midp) " (mid-p value)" else "",
                       if (conf.int) paste0("; odds ratio: ", oi$method) else "")
  res$data.name <- dname
  res$null.distribution <- data.frame(support = nd$support, prob = nd$prob)
  res$expected <- sum(nd$support * nd$prob)
  class(res) <- "htest"
  res
}

#################### INTERNAL HELPER FUNCTIONS ####################

# Relative tolerance for treating two probabilities as tied in the two-sided p-value; the
# value used by stats::fisher.test (its relErr).
.TIE_TOL <- 1e-7

#' P-values on a discrete pmf given in log scale
#'
#' Right P(T >= t), left P(T <= t) and two-sided sum of {q(s): q(s) <= q(t) (1 + tol)}
#' p-values on a pmf given as a support and log probabilities, with an optional mid-p
#' correction that subtracts half the probability of the observed point. The sums are taken
#' in log scale (.lse), so a p-value keeps its relative precision however small the
#' probabilities are, until it is exponentiated. The result is clamped to [0, 1]. The
#' observed t must be in the support.
#' @keywords internal
#' @noRd
.discrete_pvalue <- function(t, support, logprob, alternative = "two.sided", midp = FALSE) {
  idx <- match(t, support)
  if (is.na(idx)) stop("The observed count ", t, " is not in the support.")
  lp_t <- logprob[idx]
  lp <- switch(alternative,
               two.sided = .lse(logprob[logprob <= lp_t + log1p(.TIE_TOL)]),
               greater   = .lse(logprob[support >= t]),
               less      = .lse(logprob[support <= t]))
  p <- exp(lp)
  if (midp) p <- p - 0.5 * exp(lp_t)
  min(max(p, 0), 1)
}

#' Exponential tilt q(t) exp(beta t) / normalizer, from the log probabilities
#'
#' The exponent log q(t) + beta t is shifted by its maximum before exponentiation, so that
#' a large |beta| of either sign gives neither overflow nor NaN. Because the log
#' probabilities are accurate in the far tails, the tilt is accurate there too.
#' @keywords internal
#' @noRd
.tilt <- function(logprob, support, beta) {
  lw <- logprob + beta * support
  w <- exp(lw - max(lw))
  w / sum(w)
}

# The search for a root of the log odds ratio stops at +-.BETA_MAX.
.BETA_MAX <- 40
# Tolerance of the pivoted QR that selects the columns of the null design: that of glm.fit
# (min(1e-7, epsilon / 1000) at its default epsilon 1e-8), so the refit drops no column
# that the null fit keeps.
.QR_TOL <- 1e-11
# A function value within .F_TOL of zero is not taken as a sign: the bracket search expands
# past it, and at the boundary of the search range it marks an asymptote, not a root.
.F_TOL <- 1e-9

#' Root of an increasing function f of the log odds ratio by an expanding bracket
#'
#' Starts from [-1, 1] and doubles an end, capped at -beta.max and beta.max, until f is
#' below -tol_f at the lower end and above tol_f at the upper end, so that rounding near an
#' asymptote is never taken for a sign change; then returns the root by uniroot (tolerance
#' 1e-8). When the lower end reaches -beta.max with f not below -tol_f it returns -Inf, and
#' when the upper end reaches beta.max with f not above tol_f it returns Inf. The returned
#' infinity carries the attribute beyond: TRUE when f at the boundary still has the wrong
#' sign by more than tol_f (a root, if any, lies beyond the search range), FALSE when f is
#' within tol_f of zero there (f approaches zero: an asymptote). Stops with a message when f
#' is not finite.
#' @keywords internal
#' @noRd
.solve_increasing <- function(f, beta.max = .BETA_MAX, tol_f = .F_TOL) {
  g <- function(b) {
    v <- f(b)
    if (!is.finite(v)) stop("The tail probability is not finite at log odds ratio ", b, ".")
    v
  }
  # the cap makes the last bracket end exactly +-beta.max, which is evaluated (D3)
  lo <- -1; f_lo <- g(lo)
  while (f_lo >= -tol_f) {
    if (lo <= -beta.max) return(structure(-Inf, beyond = f_lo > tol_f))
    lo <- max(2 * lo, -beta.max); f_lo <- g(lo)
  }
  hi <- 1; f_hi <- g(hi)
  while (f_hi <= tol_f) {
    if (hi >= beta.max) return(structure(Inf, beyond = f_hi < -tol_f))
    hi <- min(2 * hi, beta.max); f_hi <- g(hi)
  }
  stats::uniroot(g, c(lo, hi), f.lower = f_lo, f.upper = f_hi, tol = 1e-8)$root
}

#' Estimate and confidence interval of the odds ratio by inverting the test
#'
#' For a hypothesized log odds ratio b, dist_at(b) is the null distribution of T at that
#' odds ratio. With a refit (fit NULL or a glm) the confounder slopes are refitted by
#' maximum likelihood (.logistic_refit()) on an orthonormal basis of the column space of
#' the null model's design, with offset b * x, started from the nearest refit already
#' computed, and the conditional (or Poisson-binomial) pmf is computed at the refitted
#' log-odds. The refit at b = 0 runs first, from the null fit, with the strict convergence
#' test: the maximum likelihood estimate exists at every b if and only if it exists at
#' b = 0. A failure there (separation) makes the estimate and the interval NA with a
#' warning. The tolerance of the limits' root search is min(.F_TOL, 1e-3 min(a, 1 - a)).
#' With a
#' numeric vector fit, the pmf nd at those probabilities is tilted by exp(b t).
#' The estimate solves E_b(T) = t, the lower limit P_b(T >= t) = a and the upper limit
#' P_b(T <= t) = a, with a half of 1 - conf.level (two-sided) or 1 - conf.level. At an end
#' of the support the estimate and the limit on that side are infinite (top: estimate and
#' upper limit Inf; bottom: estimate and lower limit 0). A limit whose root is beyond
#' +-.BETA_MAX is reported as 0 or Inf. For an interior t the estimate is 0 or Inf when
#' E_b(T) - t is within .F_TOL of zero at the boundary (an asymptote, as when the model
#' with the exposure is separated along it), and NA with a warning when it is clearly of
#' the wrong sign there, except with the unconditional refit: there E_b(T) = t is the score
#' equation of the logistic model with the exposure, so the estimate is that model's
#' maximum likelihood estimate (.full_mle()), finite or, under separation along x, 0 or Inf. A quantity whose search fails (a refit that does not converge
#' included) is NA, with one warning at the end; so are all three when the effect is
#' unidentified. A second warning names the quantities whose evaluated function values,
#' including those at +-.BETA_MAX, show more than one crossing, because the refitted family
#' need not be monotone in b.
#' Returns list(estimate, conf.int (odds ratio scale, attribute conf.level), method).
#' @keywords internal
#' @noRd
.odds_ratio_inference <- function(d, nd, z, fit, conditional, alternative, conf.level,
                                  unidentified = FALSE) {
  t_obs <- d$t
  lo_s <- min(nd$support); hi_s <- max(nd$support)       # the theoretical support of T
  refit <- is.null(fit) || inherits(fit, "glm")
  method <- if (refit) "slopes refitted under each hypothesized odds ratio" else
    "tilt at the supplied odds"
  na_result <- function() {
    ci <- c(NA_real_, NA_real_); attr(ci, "conf.level") <- conf.level
    list(estimate = NA_real_, conf.int = ci, method = method)
  }
  if (unidentified) return(na_result())
  cache <- new.env()                                      # b -> distribution, shared by the roots
  if (refit) {
    if (inherits(fit, "glm")) {
      if (!identical(fit$family$link, "logit"))
        stop("The odds ratio needs a logistic null model: 'fit' has the link '", fit$family$link,
             "'. Use a logit link, or set conf.int = FALSE.")
      # the refit uses the model matrix only, so weights and offsets would be dropped (D5)
      if (any(fit$prior.weights != 1))
        stop("The odds ratio refit does not support prior weights in 'fit'. ",
             "Fit the null model without weights, or set conf.int = FALSE.")
      if (!is.null(fit$offset) && any(fit$offset != 0))
        stop("The odds ratio refit does not support an offset in 'fit'. ",
             "Fit the null model without an offset, or set conf.int = FALSE.")
    }
    Zm <- .nuisance_design(d, z, fit)
    if (nrow(Zm) != d$n) stop("The model matrix of 'fit' must have one row per subject.")
    # An orthonormal basis of the design's column space (pivoted QR at the tolerance
    # 1e-11 of glm.fit, so that the refit keeps every column the null fit keeps): the
    # same fitted values as the design, and a Hessian whose conditioning does not depend
    # on the scale or location of the confounders. The coefficients serve as warm starts.
    qz <- qr(Zm, tol = .QR_TOL)
    Zr <- qr.Q(qz)[, seq_len(qz$rank), drop = FALSE]
    # The strict test at b = 0, started from the null fit (its log-odds lie in the span
    # of Zr, because offsets are refused): whether the maximum likelihood estimate exists
    # does not depend on the offset, so this one check covers every b (separation)
    g0 <- .logistic_refit(Zr, d$y, offset = numeric(d$n),
                          start = drop(crossprod(Zr, d$eta)), strict = TRUE)
    if (!g0$converged) {
      warning("The refit of the null model did not converge, which indicates ",
              "(quasi-)separation of the outcome by the confounders, for example by a ",
              "stratum with no cases or no controls. The odds ratio is therefore not ",
              "computed: the estimate and the interval are NA. A stratum without cases or ",
              "without controls carries no information on the exposure, so the test on the ",
              "data without such strata gives the odds ratio.", call. = FALSE)
      return(na_result())
    }
    starts <- new.env()                                   # b -> converged coefficients
    starts[["0"]] <- g0$coef
    compute <- function(b) {
      # warm start: the coefficients of the converged refit nearest to b (D1)
      keys <- ls(starts)
      st <- starts[[keys[which.min(abs(as.numeric(keys) - b))]]]
      g <- .logistic_refit(Zr, d$y, offset = b * d$x, start = st, strict = FALSE)
      if (!g$converged)
        stop("the logistic refit did not converge at log odds ratio ", signif(b, 4))
      starts[[sprintf("%.17g", b)]] <- g$coef
      # probability scale is enough here: only tail sums and the mean are used
      nb <- if (conditional) .cb_pmf(g$eta, d$exposed, d$r, exact = FALSE)
            else .pb_pmf(g$eta[d$exposed], exact = FALSE)
      list(support = nb$support, prob = nb$prob)
    }
  } else {
    compute <- function(b)
      list(support = nd$support, prob = .tilt(nd$logprob, nd$support, b))
  }
  dist_at <- function(b) {
    key <- sprintf("%.17g", b)
    if (is.null(cache[[key]])) cache[[key]] <- compute(b)
    cache[[key]]
  }
  right_tail <- function(b) { v <- dist_at(b); sum(v$prob[v$support >= t_obs]) }
  left_tail  <- function(b) { v <- dist_at(b); sum(v$prob[v$support <= t_obs]) }
  mean_T     <- function(b) { v <- dist_at(b); sum(v$support * v$prob) - t_obs }

  failed <- character(0); solved <- list()
  # Solve for one quantity; root NA if the search fails. The function and the root are
  # kept for the crossing check below; beyond marks an infinite root whose function still
  # had the wrong sign at the boundary of the search range.
  solve_q <- function(f, name, tol_f = .F_TOL) {
    r <- tryCatch(.solve_increasing(f, tol_f = tol_f), error = function(e) {
      failed <<- c(failed, paste0(name, " (", conditionMessage(e), ")")); NA_real_ })
    beyond <- isTRUE(attr(r, "beyond")); root <- as.numeric(r)
    if (!is.na(root)) solved[[name]] <<- list(f = f, root = root, tol = tol_f / 10)
    list(root = root, beyond = beyond)
  }
  a <- if (alternative == "two.sided") (1 - conf.level) / 2 else 1 - conf.level
  # a tail function lies in [-a, 1 - a], so its tolerance must be well below min(a, 1 - a)
  tol_lim <- min(.F_TOL, 1e-3 * min(a, 1 - a))
  # the estimate: E_b(T) = t, increasing in b. For an interior t, an asymptote at the
  # boundary gives 0 or Inf, and a root clearly beyond the search range gives NA (D3)
  beta_hat <- if (t_obs >= hi_s) Inf else if (t_obs <= lo_s) -Inf else {
    sq <- solve_q(mean_T, "the estimate")
    # unconditional refit: E_b(T) = t is the score equation of the model with the
    # exposure, so its maximum likelihood estimate decides a root beyond the range
    bf <- if (sq$beyond && refit && !conditional)
      .full_mle(Zr, d, sign(sq$root), starts[[sprintf("%.17g", sign(sq$root) * .BETA_MAX)]])
    else NA_real_
    if (sq$beyond && !is.na(bf)) {
      solved[["the estimate"]] <- NULL
      bf
    } else if (sq$beyond) {
      failed <- c(failed, paste0("the estimate (the mean of T stays ",
                                 if (sq$root > 0) "below" else "above",
                                 " the observed count over the search range exp(+-", .BETA_MAX,
                                 "), so the estimate is beyond that range or infinite)"))
      solved[["the estimate"]] <- NULL
      NA_real_
    } else sq$root
  }
  # lower limit: P_b(T >= t) = a increases in b; upper: P_b(T <= t) = a decreases in b
  lower <- if (alternative == "less" || t_obs <= lo_s) -Inf else
    solve_q(function(b) right_tail(b) - a, "the lower limit", tol_lim)$root
  upper <- if (alternative == "greater" || t_obs >= hi_s) Inf else
    solve_q(function(b) a - left_tail(b), "the upper limit", tol_lim)$root
  if (length(failed))
    warning("The root search failed for ", paste(failed, collapse = "; "),
            "; they are returned as NA.", call. = FALSE)
  # Crossing check (suspicion S1 of the review round 2026-11-04): on every b evaluated by
  # any of the searches, and at +-.BETA_MAX, an increasing f must be negative below its
  # root and positive above it. A value of the wrong sign shows a second crossing. A refit
  # that fails at +-.BETA_MAX only leaves that point out of the check.
  if (length(solved))
    for (bb in c(-.BETA_MAX, .BETA_MAX)) tryCatch(dist_at(bb), error = function(e) NULL)
  bs <- as.numeric(ls(cache))
  multiple <- character(0)
  for (name in names(solved))
    if (.second_crossing(solved[[name]]$f, solved[[name]]$root, bs, solved[[name]]$tol))
      multiple <- c(multiple, name)
  if (length(multiple))
    warning("The refitted ", if (length(multiple) > 1) "functions are" else "function is",
            " not monotone in the odds ratio for ", paste(multiple, collapse = " and "),
            ": the values evaluated cross the target more than once, so the set of odds ",
            "ratios the test does not reject may not be an interval. The value reported is ",
            "the crossing found by the bracket search.", call. = FALSE)
  ci <- exp(c(lower, upper))
  attr(ci, "conf.level") <- conf.level
  list(estimate = exp(beta_hat), conf.int = ci, method = method)
}

#' Logistic maximum likelihood with an offset, by damped Newton
#'
#' Minimizes the negative log-likelihood sum(log(1 + exp(eta)) - y eta), eta = Z coef +
#' offset, for a design Z of full column rank. A Newton step is scaled so that it changes
#' no log-odds by more than 5, and while the Newton decrement g' H^-1 g is at least 1e-8,
#' or the step is larger than that, it is halved until it decreases the objective by at
#' least a fraction 1e-4 of the decrement (Armijo backtracking). The objective therefore
#' never increases, unlike the undamped iteration of glm.fit, which can report convergence
#' at diverging coefficients when the offset is large, and a step into a region where the
#' curvature vanishes cannot be followed by an unbounded one. Below a decrement of 1e-8,
#' with a step of at most 5 in the log-odds, the full step is taken, because the objective
#' cannot resolve so small a decrease in double precision. Convergence is declared after
#' such a full step when it is negligible in the coefficients, every |step_j| at most
#' 1e-10 (1 + |coef_j|) (strict = TRUE), or, with strict = FALSE, also when it changes no
#' fitted probability by more than 1e-12 to first order (max |mu (1 - mu) delta eta| <=
#' 1e-12), which accepts a coefficient pinned only by fitted probabilities near 0 or 1.
#' Under separation the full steps keep growing the coefficients, so the strict test does
#' not pass; the relaxed one can, which is why the caller runs the strict test once, at
#' b = 0 (whether the maximum likelihood estimate exists does not depend on the offset: the
#' objective grows along a direction d at the asymptotic rate sum over cases of
#' max(-(Z d)_i, 0) plus sum over controls of max((Z d)_i, 0), which the offset does not
#' change). Starts from start, or from zero. Returns list(converged, coef, eta); converged
#' is FALSE after maxit iterations without convergence, when the line search fails, or
#' when the Hessian is singular.
#' @keywords internal
#' @noRd
.logistic_refit <- function(Z, y, offset, start = NULL, strict = TRUE, maxit = 200L) {
  nll <- function(eta) sum(pmax(eta, 0) + log1p(exp(-abs(eta))) - y * eta)  # stable form
  cf <- if (is.null(start)) numeric(ncol(Z)) else start
  eta <- drop(Z %*% cf) + offset
  f <- nll(eta)
  for (it in seq_len(maxit)) {
    mu <- stats::plogis(eta)
    w <- mu * stats::plogis(-eta)                         # mu (1 - mu), without 1 - mu
    g <- drop(crossprod(Z, mu - y))                       # gradient
    step <- tryCatch(solve(crossprod(Z, Z * w), g), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) break     # singular Hessian
    dec <- sum(g * step)                                  # Newton decrement
    deta <- drop(Z %*% step)                              # change of eta by the full step
    big <- max(abs(deta))
    if (dec < 1e-8 && big <= 5) {
      # near the optimum: the full Newton step, without the line search
      cf <- cf - step
      eta <- drop(Z %*% cf) + offset
      f <- nll(eta)
      if (all(abs(step) <= 1e-10 * (1 + abs(cf))) || (!strict && max(w * abs(deta)) <= 1e-12))
        return(list(converged = TRUE, coef = cf, eta = eta))
      next
    }
    s0 <- min(1, 5 / big)                                 # at most 5 units of log-odds
    s <- s0
    repeat {                                              # Armijo backtracking
      cf_new <- cf - s * step
      eta_new <- drop(Z %*% cf_new) + offset
      f_new <- nll(eta_new)
      if (is.finite(f_new) && f_new <= f - 1e-4 * s * dec) break
      s <- s / 2
      if (s < 1e-10 * s0) break
    }
    if (s < 1e-10 * s0) break                             # no decrease far from the optimum
    cf <- cf_new; eta <- eta_new; f <- f_new
  }
  list(converged = FALSE, coef = cf, eta = eta)
}

#' The maximum likelihood estimate of the exposure's coefficient in the logistic model
#' with design (Zr, x), for the unconditional estimate beyond the search range
#'
#' With the refitted nuisance theta(b), the score of that model in its x coefficient is
#' t - E_b(T) under the unconditional null, so the root of E_b(T) = t is its maximum
#' likelihood estimate. Returns that estimate when the strict refit converges. When it
#' does not, the model is separated along x (the null model is not separated, and x is not
#' in the span of Zr, so every separating direction has an x coefficient of one sign), and
#' the function returns Inf or -Inf by the sign of the x coefficient, provided that it
#' agrees with direction, the side where the bracket search ran out; NA otherwise.
#' @keywords internal
#' @noRd
.full_mle <- function(Zr, d, direction, start = NULL) {
  # started from the refit at the boundary of the search range, with its log odds ratio
  st <- if (is.null(start)) NULL else c(start, direction * .BETA_MAX)
  g <- .logistic_refit(cbind(Zr, d$x), d$y, offset = numeric(d$n), start = st,
                       strict = TRUE, maxit = 1000L)
  bx <- g$coef[ncol(Zr) + 1L]
  if (g$converged) return(unname(bx))
  if (is.finite(bx) && sign(bx) == direction) return(direction * Inf)
  NA_real_
}

#' Whether an increasing function f, with root found at root, shows a second crossing on
#' the points bs: a value above eps at a point more than 1e-6 below the root, or below
#' -eps at a point more than 1e-6 above it (root may be -Inf or Inf); eps is a tenth of
#' the tolerance of the root search
#' @keywords internal
#' @noRd
.second_crossing <- function(f, root, bs, eps = 1e-10) {
  v <- vapply(bs, f, numeric(1))
  any((bs < root - 1e-6 & v > eps) | (bs > root + 1e-6 & v < -eps))
}

#' Design matrix of the null model for the refit: the model matrix of a glm fit, or the
#' confounders z with an intercept
#' @keywords internal
#' @noRd
.nuisance_design <- function(d, z, fit) {
  if (inherits(fit, "glm")) {
    Zm <- tryCatch(stats::model.matrix(fit), error = function(e)
      stop("The odds ratio refit needs the model matrix of 'fit', which could not be ",
           "rebuilt (", conditionMessage(e), "). Keep the data of the glm call available, ",
           "or set conf.int = FALSE.", call. = FALSE))
    # the columns the fit used: those with a coefficient that is not aliased (NA)
    cf <- stats::coef(fit)
    if (length(cf) == ncol(Zm)) Zm <- Zm[, !is.na(cf), drop = FALSE]
    return(Zm)
  }
  .design_matrix(z, d$n)
}

#' Whether the exposure x lies in the column space of the null model's design (its
#' effect is then not identified). For a glm fit, from the QR decomposition the fit
#' stores (of the weighted model matrix; the span does not depend on positive weights),
#' so that no model matrix is rebuilt; for z, by the ranks of the design with and
#' without x; FALSE for a numeric vector fit, which has no design.
#' @keywords internal
#' @noRd
.unidentified <- function(d, z, fit) {
  if (inherits(fit, "glm")) {
    q <- fit$qr
    if (is.null(q) || NROW(q$qr) != d$n) return(FALSE)    # cannot check (zero weights)
    v <- sqrt(fit$weights) * d$x
    r <- qr.resid(q, v)
    return(sqrt(sum(r^2)) <= 1e-7 * sqrt(sum(v^2)))
  }
  if (!is.null(fit)) return(FALSE)
  .in_span(d$x, .design_matrix(z, d$n))
}

#' Whether x lies in the column space of Zm: adding x does not raise the rank
#' @keywords internal
#' @noRd
.in_span <- function(x, Zm) {
  qr(cbind(Zm, x), tol = .QR_TOL)$rank <= qr(Zm, tol = .QR_TOL)$rank
}
