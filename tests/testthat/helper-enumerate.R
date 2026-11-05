# ==============================================================================
# File: helper-enumerate.R
# Purpose: Test helpers: the exact null distribution of T by enumerating all 2^n outcome
#   vectors (a slow, obviously correct reference), and a 2 x 2 table builder.

# The null pmf of T = sum of y over `exposed`, for independent Bernoulli(mu) outcomes,
# by summing the probability of every outcome vector. If r is given, conditions on
# sum(y) = r. Returns list(support, prob) over the values with positive probability.
enumerate_null <- function(mu, exposed, r = NULL) {
  n <- length(mu)
  ys <- as.matrix(expand.grid(rep(list(0:1), n)))          # 2^n x n, all outcome vectors
  logp <- ys %*% log(mu) + (1 - ys) %*% log(1 - mu)         # log P(y)
  p <- exp(as.numeric(logp))
  tt <- rowSums(ys[, exposed, drop = FALSE])
  keep <- if (is.null(r)) rep(TRUE, nrow(ys)) else rowSums(ys) == r
  pm <- tapply(p[keep], tt[keep], sum)
  list(support = as.integer(names(pm)), prob = as.numeric(pm) / sum(pm))
}

# x and y vectors with a exposed cases, b exposed controls, c unexposed cases, d unexposed
# controls; the 2 x 2 table with rows exposed/unexposed and columns case/control.
make_xy <- function(a, b, c, d) {
  list(x = rep(c(1L, 1L, 0L, 0L), c(a, b, c, d)),
       y = rep(c(1L, 0L, 1L, 0L), c(a, b, c, d)),
       tab = matrix(c(a, c, b, d), 2))
}

# The null pmf of T at the log odds ratio b with the confounder slopes refitted (offset b * x),
# conditional or Poisson-binomial: a restatement of the construction of cfisher.test's
# estimate and interval, from the same internals.
refit_dist <- function(x, y, z, b, conditional = TRUE) {
  g <- cfisher:::.logistic_refit(cbind(1, z), y, offset = b * x)
  stopifnot(g$converged)
  exposed <- which(x == 1)
  if (conditional) cfisher:::.cb_pmf(g$eta, exposed, sum(y))
  else cfisher:::.pb_pmf(g$eta[exposed])
}

# Largest absolute difference between two vectors, counting equal infinite entries (a log
# probability of -Inf in both) as no difference: the precision tests measure the worst
# entry, not the mean relative difference that expect_equal() uses.
max_abs_diff <- function(a, b) {
  stopifnot(length(a) == length(b))
  same <- a == b                                       # covers -Inf == -Inf
  if (all(same)) return(0)
  max(abs(a[!same] - b[!same]))
}

# Log pmf of Fisher's noncentral hypergeometric distribution of T, from the closed form
# choose(m, t) choose(n - m, r - t) psi^t / normalizer, in log scale (lchoose, log-sum-exp),
# so that it is accurate in the far tails: n subjects, m exposed, r cases, odds ratio psi.
log_nchg <- function(n, m, r, psi = 1) {
  s <- max(0, r - n + m):min(m, r)
  lw <- lchoose(m, s) + lchoose(n - m, r - s) + s * log(psi)
  list(support = s, logprob = lw - (max(lw) + log(sum(exp(lw - max(lw))))))
}

# Estimate and 95% limits of the odds ratio without confounders from the hypergeometric pmf
# in log scale (dhyper(log = TRUE), log-sum-exp), by tail inversion: the reference for the
# far-tail tests. n subjects, m exposed, r cases, t exposed cases.
hyper_log_reference <- function(n, m, r, t, a = 0.025) {
  s <- max(0, r - n + m):min(m, r)
  lp <- dhyper(s, m, n - m, r, log = TRUE)
  lse <- function(v) max(v) + log(sum(exp(v - max(v))))
  logtail <- function(b, up) { lw <- lp + b * s; keep <- if (up) s <= t else s >= t
    lse(lw[keep]) - lse(lw) }
  est <- uniroot(function(b) { lw <- lp + b * s - lse(lp + b * s); sum(s * exp(lw)) - t },
                 c(-30, 30), tol = 1e-12)$root
  lo <- uniroot(function(b) exp(logtail(b, FALSE)) - a, c(-30, 30), tol = 1e-12)$root
  up <- uniroot(function(b) a - exp(logtail(b, TRUE)), c(-30, 30), tol = 1e-12)$root
  list(estimate = exp(est), lower = exp(lo), upper = exp(up))
}
