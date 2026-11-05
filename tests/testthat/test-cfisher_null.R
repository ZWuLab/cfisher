# ==============================================================================
# File: test-cfisher_null.R
# Purpose: Tests of R/cfisher_null.R: the Poisson-binomial recursion, the conditional
#   Bernoulli pmf against closed forms and against enumeration, the centering shift, the
#   fit argument, and the input checks of cfisher.null().

test_that("the Poisson-binomial recursion matches dbinom with equal probabilities and sums to one", {
  # the internals take log-odds eta = qlogis(mu)
  f <- exp(.pb_log(qlogis(rep(0.3, 20))))
  expect_equal(f, dbinom(0:20, 20, 0.3), tolerance = 1e-12)
  set.seed(1)
  expect_equal(sum(exp(.pb_log(qlogis(runif(100))))), 1, tolerance = 1e-12)
  s <- .pb_scaled(qlogis(c(0.2, 0.7)), kmax = 1)       # truncated at kmax = 1
  expect_equal(s$f * exp(-s$logscale), c(0.8 * 0.3, 0.2 * 0.3 + 0.8 * 0.7))
})

test_that("the rescaling of the probability-scale recursion reconstructs log probabilities", {
  # 2000 fair coins truncated at kmax = 5: every retained entry falls below the double
  # range (P(T = 5) is about 1e-589), so the recursion must rescale
  s <- .pb_scaled(rep(0, 2000), kmax = 5)
  expect_gt(s$logscale, 1000)                          # the rescaling branch ran
  lf <- log(s$f) - s$logscale
  expect_lt(max_abs_diff(lf, dbinom(0:5, 2000, 0.5, log = TRUE)), 1e-9)
})

test_that("the log-scale recursion matches dbinom(log = TRUE) where probability scale underflows", {
  # 2000 fair coins: P(0) = 2^-2000, far below the double range. The error is measured as
  # the largest absolute error of a log probability (about the largest relative error of a
  # probability); the rounding error of 2000 updates in log scale is near 7e-11
  lf <- .pb_log_cpp(rep(log(0.5), 2000), rep(log(0.5), 2000), 2000)
  expect_lt(max_abs_diff(lf, dbinom(0:2000, 2000, 0.5, log = TRUE)), 1e-9)
  # unequal probabilities: agrees with the probability-scale recursion where that is exact
  set.seed(2)
  eta <- qlogis(runif(300))
  lp <- plogis(eta, log.p = TRUE); lq <- plogis(-eta, log.p = TRUE)
  expect_lt(max_abs_diff(.pb_log_cpp(lp, lq, 300), log(.pb_scaled(eta)$f)), 1e-10)
  expect_identical(.pb_log_cpp(lp, lq, 40), .pb_log_cpp(lp, lq, 300)[1:41])   # truncation
})

test_that(".pb_log keeps probability scale when the needed entries are safe, and switches otherwise", {
  eta <- rep(0, 2000)                        # fair coins
  ref <- dbinom(0:2000, 2000, 0.5, log = TRUE)
  mid <- .pb_log(eta, need = 900:1100)       # needed entries near the center: probability scale
  expect_lt(max_abs_diff(mid[901:1101], ref[901:1101]), 1e-10)
  expect_equal(mid[1], -Inf)                 # the far tail it did not need has underflowed
  all <- .pb_log(eta)                        # every entry needed: the log-scale recursion
  expect_lt(max_abs_diff(all, ref), 1e-9)
})

test_that(".cb_pmf gives every probability of the support in log scale, however small", {
  # the Codex example of 2026-11-02: n = 3307, m = 302, r = 305, equal odds (hypergeometric);
  # in probability scale the support points above 272 underflowed to zero
  n <- 3307; m <- 302; r <- 305
  a <- .cb_pmf(rep(qlogis(r / n), n), 1:m, r)
  expect_equal(a$support, 0:302)
  expect_lt(max_abs_diff(a$logprob, dhyper(0:302, m, n - m, r, log = TRUE)), 1e-9)
  expect_lt(min(a$logprob), log(1e-308))     # the check reaches below the double range
  # two odds values (Fisher's noncentral hypergeometric), odds ratio 1000: the lower tail
  # falls below the double range instead
  w <- c(rep(1000, m), rep(1, n - m)) * 0.001
  b <- .cb_pmf(log(w), 1:m, r)
  ref <- log_nchg(n, m, r, psi = 1000)
  expect_equal(b$support, ref$support)
  expect_lt(max_abs_diff(b$logprob, ref$logprob), 1e-9)
  expect_lt(min(b$logprob), log(1e-308))
})

test_that(".cb_pmf keeps a centered probability near 1 apart from 1 (review round 2026-11-04, D2)", {
  # The Codex example: centering moves the two unexposed probabilities 1e-16 to about
  # 1 - 1e-17, which rounds to 1 as a probability; on the log-odds scale it does not.
  # Closed form: P(T = 0) = w2 w3 / (w1 w2 + w1 w3 + w2 w3) with w the odds
  mu <- c(0.9, 1e-16, 1e-16)
  w <- mu / (1 - mu)
  d <- .cb_pmf(qlogis(mu), 1L, 2)
  expect_equal(d$support, 0:1)
  p0 <- w[2] * w[3] / (w[1] * w[2] + w[1] * w[3] + w[2] * w[3])
  expect_equal(d$logprob[1], log(p0), tolerance = 1e-12)
  expect_warning(p <- cfisher.test(c(1, 0, 0), c(0, 1, 1), fit = mu, alternative = "less",
                                   conf.int = FALSE)$p.value, "within 1e-8")
  expect_equal(p, p0, tolerance = 1e-12)
})

test_that(".cb_pmf with equal odds is the hypergeometric pmf", {
  n <- 300; m <- 12; r <- 100
  d <- .cb_pmf(rep(qlogis(0.37), n), exposed = 1:m, r = r)
  expect_equal(d$support, max(0, r - n + m):min(m, r))
  expect_equal(d$prob, dhyper(d$support, m, n - m, r), tolerance = 1e-12)
})

test_that(".cb_pmf with two odds values is Fisher's noncentral hypergeometric pmf", {
  skip_if_not_installed("MCMCpack")
  n <- 400; m <- 15; r <- 120
  mu <- c(rep(0.5, m), rep(0.25, n - m))      # odds 1 and 1/3, odds ratio 3
  d <- .cb_pmf(qlogis(mu), exposed = 1:m, r = r)
  ref <- sapply(d$support, function(x)
    MCMCpack::dnoncenhypergeom(x, n1 = r, n2 = n - r, m1 = m, psi = 3))
  expect_equal(d$prob, ref, tolerance = 1e-10)
})

test_that("both null distributions equal the enumeration of all outcome vectors for small n", {
  set.seed(5)
  for (n in c(6, 9, 11)) {
    mu <- runif(n, 0.1, 0.9)
    exposed <- sort(sample(n, 3))
    for (r in c(1, floor(n / 2), n - 1)) {
      ref <- enumerate_null(mu, exposed, r = r)
      d <- .cb_pmf(qlogis(mu), exposed, r)
      expect_equal(d$support, ref$support)
      expect_equal(d$prob, ref$prob, tolerance = 1e-12)
    }
    ref <- enumerate_null(mu, exposed)
    expect_equal(.pb_pmf(qlogis(mu[exposed]))$prob, ref$prob, tolerance = 1e-12)
  }
})

test_that(".cb_pmf centers the total: a common log-odds shift and a far-off total change nothing", {
  set.seed(9)
  n <- 300; eta <- rnorm(n, -1); exposed <- sample(n, 12); r <- 80
  a <- .cb_pmf(eta, exposed, r)
  b <- .cb_pmf(eta + 1.7, exposed, r)                   # shifted odds, same conditional pmf
  expect_equal(a$prob, b$prob, tolerance = 1e-12)
  d <- .cb_pmf(rep(0, 5000), 1:30, 400)                 # total 2,100 below the mean
  expect_equal(d$prob, dhyper(d$support, 30, 4970, 400), tolerance = 1e-12)
  expect_equal(.cb_pmf(eta, exposed, 0), list(support = 0L, prob = 1, logprob = 0))
  expect_equal(.cb_pmf(eta, exposed, n), list(support = 12L, prob = 1, logprob = 0))
  # few subjects and a total far from sum(mu): the Newton step saturates and the bracketing
  # fallback must take over; equal odds, so the answer is the hypergeometric pmf
  e <- .cb_pmf(rep(qlogis(0.1), 10), 1:4, 9)
  expect_equal(e$prob, dhyper(e$support, 4, 6, 9), tolerance = 1e-12)
  f <- .cb_pmf(rep(qlogis(0.9), 10), 1:4, 1)
  expect_equal(f$prob, dhyper(f$support, 4, 6, 1), tolerance = 1e-12)
})

test_that("cfisher.null without confounders is hypergeometric (conditional) or binomial (unconditional)", {
  set.seed(21)
  n <- 120; x <- rbinom(n, 1, 0.2); y <- rbinom(n, 1, 0.35)
  m <- sum(x); r <- sum(y)
  d <- cfisher.null(x, y)
  expect_s3_class(d, "data.frame")
  expect_named(d, c("support", "prob"))
  expect_equal(d$prob, dhyper(d$support, m, n - m, r), tolerance = 1e-12)
  u <- cfisher.null(x, y, conditional = FALSE)
  expect_equal(u$support, 0:m)
  expect_equal(u$prob, dbinom(0:m, m, r / n), tolerance = 1e-12)
})

test_that("cfisher.null with confounders equals the enumeration at the fitted probabilities", {
  set.seed(6)
  n <- 10; z <- rnorm(n); y <- rep(0:1, c(6, 4))[sample(n)]; x <- c(1, 1, 1, rep(0, n - 3))[sample(n)]
  mu <- fitted(glm(y ~ z, family = binomial))
  exposed <- which(x == 1)
  d <- cfisher.null(x, y, z)
  ref <- enumerate_null(mu, exposed, r = sum(y))
  expect_equal(d$support, ref$support)
  expect_equal(d$prob, ref$prob, tolerance = 1e-10)
  u <- cfisher.null(x, y, z, conditional = FALSE)
  expect_equal(u$prob, enumerate_null(mu, exposed)$prob, tolerance = 1e-10)
})

test_that("fit as a glm object, fit as a vector, z as a matrix and z as a data frame give the same null distribution", {
  set.seed(3)
  n <- 800; z1 <- rnorm(n); z2 <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + z1 + z2)); x <- rbinom(n, 1, 0.02)
  fit <- glm(y ~ z1 + z2, family = binomial)
  a <- cfisher.null(x, y, z = cbind(z1, z2))
  b <- cfisher.null(x, y, fit = fit)
  c_ <- cfisher.null(x, y, fit = fitted(fit))
  d <- cfisher.null(x, y, z = data.frame(z1 = z1, z2 = factor(z2)))
  expect_equal(a, b, tolerance = 1e-10)
  expect_equal(a, c_, tolerance = 1e-10)
  expect_equal(a, d, tolerance = 1e-10)
  # fit takes precedence over z
  expect_equal(cfisher.null(x, y, z = rnorm(n), fit = fit), a, tolerance = 1e-10)
})

test_that("the input checks stop with informative messages", {
  set.seed(8)
  y <- rbinom(50, 1, 0.5); x <- rbinom(50, 1, 0.3)
  expect_error(cfisher.null(rep(0, 50), y), "constant")
  expect_error(cfisher.null(rep(1, 50), y), "constant")
  expect_error(cfisher.null(x, rep(1, 50)), "constant")
  expect_error(cfisher.null(c(x[-1], 2), y), "0/1")
  expect_error(cfisher.null(x, c(y[-1], NA)), "0/1")
  expect_error(cfisher.null(x, c(y[-1], 0.5)), "0/1")
  expect_error(cfisher.null(x[-1], y), "same length")
  expect_error(cfisher.null(x, y, fit = rep(c(0.3, 1), 25)), "strictly inside")
  expect_error(cfisher.null(x, y, fit = rep(0.3, 49)), "one fitted probability")
  expect_error(cfisher.null(x, y, fit = glm(y ~ x, family = poisson)), "binomial")
  expect_error(cfisher.null(x, y, z = rnorm(49)), "one row per subject")
  expect_error(cfisher.null(x, y, conditional = NA), "TRUE or FALSE")
  expect_error(.cb_pmf(rep(0, 5), 1:5, 2), "proper subset")
})
