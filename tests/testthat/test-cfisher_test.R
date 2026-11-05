# ==============================================================================
# File: test-cfisher_test.R
# Purpose: Tests of R/cfisher_test.R: agreement with stats::fisher.test without
#   confounders (p-values on three alternatives, estimate, interval), the identities of the
#   p-values, the tie tolerance, invariance to a common scale of the odds, the estimate and
#   the interval, the unconditional option, and the input checks.

test_that("without confounders cfisher.test reproduces fisher.test on all three alternatives", {
  # (exposed cases, exposed controls, unexposed cases, unexposed controls); the 3rd and 4th
  # have an empty cell, so the estimate is 0 or Inf
  configs <- list(c(5, 15, 10, 30), c(12, 3, 20, 40), c(0, 8, 9, 25), c(7, 0, 12, 20),
                  c(30, 50, 40, 80), c(2, 2, 3, 3), c(20, 10, 5, 40))
  for (cfg in configs) {
    dd <- make_xy(cfg[1], cfg[2], cfg[3], cfg[4])
    for (alt in c("two.sided", "greater", "less")) {
      for (lev in c(0.95, 0.9)) {
        ft <- fisher.test(dd$tab, alternative = alt, conf.level = lev)
        cf <- cfisher.test(dd$x, dd$y, alternative = alt, conf.level = lev)
        lab <- paste(paste(cfg, collapse = ","), alt, lev)
        expect_equal(cf$p.value, ft$p.value, tolerance = 1e-6, label = paste("p-value", lab))
        expect_equal(unname(cf$estimate), unname(ft$estimate), tolerance = 1e-4,
                     label = paste("estimate", lab))
        # fisher.test solves the limits with uniroot's default tolerance, so its limits miss
        # the target tail probability in the third or fourth digit (see the test that
        # checks the tails of cfisher.test's limits at 1e-6); 2e-3 absorbs that
        expect_equal(as.numeric(cf$conf.int), as.numeric(ft$conf.int), tolerance = 2e-3,
                     label = paste("interval", lab))
        expect_equal(attr(cf$conf.int, "conf.level"), lev)
      }
    }
  }
})

test_that("without confounders cfisher.test agrees with fisher.test on random tables", {
  set.seed(11)
  for (k in 1:4) {
    n <- 500; y <- rbinom(n, 1, 0.3); x <- rbinom(n, 1, c(0.04, 0.1, 0.3, 0.5)[k])
    tab <- table(factor(x, levels = 1:0), factor(y, levels = 1:0))
    for (alt in c("two.sided", "greater", "less")) {
      cf <- cfisher.test(x, y, alternative = alt)
      ft <- fisher.test(tab, alternative = alt)
      lab <- paste("table", k, alt)
      expect_equal(cf$p.value, ft$p.value, tolerance = 1e-6, label = lab)
      expect_equal(unname(cf$estimate), unname(ft$estimate), tolerance = 1e-4, label = lab)
      expect_equal(as.numeric(cf$conf.int), as.numeric(ft$conf.int), tolerance = 2e-3, label = lab)
    }
  }
})

test_that("the returned htest object has the documented components", {
  dd <- make_xy(5, 15, 10, 30)
  cf <- cfisher.test(dd$x, dd$y)
  expect_s3_class(cf, "htest")
  expect_equal(unname(cf$statistic), 5)
  expect_named(cf$statistic, "T")
  expect_equal(unname(cf$parameter), c(20, 15, 60))
  expect_named(cf$parameter, c("exposed", "cases", "n"))
  expect_named(cf$estimate, "conditional odds ratio")
  expect_equal(cf$null.value, c("odds ratio" = 1))
  expect_equal(cf$alternative, "two.sided")
  expect_match(cf$method, "Fisher's exact test controlling for confounders")
  expect_match(cf$data.name, "dd\\$x and dd\\$y")
  expect_named(cf$null.distribution, c("support", "prob"))
  expect_equal(cf$expected, sum(cf$null.distribution$support * cf$null.distribution$prob))
  expect_equal(cf$expected, 20 * 15 / 60, tolerance = 1e-10)      # hypergeometric mean
  expect_output(print(cf), "alternative hypothesis")
  expect_match(cfisher.test(dd$x, dd$y, midp = TRUE)$method, "mid-p")
  # conf.int = FALSE omits the estimate and the interval
  nc <- cfisher.test(dd$x, dd$y, conf.int = FALSE)
  expect_null(nc$conf.int)
  expect_null(nc$estimate)
  expect_equal(nc$p.value, cf$p.value)
})

test_that("the estimate solves E(T) = t under the refitted distribution at the estimate", {
  set.seed(14)
  n <- 400; z <- rnorm(n); x <- rbinom(n, 1, plogis(-2 + 0.5 * z)); y <- rbinom(n, 1, plogis(-1 + z + 0.5 * x))
  for (cond in c(TRUE, FALSE)) {
    cf <- cfisher.test(x, y, z, conditional = cond)
    nd <- refit_dist(x, y, z, log(unname(cf$estimate)), cond)
    expect_equal(sum(nd$support * nd$prob), unname(cf$statistic), tolerance = 1e-6)
  }
  expect_named(cfisher.test(x, y, z, conditional = FALSE)$estimate, "odds ratio")
})

test_that("the sides satisfy their identities, and mid-p is the defining subtraction", {
  set.seed(2)
  n <- 1000; z <- rnorm(n); y <- rbinom(n, 1, plogis(-1.5 + z)); x <- rbinom(n, 1, 0.01)
  for (cond in c(TRUE, FALSE)) {
    r <- cfisher.test(x, y, z, alternative = "greater", conditional = cond)
    nd <- r$null.distribution
    p_t <- nd$prob[match(r$statistic, nd$support)]
    pl <- cfisher.test(x, y, z, alternative = "less", conditional = cond)$p.value
    expect_equal(pl + r$p.value, 1 + p_t)                                  # p_L + p_R = 1 + q(t)
    expect_equal(cfisher.test(x, y, z, alternative = "less", midp = TRUE, conditional = cond)$p.value +
                 cfisher.test(x, y, z, alternative = "greater", midp = TRUE, conditional = cond)$p.value,
                 1)                                                         # mid-p sides sum to 1
    for (s in c("two.sided", "greater", "less"))
      expect_equal(cfisher.test(x, y, z, alternative = s, midp = TRUE, conditional = cond)$p.value,
                   cfisher.test(x, y, z, alternative = s, conditional = cond)$p.value - p_t / 2,
                   label = s)
  }
})

test_that("the mid-p value leaves the estimate and the interval unchanged", {
  dd <- make_xy(8, 12, 15, 40)
  a <- cfisher.test(dd$x, dd$y)
  b <- cfisher.test(dd$x, dd$y, midp = TRUE)
  expect_lt(b$p.value, a$p.value)
  expect_equal(b$estimate, a$estimate)
  expect_equal(b$conf.int, a$conf.int)
})

test_that("the two-sided rule counts rounding ties as ties (fisher.test's tolerance)", {
  # n = 10, m = 2, r = 3, t = 0: the hypergeometric probabilities at 0 and 1 tie in exact
  # arithmetic up to rounding; fisher.test gives 1, a literal comparison would not
  p <- dhyper(0:2, 2, 8, 3)
  expect_equal(.discrete_pvalue(0, 0:2, log(p), "two.sided"), 1)
  expect_equal(.discrete_pvalue(0, 0:2, log(p), "two.sided"), fisher.test(matrix(c(0, 2, 3, 5), 2))$p.value)
  expect_error(.discrete_pvalue(7, 0:5, log(rep(1/6, 6))), "not in the support")
  # the same table through the exported function
  dd <- make_xy(0, 2, 3, 5)
  expect_equal(cfisher.test(dd$x, dd$y)$p.value, 1)
})

test_that("the conditional test does not change when every fitted odds is multiplied by a constant", {
  set.seed(17)
  n <- 500; z <- rnorm(n); y <- rbinom(n, 1, plogis(-1 + z)); x <- rbinom(n, 1, plogis(-3 + 0.5 * z))
  mu <- fitted(glm(y ~ z, family = binomial))
  base <- cfisher.test(x, y, fit = mu)
  for (shift in c(-2, 0.7, 3)) {
    scaled <- cfisher.test(x, y, fit = plogis(qlogis(mu) + shift))        # odds times exp(shift)
    expect_equal(scaled$p.value, base$p.value, tolerance = 1e-8)
    expect_equal(scaled$conf.int, base$conf.int, tolerance = 1e-6)
    expect_equal(scaled$estimate, base$estimate, tolerance = 1e-6)
  }
})

test_that("fit as a glm, fit as a vector and z give the same p-value; a glm and z the same interval", {
  set.seed(3)
  n <- 800; z1 <- rnorm(n); z2 <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + z1 + z2)); x <- rbinom(n, 1, 0.02)
  fit <- glm(y ~ z1 + z2, family = binomial)
  a <- cfisher.test(x, y, z = cbind(z1, z2))
  b <- cfisher.test(x, y, fit = fit)
  c_ <- cfisher.test(x, y, fit = fitted(fit))
  d <- cfisher.test(x, y, z = data.frame(z1 = z1, z2 = factor(z2)))
  expect_equal(a$p.value, b$p.value, tolerance = 1e-10)
  expect_equal(a$p.value, c_$p.value, tolerance = 1e-10)
  expect_equal(a$p.value, d$p.value, tolerance = 1e-10)
  expect_equal(a$conf.int, b$conf.int, tolerance = 1e-8)
  expect_equal(unname(a$parameter), c(sum(x), sum(y), n))
  expect_equal(unname(a$statistic), sum(y[x == 1]))
})

test_that("the interval covers the estimate, shrinks with the confidence level, and has the right tails", {
  set.seed(23)
  n <- 600; z <- rnorm(n); x <- rbinom(n, 1, plogis(-2 + z)); y <- rbinom(n, 1, plogis(-1 + z + 0.6 * x))
  est <- cfisher.test(x, y, z)$estimate
  ci90 <- cfisher.test(x, y, z, conf.level = 0.90)$conf.int
  ci95 <- cfisher.test(x, y, z, conf.level = 0.95)$conf.int
  ci99 <- cfisher.test(x, y, z, conf.level = 0.99)$conf.int
  expect_true(all(ci95[1] <= est, est <= ci95[2]))
  expect_true(ci99[1] < ci95[1] && ci95[1] < ci90[1])
  expect_true(ci90[2] < ci95[2] && ci95[2] < ci99[2])
  # the interval inverts the tails of the refitted distribution: alpha / 2 at each limit
  t_obs <- unname(cfisher.test(x, y, z)$statistic)
  tl <- function(b, upper) { nd <- refit_dist(x, y, z, b, TRUE)
    if (upper) sum(nd$prob[nd$support <= t_obs]) else sum(nd$prob[nd$support >= t_obs]) }
  expect_equal(tl(log(ci95[1]), FALSE), 0.025, tolerance = 1e-6)
  expect_equal(tl(log(ci95[2]), TRUE), 0.025, tolerance = 1e-6)
  # one-sided intervals have one infinite end
  g <- cfisher.test(x, y, z, alternative = "greater")$conf.int
  l <- cfisher.test(x, y, z, alternative = "less")$conf.int
  expect_equal(g[2], Inf)
  expect_equal(l[1], 0)
  expect_equal(tl(log(g[1]), FALSE), 0.05, tolerance = 1e-6)
  expect_equal(tl(log(l[2]), TRUE), 0.05, tolerance = 1e-6)
})

test_that("the estimate and limits are at 0 or Inf when the count is at an end of the support", {
  top <- make_xy(6, 0, 4, 20)        # every exposed subject is a case
  bot <- make_xy(0, 6, 8, 20)        # no exposed subject is a case
  a <- cfisher.test(top$x, top$y); b <- cfisher.test(bot$x, bot$y)
  expect_equal(unname(a$estimate), Inf); expect_equal(a$conf.int[2], Inf); expect_gt(a$conf.int[1], 0)
  expect_equal(unname(b$estimate), 0);   expect_equal(b$conf.int[1], 0);   expect_lt(b$conf.int[2], Inf)
})

test_that("a large effect gives finite estimate and limits without overflow", {
  big <- make_xy(40, 2, 5, 400)       # 40 of 42 exposed are cases; null mean about 4
  r <- cfisher.test(big$x, big$y)
  expect_true(all(is.finite(r$conf.int)) && is.finite(r$estimate))
  expect_true(r$conf.int[1] < r$estimate && r$estimate < r$conf.int[2] && r$conf.int[1] > 1)
})

test_that("an observed count far below the double range gives p-value 0 and exact limits, silently", {
  extreme <- make_xy(300, 2, 5, 3000)  # null probability of T = 300 is below the double range
  expect_silent(r <- cfisher.test(extreme$x, extreme$y))
  expect_equal(r$p.value, 0)           # the true p-value is below the double range as well
  ref <- hyper_log_reference(3307, 302, 305, 300)
  expect_equal(unname(r$estimate), ref$estimate, tolerance = 1e-6)
  expect_equal(as.numeric(r$conf.int), c(ref$lower, ref$upper), tolerance = 1e-6)
})

test_that("the unconditional option uses the Poisson-binomial null", {
  set.seed(4)
  n <- 600; z <- rnorm(n); y <- rbinom(n, 1, plogis(-1 + z)); x <- rbinom(n, 1, 0.03)
  fit <- glm(y ~ z, family = binomial)
  r <- cfisher.test(x, y, fit = fit, conditional = FALSE)
  m <- sum(x)
  expect_equal(r$null.distribution$support, 0:m)
  # reference: the Poisson-binomial pmf by repeated convolution, one subject at a time
  ref <- Reduce(function(f, p) c(f * (1 - p), 0) + c(0, f * p), fitted(fit)[x == 1], 1)
  expect_equal(r$null.distribution$prob, ref, tolerance = 1e-12)
  expect_equal(r$expected, sum(fitted(fit)[x == 1]), tolerance = 1e-10)
  expect_match(r$method, "unconditional")
})

test_that("bad inputs stop with a message", {
  set.seed(1)
  y <- rbinom(50, 1, 0.5); x <- rbinom(50, 1, 0.4)
  expect_error(cfisher.test(rep(0, 50), y), "constant")
  expect_error(cfisher.test(rep(1, 50), y), "constant")
  expect_error(cfisher.test(c(x[-1], 2), y), "0/1")
  expect_error(cfisher.test(x[-1], y), "same length")
  expect_error(cfisher.test(x, y, fit = c(rep(0.5, 49), 1)), "strictly inside")
  expect_error(cfisher.test(x, y, fit = c(rep(0.5, 49), 0)), "strictly inside")
  expect_error(cfisher.test(x, y, conf.level = 1.2), "conf.level")
  expect_error(cfisher.test(x, y, alternative = "bigger"), "should be one of")
  expect_error(cfisher.test(x, y, midp = NA), "midp")
})

test_that("D1: the refit recovers a common exposure's odds ratio, the null-fit tilt does not", {
  # Codex round 2026-11-02, D1: x = 1 for half the subjects and a strong confounder; the null
  # fit absorbs the exposure effect, so tilting its distribution estimates about 1.9 (not 4)
  set.seed(20261102)
  x <- rep(c(0L, 1L, 0L, 1L), c(900, 100, 100, 900)); z <- rep(c(0, 0, 1, 1), c(900, 100, 100, 900))
  y <- rbinom(length(x), 1, plogis(-2 + z + log(4) * x))
  refit <- cfisher.test(x, y, z)
  null_fit <- glm(y ~ z, family = binomial)
  tilt <- cfisher.test(x, y, fit = fitted(null_fit))
  oracle <- cfisher.test(x, y, fit = plogis(-2 + z))       # true baseline probabilities
  expect_match(refit$method, "slopes refitted under each hypothesized odds ratio")
  expect_match(tilt$method, "tilt at the supplied odds")
  expect_true(abs(log(refit$estimate) - log(4)) < 0.35)
  expect_true(refit$conf.int[1] < 4 && 4 < refit$conf.int[2])
  expect_lt(unname(tilt$estimate), 2.5)                    # the null-fit tilt is biased down
  expect_true(abs(log(oracle$estimate) - log(4)) < 0.35)   # exact tilt at the true baseline
  # the unconditional null gives the same kind of answer
  expect_true(abs(log(cfisher.test(x, y, z, conditional = FALSE)$estimate) - log(4)) < 0.35)
  # the deterministic grouped data set of the review (10,000 subjects, expected counts)
  xg <- rep(c(0L, 1L, 0L, 1L), c(4500, 500, 500, 4500)); zg <- rep(c(0, 0, 1, 1), c(4500, 500, 500, 4500))
  cases <- c(536, 176, 134, 2679)
  yg <- unlist(mapply(function(k, ng) c(rep(1L, k), rep(0L, ng - k)), cases, c(4500, 500, 500, 4500)))
  rg <- cfisher.test(xg, yg, zg)
  expect_true(abs(log(rg$estimate) - log(4)) < 0.1)
  expect_true(rg$conf.int[1] < 4 && 4 < rg$conf.int[2])
})

test_that("a supplied probability vector that is the baseline under the alternative is exact", {
  # a tilt of the distribution at the baseline probabilities: compare with the refit's own
  # tilt formula by hand on a small example where the pmf is enumerated
  mu <- c(0.2, 0.3, 0.4, 0.25, 0.35, 0.5, 0.15, 0.45)
  exposed <- c(1, 2, 3, 4); y <- c(1, 0, 1, 0, 1, 1, 0, 0); x <- as.integer(1:8 %in% exposed)
  r <- cfisher.test(x, y, fit = mu)
  nd <- enumerate_null(mu, exposed, r = sum(y))
  # estimate: E_b(T) = t under the tilt
  f <- function(b) { w <- nd$prob * exp(b * nd$support); sum(nd$support * w) / sum(w) - r$statistic }
  expect_equal(unname(log(r$estimate)), uniroot(f, c(-10, 10), tol = 1e-12)$root, tolerance = 1e-6)
})

test_that("D2: the far tail is computed, and refit and tilt both match a log-scale reference", {
  # Codex round 2026-11-02, D2: n = 3307, m = 302, r = 305, no confounders. In probability
  # scale the support above 272 underflowed; in log scale every point is kept.
  n <- 3307; m <- 302; r <- 305
  mk <- function(t) list(x = c(rep(1L, m), rep(0L, n - m)),
                         y = c(rep(1L, t), rep(0L, m - t), rep(1L, r - t), rep(0L, n - m - r + t)))
  expect_equal(max(cfisher.null(mk(270)$x, mk(270)$y)$support), 302)   # support kept
  for (t in c(262, 270, 272, 290)) {
    d <- mk(t)
    ref <- hyper_log_reference(n, m, r, t)
    refv <- c(ref$estimate, ref$lower, ref$upper)
    expect_silent(res <- cfisher.test(d$x, d$y))                  # refit
    expect_equal(c(unname(res$estimate), as.numeric(res$conf.int)), refv, tolerance = 1e-6,
                 label = paste("refit", t))
    expect_silent(tl <- cfisher.test(d$x, d$y, fit = rep(r / n, n)))  # tilt of the supplied odds
    expect_equal(c(unname(tl$estimate), as.numeric(tl$conf.int)), refv, tolerance = 1e-6,
                 label = paste("tilt", t))
  }
  # the reviewer's values at T = 270, which the tilt in probability scale got wrong (957, 4839)
  ref <- hyper_log_reference(n, m, r, 270)
  expect_equal(c(ref$estimate, ref$lower, ref$upper), c(706.0552, 424.8078, 1212.8350), tolerance = 1e-6)
})

test_that("an extreme observed count with alternative = less gives 1, silently", {
  n <- 3307; m <- 302; r <- 305; t <- 300
  x <- c(rep(1L, m), rep(0L, n - m)); y <- c(rep(1L, t), rep(0L, m - t), rep(1L, r - t), rep(0L, n - m - r + t))
  expect_silent(res <- cfisher.test(x, y, alternative = "less"))
  expect_equal(res$p.value, 1)
})

test_that("a fitted null probability within 1e-8 of 0 or 1 warns", {
  z <- rep(0:1, each = 40); y <- z; x <- rep(c(1L, 0L), 40)
  w <- capture_warnings(a <- cfisher.test(x, y, z))
  expect_match(w, "within 1e-8", all = FALSE)
  # under separation the refit has no maximum likelihood estimate, so the odds ratio is NA
  expect_match(w, "did not converge, which indicates", all = FALSE)
  expect_true(is.na(a$estimate))
  expect_true(all(is.na(a$conf.int)))
  expect_warning(cfisher.null(x, y, z), "within 1e-8")
  set.seed(5)
  expect_silent(cfisher.test(rbinom(60, 1, 0.4), rbinom(60, 1, 0.5), rnorm(60)))
})

test_that("cfisher.test agrees with the reference implementation cb_confint() to 1e-6", {
  # Set CFISHER_CBTEST to the path of the reference file cbtest.R to run this test.
  src <- Sys.getenv("CFISHER_CBTEST", unset = "")
  skip_if_not(nzchar(src) && file.exists(src), "CFISHER_CBTEST (path of cbtest.R) is not set")
  env <- new.env(parent = globalenv())
  sys.source(src, envir = env)
  set.seed(20261102)
  n <- 500; z <- rnorm(n); x <- rbinom(n, 1, plogis(-0.5 + z)); y <- rbinom(n, 1, plogis(-1 + z + log(2) * x))
  xz <- rep(c(0L, 1L, 0L, 1L), c(90, 10, 10, 90)); zz <- rep(c(0, 0, 1, 1), c(90, 10, 10, 90))
  yz <- rbinom(200, 1, plogis(-2 + zz + log(4) * xz))
  cases <- list(list(x = x, y = y, z = z), list(x = xz, y = yz, z = zz),
                list(x = x, y = y, z = NULL))
  for (cs in cases) for (alt in c("two.sided", "greater", "less")) {
    a <- cfisher.test(cs$x, cs$y, cs$z, alternative = alt)
    b <- env$cb_confint(cs$y, cs$x, cs$z, alternative = alt)
    expect_equal(unname(a$estimate), b$estimate, tolerance = 1e-6, label = paste("estimate", alt))
    expect_equal(as.numeric(a$conf.int), as.numeric(b$conf.int), tolerance = 1e-6,
                 label = paste("interval", alt))
  }
  # the tilt mode, with the glm fitted values as the supplied probabilities
  fv <- fitted(glm(y ~ z, family = binomial))
  a <- cfisher.test(x, y, fit = fv); b <- env$cb_confint(y, x, fit = fv)
  expect_equal(unname(a$estimate), b$estimate, tolerance = 1e-6)
  expect_equal(as.numeric(a$conf.int), as.numeric(b$conf.int), tolerance = 1e-6)
})

#################### Review round 2026-11-04 ####################

test_that("the refit by damped Newton agrees with glm.fit on a regular problem", {
  set.seed(4)
  n <- 300; z1 <- rnorm(n); z2 <- rbinom(n, 1, 0.4); x <- rbinom(n, 1, 0.3)
  y <- rbinom(n, 1, plogis(-1 + z1 + z2 + 0.5 * x))
  Z <- cbind(1, z1, z2)
  for (b in c(-3, 0, 2)) {
    g <- cfisher:::.logistic_refit(Z, y, offset = b * x)
    h <- glm.fit(Z, y, offset = b * x, family = binomial())
    expect_true(g$converged)
    expect_equal(unname(g$coef), unname(h$coefficients), tolerance = 1e-8)
  }
})

test_that("D1: a refit that glm.fit got wrong no longer gives a finite limit", {
  # four subjects: the refitted P_b(T = 0) stays above 0.887 at every b (the reviewer's
  # independent solution of the score equations), so the upper limit is Inf
  a <- cfisher.test(c(1, 1, 0, 0), c(0, 0, 0, 1), c(3, 3, -1, 0))
  expect_equal(as.numeric(a$conf.int), c(0, Inf))
  expect_equal(unname(a$estimate), 0)
  # eight subjects: the reviewer's independent upper limit, from a direct optimization of
  # the likelihood and an enumeration of the case pairs, is 90780.69
  x <- c(1, 1, 0, 0, 0, 0, 0, 0); y <- c(0, 0, 0, 0, 0, 1, 1, 0); z <- c(3, 3, -1, 0, 0, 0, 1, -3)
  expect_silent(a <- cfisher.test(x, y, z))
  expect_equal(a$conf.int[2], 90780.69, tolerance = 1e-5)
  # and at that limit the refitted P(T <= 0) is the level 0.025, by enumeration
  d <- refit_dist(x, y, z, log(a$conf.int[2]))
  expect_equal(d$prob[d$support == 0], 0.025, tolerance = 1e-6)
})

test_that("D3: the bracket reaches +-40, and a finite estimate beyond 32 is found", {
  expect_equal(.solve_increasing(function(b) b - 35), 35, tolerance = 1e-8)
  expect_equal(.solve_increasing(function(b) b + 35), -35, tolerance = 1e-8)
  # a root beyond the range: infinite, marked beyond
  r <- .solve_increasing(function(b) b - 41)
  expect_equal(as.numeric(r), Inf); expect_true(attr(r, "beyond"))
  r <- .solve_increasing(function(b) b + 41)
  expect_equal(as.numeric(r), -Inf); expect_true(attr(r, "beyond"))
  # an asymptote, and a value of exactly 0 at the boundary, are not taken as roots (code
  # review of the fixes, finding 1)
  r <- .solve_increasing(function(b) -exp(-b))
  expect_equal(as.numeric(r), Inf); expect_false(attr(r, "beyond"))
  r <- .solve_increasing(function(b) min(b - 40, 0))
  expect_equal(as.numeric(r), Inf); expect_false(attr(r, "beyond"))
  # tilt: the weights of T = 0, 1, 2 are 1, 4w, w^2 with w = exp(b - 34), and t = 1, so the
  # estimate is exp(34) and each limit solves a quadratic in w
  x <- c(1, 1, 0, 0); y <- c(1, 0, 1, 0)
  expect_warning(a <- cfisher.test(x, y, fit = c(rep(plogis(-34), 2), 0.5, 0.5)), "within 1e-8")
  w_lo <- (-3.9 + sqrt(3.9^2 + 4 * 0.975 * 0.025)) / (2 * 0.975)   # P(T >= 1) = 0.025
  w_hi <- (3.9 + sqrt(3.9^2 + 4 * 0.025 * 0.975)) / (2 * 0.025)    # P(T <= 1) = 0.025
  expect_equal(unname(a$estimate), exp(34), tolerance = 1e-8)
  expect_equal(as.numeric(a$conf.int), c(w_lo, w_hi) * exp(34), tolerance = 1e-8)
  # an estimate beyond the search range is NA with a warning, not Inf
  w <- capture_warnings(a <- cfisher.test(x, y, fit = c(rep(plogis(-45), 2), 0.5, 0.5)))
  expect_match(w, "within 1e-8", all = FALSE)               # the near-0 probabilities
  expect_match(w, "stays below the observed count", all = FALSE)
  expect_true(is.na(a$estimate))
})

test_that("D4: an exposure in the span of the confounders gives a warning and NA", {
  x <- c(0, 0, 1, 1); y <- c(0, 1, 0, 1)
  expect_warning(a <- cfisher.test(x, y, z = x), "not identified")
  expect_true(is.na(a$estimate))
  expect_true(all(is.na(a$conf.int)))
  # a combination of confounders, and a glm fit that contains the exposure's span
  set.seed(7)
  n <- 60; z1 <- rnorm(n); x <- rbinom(n, 1, 0.4); y <- rbinom(n, 1, 0.5)
  expect_warning(cfisher.test(x, y, z = cbind(z1, x - 2 * z1)), "not identified")
  expect_warning(cfisher.test(x, y, fit = glm(y ~ I(x + z1) + z1, family = binomial)),
                 "not identified")
  # a redundant confounder column alone does not trigger it, and changes nothing
  a <- cfisher.test(x, y, z = z1)
  expect_silent(b <- cfisher.test(x, y, z = cbind(z1, 2 * z1)))
  expect_equal(a$p.value, b$p.value, tolerance = 1e-10)
  expect_equal(as.numeric(a$conf.int), as.numeric(b$conf.int), tolerance = 1e-8)
})

test_that("D5: a glm with prior weights or an offset is refused for the odds ratio", {
  x <- c(0, 0, 0, 0, 1, 1, 1, 1); y <- c(0, 0, 0, 1, 0, 1, 1, 1); o <- c(-2, -1, 0, 1, 0, 1, 2, 3)
  f <- glm(y ~ 1 + offset(o), family = binomial())
  expect_error(cfisher.test(x, y, fit = f), "offset")
  expect_equal(cfisher.test(x, y, fit = f, conf.int = FALSE)$p.value, 1)  # the p-value still works
  fw <- suppressWarnings(glm(y ~ 1, family = binomial(), weights = rep(1:2, 4)))
  expect_error(cfisher.test(x, y, fit = fw), "prior weights")
  expect_silent(cfisher.test(x, y, fit = fw, conf.int = FALSE))
})

test_that("S1: the crossing check flags a function that crosses its target twice", {
  f <- function(b) (b + 1) * (b - 1) * (b - 3)       # increasing at -1 and 3, decreasing at 1
  expect_true(.second_crossing(f, 3, c(-4, -2, 0, 2, 4)))
  expect_false(.second_crossing(function(b) b - 3, 3, c(-4, -2, 0, 2, 4)))
  expect_false(.second_crossing(function(b) -1, Inf, c(-40, 0, 40)))
  expect_true(.second_crossing(function(b) if (b == 0) 1 else -1, Inf, c(-40, 0, 40)))
})

# ---- Code review of the fixes (2026-11-04): the reviewer's inputs --------------------

test_that("an infinite estimate for an interior count is Inf, with the lower limit in closed form", {
  # Without confounders the unconditional refit has T ~ Binomial(m, p_e), and the
  # intercept solves t = m p_e + n_u p_u; the lower limit solves P(T >= t) = 0.025 for
  # p_e, which gives p_u and the odds ratio. All the n_u unexposed are controls, so the
  # model with the exposure is separated and the estimate and upper limit are Inf.
  lim <- function(m, t, nu) {
    pe <- uniroot(function(p) pbinom(t - 1, m, p, lower.tail = FALSE) - 0.025,
                  c(1e-9, 1 - 1e-9), tol = 1e-14)$root
    exp(qlogis(pe) - qlogis((t - m * pe) / nu))
  }
  # finding 1: was a finite 2.08e15
  xy <- make_xy(3, 5, 0, 10)
  expect_silent(a <- cfisher.test(xy$x, xy$y, conditional = FALSE))
  expect_equal(unname(a$estimate), Inf)
  expect_equal(as.numeric(a$conf.int), c(lim(8, 3, 10), Inf), tolerance = 1e-7)
  # finding 2: the line search used to fail at log odds ratio 16, giving NA
  xy <- make_xy(1, 7, 0, 50)
  expect_silent(a <- cfisher.test(xy$x, xy$y, conditional = FALSE))
  expect_equal(unname(a$estimate), Inf)
  expect_equal(as.numeric(a$conf.int), c(lim(8, 1, 50), Inf), tolerance = 1e-7)
})

test_that("a conditional estimate whose refitted mean stays below t is NA with a warning", {
  # Separation of the model with the exposure along x - s: as the log odds ratio grows the
  # refitted conditional mean tends to a limit below t, so E_b(T) = t has no root. The
  # unconditional mean tends to t, so that estimate is Inf.
  x <- c(1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0); y <- c(1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0)
  s <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
  expect_warning(a <- cfisher.test(x, y, s), "stays below the observed count")
  expect_true(is.na(a$estimate))
  expect_equal(a$conf.int[2], Inf)
  expect_silent(a <- cfisher.test(x, y, s, conditional = FALSE))
  expect_equal(unname(a$estimate), Inf)
})

test_that("the unconditional estimate is the logistic MLE of the exposure's coefficient", {
  # E_b(T) = t is the profile score equation of the model with the exposure (finding 1)
  set.seed(11)
  n <- 150; z <- rnorm(n); x <- rbinom(n, 1, plogis(z)); y <- rbinom(n, 1, plogis(-0.5 + z + 0.7 * x))
  g <- glm(y ~ z + x, family = binomial, control = glm.control(epsilon = 1e-14))
  a <- cfisher.test(x, y, z, conditional = FALSE)
  expect_equal(unname(a$estimate), unname(exp(coef(g)["x"])), tolerance = 1e-8)
  # and .full_mle() returns it directly, on the orthonormal basis of the null design
  d <- list(x = x, y = y, n = n)
  expect_equal(.full_mle(qr.Q(qr(cbind(1, z))), d, 1), unname(coef(g)["x"]), tolerance = 1e-8)
})

test_that("the odds ratio does not depend on the scale or location of a confounder", {
  # finding 3: solve() failed on the raw design for z * 1e8 and 1e6 + 1e3 * z
  set.seed(1)
  n <- 300; z <- rnorm(n); x <- rbinom(n, 1, plogis(-1 + z)); y <- rbinom(n, 1, plogis(-0.5 + z + 0.5 * x))
  a <- cfisher.test(x, y, z)
  for (zz in list(z * 1e8, z * 1e-8, 1e6 + 1e3 * z)) {
    expect_silent(b <- cfisher.test(x, y, zz))
    expect_equal(unname(b$estimate), unname(a$estimate), tolerance = 1e-6)
    expect_equal(as.numeric(b$conf.int), as.numeric(a$conf.int), tolerance = 1e-6)
  }
})

test_that("a refit pinned by fitted probabilities near 0 converges, and the limit is right", {
  # finding 4: the step test failed along a flat direction at log odds ratio -40
  x <- c(1L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L)
  y <- c(1L, 1L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L)
  Z <- cbind(c(0.0346, 0.4217, 0.6709, 1.203, 0.154, 1.475, 0.8507, 0.0905, -1.9904, -0.3456,
               -1.355, -0.8079), c(0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1))
  expect_silent(a <- cfisher.test(x, y, Z, conditional = FALSE, alternative = "less"))
  expect_equal(unname(a$estimate), 0)
  # at the upper limit, P(T <= t) = 0.05 by enumeration, with an independent refit by glm.fit
  b <- log(a$conf.int[2])
  h <- glm.fit(cbind(1, Z), y, offset = b * x, family = binomial(),
               control = glm.control(epsilon = 1e-14, maxit = 2000))
  e <- enumerate_null(h$fitted.values, which(x == 1))
  expect_equal(sum(e$prob[e$support <= sum(y[x == 1])]), 0.05, tolerance = 1e-7)
})

test_that("the p-value of a glm fit does not need its model matrix", {
  # finding 5: the identification check used to rebuild the model matrix
  set.seed(3)
  df <- data.frame(z1 = rnorm(80), z2 = rnorm(80)); y <- rbinom(80, 1, plogis(df$z1))
  x <- rbinom(80, 1, 0.3)
  f <- glm(y ~ z1 + z2, family = binomial, data = df, model = FALSE)
  p <- cfisher.test(x, y, fit = f, conf.int = FALSE)$p.value
  rm(df)
  expect_equal(cfisher.test(x, y, fit = f, conf.int = FALSE)$p.value, p)
  expect_error(cfisher.test(x, y, fit = f), "could not be rebuilt")
})

test_that("the null log-odds of a logit glm are its linear predictor", {
  # finding 6: qlogis(fitted) loses the complement of a fitted value near 1
  zz <- c(-3, -2, -1, 0, 1, 2, 3, 4, 40); y <- c(0, 0, 1, 0, 1, 1, 0, 1, 1)
  x <- c(0, 1, 0, 1, 0, 1, 1, 0, 1)
  f <- glm(y ~ zz, family = binomial)
  d <- suppressWarnings(.prepare(x, y, NULL, f))
  expect_identical(d$eta, unname(f$linear.predictors))
  expect_gt(max(abs(qlogis(fitted(f)) - f$linear.predictors)), 1e-10)   # what it avoids
})

# ---- Second code review of the fixes (2026-11-04) ------------------------------------

test_that("a confounder with a large offset is kept by the refit, as by the null fit", {
  set.seed(1)
  n <- 300; u <- rnorm(n); x <- rbinom(n, 1, plogis(-0.5 + 1.5 * u))
  y <- rbinom(n, 1, plogis(-0.5 + 1.5 * u + 0.4 * x))
  a <- cfisher.test(x, y, u)
  for (b in list(cfisher.test(x, y, u + 1e7),
                 cfisher.test(x, y, fit = glm(y ~ I(u + 1e8), family = binomial)))) {
    expect_equal(unname(b$estimate), unname(a$estimate), tolerance = 1e-6)
    expect_equal(as.numeric(b$conf.int), as.numeric(a$conf.int), tolerance = 1e-6)
  }
})

test_that("an outlying confounder value is not taken for separation", {
  set.seed(11)
  n <- 300; u <- rnorm(n); x <- rbinom(n, 1, plogis(-0.5 + u))
  y <- rbinom(n, 1, plogis(-0.5 + u + 0.5 * x)); z <- u; z[1] <- 3000; y[1] <- 1L
  # the fitted log-odds of subject 1 is about 2800, so the null fit warns, but the MLE exists
  a <- suppressWarnings(cfisher.test(x, y, z))
  expect_true(all(is.finite(c(a$estimate, a$conf.int))))
  b <- suppressWarnings(cfisher.test(x, y, z, conditional = FALSE))
  g <- suppressWarnings(glm(y ~ z + x, family = binomial, control = glm.control(epsilon = 1e-14)))
  expect_equal(unname(b$estimate), unname(exp(coef(g)["x"])), tolerance = 1e-6)
})

test_that("the limits follow tail targets near 0 and near 1", {
  # table (10, 5, 4, 12): the exact limits from the log-scale noncentral hypergeometric
  xy <- make_xy(10, 5, 4, 12)
  ref <- hyper_log_reference(31, 15, 14, 10, a = 5e-11)
  a <- cfisher.test(xy$x, xy$y, conf.level = 1 - 1e-10)
  expect_equal(as.numeric(a$conf.int), c(ref$lower, ref$upper), tolerance = 1e-6)
  ref <- hyper_log_reference(31, 15, 14, 10, a = 1 - 1e-9)
  a <- cfisher.test(xy$x, xy$y, conf.level = 1e-9, alternative = "greater")
  expect_equal(a$conf.int[1], ref$lower, tolerance = 1e-6)
  expect_error(cfisher.test(xy$x, xy$y, conf.level = 1 - 1e-13), "between 1e-12")
})

test_that("a stratum without cases makes the odds ratio NA, with a warning", {
  set.seed(7)
  n <- 600; f <- factor(sample(1:20, n, TRUE)); x <- rbinom(n, 1, 0.4)
  y <- rbinom(n, 1, plogis(-1 + rnorm(20, 0, 0.5)[f] + 0.5 * x))
  y[f == "3"] <- 0L                                       # stratum 3 has no cases
  w <- capture_warnings(a <- cfisher.test(x, y, data.frame(f = f)))
  expect_match(w, "did not converge, which indicates", all = FALSE)
  expect_true(is.na(a$estimate))
  expect_true(all(is.na(a$conf.int)))
  # the stratum carries no information on the exposure: without it, the p-value is the
  # same and the odds ratio is computed
  keep <- f != "3"
  b <- cfisher.test(x[keep], y[keep], data.frame(f = droplevels(f[keep])))
  expect_equal(a$p.value, b$p.value, tolerance = 1e-6)
  expect_true(all(is.finite(c(b$estimate, b$conf.int))))
})
