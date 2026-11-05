# cfisher 0.1.0

* First version of the package.
* `cfisher.test()` tests a binary exposure against a binary outcome while controlling
  for confounders. It returns an `htest` object with the p-value, the conditional
  estimate of the odds ratio and its confidence interval, for the conditional null
  (Fisher's exact test controlling for confounders) or the unconditional
  Poisson-binomial null.
* The odds ratio estimate and interval refit the confounder slopes under each hypothesized
  odds ratio (the exposure term is an offset), because tilting the distribution of the
  null model estimates the wrong family when the exposure is common. A supplied vector of
  probabilities is tilted instead, which is exact only for baseline probabilities under the
  alternative. `conf.int = FALSE` omits the estimate as well as the interval.
* Each refit is a damped Newton fit started from the nearest refit already computed, on
  the columns of the design that the null fit keeps. Separation of the outcome by the
  confounders (for example a stratum with no cases), a refit that does not converge, a
  conditional estimate beyond the search range of log odds ratios from -40 to 40, or an
  exposure that is a linear combination of the confounders gives `NA` with a warning. An estimate whose equation is met only in the
  limit is 0 or `Inf`, and an unconditional estimate beyond the search range is the
  maximum likelihood estimate of the logistic model that contains the exposure. A `glm` fit with prior weights or an offset is refused
  when the odds ratio is requested.
* No probability of the null distribution is lost to underflow, also far below the range
  of double precision: the recursion runs in probability scale when a bound on its
  underflow error allows it and in log scale otherwise. The p-values and the tilt are
  computed from the log probabilities. A warning is given when a fitted null probability
  is within 1e-8 of 0 or 1.
* `cfisher.null()` returns the null distribution of the number of exposed cases.
